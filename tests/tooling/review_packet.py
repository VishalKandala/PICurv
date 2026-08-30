#!/usr/bin/env python3
"""Assemble a bounded review packet for one documentation page or invariant contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PAGES_DIR = REPO_ROOT / "docs" / "pages"
CONTRACTS = REPO_ROOT / "tests" / "tooling" / "contract_registry.json"
FRESHNESS = REPO_ROOT / "tests" / "tooling" / "freshness_manifest.json"
FAMILIES = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
SCOPE = REPO_ROOT / "tests" / "tooling" / "capability_scope_records.json"


def resolve_page(name: str) -> Path:
    """!
    @brief Resolve a page argument, accepting a number, stem, or filename.
    @param[in] name Page identifier supplied on the command line.
    @return Path to the page.
    @throws SystemExit when no page matches.
    """
    candidates = sorted(PAGES_DIR.glob("*.md"))
    for page in candidates:
        if name in (page.name, page.stem) or page.stem.startswith(f"{name}_"):
            return page
    print(f"No page matches '{name}'.", file=sys.stderr)
    raise SystemExit(2)


def git_log(paths: list) -> str:
    """!
    @brief Recent history for a set of paths, for spotting what changed underneath a page.
    @param[in] paths Repository-relative paths.
    @return Formatted git log output.
    """
    if not paths:
        return "  (none declared)"
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "log", "--oneline", "-8", "--", *paths],
        capture_output=True, text=True, check=False,
    )
    return "\n".join(f"  {line}" for line in result.stdout.strip().splitlines()) or "  (no history)"


def contracts_for(page_stem: str) -> list:
    """!
    @brief Invariant contracts that own or depend on a page.
    @param[in] page_stem Page identifier without extension.
    @return Matching contract records.
    """
    registry = json.loads(CONTRACTS.read_text(encoding="utf-8"))
    return [
        c
        for c in registry["contracts"]
        if page_stem in c.get("canonical_documentation", []) + c.get("dependent_pages", [])
    ]


def families_for(page_stem: str) -> list:
    """!
    @brief Capability families whose entries live on a page.
    @param[in] page_stem Page identifier without extension.
    @return Matching family records.
    """
    registry = json.loads(FAMILIES.read_text(encoding="utf-8"))
    return [f for f in registry["families"] if f["family_page"] == page_stem]


def scope_records_for(page_stem: str) -> list:
    """!
    @brief Scope records that name a page as a disclosure or activation surface.
    @param[in] page_stem Page identifier without extension.
    @return Matching scope records.
    """
    if not SCOPE.is_file():
        return []
    records = json.loads(SCOPE.read_text(encoding="utf-8")).get("records", [])
    matches = []
    for record in records:
        surfaces = record.get("disclosed_at", []) + record.get("activation_condition", {}).get("surfaces", [])
        if any(page_stem in surface for surface in surfaces):
            matches.append(record)
    return matches


def freshness_for(page_stem: str) -> list:
    """!
    @brief Freshness surfaces that route review to this page.
    @param[in] page_stem Page id.
    @return List of (role, surface) pairs.
    """
    document = json.loads(FRESHNESS.read_text(encoding="utf-8"))
    found = []
    for surface in document["surfaces"]:
        if page_stem in surface.get("owning_pages", []):
            found.append(("owns", surface))
        elif page_stem in surface.get("dependent_pages", []):
            found.append(("depends on", surface))
    return found


def freshness_state(surface: dict) -> str:
    """!
    @brief Whether a surface is current, suspect, or never attested.
    @param[in] surface One freshness manifest entry.
    @return Human-readable state.
    """
    paths = ([surface["artifact"]] if surface["tier"] == "hard"
             else list(surface.get("watched_paths", [])))
    missing = [p for p in paths if not (REPO_ROOT / p).is_file()]
    if missing:
        return f"INTEGRITY FAILURE - missing {', '.join(missing)}"
    attested = surface.get("attested_digest")
    if attested in (None, "unattested"):
        return "never attested"
    return "current" if attested == _freshness_digest(paths) else "SUSPECT"


def _freshness_digest(paths: list) -> str:
    """!
    @brief Digest matching audit_freshness.digest_of.
    @param[in] paths Repository-relative paths.
    @return Hex digest with algorithm prefix.
    """
    accumulator = hashlib.sha256()
    for relative in sorted(paths):
        accumulator.update(relative.encode("utf-8"))
        accumulator.update(b"\0")
        accumulator.update((REPO_ROOT / relative).read_bytes())
        accumulator.update(b"\0")
    return f"sha256:{accumulator.hexdigest()}"


def uncommitted(paths: list) -> str:
    """!
    @brief Current uncommitted status for declared sources.
    @param[in] paths Repository-relative paths.
    @return Formatted status output.
    """
    if not paths:
        return "  (none declared)"
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--", *paths],
        capture_output=True, text=True, check=False,
    )
    body = result.stdout.strip()
    if not body:
        return "  (clean)"
    return "\n".join(f"  {line}" for line in body.splitlines())


def symbols_for(family: dict) -> list:
    """!
    @brief Named symbols a capability family depends on, rather than whole files.
    @param[in] family Family record.
    @return Symbol descriptions.
    """
    out = []
    surface = family.get("public_surface", {})
    if surface.get("symbol"):
        out.append(f"{surface.get('module', '?')}::{surface['symbol']} ({surface.get('kind')})")
    for parity in family.get("parity_sources", []):
        target = parity.get("function") or parity.get("symbol") or parity.get("prefix")
        if target:
            out.append(f"{parity['path']}::{target} ({parity['kind']})")
    return out


def evidence_sources(family: dict) -> list:
    """!
    @brief Declared evidence sources across a family's values.
    @param[in] family Family record.
    @return Sorted source identifiers.
    """
    sources = set()
    for meta in family.get("value_metadata", {}).values():
        for facet_sources in (meta.get("evidence") or {}).values():
            sources.update(facet_sources)
    return sorted(sources)


def make_targets(sources: list) -> list:
    """!
    @brief Make targets named by declared evidence sources.
    @param[in] sources Source identifiers.
    @return Verification commands.
    """
    return [f"make {s.split(':', 1)[1]}" for s in sources if s.startswith("make:")]


def contract_mode(contract_id: str) -> int:
    """!
    @brief Print a review packet for one invariant contract.
    @param[in] contract_id Contract identifier.
    @return Process status code.
    """
    registry = json.loads(CONTRACTS.read_text(encoding="utf-8"))
    contracts = {c["id"]: c for c in registry["contracts"]}
    contract = contracts.get(contract_id)
    if not contract:
        print(f"No contract '{contract_id}'. Known ids:", file=sys.stderr)
        for known in sorted(contracts):
            print(f"  {known}", file=sys.stderr)
        return 2

    print(f"REVIEW PACKET: contract {contract['id']}")
    print(f"  {contract['title']}  [{contract['status']} / {contract['enforcement']}]\n")
    if contract.get("note"):
        print("SCOPE")
        print(f"  {contract['note']}\n")

    print("AUTHORITATIVE SOURCES")
    for source in contract.get("authoritative_sources", []) or ["  (none)"]:
        print(f"  {source}")
    print("")

    print("CANONICAL DOCUMENTATION (owns the narrative)")
    for page in contract.get("canonical_documentation", []) or ["  (none)"]:
        print(f"  docs/pages/{page}.md")
    print("")

    print("DEPENDENT PAGES (review if this contract changes)")
    for page in contract.get("dependent_pages", []) or ["  (none)"]:
        print(f"  docs/pages/{page}.md")
    print("")

    if contract["id"].startswith("run.artifact"):
        topology = json.loads((REPO_ROOT / "tests" / "tooling" / "artifact_topology.json").read_text())
        print("LOGICAL IDENTITIES")
        for artifact in topology["artifacts"]:
            print(f"  {artifact['id']:26} {artifact['path_rule']}")
        for entry in topology.get("resolved_safety_history", {}).get("entries", []):
            print(f"\n  RESOLVED SAFETY DEFECT: {entry['id']} ({entry['severity']}) - "
                  f"{entry['state']}")
            print(f"      residual risk: {entry['resolution']['residual_risk'][:120]}")
        print("")

    surfaces = [s for s in json.loads(FRESHNESS.read_text(encoding="utf-8"))["surfaces"]
                if any(page in contract.get("canonical_documentation", [])
                       for page in s.get("owning_pages", []))]
    if surfaces:
        print("FRESHNESS SURFACES ON THIS CONTRACT'S PAGES")
        for surface in surfaces:
            print(f"  [{surface['tier']}] {surface['id']} - {freshness_state(surface)}")
        print("")

    print("VERIFY WITH")
    checker = contract.get("checker")
    print(f"  python3 {' '.join(checker)}" if checker else "  (no automated checker - human review)")
    print("  make audit-contracts")
    print("")

    print("UNCOMMITTED CHANGES TO DECLARED SOURCES")
    print(uncommitted(contract.get("authoritative_sources", [])))
    print("")
    print("RECENT HISTORY")
    print(git_log(contract.get("authoritative_sources", [])))
    return 0


def main() -> int:
    """!
    @brief Print a deterministic review packet for a page or an invariant contract.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(
        description="Assemble everything a reviewer needs to check one page or contract."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("page", nargs="?", help="Page number, stem, or filename (for example 44).")
    group.add_argument("--contract", help="Invariant contract id (for example run.artifact_topology).")
    args = parser.parse_args()

    if args.contract:
        return contract_mode(args.contract)

    page = resolve_page(args.page)
    stem = page.stem
    contracts = contracts_for(stem)
    families = families_for(stem)
    records = scope_records_for(stem)

    sources: list = []
    symbols: list = []
    evidence: list = []
    for contract in contracts:
        sources.extend(contract.get("authoritative_sources", []))
    for family in families:
        surface = family.get("public_surface", {})
        if surface.get("module"):
            sources.append(surface["module"].replace(".", "/") + ".py")
        for parity in family.get("parity_sources", []):
            if parity.get("path"):
                sources.append(parity["path"])
        symbols.extend(symbols_for(family))
        evidence.extend(evidence_sources(family))
    sources = sorted({s for s in sources if (REPO_ROOT / s).exists()})
    evidence = sorted(set(evidence))

    print(f"REVIEW PACKET: {page.relative_to(REPO_ROOT)}")
    print(f"  {len(page.read_text(encoding='utf-8').splitlines())} lines\n")

    print("SOURCE SYMBOLS (read these, not whole files)")
    for symbol in symbols or ["  (none declared)"]:
        print(f"  {symbol}" if symbol.strip() != "(none declared)" else symbol)
    print("")

    print("AUTHORITATIVE FILES")
    for source in sources or ["  (none declared)"]:
        print(f"  {source}" if source.strip() != "(none declared)" else source)
    print("")

    if stem == "71_Invariant_Contracts":
        registry = json.loads(CONTRACTS.read_text(encoding="utf-8"))
        print("THIS PAGE DOCUMENTS THE WHOLE REGISTRY")
        for record in sorted(registry["contracts"], key=lambda c: (c["status"], c["id"])):
            print(f"  [{record['status']:8}] {record['id']}")
        print("")

    print("INVARIANT CONTRACTS")
    for contract in contracts:
        role = "owns" if stem in contract.get("canonical_documentation", []) else "depends on"
        print(f"  [{contract['status']}] {contract['id']} ({role} this page)")
        print(f"      inspect: make review-packet CONTRACT={contract['id']}")
    if not contracts:
        print("  (none)")
    print("")

    print("CAPABILITY FAMILIES ON THIS PAGE")
    for family in families:
        mode = "enforced" if family.get("coverage_enforced") else "advisory"
        print(f"  {family['id']} ({mode}) - {len(family.get('value_metadata', {}))} declared values")
    if not families:
        print("  (none)")
    print("")

    if evidence:
        print("DECLARED EVIDENCE SOURCES")
        for source in evidence:
            print(f"  {source}")
        print("")

    if records:
        print("SCOPE RECORDS NAMING THIS PAGE")
        for record in records:
            print(f"  [{record.get('status')}] {record['id']} ({record.get('publication_policy')})")
        print("")

    surfaces = freshness_for(stem)
    print("FRESHNESS SURFACES ROUTING REVIEW HERE")
    for role, surface in surfaces:
        state = freshness_state(surface)
        blocking = surface.get("enforcement") or ("blocking" if surface["tier"] == "hard"
                                                  else "report")
        print(f"  [{surface['tier']}/{blocking}] {surface['id']} ({role}) - {state}")
        watched = ([surface["artifact"]] if surface["tier"] == "hard"
                   else surface.get("watched_paths", []))
        for path in watched:
            print(f"      {path}")
        if state != "current":
            print(f"      after review: make attest-freshness ARGS=\"{surface['id']}\"")
    if not surfaces:
        print("  (none - nothing fingerprinted routes review to this page)")
    print("")

    print("UNCOMMITTED CHANGES TO DECLARED SOURCES")
    print(uncommitted(sources))
    print("")

    print("RECENT HISTORY")
    print(git_log(sources))
    print("")

    surfaces = [s for s in json.loads(FRESHNESS.read_text(encoding="utf-8"))["surfaces"]
                if any(page in contract.get("canonical_documentation", [])
                       for page in s.get("owning_pages", []))]
    if surfaces:
        print("FRESHNESS SURFACES ON THIS CONTRACT'S PAGES")
        for surface in surfaces:
            print(f"  [{surface['tier']}] {surface['id']} - {freshness_state(surface)}")
        print("")

    print("VERIFY WITH")
    for command in sorted(set(make_targets(evidence))):
        print(f"  {command}")
    print("  make audit-capability")
    print("  make audit-contracts")
    print("  make preview-docs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
