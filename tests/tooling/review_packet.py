#!/usr/bin/env python3
"""Assemble bounded review packets from PICurv's declared routing registries."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import subprocess
import sys
from functools import lru_cache
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PAGES_DIR = REPO_ROOT / "docs" / "pages"
TOOLING_DIR = Path(__file__).resolve().parent
if str(TOOLING_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLING_DIR))

from generate_xref_index import (  # noqa: E402
    SCHEMA_VERSION as XREF_SCHEMA_VERSION,
    digest_files as xref_digest_files,
    file_digest as xref_file_digest,
    source_files as xref_source_files,
)


CONTRACTS = TOOLING_DIR / "contract_registry.json"
FRESHNESS = TOOLING_DIR / "freshness_manifest.json"
FAMILIES = TOOLING_DIR / "capability_families.json"
SUBSYSTEMS = TOOLING_DIR / "subsystem_records.json"
SCOPE = TOOLING_DIR / "capability_scope_records.json"
XREF = REPO_ROOT / "docs_build" / "xref.json"
XREF_CAVEAT = (
    "Doxygen provides source-reference edges, not a semantic call graph. Registry tables, "
    "callbacks, function pointers, macros, PETSc dispatch, and runtime-selected paths may "
    "require following intermediate symbols. No direct edge does not prove that a symbol is unused."
)
COMPLETE = "ROUTE: complete over declared registry data (not over code behavior)"
INCOMPLETE = "ROUTE: incomplete — declared identifier has unresolved routing"
UNAVAILABLE = "ROUTE: unavailable — required environment or metadata is absent"


def records(path: Path, key: str) -> list[dict]:
    """!
    @brief Read one registry list.
    @param[in] path Registry JSON path.
    @param[in] key Top-level list key.
    @return Registry records.
    """

    return json.loads(path.read_text(encoding="utf-8"))[key]


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
    print(f"ROUTE: invalid — no page matches '{name}'.", file=sys.stderr)
    print("Known pages:", file=sys.stderr)
    for page in candidates:
        print(f"  {page.stem}", file=sys.stderr)
    raise SystemExit(2)


def source_path(source: str) -> str:
    """!
    @brief Extract a repository path from a path-or-symbol source declaration.
    @param[in] source Declared source string.
    @return Repository-relative path portion.
    """

    return source.split("::", 1)[0]


def module_path(module: str) -> str:
    """!
    @brief Convert a Python module identifier to its repository path.
    @param[in] module Dotted Python module.
    @return Repository-relative Python path.
    """

    return module.replace(".", "/") + ".py"


@lru_cache(maxsize=1)
def page_aliases() -> dict[str, str]:
    """!
    @brief Map both filenames and Doxygen page identifiers to published page stems.
    @return Alias-to-page-stem mapping.
    """

    aliases: dict[str, str] = {}
    for page in PAGES_DIR.glob("*.md"):
        aliases[page.stem] = page.stem
        match = re.search(r"(?:^|\n)\s*@page\s+(\S+)", page.read_text(encoding="utf-8"))
        if match:
            aliases[match.group(1)] = page.stem
    return aliases


def canonical_page(reference: str) -> str:
    """!
    @brief Resolve a registry page reference to its filename stem when declared by Doxygen id.
    @param[in] reference Page stem or Doxygen page identifier.
    @return Canonical filename stem, or the unresolved input.
    """

    return page_aliases().get(reference, reference)


def git_log(paths: list[str]) -> str:
    """!
    @brief Recent history for a set of paths, for spotting what changed underneath a route.
    @param[in] paths Repository-relative paths.
    @return Formatted git log output.
    """

    if not paths:
        return "  (none declared)"
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "log", "--oneline", "-8", "--", *paths],
        capture_output=True,
        text=True,
        check=False,
    )
    return "\n".join(f"  {line}" for line in result.stdout.strip().splitlines()) or "  (no history)"


def contracts_for(page_stem: str) -> list[dict]:
    """!
    @brief Invariant contracts that own or depend on a page.
    @param[in] page_stem Page identifier without extension.
    @return Matching contract records.
    """

    return [
        contract
        for contract in records(CONTRACTS, "contracts")
        if page_stem
        in contract.get("canonical_documentation", []) + contract.get("dependent_pages", [])
    ]


def families_for(page_stem: str) -> list[dict]:
    """!
    @brief Capability families whose entries live on a page.
    @param[in] page_stem Page identifier without extension.
    @return Matching family records.
    """

    return [family for family in records(FAMILIES, "families") if family["family_page"] == page_stem]


def record_pages(record: dict) -> set[str]:
    """!
    @brief Collect obligation and concern pages from a subsystem record.
    @param[in] record Subsystem lifecycle record.
    @return Page stems named by the record.
    """

    pages: set[str] = set()
    for section in ("obligations", "concerns"):
        for value in record.get(section, {}).values():
            if isinstance(value, dict) and value.get("page"):
                pages.add(canonical_page(value["page"]))
    return pages


def subsystems_for_pages(page_stems: set[str]) -> list[dict]:
    """!
    @brief Find subsystem records routed through any supplied page.
    @param[in] page_stems Page stems to match.
    @return Matching subsystem records.
    """

    return [record for record in records(SUBSYSTEMS, "subsystems") if record_pages(record) & page_stems]


def scope_records_for(page_stem: str) -> list[dict]:
    """!
    @brief Scope records that name a page as a disclosure or activation surface.
    @param[in] page_stem Page identifier without extension.
    @return Matching scope records.
    """

    if not SCOPE.is_file():
        return []
    matches = []
    for record in records(SCOPE, "records"):
        surfaces = record.get("disclosed_at", []) + record.get("activation_condition", {}).get("surfaces", [])
        if any(page_stem in surface for surface in surfaces):
            matches.append(record)
    return matches


def freshness_for(page_stem: str) -> list[tuple[str, dict]]:
    """!
    @brief Freshness surfaces that route review to a page.
    @param[in] page_stem Page identifier.
    @return List of role and surface pairs.
    """

    found: list[tuple[str, dict]] = []
    for surface in records(FRESHNESS, "surfaces"):
        if page_stem in surface.get("owning_pages", []):
            found.append(("owns", surface))
        elif page_stem in surface.get("dependent_pages", []):
            found.append(("depends on", surface))
    return found


def freshness_paths(surface: dict) -> list[str]:
    """!
    @brief Return the artifact or watched paths that define a freshness surface.
    @param[in] surface Freshness manifest record.
    @return Repository-relative paths.
    """

    return [surface["artifact"]] if surface["tier"] == "hard" else list(surface.get("watched_paths", []))


def freshness_state(surface: dict) -> str:
    """!
    @brief Whether a surface is current, suspect, or never attested.
    @param[in] surface One freshness manifest entry.
    @return Human-readable state.
    """

    paths = freshness_paths(surface)
    missing = [path for path in paths if not (REPO_ROOT / path).is_file()]
    if missing:
        return f"INTEGRITY FAILURE - missing {', '.join(missing)}"
    attested = surface.get("attested_digest")
    if attested in (None, "unattested"):
        return "never attested"
    return "current" if attested == _freshness_digest(paths) else "SUSPECT"


def _freshness_digest(paths: list[str]) -> str:
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


def uncommitted(paths: list[str]) -> str:
    """!
    @brief Current uncommitted status for declared sources.
    @param[in] paths Repository-relative paths.
    @return Formatted status output.
    """

    if not paths:
        return "  (none declared)"
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--", *paths],
        capture_output=True,
        text=True,
        check=False,
    )
    body = result.stdout.strip()
    return "\n".join(f"  {line}" for line in body.splitlines()) if body else "  (clean)"


def family_source_specs(family: dict) -> list[tuple[str, str, str]]:
    """!
    @brief Resolve a capability family's declared source paths and symbols.
    @param[in] family Capability family record.
    @return Tuples of path, symbol, and source kind.
    """

    specs: list[tuple[str, str, str]] = []
    surface = family.get("public_surface", {})
    if surface.get("module") and surface.get("symbol"):
        specs.append((module_path(surface["module"]), surface["symbol"], surface.get("kind", "public")))
    for parity in family.get("parity_sources", []):
        symbol = parity.get("function") or parity.get("symbol") or parity.get("prefix")
        if parity.get("path") and symbol:
            specs.append((parity["path"], symbol, parity.get("kind", "parity")))
    return specs


def symbols_for(family: dict) -> list[str]:
    """!
    @brief Named symbols a capability family depends on, rather than whole files.
    @param[in] family Family record.
    @return Symbol descriptions.
    """

    return [f"{path}::{symbol} ({kind})" for path, symbol, kind in family_source_specs(family)]


def evidence_sources(family: dict) -> list[str]:
    """!
    @brief Declared evidence sources across a family's values.
    @param[in] family Family record.
    @return Sorted source identifiers.
    """

    sources: set[str] = set()
    for metadata in family.get("value_metadata", {}).values():
        for facet_sources in (metadata.get("evidence") or {}).values():
            sources.update(facet_sources)
    return sorted(sources)


def make_targets(sources: list[str]) -> list[str]:
    """!
    @brief Make targets named by declared evidence sources.
    @param[in] sources Source identifiers.
    @return Verification commands.
    """

    return [f"make {source.split(':', 1)[1]}" for source in sources if source.startswith("make:")]


def route_status(unresolved: list[str], unavailable: bool = False) -> int:
    """!
    @brief Print the fixed route state and any unresolved declarations.
    @param[in] unresolved Problems found while joining declared registry data.
    @param[in] unavailable Whether required environmental metadata is unavailable.
    @return Route exit status.
    """

    if unavailable:
        print(UNAVAILABLE)
        for problem in unresolved:
            print(f"  {problem}")
        return 3
    if unresolved:
        print(INCOMPLETE)
        for problem in unresolved:
            print(f"  {problem}")
        return 3
    print(COMPLETE)
    return 0


def unknown_identifier(kind: str, identifier: str, known: list[str]) -> int:
    """!
    @brief Report an invalid registry identifier with close matches and known ids.
    @param[in] kind Identifier category.
    @param[in] identifier Unknown value.
    @param[in] known Valid values.
    @return Invalid-argument status.
    """

    print(f"ROUTE: invalid — unknown {kind} '{identifier}'.", file=sys.stderr)
    close = difflib.get_close_matches(identifier, known, n=5)
    if close:
        print(f"Close matches: {', '.join(close)}", file=sys.stderr)
    print("Known ids:", file=sys.stderr)
    for value in known:
        print(f"  {value}", file=sys.stderr)
    return 2


def print_freshness(surfaces: list[tuple[str, dict]], heading: str) -> None:
    """!
    @brief Print freshness routing with current attestation state.
    @param[in] surfaces Role and surface pairs.
    @param[in] heading Section heading.
    @return None.
    """

    print(heading)
    for role, surface in surfaces:
        blocking = surface.get("enforcement") or ("blocking" if surface["tier"] == "hard" else "report")
        state = freshness_state(surface)
        print(f"  [{surface['tier']}/{blocking}] {surface['id']} ({role}) - {state}")
        for path in freshness_paths(surface):
            print(f"      {path}")
        if state != "current":
            print(f"      after review: make attest-freshness ARGS=\"{surface['id']}\"")
    if not surfaces:
        print("  (none - no declared freshness surface routes here)")
    print("")


def load_xref() -> tuple[str, dict | None]:
    """!
    @brief Load the optional cross-reference index only when its dirty-byte stamp is current.
    @return State (`current`, `missing`, or `stale`) and current index when available.
    """

    if not XREF.is_file():
        return "missing", None
    try:
        index = json.loads(XREF.read_text(encoding="utf-8"))
        current_sources = xref_source_files(REPO_ROOT)
        current_digest = xref_digest_files(REPO_ROOT, current_sources)
        doxyfile_digest = xref_file_digest(REPO_ROOT / "docs" / "Doxyfile")
    except (OSError, ValueError, KeyError, TypeError):
        return "stale", None
    if (
        index.get("schema_version") != XREF_SCHEMA_VERSION
        or index.get("source_digest") != current_digest
        or index.get("source_files") != [path.as_posix() for path in current_sources]
        or index.get("doxyfile_digest") != doxyfile_digest
        or not isinstance(index.get("symbols"), dict)
    ):
        return "stale", None
    return "current", index


def xref_matches(index: dict, path: str, symbol: str) -> list[str]:
    """!
    @brief Match a declared path and symbol to Doxygen member identifiers.
    @param[in] index Current cross-reference index.
    @param[in] path Expected definition path, or an empty string.
    @param[in] symbol Expected unqualified or qualified symbol.
    @return Matching member ids.
    """

    matches: list[str] = []
    for refid, record in index["symbols"].items():
        name_match = record.get("name") == symbol or record.get("qualified_name") == symbol
        name_match = name_match or str(record.get("qualified_name", "")).endswith(f"::{symbol}")
        name_match = name_match or str(record.get("qualified_name", "")).endswith(f".{symbol}")
        path_match = not path or record.get("definition", {}).get("path") == path
        if name_match and path_match:
            matches.append(refid)
    return sorted(matches)


def format_xref_node(index: dict, refid: str) -> str:
    """!
    @brief Format one indexed symbol with its definition location.
    @param[in] index Current cross-reference index.
    @param[in] refid Doxygen member id.
    @return Compact symbol and source location.
    """

    record = index["symbols"].get(refid, {})
    location = record.get("definition", {})
    line = f":{location.get('line')}" if location.get("line") else ""
    return f"{record.get('qualified_name', record.get('name', refid))} ({location.get('path', '?')}{line})"


def print_xref(specs: list[tuple[str, str, str]], paths: list[str] | None = None) -> None:
    """!
    @brief Print bounded direct and two-hop Doxygen reference evidence for a route.
    @param[in] specs Declared source path, symbol, and kind tuples.
    @param[in] paths Optional paths whose connected definitions should seed the view.
    @return None.
    """

    print("OPTIONAL SOURCE CROSS-REFERENCES")
    state, index = load_xref()
    if state == "missing":
        print("XREF: not built (optional)")
        print("  Registry routing below is unaffected.")
        print("  Run `make docs-xref` to add bounded source-reference evidence.")
        print(f"  {XREF_CAVEAT}\n")
        return
    if state == "stale" or index is None:
        print("XREF: stale — run make docs-xref")
        print("  Stale edges are not displayed; registry routing remains available.")
        print(f"  {XREF_CAVEAT}\n")
        return

    print("XREF: current")
    print(f"  {XREF_CAVEAT}")
    selected: list[tuple[str, str]] = []
    for path, symbol, _ in specs:
        for refid in xref_matches(index, path, symbol):
            selected.append((f"{path}::{symbol}", refid))
    if paths:
        path_set = set(paths)
        connected = [
            (f"definitions in {record.get('definition', {}).get('path')}", refid)
            for refid, record in index["symbols"].items()
            if record.get("definition", {}).get("path") in path_set
            and (record.get("incoming") or record.get("outgoing"))
        ]
        selected.extend(sorted(connected, key=lambda item: format_xref_node(index, item[1]))[:8])

    seen: set[str] = set()
    bounded = []
    for label, refid in selected:
        if refid not in seen:
            seen.add(refid)
            bounded.append((label, refid))
    if not bounded:
        print("  (no indexed definition matched the declared symbols or paths)\n")
        return
    for label, refid in bounded[:8]:
        record = index["symbols"].get(refid, {})
        print(f"  {label}")
        print(f"    definition: {format_xref_node(index, refid)}")
        incoming = [edge for edge in record.get("incoming", []) if edge in index["symbols"]]
        outgoing = [edge for edge in record.get("outgoing", []) if edge in index["symbols"]]
        for edge in incoming[:4]:
            print(f"    direct caller: {format_xref_node(index, edge)}")
        for edge in outgoing[:4]:
            print(f"    direct reference: {format_xref_node(index, edge)}")
        second_in = sorted(
            {
                second
                for edge in incoming[:8]
                for second in index["symbols"].get(edge, {}).get("incoming", [])
                if second in index["symbols"] and second != refid
            }
        )
        second_out = sorted(
            {
                second
                for edge in outgoing[:8]
                for second in index["symbols"].get(edge, {}).get("outgoing", [])
                if second in index["symbols"] and second != refid
            }
        )
        for edge in second_in[:3]:
            print(f"    two-hop caller: {format_xref_node(index, edge)}")
        for edge in second_out[:3]:
            print(f"    two-hop reference: {format_xref_node(index, edge)}")
        if len(incoming) > 4 or len(outgoing) > 4 or len(second_in) > 3 or len(second_out) > 3:
            print("    (additional edges omitted; inspect docs_build/xref.json for this symbol)")
    print("")


def contract_mode(contract_id: str) -> int:
    """!
    @brief Print a review packet for one invariant contract.
    @param[in] contract_id Contract identifier.
    @return Process status code.
    """

    known = {contract["id"]: contract for contract in records(CONTRACTS, "contracts")}
    contract = known.get(contract_id)
    if not contract:
        return unknown_identifier("contract", contract_id, sorted(known))
    sources = [source_path(source) for source in contract.get("authoritative_sources", [])]
    pages = [
        canonical_page(page)
        for page in contract.get("canonical_documentation", []) + contract.get("dependent_pages", [])
    ]
    unresolved = [f"missing source {path}" for path in sources if not (REPO_ROOT / path).exists()]
    unresolved.extend(
        f"missing page docs/pages/{page}.md" for page in pages if not (PAGES_DIR / f"{page}.md").is_file()
    )

    print(f"REVIEW PACKET: contract {contract['id']}")
    status = route_status(unresolved)
    print(f"  {contract['title']}  [{contract['status']} / {contract['enforcement']}]\n")
    if contract.get("note"):
        print("SCOPE")
        print(f"  {contract['note']}\n")
    print("AUTHORITATIVE SOURCES")
    for source in contract.get("authoritative_sources", []) or ["(none)"]:
        print(f"  {source}")
    print("\nCANONICAL DOCUMENTATION (owns the narrative)")
    for page in contract.get("canonical_documentation", []) or ["(none)"]:
        print(f"  docs/pages/{page}.md")
    print("\nDEPENDENT PAGES (review if this contract changes)")
    for page in contract.get("dependent_pages", []) or ["(none)"]:
        print(f"  docs/pages/{page}.md")
    print("")

    if contract["id"].startswith("run.artifact"):
        topology = json.loads((TOOLING_DIR / "artifact_topology.json").read_text(encoding="utf-8"))
        print("LOGICAL IDENTITIES")
        for artifact in topology["artifacts"]:
            print(f"  {artifact['id']:26} {artifact['path_rule']}")
        for entry in topology.get("resolved_safety_history", {}).get("entries", []):
            print(f"\n  RESOLVED SAFETY DEFECT: {entry['id']} ({entry['severity']}) - {entry['state']}")
            print(f"      residual risk: {entry['resolution']['residual_risk'][:120]}")
        print("")

    scoped = []
    for surface in records(FRESHNESS, "surfaces"):
        if set(surface.get("owning_pages", [])) & set(contract.get("canonical_documentation", [])):
            scoped.append(("owns", surface))
    if scoped:
        print_freshness(scoped, "FRESHNESS SURFACES ON THIS CONTRACT'S PAGES")
    print_xref([], sources)
    print("VERIFY WITH")
    checker = contract.get("checker")
    print(f"  python3 {' '.join(checker)}" if checker else "  (no automated checker - human review)")
    print("  make audit-contracts\n")
    print("UNCOMMITTED CHANGES TO DECLARED SOURCES")
    print(uncommitted(sources))
    print("\nRECENT HISTORY")
    print(git_log(sources))
    return status


def capability_mode(capability_id: str) -> int:
    """!
    @brief Join one capability family across symbols, pages, contracts, evidence, and subsystems.
    @param[in] capability_id Capability family id.
    @return Process status code.
    """

    known = {family["id"]: family for family in records(FAMILIES, "families")}
    family = known.get(capability_id)
    if not family:
        return unknown_identifier("capability", capability_id, sorted(known))
    specs = family_source_specs(family)
    page = family["family_page"]
    unresolved = [f"missing declared source {path}" for path, _, _ in specs if not (REPO_ROOT / path).is_file()]
    if not specs:
        unresolved.append("no public or parity source symbol is declared")
    if not (PAGES_DIR / f"{page}.md").is_file():
        unresolved.append(f"missing family page docs/pages/{page}.md")
    related_contracts = contracts_for(page)
    related_subsystems = [
        subsystem
        for subsystem in records(SUBSYSTEMS, "subsystems")
        if capability_id in subsystem.get("capability_families", [])
    ]
    evidence = evidence_sources(family)

    print(f"REVIEW PACKET: capability {capability_id}")
    status = route_status(unresolved)
    mode = "enforced" if family.get("coverage_enforced") else "advisory"
    print(f"  {family['title']} [{mode}; {len(family.get('value_metadata', {}))} values]")
    print(f"  selector: {family.get('selector', '(not declared)')}\n")
    print("SOURCE SYMBOLS")
    for symbol in symbols_for(family) or ["(none declared)"]:
        print(f"  {symbol}")
    print("\nFAMILY DOCUMENTATION")
    print(f"  docs/pages/{page}.md")
    print("\nINVARIANT CONTRACTS ON THAT PAGE")
    for contract in related_contracts or [None]:
        print(f"  [{contract['status']}] {contract['id']}" if contract else "  (none)")
    print("\nSUBSYSTEMS USING THIS FAMILY")
    for subsystem in related_subsystems or [None]:
        print(f"  [{subsystem['status']}] {subsystem['id']}" if subsystem else "  (none declared)")
    print("")
    print_freshness(freshness_for(page), "FRESHNESS SURFACES ROUTING TO THE FAMILY PAGE")
    print("DECLARED EVIDENCE SOURCES")
    for source in evidence or ["(none declared)"]:
        print(f"  {source}")
    print("")
    print_xref(specs)
    print("VERIFY WITH")
    for command in sorted(set(make_targets(evidence))):
        print(f"  {command}")
    print("  make audit-capability")
    print(f"  make review-packet PAGE={page}")
    return status


def subsystem_mode(subsystem_id: str) -> int:
    """!
    @brief Join one subsystem across capability families, obligation pages, contracts, and evidence.
    @param[in] subsystem_id Subsystem lifecycle id.
    @return Process status code.
    """

    known = {record["id"]: record for record in records(SUBSYSTEMS, "subsystems")}
    subsystem = known.get(subsystem_id)
    if not subsystem:
        return unknown_identifier("subsystem", subsystem_id, sorted(known))
    family_map = {family["id"]: family for family in records(FAMILIES, "families")}
    family_ids = subsystem.get("capability_families", [])
    pages = record_pages(subsystem)
    unresolved = [f"unknown capability family {family_id}" for family_id in family_ids if family_id not in family_map]
    unresolved.extend(f"missing page docs/pages/{page}.md" for page in pages if not (PAGES_DIR / f"{page}.md").is_file())
    families = [family_map[family_id] for family_id in family_ids if family_id in family_map]
    specs = [spec for family in families for spec in family_source_specs(family)]
    surfaces = [
        ("routes subsystem page", surface)
        for surface in records(FRESHNESS, "surfaces")
        if set(surface.get("owning_pages", []) + surface.get("dependent_pages", [])) & pages
    ]
    related_contracts = [
        contract
        for contract in records(CONTRACTS, "contracts")
        if set(contract.get("canonical_documentation", []) + contract.get("dependent_pages", [])) & pages
    ]
    evidence = sorted({source for family in families for source in evidence_sources(family)})

    print(f"REVIEW PACKET: subsystem {subsystem_id}")
    status = route_status(unresolved)
    print(f"  {subsystem['title']} [{subsystem['status']} / {subsystem.get('visibility', '?')}]\n")
    print("CAPABILITY FAMILIES")
    for family in families or [None]:
        print(f"  {family['id']} -> docs/pages/{family['family_page']}.md" if family else "  (none declared)")
    print("\nOBLIGATION AND CONCERN PAGES")
    for page in sorted(pages) or ["(none declared)"]:
        print(f"  docs/pages/{page}.md")
    print("\nINVARIANT CONTRACTS ON THOSE PAGES")
    for contract in related_contracts or [None]:
        print(f"  [{contract['status']}] {contract['id']}" if contract else "  (none)")
    print("")
    print_freshness(surfaces, "FRESHNESS SURFACES ROUTING TO THOSE PAGES")
    print("DECLARED EVIDENCE AND TEST TARGETS")
    for source in evidence or ["(none declared through capability families)"]:
        print(f"  {source}")
    print("")
    print_xref(specs, [path for _, surface in surfaces for path in freshness_paths(surface)])
    print("VERIFY WITH")
    for command in sorted(set(make_targets(evidence))):
        print(f"  {command}")
    print("  make audit-subsystems")
    return status


def surface_mode(surface_id: str) -> int:
    """!
    @brief Join one freshness surface across watched paths, pages, contracts, families, and subsystems.
    @param[in] surface_id Freshness surface id.
    @return Process status code.
    """

    known = {surface["id"]: surface for surface in records(FRESHNESS, "surfaces")}
    surface = known.get(surface_id)
    if not surface:
        return unknown_identifier("surface", surface_id, sorted(known))
    paths = freshness_paths(surface)
    pages = {
        canonical_page(page)
        for page in surface.get("owning_pages", []) + surface.get("dependent_pages", [])
    }
    unresolved = [f"missing watched path {path}" for path in paths if not (REPO_ROOT / path).is_file()]
    unresolved.extend(f"missing page docs/pages/{page}.md" for page in pages if not (PAGES_DIR / f"{page}.md").is_file())
    families = [family for family in records(FAMILIES, "families") if family["family_page"] in pages]
    related_contracts = [
        contract
        for contract in records(CONTRACTS, "contracts")
        if set(contract.get("canonical_documentation", []) + contract.get("dependent_pages", [])) & pages
    ]
    subsystems = subsystems_for_pages(pages)
    specs = [spec for family in families for spec in family_source_specs(family)]

    print(f"REVIEW PACKET: surface {surface_id}")
    status = route_status(unresolved)
    print(f"  {surface['title']} [{surface['tier']}] - {freshness_state(surface)}\n")
    print("WATCHED PATHS")
    for path in paths or ["(none declared)"]:
        print(f"  {path}")
    print("\nOWNING AND DEPENDENT PAGES")
    for page in sorted(pages) or ["(none declared)"]:
        role = "owns" if page in surface.get("owning_pages", []) else "depends on"
        print(f"  {role}: docs/pages/{page}.md")
    print("\nINVARIANT CONTRACTS")
    for contract in related_contracts or [None]:
        print(f"  [{contract['status']}] {contract['id']}" if contract else "  (none)")
    print("\nCAPABILITY FAMILIES")
    for family in families or [None]:
        print(f"  {family['id']}" if family else "  (none routed through these pages)")
    print("\nSUBSYSTEMS")
    for subsystem in subsystems or [None]:
        print(f"  [{subsystem['status']}] {subsystem['id']}" if subsystem else "  (none routed through these pages)")
    print("")
    print_xref(specs, paths)
    print("VERIFY WITH")
    print(f"  make review-packet PAGE={next(iter(sorted(pages)), '<page>')}")
    print(f"  make attest-freshness ARGS=\"{surface_id}\"  # only after source comparison")
    return status


def changed_paths() -> tuple[dict[str, set[str]], bool]:
    """!
    @brief Collect staged, unstaged, and untracked nonignored paths from Git.
    @return Path-to-state mapping and whether Git metadata was available.
    """

    probe = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "rev-parse", "--is-inside-work-tree"],
        capture_output=True,
        text=True,
        check=False,
    )
    if probe.returncode != 0:
        return {}, False
    commands = {
        "staged": ["git", "-C", str(REPO_ROOT), "diff", "--cached", "--name-only", "-z"],
        "unstaged": ["git", "-C", str(REPO_ROOT), "diff", "--name-only", "-z"],
        "untracked": ["git", "-C", str(REPO_ROOT), "ls-files", "--others", "--exclude-standard", "-z"],
    }
    found: dict[str, set[str]] = {}
    for state, command in commands.items():
        completed = subprocess.run(command, capture_output=True, check=False)
        if completed.returncode != 0:
            return {}, False
        for raw in completed.stdout.split(b"\0"):
            if raw:
                found.setdefault(raw.decode("utf-8", errors="replace"), set()).add(state)
    return found, True


def classify_path(path: str) -> str:
    """!
    @brief Classify a changed path for routing coverage.
    @param[in] path Repository-relative changed path.
    @return Category label excluding routed versus unrouted production state.
    """

    candidate = Path(path)
    if path.startswith("docs/generated/") or path.startswith("docs_build/"):
        return "generated"
    if (
        (candidate.parts and candidate.parts[0] == "src" and candidate.suffix in {".c", ".h"})
        or (candidate.parts and candidate.parts[0] == "include" and candidate.suffix in {".h", ".hpp"})
        or (candidate.parts and candidate.parts[0] == "picurv_cli" and candidate.suffix == ".py")
    ):
        return "production"
    routed_roots = {"docs", "config", "tests", ".github", ".agents", ".claude", "examples"}
    routed_files = {"AGENTS.md", "CLAUDE.md", "CONTRIBUTING.md", "README.md", "Makefile", ".gitignore"}
    if (candidate.parts and candidate.parts[0] in routed_roots) or path in routed_files:
        return "documentation/configuration/test/tooling"
    return "ignored/out of scope"


def routes_for_path(path: str) -> dict[str, list[str]]:
    """!
    @brief Find declared capabilities, surfaces, contracts, and subsystems touching a path.
    @param[in] path Repository-relative production path.
    @return Nonempty routing categories and sorted identifiers.
    """

    family_hits = [
        family["id"]
        for family in records(FAMILIES, "families")
        if path in {spec[0] for spec in family_source_specs(family)}
    ]
    surface_hits = [
        surface["id"] for surface in records(FRESHNESS, "surfaces") if path in freshness_paths(surface)
    ]
    contract_hits = [
        contract["id"]
        for contract in records(CONTRACTS, "contracts")
        if path in {source_path(source) for source in contract.get("authoritative_sources", [])}
    ]
    pages: set[str] = set()
    for family in records(FAMILIES, "families"):
        if family["id"] in family_hits:
            pages.add(family["family_page"])
    for surface in records(FRESHNESS, "surfaces"):
        if surface["id"] in surface_hits:
            pages.update(surface.get("owning_pages", []) + surface.get("dependent_pages", []))
    for contract in records(CONTRACTS, "contracts"):
        if contract["id"] in contract_hits:
            pages.update(contract.get("canonical_documentation", []) + contract.get("dependent_pages", []))
    subsystem_hits = [record["id"] for record in subsystems_for_pages(pages)]
    routes = {
        "capability": sorted(family_hits),
        "surface": sorted(surface_hits),
        "contract": sorted(contract_hits),
        "subsystem": sorted(subsystem_hits),
    }
    return {kind: values for kind, values in routes.items() if values}


def guide_fallback(path: str) -> str:
    """!
    @brief Return the nearest existing directory guide for an unrouted production path.
    @param[in] path Repository-relative path.
    @return Guide path or root instruction fallback.
    """

    current = (REPO_ROOT / path).parent
    while current != REPO_ROOT and REPO_ROOT in current.parents:
        guide = current / "guide.md"
        if guide.is_file():
            return guide.relative_to(REPO_ROOT).as_posix()
        current = current.parent
    return "AGENTS.md"


def changed_mode(value: str) -> int:
    """!
    @brief Report declared routing coverage for current working-tree changes.
    @param[in] value Supported selector, currently `working-tree`.
    @return Process status code.
    """

    if value != "working-tree":
        return unknown_identifier("changed-set", value, ["working-tree"])
    paths, available = changed_paths()
    routed: dict[str, dict[str, list[str]]] = {}
    unrouted: list[str] = []
    categories: dict[str, str] = {}
    for path in sorted(paths):
        category = classify_path(path)
        if category == "production":
            declared = routes_for_path(path)
            if declared:
                category = "routed production"
                routed[path] = declared
            else:
                category = "unrouted production"
                unrouted.append(path)
        categories[path] = category

    print("REVIEW PACKET: changed working-tree")
    problems = [f"unrouted production path: {path}" for path in unrouted]
    status = route_status(problems if available else ["Git working-tree metadata is unavailable"], not available)
    print("  Coverage is advisory and covers staged, unstaged, and untracked nonignored paths.")
    print("  Ignored files are intentionally excluded by Git.\n")
    if not paths and available:
        print("CHANGED PATHS\n  (working tree clean)")
        return status
    print("CHANGED PATHS")
    for path in sorted(paths):
        states = ", ".join(sorted(paths[path]))
        print(f"  [{categories[path]}] {path} ({states})")
        for kind, identifiers in routed.get(path, {}).items():
            print(f"      {kind}: {', '.join(identifiers)}")
        if categories[path] == "unrouted production":
            guide = guide_fallback(path)
            root = Path(path).parts[0]
            print(f"      fallback guide: {guide}")
            print(f"      targeted search: rg -n \"<responsibility-or-symbol>\" {root}")
    return status


def page_mode(page_name: str) -> int:
    """!
    @brief Print the original page-scoped review packet with explicit route status and xrefs.
    @param[in] page_name Page number, stem, or filename.
    @return Process status code.
    """

    page = resolve_page(page_name)
    stem = page.stem
    page_contracts = contracts_for(stem)
    families = families_for(stem)
    scope_records = scope_records_for(stem)
    sources: list[str] = []
    specs: list[tuple[str, str, str]] = []
    evidence: list[str] = []
    for contract in page_contracts:
        sources.extend(source_path(source) for source in contract.get("authoritative_sources", []))
    for family in families:
        family_specs = family_source_specs(family)
        specs.extend(family_specs)
        sources.extend(path for path, _, _ in family_specs)
        evidence.extend(evidence_sources(family))
    sources = sorted(set(sources))
    evidence = sorted(set(evidence))
    unresolved = [f"missing declared source {source}" for source in sources if not (REPO_ROOT / source).exists()]

    print(f"REVIEW PACKET: {page.relative_to(REPO_ROOT)}")
    status = route_status(unresolved)
    print(f"  {len(page.read_text(encoding='utf-8').splitlines())} lines\n")
    print("SOURCE SYMBOLS (read these, not whole files)")
    for symbol in [f"{path}::{symbol} ({kind})" for path, symbol, kind in specs] or ["(none declared)"]:
        print(f"  {symbol}")
    print("\nAUTHORITATIVE FILES")
    for source in sources or ["(none declared)"]:
        print(f"  {source}")
    print("")
    if stem == "71_Invariant_Contracts":
        print("THIS PAGE DOCUMENTS THE WHOLE REGISTRY")
        for contract in sorted(records(CONTRACTS, "contracts"), key=lambda item: (item["status"], item["id"])):
            print(f"  [{contract['status']:8}] {contract['id']}")
        print("")
    print("INVARIANT CONTRACTS")
    for contract in page_contracts:
        role = "owns" if stem in contract.get("canonical_documentation", []) else "depends on"
        print(f"  [{contract['status']}] {contract['id']} ({role} this page)")
        print(f"      inspect: make review-packet CONTRACT={contract['id']}")
    if not page_contracts:
        print("  (none)")
    print("\nCAPABILITY FAMILIES ON THIS PAGE")
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
    if scope_records:
        print("SCOPE RECORDS NAMING THIS PAGE")
        for record in scope_records:
            print(f"  [{record.get('status')}] {record['id']} ({record.get('publication_policy')})")
        print("")
    print_freshness(freshness_for(stem), "FRESHNESS SURFACES ROUTING REVIEW HERE")
    print_xref(specs, sources)
    print("UNCOMMITTED CHANGES TO DECLARED SOURCES")
    print(uncommitted(sources))
    print("\nRECENT HISTORY")
    print(git_log(sources))
    print("\nVERIFY WITH")
    for command in sorted(set(make_targets(evidence))):
        print(f"  {command}")
    print("  make audit-capability")
    print("  make audit-contracts")
    print("  make preview-docs")
    return status


def parse_args(argv: list[str]) -> argparse.Namespace:
    """!
    @brief Parse the mutually exclusive review-packet route selectors.
    @param[in] argv Arguments excluding the executable name.
    @return Parsed command-line namespace.
    """

    parser = argparse.ArgumentParser(
        description="Assemble everything a reviewer needs for one declared PICurv route."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("page", nargs="?", help="Page number, stem, or filename (for example 44).")
    group.add_argument("--contract", help="Invariant contract id.")
    group.add_argument("--capability", help="Capability family id.")
    group.add_argument("--subsystem", help="Subsystem lifecycle id.")
    group.add_argument("--surface", help="Freshness surface id.")
    group.add_argument("--changed", help="Changed-set selector; currently working-tree.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """!
    @brief Dispatch the requested bounded route and preserve differentiated exit states.
    @param[in] argv Optional arguments excluding the executable name.
    @return Zero for complete, two for invalid, and three for incomplete or unavailable.
    """

    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.contract:
        return contract_mode(args.contract)
    if args.capability:
        return capability_mode(args.capability)
    if args.subsystem:
        return subsystem_mode(args.subsystem)
    if args.surface:
        return surface_mode(args.surface)
    if args.changed:
        return changed_mode(args.changed)
    return page_mode(args.page)


if __name__ == "__main__":
    raise SystemExit(main())
