#!/usr/bin/env python3
"""Enforce capability parity across the full source chain and Tier-2 documentation coverage."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
INVENTORY_PATH = REPO_ROOT / "docs" / "generated" / "capability_inventory.json"
GENERATOR = REPO_ROOT / "tests" / "tooling" / "generate_capability_inventory.py"

CANONICAL_FIELDS = (
    "Identity",
    "What it does",
    "When to choose it",
    "Parameters it owns",
    "Interactions",
    "Diagnostics",
    "Evidence",
    "Limitations",
)
ALIAS_FIELDS = ("Identity", "Status", "Migration")


def slug(value: str) -> str:
    """!
    @brief Convert a selector value into the anchor slug its entry must use.
    @param[in] value Public selector value.
    @return Lowercase underscore slug.
    """
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def strip_non_prose(text: str) -> str:
    """!
    @brief Remove fenced code blocks and HTML comments so anchors inside examples do not count.
    @param[in] text Raw page text.
    @return Page text with code fences and comments blanked out.
    """
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    text = re.sub(r"^```.*?^```", "", text, flags=re.S | re.M)
    return text


def check_generated_current() -> list[str]:
    """!
    @brief Verify the committed inventory matches what the sources currently produce.
    @return Violation lines.
    """
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "--check"], capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        detail = result.stderr.strip().replace("\n", "\n      ")
        return [f"generated capability inventory is stale or invalid; run 'make docs-inventory'\n      {detail}"]
    return []


def parity_by_kind(family: dict) -> dict[str, dict]:
    """!
    @brief Index a family's parity records by source kind.
    @param[in] family Inventory record for one family.
    @return Mapping of source kind to its parity record.
    """
    return {record["source"]["kind"]: record for record in family["parity"]}


def check_parity(family: dict) -> list[str]:
    """!
    @brief Verify the public selector set agrees with every declared parity source.

    Every declared source kind is handled explicitly; an unrecognized kind is a
    violation rather than a silent pass, because a source that is registered but
    never compared produces a false assurance of parity.
    @param[in] family Inventory record for one family.
    @return Violation lines.
    """
    violations: list[str] = []
    public = set(family["public_values"])
    records = parity_by_kind(family)
    known = {"c_string_map", "c_token_map", "c_switch", "c_dispatch", "c_enum"}

    for record in family["parity"]:
        kind = record["source"]["kind"]
        if kind not in known:
            violations.append(
                f"{family['id']}: parity source kind '{kind}' is registered but not verified by this audit"
            )

    # Link 1: every public value must be accepted by the C parser, and vice versa.
    accepted: dict[str, str] = {}
    for kind in ("c_string_map", "c_token_map"):
        record = records.get(kind)
        if not record:
            continue
        accepted = record.get("mapping", {})
        path = record["source"]["path"]
        if kind == "c_string_map":
            expected = public
        else:
            expected = {spec.get("maps_to") for spec in family["public_values"].values()}
        for missing in sorted(expected - set(accepted)):
            violations.append(
                f"{family['id']}: '{missing}' is produced by the validator but the C parser "
                f"({path}) does not accept it"
            )
        legacy = set(record["source"].get("legacy_tokens", {}))
        for token in sorted(legacy - set(accepted)):
            violations.append(
                f"{family['id']}: '{token}' is declared a legacy token but the C parser ({path}) "
                f"no longer accepts it; remove the declaration"
            )
        for extra in sorted(set(accepted) - expected - legacy):
            violations.append(
                f"{family['id']}: the C parser ({path}) accepts '{extra}' but the validator "
                f"never produces it; either expose it or remove it"
            )

    # Link 2: every enum a value resolves to must be dispatched at runtime.
    for kind, label in (("c_switch", "factory"), ("c_dispatch", "runtime dispatch")):
        record = records.get(kind)
        if not record:
            continue
        handled = set(record["values"])
        path = record["source"]["path"]
        for token, enum in sorted(accepted.items()):
            if enum not in handled:
                violations.append(
                    f"{family['id']}: token '{token}' resolves to {enum}, which the {label} in "
                    f"{path} does not handle; selecting it would fail at runtime"
                )

    # Link 3: every resolved enum must exist in the declared enum type.
    record = records.get("c_enum")
    if record:
        declared = set(record["values"])
        path = record["source"]["path"]
        for token, enum in sorted(accepted.items()):
            if enum not in declared:
                violations.append(
                    f"{family['id']}: token '{token}' resolves to {enum}, which is not a member of "
                    f"the enum declared in {path}"
                )
    return violations


def check_metadata(family: dict, registry_entry: dict) -> list:
    """!
    @brief Verify declared value metadata is complete, well-typed, and matches the sources.

    @details Fails closed on the status field. A canonical value with no status, or with
             a typo such as "suported", previously passed every check while generation
             quietly defaulted it to supported - so a defective capability could read as
             production-ready through an omission.
    @param[in] family Inventory record for one family.
    @param[in] registry_entry Registry entry carrying value metadata.
    @return Violation lines.
    """
    # No early return on an absent metadata block: a family that exposes values but
    # declares nothing must report every value as undeclared, not pass silently.
    metadata = registry_entry.get("value_metadata", {})
    public = set(family["public_values"])
    violations = []
    for stale in sorted(set(metadata) - public):
        violations.append(
            f"{family['id']}: metadata declares '{stale}', which is no longer a public value; "
            f"remove it from capability_families.json"
        )
    for undeclared in sorted(public - set(metadata)):
        violations.append(
            f"{family['id']}: public value '{undeclared}' has no metadata entry; "
            f"declare its status and whether it is canonical"
        )
    for name, spec in sorted(metadata.items()):
        target = spec.get("alias_of") or spec.get("spelling_of")
        if spec.get("canonical") is False and not target:
            violations.append(
                f"{family['id']}: '{name}' is marked non-canonical but names no alias_of or "
                f"spelling_of target"
            )
        if target and target not in metadata:
            violations.append(f"{family['id']}: '{name}' aliases '{target}', which is not a declared value")

        if spec.get("spelling_of"):
            # A spelling has no lifecycle of its own; it inherits one.
            if "status" in spec:
                violations.append(
                    f"{family['id']}: '{name}' is a spelling of '{spec['spelling_of']}' and must not "
                    f"declare its own status; it inherits one"
                )
            continue

        status = spec.get("status")
        if status is None:
            violations.append(
                f"{family['id']}: '{name}' declares no status; every canonical value and "
                f"deprecated alias must declare one of {list(VALID_STATUSES)}"
            )
            continue
        if status not in VALID_STATUSES:
            violations.append(
                f"{family['id']}: '{name}' has status '{status}', which is not in the closed "
                f"vocabulary {list(VALID_STATUSES)}"
            )
            continue
        if spec.get("alias_of") and status != "deprecated":
            violations.append(
                f"{family['id']}: '{name}' declares alias_of '{spec['alias_of']}' but status "
                f"'{status}'; an alias is by definition deprecated"
            )
        if status in NON_SELECTABLE_STATUSES:
            # Reachability is computed separately: a value the sources declare but no
            # provider can satisfy is already marked latent, which is consistent with a
            # non-selectable lifecycle. The contradiction is a non-selectable status on
            # a value that IS reachable.
            reachable = family["public_values"].get(name, {}).get("reachable", True)
            if reachable:
                violations.append(
                    f"{family['id']}: '{name}' has status '{status}' but is publicly selectable; "
                    f"a {status} capability must not be reachable"
                )
    return violations


def check_coverage(family: dict, registry_entry: dict, all_entries: list = None) -> tuple:
    """!
    @brief Verify every public value has a Tier-2 entry carrying its required fields.
    @param[in] family Inventory record for one family.
    @param[in] registry_entry Registry entry carrying anchor prefix and enforcement flag.
    @param[in] all_entries Every registry entry, so anchors owned by a sibling family on
               the same page are not misreported as stale.
    @return Blocking violations and advisory notes.
    """
    all_entries = all_entries or [registry_entry]
    page = REPO_ROOT / "docs" / "pages" / f"{registry_entry['family_page']}.md"
    if not page.is_file():
        return ([f"{family['id']}: family page {page.name} does not exist"], [])
    prose = strip_non_prose(page.read_text(encoding="utf-8"))
    anchors = set(re.findall(r"^@anchor\s+([A-Za-z0-9_]+)\s*$", prose, re.M))
    prefix = registry_entry["entry_anchor_prefix"]
    metadata = registry_entry.get("value_metadata", {})

    # Slug collisions would make two values share one entry.
    seen: dict[str, str] = {}
    problems: list[str] = []
    for value in sorted(family["public_values"]):
        key = slug(value)
        if key in seen:
            problems.append(
                f"{family['id']}: '{value}' and '{seen[key]}' produce the same anchor slug "
                f"'{prefix}{key}'; entries would collide"
            )
        seen[key] = value

    # Latent values are declared but not selectable, so they owe no Tier-2 entry.
    # Latent values are not selectable, and accepted spellings are mere synonyms of a
    # canonical value; neither owes its own Tier-2 entry.
    selectable = {
        name: spec
        for name, spec in family["public_values"].items()
        if spec.get("reachability") != "latent" and not spec.get("spelling_of")
    }
    expected = {f"{prefix}{slug(v)}": v for v in selectable}
    for anchor, value in sorted(expected.items()):
        if anchor not in anchors:
            problems.append(
                f"{family['id']}: no Tier-2 entry for '{value}' (expected `@anchor {anchor}` in {page.name})"
            )
            continue
        body = entry_body(prose, anchor)
        spec = metadata.get(value, {})
        required = ALIAS_FIELDS if spec.get("alias_of") else CANONICAL_FIELDS
        for field in required:
            if not re.search(rf"\*\*{re.escape(field)}", body):
                problems.append(
                    f"{family['id']}: entry for '{value}' is missing the **{field}** part"
                )

    # Stale entries: an anchor with our prefix that no current value claims. A page may
    # host several families, and one prefix can be a prefix of another (`p08_cap_` vs
    # `p08_cap_conv_`), so an anchor belongs to this family only when no more specific
    # registered prefix on the same page also matches it.
    others = [
        other["entry_anchor_prefix"]
        for other in all_entries
        if other is not registry_entry
        and other["family_page"] == registry_entry["family_page"]
        and other["entry_anchor_prefix"].startswith(prefix)
        and other["entry_anchor_prefix"] != prefix
    ]
    for anchor in sorted(a for a in anchors if a.startswith(prefix)):
        if any(anchor.startswith(other) for other in others):
            continue
        if anchor not in expected:
            problems.append(
                f"{family['id']}: {page.name} carries a stale entry `@anchor {anchor}` that no "
                f"current public value claims; remove it or restore the capability"
            )

    if registry_entry.get("coverage_enforced"):
        return (problems, [])
    return ([], problems)


def entry_body(prose: str, anchor: str) -> str:
    """!
    @brief Return the text of one capability entry, from its anchor to the next entry.
    @param[in] prose Page text with code fences removed.
    @param[in] anchor Entry anchor name.
    @return Entry body text.
    """
    match = re.search(rf"^@anchor\s+{re.escape(anchor)}\s*$", prose, re.M)
    if not match:
        return ""
    rest = prose[match.end() :]
    nxt = re.search(r"^@(?:subsection|section)\s", rest, re.M)
    return rest[: nxt.start()] if nxt else rest


SCOPE_RECORDS_PATH = REPO_ROOT / "tests" / "tooling" / "capability_scope_records.json"


def source_exists(identifier: str) -> bool:
    """!
    @brief Check that one evidence source identifier names something that exists.
    @param[in] identifier Source identifier such as `make:unit-boundaries`.
    @return True when the named artifact exists.
    """
    kind, _, value = identifier.partition(":")
    if kind == "make":
        makefile = (REPO_ROOT / "Makefile").read_text(encoding="utf-8")
        return re.search(rf"^{re.escape(value)}\s*:", makefile, re.M) is not None
    if kind in {"example", "file"}:
        return (REPO_ROOT / value).exists()
    return False


def source_token(identifier: str) -> str:
    """!
    @brief The human-readable token a capability entry must cite for a declared source.
    @param[in] identifier Source identifier such as `make:unit-boundaries`.
    @return The bare target, example directory, or file path.
    """
    _, _, value = identifier.partition(":")
    return value


def check_evidence(family: dict, registry_entry: dict, facets: dict) -> list[str]:
    """!
    @brief Verify declared evidence sources exist and are cited by the capability entry.

    This checks *correspondence*, not scientific validity. It establishes that a
    declared source names something real and that the entry a reader sees cites the
    same source the registry does. It cannot establish that the source actually
    demonstrates the claimed result - that remains a human review judgement, which is
    why these are described as declared evidence sources rather than verified evidence.
    @param[in] family Inventory record for one family.
    @param[in] registry_entry Registry entry carrying value metadata.
    @param[in] facets The project-wide facet vocabulary.
    @return Violation lines.
    """
    page = REPO_ROOT / "docs" / "pages" / f"{registry_entry['family_page']}.md"
    prose = strip_non_prose(page.read_text(encoding="utf-8")) if page.is_file() else ""
    violations: list[str] = []
    for name, meta in sorted(registry_entry.get("value_metadata", {}).items()):
        if meta.get("spelling_of") or meta.get("alias_of"):
            if "evidence" in meta:
                violations.append(
                    f"{family['id']}: '{name}' is a spelling/alias and must not carry its own evidence"
                )
            continue
        evidence = meta.get("evidence")
        if evidence is None:
            violations.append(f"{family['id']}: '{name}' declares no evidence mapping (use {{}} for none)")
            continue
        for facet, sources in sorted(evidence.items()):
            if facet not in facets:
                violations.append(
                    f"{family['id']}: '{name}' claims unknown evidence facet '{facet}'; "
                    f"vocabulary is {sorted(facets)}"
                )
                continue
            if not sources:
                violations.append(f"{family['id']}: '{name}' claims '{facet}' with no source identifier")
            for source in sources:
                if not source_exists(source):
                    violations.append(
                        f"{family['id']}: '{name}' cites '{source}' for '{facet}', which does not exist"
                    )
        # Correspondence: the entry a reader sees must cite the sources the registry
        # declares, or the two records drift apart silently.
        if evidence:
            anchor = registry_entry["entry_anchor_prefix"] + slug(name)
            body = entry_body(prose, anchor)
            if body:
                section = body.split("**Evidence.**", 1)
                evidence_text = section[1] if len(section) > 1 else ""
                for facet, sources in sorted(evidence.items()):
                    for source in sources:
                        token = source_token(source)
                        if token not in evidence_text:
                            violations.append(
                                f"{family['id']}: '{name}' declares '{source}' for '{facet}', "
                                f"but its entry's Evidence part does not cite '{token}'"
                            )
        if meta.get("status") == "supported" and not evidence:
            violations.append(
                f"{family['id']}: '{name}' is marked supported but claims no evidence; "
                f"either record a facet or lower the status"
            )
    return violations


def check_scope_records() -> list[str]:
    """!
    @brief Verify every known-defective scope record is disclosed where it claims to be.

    A safety valve that nothing reads is not a safety valve. This makes the scope
    records load-bearing: a record marked known-defective must either name the pages
    that disclose it and have them actually say so, or be resolved.
    @return Violation lines.
    """
    if not SCOPE_RECORDS_PATH.is_file():
        return []
    records = json.loads(SCOPE_RECORDS_PATH.read_text(encoding="utf-8")).get("records", [])
    violations: list[str] = []
    for record in records:
        if record.get("status") != "known-defective":
            continue
        policy = record.get("publication_policy")
        if policy == "disclose-now":
            surfaces = record.get("disclosed_at", [])
            if not surfaces:
                violations.append(
                    f"scope record '{record['id']}' is known-defective with policy disclose-now "
                    f"but names no disclosure surfaces"
                )
            for surface in surfaces:
                path = REPO_ROOT / surface.split(" ")[0]
                if not path.is_file():
                    violations.append(f"scope record '{record['id']}' names missing surface {path.name}")
                elif "known-defective" not in path.read_text(encoding="utf-8"):
                    violations.append(
                        f"scope record '{record['id']}' claims disclosure in {path.name}, "
                        f"but that page contains no known-defective disclosure"
                    )
        elif policy == "scope-only":
            condition = record.get("activation_condition", {})
            if not condition.get("surfaces"):
                violations.append(
                    f"scope record '{record['id']}' is scope-only but names no activation surfaces"
                )
        else:
            violations.append(
                f"scope record '{record['id']}' has unknown publication_policy '{policy}'"
            )
    return violations


# Closed lifecycle vocabulary. A status outside this set is a violation, not an
# unrecognized value that silently loses its obligations.
VALID_STATUSES = (
    "supported", "experimental", "known-defective", "deprecated",
    "planned", "internal", "removed",
)

# Statuses that describe something not publicly selectable. A value carrying one of
# these must not appear in the public inventory at all.
NON_SELECTABLE_STATUSES = ("planned", "internal", "removed")

LIFECYCLE_REQUIREMENTS = {
    "supported": ("Evidence",),
    "experimental": ("Limitations",),
    "known-defective": ("Limitations",),
    "deprecated": ("Migration",),
}


def check_lifecycle_requirements(family: dict, registry_entry: dict) -> list:
    """!
    @brief Enforce the documentation each lifecycle status owes.

    @details Requirements grow with the status a capability claims: `supported` owes
             evidence, `experimental` and `known-defective` owe stated limitations, and
             `deprecated` owes migration guidance. The gate fires when a value claims a
             status, which is the moment the claim becomes readable.
    @param[in] family Inventory record for one family.
    @param[in] registry_entry Registry entry carrying value metadata and the page.
    @return Violation lines.
    """
    page = REPO_ROOT / "docs" / "pages" / f"{registry_entry['family_page']}.md"
    if not page.is_file():
        return []
    prose = strip_non_prose(page.read_text(encoding="utf-8"))
    prefix = registry_entry["entry_anchor_prefix"]
    violations: list = []
    for name, meta in sorted(registry_entry.get("value_metadata", {}).items()):
        if meta.get("spelling_of"):
            continue
        # Use the effective status the generated inventory displays, so the obligation
        # and the rendered claim can never disagree.
        effective = family["public_values"].get(name, {}).get("status") or meta.get("status")
        required = LIFECYCLE_REQUIREMENTS.get(effective)
        if not required:
            continue
        body = entry_body(prose, prefix + slug(name))
        if not body:
            continue  # coverage check already reports a missing entry
        for field in required:
            if not re.search(rf"\*\*{re.escape(field)}", body):
                violations.append(
                    f"{family['id']}: '{name}' is {effective} but its entry has no "
                    f"**{field}** part; that status owes it"
                )
    return violations


def main() -> int:
    """!
    @brief Fail on capability parity breaks, metadata drift, or missing entries.
    @return Process status code.
    """
    registry_doc = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    registry = {f["id"]: f for f in registry_doc["families"]}
    facets = registry_doc["evidence_facets"]
    if not INVENTORY_PATH.is_file():
        print("Capability inventory has not been generated; run 'make docs-inventory'.", file=sys.stderr)
        return 1
    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))

    blocking = check_generated_current() + check_scope_records()
    advisory: list[str] = []
    for family in inventory:
        entry = registry[family["id"]]
        blocking += check_parity(family)
        blocking += check_metadata(family, entry)
        blocking += check_evidence(family, entry, facets)
        blocking += check_lifecycle_requirements(family, entry)
        family_blocking, family_advisory = check_coverage(family, entry, list(registry.values()))
        blocking += family_blocking
        advisory += family_advisory

    if advisory:
        print("Capability documentation coverage (advisory, backfill pending):")
        for note in advisory:
            print(f"  {note}")
        print("")

    if blocking:
        print("Capability parity/coverage violations:", file=sys.stderr)
        for violation in blocking:
            print(f"  {violation}", file=sys.stderr)
        return 1

    enforced = sorted(f["id"] for f in registry.values() if f.get("coverage_enforced"))
    pending = sorted(f["id"] for f in registry.values() if not f.get("coverage_enforced"))
    selectable = alias = spelling = latent = 0
    for family in inventory:
        for spec in family["public_values"].values():
            if spec.get("reachability") == "latent":
                latent += 1
            elif spec.get("alias_of"):
                alias += 1
            elif spec.get("spelling_of"):
                spelling += 1
            else:
                selectable += 1
    print(
        f"Capability audit passed: {len(inventory)} families; "
        f"{selectable} canonical values, {spelling} accepted spelling, "
        f"{alias} deprecated alias, {latent} latent; "
        f"full-chain parity verified."
    )
    print(f"  coverage enforced: {', '.join(enforced) or 'none'}")
    if pending:
        print(f"  coverage advisory (backfill pending): {', '.join(pending)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
