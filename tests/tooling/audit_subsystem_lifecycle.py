#!/usr/bin/env python3
"""Enforce the documentation obligations each subsystem lifecycle status carries."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
RECORDS = REPO_ROOT / "tests" / "tooling" / "subsystem_records.json"
PAGE_TYPES = REPO_ROOT / "tests" / "tooling" / "page_types.json"
FAMILIES = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
PAGE_DIRS = ("docs/pages", "docs")

# The ladder is cumulative: claiming a rung owes everything below it too. This is
# the whole point of the gate - experimental work is never blocked by documentation
# it cannot honestly write yet, but it also cannot skip a rung by claiming the top.
LADDER = ("planned", "internal", "experimental", "supported")

LADDER_OBLIGATIONS = {
    "planned": ("purpose", "intended_scope", "design_owner", "not_implemented_status"),
    "internal": ("scope_boundary", "architecture_boundary", "dependencies",
                 "developer_entry_points"),
    "experimental": ("configuration", "selection_guidance", "observability", "limitations",
                     "safe_use_boundaries"),
    "supported": ("evidence", "lifecycle_and_restart", "operations", "troubleshooting",
                  "examples", "complete_reference"),
}

# Statuses off the ladder. Each owes its own set; those that describe something that
# was once built also owe the ladder up to the peak status it reached.
TERMINAL_OBLIGATIONS = {
    "known-defective": ("defect_disclosure", "defect_scope_record", "safe_use_boundaries",
                        "limitations"),
    "deprecated": ("migration_path", "replacement", "compatibility_period", "removal_policy"),
    "removed": ("history_record", "rejection_behavior"),
}
# Obligations that make sense only while the subsystem is unbuilt. A supported
# subsystem does not owe a not-implemented notice.
NON_INHERITED = ("not_implemented_status",)
NEEDS_PEAK = ("known-defective", "deprecated", "removed")
# `removed` documents an absence, so it does not owe the prose of what it replaced.
INHERITS_LADDER = ("known-defective", "deprecated")

VALID_STATUSES = tuple(LADDER) + tuple(TERMINAL_OBLIGATIONS)

# Which status may follow which. Absent from a row means the transition is invalid:
# a subsystem cannot re-enter the ladder below where it stood, and cannot reach
# `removed` without first being deprecated or declared defective.
#
# `planned -> removed` is the one deliberate exception, and it is not a removal in the
# same sense: nothing was ever built, so a cancelled design owes a history record and
# a rejection behaviour rather than a migration path. Every other route to `removed`
# passes through `deprecated` or `known-defective`.
TRANSITIONS = {
    None: set(VALID_STATUSES) - {"removed"},
    "planned": {"planned", "internal", "experimental", "supported", "removed"},
    "internal": {"internal", "experimental", "supported", "known-defective", "deprecated"},
    "experimental": {"experimental", "supported", "known-defective", "deprecated"},
    "supported": {"supported", "known-defective", "deprecated"},
    "known-defective": {"known-defective", "experimental", "supported", "deprecated", "removed"},
    "deprecated": {"deprecated", "removed"},
    "removed": {"removed"},
}

VALID_VISIBILITY = ("internal", "public")

# The concern vocabulary of 64_Documentation_Extension_Framework section 5.
VALID_CONCERNS = (
    "numerical_method", "user_selector", "units_and_nondimensionalization",
    "artifact_topology", "persistent_restart_state", "determinism_and_reproducibility",
    "mpi_distributed_execution", "external_service", "security_credentials",
    "generated_artifact", "file_format", "destructive_scope", "backward_compatibility",
    "scientific_verification",
)

RECORD_KEYS = {"id", "title", "status", "visibility", "previous_status", "peak_status",
               "proposed_status", "promotion_rationale", "capability_families",
               "obligations", "concerns", "note"}
# A reason short enough to fit here is an evasion, not a reason.
MIN_REASON_CHARS = 30
FILLER_REASONS = {"n/a", "na", "none", "not applicable", "no", "-", "tbd"}


def page_index() -> dict:
    """!
    @brief Every documentation page id with the anchors it defines.
    @return Mapping of page id to the set of `@section`, `@subsection`, and `@anchor` names.
    """
    index: dict = {}
    for directory in PAGE_DIRS:
        for markdown in sorted((REPO_ROOT / directory).glob("*.md")):
            text = markdown.read_text(encoding="utf-8")
            match = re.search(r"^@page\s+(\S+)", text, re.M)
            if not match or match.group(1).startswith("<"):
                continue
            anchors = set(re.findall(r"^@(?:section|subsection|subsubsection|anchor)\s+(\S+)",
                                     text, re.M))
            index[match.group(1)] = anchors
    return index


def required_obligations(record: dict) -> tuple:
    """!
    @brief The obligation ids a record owes.

    @details A record may declare a `proposed_status` above its claimed one. That is
             how a subsystem documented to a higher bar waits for the owner to decide
             whether it has actually earned it: obligations are checked against the
             proposal, so the writing is held to the higher standard, while the status
             the documentation publishes stays conservative until a human agrees.
    @param[in] record One subsystem record.
    @return Ordered tuple of obligation ids.
    """
    status = record.get("proposed_status") or record.get("status")
    owed: list = []
    if status in LADDER:
        for rung in LADDER[: LADDER.index(status) + 1]:
            owed.extend(LADDER_OBLIGATIONS[rung])
    else:
        if status in INHERITS_LADDER:
            peak = record.get("peak_status")
            if peak in LADDER:
                for rung in LADDER[: LADDER.index(peak) + 1]:
                    owed.extend(LADDER_OBLIGATIONS[rung])
        owed.extend(TERMINAL_OBLIGATIONS.get(status, ()))
    if status != "planned":
        owed = [key for key in owed if key not in NON_INHERITED]
    # A status may re-owe something a lower rung already required; keep it once.
    return tuple(dict.fromkeys(owed))


def check_satisfaction(context: str, key: str, spec, pages: dict, published: set,
                       public: bool) -> list:
    """!
    @brief Verify one obligation or concern is genuinely answered.

    @details An answer is a documentation reference that resolves, or a stated reason
             for non-applicability, or a literal value where the obligation is a fact
             rather than prose. An empty or filler answer is a violation.
    @param[in] context Record id, for messages.
    @param[in] key Obligation or concern id.
    @param[in] spec The declared answer.
    @param[in] pages Page index from page_index().
    @param[in] published Ids of the pages the site publishes.
    @param[in] public Whether the subsystem is publicly visible.
    @return List of violation strings.
    """
    if not isinstance(spec, dict):
        return [f"{context}: '{key}' must be an object, got {type(spec).__name__}"]
    problems = []
    if "not_applicable" in spec:
        reason = str(spec["not_applicable"]).strip()
        if reason.lower().rstrip(".") in FILLER_REASONS or len(reason) < MIN_REASON_CHARS:
            problems.append(
                f"{context}: '{key}' is declared not applicable without a stated reason "
                f"({reason!r}). A reasoned N/A is accepted; a bare one is not"
            )
        return problems
    if "value" in spec:
        if not str(spec["value"]).strip():
            problems.append(f"{context}: '{key}' declares an empty value")
        return problems
    page = spec.get("page")
    if not page:
        problems.append(
            f"{context}: '{key}' is unsatisfied - give it a page reference, a value, or a "
            f"reasoned not_applicable"
        )
        return problems
    if page not in pages:
        problems.append(f"{context}: '{key}' cites page '{page}', which does not exist")
        return problems
    anchor = spec.get("anchor")
    if anchor and anchor not in pages[page]:
        problems.append(
            f"{context}: '{key}' cites anchor '{anchor}' on page '{page}', which does not "
            f"define it"
        )
    if public and page not in published:
        problems.append(
            f"{context}: '{key}' cites '{page}', which is not a published page, but the "
            f"subsystem is publicly visible"
        )
    return problems


def validate(records: list, pages: dict, families: dict, published: set) -> list:
    """!
    @brief Validate every subsystem record against the lifecycle contract.
    @param[in] records Subsystem records.
    @param[in] pages Page index from page_index().
    @param[in] families Capability family metadata, keyed by family id.
    @param[in] published Ids of the pages the site publishes.
    @return List of violation strings; empty means the contract holds.
    """
    problems: list = []
    seen: set = set()
    for record in records:
        identifier = record.get("id")
        if not identifier:
            problems.append("a record declares no id")
            continue
        if identifier in seen:
            problems.append(f"{identifier}: declared more than once")
        seen.add(identifier)

        unknown = set(record) - RECORD_KEYS
        if unknown:
            problems.append(f"{identifier}: unknown field(s) {sorted(unknown)}")

        status = record.get("status")
        if status not in VALID_STATUSES:
            problems.append(
                f"{identifier}: status {status!r} is not one of {list(VALID_STATUSES)}"
            )
            continue

        visibility = record.get("visibility")
        if visibility not in VALID_VISIBILITY:
            problems.append(
                f"{identifier}: visibility {visibility!r} is not one of "
                f"{list(VALID_VISIBILITY)}"
            )
            continue
        if status == "internal" and visibility != "internal":
            problems.append(
                f"{identifier}: status 'internal' contradicts visibility 'public'; a "
                f"subsystem the user can reach is at least experimental"
            )

        previous = record.get("previous_status", None)
        if previous is not None and previous not in VALID_STATUSES:
            problems.append(f"{identifier}: previous_status {previous!r} is not a valid status")
        elif status not in TRANSITIONS[previous]:
            problems.append(
                f"{identifier}: {previous or 'a new record'} -> '{status}' is not a valid "
                f"lifecycle transition (allowed: {sorted(TRANSITIONS[previous])})"
            )

        proposed = record.get("proposed_status")
        if proposed is not None:
            if proposed not in LADDER:
                problems.append(
                    f"{identifier}: proposed_status {proposed!r} must name a rung of "
                    f"{list(LADDER)}"
                )
            elif status not in LADDER or LADDER.index(proposed) <= LADDER.index(status):
                problems.append(
                    f"{identifier}: proposed_status '{proposed}' is not above the claimed "
                    f"status '{status}'; a proposal that is not a promotion is noise"
                )
            elif not record.get("promotion_rationale"):
                problems.append(
                    f"{identifier}: proposes '{proposed}' but gives no promotion_rationale. "
                    f"Say what the owner is being asked to confirm"
                )

        if status in NEEDS_PEAK:
            peak = record.get("peak_status")
            if peak not in LADDER:
                problems.append(
                    f"{identifier}: status '{status}' requires peak_status naming the highest "
                    f"rung it reached, one of {list(LADDER)}"
                )
        elif record.get("peak_status") is not None:
            problems.append(
                f"{identifier}: peak_status applies only to {list(NEEDS_PEAK)}, not '{status}'"
            )

        obligations = record.get("obligations") or {}
        owed = required_obligations(record)
        for key in owed:
            if key not in obligations:
                problems.append(
                    f"{identifier}: status '{status}' owes '{key}', which is not declared"
                )
            else:
                problems.extend(
                    check_satisfaction(identifier, key, obligations[key], pages, published,
                                       visibility == "public")
                )
        for key in sorted(set(obligations) - set(owed)):
            problems.append(
                f"{identifier}: declares obligation '{key}', which status '{status}' does not "
                f"owe; remove it or claim the status that owes it"
            )

        concerns = record.get("concerns") or {}
        for key in sorted(concerns):
            if key not in VALID_CONCERNS:
                problems.append(
                    f"{identifier}: concern '{key}' is not in the concern vocabulary of "
                    f"64_Documentation_Extension_Framework"
                )
                continue
            problems.extend(
                check_satisfaction(identifier, key, concerns[key], pages, published,
                                   visibility == "public")
            )

        # A planned subsystem must not already be presented as working behavior.
        for family_id in record.get("capability_families", []):
            family = families.get(family_id)
            if family is None:
                problems.append(
                    f"{identifier}: cites capability family '{family_id}', which is not "
                    f"registered"
                )
                continue
            if status != "planned":
                continue
            live = sorted(
                name for name, meta in family.get("value_metadata", {}).items()
                if (meta or {}).get("status") == "supported"
            )
            if live:
                problems.append(
                    f"{identifier}: is 'planned' but family '{family_id}' already offers "
                    f"supported value(s) {live}; a planned subsystem must not appear as "
                    f"supported behavior"
                )
    return problems


def main() -> int:
    """!
    @brief Report subsystem lifecycle violations.
    @return Process status code.
    """
    document = json.loads(RECORDS.read_text(encoding="utf-8"))
    records = document["subsystems"]
    published = set(json.loads(PAGE_TYPES.read_text(encoding="utf-8"))["assignments"])
    families = {
        family["id"]: family
        for family in json.loads(FAMILIES.read_text(encoding="utf-8"))["families"]
    }
    problems = validate(records, page_index(), families, published)
    if problems:
        print("Subsystem lifecycle violations:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print(
            "\nObligations grow with the status a subsystem claims. See\n"
            "64_Documentation_Extension_Framework section 4. Lower the claimed status, or\n"
            "write the documentation that status owes - do not add empty prose to pass.",
            file=sys.stderr,
        )
        return 1

    from collections import Counter
    spread = Counter(record["status"] for record in records)
    summary = ", ".join(f"{count} {status}" for status, count in sorted(spread.items()))
    print(f"Subsystem lifecycle audit passed: {len(records)} subsystem(s) ({summary}).")
    proposals = [r for r in records if r.get("proposed_status")]
    if proposals:
        print(f"\n{len(proposals)} promotion(s) awaiting the owner's decision. The "
              f"documentation meets the higher bar; whether the subsystem has earned "
              f"the status is a human judgement this gate cannot make:")
        for record in proposals:
            print(f"  {record['id']}: {record['status']} -> {record['proposed_status']}")
            print(f"      {record['promotion_rationale']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
