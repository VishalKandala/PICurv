#!/usr/bin/env python3
"""Report public selector surfaces that no capability family covers."""

from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
CENSUS = REPO_ROOT / "tests" / "tooling" / "family_census.json"


#: Modules whose public constants the census reads, paired with the dotted name a
#: capability family refers to. A package is scanned file by file but reported under
#: the package, because that is the surface an importer and a registry address.
SCANNED_MODULE_ROOTS = ("picurv_cli/core.py", "picurv_cli/storage", "picurv_cli/cli.py")


def scanned_modules():
    """!
    @brief Resolve the module roots into concrete files and their public dotted names.
    @return List of `(repo-relative path, dotted module name)` pairs.
    """
    resolved = []
    for root in SCANNED_MODULE_ROOTS:
        path = REPO_ROOT / root
        if path.is_dir():
            dotted = root.replace("/", ".")
            for child in sorted(path.glob("*.py")):
                if child.name != "__init__.py":
                    resolved.append((str(child.relative_to(REPO_ROOT)), dotted))
        else:
            resolved.append((root, root.replace("/", ".")[: -len(".py")]))
    return resolved


def declared_surfaces() -> set:
    """!
    @brief Public surfaces already covered by a registered capability family.
    @return Set of `module::symbol` identifiers.
    """
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    covered = set()
    for family in registry["families"]:
        surface = family.get("public_surface", {})
        if surface.get("module") and surface.get("symbol"):
            covered.add(f"{surface['module']}::{surface['symbol']}")
    return covered


# Module-level names with these suffixes declare a closed public choice set. This is
# the discovery contract: a choice set written as an inline literal at its point of use
# is invisible to the census, so the rule is that it must be named. Adding a suffix here
# widens discovery; removing an authoritative constant narrows it, which is why the
# committed classifications are asserted against this list by the test suite.
CHOICE_SET_SUFFIXES = (
    "MAP", "SPECS", "MODES", "TYPES", "TASKS", "CHOICES", "POLICIES",
    "STRUCTURES", "MODELS", "OUTPUTS", "FORMATS", "SPELLINGS", "EXTENSIONS",
    "SYMBOLS", "NAMES", "PROFILES", "METHODS", "OPERATORS", "KINDS", "SPELLINGS",
)
_CHOICE_SET_PATTERN = re.compile(r"[A-Z0-9_]+(" + "|".join(CHOICE_SET_SUFFIXES) + r")$")


def candidate_surfaces() -> dict:
    """!
    @brief Public choice points in the CLI, found independently of the registry.

    @details A `normalize_*` function, or a module-level constant whose name ends in one
             of CHOICE_SET_SUFFIXES, is how this codebase exposes a user-selectable set.
             Finding them independently of the registry is the point: a registry that
             only checks itself cannot report a family nobody registered.
    @return Mapping of `module::symbol` to a short description.
    """
    found: dict = {}
    for module, dotted in scanned_modules():
        path = REPO_ROOT / module
        if not path.is_file():
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and re.fullmatch(r"normalize_\w+", node.name):
                found[f"{dotted}::{node.name}"] = "normalizer"
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and _CHOICE_SET_PATTERN.fullmatch(target.id):
                        found[f"{dotted}::{target.id}"] = "choice set"
    return found


# Every discovered surface must land in exactly one of these. The vocabulary is the
# result of the explicit census: a closed choice set is either documented as a
# capability family, or it is one of these five things, and saying which is a
# deliberate act rather than an omission.
VALID_CLASSIFICATIONS = (
    "public_pending",            # A public family that owes Tier-1/Tier-2 coverage.
    "parameter_of_entry",        # A knob an existing capability entry already owns.
    "spelling_alias",            # Accepted spellings resolving to a canonical value.
    "structural",                # Free-form or externally-owned; no closed PICurv set.
    "cli_inventory_sufficient",  # The generated CLI reference answers it fully.
    "not_public",                # Internal; not reachable as a user choice.
)

#: Classifications that must name what owns them, so the pointer is checkable.
NEEDS_OWNER = ("parameter_of_entry", "spelling_alias")


def classification_of(entry) -> tuple:
    """!
    @brief Normalize a census entry into a typed classification and reason.

    @details Older entries were free-text, which let a pending public family be
             summarized as "not a public family". A typed record cannot do that.
    @param[in] entry Census entry, either a mapping or a legacy string.
    @return Tuple of (classification, reason).
    """
    if isinstance(entry, dict):
        return entry.get("classification", ""), entry.get("reason", "")
    text = str(entry)
    if "not a public capability family" in text:
        return "not_public", text
    return "public_pending", text


def main() -> int:
    """!
    @brief Report uncovered public selector surfaces, typed.

    @details Fails when a surface is undiscovered, when a classification is invalid,
             when a `not_public` entry gives no reason, or when any public family
             remains pending. Publication cannot proceed while a public family is
             undocumented.
    @return Process status code.
    """
    census = json.loads(CENSUS.read_text(encoding="utf-8")) if CENSUS.is_file() else {}
    acknowledged = census.get("acknowledged_uncovered", {})
    covered = declared_surfaces()
    candidates = candidate_surfaces()

    uncovered = {k: v for k, v in candidates.items() if k not in covered}
    unacknowledged = sorted(k for k in uncovered if k not in acknowledged)
    stale = sorted(k for k in acknowledged if k not in candidates)

    pending: list = []
    not_public: list = []
    malformed: list = []
    for key, entry in sorted(acknowledged.items()):
        classification, reason = classification_of(entry)
        if classification not in VALID_CLASSIFICATIONS:
            malformed.append(f"{key}: classification '{classification}' is not one of {list(VALID_CLASSIFICATIONS)}")
            continue
        if classification != "public_pending" and not reason.strip():
            malformed.append(
                f"{key}: classified {classification} but gives no reason. Every "
                f"classification other than public_pending is a claim, and a claim owes "
                f"its reasoning"
            )
            continue
        owner = entry.get("owner") if isinstance(entry, dict) else None
        if classification in NEEDS_OWNER and not owner:
            malformed.append(
                f"{key}: classified {classification} but names no owner. Say which "
                f"capability entry or canonical value it belongs to"
            )
            continue
        (pending if classification == "public_pending" else not_public).append(key)

    problems: list = []
    if stale:
        problems += [f"census acknowledges '{k}', which no longer exists" for k in stale]
    if unacknowledged:
        problems += [f"'{k}' ({candidates[k]}) has no capability family and no census entry" for k in unacknowledged]
    problems += malformed

    if problems:
        print("Capability family census violations:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print(
            "\nRegister a family in capability_families.json, or classify the surface in\n"
            "family_census.json as one of "
            + ", ".join(VALID_CLASSIFICATIONS) + ",\nwith a reason. An unregistered "
            "family is invisible to every other gate.",
            file=sys.stderr,
        )
        return 1

    from collections import Counter
    spread = Counter(classification_of(entry)[0] for entry in acknowledged.values())
    breakdown = ", ".join(f"{count} {name}" for name, count in sorted(spread.items()))
    print(
        f"Family census: {len(covered)} covered by a family; {breakdown} "
        f"({len(candidates)} surfaces examined)."
    )
    if pending:
        print("\nPublic capability families still awaiting Tier-2 backfill:", file=sys.stderr)
        for key in pending:
            print(f"  {key}", file=sys.stderr)
        print(
            "\nThese are public choice points with no documented capability entries.\n"
            "Publication cannot proceed while any public family is pending.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
