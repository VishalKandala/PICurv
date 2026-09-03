#!/usr/bin/env python3
"""Reject public closed choices written as inline literals.

The family census discovers `normalize_*` functions and named module-level constants.
That is a complete search of what it can see - and an inline `if value not in {"a","b"}`
is invisible to it. The census could therefore call itself exhaustive while public
choices sat in literals nobody enumerated.

This closes the loop from the other side: every inline set of string choices in the CLI
package must either become a named constant, which the census then discovers and
classifies, or be listed here with the same typed classification and a reason. A
literal that is neither is a violation.
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
WAIVERS = REPO_ROOT / "tests" / "tooling" / "inline_choice_waivers.json"
#: A package is scanned file by file but reported under the package name, because that
#: is the surface a waiver and a registry address.
MODULE_ROOTS = ("picurv_cli/core.py", "picurv_cli/cli.py", "picurv_cli/storage",
                "picurv_cli/main.py")


def _resolve_modules():
    """!
    @brief Expand the module roots into `(path, dotted name)` pairs.
    @return List of scanned module pairs.
    """
    resolved = []
    for root in MODULE_ROOTS:
        path = REPO_ROOT / root
        if path.is_dir():
            dotted = root.replace("/", ".")
            for child in sorted(path.glob("*.py")):
                if child.name != "__init__.py":
                    resolved.append((str(child.relative_to(REPO_ROOT)), dotted))
        else:
            resolved.append((root, root.replace("/", ".")[: -len(".py")]))
    return resolved


MODULES = tuple(item[0] for item in _resolve_modules())
MODULE_DOTTED = dict(_resolve_modules())

# The census vocabulary, minus `public_pending`: an inline literal cannot be a pending
# family, because a family needs a named surface to point at.
VALID_CLASSIFICATIONS = (
    "parameter_of_entry", "spelling_alias", "structural",
    "cli_inventory_sufficient", "not_public",
)
MIN_REASON_CHARS = 40
#: Two-element sets are the smallest thing worth calling a choice.
MINIMUM_CHOICE_SIZE = 2


def literal_key(module: str, values: list) -> str:
    """!
    @brief Stable identity for one inline choice set.

    @details Keyed by module and sorted values rather than by line number, so the
             waiver survives edits above it and two sites sharing a set collapse to one
             entry - which is itself a signal that the set wants a name.
    @param[in] module Dotted module name.
    @param[in] values The string members.
    @return Waiver key.
    """
    return f"{module}::{{{','.join(sorted(values))}}}"


def find_inline_choices() -> dict:
    """!
    @brief Every inline set of string choices in the CLI package.

    @details Two shapes are recognised: a membership test against a literal collection,
             and an argparse `choices=` literal. Both are how a closed public set gets
             written without a name.
    @return Mapping of waiver key to the list of sites that produced it.
    """
    found: dict = {}
    for module in MODULES:
        path = REPO_ROOT / module
        if not path.is_file():
            continue
        dotted = MODULE_DOTTED.get(module, module.replace("/", ".")[: -len(".py")])
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            values = None
            shape = None
            if isinstance(node, ast.Compare) and node.ops and \
                    isinstance(node.ops[0], (ast.In, ast.NotIn)):
                collection = node.comparators[0]
                if isinstance(collection, (ast.Set, ast.List, ast.Tuple)):
                    values, shape = string_members(collection), "membership test"
            elif isinstance(node, ast.Call):
                for keyword in node.keywords:
                    if keyword.arg == "choices" and \
                            isinstance(keyword.value, (ast.Set, ast.List, ast.Tuple)):
                        values, shape = string_members(keyword.value), "argparse choices"
            if values and len(values) >= MINIMUM_CHOICE_SIZE:
                key = literal_key(dotted, values)
                found.setdefault(key, []).append(f"{module}:{node.lineno} ({shape})")
    return found


def string_members(collection) -> list:
    """!
    @brief The string members of a literal collection, or None if it is not all strings.
    @param[in] collection An ast Set, List, or Tuple node.
    @return List of string values, or None.
    """
    members = [element.value for element in collection.elts
               if isinstance(element, ast.Constant) and isinstance(element.value, str)]
    return members if len(members) == len(collection.elts) and members else None


def main() -> int:
    """!
    @brief Fail when an inline choice set is neither named nor classified.
    @return Process status code.
    """
    document = json.loads(WAIVERS.read_text(encoding="utf-8"))
    waivers = document["inline_choices"]
    found = find_inline_choices()

    problems: list = []
    for key in sorted(found):
        entry = waivers.get(key)
        if entry is None:
            problems.append(
                f"{key}\n      at {', '.join(found[key])}\n"
                f"      Give this set a name so the family census can see it, or "
                f"classify it in inline_choice_waivers.json."
            )
            continue
        classification = entry.get("classification")
        if classification not in VALID_CLASSIFICATIONS:
            problems.append(f"{key}: classification {classification!r} is not one of "
                            f"{list(VALID_CLASSIFICATIONS)}")
        elif len(str(entry.get("reason", "")).strip()) < MIN_REASON_CHARS:
            problems.append(f"{key}: classified {classification} without a stated reason")
    for stale in sorted(set(waivers) - set(found)):
        problems.append(f"{stale}: waived, but no such inline literal exists any more")

    if problems:
        print("Inline choice-set violations:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1

    from collections import Counter
    spread = Counter(waivers[key]["classification"] for key in found)
    summary = ", ".join(f"{count} {name}" for name, count in sorted(spread.items()))
    print(f"Inline choice audit passed: {len(found)} inline set(s), every one classified "
          f"({summary}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
