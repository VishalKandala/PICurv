#!/usr/bin/env python3
"""
Static ingress audit for PETSc option parsing in setup/io.

This script scans the source files named by the manifest for PetscOptionsGet*/HasName
calls, extracts option flags, and compares them against a maintained manifest.

Two kinds of option name reach C. Most are string literals. A variable-arity
configuration, such as the field-statistics window list, must construct its names
instead, and a constructed name is invisible to a literal scan. Rather than let that
become a hole, this audit also collects the format strings those names are built
from, requires each to be declared as a family, and fails on any constructed name it
cannot trace back to one. The invariant is that every option name reaching C is
either a declared literal or a declared family.

Reads that pass a private PetscOptions object rather than NULL are deliberately out
of scope: they read checkpoint state, not user configuration.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from typing import Iterable, Set


_ACCESSOR = r'PetscOptions(?:Get(?:Int|Real|Bool|String|IntArray|RealArray)|HasName)'

OPTION_RE = re.compile(_ACCESSOR + r'\s*\(\s*NULL\s*,\s*NULL\s*,\s*"(-[^"]+)"')

#: A read whose name comes from a variable rather than a literal.
CONSTRUCTED_RE = re.compile(
    _ACCESSOR + r'\s*\(\s*NULL\s*,\s*NULL\s*,\s*([A-Za-z_][A-Za-z0-9_]*)\s*,'
)

#: The start of a name-building call; the format argument is parsed separately
#: because it may be several literals spliced around a width macro.
SNPRINTF_RE = re.compile(r'PetscSNPrintf\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*,')

LITERAL_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')


def canonical_family(fragments: Iterable[str]) -> str:
    """!
    @brief Join the literal fragments of a constructed name into one stable pattern.
    @details Each substituted value collapses to a bare `%`, so a name spliced around
             `PetscInt_FMT` and one written with `%d` produce the same family.
    @param[in] fragments Literal chunks of the format expression, in source order.
    @return Value returned by `canonical_family()`.
    """
    joined = "".join(fragments)
    return re.sub(r"%[-+ #0-9.]*[a-zA-Z]?", "%", joined)


def extract_families(text: str) -> "dict[str, Set[str]]":
    """!
    @brief Map each name variable built by PetscSNPrintf to every family it carries.
    @details One buffer is normally reused for many names, so a variable maps to a
             set rather than a single pattern; collapsing it to the last assignment
             would hide every family but one.
    @param[in] text Source text of one translation unit.
    @return Value returned by `extract_families()`.
    """
    families: "dict[str, Set[str]]" = {}
    for match in SNPRINTF_RE.finditer(text):
        variable = match.group(1)
        depth = 1
        index = match.end()
        while index < len(text) and depth:
            if text[index] == "(":
                depth += 1
            elif text[index] == ")":
                depth -= 1
            index += 1
        fragments = LITERAL_RE.findall(text[match.end():index])
        if not fragments or not fragments[0].startswith("-"):
            continue
        families.setdefault(variable, set()).add(canonical_family(fragments))
    return families


def scan_option_families(paths: Iterable[pathlib.Path]) -> "tuple[Set[str], list[str]]":
    """!
    @brief Collect declared option families and any constructed name that lacks one.
    @param[in] paths Source files to scan.
    @return Value returned by `scan_option_families()`.
    """
    families: Set[str] = set()
    untraceable: "list[str]" = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        built = extract_families(text)
        for match in CONSTRUCTED_RE.finditer(text):
            variable = match.group(1)
            if variable in built:
                families.update(built[variable])
            else:
                untraceable.append(f"{path.name}: '{variable}'")
    return families, untraceable


def scan_petsc_options(paths: Iterable[pathlib.Path]) -> Set[str]:
    """!
    @brief Perform scan petsc options.
    @param[in] paths Argument passed to `scan_petsc_options()`.
    @return Value returned by `scan_petsc_options()`.
    """
    flags: Set[str] = set()
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for match in OPTION_RE.finditer(text):
            flags.add(match.group(1))
    return flags


def load_manifest(path: pathlib.Path) -> dict:
    """!
    @brief Read and validate the option names, families, and sources in the manifest.
    @param[in] path Filesystem path argument passed to `load_manifest()`.
    @return Value returned by `load_manifest()`.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    result = {}
    for key in ("known_petsc_options", "known_petsc_option_families"):
        entries = data.get(key, [])
        if not isinstance(entries, list):
            raise ValueError(f"Manifest key '{key}' must be a list.")
        bad = [opt for opt in entries if not isinstance(opt, str) or not opt.startswith("-")]
        if bad:
            raise ValueError(f"Manifest key '{key}' has invalid entries: {bad}")
        result[key] = set(entries)
    sources = data.get("sources")
    if not isinstance(sources, list) or not all(isinstance(src, str) for src in sources):
        raise ValueError("Manifest key 'sources' must be a list of paths.")
    result["sources"] = sources
    return result


def main() -> int:
    """!
    @brief Entry point for this script.
    @return Value returned by `main()`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Scan PETSc option ingress in src/setup.c and src/io.c, then compare "
            "against tests/tooling/audit_ingress_manifest.json."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 tests/tooling/audit_ingress.py\n"
            "  python3 tests/tooling/audit_ingress.py --show-scanned\n"
            "  python3 tests/tooling/audit_ingress.py --manifest tests/tooling/audit_ingress_manifest.json\n"
        ),
    )
    parser.add_argument(
        "--manifest",
        default="tests/tooling/audit_ingress_manifest.json",
        help=(
            "Manifest JSON path, relative to repository root unless absolute "
            "(default: tests/tooling/audit_ingress_manifest.json)."
        ),
    )
    parser.add_argument(
        "--show-scanned",
        action="store_true",
        help="Print discovered PETSc options before drift comparison.",
    )
    args = parser.parse_args()

    repo_root = pathlib.Path(__file__).resolve().parents[2]
    manifest_path = (repo_root / args.manifest).resolve()

    if not manifest_path.exists():
        print(f"[ERROR] Manifest not found: {manifest_path}", file=sys.stderr)
        return 2

    manifest = load_manifest(manifest_path)
    #: Scan paths come from the manifest so a new parse site cannot be added
    #: without also being declared to the audit.
    scan_paths = [repo_root / src for src in manifest["sources"]]
    absent = [str(path) for path in scan_paths if not path.exists()]
    if absent:
        print(f"[ERROR] Manifest names sources that do not exist: {absent}", file=sys.stderr)
        return 2

    scanned = scan_petsc_options(scan_paths)
    expected = manifest["known_petsc_options"]
    scanned_families, untraceable = scan_option_families(scan_paths)
    expected_families = manifest["known_petsc_option_families"]

    missing_in_manifest = sorted(scanned - expected)
    stale_in_manifest = sorted(expected - scanned)
    missing_families = sorted(scanned_families - expected_families)
    stale_families = sorted(expected_families - scanned_families)

    if args.show_scanned:
        print("[INFO] Scanned PETSc options:")
        for flag in sorted(scanned):
            print(flag)
        print("")

    print(f"[INFO] Scanned options: {len(scanned)}")
    print(f"[INFO] Manifest options: {len(expected)}")
    print(f"[INFO] Scanned option families: {len(scanned_families)}")
    print(f"[INFO] Manifest option families: {len(expected_families)}")

    if missing_in_manifest:
        print("[ERROR] New PETSc ingress options missing in manifest:")
        for flag in missing_in_manifest:
            print(f"  - {flag}")

    if stale_in_manifest:
        print("[ERROR] Manifest options no longer present in setup/io scan:")
        for flag in stale_in_manifest:
            print(f"  - {flag}")

    if missing_families:
        print("[ERROR] Constructed option families missing in manifest:")
        for family in missing_families:
            print(f"  - {family}")

    if stale_families:
        print("[ERROR] Manifest option families no longer present in the scan:")
        for family in stale_families:
            print(f"  - {family}")

    if untraceable:
        print("[ERROR] Option names built from a variable this audit cannot trace to a family:")
        for entry in untraceable:
            print(f"  - {entry}")
        print(
            "        Build the name with PetscSNPrintf from a literal format beginning with '-', "
            "so the audit can see it."
        )

    if missing_in_manifest or stale_in_manifest or missing_families or stale_families or untraceable:
        print(
            "[FAIL] Ingress drift detected. Update tests/tooling/audit_ingress_manifest.json and docs mapping.",
            file=sys.stderr,
        )
        return 1

    print("[OK] Ingress manifest matches the PETSc option scan of its declared sources.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
