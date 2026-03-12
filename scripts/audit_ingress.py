#!/usr/bin/env python3
"""
Static ingress audit for PETSc option parsing in setup/io.

This script scans src/setup.c and src/io.c for PetscOptionsGet*/HasName calls,
extracts option flags, and compares them against a maintained manifest.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from typing import Iterable, Set


OPTION_RE = re.compile(
    r'PetscOptions(?:Get(?:Int|Real|Bool|String|IntArray|RealArray)|HasName)\s*'
    r'\(\s*NULL\s*,\s*NULL\s*,\s*"(-[^"]+)"'
)


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


def load_manifest(path: pathlib.Path) -> Set[str]:
    """!
    @brief Load manifest.
    @param[in] path Filesystem path argument passed to `load_manifest()`.
    @return Value returned by `load_manifest()`.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    options = data.get("known_petsc_options")
    if not isinstance(options, list):
        raise ValueError("Manifest key 'known_petsc_options' must be a list.")
    bad = [opt for opt in options if not isinstance(opt, str) or not opt.startswith("-")]
    if bad:
        raise ValueError(f"Manifest has invalid option entries: {bad}")
    return set(options)


def main() -> int:
    """!
    @brief Entry point for this script.
    @return Value returned by `main()`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Scan PETSc option ingress in src/setup.c and src/io.c, then compare "
            "against scripts/audit_ingress_manifest.json."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/audit_ingress.py\n"
            "  python3 scripts/audit_ingress.py --show-scanned\n"
            "  python3 scripts/audit_ingress.py --manifest scripts/audit_ingress_manifest.json\n"
        ),
    )
    parser.add_argument(
        "--manifest",
        default="scripts/audit_ingress_manifest.json",
        help=(
            "Manifest JSON path, relative to repository root unless absolute "
            "(default: scripts/audit_ingress_manifest.json)."
        ),
    )
    parser.add_argument(
        "--show-scanned",
        action="store_true",
        help="Print discovered PETSc options before drift comparison.",
    )
    args = parser.parse_args()

    repo_root = pathlib.Path(__file__).resolve().parents[1]
    manifest_path = (repo_root / args.manifest).resolve()
    scan_paths = [repo_root / "src/setup.c", repo_root / "src/io.c"]

    if not manifest_path.exists():
        print(f"[ERROR] Manifest not found: {manifest_path}", file=sys.stderr)
        return 2

    scanned = scan_petsc_options(scan_paths)
    expected = load_manifest(manifest_path)

    missing_in_manifest = sorted(scanned - expected)
    stale_in_manifest = sorted(expected - scanned)

    if args.show_scanned:
        print("[INFO] Scanned PETSc options:")
        for flag in sorted(scanned):
            print(flag)
        print("")

    print(f"[INFO] Scanned options: {len(scanned)}")
    print(f"[INFO] Manifest options: {len(expected)}")

    if missing_in_manifest:
        print("[ERROR] New PETSc ingress options missing in manifest:")
        for flag in missing_in_manifest:
            print(f"  - {flag}")

    if stale_in_manifest:
        print("[ERROR] Manifest options no longer present in setup/io scan:")
        for flag in stale_in_manifest:
            print(f"  - {flag}")

    if missing_in_manifest or stale_in_manifest:
        print(
            "[FAIL] Ingress drift detected. Update scripts/audit_ingress_manifest.json and docs mapping.",
            file=sys.stderr,
        )
        return 1

    print("[OK] Ingress manifest matches setup/io PETSc option scan.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
