#!/usr/bin/env python3
"""Backward-compatible wrapper for legacy grid conversion.

This script now delegates to:
    python3 scripts/grid.gen legacy1d ...
so there is a single maintained conversion implementation.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    """Parse wrapper CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Convert a legacy 1D-axis grid payload into canonical PICGRID format.",
    )
    parser.add_argument("--input", required=True, help="Path to legacy grid input file.")
    parser.add_argument("--output", required=True, help="Path to converted PICGRID output file.")
    parser.add_argument(
        "--axis-columns",
        type=int,
        nargs=3,
        default=[0, 1, 2],
        metavar=("XCOL", "YCOL", "ZCOL"),
        help="Preferred source column index for X/Y/Z axis rows (default: 0 1 2).",
    )
    parser.add_argument(
        "--allow-trailing",
        action="store_true",
        help="Ignore trailing payload rows after expected legacy content.",
    )
    return parser.parse_args()


def main() -> int:
    """Entry point."""
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gridgen = os.path.join(script_dir, "grid.gen")
    if not os.path.isfile(gridgen):
        print(f"Error: grid.gen not found at '{gridgen}'.", file=sys.stderr)
        return 1

    cmd = [
        sys.executable,
        gridgen,
        "legacy1d",
        "--input",
        args.input,
        "--output",
        args.output,
        "--axis-columns",
        str(args.axis_columns[0]),
        str(args.axis_columns[1]),
        str(args.axis_columns[2]),
        "--no-write-vtk",
    ]
    cmd.append("--allow-trailing" if args.allow_trailing else "--strict-trailing")

    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip(), file=sys.stderr)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
