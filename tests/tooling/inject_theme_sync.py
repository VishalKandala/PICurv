#!/usr/bin/env python3
"""Inject the theme-sync script into generated Doxygen HTML so the built tree is publishable as-is."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


SCRIPT_NAME = "theme-sync.js"


def parse_args() -> argparse.Namespace:
    """!
    @brief Parse command-line arguments.
    @return Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Copy theme-sync.js into the generated HTML tree and reference it from every page.\n"
            "Idempotent: pages that already reference the script are left untouched."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root directory (default: parent of this script).",
    )
    parser.add_argument(
        "--html-dir",
        default="docs_build/html",
        help="Generated HTML directory relative to repo root (default: docs_build/html).",
    )
    return parser.parse_args()


def relative_script_path(page: Path, html_dir: Path) -> str:
    """!
    @brief Build the script reference for one page, correct for pages in subdirectories.
    @param[in] page Generated HTML page.
    @param[in] html_dir Root of the generated HTML tree.
    @return Relative href to the theme-sync script.
    """
    depth = len(page.relative_to(html_dir).parts) - 1
    return "../" * depth + SCRIPT_NAME


def inject(html_dir: Path, source: Path) -> tuple[int, int]:
    """!
    @brief Copy the script into the tree and reference it from every generated page.
    @param[in] html_dir Root of the generated HTML tree.
    @param[in] source Authoritative theme-sync script in the repository.
    @return Counts of pages injected and pages already carrying the reference.
    """
    (html_dir / SCRIPT_NAME).write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    injected = skipped = 0
    for page in sorted(html_dir.rglob("*.html")):
        markup = page.read_text(encoding="utf-8", errors="replace")
        # Match the actual script reference, not a prose mention: pages that document
        # theme-sync.js contain the bare filename and still need the tag injected.
        if re.search(rf'<script[^>]+src="[^"]*{re.escape(SCRIPT_NAME)}"', markup):
            skipped += 1
            continue
        if "</head>" not in markup:
            continue
        tag = f'<script type="text/javascript" src="{relative_script_path(page, html_dir)}"></script>\n</head>'
        page.write_text(markup.replace("</head>", tag, 1), encoding="utf-8")
        injected += 1
    return injected, skipped


def main() -> int:
    """!
    @brief Make the generated HTML tree self-contained with respect to theme syncing.
    @return Process status code.
    """
    args = parse_args()
    repo_root = args.repo_root.resolve()
    html_dir = (repo_root / args.html_dir).resolve()
    source = repo_root / "docs" / SCRIPT_NAME
    if not html_dir.is_dir():
        print(f"theme-sync injection failed: {args.html_dir} does not exist.", file=sys.stderr)
        return 1
    if not source.is_file():
        print(f"theme-sync injection failed: {source} is missing.", file=sys.stderr)
        return 1
    injected, skipped = inject(html_dir, source)
    print(f"[theme-sync] injected into {injected} page(s); {skipped} already referenced it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
