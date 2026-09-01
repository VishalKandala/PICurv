#!/usr/bin/env python3
"""Local markdown link checker across every Markdown file the current commit carries."""

from __future__ import annotations

import argparse
import subprocess
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from repo_files import enumerate_repository_files  # noqa: E402


LINK_PATTERN = re.compile(r"!\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)|\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


def parse_args() -> argparse.Namespace:
    """!
    @brief Parse command-line arguments.
    @return Value returned by `parse_args()`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Check local markdown links across every repository Markdown file.\n"
            "HTTP(S), mailto, and in-page anchor links are skipped."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 tests/tooling/check_markdown_links.py\n"
            "  python3 tests/tooling/check_markdown_links.py --repo-root .\n"
            "  python3 tests/tooling/check_markdown_links.py --no-readme"
        ),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root directory (default: parent of this script).",
    )
    parser.add_argument(
        "--no-readme",
        action="store_true",
        help="Skip README.md in scan set.",
    )
    return parser.parse_args()


SKIP_DIRS = {".git", "docs_build", "obj", "bin", "stubs", "runs", "studies", "__pycache__", ".pytest_cache"}


def git_markdown(repo_root: Path) -> list[Path] | None:
    """!
    @brief List tracked plus non-ignored untracked Markdown files, so the scan set is
           reproducible for a given commit regardless of local build or scratch output.
    @param[in] repo_root Repository root directory.
    @return Sorted Markdown paths, or None when git enumeration is unavailable.
    """
    return enumerate_repository_files(repo_root, ".md")


def iter_markdown_files(repo_root: Path, include_readme: bool):
    """!
    @brief Yield every repository Markdown file the current commit would carry.
    @param[in] repo_root Repository root directory.
    @param[in] include_readme Whether the top-level README participates in the scan.
    @return Generator of Markdown paths.
    """
    candidates = git_markdown(repo_root)
    if candidates is None:
        # Fallback for a non-git export: walk the tree and skip generated/scratch trees.
        candidates = sorted(repo_root.rglob("*.md"))
    for markdown in candidates:
        if not markdown.is_file():
            continue
        parts = markdown.relative_to(repo_root).parts
        if SKIP_DIRS.intersection(parts):
            continue
        if not include_readme and markdown == repo_root / "README.md":
            continue
        yield markdown


def should_skip_link(target: str) -> bool:
    """!
    @brief Perform should skip link.
    @param[in] target Argument passed to `should_skip_link()`.
    @return Value returned by `should_skip_link()`.
    """
    lower = target.lower()
    return (
        lower.startswith("http://")
        or lower.startswith("https://")
        or lower.startswith("mailto:")
        or lower.startswith("#")
    )


def normalize_target(raw_target: str) -> str:
    """!
    @brief Remove Markdown angle brackets, fragments, and query strings before resolving a link target.
    @param[in] raw_target Argument passed to `normalize_target()`.
    @return Value returned by `normalize_target()`.
    """
    cleaned = raw_target.strip().strip("<>").split("#", 1)[0].split("?", 1)[0]
    return cleaned


def main() -> int:
    """!
    @brief Entry point for this script.
    @return Value returned by `main()`.
    """
    args = parse_args()
    repo_root = args.repo_root.resolve()
    failures = []

    scanned = list(iter_markdown_files(repo_root, include_readme=not args.no_readme))
    for md_file in scanned:
        if not md_file.is_file():
            failures.append((str(md_file), "<file>", "Markdown file not found"))
            continue
        content = md_file.read_text(encoding="utf-8", errors="replace")
        for match in LINK_PATTERN.findall(content):
            target = match[0] or match[1]
            if should_skip_link(target):
                continue
            normalized = normalize_target(target)
            if not normalized:
                continue
            resolved = (md_file.parent / normalized).resolve()
            try:
                exists = resolved.exists()
            except OSError as exc:
                failures.append((str(md_file.relative_to(repo_root)), target, f"{resolved} ({exc})"))
                continue
            if not exists:
                failures.append((str(md_file.relative_to(repo_root)), target, str(resolved)))

    if failures:
        print("Broken markdown links detected:")
        for src, target, resolved in failures:
            print(f"  - {src}: '{target}' -> missing '{resolved}'")
        return 1

    print(f"Markdown link check passed for {len(scanned)} repository Markdown files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
