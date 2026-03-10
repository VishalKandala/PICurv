#!/usr/bin/env python3
"""Basic local markdown link checker for README, docs, and example guides."""

import argparse
import re
import sys
from pathlib import Path


LINK_PATTERN = re.compile(r"!\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)|\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Check local markdown links for README.md plus docs/examples trees.\n"
            "HTTP(S), mailto, and in-page anchor links are skipped."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/check_markdown_links.py\n"
            "  python3 scripts/check_markdown_links.py --repo-root . --docs-dir docs --examples-dir examples\n"
            "  python3 scripts/check_markdown_links.py --no-readme"
        ),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root directory (default: parent of this script).",
    )
    parser.add_argument(
        "--docs-dir",
        default="docs",
        help="Docs subtree under repo root to scan recursively for *.md (default: docs).",
    )
    parser.add_argument(
        "--examples-dir",
        default="examples",
        help="Examples subtree under repo root to scan recursively for *.md (default: examples).",
    )
    parser.add_argument(
        "--no-readme",
        action="store_true",
        help="Skip README.md in scan set.",
    )
    return parser.parse_args()


def iter_markdown_files(repo_root: Path, docs_dir: str, examples_dir: str, include_readme: bool):
    """Yield markdown files from configured roots."""
    if include_readme:
        yield repo_root / "README.md"

    docs_root = repo_root / docs_dir
    examples_root = repo_root / examples_dir

    for md in sorted(docs_root.rglob("*.md")):
        if "docs_build" in md.parts:
            continue
        yield md
    for md in sorted(examples_root.rglob("*.md")):
        yield md


def should_skip_link(target: str) -> bool:
    """Perform should skip link."""
    lower = target.lower()
    return (
        lower.startswith("http://")
        or lower.startswith("https://")
        or lower.startswith("mailto:")
        or lower.startswith("#")
    )


def normalize_target(raw_target: str) -> str:
    """Normalize target."""
    cleaned = raw_target.strip().strip("<>").split("#", 1)[0].split("?", 1)[0]
    return cleaned


def main() -> int:
    """Entry point for this script."""
    args = parse_args()
    repo_root = args.repo_root.resolve()
    failures = []

    for md_file in iter_markdown_files(
        repo_root,
        docs_dir=args.docs_dir,
        examples_dir=args.examples_dir,
        include_readme=not args.no_readme,
    ):
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

    scope = []
    if not args.no_readme:
        scope.append("README.md")
    scope.append(f"{args.docs_dir}/**/*.md")
    scope.append(f"{args.examples_dir}/**/*.md")
    print(f"Markdown link check passed for {', '.join(scope)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
