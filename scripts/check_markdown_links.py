#!/usr/bin/env python3
"""Basic local markdown link checker for README + authored docs pages."""

import re
import sys
from pathlib import Path


LINK_PATTERN = re.compile(r"!\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)|\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


def iter_markdown_files(repo_root: Path):
    yield repo_root / "README.md"
    for md in sorted((repo_root / "docs").rglob("*.md")):
        if "docs_build" in md.parts:
            continue
        yield md


def should_skip_link(target: str) -> bool:
    lower = target.lower()
    return (
        lower.startswith("http://")
        or lower.startswith("https://")
        or lower.startswith("mailto:")
        or lower.startswith("#")
    )


def normalize_target(raw_target: str) -> str:
    cleaned = raw_target.strip().strip("<>").split("#", 1)[0].split("?", 1)[0]
    return cleaned


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    failures = []

    for md_file in iter_markdown_files(repo_root):
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
            if not resolved.exists():
                failures.append((str(md_file.relative_to(repo_root)), target, str(resolved)))

    if failures:
        print("Broken markdown links detected:")
        for src, target, resolved in failures:
            print(f"  - {src}: '{target}' -> missing '{resolved}'")
        return 1

    print("Markdown link check passed for README.md and docs/**/*.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
