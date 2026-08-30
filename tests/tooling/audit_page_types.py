#!/usr/bin/env python3
"""Enforce that every published documentation page has a declared page type."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY = REPO_ROOT / "tests" / "tooling" / "page_types.json"
HTML_DIR = REPO_ROOT / "docs_build" / "html"
PAGE_DIRS = ("docs/pages", "docs")
# Authoring templates use a literal placeholder rather than a real page id.
PLACEHOLDER_IDS = {"<ID>", "<id>"}


def declared_pages() -> dict:
    """!
    @brief Every `@page` id in the documentation sources, with its file.
    @return Mapping of page id to source path.
    """
    pages: dict = {}
    for directory in PAGE_DIRS:
        for markdown in sorted((REPO_ROOT / directory).glob("*.md")):
            match = re.search(r"^@page\s+(\S+)", markdown.read_text(encoding="utf-8"), re.M)
            if match and match.group(1) not in PLACEHOLDER_IDS:
                pages.setdefault(match.group(1), markdown)
    return pages


def inline_type(path: Path) -> str:
    """!
    @brief The page type a page declares inline, if any.
    @param[in] path Page source.
    @return Declared type, or an empty string.
    """
    match = re.search(r"@pagemeta\{([^,}]+)", path.read_text(encoding="utf-8"))
    return match.group(1).strip() if match else ""


def main() -> int:
    """!
    @brief Fail when a published page is untyped or its declarations disagree.

    @details Coverage is enforced against the pages the build actually publishes. A
             central registry types all of them without adding repetitive chrome to each
             page; where a page also declares its type inline, the two must agree.
    @return Process status code.
    """
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    valid = set(registry["valid_types"])
    assignments = registry["assignments"]
    pages = declared_pages()

    published = {
        page_id for page_id in pages if (HTML_DIR / f"{page_id}.html").is_file()
    } or set(pages)

    problems: list = []
    for page_id in sorted(published):
        assigned = assignments.get(page_id)
        if not assigned:
            problems.append(f"{pages[page_id].name}: page '{page_id}' has no type assignment")
            continue
        if assigned not in valid:
            problems.append(f"{page_id}: type '{assigned}' is not one of {sorted(valid)}")
            continue
        declared = inline_type(pages[page_id])
        if declared and declared != assigned:
            problems.append(
                f"{pages[page_id].name}: inline @pagemeta says '{declared}' but the registry "
                f"says '{assigned}'"
            )
    for stale in sorted(set(assignments) - set(pages)):
        problems.append(f"registry assigns a type to '{stale}', which is not a declared page")

    if problems:
        print("Page-type coverage violations:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print(
            "\nAssign a type in tests/tooling/page_types.json. See 63_Page_Type_Contract for\n"
            "what each type owes the reader.",
            file=sys.stderr,
        )
        return 1

    from collections import Counter
    spread = Counter(assignments[p] for p in published)
    summary = ", ".join(f"{count} {kind}" for kind, count in sorted(spread.items()))
    print(f"Page-type audit passed: {len(published)} published pages typed ({summary}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
