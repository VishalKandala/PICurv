#!/usr/bin/env python3
"""Reject generic documentation-expansion debris anywhere in the repository Markdown corpus."""

from __future__ import annotations

import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "tests" / "tooling" / "generic_expansion_contract.json"


def load_contract() -> dict:
    """!
    @brief Load the forbidden-signature contract.
    @return Parsed contract mapping.
    """
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def iter_markdown(contract: dict):
    """!
    @brief Yield every repository Markdown file inside the contract's scan scope.
    @param[in] contract Parsed contract mapping.
    @return Generator of Markdown paths.
    """
    excluded = set(contract["scan"]["exclude_dirs"])
    for path in sorted(REPO_ROOT.rglob("*.md")):
        if excluded.intersection(path.relative_to(REPO_ROOT).parts):
            continue
        yield path


def normalize(text: str) -> str:
    """!
    @brief Collapse whitespace so formatting changes cannot hide an exact template copy.
    @param[in] text Raw file text.
    @return Whitespace-normalized text.
    """
    return " ".join(text.split())


def scan(contract: dict) -> list[str]:
    """!
    @brief Collect every forbidden-signature hit in the Markdown corpus.
    @param[in] contract Parsed contract mapping.
    @return Human-readable violation lines.
    """
    signatures = [
        ("marker", contract["forbidden_markers"]),
        ("heading", contract["forbidden_headings"]),
        ("fragment", contract["forbidden_fragments"]),
    ]
    violations: list[str] = []
    for path in iter_markdown(contract):
        text = path.read_text(encoding="utf-8", errors="replace")
        collapsed = normalize(text)
        relative = path.relative_to(REPO_ROOT)
        for kind, patterns in signatures:
            for pattern in patterns:
                if pattern in text or normalize(pattern) in collapsed:
                    violations.append(f"{relative}: forbidden {kind}: {pattern[:72]}")
    return violations


def main() -> int:
    """!
    @brief Fail when generic expansion debris is present in the Markdown corpus.
    @return Process status code.
    """
    contract = load_contract()
    violations = scan(contract)
    if violations:
        print("Generic documentation-expansion debris detected:", file=sys.stderr)
        for violation in violations:
            print(f"  {violation}", file=sys.stderr)
        print(
            "\nThese signatures mark uniform agent-generated filler. Replace them with "
            "guidance specific to this page's subject, or remove the section.",
            file=sys.stderr,
        )
        return 1
    print(f"Generic-expansion audit passed: {sum(1 for _ in iter_markdown(contract))} Markdown files clean.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
