#!/usr/bin/env python3
"""Run the invariant-contract registry: verify enforced contracts and report tracked ones."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = REPO_ROOT / "tests" / "tooling" / "contract_registry.json"
PAGES_DIR = REPO_ROOT / "docs" / "pages"


def load_registry() -> dict:
    """!
    @brief Load the invariant-contract registry.
    @return Parsed registry mapping.
    """
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


def validate_records(contracts: list, registry: dict) -> list[str]:
    """!
    @brief Verify each contract record is complete, uses the closed vocabularies, and
           points at things that exist.

    This fails closed on purpose. A misspelled status such as "enfroced" previously
    dropped a contract out of enforcement *and* out of every report, so the registry
    silently shrank while still looking healthy.
    @param[in] contracts Contract records.
    @param[in] registry Full registry, carrying the vocabularies and required fields.
    @return Violation lines.
    """
    vocab = registry["vocabularies"]
    required = registry["required_fields"]
    violations: list[str] = []
    seen: set = set()
    for index, contract in enumerate(contracts):
        cid = contract.get("id") or f"<record {index}>"
        if not contract.get("id"):
            violations.append(f"{cid}: record has no id")
        elif cid in seen:
            violations.append(f"{cid}: duplicate contract id")
        seen.add(cid)

        for field in required:
            if field not in contract:
                violations.append(f"{cid}: missing required field '{field}'")

        for field in ("status", "kind", "enforcement"):
            value = contract.get(field)
            if value is not None and value not in vocab[field]:
                violations.append(
                    f"{cid}: {field} '{value}' is not in the closed vocabulary {vocab[field]}"
                )

        status = contract.get("status")
        enforcement = contract.get("enforcement")
        if status == "enforced" and enforcement != "blocking":
            violations.append(f"{cid}: status enforced requires enforcement 'blocking', got '{enforcement}'")
        if status != "enforced" and enforcement == "blocking":
            violations.append(f"{cid}: enforcement 'blocking' requires status 'enforced', got '{status}'")

        checker = contract.get("checker")
        if status == "enforced" and not checker:
            violations.append(f"{cid}: status is enforced but no checker is registered")
        if checker is not None:
            if not isinstance(checker, list) or not checker or not isinstance(checker[0], str):
                violations.append(f"{cid}: checker must be a non-empty list whose first element is a script path")
            elif not (REPO_ROOT / checker[0]).is_file():
                violations.append(f"{cid}: checker script '{checker[0]}' does not exist")
            elif not all(isinstance(part, str) for part in checker):
                violations.append(f"{cid}: checker command contains a non-string argument")

        if status != "planned" and not contract.get("authoritative_sources"):
            violations.append(f"{cid}: names no authoritative sources")
        for source in contract.get("authoritative_sources", []):
            if not (REPO_ROOT / source).exists():
                violations.append(f"{cid}: authoritative source '{source}' does not exist")
        for page in contract.get("canonical_documentation", []) + contract.get("dependent_pages", []):
            if not (PAGES_DIR / f"{page}.md").is_file():
                violations.append(f"{cid}: references page '{page}', which does not exist")
    return violations


def run_checkers(contracts: list) -> tuple:
    """!
    @brief Execute each enforced contract's checker, honouring declared prerequisites.

    A checker that reads the generated HTML cannot run on a clean checkout. Skipping it
    with an explicit note is honest; running it and passing because a stale build
    happens to exist is not.
    @param[in] contracts Contract records.
    @return Failure lines, ids that passed, and skip notes.
    """
    failures: list = []
    passed: list = []
    skipped: list = []
    html_ready = (REPO_ROOT / "docs_build" / "html" / "index.html").is_file()
    for contract in contracts:
        checker = contract.get("checker")
        if contract["status"] != "enforced" or not checker:
            continue
        if contract.get("requires_built_docs") and not html_ready:
            skipped.append(
                f"{contract['id']}: requires the built site; run 'make build-docs' first"
            )
            continue
        command = [sys.executable, str(REPO_ROOT / checker[0]), *checker[1:]]
        result = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            detail = (result.stderr or result.stdout).strip().replace("\n", "\n        ")
            failures.append(f"{contract['id']}: checker failed\n        {detail}")
        else:
            passed.append(contract["id"])
    return failures, passed, skipped


def main() -> int:
    """!
    @brief Verify enforced invariant contracts and report the tracked ones.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(description="Run the invariant-contract registry.")
    parser.add_argument("--list", action="store_true", help="List contracts without running checkers.")
    args = parser.parse_args()

    registry = load_registry()
    contracts = registry["contracts"]

    violations = validate_records(contracts, registry)
    if args.list:
        for contract in sorted(contracts, key=lambda c: (c["status"], c["id"])):
            print(f"  {contract['status']:10} {contract['kind']:18} {contract['id']}")
        return 1 if violations else 0

    failures, passed, skipped = ([], [], []) if violations else run_checkers(contracts)

    tracked = [c for c in contracts if c["status"] == "tracked"]
    planned = [c for c in contracts if c["status"] == "planned"]
    if tracked or planned:
        print("Invariant contracts without an automated checker (locatable, not verified):")
        for contract in tracked + planned:
            print(f"  {contract['status']:8} {contract['id']} -> owns {', '.join(contract['canonical_documentation']) or 'nothing'}")
        print("")

    if skipped:
        print("Enforced contracts skipped (prerequisite missing):")
        for note in skipped:
            print(f"  {note}")
        print("")

    if violations or failures:
        print("Invariant contract violations:", file=sys.stderr)
        for line in violations + failures:
            print(f"  {line}", file=sys.stderr)
        return 1

    enforced_total = len([c for c in contracts if c["status"] == "enforced"])
    print(
        f"Invariant contract audit passed: {len(passed)}/{enforced_total} enforced contract(s) "
        f"verified"
        + (f", {len(skipped)} skipped" if skipped else "")
        + f"; {len(tracked)} tracked, {len(planned)} planned (not verified)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
