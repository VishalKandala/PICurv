#!/usr/bin/env python3
"""Tiered documentation freshness: hard suspicion, soft suspicion, integrity failure.

A fingerprint answers one question: has the thing a page describes changed since a
human or agent last compared the page against it? It never claims the page is wrong,
and matching digests never claim the page is right. It routes review.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST = REPO_ROOT / "tests" / "tooling" / "freshness_manifest.json"
PAGE_TYPES = REPO_ROOT / "tests" / "tooling" / "page_types.json"

VALID_TIERS = ("hard", "soft")
VALID_ENFORCEMENT = ("blocking", "report")
SURFACE_KEYS = {"id", "title", "tier", "enforcement", "promotion_reason", "artifact",
                "regenerate", "watched_paths", "owning_pages", "dependent_pages",
                "attested_digest", "attested_scope", "note"}

UNATTESTED = "unattested"


def digest_of(paths: list) -> str:
    """!
    @brief Deterministic digest over an ordered list of files.

    @details The path is hashed alongside the bytes so that moving content between
             watched files is itself a change.
    @param[in] paths Repository-relative paths.
    @return Hex digest prefixed with its algorithm.
    """
    accumulator = hashlib.sha256()
    for relative in sorted(paths):
        accumulator.update(relative.encode("utf-8"))
        accumulator.update(b"\0")
        accumulator.update((REPO_ROOT / relative).read_bytes())
        accumulator.update(b"\0")
    return f"sha256:{accumulator.hexdigest()}"


def surface_paths(surface: dict) -> list:
    """!
    @brief The files a surface fingerprints.
    @param[in] surface One manifest entry.
    @return Repository-relative paths.
    """
    if surface["tier"] == "hard":
        return [surface["artifact"]]
    return list(surface.get("watched_paths", []))


def enforcement_of(surface: dict) -> str:
    """!
    @brief Whether a suspicion on this surface blocks.
    @param[in] surface One manifest entry.
    @return "blocking" or "report".
    """
    return surface.get("enforcement") or ("blocking" if surface["tier"] == "hard" else "report")


def validate_manifest(surfaces: list, published: set) -> list:
    """!
    @brief Structural checks on the manifest itself.

    @details A malformed manifest is an integrity failure, not ordinary staleness: it
             means coverage is broken and the absence of suspicion proves nothing.
    @param[in] surfaces Manifest entries.
    @param[in] published Ids of published pages.
    @return List of integrity failures.
    """
    failures: list = []
    seen: set = set()
    for surface in surfaces:
        identifier = surface.get("id")
        if not identifier:
            failures.append("a surface declares no id")
            continue
        if identifier in seen:
            failures.append(f"{identifier}: declared more than once")
        seen.add(identifier)
        unknown = set(surface) - SURFACE_KEYS
        if unknown:
            failures.append(f"{identifier}: unknown field(s) {sorted(unknown)}")
        if surface.get("tier") not in VALID_TIERS:
            failures.append(f"{identifier}: tier {surface.get('tier')!r} is not hard or soft")
            continue
        if enforcement_of(surface) not in VALID_ENFORCEMENT:
            failures.append(f"{identifier}: enforcement must be blocking or report")
        if surface["tier"] == "hard":
            if not surface.get("artifact"):
                failures.append(
                    f"{identifier}: a hard surface must name the normalized artifact it "
                    f"fingerprints"
                )
            if not surface.get("regenerate"):
                failures.append(
                    f"{identifier}: a hard surface must name the command that regenerates it"
                )
        else:
            if not surface.get("watched_paths"):
                failures.append(f"{identifier}: a soft surface must watch at least one path")
            if enforcement_of(surface) == "blocking" and not surface.get("promotion_reason"):
                failures.append(
                    f"{identifier}: a soft surface promoted to blocking must state why"
                )
        if not surface.get("owning_pages"):
            failures.append(
                f"{identifier}: no owning page; a fingerprint with nowhere to route its "
                f"suspicion cannot ask anyone to review anything"
            )
        for field in ("owning_pages", "dependent_pages"):
            for page in surface.get(field, []):
                if page not in published:
                    failures.append(f"{identifier}: {field} names '{page}', which is not published")
        for relative in surface_paths(surface):
            if not (REPO_ROOT / relative).is_file():
                failures.append(
                    f"{identifier}: watched path '{relative}' does not exist; coverage is "
                    f"broken, not merely stale"
                )
    return failures


def evaluate(surfaces: list, published: set) -> dict:
    """!
    @brief Classify every surface as current, suspect, or an integrity failure.
    @param[in] surfaces Manifest entries.
    @param[in] published Ids of published pages.
    @return Report with the four classifications and the pages needing review.
    """
    integrity = validate_manifest(surfaces, published)
    report = {"integrity_failures": integrity, "hard_current": [], "hard_suspect": [],
              "soft_current": [], "soft_suspect": [], "unattested": []}
    if integrity:
        return report
    for surface in surfaces:
        current = digest_of(surface_paths(surface))
        attested = surface.get("attested_digest")
        entry = {
            "id": surface["id"],
            "title": surface.get("title", surface["id"]),
            "tier": surface["tier"],
            "enforcement": enforcement_of(surface),
            "current_digest": current,
            "attested_digest": attested,
            "pages": list(surface.get("owning_pages", [])) +
                     list(surface.get("dependent_pages", [])),
            "regenerate": surface.get("regenerate"),
        }
        if attested in (None, UNATTESTED):
            report["unattested"].append(entry)
        elif attested == current:
            report[f"{surface['tier']}_current"].append(entry)
        else:
            report[f"{surface['tier']}_suspect"].append(entry)
    return report


def blocking_entries(report: dict) -> list:
    """!
    @brief The suspicions that must fail the build.
    @param[in] report Output of evaluate().
    @return Entries whose enforcement is blocking.
    """
    suspects = report["hard_suspect"] + report["soft_suspect"] + report["unattested"]
    return [entry for entry in suspects if entry["enforcement"] == "blocking"]


def attest(surfaces: list, wanted: list) -> list:
    """!
    @brief Record the current digests as reviewed.
    @param[in] surfaces Manifest entries, mutated in place.
    @param[in] wanted Surface ids, or ["all"].
    @return Ids that were re-attested.
    """
    changed = []
    for surface in surfaces:
        if wanted != ["all"] and surface["id"] not in wanted:
            continue
        current = digest_of(surface_paths(surface))
        if surface.get("attested_digest") != current:
            surface["attested_digest"] = current
            changed.append(surface["id"])
    return changed


def main() -> int:
    """!
    @brief Report freshness, or re-attest reviewed surfaces.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--attest", nargs="+", metavar="ID",
                        help="Record current digests as reviewed. Pass 'all' for every surface.")
    parser.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    parser.add_argument("--page", metavar="PAGE_ID",
                        help="Report only the surfaces that route review to this page.")
    arguments = parser.parse_args()

    document = json.loads(MANIFEST.read_text(encoding="utf-8"))
    surfaces = document["surfaces"]
    published = set(json.loads(PAGE_TYPES.read_text(encoding="utf-8"))["assignments"])

    if arguments.attest:
        failures = validate_manifest(surfaces, published)
        if failures:
            print("Cannot attest against a manifest with integrity failures:", file=sys.stderr)
            for failure in failures:
                print(f"  {failure}", file=sys.stderr)
            return 2
        changed = attest(surfaces, arguments.attest)
        MANIFEST.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
        if changed:
            print("Re-attested: " + ", ".join(changed))
            print("This records that the owning pages were compared against these sources.\n"
                  "It is not a claim that the comparison found nothing - only that it happened.")
        else:
            print("Nothing to re-attest; every requested surface already matches.")
        return 0

    report = evaluate(surfaces, published)
    if arguments.page:
        for key in ("hard_current", "hard_suspect", "soft_current", "soft_suspect", "unattested"):
            report[key] = [e for e in report[key] if arguments.page in e["pages"]]
    if arguments.json:
        print(json.dumps(report, indent=2))
        return 1 if report["integrity_failures"] or blocking_entries(report) else 0

    if report["integrity_failures"]:
        print("Freshness integrity failures (coverage is broken):", file=sys.stderr)
        for failure in report["integrity_failures"]:
            print(f"  {failure}", file=sys.stderr)
        return 2

    for entry in report["hard_suspect"]:
        print(f"HARD SUSPECT  {entry['id']} - {entry['title']}", file=sys.stderr)
        print(f"    the normalized artifact changed since it was last reviewed", file=sys.stderr)
        print(f"    review: {', '.join(entry['pages'])}", file=sys.stderr)
        print(f"    then:   python3 tests/tooling/audit_freshness.py --attest {entry['id']}",
              file=sys.stderr)
    for entry in report["soft_suspect"]:
        stream = sys.stderr if entry["enforcement"] == "blocking" else sys.stdout
        label = "SOFT SUSPECT (blocking)" if entry["enforcement"] == "blocking" else "soft suspect"
        print(f"{label}  {entry['id']} - {entry['title']}", file=stream)
        print(f"    watched sources changed; no semantic extractor exists for this surface",
              file=stream)
        print(f"    review: {', '.join(entry['pages'])}", file=stream)
    for entry in report["unattested"]:
        blocking_entry = entry["enforcement"] == "blocking"
        stream = sys.stderr if blocking_entry else sys.stdout
        label = "UNATTESTED (blocking)" if blocking_entry else "unattested"
        print(f"{label}  {entry['id']} - {entry['title']}: never reviewed against its "
              f"sources", file=stream)
        print(f"    review: {', '.join(entry['pages'])}", file=stream)

    counts = (f"{len(report['hard_current'])} hard-current, {len(report['hard_suspect'])} "
              f"hard-suspect, {len(report['soft_current'])} soft-current, "
              f"{len(report['soft_suspect'])} soft-suspect, "
              f"{len(report['unattested'])} never attested")
    blocking = blocking_entries(report)
    if blocking:
        print(f"\nFreshness: {counts}. Blocking on {len(blocking)} surface(s).", file=sys.stderr)
        print("A suspicion is a request to review, not a claim the page is wrong. Compare the\n"
              "page with its sources, correct what drifted, then re-attest.", file=sys.stderr)
        return 1
    print(f"Freshness: {counts}.")
    if report["soft_suspect"] or report["unattested"]:
        print("Advisory entries above do not block. A never-attested surface means no page has\n"
              "yet been compared against it - an honest gap, not a failure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
