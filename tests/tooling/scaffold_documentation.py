#!/usr/bin/env python3
"""Generate documentation skeletons that satisfy the repository's own contracts."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PAGES_DIR = REPO_ROOT / "docs" / "pages"
FAMILIES_PATH = REPO_ROOT / "tests" / "tooling" / "capability_families.json"

PAGE_TYPES = {
    "tutorial": "Tutorial",
    "how-to": "How-to",
    "reference": "Reference",
    "explanation": "Explanation",
    "hub": "Hub",
}


# Every scaffold carries this until an author removes it. Structural correctness and
# completion are separate conditions: a skeleton is deliberately shaped to satisfy the
# structural audits, so without an explicit marker an untouched scaffold could be
# published as if it were finished.
INCOMPLETE_MARKER = "@par SCAFFOLD-INCOMPLETE"
INCOMPLETE_NOTE = (
    f"{INCOMPLETE_MARKER}\n"
    "This was generated from a scaffold and is not finished. Replace every TODO, then\n"
    "delete this block. Publication gates reject the marker, so an untouched scaffold\n"
    "cannot ship."
)


def identifier_slug(value: str) -> str:
    """!
    @brief Build a safe Doxygen/file identifier from a free-form title.

    @details Titles may contain characters that are illegal in a Doxygen identifier or
             that would imply a nested path in a filename. The display title is kept
             intact; only the identifier is sanitized.
    @param[in] value Free-form title.
    @return Safe identifier.
    @throws SystemExit when nothing usable remains.
    """
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    cleaned = re.sub(r"_{2,}", "_", cleaned)
    if not cleaned:
        raise SystemExit(
            f"Title {value!r} contains no characters usable in an identifier; "
            f"choose a title with letters or digits."
        )
    if cleaned[0].isdigit():
        cleaned = f"p_{cleaned}"
    return cleaned


def slug(value: str) -> str:
    """!
    @brief Anchor slug for a selector value.
    @param[in] value Public selector value.
    @return Lowercase underscore slug.
    """
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def next_page_number() -> int:
    """!
    @brief Lowest unused page number in docs/pages.
    @return Next available page number.
    """
    used = {
        int(match.group(1))
        for path in PAGES_DIR.glob("*.md")
        if (match := re.match(r"(\d+)_", path.name))
    }
    return max(used) + 1 if used else 1


def capability_entry(family_id: str, value: str) -> str:
    """!
    @brief Render a Tier-2 capability entry carrying every required contract part.
    @param[in] family_id Capability family identifier.
    @param[in] value Public selector value.
    @return Markdown text for the entry.
    """
    families = {f["id"]: f for f in json.loads(FAMILIES_PATH.read_text(encoding="utf-8"))["families"]}
    family = families.get(family_id)
    if family is None:
        raise SystemExit(
            f"Unknown family '{family_id}'. Known: {', '.join(sorted(families))}"
        )
    anchor = family["entry_anchor_prefix"] + slug(value)
    return f"""@subsection {anchor}_sub {value}

@anchor {anchor}

{INCOMPLETE_NOTE}

**Identity.** `{family['selector']}` = `{value}` -> TODO: generated flag -> TODO: C enum
-> @ref TODO_implementing_function.

**What it does.** TODO: the behavior, in a short paragraph.

**When to choose it.** TODO: compare explicitly against its siblings in this family.
An entry that does not compare is a description, not guidance.

**Parameters it owns.** TODO: options exclusive to this value, required and optional
separated. Options shared across the family belong on the family page.

**Interactions.** TODO: what it requires, conflicts with, or implies, including
anything enforced by validation.

**Diagnostics.** TODO: what it logs, and how a reader tells it is working. Name the
failure signature, not only the success one.

**Evidence.** TODO: facets with sources, or state "Implemented only." Record the same
sources in `value_metadata["{value}"]["evidence"]` in
`tests/tooling/capability_families.json` - the audit checks they match.

**Limitations.** TODO: what it cannot do.
"""


def alias_stub(family_id: str, value: str, canonical: str) -> str:
    """!
    @brief Render a deprecated-alias stub, which owes three parts rather than eight.
    @param[in] family_id Capability family identifier.
    @param[in] value Deprecated selector value.
    @param[in] canonical Canonical value it resolves to.
    @return Markdown text for the stub.
    """
    families = {f["id"]: f for f in json.loads(FAMILIES_PATH.read_text(encoding="utf-8"))["families"]}
    family = families.get(family_id)
    if family is None:
        raise SystemExit(f"Unknown family '{family_id}'. Known: {', '.join(sorted(families))}")
    known = set(family.get("value_metadata", {}))
    if known and canonical not in known:
        raise SystemExit(
            f"Canonical target '{canonical}' is not a declared value of '{family_id}'. "
            f"Known: {', '.join(sorted(known))}"
        )
    anchor = family["entry_anchor_prefix"] + slug(value)
    return f"""@subsection {anchor}_sub {value} (deprecated)

@anchor {anchor}

{INCOMPLETE_NOTE}

**Identity.** `{value}` - a **deprecated alias** that normalizes to `{canonical}`.

**Status.** Deprecated. TODO: why it still parses.

**Migration.** TODO: the exact edits to move to `{canonical}`.
"""


PAGE_SECTIONS = {
    "tutorial": [
        ("Outcome", "TODO: what the reader will have produced by the end. State it concretely."),
        ("Prerequisites", "TODO: software, build state, inputs, and the starting directory."),
        ("Steps", "TODO: the ordered path. Every command must run as written, in this order."),
        ("Confirm It Worked", "TODO: the exact files, log lines, or values that mean success."),
        ("If It Failed", "TODO: the most common failure for this path, and its fix."),
        ("Next", "TODO: where the reader should go now."),
    ],
    "how-to": [
        ("Goal", "TODO: the one task this page completes."),
        ("Before You Start", "TODO: what must already be true. Assume the reader can run PICurv."),
        ("Procedure", "TODO: the steps. Branch only where the situation genuinely differs."),
        ("Verify", "TODO: how the reader confirms the task is done."),
    ],
    "reference": [
        ("Scope", "TODO: exactly what this page covers, and what it does not."),
        ("Reference", "TODO: the structured inventory, schema, or option table. Prefer a generated\nfragment over a hand-maintained one."),
        ("Constraints", "TODO: validation rules, mutual exclusions, and required combinations."),
        ("Related Surfaces", "TODO: the other pages that own adjacent contracts."),
    ],
    "explanation": [
        ("Concepts", "TODO: the ideas a reader needs before the rationale makes sense."),
        ("Why It Works This Way", "TODO: the reasoning behind the current design."),
        ("Alternatives Considered", "TODO: what else was possible, and why it was not chosen."),
        ("Limitations", "TODO: what this design cannot do."),
    ],
    "hub": [
        ("Routes", "TODO: the pages this hub owns, adopted with @subpage. Routing only - a hub\nthat restates reference prose has become a duplicate of what it routes to."),
    ],
}


def suggested_filename(title: str, number: int) -> str:
    """!
    @brief Filename for a scaffolded page, using the same sanitized identifier as the page ID.
    @param[in] title Free-form page title.
    @param[in] number Page number.
    @return Safe filename.
    """
    return f"{number:02d}_{identifier_slug(title)}.md"


def page(kind: str, title: str, number: int) -> str:
    """!
    @brief Render a page skeleton shaped to its document type's contract.

    @details Each type gets the sections it actually owes, not a generic outline. A
             tutorial that does not state its outcome, or a reference that does not
             declare its scope, is not the type it claims to be.
    @param[in] kind Page type key.
    @param[in] title Human-readable page title.
    @param[in] number Page number.
    @return Markdown text for the page.
    """
    safe = identifier_slug(title)
    identifier = f"{number:02d}_{safe}"
    owes = {
        "tutorial": "a working outcome: every command runs as written, from a stated starting state",
        "how-to": "a completed task, assuming the reader can already run PICurv",
        "reference": "accuracy and completeness within its declared scope",
        "explanation": "understanding: why it is like this, and what the alternatives were",
        "hub": "complete and current routing, and nothing else",
    }[kind]

    body = [
        f"@page {identifier} {title}",
        "",
        f"@anchor _{safe}",
        "",
        f"@pagemeta{{{PAGE_TYPES[kind]}, TODO audience, TODO status}}",
        "",
        INCOMPLETE_NOTE,
        "",
        "TODO: one paragraph saying what this page is for.",
        "",
        f"This page owes the reader {owes}. See **@subpage 63_Page_Type_Contract** for what a",
        f"{PAGE_TYPES[kind]} must refuse to do.",
        "",
        "@tableofcontents",
        "",
    ]
    for index, (heading, guidance) in enumerate(PAGE_SECTIONS[kind], start=1):
        section = identifier_slug(heading).lower()
        body += [f"@section p{number}_{section}_sec {index}. {heading}", "", guidance, ""]
    last = len(PAGE_SECTIONS[kind]) + 1
    body += [
        f"@section p{number}_related_sec {last}. Related Documentation",
        "",
        "- **@subpage 47_Documentation_Catalog**",
        "",
    ]
    return "\n".join(body)


def subsystem_record(identifier: str, title: str) -> str:
    """!
    @brief Render the machine-readable lifecycle record a `planned` subsystem owes.

    @details The charter is prose; this is the record `make audit-subsystems` reads.
             A planned subsystem owes only four obligations, so the stub is short by
             design - the gate grows with the status, not with the scaffold.
    @param[in] identifier Subsystem id.
    @param[in] title Human-readable subsystem name.
    @return JSON text to paste into tests/tooling/subsystem_records.json.
    """
    record = {
        "id": identifier,
        "title": title,
        "status": "planned",
        "visibility": "public",
        "previous_status": None,
        "capability_families": [],
        "obligations": {
            "purpose": {"page": "TODO_page_id", "anchor": "TODO_anchor"},
            "intended_scope": {"page": "TODO_page_id", "anchor": "TODO_anchor"},
            "design_owner": {"value": "TODO: name"},
            "not_implemented_status": {"page": "TODO_page_id", "anchor": "TODO_anchor"},
        },
        "concerns": {},
    }
    return json.dumps(record, indent=2)


def subsystem_charter(identifier: str, title: str) -> str:
    """!
    @brief Render a subsystem charter at the `planned` lifecycle stage.

    @details Requirements grow with the status a subsystem claims, so this scaffold
             asks only what `planned` owes. The later stages are listed so the author
             can see what each promotion will require.
    @param[in] identifier Subsystem identifier.
    @param[in] title Human-readable subsystem name.
    @return Markdown text for the charter.
    """
    return f"""# Subsystem charter: {title}

{INCOMPLETE_MARKER} - replace every TODO, then delete this line.

Identifier: `{identifier}`
Lifecycle status: **planned**

## Purpose

TODO: what problem this solves, and for whom.

## Intended scope

TODO: what it will cover.

## Explicitly out of scope

TODO: what it will not cover. This is the most useful section to write early.

## Design owner

TODO: name.

## Not implemented

This subsystem is `planned`. Nothing here describes working behavior, and it must not
appear in user-facing reference until it reaches `public experimental`.

## Concern modules

Declare which apply, per **64_Documentation_Extension_Framework**. A concern may be
marked not applicable **with a stated reason**; an empty heading is not an answer.

- [ ] Numerical method
- [ ] User selector
- [ ] Units and nondimensionalization
- [ ] Artifact topology and storage lifecycle
- [ ] Persistent/restart state
- [ ] Determinism and reproducibility
- [ ] MPI/distributed execution
- [ ] External service
- [ ] Security/credentials
- [ ] Generated artifact
- [ ] File format
- [ ] Concurrency, permissions, destructive scope
- [ ] Backward compatibility
- [ ] Scientific verification and validation

## What each promotion will require

| To reach | Add |
|---|---|
| `internal` | Scope boundary, architecture boundary, dependencies, developer entry points |
| `public experimental` | Configuration, selection guidance, observability, limitations, safe-use boundaries |
| `supported` | Evidence, lifecycle and restart behavior, operations, troubleshooting, examples, complete reference |
"""


def main() -> int:
    """!
    @brief Emit a documentation skeleton for the requested kind.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(
        description="Generate documentation skeletons that satisfy the repository's contracts.",
        epilog=(
            "Examples:\n"
            "  scaffold_documentation.py capability --family boundary.handler --value slip_wall\n"
            "  scaffold_documentation.py alias --family momentum.solver --value 'Old Name' "
            "--canonical 'New Name'\n"
            "  scaffold_documentation.py page --type how-to --title 'Running On A Cluster'\n"
            "  scaffold_documentation.py subsystem --id ibm --title 'Immersed Boundaries'\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    sub = parser.add_subparsers(dest="kind", required=True)

    p_cap = sub.add_parser("capability", help="A Tier-2 entry for a new selector value.")
    p_cap.add_argument("--family", required=True)
    p_cap.add_argument("--value", required=True)

    p_alias = sub.add_parser("alias", help="A deprecated-alias stub.")
    p_alias.add_argument("--family", required=True)
    p_alias.add_argument("--value", required=True)
    p_alias.add_argument("--canonical", required=True)

    p_page = sub.add_parser("page", help="A new documentation page of a declared type.")
    p_page.add_argument("--type", required=True, choices=sorted(PAGE_TYPES))
    p_page.add_argument("--title", required=True)

    p_sub = sub.add_parser("subsystem", help="A subsystem charter at the planned stage.")
    p_sub.add_argument("--id", required=True)
    p_sub.add_argument("--title", required=True)

    args = parser.parse_args()

    if args.kind == "capability":
        print(capability_entry(args.family, args.value))
        print("# Next: paste into the family page, then run 'make audit-capability'.", file=sys.stderr)
    elif args.kind == "alias":
        print(alias_stub(args.family, args.value, args.canonical))
        print("# Next: declare alias_of in capability_families.json.", file=sys.stderr)
    elif args.kind == "page":
        number = next_page_number()
        print(page(args.type, args.title, number))
        print(
            f"# Next: save as docs/pages/{suggested_filename(args.title, number)}, "
            f"then adopt it with @subpage from a hub page or it will be reported as orphaned.",
            file=sys.stderr,
        )
    else:
        print(subsystem_charter(args.id, args.title))
        print("\n<!-- Paste into the \"subsystems\" list of "
              "tests/tooling/subsystem_records.json:\n")
        print(subsystem_record(args.id, args.title))
        print("-->")
        print("# Next: add the record above to tests/tooling/subsystem_records.json and run\n"
              "#       make audit-subsystems. Record it in the contract registry too if it\n"
              "#       introduces invariants.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
