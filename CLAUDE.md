# PICurv — working agreement for AI coding sessions

Applies to Claude Code, Codex, and any other agent working in this repo.
`AGENTS.md` is a symlink to this file, so both toolchains read the same rules.

## Never create new top-level directories

The repository root is a fixed, curated set. Do not add directories to it —
not for scratch work, not for results, not for notes, not "temporarily".
If you believe a new top-level directory is genuinely required, stop and ask
first; do not create it and mention it afterwards.

## Where work goes

| Kind of file | Location |
|---|---|
| Transient files with no lasting value (intermediate output, one-off scripts, parsing buffers) | the session scratchpad **outside** the repo |
| In-repo scratch: prototypes, exploratory notes, trial configs | `sandbox/` (see `sandbox/guide.md`) |
| Generated transient files that must sit in-repo | `sandbox/tmp/` |
| Solver run output | `runs/` — already gitignored |
| Sweep/study output | `studies/<study_id>/` — created by `picurv sweep`/`submit`, already gitignored |
| Solver and tooling logs | `<repo>/logs/` — already gitignored |

Prefer the session scratchpad. Only put scratch in `sandbox/` when it needs to
outlive the session or be shared with a human.

`runs/`, `studies/`, and `<repo>/logs/` are the product's own output paths, owned by
`picurv_cli`. Let the tooling create them. Do not hand-create them, write
unrelated files into them, or invent parallel alternatives.

## Sandbox hygiene

`sandbox/` is scratch, not a staging area for deliverables. Its contents are
gitignored except `guide.md`, so nothing there reaches a commit by accident —
which also means nothing there is backed up. Anything that matters gets
promoted out per the promotion rules in `sandbox/guide.md`, or it is lost.

Clean up after yourself: delete scratch you created once it has served its
purpose, rather than leaving it for the next session to inherit.

## Committing

Do not commit scratch, temp files, run output, or generated artifacts. Before
committing, check `git status` and confirm every path is one you intend to
track. Commit or push only when explicitly asked.

## Documentation: what to trust, and what you owe

### Trust hierarchy

When sources disagree, believe them in this order. Do not "fix" a lower tier to
match a higher one without checking which is actually wrong.

1. **Runtime and code behaviour** — what the program does.
2. **Enforced contracts** — schemas, manifests, and the registries under
   `tests/tooling/`. These are machine-checked, so they are rarely wrong for long.
3. **Generated documentation** — anything under `docs/generated/`. Never edit these
   by hand; regenerate with `make docs-inventory` or `make docs-topology`.
4. **Directory `guide.md` files** — local routing. Useful, not authoritative.
5. **Explanatory prose** — long-form pages. Most likely to have drifted.
6. **Proposals and history** — `57_Future_Architecture_Specifications`, changelog
   entries. These describe intent or the past, never current behaviour.

A page being well-written is not evidence that it is current. Verify any enumeration,
default, or ordering against tier 1 or 2 before repeating it.

### Before changing documentation

Run `make review-packet PAGE=<n>` or `make review-packet CONTRACT=<id>`. It assembles
the declared sources, owning and dependent contracts, capability families, evidence
sources, uncommitted changes, and how to verify — so you read the right few hundred
lines instead of grepping.

### If you add or change a subsystem's status

Statuses are `planned`, `internal`, `experimental`, `supported`, `known-defective`,
`deprecated`, `removed`. Each owes more documentation than the one below it, and the
ladder is cumulative. Update the record in `tests/tooling/subsystem_records.json` and
run `make audit-subsystems`; it will name what the new status owes and refuse a
transition the lifecycle does not allow. A reasoned "not applicable" is a real answer.
Empty prose written to pass the gate is not.

### If you add or change a user-selectable option

The documentation obligation is enforced, not advisory. See
`docs/pages/50_Modular_Selector_Extension_Guide.md` §1.1. In short: declare it in
`tests/tooling/capability_families.json`, scaffold its entry with
`tests/tooling/scaffold_documentation.py`, then `make docs-inventory` and
`make audit-capability`.

### Verification

```bash
make audit-capability      # capability parity, coverage, evidence, value lifecycle
make audit-subsystems      # what each subsystem's claimed status owes
make audit-page-types      # every published page is typed
make audit-freshness       # hard (blocking) and soft (advisory) staleness suspicion
make audit-contracts       # invariant contracts (run after build-docs)
make audit-path-literals   # no hardcoded run paths in prose
make preview-docs          # build, gate, and serve the site
make certify-docs          # everything, for a clean tree
```

`make audit-freshness` reports which pages need review because something they describe
changed. After you have actually compared a page with its sources, record it:

```bash
make attest-freshness ARGS="<surface-id>"
```

Attesting records that the comparison happened. It is not a claim that the page is
correct, and it is never a substitute for making the page correct.

### Two standing rules

**Never claim more than was checked.** "Built from commit X" is not "verified". A
declared evidence source is not a verified result. If automation cannot establish
something, say so rather than implying it did.

**Never write filler.** Guidance that would read the same on any page is not guidance.
`make audit-docs-expansion` rejects the specific filler that was removed from 51 pages;
the rule is broader than the check. A section you cannot fill with something specific
should not exist.
