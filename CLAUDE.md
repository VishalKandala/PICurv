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
| Solver and tooling logs | `logs/` — already gitignored |

Prefer the session scratchpad. Only put scratch in `sandbox/` when it needs to
outlive the session or be shared with a human.

`runs/`, `studies/`, and `logs/` are the product's own output paths, owned by
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
