# Contributing to PICurv

PICurv welcomes focused bug reports, documentation corrections, tests, and code
changes. The fastest contributions are narrow, reproducible, and explicit about what
was verified.

## Before opening an issue

For a solver, MPI, restart, or numerical problem, include enough provenance to
distinguish configuration, environment, and implementation effects:

- the PICurv commit and whether the worktree was dirty;
- the copied/effective run configuration, not only the source YAML;
- rank count and decomposition, plus PETSc, MPI, compiler, and launcher versions;
- whether the run was fresh or restarted and the exact restart segment/state;
- the smallest relevant monitor, solver, or runtime logs; and
- a minimal reproducer when practical.

Remove private paths, credentials, and cluster-account details before sharing output.
The [troubleshooting guide](docs/pages/67_Troubleshooting.md) describes the existing
diagnostics and safe first checks.

## Preparing a change

Start with the [Developer Portal](docs/pages/43_Developer_Portal_Index.md), then use
the nearest directory `guide.md` as a map. The implementation and runtime behavior are
the source of truth; generated inventories and registries enforce only the contracts
they explicitly declare.

Keep each pull request focused. Avoid combining behavior changes with unrelated file
moves, formatting, or refactoring. Before adding a helper, branch, state field, or
production file, inspect the current owner and its callers. Reuse or coherently
generalize existing infrastructure when its semantics and performance fit. Explain any
new abstraction or file and name the alternatives inspected.

Update tests and documentation in the same change when behavior or a public contract
changes. The [documentation extension framework](docs/pages/64_Documentation_Extension_Framework.md)
routes selector and subsystem obligations; the [testing guide](docs/pages/40_Testing_and_Quality_Guide.md)
and [test target map](tests/guide.md) route verification.

## Verification

Run the narrowest checks that directly exercise the changed layers, then broaden in
proportion to risk. Report every command run and its result. Also list applicable checks
that were not run or could not run. A passing registry, schema, or documentation audit
does not establish runtime or numerical correctness; MPI, restart, and performance
claims need corresponding evidence.

Use `make certify-docs-fast` or `make certify-docs` only on a clean commit, because the
certificate is revision-scoped. During development, run the constituent targets named
by the testing and documentation guides.

## Pull requests

A reviewable pull request should state:

- the problem and deliberately bounded scope;
- the change class and affected runtime/configuration/documentation contracts;
- reuse candidates inspected and the reason for any new helper or file;
- tests and audits run, including rank/restart/performance coverage where relevant;
- checks not run and behavior not verified; and
- documentation or freshness records changed.

Review findings should point to files and lines. Keep demonstrated defects separate
from plausible risks, and never describe an untested surface as verified.

## Agent-assisted development

Agents are optional. Repository-aware agents should read `AGENTS.md`; Claude Code also
loads `CLAUDE.md`, which imports the same contract. Reusable skills live canonically
under `.agents/skills/` and are materialized under `.claude/skills/` for portable
discovery. Run `make audit-agent-setup` to check parity or
`make sync-agent-skills` after intentionally editing a canonical skill.

Whether a change is written by a person or an agent, the same scope, reuse,
verification, and evidence rules apply.
