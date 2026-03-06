# Smoke Case Fixtures

The smoke harness materializes tiny runtime cases dynamically rather than relying on static prebuilt artifacts. This keeps fixture logic aligned with current templates and catches regressions in initialization behavior directly.

## Initialization Sources

- `picurv init flat_channel`
- `picurv init bent_channel`
- `picurv init brownian_motion`

Each initialized case is then validated and dry-run checked as part of a template-matrix coverage phase before runtime rewrites are applied.

## Why Dynamic Case Materialization Is Used

- It exercises real `init` logic, not a stale fixture snapshot.
- It validates template copy behavior and origin metadata generation.
- It ensures smoke tests evolve with template contracts.
- It keeps runtime artifacts deterministic and small.

## Runtime Rewrite Scope

Selected files are rewritten in-place for tiny real executions that cover:

- short solve+post runs,
- restart branch checks,
- restart-equivalence checks (continuous vs split restart),
- multi-rank particle restart branch checks (`load` and `init`).

## Maintenance Notes

- These fixtures are workflow-regression assets, not long-run physics validation cases.
- If template contracts change, synchronize this fixture layer and smoke assertions in the same commit.
- Keep runtime durations minimal to preserve fast local/CI feedback.
