# Smoke Solver Fixtures

Solver profiles used in smoke tests are derived from initialized examples and then rewritten in `tests/smoke/run_smoke.sh` to exercise critical runtime branches with minimal wall-clock cost.

## Scenario Coverage

- template-matrix dry-run planning for flat/bent/brownian examples,
- `solve` mode on tiny Eulerian and particle-enabled runs,
- `load` mode for restart branch coverage,
- restart-equivalence split runs,
- MPI particle restart behavior across `load` and `init` branches,
- `analytical` mode (`ZERO_FLOW`) for Brownian verification.

## Design Rationale

The smoke layer prioritizes branch/path coverage and orchestration stability over strict numerical fidelity. These tests should fail when runtime contracts break, even if the breakage appears in rarely used modes.

## Contributor Guidance

- Keep rewrites explicit and localized so scenario intent stays readable.
- When adding a new runtime branch in solver orchestration, add a smoke path that exercises it.
- Pair solver fixture changes with clear assertions in smoke output checks.
