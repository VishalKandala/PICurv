# Smoke Solver Fixtures

Solver profiles are derived from initialized example templates and rewritten
inside `tests/smoke/run_smoke.sh` to exercise:

- template-matrix dry-run planning for flat/bent/brownian examples,
- `solve` mode on tiny Eulerian and particle-enabled runs,
- `load` mode for restart branch coverage, restart-equivalence split runs, and MPI particle restarts,
- `analytical` mode (`ZERO_FLOW`) for Brownian verification.
