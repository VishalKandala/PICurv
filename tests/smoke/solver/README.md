# Smoke Solver Fixtures

Solver profiles are derived from initialized example templates and rewritten
inside `tests/smoke/run_smoke.sh` to exercise:

- `solve` mode on tiny Eulerian and particle-enabled runs,
- `load` mode for restart branch coverage,
- `analytical` mode (`ZERO_FLOW`) for Brownian verification.
