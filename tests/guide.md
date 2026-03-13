# Tests Guide

PICurv testing is intentionally layered: Python control-plane validation, PETSc installation checks, focused C unit suites, executable smoke runs, MPI variants, and coverage gates. Choose the narrowest layer that answers your current question.

Function-level docs in `tests/c/` are part of the repository contract. For new
or touched C tests/helpers, keep the Doxygen blocks concise but descriptive:
the `@brief` should state what the routine verifies or sets up, not just that it
is a test-local routine.

## Canonical Targets

- repository/doc contract:
  - `python3 scripts/audit_function_docs.py`
- Python and coverage:
  - `make test-python` (`make test` alias)
  - `make coverage-python`
  - `make coverage-c`
  - `make coverage`
- installation/toolchain:
  - `make doctor` (`make install-check` alias)
- C unit suites:
  - `make unit`
  - `make unit-geometry`
  - `make unit-setup`
  - `make unit-solver`
  - `make unit-particles`
  - `make unit-io`
  - `make unit-logging`
  - `make unit-post`
  - `make unit-grid`
  - `make unit-metric`
  - `make unit-boundaries`
  - `make unit-poisson-rhs`
  - `make unit-runtime`
  - `make unit-simulation`
  - `make unit-mpi`
  - `make unit-periodic-dev`
- smoke/integration:
  - `make smoke`
  - `make smoke-mpi`
  - `make smoke-mpi-matrix`
  - `make smoke-stress`
  - `make smoke-periodic-dev`
- aggregate gates:
  - `make check`
  - `make check-mpi`
  - `make check-mpi-matrix`
  - `make check-full`
  - `make check-stress`

## Python Test Files (`tests/test_*.py`)

- `test_cli_smoke.py`
  - CLI help and argument contract checks
  - dry-run plan schema checks (text/json)
  - staged Slurm workflow coverage for `submit`, `cancel`, and `sweep`
  - run summary coverage for `summarize` JSON/text output plus failure paths
  - restart path resolution checks
  - cluster no-submit manifest/script checks
  - grid-gen/PICGRID header and node-count translation checks
  - case-local binary preference behavior for copied/symlinked `picurv`
- `test_case_maintenance.py`
  - `init` origin metadata behavior
  - real CLI wrapper coverage for `build`, `sync-binaries`, `sync-config`, `status-source`, and `pull-source`
  - source-root resolution for build/sync/pull commands
  - template sync behavior (`overwrite`, `prune`)
  - source/case drift reporting (`status-source`)
- `test_config_regressions.py`
  - ingress-manifest drift checks
  - post recipe alias compatibility
  - post validation guards and statistics artifact pathing
- `test_repo_consistency.py`
  - validates example bundles and study bundles via `picurv validate`
  - scans docs/examples/tests for stale/forbidden contract literals
  - wraps the repository-wide function documentation audit script

## C Unit Files (`tests/c/test_*.c`)

- `test_install_check.c`: PETSc environment and basic object viability (`doctor`: `environment-visible`, `basic-petsc-objects`)
- `test_geometry.c`: interpolation and geometric signed-distance helpers
- `test_setup_lifecycle.c`: setup/cleanup lifecycle, RNG, and initialized particle-settlement contracts
- `test_solver_kernels.c`: analytical geometry/particle dispatch, LES filter/eddy-viscosity, FlowSolver guardrails, and analytical source helpers
- `test_particle_kernels.c`: walking-search helper kernels
- `test_io.c`: I/O path checks, parser helpers, scaling-ingestion contracts, and startup-banner summary contracts
- `test_logging.c`: log-level, allow-list, string-conversion, continuity/min-max/interpolation diagnostics, profiling, and snapshot-cadence contracts
- `test_postprocessing.c`: post-processing kernel contracts (specific-KE, displacement, nodal average, normalization, dimensionalization, Q-criterion)
- `test_vtk_io.c`: VTK writer and data-preparation contracts (coordinates, field gather/subsampling, particle prep)
- `test_postprocessor.c`: postprocessing orchestration contracts (swarm setup, pipeline dispatch, eulerian/particle output, statistics dispatch)
- `test_statistics.c`: statistics kernel contracts (MSD CSV output and empty-swarm behavior)
- `test_grid.c`: local/global bounding-box helpers
- `test_metric.c`: metric inversion, contravariant velocity, face geometry helpers
- `test_boundaries.c`: boundary factory plus direct handler-behavior checks
- `test_periodic_dev.c`: non-gating periodic geometric/driven boundary harnesses
- `test_poisson_rhs.c`: pressure update, RHS, projection, body-force and diffusivity helpers
- `test_runtime_kernels.c`: setup/runloop/particle/interpolation/scatter/wall/walltime-guard/LES helper contracts
- `test_mpi_kernels.c`: multi-rank particle distribution, bbox collectives, and restart migration behavior
- shared fixture layer:
  - `test_support.c`
  - `test_support.h`
  - exposes both a fast minimal fixture and a richer tiny-runtime fixture built through the real setup path
  - mirrors the production `da/fda/swarm` contract (`da = IM+1/JM+1/KM+1`, coordinate-DM `fda`, production swarm fields)

## Smoke Harness (`tests/smoke/run_smoke.sh`)

Single-rank smoke (`make smoke`) verifies:

- binary `-help` launch viability (`simulator`, `postprocessor`)
  - banner presence is authoritative; local PETSc may exit with code `62` (`PETSC_ERR_ARG_WRONG`) and still satisfy the smoke contract
- `picurv init` self-contained case creation and metadata
- template matrix `init + validate + dry-run` checks (`flat_channel`, `bent_channel`, `brownian_motion`)
- dry-run plan schema and restart-source resolution
- tiny real solve+post for flat and bent channels
- tiny particle solve+post and restart branches (`load`, `init`)
- restart-equivalence continuity check
- tiny analytical Brownian run with VTP + MSD CSV checks

Opt-in stress smoke (`make smoke-stress`) additionally verifies:

- longer particle-cycle runtime sequences
- chained restart workflows over more than one split point
- BC-focused runtime variants for the stable path (parabolic inlet) plus periodic constant-flux validate/dry-run coverage
- an extra-rank MPI particle runtime stress case

Periodic development smoke (`make smoke-periodic-dev`) is separate and non-gating:

- real runtime solve+post coverage for the current periodic constant-flux development path
- allowed to fail honestly while periodic BC support is still in development
- never part of `make check`, `make check-full`, or `make coverage-c`

Multi-rank smoke (`make smoke-mpi`, `make smoke-mpi-matrix`) additionally verifies:

- rank-dependent runtime launch behavior
- flat/bent multi-rank tiny solves
- particle restart branches under multi-rank execution

Useful env knobs:

- `TEST_MPI_NPROCS` for `unit-mpi`
- `SMOKE_MPI_NPROCS` for `smoke-mpi`
- `SMOKE_MPI_MATRIX_NPROCS` for `smoke-mpi-matrix`
- `KEEP_SMOKE_TMP=1` to preserve smoke temp workspace for debugging

## Suggested Command Cadence

- editing `scripts/picurv` or YAML contracts:
  - `make test-python`
- editing C/Python functions or helper/test docstrings:
  - `python3 scripts/audit_function_docs.py`
- editing one C subsystem:
  - targeted `make unit-<area>`
- editing setup/teardown lifecycle code:
  - `make unit-setup`
- editing core simulation numerics or particle orchestration:
  - `make unit-simulation`
- editing periodic BC code under development:
  - `make unit-periodic-dev`
  - `make smoke-periodic-dev`
- editing runtime orchestration, restart, or output contracts:
  - `make smoke` plus MPI variant if rank behavior is involved
- running the opt-in medium-budget extension tier:
  - `make smoke-stress`
- pre-merge:
  - `make check` (or `make check-mpi`)
- pre-release:
  - `make check-full`
  - `make coverage`

## Current targeted backlog

- walking search:
  - add direct unit coverage for `LocateParticleOrFindMigrationTarget` boundary-clamp, ghost-handoff, tie-breaker, `LOST`, and `MIGRATING_OUT` branches
  - add explicit direction-complete and failure-path checks for the `GuessParticleOwnerWithBBox` heuristic
- particle migration:
  - add non-restart MPI migration tests for multi-pass handoff, newcomer flagging, and particle-count conservation
- momentum:
  - add direct positive-path tests for `MomentumSolver_Explicit_RungeKutta4`
  - add one small direct invariant harness for `MomentumSolver_DualTime_Picard_RK4`
- pressure/Poisson:
  - add deeper `PoissonSolver_MG` and periodic/IBM stencil checks beyond the current `Projection`/`PoissonLHSNew` helper surface
- grid/metrics/setup:
  - broaden the richer runtime fixture to more geometry/topology variants; current setup coverage is production-faithful but still mostly tiny Cartesian cases
- periodic BC:
  - keep routing new periodic work into `unit-periodic-dev` and `smoke-periodic-dev` until the product runtime path is stable enough for the default gate

## Notes

- Python tests do not require PETSc.
- GitHub Actions quality CI runs `python scripts/audit_function_docs.py`, then `pytest -q`, then markdown link checks.
- `doctor`, `unit-*`, `smoke*`, `check*`, and `coverage-c` require PETSc/MPI tooling.
- `check-full` is the single-command comprehensive gate (`check` + `unit-mpi` + `smoke-mpi` + `smoke-mpi-matrix`).
- `check-stress` extends `check-full` with the opt-in `smoke-stress` layer.
- periodic development harnesses are intentionally non-gating.
- compatibility aliases (`install-check`, `ctest-*`) still exist, but canonical names are preferred in docs and CI.

## Authoritative Docs

- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/picurv-docs/51_C_Test_Suite_Developer_Guide.html
