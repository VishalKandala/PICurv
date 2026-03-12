# Tests Guide

PICurv testing is intentionally layered: Python control-plane validation, PETSc installation checks, focused C unit suites, executable smoke runs, MPI variants, and coverage gates. Choose the narrowest layer that answers your current question.

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
  - `make unit-mpi`
- smoke/integration:
  - `make smoke`
  - `make smoke-mpi`
  - `make smoke-mpi-matrix`
- aggregate gates:
  - `make check`
  - `make check-mpi`
  - `make check-mpi-matrix`
  - `make check-full`

## Python Test Files (`tests/test_*.py`)

- `test_cli_smoke.py`
  - CLI help and argument contract checks
  - dry-run plan schema checks (text/json)
  - restart path resolution checks
  - cluster no-submit manifest/script checks
  - grid-gen/PICGRID header and node-count translation checks
  - case-local binary preference behavior for copied/symlinked `picurv`
- `test_case_maintenance.py`
  - `init` origin metadata behavior
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

- `test_install_check.c`: PETSc environment and basic object viability (`doctor`)
- `test_geometry.c`: interpolation and geometric signed-distance helpers
- `test_solver_kernels.c`: LES filter/analytical source helpers
- `test_particle_kernels.c`: walking-search helper kernels
- `test_io.c`: I/O path checks, parser helpers, and scaling-ingestion contracts
- `test_logging.c`: log-level/allow-list/snapshot-cadence contracts
- `test_postprocessing.c`: post-processing kernel contracts (specific-KE, displacement, nodal average, normalization, dimensionalization, Q-criterion)
- `test_vtk_io.c`: VTK writer and data-preparation contracts (coordinates, field gather/subsampling, particle prep)
- `test_postprocessor.c`: postprocessing orchestration contracts (swarm setup, pipeline dispatch, eulerian/particle output, statistics dispatch)
- `test_statistics.c`: statistics kernel contracts (MSD CSV output and empty-swarm behavior)
- `test_grid.c`: local/global bounding-box helpers
- `test_metric.c`: metric inversion, contravariant velocity, face geometry helpers
- `test_boundaries.c`: boundary factory and face-service matrix checks
- `test_poisson_rhs.c`: pressure update, RHS, body-force and diffusivity helpers
- `test_runtime_kernels.c`: setup/runloop/particle/wall/LES helper contracts
- `test_mpi_kernels.c`: multi-rank particle distribution and bbox collective behavior
- shared fixture layer:
  - `test_support.c`
  - `test_support.h`

## Smoke Harness (`tests/smoke/run_smoke.sh`)

Single-rank smoke (`make smoke`) verifies:

- binary `-help` launch viability (`simulator`, `postprocessor`)
- `picurv init` self-contained case creation and metadata
- template matrix `init + validate + dry-run` checks (`flat_channel`, `bent_channel`, `brownian_motion`)
- dry-run plan schema and restart-source resolution
- tiny real solve+post for flat and bent channels
- tiny particle solve+post and restart branches (`load`, `init`)
- restart-equivalence continuity check
- tiny analytical Brownian run with VTP + MSD CSV checks

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
- editing runtime orchestration, restart, or output contracts:
  - `make smoke` plus MPI variant if rank behavior is involved
- pre-merge:
  - `make check` (or `make check-mpi`)
- pre-release:
  - `make check-full`
  - `make coverage`

## Notes

- Python tests do not require PETSc.
- GitHub Actions quality CI runs `python scripts/audit_function_docs.py`, then `pytest -q`, then markdown link checks.
- `doctor`, `unit-*`, `smoke*`, `check*`, and `coverage-c` require PETSc/MPI tooling.
- `check-full` is the single-command comprehensive gate (`check` + `unit-mpi` + `smoke-mpi` + `smoke-mpi-matrix`).
- compatibility aliases (`install-check`, `ctest-*`) still exist, but canonical names are preferred in docs and CI.

## Authoritative Docs

- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/picurv-docs/51_C_Test_Suite_Developer_Guide.html
