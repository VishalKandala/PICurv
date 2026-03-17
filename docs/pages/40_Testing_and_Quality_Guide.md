@page 40_Testing_and_Quality_Guide Testing and Validation Guide

@anchor _Testing_and_Quality_Guide

This page documents PICurv's local testing model. The suite is intentionally split by intent so users can choose the smallest command that answers the question they actually have.

@tableofcontents

@section p40_taxonomy_sec 1. Testing Taxonomy

PICurv exposes four validation layers plus aggregate gates:

1. Python/control-plane validation
2. installation and PETSc provisioning validation
3. isolated C unit/component validation
4. executable entrypoint smoke validation
5. full local validation sweep
6. comprehensive MPI-inclusive validation sweep

Canonical commands:

- `make test-python`
- `make coverage-python`
- `make coverage-c`
- `make coverage`
- `make doctor`
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
- `make smoke`
- `make smoke-mpi`
- `make smoke-mpi-matrix`
- `make smoke-stress`
- `make smoke-periodic-dev`
- `make check`
- `make check-mpi`
- `make check-mpi-matrix`
- `make check-full`
- `make check-stress`

Compatibility aliases remain available:

- `make test` -> `make test-python`
- `make install-check` -> `make doctor`
- `make ctest*` -> `make unit*`

@section p40_quickstart_sec 2. Quick Start

Fast local checks:

```bash
python3 scripts/audit_function_docs.py
make test
make coverage-python
make doctor
make unit-setup
make unit-simulation
make unit-io
make unit-runtime
make unit-mpi
make smoke
make smoke-periodic-dev
make check
make check-full
```

Guidance:

- Use `python3 scripts/audit_function_docs.py` when changing C/Python function signatures, docstrings, or test helpers.
- Use `make test` when working on `scripts/picurv`, schemas, or repository metadata.
- Use `python3 scripts/check_markdown_links.py` when changing docs/examples.
- Use `make doctor` after provisioning PETSc on a new machine.
- Use `make unit-setup` when changing setup, teardown, initialization, or rank-info lifecycle code.
- Use `make unit-simulation` for the normal simulation-core debugging loop (`unit-boundaries + unit-solver + unit-poisson-rhs + unit-runtime + unit-particles`).
- Use `make unit-<area>` while changing a subsystem in isolation.
- Use `make smoke` after building binaries to execute tiny real solve/post/restart workflows.
- Use `make smoke-mpi-matrix` when you need a rank-sweep MPI runtime sanity check.
- Use `make smoke-stress` when you want the opt-in medium-budget runtime extension tier.
- Use `make unit-periodic-dev` and `make smoke-periodic-dev` only while working on the in-development periodic BC path; they are intentionally non-gating.
- Use `make coverage` to enforce line-coverage floors for core Python scripts and C sources.
- Use `make check` as the pre-merge gate at the end of a development cycle.
- Use `make check-mpi` when multi-rank MPI behavior is in scope.
- Use `make check-full` for comprehensive branch/CI/release validation that must cover all MPI layers.
- Use `make check-stress` when you want the full default gate plus the stress tier in one pass.

@section p40_python_sec 3. Python Suite (`test-python`)

Purpose:

- validate `scripts/picurv`
- validate config translation and contract behavior
- validate repository-level regression checks

This suite is implemented with `pytest` and does not require PETSc.

Current coverage includes:

- CLI help/validate/dry-run smoke
- config regression checks
- case maintenance regressions
- repository consistency scans
- function-level documentation contract checks for C and Python executable code

For newly added or modified C test helpers, keep Doxygen blocks concise but
descriptive: the `@brief` should state what the test or helper verifies rather
than using a placeholder label.

Run locally:

```bash
make test-python
```

CI note:

- `.github/workflows/quality.yml` runs `python scripts/audit_function_docs.py` explicitly before `pytest -q`, then runs markdown link checks on pull requests and pushes to `main`.

@subsection p40_python_files_ssec 3.1 Python Test File Matrix

Current Python files and primary responsibilities:

- `tests/test_cli_smoke.py`
  - CLI help and argument contract checks
  - dry-run JSON schema assertions
  - staged Slurm workflow coverage for `submit`, `cancel`, and `sweep`
  - `summarize` JSON/text output and failure-path checks
  - restart resolution/pathing checks
  - cluster no-submit manifest and sbatch-rank checks
  - grid-gen header/node-count checks
- `tests/test_case_maintenance.py`
  - case-origin metadata behavior
  - real CLI wrapper coverage for `build`, `sync-binaries`, `sync-config`, `status-source`, and `pull-source`
  - sync/pull/build source-root resolution behavior
  - template sync behavior (`overwrite`, `prune`)
  - status drift report behavior
- `tests/test_config_regressions.py`
  - ingress-manifest drift checks
  - post-recipe alias compatibility checks
  - post-task validation guards and stats artifact prediction
- `tests/test_repo_consistency.py`
  - validates canonical example bundles with `picurv validate`
  - scans docs/examples/tests for forbidden stale contract terms
  - runs the repository-wide function documentation audit (`scripts/audit_function_docs.py`)

These files are intentionally role-specific to keep failures actionable.

@section p40_doctor_sec 4. Installation Validation (`doctor`)

Purpose:

- confirm `PETSC_DIR` is present and usable
- confirm the configured toolchain can compile against PETSc
- confirm a minimal PETSc-backed binary can initialize and create core PETSc objects

The `doctor` target builds and runs a small C smoke binary under `tests/c/test_install_check.c`.

It validates:

- named case `environment-visible` for `PETSC_DIR` visibility
- named case `basic-petsc-objects` for `PetscInitialize`, `DMDA`, `Vec`, `DMSwarm`, and `PetscFinalize`

This target answers: "Is this machine set up correctly for PICurv development?"

Run locally:

```bash
make doctor
```

@section p40_unit_sec 5. Isolated C Unit Suites (`unit-*`)

Purpose:

- support kernel and component development in isolation
- reduce the blast radius while changing numerical code

Current suites:

- `make unit-geometry`
- `make unit-setup`
- `make unit-solver`
- `make unit-particles`
- `make unit-io`
- `make unit-post`
- `make unit-grid`
- `make unit-metric`
- `make unit-boundaries`
- `make unit-poisson-rhs`
- `make unit-runtime`
- `make unit-mpi`
- `make unit` (all of the above)

Suite focus areas:

- setup, initialization, and cleanup lifecycle coverage
- geometry and interpolation invariants
- solver utility kernels
- particle location helpers
- I/O pathing and field round-trips
- logging contracts (level/allow-list/console snapshot cadence)
- post-processing/statistics kernels and orchestration
- grid decomposition and bounding-box exchange
- metric tensor/face metric kernels
- boundary factory/face-service checks
- Poisson/RHS and diffusivity kernels
- runtime helper kernels (initialization, particle helpers, wall models)

These tests use real PETSc objects and run single-rank by default.
`unit-mpi` is the dedicated multi-rank suite (default `TEST_MPI_NPROCS=2`).

@subsection p40_c_file_matrix_ssec 5.1 C Test File Matrix

Current C test files and their main purpose:

- `tests/c/test_install_check.c`: PETSc installation viability (`doctor`: `environment-visible`, `basic-petsc-objects`)
- `tests/c/test_geometry.c`: interpolation/signed-distance helpers
- `tests/c/test_setup_lifecycle.c`: setup/cleanup lifecycle, RNG, and initialized particle-settlement checks
- `tests/c/test_solver_kernels.c`: analytical geometry/particle dispatch, LES filter/eddy-viscosity, and FlowSolver guardrails
- `tests/c/test_particle_kernels.c`: walking-search helper kernels
- `tests/c/test_io.c`: I/O path, parser, scaling-ingestion, and startup-banner checks
- `tests/c/test_logging.c`: logging-contract checks (log level, allow-list, continuity/min-max/interpolation diagnostics, string conversion, profiling, snapshot cadence)
- `tests/c/test_postprocessing.c`: post kernel checks (specific KE, displacement, nodal averaging, normalization, dimensionalization, Q-criterion)
- `tests/c/test_vtk_io.c`: VTK prep/writer checks (Eulerian + particle data shaping)
- `tests/c/test_postprocessor.c`: postprocessor pipeline/orchestration checks
- `tests/c/test_statistics.c`: statistics-kernel checks (MSD CSV behavior)
- `tests/c/test_grid.c`: bounding-box exchange checks
- `tests/c/test_metric.c`: metric and face-geometry checks
- `tests/c/test_boundaries.c`: boundary factory plus direct handler-behavior checks
- `tests/c/test_poisson_rhs.c`: pressure/rhs/projection/diffusivity helper checks
- `tests/c/test_runtime_kernels.c`: runtime orchestration, interpolation/scatter, particle lifecycle, wall, and walltime-guard checks
- `tests/c/test_mpi_kernels.c`: dedicated multi-rank consistency and restart-migration checks
- `tests/c/test_support.{c,h}`: shared PETSc fixture and assertion layer

@section p40_smoke_sec 6. Executable Smoke (`smoke`)

Purpose:

- confirm the compiled binaries still launch successfully
- catch integration breakage that isolated unit tests may miss

The current smoke runner verifies:

- `bin/simulator` launches and responds to `-help`
- `bin/postprocessor` launches and responds to `-help`
- `bin/picurv init` creates a self-contained case with copied binaries + origin metadata
- template matrix init/validate/dry-run checks across `flat_channel`, `bent_channel`, and `brownian_motion`
- `picurv run --dry-run --format json` emits a valid solve/post execution plan
- restart dry-run planning resolves `--restart-from` into the expected restart source directory
- tiny end-to-end solve + post run (flat channel)
- tiny end-to-end solve + post run (bent channel)
- tiny end-to-end particle solve + post run (flat channel, default Trilinear interpolation)
- tiny end-to-end particle solve + post run with corner-averaged interpolation (legacy regression)
- restart execution with both particle modes (`load` and `init`)
- restart-equivalence check for flat channel (continuous tiny run vs split restart tiny run continuity metric agreement)
- tiny end-to-end analytical Brownian run with particle VTP + MSD CSV output
- multi-rank tiny solve + post runs for flat and bent channels (`make smoke-mpi`)
- multi-rank flat particle runtime + restart (`load`/`init`) runs (`make smoke-mpi`)
- rank-matrix MPI runtime sweep across flat+bent+flat-particle-restart (`make smoke-mpi-matrix`)
- opt-in stress extensions for longer particle cycling, chained restarts, parabolic-inlet runtime coverage, periodic constant-flux validate/dry-run coverage, and extra-rank MPI particle runs (`make smoke-stress`)
- periodic runtime development harness for the in-progress periodic constant-flux path (`make smoke-periodic-dev`, non-gating)

These checks are intentionally tiny but execute real solver/postprocessor runtime paths.
For the `-help` smoke checks, banner presence is authoritative: the local PETSc
debug build may exit with code `62` (`PETSC_ERR_ARG_WRONG`) after printing help,
and that still counts as a pass.

Run locally:

```bash
make smoke
make smoke-mpi
make smoke-mpi-matrix
make smoke-stress
make smoke-periodic-dev
```

@subsection p40_smoke_knobs_ssec 6.1 Useful Smoke Knobs

Environment controls used by the smoke layer:

- `KEEP_SMOKE_TMP=1`:
  - preserves the smoke temporary workspace for post-failure inspection.
- `SMOKE_MPI_NPROCS=<n>`:
  - rank count for `make smoke-mpi`.
- `SMOKE_MPI_MATRIX_NPROCS="<n1> <n2> ..."`:
  - rank matrix for `make smoke-mpi-matrix`.

Smoke orchestration lives in `tests/smoke/run_smoke.sh`; runtime profile mutation helpers there are part of the tested contract.

@section p40_aggregate_sec 7. Aggregate Validation Gates (`check*`)

`make check` is the top-level local validation sweep.

It runs, in order:

1. `make test-python`
2. `make doctor`
3. `make unit`
4. `make smoke`

Use it when you want maximum local confidence before ending a development cycle.

`make check-mpi` extends this sweep by running `make smoke-mpi` and `make unit-mpi` afterward.

`make check-mpi-matrix` extends this sweep by running `make smoke-mpi-matrix` and `make unit-mpi`.

`make check-full` is the comprehensive MPI-inclusive gate. It runs `make check`, then `make unit-mpi`, `make smoke-mpi`, and `make smoke-mpi-matrix`.

`make check-stress` extends `make check-full` by running `make smoke-stress` afterward.
`make smoke-periodic-dev` is intentionally outside the default gates and may fail honestly while periodic BC work is still in development.

Use `check-full` for release candidates or CI workflows where both focused multi-rank tests and rank-matrix runtime checks are required in a single pass.

@section p40_coverage_impl_sec 8. Coverage Gate Implementation Notes

Coverage gates are script-backed and checked in-repo:

- `make coverage-python`:
  - runs `scripts/python_coverage_gate.py`
  - uses stdlib `trace` and currently focuses on core runtime scripts (default: `scripts/picurv`)
- `make coverage-c`:
  - rebuilds with coverage flags, runs `unit + smoke`, then executes `scripts/c_coverage_gate.py`
  - parses gcov output for weighted line coverage over `src/*.c`

Coverage artifacts:

- Python summary: `coverage/python/summary.txt`
- C summary: `coverage/c/summary.txt`

@section p40_troubleshooting_sec 9. Common Failure Modes

Common issues:

- `PETSC_DIR` is unset:
  - `make doctor`, `make unit*`, `make smoke*`, and `make check*` will fail.
- `PETSC_ARCH` points to the wrong build:
  - compilation may fail while including PETSc make variables.
- `mpicc` or MPI runtime wrappers are missing:
  - PETSc-backed builds may fail even if Python tests still pass.
- a unit test fails after writing temporary files:
  - inspect `/tmp/picurv-test-*` for preserved intermediate artifacts.
- `pytest` is unavailable:
  - only the Python suite is blocked; PETSc-backed C tests are separate.
- `gcov` is unavailable:
  - `make coverage-c` cannot produce C line-coverage reports.

@section p40_extend_sec 10. Extending The Test Suite

When adding new tests:

1. choose the narrowest valid layer (`test-python`, `doctor`, `unit-*`, or `smoke`)
2. reuse `tests/c/test_support.*` for PETSc-backed fixtures; it now provides both a fast minimal fixture and a richer tiny-runtime fixture while mirroring the production `da/fda/swarm` contract instead of a same-size synthetic DM
3. keep temporary artifacts under `/tmp`
4. document new user-facing targets or workflows in the README and this guide

For detailed C test maintenance guidance, see **@subpage 51_C_Test_Suite_Developer_Guide**.

@section p40_runtime_coverage_sec 11. Runtime Coverage Map

The smoke suite uses these named runtime sequences:

- `S0`: template matrix init/validate/dry-run across `flat_channel`, `bent_channel`, and `brownian_motion`
- `S1`: tiny flat-channel solve+post run
- `S1b`: tiny bent-channel solve+post run
- `S2`: tiny flat-channel solve+post with particles enabled (default Trilinear interpolation)
- `S2b`: tiny flat-channel particle solve+post with `CornerAveraged` interpolation (legacy path regression)
- `S3`: tiny restart runs from `S2` with `particle_restart_mode=load` and `particle_restart_mode=init`
- `S4`: tiny Brownian analytical solve+post with particle outputs and MSD statistics
- `S5`: multi-rank tiny solve+post runs for flat and bent channels, plus flat particle base/restart (`load` and `init`) branches
- `S6`: restart-equivalence run for flat channel (continuous vs split restart continuity metric agreement)
- `S7`: periodic constant-flux validate + dry-run contract coverage in the opt-in stress tier
- `S8`: periodic constant-flux real runtime development harness (`make smoke-periodic-dev`, non-gating)

Runtime file coverage map (unit targets + runtime sequences):

- `src/AnalyticalSolutions.c`: `unit-solver`, `S4`
- `src/BC_Handlers.c`: `unit-boundaries`, `unit-periodic-dev`, `unit-runtime`, `S1`, `S1b`, `S2`, `S3`, `S5`, `S6`, `S7`, `S8`
- `src/BodyForces.c`: `unit-solver`, `unit-poisson-rhs`, `S1`, `S2`
- `src/Boundaries.c`: `unit-boundaries`, `unit-periodic-dev`, `S1`, `S1b`, `S2`, `S3`, `S5`, `S6`, `S7`, `S8`
- `src/Filter.c`: `unit-solver`
- `src/Metric.c`: `unit-metric`, `unit-grid`, `S1`, `S1b`, `S2`, `S5`, `S6`
- `src/ParticleMotion.c`: `unit-runtime`, `S2`, `S3`, `S4`, `S5`
- `src/ParticlePhysics.c`: `unit-runtime`, `S2`, `S3`, `S4`, `S5`
- `src/ParticleSwarm.c`: `unit-runtime`, `S2`, `S3`, `S4`, `S5`
- `src/grid.c`: `unit-grid`, `unit-setup`, `S1`, `S1b`, `S2`, `S4`, `S5`, `S6`, `S8`
- `src/initialcondition.c`: `unit-runtime`, `S1`, `S1b`, `S2`, `S4`, `S5`, `S6`
- `src/interpolation.c`: `unit-geometry`, `unit-runtime`, `unit-particles`, `S2`, `S2b`, `S3`, `S4`, `S5`
- `src/io.c`: `unit-io`, `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`
- `src/les.c`: `unit-solver`, `unit-runtime`
- `src/logging.c`: `unit-logging`, `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`, `S8`
- `src/momentumsolvers.c`: `S1`, `S1b`, `S2`, `S5`, `S6`
- `src/particle_statistics.c`: `unit-post`, `S4`
- `src/poisson.c`: `unit-poisson-rhs`, `S1`, `S1b`, `S2`, `S5`, `S6`
- `src/postprocessing_kernels.c`: `unit-post`, `S1`, `S1b`, `S2`, `S4`, `S5`, `S6`
- `src/postprocessor.c`: `unit-post`, `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`
- `src/rhs.c`: `unit-poisson-rhs`, `S1`, `S1b`, `S2`, `S5`, `S6`
- `src/runloop.c`: `unit-runtime`, `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`
- `src/setup.c`: `unit-setup`, `unit-runtime`, `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`, `S8`
- `src/simulator.c`: `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`
- `src/solvers.c`: `unit-solver`, `S1`, `S1b`, `S2`, `S5`, `S6`
- `src/vtk_io.c`: `S1`, `S1b`, `S2`, `S3`, `S4`, `S5`, `S6`
- `src/walkingsearch.c`: `unit-geometry`, `unit-particles`, `S2`, `S3`, `S4`
- `src/wallfunction.c`: `unit-runtime`

@section p40_exhaustive_backlog_sec 12. Exhaustive-Readiness Backlog

P0 (implemented):

- coverage gates: `make coverage-python`, `make coverage-c`, `make coverage`
- restart-equivalence runtime smoke (`S6`)
- rank-sweep MPI runtime smoke (`make smoke-mpi-matrix`)

P1 (implemented):

- broadened MPI matrix to particle/restart branches
- dedicated setup/runloop orchestration branch tests in `unit-runtime`
- boundary-condition matrix expansion for handler/face/orientation combinations in `unit-boundaries`

P1 (next):

- direct walking-search branch pinning for `LocateParticleOrFindMigrationTarget`:
  - boundary clamp
  - ghost-region handoff
  - tie-breaker
  - `LOST` and `MIGRATING_OUT` outcomes
- explicit direction-complete and failure-path coverage for the `GuessParticleOwnerWithBBox` heuristic
- non-restart MPI particle-migration tests for multi-pass handoff, newcomer flagging, and count conservation
- direct positive-path momentum harnesses:
  - `MomentumSolver_Explicit_RungeKutta4`
  - one small invariant case for `MomentumSolver_DualTime_Picard_RK4`
- deeper `PoissonSolver_MG` and periodic/IBM stencil behavior checks beyond the current helper-level `unit-poisson-rhs` surface
- broader richer-runtime fixture variants so grid/setup/metric tests cover more than the tiny Cartesian baseline
- broaden the MPI rank matrix to larger optional decompositions (for example `SMOKE_MPI_MATRIX_NPROCS="2 3 4 6"`) in CI/nightly profiles

P2 (deeper hardening):

- long-duration nightly/weekly stability suites (beyond tiny smoke budgets)
- numerical oracle/golden-output tolerance checks for more physical scenarios
- performance regression gates (runtime/memory envelopes)

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Testing and Validation Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
