@page 40_Testing_and_Quality_Guide Testing and Validation Guide

This page documents PICurv's local testing model. The suite is intentionally split by intent so users can choose the smallest command that answers the question they actually have.

@tableofcontents

@section taxonomy_sec 1. Testing Taxonomy

PICurv exposes four validation layers plus one aggregate target:

1. Python/control-plane validation
2. installation and PETSc provisioning validation
3. isolated C unit/component validation
4. executable entrypoint smoke validation
5. full local validation sweep

Canonical commands:

- `make test-python`
- `make doctor`
- `make unit`
- `make unit-geometry`
- `make unit-solver`
- `make unit-particles`
- `make unit-io`
- `make unit-post`
- `make unit-grid`
- `make unit-metric`
- `make unit-boundaries`
- `make unit-poisson-rhs`
- `make unit-runtime`
- `make smoke`
- `make check`

Compatibility aliases remain available:

- `make test` -> `make test-python`
- `make install-check` -> `make doctor`
- `make ctest*` -> `make unit*`

@section quickstart_sec 2. Quick Start

Fast local checks:

```bash
make test
make doctor
make unit-io
make unit-runtime
make smoke
make check
```

Guidance:

- Use `make test` when working on `scripts/picurv`, schemas, or repository metadata.
- Use `make doctor` after provisioning PETSc on a new machine.
- Use `make unit-<area>` while changing a subsystem in isolation.
- Use `make smoke` after building binaries to execute tiny real solve/post/restart workflows.
- Use `make check` as the pre-merge gate at the end of a development cycle.

@section python_sec 3. Python Suite (`test-python`)

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

Run locally:

```bash
make test-python
```

@section doctor_sec 4. Installation Validation (`doctor`)

Purpose:

- confirm `PETSC_DIR` is present and usable
- confirm the configured toolchain can compile against PETSc
- confirm a minimal PETSc-backed binary can initialize and create core PETSc objects

The `doctor` target builds and runs a small C smoke binary under `tests/c/test_install_check.c`.

It validates:

- environment visibility for `PETSC_DIR`
- `PetscInitialize`
- `DMDA` creation
- `Vec` creation
- `DMSwarm` creation
- `PetscFinalize`

This target answers: "Is this machine set up correctly for PICurv development?"

Run locally:

```bash
make doctor
```

@section unit_sec 5. Isolated C Unit Suites (`unit-*`)

Purpose:

- support kernel and component development in isolation
- reduce the blast radius while changing numerical code

Current suites:

- `make unit-geometry`
- `make unit-solver`
- `make unit-particles`
- `make unit-io`
- `make unit-post`
- `make unit-grid`
- `make unit-metric`
- `make unit-boundaries`
- `make unit-poisson-rhs`
- `make unit-runtime`
- `make unit` (all of the above)

Suite focus areas:

- geometry and interpolation invariants
- solver utility kernels
- particle location helpers
- I/O pathing and field round-trips
- post-processing/statistics kernels
- grid decomposition and bounding-box exchange
- metric tensor/face metric kernels
- boundary factory/face-service checks
- Poisson/RHS and diffusivity kernels
- runtime helper kernels (initialization, particle helpers, wall models)

These tests use real PETSc objects and run single-rank by default.

@section smoke_sec 6. Executable Smoke (`smoke`)

Purpose:

- confirm the compiled binaries still launch successfully
- catch integration breakage that isolated unit tests may miss

The current smoke runner verifies:

- `bin/simulator` launches and responds to `-help`
- `bin/postprocessor` launches and responds to `-help`
- `bin/picurv init` creates a self-contained case with copied binaries + origin metadata
- `picurv run --dry-run --format json` emits a valid solve/post execution plan
- restart dry-run planning resolves `run_control.restart_from_run_dir` into the expected restart source directory
- tiny end-to-end solve + post run (flat channel)
- tiny end-to-end particle solve + post run (flat channel)
- restart execution with both particle modes (`load` and `init`)
- tiny end-to-end analytical Brownian run with particle VTP + MSD CSV output

These checks are intentionally tiny but execute real solver/postprocessor runtime paths.

Run locally:

```bash
make smoke
```

@section aggregate_sec 7. Full Validation (`check`)

`make check` is the top-level local validation sweep.

It runs, in order:

1. `make test-python`
2. `make doctor`
3. `make unit`
4. `make smoke`

Use it when you want maximum local confidence before ending a development cycle.

@section troubleshooting_sec 8. Common Failure Modes

Common issues:

- `PETSC_DIR` is unset:
  - `make doctor`, `make unit*`, `make smoke`, and `make check` will fail.
- `PETSC_ARCH` points to the wrong build:
  - compilation may fail while including PETSc make variables.
- `mpicc` or MPI runtime wrappers are missing:
  - PETSc-backed builds may fail even if Python tests still pass.
- a unit test fails after writing temporary files:
  - inspect `/tmp/picurv-test-*` for preserved intermediate artifacts.
- `pytest` is unavailable:
  - only the Python suite is blocked; PETSc-backed C tests are separate.

@section extend_sec 9. Extending The Test Suite

When adding new tests:

1. choose the narrowest valid layer (`test-python`, `doctor`, `unit-*`, or `smoke`)
2. reuse `tests/c/test_support.*` for PETSc-backed fixtures
3. keep temporary artifacts under `/tmp`
4. document new user-facing targets or workflows in the README and this guide

For detailed C test maintenance guidance, see **@subpage 51_C_Test_Suite_Developer_Guide**.

@section runtime_coverage_sec 10. Runtime Coverage Map

The smoke suite uses these named runtime sequences:

- `S1`: tiny flat-channel solve+post run
- `S2`: tiny flat-channel solve+post with particles enabled
- `S3`: tiny restart runs from `S2` with `particle_restart_mode=load` and `particle_restart_mode=init`
- `S4`: tiny Brownian analytical solve+post with particle outputs and MSD statistics

Runtime file coverage map (unit targets + runtime sequences):

- `src/AnalyticalSolutions.c`: `unit-solver`, `S4`
- `src/BC_Handlers.c`: `unit-boundaries`, `unit-runtime`, `S1`, `S2`, `S3`
- `src/BodyForces.c`: `unit-solver`, `unit-poisson-rhs`, `S1`, `S2`
- `src/Boundaries.c`: `unit-boundaries`, `S1`, `S2`, `S3`
- `src/Filter.c`: `unit-solver`
- `src/Metric.c`: `unit-metric`, `unit-grid`, `S1`, `S2`
- `src/ParticleMotion.c`: `unit-runtime`, `S2`, `S3`, `S4`
- `src/ParticlePhysics.c`: `unit-runtime`, `S2`, `S3`, `S4`
- `src/ParticleSwarm.c`: `unit-runtime`, `S2`, `S3`, `S4`
- `src/grid.c`: `unit-grid`, `S1`, `S2`, `S4`
- `src/initialcondition.c`: `unit-runtime`, `S1`, `S2`, `S4`
- `src/interpolation.c`: `unit-geometry`, `unit-particles`, `S2`, `S3`, `S4`
- `src/io.c`: `unit-io`, `S1`, `S2`, `S3`, `S4`
- `src/les.c`: `unit-runtime`
- `src/logging.c`: `S1`, `S2`, `S3`, `S4`
- `src/momentumsolvers.c`: `S1`, `S2`
- `src/particle_statistics.c`: `unit-post`, `S4`
- `src/poisson.c`: `unit-poisson-rhs`, `S1`, `S2`
- `src/postprocessing_kernels.c`: `unit-post`, `S1`, `S2`, `S4`
- `src/postprocessor.c`: `S1`, `S2`, `S3`, `S4`
- `src/rhs.c`: `unit-poisson-rhs`, `S1`, `S2`
- `src/runloop.c`: `S1`, `S2`, `S3`, `S4`
- `src/setup.c`: `S1`, `S2`, `S3`, `S4`
- `src/simulator.c`: `S1`, `S2`, `S3`, `S4`
- `src/solvers.c`: `S1`, `S2`
- `src/vtk_io.c`: `S1`, `S2`, `S3`, `S4`
- `src/walkingsearch.c`: `unit-geometry`, `unit-particles`, `S2`, `S3`, `S4`
- `src/wallfunction.c`: `unit-runtime`
