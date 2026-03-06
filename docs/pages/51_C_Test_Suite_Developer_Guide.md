@page 51_C_Test_Suite_Developer_Guide C Test Suite Developer Guide

@anchor _C_Test_Suite_Developer_Guide

This page documents how the PETSc-backed C testing layer is structured and how to extend it safely.

@tableofcontents

@section p51_layout_sec 1. Test Layout

Key directories:

- `tests/c/`: C test sources and shared fixture helpers
- `tests/smoke/`: executable smoke runner for tiny end-to-end flat/bent/brownian solve/post/restart/analytical workflows (single-rank and multi-rank, including particle restart branches)
- `tests/fixtures/`: Python/control-plane YAML fixtures

Current C files:

- `test_support.h`
- `test_support.c`
- `test_install_check.c`
- `test_geometry.c`
- `test_solver_kernels.c`
- `test_particle_kernels.c`
- `test_io.c`
- `test_logging.c`
- `test_postprocessing.c`
- `test_vtk_io.c`
- `test_postprocessor.c`
- `test_statistics.c`
- `test_grid.c`
- `test_metric.c`
- `test_boundaries.c`
- `test_poisson_rhs.c`
- `test_runtime_kernels.c`
- `test_mpi_kernels.c`

@subsection p51_case_matrix_ssec 1.1 Named C Test Case Matrix

Representative named cases by suite (exact strings used in `PicurvTestCase` arrays):

- `test_geometry.c`:
  - `scalar-interpolation-constant`
  - `vector-interpolation-constant`
  - `signed-distance-and-classification`
- `test_solver_kernels.c`:
  - `les-filter-paths`
  - `analytical-geometry-selection`
  - `driven-channel-flow-source`
- `test_particle_kernels.c`:
  - `check-cell-within-local-grid`
  - `initialize-traversal-parameters`
  - `retrieve-current-cell`
- `test_io.c`:
  - `should-write-data-output`
  - `verify-path-existence`
  - `write-and-read-simulation-fields`
  - `parse-post-processing-settings`
  - `trim-whitespace`
  - `bc-string-parsers`
  - `validate-bc-handler-for-type`
  - `parse-scaling-information`
- `test_logging.c`:
  - `get-log-level-from-environment`
  - `allowed-functions-filter`
  - `particle-console-snapshot-cadence`
- `test_postprocessing.c`:
  - `compute-specific-ke`
  - `compute-displacement`
  - `compute-nodal-average-scalar`
  - `normalize-relative-field`
  - `dimensionalize-pressure-field`
  - `compute-qcriterion-zero-flow`
- `test_vtk_io.c`:
  - `prepare-output-coordinates-subsamples-interior-grid`
  - `prepare-output-eulerian-field-data-subsamples-scalar`
  - `prepare-output-particle-data-subsampling`
  - `create-vtk-file-from-metadata-writes-structured-grid`
- `test_postprocessor.c`:
  - `setup-postprocess-swarm-registers-pipeline-fields`
  - `eulerian-data-processing-pipeline-runs-configured-kernels`
  - `particle-data-processing-pipeline-computes-specific-ke`
  - `global-statistics-pipeline-writes-msd-csv`
  - `write-eulerian-file-writes-vts`
  - `write-particle-file-writes-vtp`
- `test_statistics.c`:
  - `compute-particle-msd-writes-csv`
  - `compute-particle-msd-empty-swarm-no-output`
- `test_grid.c`:
  - `compute-local-bbox-uniform-grid`
  - `gather-and-broadcast-bboxes`
- `test_metric.c`:
  - `invert-covariant-metric-diagonal`
  - `metric-velocity-contravariant-identity`
  - `metric-velocity-contravariant-scaled-axes`
  - `face-normal-and-area-axis-aligned`
  - `characteristic-length-scale-axis-aligned`
- `test_boundaries.c`:
  - `can-rank-service-face-matches-inlet-when-defined`
  - `can-rank-service-inlet-face-requires-definition`
  - `boundary-condition-factory-assignments`
  - `boundary-condition-factory-implemented-handler-matrix`
  - `boundary-condition-factory-rejects-unsupported-handler`
  - `deterministic-face-grid-location-matrix`
  - `random-inlet-face-location-matrix`
- `test_poisson_rhs.c`:
  - `update-pressure-adds-phi`
  - `poisson-rhs-zero-divergence`
  - `compute-body-forces-dispatcher`
  - `compute-eulerian-diffusivity-molecular-only`
- `test_runtime_kernels.c`:
  - `distribute-particles-remainder-handling`
  - `is-particle-inside-bbox-basic-cases`
  - `update-particle-weights-computes-expected-ratios`
  - `update-particle-position-without-brownian-contribution`
  - `update-particle-field-iem-relaxation`
  - `set-initial-interior-field-ignores-non-ucont-request`
  - `set-initial-interior-field-constant-profile-on-z-inlet`
  - `wall-noslip-and-freeslip-helpers`
  - `wall-model-scalar-helpers`
  - `validate-driven-flow-configuration-no-driven-handlers`
  - `compute-smagorinsky-constant-constant-model`
  - `update-solver-history-vectors-shifts-states`
  - `get-owned-cell-range-single-rank-accounting`
  - `compute-and-store-neighbor-ranks-single-rank`
- `test_mpi_kernels.c`:
  - `distribute-particles-collective-consistency`
  - `bounding-box-collectives-multi-rank`

@section p51_targets_sec 2. Canonical Targets

Canonical user-facing targets:

- `make doctor`
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
- `make smoke`
- `make smoke-mpi`
- `make smoke-mpi-matrix`
- `make coverage-c`
- `make check`
- `make check-mpi`
- `make check-mpi-matrix`
- `make check-full`

Compatibility aliases (`install-check`, `ctest*`) should be preserved for continuity, but documentation and new contributor guidance should always prefer the canonical names.

@section p51_choosing_sec 3. Where A Test Belongs

Use this rule set:

- `doctor`:
  - only for installation/toolchain/PETSc viability checks
- `unit-*`:
  - for isolated functions, kernels, and small component seams
- `smoke`:
  - for real executable launch or tiny end-to-end entrypoint checks
- `test-python`:
  - for `scripts/picurv`, schema, and repository metadata behavior

If a test can run against a small local fixture without invoking the full executables, it almost always belongs in `unit-*`.

Aggregate gate guidance for C contributors:

- use `make check` for default local branch gating
- use `make check-mpi` or `make check-mpi-matrix` for targeted MPI validation
- use `make check-full` when you need one command that runs all MPI layers (`unit-mpi`, `smoke-mpi`, and `smoke-mpi-matrix`) after baseline `check`

@section p51_smoke_bridge_sec 4. Smoke Harness Relationship

`tests/smoke/run_smoke.sh` is not a C unit file, but it is the integration bridge that verifies how C executables behave under real launcher orchestration.

Smoke checks include:

- executable `-help` viability
- tiny solve/post end-to-end runs
- particle restart branches (`load` and `init`)
- restart-equivalence continuity check
- analytical Brownian branch
- MPI rank variants and rank-matrix sweeps

Environment knobs relevant to C contributors:

- `TEST_MPI_NPROCS` (for `unit-mpi`)
- `SMOKE_MPI_NPROCS` (for `smoke-mpi`)
- `SMOKE_MPI_MATRIX_NPROCS` (for `smoke-mpi-matrix`)
- `KEEP_SMOKE_TMP=1` (retain smoke workspace)

@section p51_fixtures_sec 5. Fixture Helpers

`tests/c/test_support.*` provides the shared PETSc-backed fixture layer.

Current responsibilities:

- create minimal `SimCtx` and `UserCtx` fixtures
- create deterministic `DMDA` objects
- create deterministic `DMSwarm` objects
- seed identity metric fields
- create unique temporary directories under `/tmp`
- provide assertion helpers for reals, ints, vectors, and file existence

Future C tests should reuse these helpers instead of rebuilding fixtures inline unless a suite truly needs a specialized setup.

@section p51_style_sec 6. Conventions

Naming:

- keep file names aligned with the target name (`test_io.c` -> `make unit-io`)
- use function names that state the behavior under test

Assertions:

- fail fast
- print the failing context plus expected and actual values

Tolerances:

- `1e-12` for algebraic invariants
- `1e-10` to `1e-9` for interpolation and I/O round-trips
- exact equality for integer/state checks

@section p51_fs_sec 7. Temporary Files

Rules:

- write temporary artifacts under `/tmp`
- never write repo-tracked test outputs
- keep paths easy to inspect (`/tmp/picurv-test-*`) so a failed test can be debugged manually

@section p51_extending_sec 8. Adding A New C Unit Suite

Recommended workflow:

1. add or reuse a `tests/c/test_<area>.c` file
2. build on `test_support.*`
3. add a dedicated executable in the `Makefile`
4. expose a canonical `unit-<area>` target
5. add a compatibility alias only if it matches an existing naming family
6. document the new target in the README and testing guide if it is user-facing

@section p51_refs_sec 9. Related Pages

- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 01_Installation**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **C Test Suite Developer Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
