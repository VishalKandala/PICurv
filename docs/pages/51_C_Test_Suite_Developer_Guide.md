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
- `test_setup_lifecycle.c`
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
- `test_periodic_dev.c`
- `test_poisson_rhs.c`
- `test_runtime_kernels.c`
- `test_mpi_kernels.c`

@subsection p51_case_matrix_ssec 1.1 Named C Test Case Matrix

Representative named cases by suite (exact strings used in `PicurvTestCase` arrays):

- `test_install_check.c`:
  - `environment-visible`
  - `basic-petsc-objects`
- `test_geometry.c`:
  - `scalar-interpolation-constant`
  - `vector-interpolation-constant`
  - `signed-distance-and-classification`
- `test_setup_lifecycle.c`:
  - `shared-runtime-fixture-contracts`
  - `setup-lifecycle-core-solver-setup`
  - `setup-lifecycle-particle-initialization`
  - `setup-lifecycle-random-generators-and-cleanup`
  - `setup-lifecycle-cleanup-across-initialization-states`
- `test_solver_kernels.c`:
  - `les-filter-paths`
  - `analytical-geometry-selection`
  - `analytical-solution-engine-dispatch`
  - `analytical-particle-dispatch`
  - `compute-eddy-viscosity-les-deterministic-field`
  - `analytical-solution-engine-taylor-green-samples`
  - `flow-solver-rejects-unsupported-momentum-solver-type`
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
  - `display-banner-startup-summary`
- `test_logging.c`:
  - `string-conversion-helpers`
  - `get-log-level-from-environment`
  - `allowed-functions-filter`
  - `particle-console-snapshot-cadence`
  - `logging-file-parsing-and-formatting-helpers`
  - `logging-continuity-and-field-diagnostics`
  - `interpolation-error-logging`
  - `particle-field-table-logging`
  - `particle-console-snapshot-logging`
  - `particle-metrics-logging`
  - `field-anatomy-logging`
  - `profiling-lifecycle-helpers`
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
  - `inlet-constant-velocity-handler-behavior`
  - `inlet-parabolic-profile-handler-behavior`
  - `outlet-conservation-handler-behavior`
- `test_periodic_dev.c`:
  - `periodic-geometric-factory-assignment`
  - `transfer-periodic-face-field-copies-x-faces`
  - `apply-metrics-periodic-bcs-copies-aj-faces`
  - `periodic-driven-constant-handler-behavior`
  - `periodic-driven-constant-rejects-non-periodic-face`
- `test_poisson_rhs.c`:
  - `update-pressure-adds-phi`
  - `poisson-rhs-zero-divergence`
  - `compute-body-forces-dispatcher`
  - `compute-eulerian-diffusivity-molecular-only`
  - `convection-zero-field`
  - `viscous-uniform-field`
  - `compute-rhs-zero-field-no-forcing`
  - `compute-eulerian-diffusivity-gradient-constant-field`
  - `poisson-null-space-function-removes-mean`
  - `poisson-lhsnew-assembles-operator`
  - `projection-zero-phi-leaves-velocity-unchanged`
  - `projection-linear-phi-corrects-velocity`
- `test_runtime_kernels.c`:
  - `distribute-particles-remainder-handling`
  - `is-particle-inside-bbox-basic-cases`
  - `update-particle-weights-computes-expected-ratios`
  - `update-particle-position-without-brownian-contribution`
  - `update-particle-field-iem-relaxation`
  - `set-initial-interior-field-ignores-non-ucont-request`
  - `set-initial-interior-field-constant-profile-on-z-inlet`
  - `interpolate-all-fields-to-swarm-constant-fields`
  - `scatter-all-particle-fields-to-euler-fields-averages-psi`
  - `calculate-particle-count-per-cell-counts-global-cell-ids`
  - `reset-all-particle-statuses-leaves-lost-particles-untouched`
  - `check-and-remove-out-of-bounds-particles-removes-escaped-particle`
  - `check-and-remove-lost-particles-removes-lost-entries`
  - `calculate-brownian-displacement-deterministic-seed`
  - `update-all-particle-positions-moves-swarm-entries`
  - `locate-all-particles-in-grid-prior-cell-fast-path`
  - `locate-all-particles-in-grid-guess-path-resolves-local-particle`
  - `wall-noslip-and-freeslip-helpers`
  - `wall-model-scalar-helpers`
  - `wall-model-velocity-helpers`
  - `wall-function-vector-wrappers`
  - `validate-driven-flow-configuration-no-driven-handlers`
  - `compute-smagorinsky-constant-constant-model`
  - `update-solver-history-vectors-shifts-states`
  - `get-owned-cell-range-single-rank-accounting`
  - `compute-and-store-neighbor-ranks-single-rank`
  - `runtime-walltime-guard-parses-positive-seconds`
  - `runtime-walltime-guard-estimator-helpers`
  - `runtime-walltime-guard-trigger-decision`
- `test_mpi_kernels.c`:
  - `distribute-particles-collective-consistency`
  - `bounding-box-collectives-multi-rank`
  - `restart-cellid-migration-moves-particle-to-owner`

@section p51_targets_sec 2. Canonical Targets

Canonical user-facing targets:

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
- `make coverage-c`
- `make check`
- `make check-mpi`
- `make check-mpi-matrix`
- `make check-full`
- `make check-stress`

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
  - banner presence is authoritative; local PETSc may exit with code `62` (`PETSC_ERR_ARG_WRONG`) after printing help
- tiny solve/post end-to-end runs
- particle restart branches (`load` and `init`)
- restart-equivalence continuity check
- analytical Brownian branch
- MPI rank variants and rank-matrix sweeps
- opt-in stress extensions (`make smoke-stress`) for longer particle cycling, chained restarts, parabolic-inlet runtime coverage, periodic constant-flux validate/dry-run coverage, and extra-rank MPI particle runs
- periodic runtime development harness (`make smoke-periodic-dev`) for the in-progress periodic constant-flux path; not part of the default gate

Environment knobs relevant to C contributors:

- `TEST_MPI_NPROCS` (for `unit-mpi`)
- `SMOKE_MPI_NPROCS` (for `smoke-mpi`)
- `SMOKE_MPI_MATRIX_NPROCS` (for `smoke-mpi-matrix`)
- `KEEP_SMOKE_TMP=1` (retain smoke workspace)

@section p51_fixtures_sec 5. Fixture Helpers

`tests/c/test_support.*` provides the shared PETSc-backed fixture layer.

Current responsibilities:

- create minimal `SimCtx` and `UserCtx` fixtures for fast kernel-level tests
- create a richer tiny-runtime fixture through the real setup path for behavior/orchestrator tests
- create deterministic `DMDA` objects using the production node contract (`da` sized as `IM+1/JM+1/KM+1`)
- derive `fda` from the `da` coordinate-DM path and keep `user->info` aligned with `da`
- create deterministic `DMSwarm` objects with the production solver field set (`position`, `velocity`, `DMSwarm_CellID`, `weight`, `Diffusivity`, `DiffusivityGradient`, `Psi`, `DMSwarm_location_status`)
- seed identity metric fields
- normalize DM coordinates to the same `0..1` convention used by the main grid path
- create unique temporary directories under `/tmp`
- provide assertion helpers for reals, ints, vectors, and file existence

Future C tests should reuse these helpers instead of rebuilding fixtures inline unless a suite truly needs a specialized setup. Periodic boundary harnesses belong in the dedicated development targets until the product implementation is stable enough for the default gate.

@subsection p51_gap_ssec 5.1 Current next-gap priorities

The current highest-value additions are:

- direct walking-search branch tests for `LocateParticleOrFindMigrationTarget`
  - boundary clamp
  - ghost-region handoff
  - tie-breaker
  - `LOST` and `MIGRATING_OUT` outcomes
- explicit direction-complete and failure-path coverage for the `GuessParticleOwnerWithBBox` heuristic
- non-restart MPI particle migration tests covering multi-pass handoff, newcomer flagging, and count conservation
- direct positive-path momentum harnesses for `MomentumSolver_Explicit_RungeKutta4` and one small invariant case for `MomentumSolver_DualTime_Picard_RK4`
- deeper `PoissonSolver_MG` and periodic/IBM stencil behavior checks beyond the current helper-level `unit-poisson-rhs` surface
- richer-runtime fixture variants for more geometry/topology cases than the tiny Cartesian baseline

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
