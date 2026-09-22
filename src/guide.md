# Source Guide

This directory contains the C implementation for the solver, postprocessor, and shared runtime subsystems. The source tree is organized by functional responsibilities rather than by executable boundary, so understanding module interaction is essential for safe changes.

## Runtime Entry Points

- `simulator.c`: main entry for solver execution.
- `postprocessor.c`: main entry for postprocessing execution.

## Module Map (File -> Responsibility -> Key Public APIs)

- startup/context:
  - files: `setup.c`, `field_catalog.c`, `particle_field_catalog.c`, `simulator.c`
  - APIs: `CreateSimulationContext`, `SetupSimulationEnvironment`, `SetupGridAndSolvers`, `FieldGetView`, `FinalizeSimulation`
- runtime loop:
  - files: `runloop.c`, `solvers.c`
  - APIs: `AdvanceSimulation`, `FlowSolver`, `UpdateSolverHistoryVectors`
- momentum/rhs/pressure:
  - files: `momentumsolvers.c`, `momentum_newton_krylov.c`, `rhs.c`, `poisson.c`, `BodyForces.c`, `Filter.c`, `les.c`
  - APIs: `MomentumSolver_DualTime_Picard_JamesonRK`, `MomentumSolver_NewtonKrylov`, `ComputeRHS`, `PoissonSolver_MG`, `Projection`, `ComputeSmagorinskyConstant`
- grid/metrics:
  - files: `grid.c`, `Metric.c`
  - APIs: `DefineAllGridDimensions`, `InitializeAllGridDMs`, `AssignAllGridCoordinates`, `CalculateAllGridMetrics`
- boundary system:
  - files: `Boundaries.c`, `BC_Handlers.c`, `wallfunction.c`
  - APIs: `BoundarySystem_Initialize`, `BoundaryCondition_Create`, `ApplyBoundaryConditions`, `Validate_DrivenFlowConfiguration`
- particle transport/coupling:
  - files: `ParticleSwarm.c`, `ParticleMotion.c`, `ParticlePhysics.c`, `walkingsearch.c`, `interpolation.c`
  - APIs: `InitializeParticleSwarm`, `LocateAllParticlesInGrid`, `PerformMigration`, `UpdateAllParticleFields`, `InterpolateAllFieldsToSwarm`
- I/O and post:
  - files: `io.c`, `postprocessor.c`, `postprocessing_kernels.c`, `particle_statistics.c`, `vtk_io.c`
  - APIs: `ReadSimulationFields`, `WriteSimulationFields`, `ParsePostProcessingSettings`, `EulerianDataProcessingPipeline`, `GlobalStatisticsPipeline`
- field statistics:
  - files: `statistics_moments.c`, `statistics_window.c`, `statistics_accumulator.c`, `statistics_target.c`, `statistics_config.c`
  - APIs: `ParseFieldStatisticsConfig`, `FieldStatisticsIsActive`, `FieldStatisticsUpdateWindows`, `PicurvWindowOfferState`, `PicurvWindowAccumulate`, `FieldStatisticsPipeline`
- observability:
  - files: `logging.c`
  - APIs: `LOG_ALLOW`, `get_log_level`, `is_function_allowed`, `PicurvOpenDiagnosticsCsv`, `LOG_CONTINUITY_METRICS`, `LOG_SEARCH_METRICS`, `EmitStatisticsConsoleSnapshot`
- analytical/initialization:
  - files: `AnalyticalSolutions.c`, `initialcondition.c`
  - APIs: `AnalyticalSolutionEngine`, `InitializeEulerianState`, `SetInitialInteriorField`

## How To Navigate During Development

1. Start from the top-level orchestrator (`runloop.c` and `solvers.c`) to identify call sequence.
2. Move into one subsystem module at a time.
3. Track shared state usage through `SimCtx`/`UserCtx` and supporting structs.
4. Confirm header contract alignment in `include/` before finalizing changes.

Practical tip:

- when tracing a YAML key into C, start at generated `<run.config>/*.control` and then follow `setup.c`/`io.c` option ingestion before jumping into physics kernels.

## Development Notes

- Shared state is centralized in `SimCtx`/`UserCtx` from `variables.h`.
- Persistent Eulerian field identity and DM/Vec binding metadata is centralized
  in `field_catalog.c`; storage allocation remains explicit in `setup.c`.
- Persistent solver-particle identity and DMSwarm metadata is centralized in
  `particle_field_catalog.c`; dynamic postprocessor fields remain recipe-owned.
- Prefer extending existing module boundaries before creating new files.
- Keep high-level orchestration logic in orchestrator modules; keep math kernels in subsystem files.
- Avoid hidden behavior in utility helpers that bypass main execution flow assumptions.

## Removed Legacy Runtime Flags

Internal record; not published. On 2026-09-18 these PETSc options were deleted from
`CreateSimulationContext()` in `setup.c`, together with the `SimCtx` fields they wrote
(`include/variables.h`) and their entries in `tests/tooling/audit_ingress_manifest.json`.
Each was read into a field that nothing outside `setup.c` ever used - the only other
references were default assignments or code that was already commented out - so no run
could depend on them. The conductor never emitted any of them; they were reachable only
through a PETSc passthrough, where they now surface as PETSc "unused option" warnings.

| Removed option | Former `SimCtx` field | What it was |
| --- | --- | --- |
| `-sediment`, `-rheology` | `sediment`, `rheology` | legacy physics switches, never implemented here |
| `-thin`, `-blk` | `thin`, `blank` | legacy geometry switches |
| `-dgf_x/-dgf_y/-dgf_z/-dgf_ax/-dgf_ay/-dgf_az` | same names | legacy body degree-of-freedom switches |
| `-cop`, `-fish`, `-cstart`, `-fishcyl`, `-eel`, `-pizza`, `-turbine`, `-wing`, `-hydro`, `-Pipe` | `cop`, `fish`, `fish_c`, `fishcyl`, `eel`, `pizza`, `turbine`, `wing`, `hydro`, `Pipe` | legacy case-specific switches |
| `-Turbulent_Channel_z` | `channelz` | legacy case-specific switch |
| `-mg_max_it`, `-mg_idx` | `mg_MAX_IT`, `mg_idx` | legacy multigrid controls; the Poisson solve reads `-ps_*` options |
| `-str` | `STRONG_COUPLING` | legacy FSI coupling switch |
| `-grid1d`, `-Ogrid` | `grid1d`, `Ogrid` | legacy grid switches; the O-grid force branch was commented out |
| `-pbc_domain` | `blkpbc` | legacy periodic-block switch |
| `-grid_rotation_angle`, `-Croty`, `-Crotz` | same names | legacy grid-rotation parameters |
| `-U_bc` | `U_bc` | legacy boundary velocity |
| `-read_fields` | `readFields` | legacy field-read switch; `-euler_field_source load` owns this |
| `-rs_fsi`, `-duplicate` | `rstart_fsi`, `duplicate` | legacy FSI restart switches |
| `-no_of_bodies` | read removed; `NumberOfBodies` kept at 1 | the dormant immersed-body flux routine in `poisson.c` still reads the field |
| `-poisson_tol` | `poisson_tol` | a Poisson tolerance nothing read; removed with the `poisson_solver.tolerance` YAML key, which validation now refuses in favour of `absolute_tolerance`/`relative_tolerance`, and the conductor no longer emits it beside `-ps_ksp_atol` |

Recover any of them from history (the parent of the commit that removes this table's
source lines) if an implementation ever needs one; do not re-add a read without a
consumer.

## Removed RANS Hooks

Internal record; not published. On 2026-09-22 the k-omega RANS closure was removed and
the `turbulence.rans` subsystem returned to `planned`. Nothing was ever built: setup
never allocated the fields, and `FlowSolver`'s transport update was commented out, so a
case that enabled it aborted on a null vector at the end of step one. What went with it:

| Removed | Where it lived | What it was |
| --- | --- | --- |
| `-rans` option and `SimCtx::rans` | `setup.c`, `variables.h` | the selector; `-rans` is now refused as planned, not implemented |
| `K_Omega`, `lK_Omega`, `K_Omega_o`, `lK_Omega_o` | `variables.h`, `field_catalog.c` | the k-omega vectors and their two catalog entries; never allocated |
| `FIELD_AVAILABILITY_RANS`, `FIELD_DM_FDA2`, `UserCtx::fda2` | `field_catalog.h`, `variables.h`, `setup.c` | the availability flag and the two-component DM only those fields used; `fda2` was destroyed but never created |
| the commented k-omega update | `solvers.c` | calls to `K_Omega_Set_Constant`, `K_Omega_IC` and `Solve_K_Omega`, none of which exist |
| `simCtx->rans` guards | `rhs.c`, `runloop.c`, `io.c`, `momentumsolvers.c`, `momentum_newton_krylov.c` | eddy-viscosity, history-copy and Newton-Krylov scope branches that could never be taken |
| `-checkpoint_rans` | `io.c` | checkpoint metadata key, written and required; a checkpoint written now is not readable by a pre-removal binary, which still requires it. The format version stays 1, so existing checkpoints remain readable |
| `models.physics.turbulence.rans`, `normalize_rans_model()` | `picurv_cli/core.py` | the YAML block and its translation; the block is now refused with the planned-status message |

The wall-function treatment is unaffected: it is configured independently of any
closure, and only its two RANS-specific pairing refusals (`cabot` and `werner` under
RANS) went with the selector.

A future implementation starts from the charter on page 57, not from this history.

## Change-Safety Checklist

- Update matching headers for every public symbol change.
- Add/adjust C tests when modifying numerics, BC handling, or data exchange paths.
- Re-run relevant `make unit-*` and smoke targets before merge.
- For control-plane touching changes, also run `make test-python` because `picurv_cli/core.py` controls many runtime defaults and mode gates.

## Coverage Lens

The source tree is broadly touched by current tests, but function-level direct documentation remains thinner than module-level coverage.

- highest-ingestion modules (`setup.c`, `io.c`) are heavily documented and tested
- runtime-heavy modules (`runloop.c`, `solvers.c`, `poisson.c`, particle stack) are strongly smoke-covered
- utility-heavy modules (`wallfunction.c`, `walkingsearch.c`, `Filter.c`) rely more on focused unit tests and have less API prose density

Current next-gap priorities:

- direct walking-search branch coverage for locate/migrate edge cases
- unit-level single-step momentum harnesses: Explicit RK4 and Picard are gated end to end by `make smoke` (RK4's stable step and stability-limit stop) and measured for order, but no unit test drives one step on a fixture
- periodic Poisson/multigrid stencil branches (`PoissonSolver_MG` itself is covered by `make unit-poisson-rhs`)
- non-restart MPI migration and multi-pass particle handoff coverage
- richer-runtime fixture variants beyond the tiny Cartesian baseline

When expanding docs, prefer module-level execution narratives plus targeted API notes for high-risk entry points rather than trying to document every internal helper in one pass.

## Related Docs

- `include/guide.md`
- https://vishalkandala.me/docs/picurv/13_Code_Architecture.html
- https://vishalkandala.me/docs/picurv/21_Methods_Overview.html
