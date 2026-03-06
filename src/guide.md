# Source Guide

This directory contains the C implementation for the solver, postprocessor, and shared runtime subsystems. The source tree is organized by functional responsibilities rather than by executable boundary, so understanding module interaction is essential for safe changes.

## Runtime Entry Points

- `simulator.c`: main entry for solver execution.
- `postprocessor.c`: main entry for postprocessing execution.

## Module Map (File -> Responsibility -> Key Public APIs)

- startup/context:
  - files: `setup.c`, `simulator.c`
  - APIs: `CreateSimulationContext`, `SetupSimulationEnvironment`, `SetupGridAndSolvers`, `FinalizeSimulation`
- runtime loop:
  - files: `runloop.c`, `solvers.c`
  - APIs: `AdvanceSimulation`, `FlowSolver`, `UpdateSolverHistoryVectors`
- momentum/rhs/pressure:
  - files: `momentumsolvers.c`, `rhs.c`, `poisson.c`, `BodyForces.c`, `Filter.c`, `les.c`
  - APIs: `MomentumSolver_DualTime_Picard_RK4`, `ComputeRHS`, `PoissonSolver_MG`, `Projection`, `ComputeSmagorinskyConstant`
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
- analytical/initialization:
  - files: `AnalyticalSolutions.c`, `initialcondition.c`
  - APIs: `AnalyticalSolutionEngine`, `InitializeEulerianState`, `SetInitialInteriorField`

## How To Navigate During Development

1. Start from the top-level orchestrator (`runloop.c` and `solvers.c`) to identify call sequence.
2. Move into one subsystem module at a time.
3. Track shared state usage through `SimCtx`/`UserCtx` and supporting structs.
4. Confirm header contract alignment in `include/` before finalizing changes.

Practical tip:

- when tracing a YAML key into C, start at generated `runs/<run_id>/config/*.control` and then follow `setup.c`/`io.c` option ingestion before jumping into physics kernels.

## Development Notes

- Shared state is centralized in `SimCtx`/`UserCtx` from `variables.h`.
- Prefer extending existing module boundaries before creating new files.
- Keep high-level orchestration logic in orchestrator modules; keep math kernels in subsystem files.
- Avoid hidden behavior in utility helpers that bypass main execution flow assumptions.

## Change-Safety Checklist

- Update matching headers for every public symbol change.
- Add/adjust C tests when modifying numerics, BC handling, or data exchange paths.
- Re-run relevant `make unit-*` and smoke targets before merge.
- For control-plane touching changes, also run `make test-python` because `scripts/picurv` controls many runtime defaults and mode gates.

## Coverage Lens

The source tree is broadly touched by current tests, but function-level direct documentation remains thinner than module-level coverage.

- highest-ingestion modules (`setup.c`, `io.c`) are heavily documented and tested
- runtime-heavy modules (`runloop.c`, `solvers.c`, `poisson.c`, particle stack) are strongly smoke-covered
- utility-heavy modules (`wallfunction.c`, `walkingsearch.c`, `Filter.c`) rely more on focused unit tests and have less API prose density

When expanding docs, prefer module-level execution narratives plus targeted API notes for high-risk entry points rather than trying to document every internal helper in one pass.

## Related Docs

- `include/guide.md`
- https://vishalkandala.me/picurv-docs/13_Code_Architecture.html
- https://vishalkandala.me/picurv-docs/21_Methods_Overview.html
