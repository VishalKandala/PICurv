@page 13_Code_Architecture Code Architecture

@anchor _Code_Architecture

This page is the developer-oriented map of the current PICurv C codebase.

@tableofcontents

@section p13_entry_sec 1. Executable Entry Points

- Solver executable entry: `src/simulator.c`
- Postprocessor executable entry: `src/postprocessor.c`

Both rely on shared setup/context infrastructure from `setup.c`, `io.c`, and `variables.h`.

@section p13_solver_flow_sec 2. Solver Runtime Flow (simulator.c)

High-level stages:

1. `PetscInitialize`
2. `CreateSimulationContext` (`setup.c`) parses control/options and initializes defaults
3. `SetupSimulationEnvironment` configures run directories and environment-dependent logging setup
4. Setup stack:
   - `SetupGridAndSolvers`
   - `SetupBoundaryConditions`
   - `SetupDomainRankInfo`
   - `InitializeEulerianState`
   - `InitializeParticleSwarm` (if particles enabled)
5. Time integration via `AdvanceSimulation`
6. Finalization via profiling teardown + `FinalizeSimulation` + `PetscFinalize`

@section p13_post_flow_sec 3. Postprocessor Runtime Flow (postprocessor.c)

- Parse post recipe (`post.run`) into `PostProcessParams`
- Load requested Eulerian/particle fields by timestep
- Execute configured Eulerian/Lagrangian/statistics pipelines
- Write VTK outputs (`.vts`, `.vtp`) and statistics CSV outputs

@section p13_contexts_sec 4. Core Context Objects

@subsection p13_simctx_ssec 4.1 SimCtx

- Declared in `include/variables.h`
- Holds global run configuration and top-level handles
- Populated mainly in `CreateSimulationContext`

@subsection p13_userctx_ssec 4.2 UserCtx

- Block/grid-level state container (`DM`, vectors, metrics, block-local geometry)
- Used throughout grid, solver, BC, and post kernels
- Finest-level block arrays are central runtime work objects

@section p13_modules_sec 5. Module Responsibilities

- Setup/config ingestion:
  - `src/setup.c`
  - `src/io.c`
- Grid and metrics:
  - `src/grid.c`
  - `src/Metric.c`
- Momentum/pressure/rhs:
  - `src/momentumsolvers.c`
  - `src/poisson.c`
  - `src/rhs.c`
  - `src/solvers.c`
- BC system:
  - `src/Boundaries.c`
  - `src/BC_Handlers.c`
- Particles and coupling:
  - `src/ParticleSwarm.c`
  - `src/ParticleMotion.c`
  - `src/ParticlePhysics.c`
  - `src/interpolation.c`
  - `src/particle_statistics.c`
- Postprocessing and VTK:
  - `src/postprocessor.c`
  - `src/postprocessing_kernels.c`
  - `src/vtk_io.c`

@subsection p13_module_api_matrix_ssec 5.1 File-Level API Entry Points

For contributor orientation, the table below lists high-value public entry points per subsystem.
Function names come from `include/*.h` and represent the safest integration seams.

- startup/context:
  - files: `setup.c`, `simulator.c`
  - APIs: @ref CreateSimulationContext, @ref SetupSimulationEnvironment, @ref SetupGridAndSolvers, @ref FinalizeSimulation
- runtime loop:
  - files: `runloop.c`, `solvers.c`
  - APIs: @ref AdvanceSimulation, @ref FlowSolver, @ref UpdateSolverHistoryVectors
- pressure/projection/rhs:
  - files: `poisson.c`, `rhs.c`, `momentumsolvers.c`, `BodyForces.c`
  - APIs: @ref PoissonSolver_MG, @ref Projection, @ref ComputeRHS, @ref MomentumSolver_DualTime_Picard_RK4
- grid/metrics:
  - files: `grid.c`, `Metric.c`
  - APIs: @ref DefineAllGridDimensions, @ref InitializeAllGridDMs, @ref AssignAllGridCoordinates, @ref CalculateAllGridMetrics
- boundary system:
  - files: `Boundaries.c`, `BC_Handlers.c`, `wallfunction.c`
  - APIs: @ref BoundarySystem_Initialize, @ref BoundaryCondition_Create, @ref ApplyBoundaryConditions, @ref Validate_DrivenFlowConfiguration
- particle stack:
  - files: `ParticleSwarm.c`, `ParticleMotion.c`, `ParticlePhysics.c`, `walkingsearch.c`, `interpolation.c`
  - APIs: @ref InitializeParticleSwarm, @ref LocateAllParticlesInGrid, @ref PerformMigration, @ref UpdateAllParticleFields, @ref InterpolateAllFieldsToSwarm
- I/O and post:
  - files: `io.c`, `postprocessor.c`, `postprocessing_kernels.c`, `particle_statistics.c`, `vtk_io.c`
  - APIs: @ref ReadSimulationFields, @ref WriteSimulationFields, @ref ParsePostProcessingSettings, @ref EulerianDataProcessingPipeline, @ref GlobalStatisticsPipeline
- analytical/init:
  - files: `AnalyticalSolutions.c`, `initialcondition.c`
  - APIs: @ref AnalyticalSolutionEngine, @ref InitializeEulerianState, @ref SetInitialInteriorField

@subsection p13_source_to_tests_ssec 5.2 Source-to-Test Coverage Lens

Current tests cover all `src/*.c` files at least at module level (unit suites and/or smoke).
Coverage depth is intentionally uneven:

- orchestration-heavy files (`setup.c`, `io.c`, `runloop.c`, `postprocessor.c`) are mainly covered by smoke and Python orchestration tests.
- math-kernel-heavy files (`Metric.c`, `rhs.c`, `poisson.c`, post kernels) are covered by targeted C unit suites plus smoke.
- utility-heavy files (`walkingsearch.c`, `wallfunction.c`, `Filter.c`) are covered by specific unit cases but still benefit from deeper function-level API prose.

When adding docs, prioritize:

1. explicit call-sequence narratives for orchestration modules
2. contract/units/assumptions for numerical kernels
3. failure signals and diagnostics for runtime-heavy paths

@section p13_ingestion_sec 6. Configuration Ingestion Boundaries

Primary ingestion sites:
- `setup.c`: PETSc option parsing for solver/post shared runtime controls
- `io.c`: grid read/generation inputs, restart/data IO, post recipe parsing
- logging path includes environment variable ingress (`LOG_LEVEL`)

Not all option consumption is explicit `PetscOptionsGet*`; PETSc dynamic ingestion also occurs through calls like `KSPSetFromOptions` in `poisson.c`.

@section p13_extension_sec 7. Where to Extend

- New user-visible config key:
  1. add schema/template key
  2. validate/map in `scripts/picurv`
  3. parse in `setup.c` or `io.c`
  4. wire runtime consumer
  5. update docs pages 14/15/16
- New post kernel:
  - implement in `postprocessing_kernels.c`
  - expose task mapping in `picurv`
  - document in post reference
- New particle model:
  - extend `ParticlePhysics.c`/`ParticleMotion.c` interface
  - add structured schema and ingestion wiring

@section p13_next_steps_sec 8. Next Steps

- Config contract: **@subpage 14_Config_Contract**
- Ingestion map: **@subpage 15_Config_Ingestion_Map**
- Extension workflow: **@subpage 16_Config_Extension_Playbook**
- Workflow growth paths: **@subpage 17_Workflow_Extensibility**
- Numerical-method overviews: **@subpage 21_Methods_Overview**
- Repository map: **@subpage 30_Repository_Navigation**
- Momentum solver index: **@subpage 31_Momentum_Solvers**
- Analytical solution index: **@subpage 32_Analytical_Solutions**
- Initial-condition index: **@subpage 33_Initial_Conditions**
- Particle-model index: **@subpage 34_Particle_Model_Overview**
- API docs quality status: **@subpage 35_API_Documentation_Status**
- Low-priority fix queue: **@subpage 29_Maintenance_Backlog**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Code Architecture** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
