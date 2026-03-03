@page 13_Code_Architecture Code Architecture

This page is the developer-oriented map of the current PICurv C codebase.

@tableofcontents

@section entry_sec 1. Executable Entry Points

- Solver executable entry: `src/simulator.c`
- Postprocessor executable entry: `src/postprocessor.c`

Both rely on shared setup/context infrastructure from `setup.c`, `io.c`, and `variables.h`.

@section solver_flow_sec 2. Solver Runtime Flow (simulator.c)

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

@section post_flow_sec 3. Postprocessor Runtime Flow (postprocessor.c)

- Parse post recipe (`post.run`) into `PostProcessParams`
- Load requested Eulerian/particle fields by timestep
- Execute configured Eulerian/Lagrangian/statistics pipelines
- Write VTK outputs (`.vts`, `.vtp`) and statistics CSV outputs

@section contexts_sec 4. Core Context Objects

@subsection simctx_ssec 4.1 SimCtx

- Declared in `include/variables.h`
- Holds global run configuration and top-level handles
- Populated mainly in `CreateSimulationContext`

@subsection userctx_ssec 4.2 UserCtx

- Block/grid-level state container (`DM`, vectors, metrics, block-local geometry)
- Used throughout grid, solver, BC, and post kernels
- Finest-level block arrays are central runtime work objects

@section modules_sec 5. Module Responsibilities

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

@section ingestion_sec 6. Configuration Ingestion Boundaries

Primary ingestion sites:
- `setup.c`: PETSc option parsing for solver/post shared runtime controls
- `io.c`: grid read/generation inputs, restart/data IO, post recipe parsing
- logging path includes environment variable ingress (`LOG_LEVEL`)

Not all option consumption is explicit `PetscOptionsGet*`; PETSc dynamic ingestion also occurs through calls like `KSPSetFromOptions` in `poisson.c`.

@section extension_sec 7. Where to Extend

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

@section next_steps_sec 8. Next Steps

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
