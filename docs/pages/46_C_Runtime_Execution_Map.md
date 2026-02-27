@page 46_C_Runtime_Execution_Map C Runtime Execution Map

This page is a contributor-oriented map of how the C solver executes from process start to timestep loop.
Use it as a practical companion to architecture and methods pages when modifying solver behavior.

@tableofcontents

@section startup_sec 1. Solver Startup Order (`src/picsolver.c`)

High-level startup sequence:

1. `PetscInitialize`
2. @ref CreateSimulationContext
3. `SetupSimulationEnvironment`
4. `SetupGridAndSolvers`
5. `SetupBoundaryConditions`
6. `SetupDomainRankInfo`
7. @ref InitializeEulerianState
8. @ref InitializeParticleSwarm (if `np > 0`)
9. `DisplayBanner`
10. initial settlement/restart finalization
11. @ref AdvanceSimulation

@section ingestion_sec 2. Python-to-C Configuration Boundary

`scripts/pic.flow` is the control-plane generator.
It writes normalized runtime artifacts under `runs/<run_id>/config/` and launches C binaries with `-control_file`.

Core generated files consumed by C:

- `*.control` (solver flags),
- `bcs.run` (boundary face/type/handler + params),
- `whitelist.run` and `profile.run`,
- `grid.run` (for file/grid_gen paths),
- `post.run` (postprocessor path).

@section core_structs_sec 3. Core Runtime Structs

Most solver-wide state flows through:

- `SimCtx`: global run configuration, solver controls, pointers to hierarchy and shared runtime services.
- `UserCtx`: per-block/per-level field ownership, DM/Vec handles, boundary configs, and local coupling context.
- `BoundaryFaceConfig`: per-face mathematical type + handler + param list.
- swarm particle records and DMSwarm fields (`position`, `DMSwarm_CellID`, status fields, etc.).

@section init_branches_sec 4. Initialization Branches

Eulerian:

- @ref InitializeEulerianState selects fresh/load/analytical branch and then seeds history vectors.

Lagrangian:

- @ref InitializeParticleSwarm selects initialize-vs-load based on `StartStep` and `particle_restart_mode`.
- fresh path uses random/deterministic placement logic depending on `-pinit`.
- load path may use @ref MigrateRestartParticlesUsingCellID for direct ownership migration.

@section loop_sec 5. Timestep Loop (`AdvanceSimulation`)

Per-step sequence:

1. update `step` and `time`,
2. reset particle statuses (if particles enabled),
3. Eulerian step:
   - `FlowSolver` for solve mode, or load/analytical branch,
4. Lagrangian step:
   - position update,
   - location/migration,
   - lost-particle handling,
   - interpolation of new Eulerian fields to particles,
   - particle-field updates,
   - scatter particle fields back to Eulerian grid,
5. update history vectors,
6. write outputs on configured cadence.

@section boundaries_sec 6. Boundary System Runtime Hooks

Boundary lifecycle is object-style (function pointers per handler):

1. parse `bcs.run`,
2. create handler objects via factory,
3. run `Initialize` once,
4. on each step, run handler phases (`PreStep`, `Apply`, optional `PostStep`) in priority order.

Useful entry points:

- @ref BoundarySystem_Initialize
- @ref BoundaryCondition_Create
- @ref ApplyBoundaryConditions
- @ref DeterminePeriodicity

@section extension_sec 7. Safe Extension Workflow (C Side)

When adding or changing physics behavior:

1. identify owning module (`solvers.c`, `rhs.c`, `poisson.c`, `Particle*`, `Boundaries*`),
2. add/modify fields in `SimCtx` or `UserCtx` only when ownership is clear,
3. keep history-vector and ghost-update contracts intact,
4. update logging labels and diagnostics,
5. update Python ingestion path so YAML and C flags remain consistent.

@section debug_sec 8. Debugging Entry Points

High-value checks during development:

- `DisplayBanner` summary (BCs, modes, solver selection),
- per-step profile summaries (`ProfilingLogTimestepSummary`),
- particle metrics and location/migration counters,
- strict YAML validation through `pic.flow validate` before solver execution.

@section refs_sec 9. Related Pages

- **@subpage 13_Code_Architecture**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 21_Methods_Overview**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**
