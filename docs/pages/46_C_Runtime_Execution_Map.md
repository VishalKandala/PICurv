@page 46_C_Runtime_Execution_Map C Runtime Execution Map

@anchor _C_Runtime_Execution_Map

This page is a contributor-oriented map of how the C solver executes from process start to timestep loop.
Use it as a practical companion to architecture and methods pages when modifying solver behavior.

@tableofcontents

@section p46_startup_sec 1. Solver Startup Order (`src/simulator.c`)

High-level startup sequence:

1. `PetscInitialize`
2. @ref CreateSimulationContext
3. @ref SetupSimulationEnvironment
4. @ref SetupGridAndSolvers
5. @ref SetupBoundaryConditions
6. @ref SetupDomainRankInfo
7. @ref InitializeEulerianState
8. @ref InitializeParticleSwarm (if `np > 0`)
9. @ref DisplayBanner
10. initial settlement/restart finalization
11. @ref AdvanceSimulation

Startup branch details:

- if `StartStep == 0`:
  - particles enabled: @ref PerformInitializedParticleSetup
  - particles disabled: initial Eulerian field write path
- if `StartStep > 0`:
  - @ref FinalizeRestartState dispatches by particle restart mode:
    - `load` -> @ref PerformLoadedParticleSetup
    - `init` -> @ref PerformInitializedParticleSetup

@section p46_ingestion_sec 2. Python-to-C Configuration Boundary

`picurv_cli/core.py` is the control-plane generator.
It writes normalized runtime artifacts under `<run.config>/` and launches C binaries with `-control_file`.

Core generated files consumed by C:

- `*.control` (solver flags),
- `bcs.run` (boundary face/type/handler + params),
- `whitelist.run` and `profile.run`,
- `grid.run` (for file/grid_gen paths),
- `post.run` (postprocessor path).

Physical-solution monitoring is normalized by `picurv` and written directly
into `*.control` as `-solution_convergence_*` options. There is no separate
observation-plan sidecar. Scientific field-statistics ingress will only be
added when the runtime, checkpoint, and postprocessor implementation is ready
as one usable feature.

@section p46_core_structs_sec 3. Core Runtime Structs

Most solver-wide state flows through:

- `SimCtx`: global run configuration, solver controls, pointers to hierarchy and shared runtime services.
- `UserCtx`: per-block/per-level field ownership, DM/Vec handles, boundary configs, and local coupling context.
- `BoundaryFaceConfig`: per-face mathematical type + handler + param list.
- swarm particle records and DMSwarm fields (`position`, `DMSwarm_CellID`, status fields, etc.).

Related ownership note:

- most numerics run at finest level `simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user[*]`
- multiblock loops (`for bi in block_number`) are common in setup, IO, and update stages

@section p46_init_branches_sec 4. Initialization Branches

Eulerian:

- @ref InitializeEulerianState selects fresh/load/analytical branch and then seeds history vectors.

Lagrangian:

- @ref InitializeParticleSwarm selects initialize-vs-load based on `StartStep` and `particle_restart_mode`.
- fresh path uses random/deterministic placement logic depending on `-pinit`.
- load path may use @ref MigrateRestartParticlesUsingCellID for direct ownership migration.

@section p46_loop_sec 5. Timestep Loop (`AdvanceSimulation`)

Per-step sequence:

1. update `step` and `time`,
2. reset particle statuses (if particles enabled),
3. Eulerian step:
   - @ref FlowSolver for solve mode, or load/analytical branch,
4. Lagrangian step:
   - @ref UpdateAllParticlePositions,
   - @ref LocateAllParticlesInGrid and migration,
   - @ref CheckAndRemoveLostParticles,
   - @ref InterpolateAllFieldsToSwarm,
   - @ref UpdateAllParticleFields,
   - @ref CalculateParticleCountPerCell and @ref ScatterAllParticleFieldsToEulerFields,
5. write the configured physical-solution convergence observation through the
   existing logger (unless disabled),
6. update history vectors,
7. commit a checkpoint bundle on configured cadence.

The loop calls @ref WriteCheckpointBundle. That thin coordinator delegates
payload I/O to:

- @ref WriteSimulationFields
- @ref WriteAllSwarmFields

@ref ProfilingLogTimestepSummary and all other persistent logging remain
separate loop services and are not checkpoint payloads.

Loop-time branch notes:

- Eulerian source mode is selected once per step (`solve`, `load`, or `analytical`), then particle coupling follows.
- particle migration is iterative with global settlement passes and explicit lost-particle handling.
- periodic particle console snapshots are controlled by monitor/profiling settings and rank-aware logging helpers.

@section p46_boundaries_sec 6. Boundary System Runtime Hooks

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

@subsection p46_residual_purity_sec 6a. Residual-Purity Invariant (Matrix-Free Solvers)

Any field read by the momentum residual or its boundary handlers must be
reconstructed or synchronized from the current SNES trial vector **before** it can
influence `F(X)`. A matrix-free residual (@ref MomentumNewtonKrylov_FormResidual)
that reads stale state is not a deterministic function of `X`, which breaks the
finite-difference Jacobian action and the line search.

Concretely, the Cartesian velocity state must be seeded in dependency order at the
top of each residual evaluation, before @ref ApplyBoundaryConditions runs (the
conservation-outlet handler reads `lUcat` on its first pass):

```
X / Ucont -> lUcont -> Ucat -> periodic Ucat finalization -> lUcat -> boundary application
```

`Contra2Cart()` alone is insufficient (it does not refresh `lUcat`), and the
reconstruction inside `ApplyBoundaryConditions()` runs after the first handler
sweep, so it cannot prepare that sweep. When adding a new field read by the
residual or a boundary handler, seed it from `X` here or expect residual-purity
regressions. See @ref 55_Newton_Krylov_Momentum_Solver, Section 5.

@section p46_extension_sec 7. Safe Extension Workflow (C Side)

When adding or changing physics behavior:

1. identify owning module (`solvers.c`, `rhs.c`, `poisson.c`, `Particle*`, `Boundaries*`),
2. add/modify fields in `SimCtx` or `UserCtx` only when ownership is clear,
3. keep history-vector and ghost-update contracts intact,
4. update logging labels and diagnostics,
5. update Python ingestion path so YAML and C flags remain consistent.

Cross-layer reminder:

- changes in `setup.c` / `io.c` option ingestion should also update:
  - `tests/tooling/audit_ingress_manifest.json`
  - `docs/pages/15_Config_Ingestion_Map.md`
  - relevant Python regression tests (`test_config_regressions.py` / `test_cli_smoke.py`)

@section p46_debug_sec 8. Debugging Entry Points

High-value checks during development:

- `DisplayBanner` summary (BCs, modes, solver selection),
- per-step profile summaries (`ProfilingLogTimestepSummary`),
- particle metrics and location/migration counters,
- strict YAML validation through `picurv validate` before solver execution.

Additional fast diagnostics:

- compare `manifest.json` + generated `<run.config>/*.control` against expected YAML role settings
- run `picurv run --dry-run --format json` to verify launch-mode/rank assumptions before expensive runs
- for restart issues, verify `--restart-from` / `--continue` resolution and the generated `-restart_dir` in control artifacts

@section p46_refs_sec 9. Related Pages

- **@subpage 13_Code_Architecture**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 21_Methods_Overview**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 56_Field_Identity_and_Layout_Catalog**
- **@subpage 58_Field_Statistics** (per-step accumulation hook and checkpoint state)
