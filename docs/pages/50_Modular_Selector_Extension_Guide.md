@page 50_Modular_Selector_Extension_Guide Modular Selector Extension Guide

@anchor _Modular_Selector_Extension_Guide

This page is the contributor checklist for adding new selector-driven options
without drifting the YAML -> launcher -> C runtime contract.

@tableofcontents

@section p50_general_sec 1. Standard Checklist For Any Selector

For every new selector value:

1. define the canonical user-facing value in YAML docs/templates,
2. add or update Python normalization/validation in `scripts/picurv`,
3. emit the correct generated control/post artifact mapping,
4. add any C enum/storage and parser wiring,
5. add the runtime dispatch/consumer branch,
6. update tests,
7. update the relevant reference pages.

Use the canonical value only. Do not add placeholder enum values or compatibility aliases for unimplemented options.

@section p50_momentum_sec 2. Momentum Solver Selector

- Schema home: `solver.yml -> strategy.momentum_solver`
- Canonical values:
  - `Explicit RK4`
  - `Dual Time Picard RK4`
- Python hook:
  - `normalize_momentum_solver_type()` in `scripts/picurv`
- Generated mapping:
  - `strategy.momentum_solver` -> `-mom_solver_type`
- C enum/storage:
  - `MomentumSolverType` in `include/variables.h`
- C parser:
  - `-mom_solver_type` in @ref CreateSimulationContext (`src/setup.c`)
- Runtime dispatch:
  - @ref FlowSolver in `src/solvers.c`
- String/log display:
  - `MomentumSolverTypeToString()` in `src/logging.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 08_Solver_Reference**
  - **@subpage 24_Dual_Time_Picard_RK4**
  - **@subpage 31_Momentum_Solvers**

@section p50_bc_sec 3. Boundary Condition Type Or Handler

- Schema home: `case.yml -> boundary_conditions`
- Canonical values:
  - see the supported type/handler matrix in **@subpage 44_Boundary_Conditions_Guide**
- Python hooks:
  - `BC_TYPE_MAP`
  - `BC_HANDLER_SPECS`
  - `validate_and_prepare_boundary_conditions()` in `scripts/picurv`
- Generated artifact:
  - `bcs.run`
- C parse/load:
  - @ref ParseAllBoundaryConditions in `src/io.c`
- C validation/factory:
  - `ValidateBCHandlerForBCType`
  - @ref BoundaryCondition_Create
- Runtime application:
  - `src/Boundaries.c`
  - `src/BC_Handlers.c`
- Canonical runtime representation:
  - `UserCtx.boundary_faces`
- Compatibility policy:
  - add any temporary adapter only when a concrete consumer exists; do not
    maintain parallel legacy BC arrays in `UserCtx`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 07_Case_Reference**
  - **@subpage 44_Boundary_Conditions_Guide**

@section p50_particle_init_sec 4. Particle Initialization Mode

- Schema home: `case.yml -> models.physics.particles.init_mode`
- Canonical values:
  - `Surface`
  - `Volume`
  - `PointSource`
  - `SurfaceEdges`
- Python hook:
  - `normalize_particle_init_mode()` in `scripts/picurv`
- Generated mapping:
  - `init_mode` -> `-pinit`
- C enum/storage:
  - `ParticleInitializationType` in `include/variables.h`
- C parser:
  - `-pinit` in @ref CreateSimulationContext
- Runtime:
  - `src/ParticleSwarm.c`
  - `src/ParticleMotion.c`
  - `src/runloop.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 45_Particle_Initialization_and_Restart**

@section p50_interp_sec 5. Interpolation Method

- Schema home: `solver.yml -> interpolation.method`
- Canonical values:
  - `Trilinear`
  - `CornerAveraged`
- Python hook:
  - `normalize_interpolation_method()` in `scripts/picurv`
- Generated mapping:
  - `interpolation_method` -> `-interpolation_method`
- C enum/storage:
  - `InterpolationMethod` in `include/variables.h`
- C parser:
  - `-interpolation_method` in @ref CreateSimulationContext (`src/setup.c`)
- Runtime dispatch:
  - @ref InterpolateEulerFieldToSwarm in `src/interpolation.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 27_Trilinear_Interpolation_and_Projection**
  - **@subpage 45_Particle_Initialization_and_Restart**

@section p50_field_init_sec 6. Field Initialization Mode

- Schema home: `case.yml -> properties.initial_conditions.mode`
- Canonical values:
  - `Zero`
  - `Constant`
  - `Poiseuille`
- Python hooks:
  - `normalize_field_init_mode()`
  - `resolve_initial_velocity_components()`
- Generated mapping:
  - `mode` -> `-finit`
  - `u_physical/v_physical/w_physical` -> `-ucont_x/-ucont_y/-ucont_z`
- C storage/parser:
  - `FieldInitialization` in `SimCtx`
  - `-finit` in @ref CreateSimulationContext
- Runtime:
  - `src/initialcondition.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 33_Initial_Conditions**

@section p50_analytical_sec 7. Analytical Type

- Schema home: `solver.yml -> operation_mode.analytical_type`
- Canonical values:
  - `TGV3D`
  - `ZERO_FLOW`
- Python hook:
  - `normalize_analytical_type()` in `scripts/picurv`
- Generated mapping:
  - `analytical_type` -> `-analytical_type`
- C storage/parser:
  - `AnalyticalSolutionType` in `SimCtx`
  - `-analytical_type` in @ref CreateSimulationContext
- Runtime:
  - `src/AnalyticalSolutions.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 08_Solver_Reference**
  - **@subpage 32_Analytical_Solutions**

@section p50_grid_sec 8. Grid Selector / Generator Selector

- Schema homes:
  - `case.yml -> grid.mode`
  - `case.yml -> grid.generator.grid_type`
- Canonical values:
  - `programmatic_c`
  - `file`
  - `grid_gen`
- Python orchestration:
  - `scripts/picurv` grid validation and run-dir staging
- Generator implementation:
  - `scripts/grid.gen`
- Generated mapping:
  - `programmatic_c` emits `-im/-jm/-km/...`
  - `file` emits `-grid_file`
  - `grid_gen` shells out to `grid.gen`
- C consumption:
  - `src/grid.c`
  - `src/io.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 07_Case_Reference**
  - **@subpage 48_Grid_Generator_Guide**

@section p50_profiling_sec 9. Profiling Selector

- Schema home: `monitor.yml -> profiling.timestep_output.mode`
- Canonical values:
  - `off`
  - `selected`
  - `all`
- Python hooks:
  - `resolve_profiling_config()`
  - `prepare_monitor_files()` in `scripts/picurv`
- Generated mapping:
  - `timestep_output.mode` -> `-profiling_timestep_mode`
  - `timestep_output.file` -> `-profiling_timestep_file`
  - `final_summary.enabled` -> `-profiling_final_summary`
  - selected function list -> `profile.run` + `-profile_config_file`
- C storage/parser:
  - profiling fields in `SimCtx`
  - option parsing in @ref CreateSimulationContext
- Runtime:
  - `src/logging.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 09_Monitor_Reference**

@section p50_post_sec 10. Postprocessing / Statistics Tasks

- Schema home: `post.yml -> statistics_pipeline`
- Canonical values:
  - `msd`
- Python hook:
  - `normalize_statistics_task()` in `scripts/picurv`
- Generated mapping:
  - `msd` -> `ComputeMSD` token in generated `post.run` (dispatched to `ComputeParticleMSD`)
- C dispatch:
  - `src/postprocessor.c`
- Kernel declaration/implementation:
  - `include/particle_statistics.h`
  - `src/particle_statistics.c`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 10_Post_Processing_Reference**

@section p50_related_sec 11. Related Pages

- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 31_Momentum_Solvers**
- **@subpage 44_Boundary_Conditions_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Modular Selector Extension Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

