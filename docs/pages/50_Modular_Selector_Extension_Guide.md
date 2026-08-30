@page 50_Modular_Selector_Extension_Guide Modular Selector Extension Guide

@anchor _Modular_Selector_Extension_Guide

This page is the contributor checklist for adding new selector-driven options
without drifting the YAML -> launcher -> C runtime contract.

@tableofcontents

@section p50_general_sec 1. Standard Checklist For Any Selector

For every new selector value:

1. define the canonical user-facing value in YAML docs/templates,
2. add or update Python normalization/validation in `picurv_cli/core.py`,
3. emit the correct generated control/post artifact mapping,
4. add any C enum/storage and parser wiring,
5. add the runtime dispatch/consumer branch,
6. add tests,
7. **declare the value in the capability registry**, and
8. **write its Tier-2 capability entry**.

Use the canonical value only. Do not add placeholder enum values or compatibility aliases for unimplemented options.

@subsection p50_obligations_ssec 1.1 The Documentation Obligation Is Enforced

Steps 7 and 8 are not advisory. Coverage is enforced for every registered family, so
`make audit-capability` fails until they are done, and the full source chain -
validator, C parser, enum, runtime dispatch - must agree. If the value belongs to a
family nothing yet covers, `make audit-family-census` will have already refused the
new public surface.

Concretely, adding a selector value means:

```bash
# 1. Declare it: add an entry to value_metadata in
#    tests/tooling/capability_families.json with its status, canonical/alias
#    relationship, and evidence sources.

# 2. Scaffold the entry with every required part already present:
python3 tests/tooling/scaffold_documentation.py capability \
    --family <family.id> --value '<Your Value>'

# 3. Paste it into the family page, fill in the TODOs, then:
make docs-inventory      # regenerate the inventory and its tables
make audit-capability    # parity, coverage, evidence, lifecycle
```

The audit checks four things a reviewer would otherwise have to check by hand:

- **Parity** - the value is accepted end to end, and nothing accepts a value the
  validator does not expose.
- **Coverage** - an entry exists, anchored outside code blocks, with all eight
  contract parts (three for a deprecated alias).
- **Evidence** - each declared source exists *and* is cited in the entry's Evidence
  part, so the registry and the prose cannot drift.
- **Lifecycle** - a `supported` value has evidence, an `experimental` or
  `known-defective` value states limitations, a `deprecated` value gives migration.

If the value is a synonym rather than a new capability, declare `spelling_of` instead
and it inherits the canonical value's status without needing its own entry. If it is a
retired name, declare `alias_of` and scaffold the three-part stub with
`scaffold_documentation.py alias`.

Full contract: **@subpage 64_Documentation_Extension_Framework**.

@section p50_momentum_sec 2. Momentum Solver Selector

- Schema home: `solver.yml -> strategy.momentum_solver`
- Canonical values:
  - `Explicit RK4`
  - `Dual Time Picard Jameson RK`
  - `Newton Krylov`
- Deprecated alias:
  - `Dual Time Picard RK4` normalizes to `Dual Time Picard Jameson RK`; the C enum
    keeps `MOMENTUM_SOLVER_DUALTIME_PICARD_RK4` as an alias of the Jameson value.
- Python hook:
  - `normalize_momentum_solver_type()` in `picurv_cli/core.py`
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
  - **@subpage 24_Dual_Time_Picard_Jameson_RK**
  - **@subpage 31_Momentum_Solvers**
  - **@subpage 55_Newton_Krylov_Momentum_Solver**

@section p50_bc_sec 3. Boundary Condition Type Or Handler

- Schema home: `case.yml -> boundary_conditions`
- Canonical values:
  - see the supported type/handler matrix in **@subpage 44_Boundary_Conditions_Guide**
- Python hooks:
  - `BC_TYPE_MAP`
  - `BC_HANDLER_SPECS`
  - `validate_and_prepare_boundary_conditions()` in `picurv_cli/core.py`
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
  - `normalize_particle_init_mode()` in `picurv_cli/core.py`
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
  - `normalize_interpolation_method()` in `picurv_cli/core.py`
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
  - `generated` with `generator: zero|constant|streamwise_constant|poiseuille|ic_gen`
  - `file` with `field: Ucat|Ucont`
- Python hooks:
  - `resolve_initial_condition_config()`
  - `stage_initial_condition_file()`
  - `run_initial_condition_generator()`
  - `normalize_flow_direction_token()`
  - `resolve_ic_cli_params()`
- Generated mapping:
  - built-in generator -> matching `-finit` and parameter flags
  - file or `ic_gen` -> `-finit 4 -ic_field ... -ic_dir ...`
- C storage/parser:
  - `initialConditionMode`, `initialConditionField`, `initialConditionDirectory`,
    `icVelocityPhysical`, and `flowDirection` in `SimCtx`
- Runtime:
  - `src/initialcondition.c` — `PopulateInitialUcont`, `SetInitialInteriorField`
  - `src/setup.c` — `Cart2Contra`, `UniformCart2Contra`
- Tests/docs to update:
  - `tests/test_cli_smoke.py`
  - **@subpage 33_Initial_Conditions**

@section p50_analytical_sec 7. Analytical Type

- Schema home: `solver.yml -> operation_mode.analytical_type`
- Canonical values:
  - `TGV3D`
  - `ZERO_FLOW`
  - `UNIFORM_FLOW`
- Python hook:
  - `normalize_analytical_type()` in `picurv_cli/core.py`
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

@section p50_verification_sources_sec 8. Verification Source Overrides

- Schema home: `solver.yml -> verification.sources.*`
- Design intent:
  - reserved for verification-only source injections when no ordinary end-to-end path is sufficient
  - keep runtime call sites thin and centralize override implementations in `verification_sources.*`
- Current support:
  - `verification.sources.diffusivity`
  - `verification.sources.scalar`
  - diffusivity profile `LINEAR_X`
  - scalar profiles `CONSTANT`, `LINEAR_X`, `SIN_PRODUCT`
- Python hooks:
  - validation and flag emission in `picurv_cli/core.py`
- Generated mapping:
  - `verification.sources.diffusivity.mode` -> `-verification_diffusivity_mode`
  - `verification.sources.diffusivity.profile` -> `-verification_diffusivity_profile`
  - `verification.sources.diffusivity.gamma0` -> `-verification_diffusivity_gamma0`
  - `verification.sources.diffusivity.slope_x` -> `-verification_diffusivity_slope_x`
  - `verification.sources.scalar.mode` -> `-verification_scalar_mode`
  - `verification.sources.scalar.profile` -> `-verification_scalar_profile`
  - `verification.sources.scalar.value` -> `-verification_scalar_value`
  - `verification.sources.scalar.phi0` -> `-verification_scalar_phi0`
  - `verification.sources.scalar.slope_x` -> `-verification_scalar_slope_x`
  - `verification.sources.scalar.amplitude` -> `-verification_scalar_amplitude`
  - `verification.sources.scalar.kx/ky/kz` -> `-verification_scalar_kx/-verification_scalar_ky/-verification_scalar_kz`
- C storage/parser:
  - `verificationDiffusivity` and `verificationScalar` in `SimCtx`
  - option parsing in @ref CreateSimulationContext
- Runtime:
  - `src/verification_sources.c`
  - thin delegation from `src/rhs.c` for diffusivity
  - analytical scalar truth in `src/AnalyticalSolutions.c`
  - gated particle-physics bypass in `src/ParticlePhysics.c`
  - runtime metric writer in `src/logging.c` (`<run.runtime_logs>/scatter_metrics.csv`)
- Design boundary:
  - verification sources exist only for otherwise-unreachable verification scenarios
  - analytical truth definitions belong in `AnalyticalSolutions`, not in model-evolution code such as `ParticlePhysics`
- Tests/docs to update:
  - `tests/test_config_regressions.py`
  - `tests/c/test_poisson_rhs.c`
  - `tests/c/test_solver_kernels.c`
  - `tests/c/test_logging.c`
  - `tests/c/test_setup_lifecycle.c`
  - **@subpage 08_Solver_Reference**
  - **@subpage 32_Analytical_Solutions**
  - **@subpage 16_Config_Extension_Playbook**

@section p50_grid_sec 9. Grid Selector / Generator Selector

- Schema homes:
  - `case.yml -> grid.mode`
  - `case.yml -> grid.generator.grid_type`
- Canonical values:
  - `programmatic_c`
  - `file`
  - `grid_gen`
- Python orchestration:
  - `picurv_cli/core.py` grid validation and run-dir staging
- Generator implementation:
  - `generators/grid.gen`
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

@section p50_profiling_sec 10. Profiling Selector

- Schema home: `monitor.yml -> profiling.timestep_output.mode`
- Canonical values:
  - `off`
  - `selected`
  - `all`
- Python hooks:
  - `resolve_profiling_config()`
  - `prepare_monitor_files()` in `picurv_cli/core.py`
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

@section p50_post_sec 11. Postprocessing / Statistics Tasks

- Schema home: `post.yml -> statistics_pipeline`
- Canonical values:
  - `msd`
- Python hook:
  - `normalize_statistics_task()` in `picurv_cli/core.py`
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
