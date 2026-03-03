@page 08_Solver_Reference Configuration Reference: Solver YAML

For the full commented template, see:

@verbinclude master_template/master_solver.yml

`solver.yml` controls numerical strategy and solver internals.

@tableofcontents

@section opmode_sec 1. operation_mode

```yaml
operation_mode:
  eulerian_field_source: "solve"
  analytical_type: "TGV3D"
```

Mappings:
- `eulerian_field_source` -> `-euler_field_source` (`solve`, `load`, `analytical`)
- `analytical_type` -> `-analytical_type`

@section strategy_sec 2. strategy

```yaml
strategy:
  momentum_solver: "Dual Time Picard RK4"
  central_diff: false
```

Mappings:
- `momentum_solver` -> `-mom_solver_type` (`picurv` accepts `Explicit RK4` or `Dual Time Picard RK4`)
- `central_diff` -> `-central`

Older boolean toggles are not supported; use `strategy.momentum_solver`.
Only implemented momentum solver values are accepted by `picurv` and the C runtime.

@section tol_sec 3. tolerances

```yaml
tolerances:
  max_iterations: 50
  absolute_tol: 1.0e-7
  relative_tol: 1.0e-4
  step_tol: 1.0e-8
```

Mappings:
- `max_iterations` -> `-mom_max_pseudo_steps`
- `absolute_tol` -> `-mom_atol`
- `relative_tol` -> `-mom_rtol`
- `step_tol` -> `-imp_stol`

@section msolver_sec 4. momentum_solver (Solver-Specific Block)

```yaml
momentum_solver:
  dual_time_picard_rk4:
    max_pseudo_steps: 50
    absolute_tol: 1.0e-8
    relative_tol: 1.0e-5
    step_tol: 1.0e-8
    pseudo_cfl:
      initial: 0.1
      minimum: 0.001
      maximum: 1.0
      growth_factor: 1.0
      reduction_factor: 1.0
    rk4_residual_noise_allowance_factor: 1.05
```

Mappings include:
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- `-pseudo_cfl_growth_factor`, `-pseudo_cfl_reduction_factor`
- `-mom_dt_rk4_residual_norm_noise_allowance_factor`

Rule: solver-specific blocks must match selected momentum solver type.

@section pressure_sec 5. pressure_solver

```yaml
pressure_solver:
  tolerance: 5.0e-9
  multigrid:
    levels: 3
    pre_sweeps: 1
    post_sweeps: 1
```

Mappings:
- `tolerance` -> `-poisson_tol`
- `multigrid.levels` -> `-mg_level`
- `multigrid.pre_sweeps` -> `-mg_pre_it`
- `multigrid.post_sweeps` -> `-mg_post_it`
- `multigrid.semi_coarsening.i/j/k` -> `-mg_i_semi/-mg_j_semi/-mg_k_semi`

@section petsc_sec 6. petsc_passthrough_options

Advanced escape hatch for raw PETSc flags:

```yaml
petsc_passthrough_options:
  -ps_ksp_type: "gmres"
  -ps_pc_type: "ilu"
```

These are passed into PETSc options DB and consumed by runtime calls like `KSPSetFromOptions`.

@section next_steps_sec 7. Next Steps

Proceed to **@subpage 09_Monitor_Reference**.

For mapping and extension workflows:
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 31_Momentum_Solvers**
- **@subpage 24_Dual_Time_Picard_RK4**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
- **@subpage 50_Modular_Selector_Extension_Guide**
