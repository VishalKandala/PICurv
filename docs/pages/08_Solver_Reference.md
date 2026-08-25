@page 08_Solver_Reference Configuration Reference: Solver YAML

@anchor _Solver_Reference

For the full commented template, see:

@verbinclude master_template/master_solver.yml

`solver.yml` controls numerical strategy and solver internals.

@tableofcontents

@section p08_opmode_sec 1. operation_mode

```yaml
operation_mode:
  eulerian_field_source: "solve"
  analytical_type: "TGV3D"
  uniform_flow:
    u: 0.0
    v: 0.0
    w: 0.0
```

Mappings:
- `eulerian_field_source` -> `-euler_field_source` (`solve`, `load`, `analytical`)
- `analytical_type` -> `-analytical_type`
- `uniform_flow.u/v/w` -> `-analytical_uniform_u/-analytical_uniform_v/-analytical_uniform_w` when `analytical_type: "UNIFORM_FLOW"`

`uniform_flow` is only valid when `analytical_type: "UNIFORM_FLOW"`.

@section p08_strategy_sec 2. strategy

```yaml
strategy:
  momentum_solver: "Dual Time Picard Jameson RK"
  central_diff: false
```

Mappings:
- `momentum_solver` -> `-mom_solver_type` (`picurv` accepts `Explicit RK4`, `Dual Time Picard Jameson RK`, or `Newton Krylov`)
- `central_diff` -> `-central`

Older boolean toggles are not supported; use `strategy.momentum_solver`.
Only implemented momentum solver values are accepted by `picurv` and the C runtime.
The deprecated `Dual Time Picard RK4` display name, `dual_time_picard_rk4`
solver block, and `rk4_residual_noise_allowance_factor` key remain readable
compatibility aliases; generated controls always use the canonical Jameson names.

@section p08_tol_sec 3. tolerances

```yaml
tolerances:
  max_iterations: 50
  absolute_tol: 1.0e-7
  relative_tol: 1.0e-4
  residual_absolute_tol: 0.0
  residual_relative_tol: 1.0e-3
```

Mappings:
- `max_iterations` -> `-mom_max_pseudo_steps`
- `absolute_tol` -> `-mom_atol`
- `relative_tol` -> `-mom_rtol`
- `residual_absolute_tol` -> `-mom_resid_atol`
- `residual_relative_tol` -> `-mom_resid_rtol`

When both residual tolerances are non-positive, the solver preserves its
existing update-only convergence criterion. `step_tol`/`-imp_stol` remains
accepted as a deprecated compatibility option but is unused by active
momentum solvers.

@section p08_msolver_sec 4. momentum_solver (Solver-Specific Block)

```yaml
momentum_solver:
  dual_time_picard_jameson_rk:
    max_pseudo_steps: 50
    absolute_tol: 1.0e-8
    relative_tol: 1.0e-5
    pseudo_cfl:           # dimensionless Courant number: dtau = pseudo_cfl / lambda_max
      initial: 0.5
      minimum: 0.001
      maximum: 2.0
      growth_factor: 1.1
      reduction_factor: 0.75
    jameson_residual_noise_allowance_factor: 1.1
```

Mappings include:
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- `-pseudo_cfl_growth_factor`, `-pseudo_cfl_reduction_factor`
- `-mom_dt_jameson_residual_norm_noise_allowance_factor`

Rule: solver-specific blocks must match selected momentum solver type.
Do not set canonical Jameson keys and their deprecated RK4 aliases together.
Pseudo-CFL is exclusively a Dual Time Picard--Jameson RK control: it is neither
accepted from the structured Newton--Krylov block nor shown in a Newton--Krylov
startup banner.

For `strategy.momentum_solver: "Newton Krylov"`, the structured PETSc controls are:

```yaml
momentum_solver:
  newton_krylov:
    jacobian:
      type: finite_difference
      finite_difference:
        mode: matrix_free
    preconditioner:
      model: none
    nonlinear_solver:
      method: newtonls
      absolute_tolerance: 1.0e-10
      relative_tolerance: 1.0e-8
      step_tolerance: 1.0e-12
      max_iterations: 12
      line_search:
        type: bt
    linear_solver:
      method: gmres
      absolute_tolerance: 1.0e-10
      relative_tolerance: 1.0e-6
      max_iterations: 400
      gmres:
        restart: 80
```

Mappings are `nonlinear_solver.method/absolute_tolerance/relative_tolerance/step_tolerance/max_iterations`
to `-mom_nk_snes_type/-mom_nk_snes_atol/-mom_nk_snes_rtol/-mom_nk_snes_stol/-mom_nk_snes_max_it`,
`line_search.type` to `-mom_nk_snes_linesearch_type`, and the corresponding
`linear_solver` fields to `-mom_nk_ksp_type/-mom_nk_ksp_atol/-mom_nk_ksp_rtol/-mom_nk_ksp_max_it`.
`gmres.restart` maps to `-mom_nk_ksp_gmres_restart`. The Jacobian fields map to
`-mom_nk_jacobian_type/-mom_nk_jacobian_fd_mode`; the preconditioner
fields map to `-mom_nk_preconditioner_model/-mom_nk_preconditioner_structure`.

Newton tolerances are nonnegative, iteration/restart counts are positive integers,
and GMRES restart is valid only for `gmres`, `fgmres`, or `lgmres`. Supported
combinations are finite difference/matrix free with either no
preconditioner or frozen momentum Jacobian/point block. The released
`linear_solver.preconditioner.type: none` remains a deprecated alias for the
no-preconditioner model. Raw `petsc_passthrough_options` are applied last, but an
incompatible raw `-mom_nk_pc_type` override is rejected by the runtime.
The Jacobian block is a strict discriminated configuration: an explicit
`type: finite_difference` requires `finite_difference.mode: matrix_free`.
`colored_sparse`, `frozen_momentum_approximation`, and irrelevant sibling
configuration are rejected because their implementations are not present.

@section p08_poisson_sec 5. poisson_solver

```yaml
poisson_solver:
  method: "fgmres"
  absolute_tolerance: 1.0e-5
  relative_tolerance: 1.0e-11
  max_iterations: 50
  gmres:
    restart: 20
  preconditioner:
    type: "multigrid"
  multigrid:
    levels: 3
    pre_sweeps: 2
    post_sweeps: 2
    cycle: "v"
    mode: "multiplicative"
    semi_coarsening:
      i: false
      j: false
      k: true
    level_solvers:
      level_0:
        method: "preonly"
        preconditioner: "redundant"
      level_1:
        method: "richardson"
        preconditioner: "bjacobi"
      level_2:
        method: "richardson"
        preconditioner: "bjacobi"
```

Mappings:
- `method` -> `-ps_ksp_type`
- `absolute_tolerance` -> `-ps_ksp_atol` and legacy `-poisson_tol`
- `relative_tolerance` -> `-ps_ksp_rtol`
- `max_iterations` -> `-ps_ksp_max_it`
- `gmres.restart` -> `-ps_ksp_gmres_restart`; valid only for `gmres`, `fgmres`, or `lgmres`
- `preconditioner.type` -> `-ps_pc_type`; currently only `multigrid` is supported
- `multigrid.levels` -> `-mg_level`
- `multigrid.pre_sweeps` -> `-mg_pre_it`
- `multigrid.post_sweeps` -> `-mg_post_it`
- `multigrid.semi_coarsening.i/j/k` -> `-mg_i_semi/-mg_j_semi/-mg_k_semi`
- `multigrid.level_solvers.level_N.method` -> `-ps_mg_levels_N_ksp_type` for `N > 0`
- `multigrid.level_solvers.level_N.preconditioner` -> `-ps_mg_levels_N_pc_type` for `N > 0`
- `multigrid.level_solvers.level_0.*` -> `-ps_mg_coarse_*`; PETSc names the coarsest
  solver separately from the positive levels
- `multigrid.cycle` and `multigrid.mode` are validated structured keys; current supported values are `v` and `multiplicative`.

Rules:
- `pressure_solver` is accepted as a legacy alias, but `poisson_solver` is preferred because the linear solve computes pressure correction `Phi`.
- MG level numbering follows PETSc/PICurv convention: `level_0` is the coarsest level and larger numbers are finer.
- `level_0` is the multigrid **coarse solve**, not a smoother, and the naming hides that. `level_1..N` are
  smoothers; `level_0` sits at the base of the V-cycle and removes the smooth error the smoothers cannot see.
  Because multigrid is used here as a preconditioner, `level_0` must be a **fixed linear operator**. A Krylov
  method there (`gmres`, `fgmres`, `cg`, `bcgs`, ...) makes the whole preconditioner nonlinear, which decouples
  the outer FGMRES tracked residual from the true residual `b - Ax`: the solver then reports convergence on a
  number that no longer describes the constraint. Pinning `ksp_max_it` does not fix it. PICurv logs a startup
  warning if you configure one. Use `{method: preonly, preconditioner: redundant}` for coarse grids up to
  roughly 1e4 unknowns, `telescope` above that. Full discussion, the worked tracked-vs-true residual table, and
  the level-count sizing rule: @ref 25_Pressure_Poisson_GMRES_Multigrid.
- `multigrid.levels` is bounded by the MPI decomposition, not chosen freely: every level must leave each rank
  at least `stencil_width` nodes per axis (3 when any axis is periodic, 2 otherwise). Exceeding it aborts during
  DM creation with `Local x-width of domain ... is smaller than stencil width`, and `picurv validate` cannot
  catch it because it does not see the rank layout. Formula and worked maxima:
  @ref 25_Pressure_Poisson_GMRES_Multigrid.
- The outer Poisson preconditioner is multigrid-only in the current runtime.
- The current PETSc binding applies one MG smoother count; when `pre_sweeps` and `post_sweeps` differ, PICurv uses the larger value and logs a warning.
- Advanced PETSc tuning remains available through `petsc_passthrough_options`; common examples include `-ps_mg_levels_N_pc_sor_omega` for SOR and `-ps_mg_levels_N_pc_factor_shift_amount` / `-ps_mg_levels_N_pc_factor_levels` for factor PCs.

@section p08_solution_conv_sec 6. Physical-solution convergence monitoring

Convergence monitoring is observation policy rather than a numerical solver
selection, so it is configured in `monitor.yml -> solution_monitoring.convergence`
rather than here. See **@subpage 09_Monitor_Reference** for modes, mappings, and
defaults. The monitor records every completed timestep; cadence is not a user
setting.

`solver.yml` accepts no `solution_convergence` key. A file carrying one is
rejected by validation with the location to move it to.

Scientific field statistics are likewise a monitor concern, configured at
`monitor.yml -> field_statistics`; see **@subpage 58_Field_Statistics**.

@section p08_interp_sec 7. interpolation

```yaml
interpolation:
  method: "Trilinear"     # default; or "CornerAveraged"
```

Mappings:
- `method` -> `-interpolation_method` (`Trilinear` = `0`, `CornerAveraged` = `1`)

The **Trilinear** method (default) performs direct trilinear interpolation from the 8 nearest cell centers, providing second-order accuracy on both uniform and curvilinear grids. The **CornerAveraged** method is the legacy two-stage path (center-to-corner average, then trilinear from corners), which is second-order only on uniform Cartesian grids.

See **@subpage 27_Trilinear_Interpolation_and_Projection** for algorithmic details.

@section p08_scalar_transport_sec 8. scalar_transport

```yaml
scalar_transport:
  schmidt_number: 1.0
  turbulent_schmidt_number: 0.7
```

Mappings:
- `schmidt_number` -> `-schmidt_number`
- `turbulent_schmidt_number` -> `-turb_schmidt_number`

Rules:
- values must be positive numbers
- omitted values use the C runtime defaults: `schmidt_number = 1.0` and `turbulent_schmidt_number = 0.7`
- use this structured block for ordinary scalar/Brownian transport tuning; reserve `petsc_passthrough_options` for flags without a YAML schema

@section p08_verification_sec 9. verification

```yaml
verification:
  sources:
    diffusivity:
      mode: "analytical"
      profile: "LINEAR_X"
      gamma0: 1.0e-3
      slope_x: 2.0e-4

    scalar:
      mode: "analytical"
      profile: "CONSTANT"  # or LINEAR_X, SIN_PRODUCT
      value: 1.0
      # phi0: 0.0
      # slope_x: 1.0
      # amplitude: 1.0
      # kx: 3.141592653589793
      # ky: 3.141592653589793
      # kz: 3.141592653589793
```

Mappings:
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

Rules:
- this path is verification-only and should be used only when no ordinary end-to-end workflow can expose the behavior under test
- it is only valid with `operation_mode.eulerian_field_source: "analytical"`
- `verification.sources.scalar` prescribes particle `Psi` from analytical truth and enables the runtime diagnostic `logs/scatter_metrics.csv`
- scalar profiles currently supported are `CONSTANT`, `LINEAR_X`, and `SIN_PRODUCT`
- new verification source overrides must be implemented in `include/verification_sources.h` and `src/verification_sources.c`

@section p08_petsc_sec 10. petsc_passthrough_options

Advanced escape hatch for raw PETSc flags:

```yaml
petsc_passthrough_options:
  -ps_ksp_type: "gmres"
  -ps_pc_type: "ilu"
```

These are passed into PETSc options DB and consumed by runtime calls like `KSPSetFromOptions`.

@section p08_next_steps_sec 11. Next Steps

Proceed to **@subpage 09_Monitor_Reference**.

For mapping and extension workflows:
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 31_Momentum_Solvers**
- **@subpage 24_Dual_Time_Picard_Jameson_RK**
- **@ref 55_Newton_Krylov_Momentum_Solver** — matrix-free Newton--Krylov momentum solver.
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
- **@subpage 50_Modular_Selector_Extension_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Configuration Reference: Solver YAML** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
