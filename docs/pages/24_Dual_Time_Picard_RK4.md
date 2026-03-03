@page 24_Dual_Time_Picard_RK4 Dual-Time Picard RK4 Momentum Solver

This page documents the currently active implicit-style momentum strategy in PICurv: dual-time Picard iteration with RK4 pseudo-time smoothing.

@tableofcontents

@section model_sec 1. Algorithmic Model

The solver enforces the physical-time momentum equation by iterating a pseudo-time equation:

\f[
\frac{\partial \mathbf{U}}{\partial \tau} = -\Big(R_{spatial}(\mathbf{U}) + R_{time}(\mathbf{U})\Big).
\f]

Implementation details from `ComputeTotalResidual`:

- `R_spatial` comes from @ref ComputeRHS
- `R_time` uses BDF1/BDF2-style terms from `Ucont`, `Ucont_o`, and `Ucont_rm1`,
- RK4 stages use coefficients \f$\{1/4,1/3,1/2,1\}\f$.

@section convergence_sec 2. Convergence and Backtracking

Per pseudo-iteration, the solver tracks:

- \f$\|\Delta U\|_\infty\f$,
- relative solution update \f$\|\Delta U_k\|/\|\Delta U_0\|\f$,
- residual norms and growth ratios.

Backtracking is triggered when residual and update growth indicate divergence (or NaN), then:

1. state rolls back to the last accepted pseudo-state,
2. pseudo-CFL is reduced,
3. iteration retries from the restored state.

Pseudo-CFL is also adaptively ramped on successful steps and clamped by configured min/max bounds.

@section config_sec 3. YAML -> Runtime Controls

User-facing configuration (`solver.yml`) maps to:

- `strategy.momentum_solver` -> `-mom_solver_type`
- `tolerances.max_iterations` -> `-mom_max_pseudo_steps`
- `tolerances.absolute_tol` -> `-mom_atol`
- `tolerances.relative_tol` -> `-mom_rtol`
- `tolerances.step_tol` -> `-imp_stol`
- `momentum_solver.dual_time_picard_rk4.pseudo_cfl.*` -> pseudo-CFL flags
- `rk4_residual_noise_allowance_factor` -> `-mom_dt_rk4_residual_norm_noise_allowance_factor`

Parsing and normalization are performed in `scripts/picurv`, with final option ingestion in function @ref CreateSimulationContext during setup.
Only the currently implemented momentum solver values are exposed; add new ones
only when the parser and dispatcher are extended in the same change.

@section touchpoints_sec 4. Core Code Touchpoints

- implementation: @ref MomentumSolver_DualTime_Picard_RK4
- explicit comparator path: @ref MomentumSolver_Explicit_RungeKutta4
- runtime dispatch: @ref FlowSolver
- options ingestion: @ref CreateSimulationContext

@section practical_sec 5. Practical Tuning Guidance

Common stability tuning order:

1. reduce initial pseudo-CFL,
2. tighten/loosen residual-noise allowance slightly,
3. adjust pseudo-CFL growth/reduction factors,
4. revisit grid quality and boundary consistency if instability persists.

For many cases, robust Poisson settings and sane initialization matter as much as dual-time tolerances.

For contributor extension steps, see **@subpage 50_Modular_Selector_Extension_Guide**.
