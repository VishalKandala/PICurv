@page 24_Dual_Time_Picard_RK4 Dual-Time Picard RK4 Momentum Solver

@anchor _Dual_Time_Picard_RK4

This page documents the currently active implicit-style momentum strategy in PICurv: dual-time Picard iteration with RK4 pseudo-time smoothing.

@tableofcontents

@section p24_model_sec 1. Algorithmic Model

The solver enforces the physical-time momentum equation by iterating a pseudo-time equation:

\f[
\frac{\partial \mathbf{U}}{\partial \tau} = -\Big(R_{spatial}(\mathbf{U}) + R_{time}(\mathbf{U})\Big).
\f]

Implementation details from `ComputeTotalResidual`:

- `R_spatial` comes from @ref ComputeRHS
- `R_time` uses BDF1/BDF2-style terms from `Ucont`, `Ucont_o`, and `Ucont_rm1`,
- RK4 stages use coefficients \f$\{1/4,1/3,1/2,1\}\f$.

@section p24_convergence_sec 2. Convergence and Backtracking

Per pseudo-iteration, the solver tracks:

- \f$\|\Delta U\|_\infty\f$,
- relative solution update \f$\|\Delta U_k\|/\|\Delta U_0\|\f$,
- residual norms and growth ratios.

Backtracking is triggered when residual and update growth indicate divergence (or NaN), then:

1. state rolls back to the last accepted pseudo-state,
2. pseudo-CFL is reduced,
3. iteration retries from the restored state.

Pseudo-CFL is also adaptively ramped on successful steps and clamped by configured min/max bounds.

@section p24_config_sec 3. YAML -> Runtime Controls

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

@section p24_touchpoints_sec 4. Core Code Touchpoints

- implementation: @ref MomentumSolver_DualTime_Picard_RK4
- explicit comparator path: @ref MomentumSolver_Explicit_RungeKutta4
- runtime dispatch: @ref FlowSolver
- options ingestion: @ref CreateSimulationContext

@section p24_practical_sec 5. Practical Tuning Guidance

Common stability tuning order:

1. reduce initial pseudo-CFL,
2. tighten/loosen residual-noise allowance slightly,
3. adjust pseudo-CFL growth/reduction factors,
4. revisit grid quality and boundary consistency if instability persists.

For many cases, robust Poisson settings and sane initialization matter as much as dual-time tolerances.

For contributor extension steps, see **@subpage 50_Modular_Selector_Extension_Guide**.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Dual-Time Picard RK4 Momentum Solver** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

