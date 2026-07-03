@page 24_Dual_Time_Picard_Jameson_RK Dual-Time Picard Jameson RK Momentum Solver

@anchor _Dual_Time_Picard_Jameson_RK

This page documents one of PICurv's two implicit-in-physical-time momentum
strategies: the **dual-time Picard fixed-point / pseudo-time** solver, using
Jameson RK pseudo-time smoothing. It is the established, broadly exercised default.
For the alternative matrix-free Newton solver and for selection guidance between
the two, see @ref 55_Newton_Krylov_Momentum_Solver and @ref 31_Momentum_Solvers.

@tableofcontents

@section p24_when_to_use_sec 1. When to Use This Solver

This solver implements global pseudo-time explicit smoothing on top of an implicit BDF2 physical-time discretization. It is the recommended momentum solver when:

- The flow is **laminar or low-to-moderate Re (Re ≤ ~2000)** on body-fitted curvilinear grids.
- A **larger physical timestep** is required than explicit stability allows — typically 5–20× the explicit CFL limit.
- The flow geometry involves **curves, bends, or bodies** where curvilinear coordinates are essential and the explicit RK4 solver would need prohibitively small `dt`.
- **5–30 pseudo-iterations per physical step** is an acceptable overhead (typical for smooth, slowly-varying flows).
- **Robustness over throughput** is the priority: automatic CFL rollback and EMA-smoothed rejection reduce manual tuning for stable convergence.

Use the **Explicit RK4** solver instead when:

- The explicit stability limit is affordable (coarse grid, short verification, or low-Re 1D-like flows).
- No pseudo-time truncation error is acceptable (verification against analytical solutions).
- Per-step throughput dominates cost and no implicit stability is needed.

**Pseudo-CFL semantics (Phase 3+):** the pseudo-time step is `dtau = pseudo_cfl / lambda_max`, where `lambda_max` is the global maximum convective spectral radius computed from the velocity field at the start of each physical timestep. This makes `pseudo_cfl` a true dimensionless Courant number independent of `dt`, grid size, or flow speed — `pseudo_cfl = 1.0` has the same physical meaning across all grids and all `dt`. The stable range for the 4-stage Jameson RK smoother is `pseudo_cfl ≈ 0–2.83`; default `initial: 0.5` gives a comfortable margin. See Section 2 for the spectral-radius estimation details.

@section p24_model_sec 2. Algorithmic Model

The solver enforces the physical-time momentum equation by iterating a pseudo-time equation:

\f[
\frac{\partial \mathbf{U}}{\partial \tau} = -\Big(R_{spatial}(\mathbf{U}) + R_{time}(\mathbf{U})\Big).
\f]

Implementation details from `ComputeTotalResidual`:

- `R_spatial` comes from @ref ComputeRHS
- `R_time` uses BDF2-style terms from `Ucont`, `Ucont_o`, and `Ucont_rm1` (BDF1 on the first physical step)
- Jameson RK smoothing stages use coefficients \f$\{1/4,1/3,1/2,1\}\f$
- Each stage evaluates a fresh residual at the previous stage state, then forms the next stage from the fixed pseudo-iteration base state — this is a 4-stage dissipation smoother, not classical 4th-order RK time integration
- The pseudo-time step is `dtau = pseudo_cfl / lambda_max`, where `lambda_max` is the global maximum convective spectral radius; stability requires `pseudo_cfl < C_RK ≈ 2.83`

**Spectral radius estimation** (`ComputeGlobalSpectralRadiusEstimate`): at the start of each physical timestep (after BCs and ghost synchronization), the solver computes per-cell:

\f[
\lambda_{cell} = \left(|\tilde{U}_{x}| + |\tilde{U}_{y}| + |\tilde{U}_{z}|\right) \cdot J^{-1}
\f]

where \f$\tilde{U}\f$ components are the contravariant volume fluxes (`ucont`, units m³/s) and \f$J^{-1} = 1/\text{volume}\f$ (`lAj`). The product yields units [1/s]. A global MPI `MAX` reduction gives `lambda_max`. A lower bound `PetscMax(lambda_max, 1.5/dt)` (the BDF2 accuracy coefficient divided by `dt`) prevents division by zero at startup or in zero-flow regions, falling back to `dtau ≈ pseudo_cfl × dt / 1.5`.

@section p24_convergence_sec 3. Convergence and Adaptive Rollback

Per pseudo-iteration, the solver tracks:

- \f$\|\Delta U\|_\infty\f$ (solution update norm),
- \f$\|\Delta U_k\|/\|\Delta U_0\|\f$ (relative update),
- residual norms and their step-to-step ratio (EMA-smoothed).

Adaptive pseudo-CFL rollback is triggered when the EMA-smoothed step-to-step residual ratio exceeds the configured noise allowance (`jameson_residual_noise_allowance_factor`, default 1.1), or a non-finite trial is detected:

1. State rolls back to the last accepted pseudo-state.
2. Pseudo-CFL is multiplied by `reduction_factor`.
3. Iteration retries from the restored state without consuming the accepted-iteration budget.

Acceptance and rollback are global across blocks and MPI ranks. The `max_iterations` parameter bounds **accepted** pseudo-iterations. A separate hard cap of `3 × max_iterations` bounds total attempts (accepted plus rejected) to prevent infinite rejection loops. A finite solve that exhausts its accepted-iteration budget exits with the last accepted finite state.

Pseudo-CFL is adaptively ramped on successful trials (ratio < 0.90: immediate growth; 0.90–1.0: growth after 3 consecutive clean trials), reduced on noisy accepted trials or rejection, and clamped by configured min/max bounds. The controller-selected next CFL carries directly into the next physical timestep.

**Convergence criteria** (both must pass simultaneously):

- **Update pass**: `|ΔU| ≤ absolute_tol` **AND** `|ΔU|/|ΔU₀| ≤ relative_tol`
- **Residual pass** (only checked when either residual tolerance is set to a positive value): at least one of `|R| ≤ mom_resid_atol` or `|R|/|R₀| ≤ mom_resid_rtol` must hold

Both update and residual passes must be true simultaneously for the solver to declare convergence.

@section p24_config_sec 4. YAML → Runtime Controls

User-facing configuration (`solver.yml`) maps to:

- `strategy.momentum_solver` → `-mom_solver_type`
- `tolerances.max_iterations` → `-mom_max_pseudo_steps`
- `tolerances.absolute_tol` → `-mom_atol`
- `tolerances.relative_tol` → `-mom_rtol`
- `tolerances.residual_absolute_tol` → `-mom_resid_atol`
- `tolerances.residual_relative_tol` → `-mom_resid_rtol`
- `momentum_solver.dual_time_picard_jameson_rk.pseudo_cfl.*` → pseudo-CFL flags
- `jameson_residual_noise_allowance_factor` → `-mom_dt_jameson_residual_norm_noise_allowance_factor` (default: 1.1)
- `ratio_ema_alpha` → `-mom_ratio_ema_alpha` (default: 0.3; range [0, 1])

The `ratio_ema_alpha` parameter controls EMA smoothing of the step-to-step residual ratio before the rejection decision:
```
smoothed = alpha × raw_ratio + (1 − alpha) × smoothed_prev
```
`alpha = 1.0` recovers the original raw-ratio behavior (most aggressive rejection). `alpha = 0.3` (default) requires approximately 3–4 consecutive bad trials before triggering rejection, tolerating transient residual bumps common in convection-dominated flows.

The former `Dual Time Picard RK4`, `dual_time_picard_rk4`,
`rk4_residual_noise_allowance_factor`, `DUALTIME_PICARD_RK4`, and
`-mom_dt_rk4_residual_norm_noise_allowance_factor` spellings remain deprecated
compatibility aliases. Canonical configuration and generated controls use the
Jameson names.

Parsing and normalization are performed in `picurv_cli/core.py`, with final option ingestion in function @ref CreateSimulationContext during setup.
Only the currently implemented momentum solver values are exposed; add new ones
only when the parser and dispatcher are extended in the same change.

@section p24_logging_sec 5. Logging and Diagnostics

The persistent momentum convergence-history format (`logs/Momentum_Solver_Convergence_History_Block_N.log`) includes per-trial fields:

- `PseudoIter(k)`: total attempted trial index (includes rejected trials)
- `dtau`: physical-time pseudo-step used for this trial [s] — equals `pseudo_cfl / lambda_max`
- `cfl_eff`: effective dimensionless Courant number for this trial — equals `dtau × lambda_max`; controlled by `pseudo_cfl.*` YAML keys
- `|dUk|`, `|dUk|/|dU0|`: solution update norms
- `|Rk|`, `|Rk|/|R0|`: residual norms
- `trial_ratio`: raw step-to-step residual ratio
- `smoothed_ratio`: EMA-smoothed ratio used for the rejection decision
- `status`: `accepted` or `rejected`
- `dtau_after`: physical-time pseudo-step selected for the next trial [s]
- `cfl_eff_after`: corresponding dimensionless Courant number for the next trial

A startup INFO log line prints the active CFL bounds, rejection threshold, EMA alpha, growth/reduction factors, and iteration budget. Non-convergence (exhausted accepted-iteration budget) is reported unconditionally via `PetscPrintf`. Internal ratios, rollback decisions, and CFL changes are logged at DEBUG.

@section p24_touchpoints_sec 6. Core Code Touchpoints

- implementation: @ref MomentumSolver_DualTime_Picard_JamesonRK
- explicit comparator path: @ref MomentumSolver_Explicit_RungeKutta4
- runtime dispatch: @ref FlowSolver
- options ingestion: @ref CreateSimulationContext

@section p24_practical_sec 7. Practical Tuning Guidance

Common stability tuning order:

1. Start from the shipped defaults: `initial: 0.5`, `maximum: 2.0`, `growth_factor: 1.1`, `reduction_factor: 0.75`, `jameson_residual_noise_allowance_factor: 1.1`, `ratio_ema_alpha: 0.3`. `pseudo_cfl` is now a dimensionless Courant number (Phase 3+); `initial: 0.5` sits at ~18% of the 4-stage Jameson stability limit (2.83) and is a safe universal starting point regardless of `dt` or grid size.
2. If trials repeatedly reject: reduce `pseudo_cfl.maximum` and/or `pseudo_cfl.initial`. For flows near the stability limit, try `maximum: 1.5` first; `pseudo_cfl = 2.83` is the theoretical convection-stability limit, so practical `maximum` should not exceed 2.5.
3. If the residual history is non-monotonic (common for convection-dominated flows): raise `ratio_ema_alpha` toward 0.5–0.7 to make the EMA respond faster, or raise `jameson_residual_noise_allowance_factor` to 1.2–1.3.
4. Use `residual_relative_tol: 1.0e-3` for robust production runs or `1.0e-2` for exploratory LES where looser inner convergence is acceptable.
5. Tighten/loosen the residual-noise allowance only after examining the `smoothed_ratio` column in the convergence log.
6. If instability persists after CFL tuning: revisit physical timestep, grid quality near bends or walls, and boundary condition consistency.

For many cases, robust Poisson settings and sane initialization matter as much as dual-time tolerances.

For contributor extension steps, see **@subpage 50_Modular_Selector_Extension_Guide**.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Dual-Time Picard Jameson RK Momentum Solver** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
