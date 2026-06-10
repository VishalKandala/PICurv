@page 24_Dual_Time_Picard_Jameson_RK Dual-Time Picard Jameson RK Momentum Solver

@anchor _Dual_Time_Picard_Jameson_RK

This page documents the currently active implicit-style momentum strategy in PICurv: dual-time Picard iteration with Jameson RK pseudo-time smoothing.

@tableofcontents

@section p24_model_sec 1. Algorithmic Model

The solver enforces the physical-time momentum equation by iterating a pseudo-time equation:

\f[
\frac{\partial \mathbf{U}}{\partial \tau} = -\Big(R_{spatial}(\mathbf{U}) + R_{time}(\mathbf{U})\Big).
\f]

Implementation details from `ComputeTotalResidual`:

- `R_spatial` comes from @ref ComputeRHS
- `R_time` uses BDF1/BDF2-style terms from `Ucont`, `Ucont_o`, and `Ucont_rm1`,
- Jameson RK smoothing stages use coefficients \f$\{1/4,1/3,1/2,1\}\f$.
- Each stage evaluates a fresh residual at the previous stage state, then forms
  the next stage from the fixed pseudo-iteration base state. This is not
  classical fourth-order Runge-Kutta time integration.

@section p24_convergence_sec 2. Convergence and Adaptive Rollback

Per pseudo-iteration, the solver tracks:

- \f$\|\Delta U\|_\infty\f$,
- relative solution update \f$\|\Delta U_k\|/\|\Delta U_0\|\f$,
- residual norms and growth ratios.

Adaptive pseudo-CFL rollback is triggered when residual growth exceeds the
configured noise allowance or a non-finite trial is detected, then:

1. state rolls back to the last accepted pseudo-state,
2. pseudo-CFL is reduced,
3. iteration retries from the restored state.

Acceptance and rollback are global across blocks and MPI ranks. Rejected
trials do not update accepted residual history, and `max_iterations` bounds
total attempted trials rather than only accepted trials. A finite solve that
exhausts its attempts exits with the last accepted finite state.

Pseudo-CFL is adaptively ramped on successful trials, reduced on noisy accepted
trials or rejection, and clamped by configured min/max bounds. The
controller-selected next CFL carries directly into the next physical timestep;
there is no separate end-of-step rebound.

Convergence has a compatibility switch:

- both residual tolerances disabled: update absolute **and** relative criteria,
- either residual tolerance enabled: update absolute **or** relative criterion,
  combined with the enabled residual absolute/relative criteria.

@section p24_config_sec 3. YAML -> Runtime Controls

User-facing configuration (`solver.yml`) maps to:

- `strategy.momentum_solver` -> `-mom_solver_type`
- `tolerances.max_iterations` -> `-mom_max_pseudo_steps`
- `tolerances.absolute_tol` -> `-mom_atol`
- `tolerances.relative_tol` -> `-mom_rtol`
- `tolerances.residual_absolute_tol` -> `-mom_resid_atol`
- `tolerances.residual_relative_tol` -> `-mom_resid_rtol`
- `momentum_solver.dual_time_picard_jameson_rk.pseudo_cfl.*` -> pseudo-CFL flags
- `jameson_residual_noise_allowance_factor` -> `-mom_dt_jameson_residual_norm_noise_allowance_factor`

The former `Dual Time Picard RK4`, `dual_time_picard_rk4`,
`rk4_residual_noise_allowance_factor`, `DUALTIME_PICARD_RK4`, and
`-mom_dt_rk4_residual_norm_noise_allowance_factor` spellings remain deprecated
compatibility aliases. Canonical configuration and generated controls use the
Jameson names.

Parsing and normalization are performed in `scripts/picurv`, with final option ingestion in function @ref CreateSimulationContext during setup.
Only the currently implemented momentum solver values are exposed; add new ones
only when the parser and dispatcher are extended in the same change.

@section p24_logging_sec 4. Logging and Diagnostics

The persistent momentum convergence-history format is unchanged.
`PseudoIter(k)` now denotes total attempted trials, including rejected trials.
Controller ratios, rollback decisions, CFL changes, and internal counters are
available at `DEBUG`; finite nonconvergence is reported at `WARNING`; and each
solve emits a concise final summary at `INFO`.

@section p24_touchpoints_sec 5. Core Code Touchpoints

- implementation: @ref MomentumSolver_DualTime_Picard_JamesonRK
- explicit comparator path: @ref MomentumSolver_Explicit_RungeKutta4
- runtime dispatch: @ref FlowSolver
- options ingestion: @ref CreateSimulationContext

@section p24_practical_sec 6. Practical Tuning Guidance

Common stability tuning order:

1. start from the shipped `growth_factor: 1.05` and `reduction_factor: 0.75`,
2. reduce initial or maximum pseudo-CFL if trials repeatedly reject,
3. use `residual_relative_tol: 1.0e-3` for robust profiles or `1.0e-2` for
   exploratory LES where looser inner convergence is acceptable,
4. tighten/loosen residual-noise allowance slightly only after examining logs,
5. revisit physical timestep, grid quality, and boundary consistency if
   instability persists.

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
