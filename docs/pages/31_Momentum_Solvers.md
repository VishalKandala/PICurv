@page 31_Momentum_Solvers Momentum Solver Implementations

@anchor _Momentum_Solvers

This page tracks momentum-solver options accepted by the current configuration and their runtime implementation status.

@tableofcontents

@section p31_selection_sec 1. Selection and Dispatch

Runtime selection is controlled by `-mom_solver_type`, produced from `solver.yml` (`strategy.momentum_solver`).
Dispatch currently happens in function @ref FlowSolver within the step orchestrator.

Accepted YAML values:

- `Explicit RK4` -> `EXPLICIT_RK`
- `Dual Time Picard Jameson RK` -> `DUALTIME_PICARD_JAMESON_RK`

Only implemented values are exposed in the enum, parser, and dispatcher. New
solver values should be added only with a real implementation plus matching
parser, docs, and test updates.

For compatibility, the former `Dual Time Picard RK4` YAML display name and
`DUALTIME_PICARD_RK4` C CLI value still select the Jameson solver. New
configuration and code must use the Jameson names.

@section p31_status_sec 2. Implementation Status Matrix

- `EXPLICIT_RK`: implemented by @ref MomentumSolver_Explicit_RungeKutta4
- `DUALTIME_PICARD_JAMESON_RK`: implemented by @ref MomentumSolver_DualTime_Picard_JamesonRK

@section p31_controls_sec 3. Numerical Controls In Use

Main controls consumed by implemented solvers:

- `-mom_max_pseudo_steps`
- `-mom_atol`
- `-mom_rtol`
- `-mom_resid_atol`, `-mom_resid_rtol`
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- `-pseudo_cfl_growth_factor`, `-pseudo_cfl_reduction_factor`
- `-mom_dt_jameson_residual_norm_noise_allowance_factor`
- `-mom_ratio_ema_alpha`

Defaults and final option ingestion are in function @ref CreateSimulationContext during startup parsing.

For the dual-time Jameson solver, `max_iterations` bounds **accepted** pseudo-iterations. A separate
hard cap of `3 × max_iterations` limits total attempts (accepted plus rejected) to prevent infinite
rejection loops. Convergence requires both the update pass (`|ΔU| ≤ atol` **AND** `|ΔU|/|ΔU₀| ≤ rtol`)
and, when either residual tolerance is positive, at least one enabled residual criterion to hold.

The dual-time controller uses one global pseudo-CFL and globally accepts or
rolls back a complete four-stage trial. The selected next pseudo-CFL is carried
directly into the next physical timestep. `step_tol`/`-imp_stol` remains
accepted only as a deprecated compatibility input and is unused by active
momentum solvers.

`pseudo_cfl.*` values are **dimensionless Courant numbers** (Phase 3+), not fractions of the physical
timestep `dt`. The solver computes `dtau = pseudo_cfl / lambda_max` where `lambda_max` is the global
maximum convective spectral radius of the current velocity field. This makes `pseudo_cfl` independent
of `dt`, grid size, and flow speed. The stable range for the 4-stage Jameson RK smoother is
`pseudo_cfl ≈ 0–2.83`; the shipped defaults are `initial: 0.5`, `maximum: 2.0`.

@section p31_testing_sec 4. Current test status

Current testing is uneven by solver path:

- dispatch and guardrails are directly covered through `FlowSolver`-side unit tests
- `MomentumSolver_DualTime_Picard_JamesonRK` is exercised mainly through smoke and runtime orchestration
- `MomentumSolver_Explicit_RungeKutta4` still needs a direct positive-path harness

That means the momentum stack is currently a stronger regression gate than bespoke debugging surface.

@section p31_extension_sec 5. Adding A New Momentum Solver

Required steps:

1. define solver implementation function in `src/momentumsolvers.c`,
2. ensure enum and parser mapping are present (`variables.h`, `setup.c`, `picurv_cli/core.py`),
3. add dispatch branch in function @ref FlowSolver for the new enum value,
4. expose and document solver-specific YAML options,
5. add smoke tests and docs updates.

For user-facing contract updates, also update:

- **@subpage 08_Solver_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 50_Modular_Selector_Extension_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Momentum Solver Implementations** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
