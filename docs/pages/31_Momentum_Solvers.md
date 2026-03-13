@page 31_Momentum_Solvers Momentum Solver Implementations

@anchor _Momentum_Solvers

This page tracks momentum-solver options accepted by the current configuration and their runtime implementation status.

@tableofcontents

@section p31_selection_sec 1. Selection and Dispatch

Runtime selection is controlled by `-mom_solver_type`, produced from `solver.yml` (`strategy.momentum_solver`).
Dispatch currently happens in function @ref FlowSolver within the step orchestrator.

Accepted YAML values:

- `Explicit RK4` -> `EXPLICIT_RK`
- `Dual Time Picard RK4` -> `DUALTIME_PICARD_RK4`

Only implemented values are exposed in the enum, parser, and dispatcher. New
solver values should be added only with a real implementation plus matching
parser, docs, and test updates.

@section p31_status_sec 2. Implementation Status Matrix

- `EXPLICIT_RK`: implemented by @ref MomentumSolver_Explicit_RungeKutta4
- `DUALTIME_PICARD_RK4`: implemented by @ref MomentumSolver_DualTime_Picard_RK4

@section p31_controls_sec 3. Numerical Controls In Use

Main controls consumed by implemented solvers:

- `-mom_max_pseudo_steps`
- `-mom_atol`
- `-mom_rtol`
- `-imp_stol`
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- `-pseudo_cfl_growth_factor`, `-pseudo_cfl_reduction_factor`
- `-mom_dt_rk4_residual_norm_noise_allowance_factor`

Defaults and final option ingestion are in function @ref CreateSimulationContext during startup parsing.

@section p31_testing_sec 4. Current test status

Current testing is uneven by solver path:

- dispatch and guardrails are directly covered through `FlowSolver`-side unit tests
- `MomentumSolver_DualTime_Picard_RK4` is exercised mainly through smoke and runtime orchestration
- `MomentumSolver_Explicit_RungeKutta4` still needs a direct positive-path harness

That means the momentum stack is currently a stronger regression gate than bespoke debugging surface.

@section p31_extension_sec 5. Adding A New Momentum Solver

Required steps:

1. define solver implementation function in `src/momentumsolvers.c`,
2. ensure enum and parser mapping are present (`variables.h`, `setup.c`, `scripts/picurv`),
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
