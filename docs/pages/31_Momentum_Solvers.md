@page 31_Momentum_Solvers Momentum Solver Implementations

This page tracks momentum-solver options accepted by the current configuration and their runtime implementation status.

@tableofcontents

@section selection_sec 1. Selection and Dispatch

Runtime selection is controlled by `-mom_solver_type`, produced from `solver.yml` (`strategy.momentum_solver`).
Dispatch currently happens in function @ref FlowSolver within the step orchestrator.

Accepted YAML values:

- `Explicit RK4` -> `EXPLICIT_RK`
- `Dual Time Picard RK4` -> `DUALTIME_PICARD_RK4`

Only implemented values are exposed in the enum, parser, and dispatcher. New
solver values should be added only with a real implementation plus matching
parser, docs, and test updates.

@section status_sec 2. Implementation Status Matrix

- `EXPLICIT_RK`: implemented by @ref MomentumSolver_Explicit_RungeKutta4
- `DUALTIME_PICARD_RK4`: implemented by @ref MomentumSolver_DualTime_Picard_RK4

@section controls_sec 3. Numerical Controls In Use

Main controls consumed by implemented solvers:

- `-mom_max_pseudo_steps`
- `-mom_atol`
- `-mom_rtol`
- `-imp_stol`
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- `-pseudo_cfl_growth_factor`, `-pseudo_cfl_reduction_factor`
- `-mom_dt_rk4_residual_norm_noise_allowance_factor`

Defaults and final option ingestion are in function @ref CreateSimulationContext during startup parsing.

@section extension_sec 4. Adding A New Momentum Solver

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
