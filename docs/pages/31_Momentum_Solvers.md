@page 31_Momentum_Solvers Momentum Solver Implementations

This page tracks momentum-solver options accepted by the current configuration and their runtime implementation status.

@tableofcontents

@section selection_sec 1. Selection and Dispatch

Runtime selection is controlled by `-mom_solver_type`, produced from `solver.yml` (`strategy.momentum_solver` and/or `momentum_solver.type`).
Dispatch currently happens in function @ref FlowSolver within the step orchestrator.

Accepted canonical values:

- `EXPLICIT_RK`
- `DUALTIME_PICARD_RK4`
- `DUALTIME_NK_ARNOLDI`
- `DUALTIME_NK_ANALYTICAL_JACOBIAN`

Normalization of user aliases is handled in `scripts/pic.flow` (`normalize_momentum_solver_type`).

@section status_sec 2. Implementation Status Matrix

- `EXPLICIT_RK`: implemented by @ref MomentumSolver_Explicit_RungeKutta4
- `DUALTIME_PICARD_RK4`: implemented by @ref MomentumSolver_DualTime_Picard_RK4
- `DUALTIME_NK_ARNOLDI`: enum/parser path exists; runtime implementation is placeholder in current dispatch path
- `DUALTIME_NK_ANALYTICAL_JACOBIAN`: enum/parser path exists; runtime implementation is placeholder in current dispatch path

If an unimplemented enum is selected, dispatch currently logs warning behavior.

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
2. ensure enum and parser mapping are present (`variables.h`, `setup.c`, `scripts/pic.flow`),
3. add dispatch branch in function @ref FlowSolver for the new enum value,
4. expose and document solver-specific YAML options,
5. add smoke tests and docs updates.

For user-facing contract updates, also update:

- **@subpage 08_Solver_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 40_Testing_and_Quality_Guide**
