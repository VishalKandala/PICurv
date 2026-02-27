@page 24_Dual_Time_Picard_RK4 Dual-Time Picard RK4 Momentum Solver

PICurv’s current named momentum strategy is a dual-time fixed-point (Picard-style) iteration with RK4 pseudo-time smoothing.

@tableofcontents

@section algorithm_sec 1. Core Idea

For each physical timestep, solve momentum implicitly in pseudo-time:

- Iterate pseudo-steps until convergence criteria are met.
- Use RK4 staged updates to stabilize pseudo-time evolution.
- Check residual/update norms against absolute/relative/step thresholds.

Pseudo-CFL controls adapt iteration aggressiveness and robustness.

@section cfg_sec 2. User-Facing Controls

From `solver.yml`:
- `strategy.momentum_solver`
- `tolerances.*`
- `momentum_solver.dual_time_picard_rk4.*`

Mapped to C flags such as:
- `-mom_max_pseudo_steps`, `-mom_atol`, `-mom_rtol`, `-imp_stol`
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- growth/reduction factors and RK4 residual-noise allowance.

@section code_sec 3. Code Touchpoints

- Implemented dual-time solver: @ref MomentumSolver_DualTime_Picard_RK4
- Explicit RK fallback solver: @ref MomentumSolver_Explicit_RungeKutta4
- Runtime dispatch: @ref FlowSolver
- Option ingestion: @ref CreateSimulationContext (via `-mom_solver_type` and related controls)

@section refs_sec 4. References

- Classical RK formulations: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
- PETSc TS/KSP ecosystem (general): https://petsc.org/release/docs/manual/
