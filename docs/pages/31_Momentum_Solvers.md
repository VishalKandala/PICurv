@page 31_Momentum_Solvers Momentum Solver Implementations

This page tracks momentum-solver strategies, current implementation status, and extension points.

@tableofcontents

@section current_impl_sec 1. Current Runtime Selection

Selection is configured by `-mom_solver_type` and dispatched in @ref FlowSolver.

Supported names at ingestion:
- `EXPLICIT_RK`
- `DUALTIME_PICARD_RK4`
- `DUALTIME_NK_ARNOLDI` (enum recognized, runtime implementation placeholder)
- `DUALTIME_NK_ANALYTICAL_JACOBIAN` (enum recognized, runtime implementation placeholder)

@section implemented_sec 2. Implemented Solvers

- @ref MomentumSolver_Explicit_RungeKutta4
- @ref MomentumSolver_DualTime_Picard_RK4

The dual-time path uses additional controls (`-mom_max_pseudo_steps`, `-mom_atol`, `-mom_rtol`, pseudo-CFL controls).

@section touchpoints_sec 3. Key Code Touchpoints

- Selection and orchestration: @ref FlowSolver
- Solver implementations: `src/momentumsolvers.c`
- Option ingestion: @ref CreateSimulationContext
- Human-readable summaries: @ref MomentumSolverTypeToString

@section extension_sec 4. Adding a New Momentum Solver

1. Add enum + string mapping.
2. Add parser mapping for `-mom_solver_type`.
3. Implement solver function in `src/momentumsolvers.c`.
4. Add dispatch branch in @ref FlowSolver.
5. Expose structured YAML mapping in `scripts/pic.flow` and templates.
