@page 23_Fractional_Step_Method Fractional-Step (Projection) Method

PICurv uses a predictor-corrector style incompressible update: predict momentum, solve pressure-Poisson, project velocity to divergence-free space.

@tableofcontents

@section model_sec 1. High-Level Algorithm

Per physical timestep:

1. Compute momentum predictor (provisional velocity).
2. Solve pressure-Poisson equation from continuity constraint.
3. Correct velocity with pressure gradient (projection step).
4. Apply boundary/refresh operations and proceed.

This enforces incompressibility through the projection/correction stage.

@section code_sec 2. Code Touchpoints

- High-level step orchestration: @ref FlowSolver
- Momentum RHS assembly: @ref ComputeRHS
- Pressure solve stage: @ref PoissonSolver_MG
- Pressure/velocity correction: @ref UpdatePressure, @ref Projection

@section refs_sec 3. References

- Chorin projection-method background: https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
- PETSc KSP manual: https://petsc.org/release/manual/ksp/
