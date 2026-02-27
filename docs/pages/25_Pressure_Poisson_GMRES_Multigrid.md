@page 25_Pressure_Poisson_GMRES_Multigrid Pressure-Poisson, GMRES, and Multigrid

The pressure solve is the key linear system in the projection step. PICurv configures this through PETSc `KSP` with geometric multigrid support.

@tableofcontents

@section stack_sec 1. Solver Stack

- Build pressure linear system from continuity/projection formulation.
- Use PETSc `KSP` for linear iterations.
- Apply multigrid preconditioning (`PCMG`) across geometric levels.
- Use level-specific smoothers and coarse-level strategy through PETSc options.

@section cfg_sec 2. YAML Controls

From `solver.yml`:
- `pressure_solver.tolerance`
- `pressure_solver.multigrid.levels/pre_sweeps/post_sweeps`
- `pressure_solver.multigrid.semi_coarsening.*`
- `petsc_passthrough_options` (advanced PETSc tuning)

@section code_sec 3. Code Touchpoints

- Multigrid pressure solve entry: @ref PoissonSolver_MG
- Matrix/RHS assembly: @ref PoissonLHSNew, @ref PoissonRHS
- Correction/projection: @ref UpdatePressure, @ref Projection
- MG transfer helpers: @ref MyRestriction, @ref MyInterpolation
- Config ingestion: @ref CreateSimulationContext

@section notes_sec 4. Practical Notes

- The default profile uses robust settings; advanced PETSc tuning is optional.
- Coarse/fine smoother behavior is one of the most impactful performance knobs.

@section refs_sec 5. References

- PETSc KSP guide: https://petsc.org/release/manual/ksp/
- PETSc multigrid notes: https://petsc.org/release/manual/ksp/#multigrid-methods
