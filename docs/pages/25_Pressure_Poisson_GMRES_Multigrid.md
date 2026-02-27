@page 25_Pressure_Poisson_GMRES_Multigrid Pressure-Poisson, GMRES, and Multigrid

This page describes the pressure-correction solve path used by projection in PICurv.

@tableofcontents

@section equation_sec 1. Pressure-Correction Equation

The correction solve enforces incompressibility through:

\f[
\nabla^2 \phi = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*.
\f]

In code terms:

- RHS assembly: divergence source in @ref PoissonRHS
- LHS assembly: metric-aware operator in @ref PoissonLHSNew
- solve orchestration: @ref PoissonSolver_MG

Null-space handling is explicitly configured for Neumann-like pressure systems via function @ref PoissonNullSpaceFunction in the Poisson module.

@section mg_sec 2. Multigrid/KSP Stack In Code

@ref PoissonSolver_MG currently:

1. assembles per-level operators,
2. configures `KSP` + `PCMG`,
3. sets restriction/interpolation operators (@ref MyRestriction and @ref MyInterpolation plus solid-aware variants),
4. applies level smoothers/coarse solve,
5. solves finest-level system for `Phi`.

After Poisson solve:

- pressure is updated by @ref UpdatePressure
- velocity is projected by @ref Projection

@section config_sec 3. YAML Mapping and PETSc Options

From `solver.yml` via `scripts/pic.flow`:

- `pressure_solver.tolerance` -> `-poisson_tol`
- `pressure_solver.multigrid.levels` -> `-mg_level`
- `pressure_solver.multigrid.pre_sweeps` -> `-mg_pre_it`
- `pressure_solver.multigrid.post_sweeps` -> `-mg_post_it`
- `pressure_solver.multigrid.semi_coarsening.{i,j,k}` -> `-mg_i_semi`, `-mg_j_semi`, `-mg_k_semi`
- optional level solver keys -> `-ps_mg_levels_*`
- `petsc_passthrough_options` -> raw PETSc flags

Final option parsing happens in function @ref CreateSimulationContext during context creation.

@section robustness_sec 4. Robustness Characteristics

Current implementation includes:

- periodic-boundary pressure synchronization,
- immersed-boundary-aware treatment paths (`Nvert`/solid checks),
- optional Poisson monitor logging (`-ps_ksp_pic_monitor_true_residual`).

If pressure solve quality degrades, check first:

1. BC consistency,
2. MG level/smoother settings,
3. grid metrics/orientation quality,
4. solver tolerances vs timestep.

@section refs_sec 5. Related Pages

- **@subpage 23_Fractional_Step_Method**
- **@subpage 24_Dual_Time_Picard_RK4**
- **@subpage 31_Momentum_Solvers**
