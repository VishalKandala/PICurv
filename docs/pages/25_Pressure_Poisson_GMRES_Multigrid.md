@page 25_Pressure_Poisson_GMRES_Multigrid Pressure-Poisson, GMRES, and Multigrid

@anchor _Pressure_Poisson_GMRES_Multigrid

This page describes the pressure-correction solve path used by projection in PICurv.

@tableofcontents

@section p25_equation_sec 1. Pressure-Correction Equation

The correction solve enforces incompressibility through:

\f[
\nabla^2 \phi = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*.
\f]

In code terms:

- RHS assembly: divergence source in @ref PoissonRHS
- LHS assembly: metric-aware operator in @ref PoissonLHSNew
- solve orchestration: @ref PoissonSolver_MG

Null-space handling is explicitly configured for Neumann-like pressure systems via function @ref PoissonNullSpaceFunction in the Poisson module.

@section p25_mg_sec 2. Multigrid/KSP Stack In Code

@ref PoissonSolver_MG currently:

1. assembles per-level operators,
2. configures `KSP` + `PCMG`,
3. sets restriction/interpolation operators (@ref MyRestriction and @ref MyInterpolation plus solid-aware variants),
4. applies level smoothers/coarse solve,
5. solves finest-level system for `Phi`.

After Poisson solve:

- pressure is updated by @ref UpdatePressure
- velocity is projected by @ref Projection

@section p25_config_sec 3. YAML Mapping and PETSc Options

From `solver.yml` via `scripts/picurv`:

- `pressure_solver.tolerance` -> `-poisson_tol`
- `pressure_solver.multigrid.levels` -> `-mg_level`
- `pressure_solver.multigrid.pre_sweeps` -> `-mg_pre_it`
- `pressure_solver.multigrid.post_sweeps` -> `-mg_post_it`
- `pressure_solver.multigrid.semi_coarsening.{i,j,k}` -> `-mg_i_semi`, `-mg_j_semi`, `-mg_k_semi`
- optional level solver keys -> `-ps_mg_levels_*`
- `petsc_passthrough_options` -> raw PETSc flags

Final option parsing happens in function @ref CreateSimulationContext during context creation.

@section p25_robustness_sec 4. Robustness Characteristics

Current implementation includes:

- periodic-boundary pressure synchronization,
- immersed-boundary-aware treatment paths (`Nvert`/solid checks),
- optional Poisson monitor logging (`-ps_ksp_pic_monitor_true_residual`).

If pressure solve quality degrades, check first:

1. BC consistency,
2. MG level/smoother settings,
3. grid metrics/orientation quality,
4. solver tolerances vs timestep.

@section p25_testing_sec 5. Current test status

Current direct tests are strongest for helper and invariant behavior:

- `PoissonLHSNew`
- `Projection`
- `PoissonNullSpaceFunction`
- RHS-related helpers used by `ComputeRHS`

The main remaining gap is `PoissonSolver_MG`: it is exercised in runtime smoke, but still lacks equivalent direct bespoke coverage for debugging. Periodic and immersed-boundary stencil branches also remain thinner than the core Cartesian helper surface.

@section p25_refs_sec 6. Related Pages

- **@subpage 23_Fractional_Step_Method**
- **@subpage 24_Dual_Time_Picard_RK4**
- **@subpage 31_Momentum_Solvers**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Pressure-Poisson, GMRES, and Multigrid** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
