@page 08_Solver_Reference Configuration Reference: Solver Profiles (`solver.yml`)

For a complete, heavily commented reference file showing every possible option, please see the master template:

@verbinclude master_template/master_solver.yml

A **Solver Profile** is a `.yml` file that defines the complete numerical strategy for the C-solver. It controls everything from the core time-stepping scheme to the fine-grained tolerances of the linear solvers.

The PICurv platform is designed to be modular: you can combine any `case.yml` with any `solver.yml` at runtime. This allows you to test different numerical methods on the same physical problem without ever changing the case definition.

This document serves as a reference for all available sections and parameters within a `solver.yml` file.

@tableofcontents

@section operation_mode_sec 1. The `operation_mode` Section

This section controls how Eulerian fields are sourced and, for analytical mode, which analytical family is requested.

```yaml
operation_mode:
  eulerian_field_source: "solve"
  analytical_type: "TGV3D"
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `eulerian_field_source` | String | Selects Eulerian source mode. Common values: `"solve"`, `"load"`, `"analytical"`. | `-euler_field_source` |
| `analytical_type` | String | Analytical solution selector used when `eulerian_field_source` is analytical. | `-analytical_type` |

@section strategy_sec 2. The `strategy` Section

This section controls the high-level numerical scheme for the momentum equations.

```yaml
strategy:
  momentum_solver: "Dual Time Picard RK4"
  central_diff: false
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `momentum_solver` | String | Preferred momentum solver selector. User-friendly values like `"Explicit RK4"` and `"Dual Time Picard RK4"` are translated by `pic-flow` to C enum strings. | `-mom_solver_type` |
| `central_diff` | Boolean | If `true`, uses a central differencing scheme for the convective term. If `false` (default), uses a higher-order QUICK scheme. | `-central` |

@section tolerances_sec 3. The `tolerances` Section

This section configures momentum-solver convergence controls (primarily used by dual-time solvers).

```yaml
tolerances:
  max_iterations: 50
  absolute_tol: 1.0e-7
  relative_tol: 1.0e-4
  step_tol: 1.0e-8
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `max_iterations` | Integer | The maximum number of pseudo-time iterations for the momentum solver per time step. | `-mom_max_pseudo_steps` |
| `absolute_tol` | Real | The absolute tolerance for the momentum update norm. | `-mom_atol` |
| `relative_tol` | Real | The relative tolerance for the momentum update norm. | `-mom_rtol` |
| `step_tol` | Real | The step tolerance. The solver converges if the norm of the solution update is smaller than this value. | `-imp_stol` |

@section momentum_solver_sec 3.1. The `momentum_solver` Section (Advanced)

This optional section exposes solver-specific controls directly in YAML.

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `type` | String | Optional canonical enum fallback (`DUALTIME_PICARD_RK4`, `EXPLICIT_RK`, etc.) if `strategy.momentum_solver` is omitted. | `-mom_solver_type` |
| `dual_time_picard_rk4.max_pseudo_steps` | Integer | Maximum pseudo-time iterations per physical step. | `-mom_max_pseudo_steps` |
| `dual_time_picard_rk4.absolute_tol` | Real | Absolute pseudo-time convergence threshold. | `-mom_atol` |
| `dual_time_picard_rk4.relative_tol` | Real | Relative pseudo-time convergence threshold. | `-mom_rtol` |
| `dual_time_picard_rk4.step_tol` | Real | Step-size convergence threshold. | `-imp_stol` |
| `dual_time_picard_rk4.pseudo_cfl.initial` | Real | Initial pseudo-CFL at the start of each physical step. | `-pseudo_cfl` |
| `dual_time_picard_rk4.pseudo_cfl.minimum` | Real | Lower bound for adaptive pseudo-CFL. | `-min_pseudo_cfl` |
| `dual_time_picard_rk4.pseudo_cfl.maximum` | Real | Upper bound for adaptive pseudo-CFL. | `-max_pseudo_cfl` |
| `dual_time_picard_rk4.pseudo_cfl.growth_factor` | Real | Growth factor used when convergence is healthy. | `-pseudo_cfl_growth_factor` |
| `dual_time_picard_rk4.pseudo_cfl.reduction_factor` | Real | Reduction factor used when convergence degrades. | `-pseudo_cfl_reduction_factor` |
| `dual_time_picard_rk4.rk4_residual_noise_allowance_factor` | Real | Noise threshold for RK4 residual-based divergence detection. | `-mom_dt_rk4_residual_norm_noise_allowance_factor` |

@section pressure_solver_sec 4. The `pressure_solver` Section

This section provides detailed control over the linear solver used for the pressure-Poisson equation, which is often the most computationally expensive part of the simulation.

```yaml
pressure_solver:
  tolerance: 5.0e-9
  multigrid:
    levels: 3
    pre_sweeps: 1
    post_sweeps: 1
    # ... etc. ...
```

@subsection p_main_ssec 3.1. Main Settings

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `tolerance` | Real | The relative convergence tolerance for the pressure-Poisson linear solve. | `-poisson_tol` |

@subsection p_mg_ssec 3.2. `multigrid` Settings

These parameters configure the geometric multigrid (PCMG) preconditioner.

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `levels` | Integer | The total number of grid levels to use in the V-cycle (including the finest). | `-mg_level` |
| `pre_sweeps` | Integer | The number of smoothing iterations to perform before restricting to a coarser grid. | `-mg_pre_it` |
| `post_sweeps` | Integer | The number of smoothing iterations to perform after interpolating from a coarser grid. | `-mg_post_it` |
| `semi_coarsening.i`| Bool Array | An array of booleans (one per block) specifying whether to disable coarsening in the i-direction. | `-mg_i_semi` |
| `semi_coarsening.j`| Bool Array | An array of booleans (one per block) specifying whether to disable coarsening in the j-direction. | `-mg_j_semi` |
| `semi_coarsening.k`| Bool Array | An array of booleans (one per block) specifying whether to disable coarsening in the k-direction. | `-mg_k_semi` |

@section petsc_sec 5. The `petsc_passthrough_options` Section

This is an advanced section that allows you to pass any command-line option directly to the underlying PETSc library. This gives you complete control over the linear solvers (KSP) and preconditioners (PC).

The key is the PETSc option prefix. For the main pressure solver, the prefix is **`-ps_`**.

**Example:** To change the main KSP solver to GMRES and its preconditioner to ILU on all multigrid levels:
```yaml
petsc_passthrough_options:
  -ps_ksp_type: "gmres"
  -ps_pc_type: "ilu"
```

**Example:** To configure the solver on a specific multigrid level, use the `-ps_mg_levels_<level_num>_` prefix. This is useful for setting a direct solver (LU) on the coarsest grid:
```yaml
petsc_passthrough_options:
  # Set the coarsest level (level 0) to use a direct LU solver
  -ps_mg_levels_0_ksp_type: "preonly"
  -ps_mg_levels_0_pc_type: "lu"
  # Set the smoother on all finer levels to be 2 iterations of SOR
  -ps_mg_levels_ksp_type: "richardson"
  -ps_mg_levels_ksp_max_it: 2
  -ps_mg_levels_pc_type: "sor"
```
For a full list of PETSc options, please refer to the [PETSc Users Manual](https://petsc.org/release/docs/manual/).

@section next_steps_sec 6. Next Steps

You now know how to define the physics of a problem (`case.yml`) and control the numerical strategy used to solve it (`solver.yml`). The next step is to learn how to manage the simulation's output and runtime monitoring.

Proceed to the **@subpage 09_Monitor_Reference** to learn about all the parameters available in the `monitor.yml` file.

For end-to-end mapping and extension workflow, see **@subpage 15_Config_Ingestion_Map** and **@subpage 16_Config_Extension_Playbook**.
