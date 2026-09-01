@page 21_Methods_Overview Methods and Models Overview

@anchor _Methods_Overview

@pagemeta{Explanation, Readers evaluating the numerics, Current behavior}

This section maps PICurv's numerical methods to the code paths that execute each step.
It is intended as the bridge between theory-level terminology and what the current codebase actually does.

@tableofcontents

@section p21_governing_sec 1. Governing Model Snapshot

PICurv advances incompressible flow in non-dimensional form on curvilinear grids:

\f[
\frac{\partial \mathbf{u}}{\partial t} + \nabla\cdot(\mathbf{u}\mathbf{u}) = -\nabla p + \frac{1}{Re}\nabla^2\mathbf{u} + \mathbf{f},
\qquad
\nabla\cdot\mathbf{u}=0.
\f]

Operationally, the solver uses a projection workflow:

1. build momentum residuals,
2. advance momentum with selected strategy,
3. solve Poisson equation for pressure correction,
4. project velocity to divergence-free space,
5. execute particle coupling (if enabled).

@section p21_runtime_map_sec 2. Runtime Execution Order

At runtime, the top-level sequence is:

1. @ref CreateSimulationContext parses generated control files.
2. @ref InitializeEulerianState sets initial Eulerian fields.
3. @ref FlowSolver advances one fluid step (or @ref AnalyticalSolutionEngine for analytical mode).
4. Particle stage (when active): interpolation -> motion -> relocation -> physics -> scatter.

The method pages below document each major stage in detail.

@section p21_method_map_sec 3. Method Map

- **@subpage 22_CURVIB_Method**: curvilinear grid/metric framework and immersed-boundary context.
- **@subpage 23_Fractional_Step_Method**: predictor/projection incompressible update.
- **@subpage 24_Dual_Time_Picard_Jameson_RK**: implicit-in-physical-time momentum iteration.
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**: Poisson assembly and PETSc multigrid/KSP path.
- **@subpage 26_Walking_Search_Method**: particle cell location and migration orchestration.
- **@subpage 27_Trilinear_Interpolation_and_Projection**: Eulerian-Lagrangian field exchange.
- **@subpage 28_IEM_and_Statistical_Averaging**: particle micromixing and post statistics kernels.
- **@subpage 72_LES_Turbulence_Closure**: subgrid closure, the Germano-Lilly dynamic procedure, coefficient averaging, and limiting.
- **@subpage 31_Momentum_Solvers**: supported momentum-solver options and dispatch status.
- **@subpage 32_Analytical_Solutions**: analytical Eulerian modes and geometry policies.
- **@subpage 33_Initial_Conditions**: Eulerian and particle initialization modes.
- **@subpage 55_Newton_Krylov_Momentum_Solver**: matrix-free PETSc SNES/GMRES momentum solver and its deterministic-residual invariant.
- **@subpage 34_Particle_Model_Overview**: end-to-end particle model pipeline.
- **@subpage 44_Boundary_Conditions_Guide**: boundary handler combinations and runtime enforcement path.
- **@subpage 45_Particle_Initialization_and_Restart**: detailed particle seeding/restart/migration behavior.
- **@subpage 46_C_Runtime_Execution_Map**: code-level startup and timestep execution map.

@section p21_references_sec 4. External Method References

For readers connecting PICurv implementation to the CURVIB literature, these are the primary starting references:

- Borazjani I, Ge L, Sotiropoulos F. "Curvilinear immersed boundary method for simulating fluid structure interaction with complex 3D rigid bodies." *Journal of Computational Physics* 227(16), 7587-7620 (2008). DOI: `10.1016/j.jcp.2008.04.024`.
- Borazjani I, Di Achille P, D'Souza RM, et al. "The functional role of left atrial flow in ventricular filling and flow evolution in the left ventricle." *Annals of Biomedical Engineering* 41(6), 1265-1275 (2013). DOI: `10.1007/s10439-013-0758-9`.

Project context:

- PICurv project page: `https://vishalkandala.me/projects/Picurv/`

@section p21_usage_sec 5. How To Read These Pages

Use each page with two goals:

1. theory alignment (which equation/model is being approximated),
2. implementation alignment (which function actually executes it now).

When behavior differs from classical textbook formulations, the implementation notes take precedence for this repository.

@section p21_authority_sec 6. What Makes a Method Page Authoritative

Each method page is a contract between users and maintainers, not a literature
summary. Read it in this order:

1. **Scope and assumptions** — supported grid, boundary, solver, and restart
   conditions; unsupported combinations must be stated explicitly.
2. **Discrete/runtime behavior** — equation or algorithm, field locations,
   update order, and the C functions that implement it.
3. **User controls** — the YAML fields, defaults, generated PETSc options, and
   logs that show whether the method is behaving as configured.
4. **Validation evidence** — the unit, smoke, or regression checks that protect
   the stated behavior, including any limits on what has been exercised.
5. **Extension boundary** — which parser, data structure, runtime dispatch, and
   tests must change together when adding a model or selector.

A numerical claim is authoritative only when it can be traced to the current
runtime path and a named validation surface. When this page and an older
example, paper, or legacy option disagree, use the current configuration
reference and runtime-specific method page, then correct the stale material.
