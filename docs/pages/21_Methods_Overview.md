@page 21_Methods_Overview Methods and Models Overview

This section maps PICurv's numerical methods to the code paths that execute each step.
It is intended as the bridge between theory-level terminology and what the current codebase actually does.

@tableofcontents

@section governing_sec 1. Governing Model Snapshot

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

@section runtime_map_sec 2. Runtime Execution Order

At runtime, the top-level sequence is:

1. @ref CreateSimulationContext parses generated control files.
2. @ref InitializeEulerianState sets initial Eulerian fields.
3. @ref FlowSolver advances one fluid step (or @ref AnalyticalSolutionEngine for analytical mode).
4. Particle stage (when active): interpolation -> motion -> relocation -> physics -> scatter.

The method pages below document each major stage in detail.

@section method_map_sec 3. Method Map

- **@subpage 22_CURVIB_Method**: curvilinear grid/metric framework and immersed-boundary context.
- **@subpage 23_Fractional_Step_Method**: predictor/projection incompressible update.
- **@subpage 24_Dual_Time_Picard_RK4**: implicit-in-physical-time momentum iteration.
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**: Poisson assembly and PETSc multigrid/KSP path.
- **@subpage 26_Walking_Search_Method**: particle cell location and migration orchestration.
- **@subpage 27_Trilinear_Interpolation_and_Projection**: Eulerian-Lagrangian field exchange.
- **@subpage 28_IEM_and_Statistical_Averaging**: particle micromixing and post statistics kernels.
- **@subpage 31_Momentum_Solvers**: supported momentum-solver options and dispatch status.
- **@subpage 32_Analytical_Solutions**: analytical Eulerian modes and geometry policies.
- **@subpage 33_Initial_Conditions**: Eulerian and particle initialization modes.
- **@subpage 34_Particle_Model_Overview**: end-to-end particle model pipeline.
- **@subpage 44_Boundary_Conditions_Guide**: boundary handler combinations and runtime enforcement path.
- **@subpage 45_Particle_Initialization_and_Restart**: detailed particle seeding/restart/migration behavior.
- **@subpage 46_C_Runtime_Execution_Map**: code-level startup and timestep execution map.

@section usage_sec 4. How To Read These Pages

Use each page with two goals:

1. theory alignment (which equation/model is being approximated),
2. implementation alignment (which function actually executes it now).

When behavior differs from classical textbook formulations, the implementation notes take precedence for this repository.
