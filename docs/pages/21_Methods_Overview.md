@page 21_Methods_Overview Methods and Models Overview

@anchor _Methods_Overview

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

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Methods and Models Overview** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
