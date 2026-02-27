@page 22_CURVIB_Method CurvIB Method Overview

PICurv solves flow on structured curvilinear grids and applies immersed-boundary-aware logic in metric, pressure, and projection stages.

@tableofcontents

@section formulation_sec 1. Curvilinear Formulation Context

The code evolves contravariant flux-like velocity components (`Ucont`) and derives Cartesian velocity (`Ucat`) as needed.
Curvilinear metrics map physical derivatives to computational coordinates:

\f[
\nabla \phi = \phi_{\xi}\,\mathbf{\xi} + \phi_{\eta}\,\mathbf{\eta} + \phi_{\zeta}\,\mathbf{\zeta},
\qquad
J^{-1} = \frac{\partial(\xi,\eta,\zeta)}{\partial(x,y,z)}.
\f]

Face and cell metric tensors are precomputed and reused by RHS, Poisson, and projection stages.

@section setup_sec 2. Grid and Metric Build Pipeline

Main setup touchpoints:

- grid dimensions and decomposition: @ref DefineAllGridDimensions
- DM creation per block and level: @ref InitializeAllGridDMs
- coordinate assignment: @ref AssignAllGridCoordinates
- face metrics: @ref ComputeFaceMetrics
- cell-centered Jacobian inverse: @ref ComputeCellCenteredJacobianInverse
- orientation checks/fixups: @ref CheckAndFixGridOrientation

Useful geometric helper for BC and flux logic:

- @ref CalculateFaceCenterAndArea

@section ibm_sec 3. Immersed-Boundary Role In Current Code

The current branch keeps immersed-boundary hooks in solver paths, while several IBM-specific calls are conditional or currently inactive in default runs.
In practice:

- geometry masking uses `Nvert`/solid markers in key kernels,
- Poisson and projection stencils include boundary-aware logic,
- IBM-specific interpolation hooks remain extension points for fully active IBM runs.

@section practical_sec 4. What This Means For Users

You mainly control CurvIB behavior through:

- grid generation/ingestion choice,
- boundary-condition handlers,
- solver settings that influence projection/Poisson robustness,
- case geometry and decomposition settings.

Even when users never interact with metric tensors directly, they govern stability, pressure correction quality, and particle coupling accuracy.

@section refs_sec 5. Related Pages

- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@subpage 23_Fractional_Step_Method**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
