@page 22_CURVIB_Method CurvIB Method Overview

@anchor _CURVIB_Method

PICurv solves flow on structured curvilinear grids using the curvilinear formulation of CURVIB. It has no immersed-boundary capability.

@tableofcontents

@section p22_formulation_sec 1. Curvilinear Formulation Context

The code evolves contravariant flux-like velocity components (`Ucont`) and derives Cartesian velocity (`Ucat`) as needed.
Curvilinear metrics map physical derivatives to computational coordinates:

\f[
\nabla \phi = \phi_{\xi}\,\mathbf{\xi} + \phi_{\eta}\,\mathbf{\eta} + \phi_{\zeta}\,\mathbf{\zeta},
\qquad
J^{-1} = \frac{\partial(\xi,\eta,\zeta)}{\partial(x,y,z)}.
\f]

Face and cell metric tensors are precomputed and reused by RHS, Poisson, and projection stages.

@section p22_setup_sec 2. Grid and Metric Build Pipeline

Main setup touchpoints:

- grid dimensions and decomposition: @ref DefineAllGridDimensions
- DM creation per block and level: @ref InitializeAllGridDMs
- coordinate assignment: @ref AssignAllGridCoordinates
- face metrics: @ref ComputeFaceMetrics
- cell-centered Jacobian inverse: @ref ComputeCellCenteredJacobianInverse
- orientation check, which refuses a uniformly left-handed grid rather than repairing it: @ref CheckAndFixGridOrientation

Useful geometric helper for BC and flux logic:

- @ref CalculateFaceCenterAndArea

Per-cell geometry for the subgrid models, computed from a cell's face-area vectors and
inverse volume alone, so that it is the same for an identical cell wherever and however
it sits in space:

- @ref ComputeCellDirectionalExtents returns the cell's extent across each face pair (volume
  over that pair's area), from which the `geometric_mean`, `max_edge` and `scotti` filter
  widths are built;
- @ref ComputeCellEdgeVectors returns the covariant edge vector along each grid direction,
  direction and length together, which the Vreman model projects the velocity gradient
  onto.

@ref ComputeCellCharacteristicLengthScale returns the Cartesian components of the cell
diagonal instead. Those change when the same cell is rotated, which is why no subgrid
path uses it; it has no production caller.

@section p22_literature_sec 3. Literature Anchors

Relevant CURVIB references for this code path:

- Borazjani I, Ge L, Sotiropoulos F. "Curvilinear immersed boundary method for simulating fluid structure interaction with complex 3D rigid bodies." *Journal of Computational Physics* 227(16), 7587-7620 (2008). DOI: `10.1016/j.jcp.2008.04.024`.
- Borazjani I, Di Achille P, D'Souza RM, et al. "The functional role of left atrial flow in ventricular filling and flow evolution in the left ventricle." *Annals of Biomedical Engineering* 41(6), 1265-1275 (2013). DOI: `10.1007/s10439-013-0758-9`.

PICurv implementation notes:

- PICurv follows the same curvilinear metric-centered philosophy (precomputed geometric tensors and metric-aware operators).
- PICurv implements the curvilinear formulation these references build on, together with the particle workflows; it does not implement their immersed-boundary or fluid-structure components.

@section p22_practical_sec 4. What This Means For Users

You mainly control CurvIB behavior through:

- grid generation/ingestion choice,
- boundary-condition handlers,
- solver settings that influence projection/Poisson robustness,
- case geometry and decomposition settings.

Even when users never interact with metric tensors directly, they govern stability, pressure correction quality, and particle coupling accuracy.

@section p22_refs_sec 5. Related Pages

- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@subpage 23_Fractional_Step_Method**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
