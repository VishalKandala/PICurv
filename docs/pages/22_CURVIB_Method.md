@page 22_CURVIB_Method CurvIB Method Overview

@anchor _CURVIB_Method

PICurv solves flow on structured curvilinear grids and applies immersed-boundary-aware logic in metric, pressure, and projection stages.

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
- orientation checks/fixups: @ref CheckAndFixGridOrientation

Useful geometric helper for BC and flux logic:

- @ref CalculateFaceCenterAndArea

@section p22_ibm_sec 3. Immersed-Boundary Role In Current Code

The current branch keeps immersed-boundary hooks in solver paths, while several IBM-specific calls are conditional or currently inactive in default runs.
In practice:

- geometry masking uses `Nvert`/solid markers in key kernels,
- Poisson and projection stencils include boundary-aware logic,
- IBM-specific interpolation hooks remain extension points for fully active IBM runs.

@section p22_literature_sec 4. Literature Anchors

Relevant CURVIB references for this code path:

- Borazjani I, Ge L, Sotiropoulos F. "Curvilinear immersed boundary method for simulating fluid structure interaction with complex 3D rigid bodies." *Journal of Computational Physics* 227(16), 7587-7620 (2008). DOI: `10.1016/j.jcp.2008.04.024`.
- Borazjani I, Di Achille P, D'Souza RM, et al. "The functional role of left atrial flow in ventricular filling and flow evolution in the left ventricle." *Annals of Biomedical Engineering* 41(6), 1265-1275 (2013). DOI: `10.1007/s10439-013-0758-9`.

PICurv implementation notes:

- PICurv follows the same curvilinear metric-centered philosophy (precomputed geometric tensors and metric-aware operators).
- the current branch emphasizes stable production paths for structured curvilinear flow + particle workflows, while preserving immersed-boundary extension hooks for deeper IBM/FSI development.

@section p22_practical_sec 5. What This Means For Users

You mainly control CurvIB behavior through:

- grid generation/ingestion choice,
- boundary-condition handlers,
- solver settings that influence projection/Poisson robustness,
- case geometry and decomposition settings.

Even when users never interact with metric tensors directly, they govern stability, pressure correction quality, and particle coupling accuracy.

@section p22_refs_sec 6. Related Pages

- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@subpage 23_Fractional_Step_Method**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **CurvIB Method Overview** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
