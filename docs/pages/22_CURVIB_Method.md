@page 22_CURVIB_Method CurvIB Method Overview

PICurv follows a curvilinear immersed-boundary (CurvIB-style) workflow where flow variables evolve on structured curvilinear grids and immersed-body effects are applied through dedicated IBM hooks.

@tableofcontents

@section concepts_sec 1. Core Concepts

- Structured curvilinear geometry and metrics are managed per block (`DMDA` + metric vectors/Jacobians).
- Momentum and pressure updates occur in grid-aligned computational space.
- Immersed-boundary corrections are injected conditionally when IBM is active.

@section code_sec 2. Main Code Touchpoints

- Grid dimensions and DM setup: @ref DefineAllGridDimensions, @ref InitializeAllGridDMs, @ref AssignAllGridCoordinates
- Metric construction: @ref ComputeFaceMetrics, @ref ComputeCellCenteredJacobianInverse, @ref CheckAndFixGridOrientation
- Flow-step orchestration: @ref FlowSolver
- Boundary geometry helpers used by setup/BC layers: @ref CalculateFaceCenterAndArea

@section practical_sec 3. Practical Interpretation

For users, this mostly appears as:
- grid choice (`programmatic_c`, `file`, `grid_gen`),
- boundary-condition/physics configuration,
- optional immersed-mode flags in advanced settings.

@section refs_sec 4. References

- Curvilinear immersed-boundary literature surveys (general background)
- PETSc DMDA overview: https://petsc.org/release/manualpages/DMDA/
