@page 21_Methods_Overview Methods and Models Overview

This section gives concise method-level summaries of the major numerical and modeling components used in PICurv.

@tableofcontents

@section map_sec 1. Method Map

- **@subpage 22_CURVIB_Method**: Curvilinear immersed-boundary framing and data layout.
- **@subpage 23_Fractional_Step_Method**: Predictor-projection incompressible flow update.
- **@subpage 24_Dual_Time_Picard_RK4**: Dual-time momentum solver with RK4 smoothing.
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**: Pressure solve stack and PETSc controls.
- **@subpage 26_Walking_Search_Method**: Particle cell location and migration logic.
- **@subpage 27_Trilinear_Interpolation_and_Projection**: Grid-particle interpolation and scatter averaging.
- **@subpage 28_IEM_and_Statistical_Averaging**: Particle mixing model and flow/statistics averaging paths.
- **@subpage 31_Momentum_Solvers**: Named momentum solver implementations and dispatch points.
- **@subpage 32_Analytical_Solutions**: Analytical flow modes and geometry policy.
- **@subpage 33_Initial_Conditions**: Eulerian and particle IC options.
- **@subpage 34_Particle_Model_Overview**: End-to-end particle model/coupling pipeline.

@section scope_sec 2. Scope Note

These are implementation-aligned overviews, not derivation-heavy treatises.
Each page includes code touchpoints and a few external references for deeper reading.
