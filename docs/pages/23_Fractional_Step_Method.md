@page 23_Fractional_Step_Method Fractional-Step (Projection) Method

PICurv uses a projection-style incompressible update: momentum predictor, pressure correction, velocity projection.

@tableofcontents

@section equations_sec 1. Discrete Method Skeleton

Given velocity \f$\mathbf{u}^n\f$, the solver computes a provisional field \f$\mathbf{u}^*\f$ from momentum terms, then solves pressure correction:

\f[
\nabla^2 \phi = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*,
\f]

followed by:

\f[
\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\nabla\phi,
\qquad
p^{n+1}=p^n+\phi.
\f]

In curvilinear form, gradient components use metric terms (`ICsi/IEta/IZet`, etc.) and Jacobian inverses in the projection kernel.

@section implementation_sec 2. Code Path In PICurv

- stage orchestrator: @ref FlowSolver
- momentum residual assembly: @ref ComputeRHS
- pressure solve (MG/KSP): @ref PoissonSolver_MG
- pressure update: @ref UpdatePressure
- velocity projection: @ref Projection

`@ref Projection` explicitly uses boundary-aware stencils and periodic-edge correction logic before final ghost updates and `Contra2Cart` conversion.

@section boundary_sec 3. Boundary and Geometry Handling

Key implementation details:

- one-sided derivative stencils are used when centered stencils touch solids (`Nvert` checks),
- periodic edges are corrected in dedicated boundary loops,
- driven-periodic channel cases call the flux-profile corrector in projection path.

These are part of why projection behavior in PICurv differs from simplified textbook pseudocode.

@section runtime_sec 4. Diagnostics To Watch

For projection health, monitor:

- continuity/divergence diagnostics after correction,
- Poisson KSP convergence history logs,
- step-to-step stability of pressure and corrected velocity.

Persistent divergence or noisy pressure correction usually indicates a Poisson tuning, grid-quality, or BC consistency issue.

@section refs_sec 5. Related Pages

- **@subpage 24_Dual_Time_Picard_RK4**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
- **@subpage 39_Common_Fatal_Errors**
