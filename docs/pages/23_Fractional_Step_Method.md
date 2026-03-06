@page 23_Fractional_Step_Method Fractional-Step (Projection) Method

@anchor _Fractional_Step_Method

PICurv uses a projection-style incompressible update: momentum predictor, pressure correction, velocity projection.

@tableofcontents

@section p23_equations_sec 1. Discrete Method Skeleton

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

@section p23_implementation_sec 2. Code Path In PICurv

- stage orchestrator: @ref FlowSolver
- momentum residual assembly: @ref ComputeRHS
- pressure solve (MG/KSP): @ref PoissonSolver_MG
- pressure update: @ref UpdatePressure
- velocity projection: @ref Projection

`@ref Projection` explicitly uses boundary-aware stencils and periodic-edge correction logic before final ghost updates and `Contra2Cart` conversion.

@section p23_boundary_sec 3. Boundary and Geometry Handling

Key implementation details:

- one-sided derivative stencils are used when centered stencils touch solids (`Nvert` checks),
- periodic edges are corrected in dedicated boundary loops,
- driven-periodic channel cases call the flux-profile corrector in projection path.

These are part of why projection behavior in PICurv differs from simplified textbook pseudocode.

@section p23_runtime_sec 4. Diagnostics To Watch

For projection health, monitor:

- continuity/divergence diagnostics after correction,
- Poisson KSP convergence history logs,
- step-to-step stability of pressure and corrected velocity.

Persistent divergence or noisy pressure correction usually indicates a Poisson tuning, grid-quality, or BC consistency issue.

@section p23_refs_sec 5. Related Pages

- **@subpage 24_Dual_Time_Picard_RK4**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
- **@subpage 39_Common_Fatal_Errors**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Fractional-Step (Projection) Method** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

