@page 32_Analytical_Solutions Analytical Fields

@anchor _Analytical_Solutions

With `operation_mode.eulerian_field_source: analytical` nothing is solved. The Eulerian field is the
case's `properties.initial_conditions.source`, the definition an initial condition and an inlet take
(@ref p33_source), evaluated at each step's time.

@tableofcontents

@section p32_activation_sec 1. Activation Path

- `solver.yml -> operation_mode.eulerian_field_source: analytical` -> `-euler_field_source analytical`.
- `case.yml -> properties.initial_conditions.source` -> `-ic_source_*`, as for any initial condition.
- @ref InitializeEulerianState and the run loop call @ref AnalyticalSolutionEngine.

Each evaluation sets the fluxes on every face, including the boundary faces, `Ubcs` on every boundary
face, and at the cell centres `Ucat` and `P`, then finalizes the dummy cells as a solved flow is. As on
a time-dependent inlet, a source that does not read `t` is evaluated once.

@section p32_types_sec 2. What Can Be Prescribed

Anything an expression can write, in `cartesian` or `curvilinear` components, with these limits:

- a generator or file term is refused, because the field is re-evaluated, and so is `project`;
- `grid.mode` is `programmatic_c` or `file`;
- a `curvilinear` source derives `Ucat` from its fluxes; a `cartesian` one sets it exactly at the cell
  centres;
- a `pressure` expression [Pa] sets `P` at every evaluation. Without it `P` is zero.

The domain is the case's grid: no field assigns its own.

Zero flow:

```yaml
properties:
  initial_conditions:
    source: {components: cartesian, expression: {u: 0.0, v: 0.0, w: 0.0}}
```

A uniform stream, with its speed as a constant a study can sweep
(`case.properties.initial_conditions.source.constants.U`):

```yaml
properties:
  initial_conditions:
    source: {components: cartesian, constants: {U: 0.1}, expression: {u: U}}
```

@section p32_tgv_sec 3. Taylor-Green Vortex

```yaml
properties:
  initial_conditions:
    source:
      components: cartesian
      constants: {nu: 0.01}
      expression:
        u: "sin(x)*cos(y)*cos(z)*exp(-2*nu*t)"
        v: "-cos(x)*sin(y)*cos(z)*exp(-2*nu*t)"
      pressure: "0.25*(cos(2*x) + cos(2*y))*exp(-4*nu*t)"
```

with `L = U = rho = 1` and `nu = 1/Re`. `examples/interpolation_test` runs it on
`[0, 2 pi]^2 x [0, 0.2 pi]`, which its grid states.

@section p32_particles_sec 4. Particle Consistency

When the source is `cartesian`, reads only `x y z t` and varies in space, the interpolated particle
velocities are compared with it at the particle positions, and `interpolation_error.csv` is written.
A uniform field has no interpolation error to measure, and a source in the grid-line fractions cannot
be evaluated at an arbitrary point. `FaceFlux` interpolation reads the fluxes every analytical field
sets.

@section p32_evidence_sec 5. Evidence

Production exercised in `examples/interpolation_test`, `examples/drift_uniform_flow`,
`examples/drift_diffusivity_gradient` and `examples/brownian_motion`. Measured when these fields were
selected by name, and reproduced from their sources since: the Brownian cloud of
`brownian-msd-2026-09-18` in zero flow, the face fluxes on a curved file grid of
`curvilinear-gcl-2026-09-18` in a uniform stream, and the staged velocity and pressure of
`solution-monitoring-tgv3d-2026-09-18` for the vortex. On reduced grids over three steps, the source
forms match those fields to 2.2e-16 in `Ucat`, and for the vortex exactly in `Ucat`, `P`, particle
positions and velocities and the interpolation error. A vortex's `FaceFlux` interpolation and its
fluxes have not been measured against an exact solution.

@section p32_extension_sec 6. Adding A Field

Write its expressions; no code is needed. A field an expression cannot write belongs in a generator
as an initial condition, which is solved from, not prescribed.
