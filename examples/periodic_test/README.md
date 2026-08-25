# Periodic Examples

Two families live here: plain geometric-periodic smoke cases, and driven
periodic flows that sustain a bulk flux against wall drag.

## Geometric periodic

`constant_uniform_flow/` is the supported Eulerian geometric-periodic smoke
case. Its three opposite face pairs are declared `PERIODIC/geometric`; no
separate periodic flags are configured.

The runtime derives periodic axes from those BC pairs before creating the
DMDA. Each paired grid surface must match pointwise under one nonzero constant
Cartesian translation, and every periodic axis must have at least four
physical nodes. The startup banner reports the derived axes and validated
translations.

The node minimum applies on every multigrid level, so this example uses enough
cells and only two levels to keep its coarsest grid valid.

Run focused coverage with:

```bash
make unit-periodic
make smoke-periodic
```

Geometric periodicity currently applies to Eulerian fields only. Particle
positions are not wrapped across periodic boundaries.

## Driven periodic

A streamwise-periodic domain has no inlet or outlet, so nothing sustains a bulk
flow. A driven handler measures the volumetric flux, compares it against a
target, and feeds the difference back as a body force.

- `driven_channel/` — plane channel: an exact laminar verification case, DNS at
  `Re_tau = 180` and `395` against Moser, Kim & Mansour (1999), and a
  constant-Smagorinsky LES repeat. See `driven_channel/driven_channel.md`.
- `driven_duct/` — square duct at `Re_b = 4410`, which sustains the
  turbulence-driven secondary flow of the second kind a plane channel lacks.
  See `driven_duct/driven_duct.md`.

Both directories ship a `case_constant_flux.yml` and (where useful) a
`case_initial_flux.yml`, exercising the two driven handlers: one drives to a
flux you prescribe, the other to the flux it measures from the initial
condition. Grids are generated from checked-in `.cfg` files under
`config/grids/`; no `.picgrid` binaries are shipped.

Run focused coverage with:

```bash
make smoke-driven-periodic
```

That asserts the multigrid coarse-solve residual invariant at 4 and 10 ranks
and the `initial_flux` latch-and-restart contract.

Semantics, control law, restart behaviour and known limitations:
`docs/pages/54_Geometric_Periodic_Boundaries.md`, section 5.
