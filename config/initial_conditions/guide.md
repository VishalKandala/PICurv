# Initial-Condition Generation

`generators/ic.gen` generates one complete PETSc binary `Ucat` or `Ucont` vector
from expressions evaluated on a staged nondimensional `PICGRID`. The staged grid
may come from `file`, `grid_gen`, or single-block `programmatic_c` cases.

```ini
[expression]
u = sin(2*pi*x) * cos(2*pi*y)
v = -cos(2*pi*x) * sin(2*pi*y)
w = 0.0
```

Select it from `case.yml`:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: ic_gen
    params:
      field: Ucat
      config_file: config/initial_conditions/expression.cfg
```

For `Ucat`, `u`, `v`, and `w` are Cartesian velocity components evaluated at
actual cell centers, including extrapolated dummy-cell centers. For `Ucont`,
define `u_xi`, `u_eta`, and `u_zeta`; each component is evaluated at its own
geometric face centers and represents a contravariant flux.

Available values are `x`, `y`, `z`, `xi`, `eta`, `zeta`, `i`, `j`, `k`, and
`pi`. Supported functions are `sin`, `cos`, `tan`, `exp`, `sqrt`, `abs`,
`minimum`, `maximum`, and `where`.

Use `picurv precompute --case case.yml` to inspect the generated and staged
PETSc vectors before launching a solve.

## Spectral random velocity

The same `ic.gen` entry point provides a composable spectral generator:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: spectral_random_velocity
    params:
      field: Ucat
      seed: 12345
      random: {distribution: gaussian, mean: [0.0, 0.0, 0.0]}
      spectrum: {type: k4_exponential, k0: 4.0, k_cut: 20.0}
      projection: {type: solenoidal, operator: picurv_discrete}
      normalization: {type: component_rms, target: 1.0}
      remove_mean: true
```

Defaults are `field: Ucat`, `seed: 12345`, Gaussian sampling with zero requested
mean, `spectrum.type: white`, `projection.type: none`,
`normalization.type: none`, and `remove_mean: true`. Selectors are strict; no
hyphenated or case-specific compatibility aliases are accepted.

`spectrum.type` may be `white` or `k4_exponential`. White applies no spectral
shaping. The latter gives `E(k) proportional to k^4 exp[-2(k/k0)^2]`, requires
positive `k0` and `k_cut`, and distributes each target shell's energy over the
finite number of active lattice modes in that shell. Its cutoff may not exceed
the grid's two-thirds-resolution limit. `normalization.type: component_rms` sets the fluctuation measure
`sqrt(<|u'|^2>/3)`; `none` leaves the filtered realization unscaled.

The target and `resolved_kinetic_energy` describe the staged cell-centered
`Ucat` field. On a uniform Cartesian grid, PICurv's native
`Ucat -> Cart2Contra -> Contra2Cart` reconstruction multiplies component
`alpha` by `cos^2(k_alpha*Delta_alpha/2)`. The summary therefore also reports
`picurv_reconstructed_component_rms` and
`picurv_reconstructed_kinetic_energy`; keep the energetic spectrum well below
the grid cutoff when the reconstructed and staged energy levels must be close.
These reconstructed diagnostics are predictions for planning and consistency
checks. For decay studies, the solver's actual step-zero output is the
authoritative measured initial energy.

Set `projection.type: none` for an unconstrained field. A solenoidal projection
must explicitly select one operator: `continuum` projects against `k`, while
`picurv_discrete` projects against `sin(k dx)/dx`, the uniform-grid symbol
associated with PICurv's centered flux divergence. Selecting one operator keeps
the usual two transverse polarizations. Enforcing both simultaneously would
normally retain only the one-dimensional intersection of their nullspaces.
Consequently only the selected divergence measure is expected at roundoff; the
summary reports both measures.
An unprojected random velocity is generally compressible and may produce a
large initial pressure-projection transient; it is intended for explicit
diagnostic and perturbation experiments.

Generation order is sampling, spectral shaping, projection, optional sampled
mean removal, normalization, then addition of `random.mean`. Thus
`remove_mean: true` removes the finite-sample fluctuation mean before
normalization, while `false` retains it. The configured mean is added afterward
and is never removed accidentally. The generator requires a fresh, single-block, 3-D,
uniform Cartesian grid with geometric-periodic faces and `Ucat` output.
Preparation writes `diagnostics/initial_condition_summary.json` and
`diagnostics/initial_condition_spectrum.csv`.
