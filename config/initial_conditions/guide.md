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
