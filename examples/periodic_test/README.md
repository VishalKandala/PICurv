# Geometric-Periodic Example

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
