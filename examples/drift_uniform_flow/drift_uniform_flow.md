# Uniform Flow Drift Verification

## Purpose

Verify the deterministic particle advection term in `UpdateParticlePosition` by
advecting a point-source cloud in a constant analytical carrier flow. This case
targets the `vel * dt` contribution only.

## Analytical Reference

For a constant imposed velocity `U = (u0, v0, w0)`, the particle cloud centre
of mass should move linearly:

- `com_x(t) = x0 + u0 * t`
- `com_y(t) = y0 + v0 * t`
- `com_z(t) = z0 + w0 * t`

This profile uses `U = (0.1, 0.0, 0.0)` and suppresses Brownian motion by
setting an effectively infinite Schmidt number.

## Run

```bash
./bin/picurv run --solve --post-process -n 4 \
  --case  examples/drift_uniform_flow/drift_uniform_flow.yml \
  --solver examples/drift_uniform_flow/Analytical-UniformFlow.yml \
  --monitor examples/drift_uniform_flow/Standard_Output.yml \
  --post examples/drift_uniform_flow/uniform_flow_analysis.yml
```

Primary output:

- `<run_dir>/UniformFlowDrift_msd.csv`

## Verification

Fit `com_x` vs `t` and confirm the slope is `0.1`. `com_y` and `com_z` should
remain near zero. `MSD_total` should remain negligible because diffusivity is
suppressed for this verification case.
