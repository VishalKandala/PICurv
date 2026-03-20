# Diffusivity-Gradient Drift Verification

## Purpose

Verify the deterministic `diffusivitygradient * dt` term in
`UpdateParticlePosition` using a full analytical-runtime case. The carrier flow
is `ZERO_FLOW`; the only deterministic motion comes from the verification-only
linear diffusivity profile.

## Analytical Reference

The verification override imposes

- `Gamma(x) = gamma0 + a x`

with `gamma0 = 0.002` and `a = 0.001`, so

- `grad(Gamma) = (a, 0, 0)`

The cloud centre of mass should therefore drift as

- `com_x(t) = x0 + a * t`

while `com_y` and `com_z` remain near zero. Brownian spreading is still present,
but it is centered around that drifting mean.

## Run

```bash
./bin/picurv run --solve --post-process -n 4 \
  --case  examples/drift_diffusivity_gradient/case.yml \
  --solver examples/drift_diffusivity_gradient/solver.yml \
  --monitor examples/drift_diffusivity_gradient/monitor.yml \
  --post examples/drift_diffusivity_gradient/post.yml
```

Primary output:

- `<run_dir>/DiffusivityGradientDrift_msd.csv`

## Verification

Fit `com_x` vs `t` and confirm the slope matches `0.001`. `com_y` and `com_z`
should stay near zero. This exercises:

- verification-source injection into Eulerian diffusivity
- `ComputeEulerianDiffusivityGradient`
- interpolation of `DiffusivityGradient` to particles
- deterministic drift inside `UpdateParticlePosition`
