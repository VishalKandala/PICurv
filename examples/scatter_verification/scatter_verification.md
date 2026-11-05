# Scatter Verification

## Purpose

This example verifies the particle-to-grid scatter path for the particle scalar
`Psi` using a prescribed analytical truth field. The run is intentionally
**deposition-only**:

- particles are initialized and settled
- `verification.sources.scalar` overwrites particle `Psi`
- the existing runtime `Psi -> user->Psi` scatter path is reused
- the scattered Eulerian `Psi` field is compared against analytical truth at
  **physical cell centers**
- the runtime writes `<run.analysis.metrics>/scatter_metrics.csv`

This example exists because scatter verification needs an artificial particle
truth field. Unlike the interpolation test, there is no ordinary end-to-end
runtime path that naturally produces an exact analytical `Psi` reference.

## Why This Uses The Verification Path

The verification pathway is used here for one reason only: we need to prescribe
particle `Psi` deliberately.

Normal production behavior is preserved:

- ordinary `Psi` model evolution in `ParticlePhysics` is untouched
- the verification override is inactive unless `verification.sources.scalar` is
  configured
- the existing scatter operator is still the operator under test
- no `statistics_pipeline` extension is used for the core metric

When scalar verification is active, only the `Psi` update is bypassed. Carrier
flow selection remains independent, so the same scalar truth machinery can later
be paired with zero flow, a uniform stream, or the Taylor-Green vortex.

## What `scatter_metrics.csv` Measures

The runtime diagnostic writes one row per output point to:

- `<run.analysis.metrics>/scatter_metrics.csv`

Columns include:

- `relative_L2_error`
- `Linf_error`
- `occupancy_fraction`
- `mean_particles_per_occupied_cell`
- `particle_integral`
- `grid_integral`
- `conservation_error_abs`

Reference definition used here:

- analytical truth is evaluated at each **physical cell center** from `lCent`
- the scattered Eulerian `Psi` field is treated as the solver's cell-centered
  stored field
- the metric therefore checks how well the deposited cell-centered field agrees
  with the analytical cell-centered truth the solver would use

This is deliberate. The current verification feature is about the solver's
cell-centered scalar representation, not about reconstructing an exact cell
integral.

## Profiles And Physical Meaning

The scalar is one expression in `x y z` and `t`. Useful choices:

- a constant: conservation and preservation checks
- linear in `x`: affine and gradient-sensitive checks
- `sin(pi*x)*sin(pi*y)*sin(pi*z)`: a smooth nontrivial field for grid sensitivity studies

This example uses the sine product because it creates a meaningful grid-size trend for the
cell-centered scatter metric.

## Files In This Example

- `scatter_verification.yml`: deposition-only case definition
- `Analytical-Verification-Scalar.yml`: analytical solver plus
  scalar verification source
- `Standard_Output.yml`: monitor profile that keeps `scatter_metrics.csv`
  quiet on console unless `LOG_SCATTER_METRICS` is whitelisted at an enabled log
  level
- `scatter_verification_analysis.yml`: lightweight optional visualization recipe
- `fixed_total_particles_grid_study.yml`: vary grid size at fixed total
  particle count
- `fixed_ppc_grid_study.yml`: vary grid size while keeping particles per cell
  approximately constant

## Run

```bash
./bin/picurv run --solve --post-process -n 4   --case  examples/scatter_verification/scatter_verification.yml   --solver examples/scatter_verification/Analytical-Verification-Scalar.yml   --monitor examples/scatter_verification/Standard_Output.yml   --post examples/scatter_verification/scatter_verification_analysis.yml
```

Primary output:

- `<run.analysis.metrics>/scatter_metrics.csv`

Optional visualization output:

- `<run.visualization>/Field_*.vts`
- `<run.visualization>/Particle_*.vtp`

## Studies

### Fixed Total Particles

`fixed_total_particles_grid_study.yml` holds the total particle count fixed at
32768 while sweeping `8^3`, `16^3`, `32^3`, and `64^3` grids.

Interpretation:

- occupancy usually drops as the grid refines
- `mean_particles_per_occupied_cell` drops with refinement
- `relative_L2_error` captures the practical loss of scatter fidelity when the
  same particle budget is spread over more cells

### Fixed Particles Per Cell

`fixed_ppc_grid_study.yml` keeps approximately 8 particles per cell while
sweeping the same grid sizes.

Interpretation:

- occupancy should remain comparatively healthy
- this isolates how the grid itself changes the cell-centered scatter metric
- it also exercises the new general `parameter_sets` study capability, which is
  useful any time multiple config keys must move together instead of forming the full
  cross-product

Stage either study without submitting jobs:

```bash
./bin/picurv sweep   --study examples/scatter_verification/fixed_total_particles_grid_study.yml   --cluster config/schedulers/slurm_default.yml   --no-submit
```

```bash
./bin/picurv sweep   --study examples/scatter_verification/fixed_ppc_grid_study.yml   --cluster config/schedulers/slurm_default.yml   --no-submit
```

## Other Checks The Same Pathway Supports

Prescribed scalar truth injection plus runtime scatter metrics is not specific to
this study. The same configuration supports:

- static deposition verification with zero flow
- moving-cloud verification with a uniform stream
- coupled flow-plus-scatter verification with the Taylor-Green vortex
- conservation diagnostics for deposited scalar fields
- grid and curvilinear-geometry sensitivity studies

That is why the analytical scalar truth lives in `AnalyticalSolutions` and the
runtime override lives in `verification_sources`, rather than in a study-specific
hook.
