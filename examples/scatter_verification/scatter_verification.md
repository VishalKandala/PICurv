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
- the runtime writes `logs/scatter_metrics.csv`

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
be paired with `ZERO_FLOW`, `UNIFORM_FLOW`, or `TGV3D`.

## What `scatter_metrics.csv` Measures

The runtime diagnostic writes one row per output point to:

- `<run_dir>/logs/scatter_metrics.csv`

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

The reusable verification feature supports three scalar profiles:

- `CONSTANT`: useful for conservation / preservation checks
- `LINEAR_X`: useful for affine and gradient-sensitive checks
- `SIN_PRODUCT`: useful as a smooth nontrivial verification profile for grid
  sensitivity studies

This example uses `SIN_PRODUCT` because it creates a meaningful grid-size trend
for the cell-centered scatter metric. If you want a more directly physical
sanity check, `CONSTANT` and `LINEAR_X` are already available through the same
feature.

## Files In This Example

- `scatter_verification.yml`: deposition-only case definition
- `Analytical-Zero-Verification-Scalar.yml`: analytical zero-flow solver plus
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
./bin/picurv run --solve --post-process -n 4   --case  examples/scatter_verification/scatter_verification.yml   --solver examples/scatter_verification/Analytical-Zero-Verification-Scalar.yml   --monitor examples/scatter_verification/Standard_Output.yml   --post examples/scatter_verification/scatter_verification_analysis.yml
```

Primary output:

- `<run_dir>/logs/scatter_metrics.csv`

Optional visualization output:

- `<run_dir>/viz/Field_*.vts`
- `<run_dir>/viz/Particle_*.vtp`

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

## Future Verification And Diagnostic Pathways Enabled

This feature is intentionally broader than one study. Prescribed scalar truth
injection plus runtime scatter metrics opens up several future pathways:

- static deposition verification with `ZERO_FLOW`
- moving-cloud verification with `UNIFORM_FLOW`
- coupled flow-plus-scatter verification with `TGV3D`
- conservation diagnostics for deposited scalar fields
- grid and curvilinear-geometry sensitivity studies
- future scalar-transport or scalar-mixing verification reuse without adding a
  one-off test hook each time

That is the main reason the analytical scalar truth lives in
`AnalyticalSolutions` and the runtime override lives in `verification_sources`:
this is meant to be reusable verification infrastructure, not a one-off study
patch.
