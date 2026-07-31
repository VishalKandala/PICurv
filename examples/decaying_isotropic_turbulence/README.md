# Decaying isotropic turbulence

This is a reproducible 64-cubed-cell LES decay benchmark in a triply periodic
`[0,2*pi]^3` box. It is an LES demonstration, not a DNS reference.

Create and validate it with:

```bash
picurv init decaying_isotropic_turbulence --dest dit_case
picurv validate --case dit_case/case.yml --solver dit_case/solver.yml --monitor dit_case/monitor.yml --post dit_case/post.yml
picurv precompute --case dit_case/case.yml --output-dir /tmp/dit-precompute
picurv run --solve -n 8 --case dit_case/case.yml --solver dit_case/solver.yml --monitor dit_case/monitor.yml
```

This example configures the reusable `spectral_random_velocity` provider in
`generators/ic.gen`; the production implementation contains no case-specific
DIT generator. Its `k4_exponential` envelope gives
`E(k) proportional to k^4 exp[-2(k/k0)^2]`, and `component_rms: 1` means
`sqrt(<|u'|^2>/3) = 1`. The viscosity used for dissipation and Reynolds-scale
diagnostics comes from `properties.fluid` through the same nondimensional
scaling used by the solver, rather than from IC parameters.

The normalization target describes the staged `Ucat` field. The diagnostic
summary also predicts the cell velocity after PICurv's native face/cell
reconstruction in `picurv_reconstructed_kinetic_energy`. Treat that prediction
as a planning and consistency-check value: record the solver's actual step-zero
kinetic-energy output as the measured `K(0)` used to normalize and interpret the
decay. Use the staged value only when checking the provider contract.

The example selects `projection.operator: picurv_discrete`, so PICurv's
uniform-grid centered-divergence symbol `sin(k dx)/dx` is roundoff-level. The
continuum divergence is reported but is not expected to vanish. A continuum
projection can be selected instead; either single-operator projector retains
two transverse polarizations. The former common-nullspace construction usually
retained only one and is intentionally not used.

Preparation writes `diagnostics/initial_condition_summary.json` and
`diagnostics/initial_condition_spectrum.csv` under the selected output/run
directory. For a quick smoke run, copy the case, set `im/jm/km` to 16,
`k_cut` to 5, multigrid levels to 2, and `total_steps` to 1 before running the
full 64-cubed case.
