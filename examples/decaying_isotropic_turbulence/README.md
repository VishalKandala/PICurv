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

## Turbulence statistics

`monitor.yml` accumulates two field-statistics windows while the solver runs, so
the turbulence statistics do not depend on how many instantaneous states you
happened to save. Each window keeps its own weighted centered moments of `Ucat`,
`P`, and the LES eddy viscosity `Nu_t`, plus the `Ucat`-`P` co-moment:

- `early_decay` covers `t` in `[0.2, 0.8]`,
- `late_decay` covers `t` in `[1.2, 2.0]`.

Both start after the synthetic initial field has adjusted to the discrete
operators, and both weight each accepted state by the physical time it
represents, so the averages stay correct if the timestep ever varies.

`post.yml` derives them into Reynolds stresses, RMS, turbulent kinetic energy,
and the pressure-velocity flux:

```bash
picurv run --post-process --run-dir runs/<run_id> --post dit_case/post.yml
```

That writes, per window, a `.vts` per processed step under
`visualization/dit/` and one convergence-history `.csv` recording sample count,
accumulated weight, represented time, mask coverage, and the domain-mean TKE.

### What to check

Isotropy makes this case self-checking. Within a converged window:

- the three normal Reynolds stresses `<u'u'>`, `<v'v'>`, `<w'w'>` should agree to
  within the sampling error of the window, and
- the shear components `<u'v'>`, `<u'w'>`, `<v'w'>` should be small compared with
  the normal components.

A window that has not yet gathered enough states shows this as a spread in the
normal stresses rather than as an obvious error, so read the sample count in the
CSV before drawing a conclusion from a short run. Across windows, `late_decay`
should carry the lower TKE; the ratio between the two window means is the decay
measured over the interval separating them.

Because the box is triply periodic, no plane is counted twice and every cell is
sampled, so the mask coverage column stays at `1.00`. A coverage below one here
would indicate a masking or target-plan problem rather than a physical one.

Both windows are bounded, so they complete before the run ends and their saved
state stops changing. Resuming the run with `--continue` carries completed
windows forward untouched and keeps accumulating any window still open.
