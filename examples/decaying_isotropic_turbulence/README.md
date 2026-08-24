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

## Energy spectra

The case seeds a prescribed `k4_exponential` envelope, so measuring the realized
spectrum back is the most direct statement of what the solver did to it.
`post.yml` requests one:

```bash
picurv run --post-process --only spectra --run-dir runs/<run_id> --post dit_case/post.yml
picurv summarize --run-dir runs/<run_id> --plot-spectrum
```

`--only spectra` skips the field post-processor, which makes re-measuring cheap.
The report plot draws up to six evenly spaced measured states (always including
the first and last) over `diagnostics/initial_condition_spectrum.csv`; the CSV
still retains every processed state. The per-step scalars
(`spectra.resolved_kinetic_energy`, `spectra.integral_length_scale`,
`spectra.taylor_microscale`, `spectra.dissipation_over_viscosity`) are ordinary
series that `picurv summarize --plot` draws against physical time.

### What to expect

- **The staged initial condition and step zero do not coincide.** The solver
  round-trips `Ucat` through the contravariant fluxes, which applies
  `cos^2(k dx/2)` per component. The gap is concentrated at high `k` and matches
  `picurv_reconstructed_kinetic_energy` in the initial-condition summary; treat
  that value, not the staged one, as `K(0)`.
- **Modes above `k_cut` fill in.** The initial field carries no energy there at
  all, so anything that appears is the nonlinear term building the cascade. This
  is the clearest sign the measurement is physical.
- **Dissipation rises before it falls.** A synthetic field starts with no
  small-scale content, so the cascade spins up first: `dissipation_over_viscosity`
  increases and the Taylor microscale shrinks. Self-similar decay, where both
  reverse, comes afterwards. Do not fit a decay exponent through the spin-up.
- **`parseval_residual` should stay at round-off.** It is a self-check: summed
  shell energy must equal the resolved kinetic energy of the transformed field.

Note the two-thirds rule when changing resolution: `k_cut` must stay below
`(2/3) * pi / dx`, so a 32-cubed variant of this case needs `k_cut` at or below
about 10.6.

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
