# Turbulent plane channel (wall-modelled LES)

A driven periodic channel at `Re_tau ~ 1000`. Where
`decaying_isotropic_turbulence` exercises the subgrid closure in a flow with no
wall, this case exercises everything that only exists because there is one: the
wall model and its three laws, wall-normal-inhomogeneous coefficient averaging,
the driven periodic pair, and the near-wall diagnostics.

## What it is set up to be

| Quantity | Value | Why |
|---|---|---|
| `Re = U_b delta / nu` | 20000 | with `Re_b = 2 U_b delta / nu = 40000` |
| `Re_tau` (Dean) | ~1016 | `C_f = 0.073 Re_b^(-1/4)` gives `u_tau/U_b = 0.0508` |
| Domain | `pi x 2 x 2pi` | spanwise x wall-normal x streamwise, `delta = 1` |
| Grid | `33 x 21 x 49` nodes | uniform; `dz+ = 100`, first cell `y+ = 51`, `dx+ = 133` |
| Driving | `constant_flux` | `target_flux = 6.2832`, the volumetric flux at `U_b = 1` |

The wall-normal grid is uniform on purpose. A wall model exists so the boundary
layer need not be resolved, so its first cell belongs in the logarithmic region,
not packed against the wall. That is also why this mesh cannot be reused with the
wall model switched off: a wall-resolved LES at this Reynolds number needs a first
cell near `y+ = 1`, roughly fifty times finer, and clustered at both walls. The
programmatic grid generator applies a one-sided geometric stretch, so a
symmetrically clustered channel needs a grid file rather than
`mode: programmatic_c`.

## What it is for

Toggle one thing at a time and compare against the literature:

- `wall_function.model` between `log_law`, `werner`, and `cabot`;
- `wall_function.enabled` off, to see an under-resolved run on a mesh built for a
  wall model - informative, not a wall-resolved reference;
- `les.averaging.mode` between `homogeneous` and `local`, and
  `les.clipping.mode` between `clamp` and `none`, to see the coefficient go from
  a wall-normal profile to a signed field;
- `les.gradient_model.enabled` for the mixed model.

## What to compare against

- **The log law.** Above `y+ ~ 30` the mean profile should follow
  `u+ = (1/0.41) ln(y+) + 5.2`. This is the first thing to check and the one that
  fails most visibly if the wall stress is wrong.
- **Dean's correlation.** `C_f = 0.073 Re_b^(-1/4)`, here `5.16e-3`, equivalently
  `u_tau/U_b = 0.0508`. `u_tau_mean` in `<run.runtime_logs>/wall_model.csv` reports
  the model's own value directly.
- **Lee & Moser (2015)** channel DNS at `Re_tau = 1000` for mean and Reynolds
  stress profiles. The `developed` statistics window collects the first and second
  moments of `Ucat` that these compare to.

Read `y_plus_mean` and `nu_wall_over_nu_mean` in `wall_model.csv` first. The first
says whether the mesh puts the first cell where the selected law is valid; the
second says whether the modelled stress is reaching the momentum equation, and
should be near `y+/u+ - 1`.

## Known gap: the initial condition

**This case as shipped will not become turbulent.** It starts from the laminar
Poiseuille profile that carries the target flux, which is the correct mean and the
wrong state - nothing in a well-behaved implicit solve makes a one-dimensional
laminar channel three-dimensional. A perturbed field is needed, and this
repository has no generator that produces one for a wall-bounded domain:
`generators/ic.gen` builds a solenoidal field by a spectral projection that
assumes periodicity in all three directions, which a channel does not have.

Until such a generator exists, use this case either as a configuration and
plumbing exercise, which is what it is verified as, or restart it from a perturbed
field produced elsewhere. Treat any statistics gathered from the shipped initial
condition as laminar.

## Cost

Measured at 3.5 s/step on four ranks for the shipped grid. Converged statistics
need order fifty flow-through times, `t ~ 300` or 15000 steps, which is hours -
a cluster job, not a session. `total_steps` ships at 1000, a spin-up length.

```bash
picurv init turbulent_channel --dest my_channel
cd my_channel
picurv validate --case case.yml --solver solver.yml --monitor monitor.yml
picurv run --solve --case case.yml --solver solver.yml --monitor monitor.yml -n 4
picurv summarize --run-dir runs/<run> --plot wall_model.y_plus_mean
```
