# Temporal Mixing Layer (Colucci et al. 1998)

The two-dimensional, non-reacting temporal mixing layer of Colucci, Jaberi, Givi and
Pope, *Filtered density function for large eddy simulation of turbulent reacting
flows*, Phys. Fluids 10:499 (1998): the benchmark for the particle filtered density
function (FDF) of a conserved scalar.

Two streams at `+1` and `-1` roll up into two Kelvin-Helmholtz vortices that pair once.
The scalar `A` (PICurv's `Psi`) is 1 in the upper stream and 0 in the lower; Monte Carlo
particles carry it, mix it through IEM, and diffuse with `Gamma + Gamma_t`. At Re = 500 in
2D the flow is transitional, not turbulent: the benchmark tests the FDF machinery
against filtered DNS on an under-resolved grid, not a turbulence closure.

## Run

```bash
./bin/picurv init temporal_mixing_layer --dest my_tml
cd my_tml
picurv run --solve --post-process -n 4 --case config/case.yml --solver config/solver.yml \
  --monitor config/monitor.yml --post config/post.yml
```

`solver_newton_krylov.yml`, beside `solver.yml` in the initialized case, runs the same
case with the Newton Krylov momentum solver. Measured on 4 ranks (one workstation, 2026-10-22): about 11 s per step for both
solvers, so the 100 steps to t = 50 take about 18 minutes; particle location and
interpolation are half of it.

## How the paper maps onto PICurv

| Colucci et al. | Here |
| --- | --- |
| units: U_r, L_r with delta_v0 / L_r = 2.83 | `velocity_ref = length_ref = 1` |
| Re = 500 on 2 U_r and delta_v0 | `viscosity: 0.01132` |
| x periodic, 2 most-unstable wavelengths | k periodic, `Lx = 39.97` (`<repo>/config/grids/temporal_mixing_layer_2d.cfg`) |
| mirror symmetry at the top and bottom | `-Eta`/`+Eta` `WALL` `free_slip`; particles reflect |
| LES grid 37 x 49, equal spacing | 36 x 48 cells, `Delta = 1.110`, plus 4 periodic cells in i for 2D |
| Smagorinsky `Cs = 0.014`, `Delta_G = 2 Delta` | `constant_smagorinsky`, `constant_cs: 0.0396` |
| Sc = 1, Sc_t = 0.7 | `schmidt_number`, `turbulent_schmidt_number` |
| IEM `C_Omega = 3` | `iem_constant: 0.75` |
| ensemble domain `Delta_E = Delta` | one cell |
| NPG = 50 | `count: 345600` (50 per cubic cell) |
| CFL 0.5 | `dt_physical: 0.5` (0.45) |
| tanh layer + most unstable mode + subharmonic | `stability_mode_velocity`, harmonics 2 and 1 at 1% |

The initial-condition summary's growth-rate scan confirms the box: harmonic 2 grows
fastest (0.133). The scalar starts as an error function with the velocity layer's
thickness. Initial details the paper takes from Metcalfe et al. (1988) are not in the
paper itself: the 1% mode amplitudes, the 2% planar noise and the scalar thickness are
choices here, not the paper's values.

## Outputs and the paper's figures

`post.yml` writes, every t = 2, cross-stream profiles averaged over the homogeneous i and
k directions, and `profile_series` histories of them:

- `tml_profile_fdf_<step>.csv`: `Psi` mean (the paper's Reynolds-averaged `<A>_L`,
  Fig. 2a), `point_variance` (the generalized subgrid variance sigma_A averaged cell by
  cell, Figs. 2b-d, 3), and the pooled `rms`, which also counts the spread of the
  filtered values and is not the paper's quantity.
- `tml_profile_fdf_history.csv`: `Psi_point_variance.integral`, the integrated subgrid
  variance of Fig. 10.
- `tml_profile_resolved_history.csv`: `Ucat_mean_z.max_gradient_thickness`, the vorticity
  thickness delta_v of Figs. 9 and 13.

The thickness `max_gradient_thickness` takes the slope as a central difference over two
cells, `(U(y+Delta) - U(y-Delta)) / 2 Delta`, which is the exact slope of the profile filtered
with a top-hat of width `2 Delta`, the paper's `Delta_G`. That is why the paper's curves start
near 3.4 rather than at the unfiltered 2.83: `2 / (2 tanh(1.11/1.415) / 2.22) = 3.38`. The grid
is shifted by `Delta/2` in y so that a cell centre lies on the layer centre, where that slope
peaks, as a grid point does in the paper; with a face there the reading would be 3.70.
Digitized curves of Figs. 2, 3, 9(a), 10 and 13 are kept with the reference data, outside
the repository.

## Limitations

- Non-reacting only: PICurv has no particle reaction source, so the A + B -> P figures are
  out of reach.
- Smagorinsky only: the MKEV closure is not implemented.
- The paper's `Delta_E = 2 Delta` and `Delta/2` ensembles are not available.
- The paper re-enters a particle that leaves through the top at the bottom with A and B
  swapped; here it reflects. The two agree while the free streams stay unmixed.
- The paper's solver is compressible at M = 0.3; PICurv's is incompressible.
