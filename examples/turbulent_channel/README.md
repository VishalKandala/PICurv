# Seeded turbulent plane channel

A wall-resolved channel startup at nominal `Re_tau ~ 180`, using a generated
three-dimensional perturbation. This replaces the former coarse wall-modelled LES
setup. It is a candidate DNS configuration, not a validated turbulent benchmark.

The Cartesian axes are spanwise `i`, wall-normal `j`, and streamwise `k`.
`grid.gen` reads `grid.cfg`: 128 x 128 x 256 cells on a 2*pi x 2 x 4*pi domain,
with symmetric wall clustering and three-level multigrid-compatible node counts.
At nominal Re_tau=180, both periodic spacings are about 8.84 wall units. The first
full wall cell is 1 wall unit wide, placing its center at y+=0.5. These nominal
spacings must be recomputed from the achieved wall stress; grid/time convergence
is still required. The viscosity gives U_b*h/nu=2800. Re_tau is an outcome of
constant-flux driving, not an input imposed by the Reynolds-number setting.

`channel_spectral_velocity` generates a discrete-curl perturbation using seeded
random Fourier content in the periodic directions and smooth basis functions in
the wall-normal direction. Its volume-weighted component-equivalent RMS is 0.1 U_b.
The parabolic mean is normalized to discrete bulk velocity 1; the driven target
flux is 4*pi. There is no recurring perturbation injection. Transition is not
guaranteed by this amplitude; monitor fluctuation energy and Reynolds shear stress.

## Generate and inspect before solving

```bash
picurv init turbulent_channel --dest my_channel
cd my_channel
picurv validate --case config/case.yml --solver config/solver.yml --monitor config/monitor.yml --post config/post.yml
picurv precompute --case config/case.yml
```

Precompute publishes the grid and IC as workspace assets. Their payload retains the canonical grid, initial-condition, metrics, and spectra
locations; staging a run materializes those existing assets.
The IC summary reports realized bulk velocity, perturbation RMS, and discrete
divergence. Initial plane spectra sample zero-based physical wall-normal cells
8, 32, and 63. Each has a CSV and a JSON summary recording the actual location,
mean subtraction, normalization and Parseval residual. Indices exclude dummy cells.

Spectra contain signed streamwise/spanwise wavenumbers and per-component modal
energies. Summing all signed modes recovers sample kinetic energy. They measure
the seed, before runtime velocity reconstruction; they are not an expected
fully developed turbulence spectrum. For a real line instead, use e.g.
`{task: line_spectrum, axes: [k], fixed_indices: {i: 32, j: 32}, subtract_mean: sample}`.
No neighboring planes or lines are averaged.

## Cluster experiment

The shipped 1000 steps cover only startup. Measure a short pilot at the intended
grid and MPI layout and extrapolate before choosing campaign length and walltime.
Copy the repository's `examples/master_template/master_cluster.yml` to your case
as `cluster.yml`; replace account, partition, node/rank counts, walltime, module
setup, and notification settings with your site's values. Then stage with:

```bash
picurv run --solve --case config/case.yml --solver config/solver.yml --monitor config/monitor.yml --cluster cluster.yml --no-submit
picurv submit --run-dir runs/<run_id>
picurv run --post-process --only spectra --run-dir runs/<run_id> --post config/post.yml
```

The spectra-only stage uses Python and reads committed raw Ucat checkpoints,
including step zero; it does not launch the C field postprocessor. Its files live
under the run's recipe-specific spectra directory. The initial plane and line
products use the same computation as checkpoint products.

Field statistics are disabled during startup. After wall stress, bulk velocity,
and fluctuation energy become stationary, choose a developed averaging window,
enable it in the monitor, and add `field_statistics` with that window to the post
recipe. Compare mean velocity, all velocity RMS components, Reynolds shear stress,
wall friction and total-shear balance against matching-Reynolds-number channel
reference data. Account for sampling uncertainty and resolution sensitivity.
A visually fluctuating field or a successful short solve is not that validation.
