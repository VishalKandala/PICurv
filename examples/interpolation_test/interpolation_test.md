# Interpolation Accuracy Test

## Purpose

Quantify the accuracy of grid-to-particle velocity interpolation by advecting
particles in the 3D Taylor-Green Vortex (TGV3D) analytical flow field. Because
the analytical velocity is known exactly at every point, interpolated particle
velocities can be compared against the true solution to measure L2 and Linf
interpolation error.

---

## Analytical Reference

The TGV3D velocity field at position (x, y, z) and time t is:

    u =  sin(x)*cos(y)*cos(z) * exp(-2*nu*t)
    v = -cos(x)*sin(y)*cos(z) * exp(-2*nu*t)
    w =  0

With the test parameters:

    Re  = 100
    nu  = 1/Re = 0.01
    dt* = 0.01,  N_steps = 200  -->  t_final* = 2.0

Velocity decay factor at t* = 2.0:

    exp(-2 * nu * k^2 * t) = exp(-0.04) ~ 0.961

The flow retains 96% of its initial strength throughout the run, so
interpolation error is measured against a near-constant-amplitude field.

Domain: [0, 2*pi]^2 x [0, 0.2*pi], with 32^3 uniform cells.
Cell sizes: dx = dy ~ 0.196, dz ~ 0.020.

---

## Files in this Example

- `interpolation_test.yml` -- case definition (physics, grid, BCs, run control)
- `Analytical-TGV.yml` -- solver profile (analytical mode, TGV3D field, Trilinear interpolation)
- `Standard_Output.yml` -- monitor profile
- `interpolation_analysis.yml` -- post-processor recipe (particle VTP output)
- `execution.example.yml` -- optional site launcher example

---

## Run Command

```bash
# From project root
./bin/picurv run --solve --post-process -n 4 \
  --case  examples/interpolation_test/interpolation_test.yml \
  --solver examples/interpolation_test/Analytical-TGV.yml \
  --monitor examples/interpolation_test/Standard_Output.yml \
  --post examples/interpolation_test/interpolation_analysis.yml
```

Output files:
- `<run_dir>/logs/interpolation_error.csv` -- per-step error statistics
- `<run_dir>/viz/Particle_*.vtp` -- particle positions and velocities for ParaView

---

## Particle Count Study

This example also supports a particle-count sensitivity study where only
`models.physics.particles.count` changes while the grid, timestep, analytical
field, and interpolation method remain fixed.

Recommended counts:

- 50,000
- 500,000
- 5,000,000

### Option A: Local Parameter Study

For a workstation or quick debugging workflow, repeat the same run with
different case variants. Keep:

- `Analytical-TGV.yml`
- `Standard_Output.yml`
- `interpolation_analysis.yml`

unchanged, and change only `models.physics.particles.count` in copied case
files such as:

- `interpolation_test_np_50000.yml`
- `interpolation_test_np_500000.yml`
- `interpolation_test_np_5000000.yml`

Example local command:

```bash
./bin/picurv run --solve --post-process -n 4 \
  --case  examples/interpolation_test/interpolation_test_np_50000.yml \
  --solver examples/interpolation_test/Analytical-TGV.yml \
  --monitor examples/interpolation_test/Standard_Output.yml \
  --post examples/interpolation_test/interpolation_analysis.yml
```

After each run, compare:

- `<run_dir>/logs/interpolation_error.csv`
- max `L2_u`, `L2_v`, `Linf_u`, `Linf_v`
- optional wall-clock/runtime usage from solver logs or scheduler accounting

Use this path when you want a quick manual comparison without Slurm study
artifacts.

### Option B: Job-Based Study with `picurv sweep`

For a cluster or any workflow where you want generated study artifacts and one
aggregate comparison table, use the provided:

- `particle_count_study.yml`

This study uses explicit CSV metrics from `logs/interpolation_error.csv` rather
than the default `msd_final` shorthand.

Stage the study without submitting jobs:

```bash
./bin/picurv sweep \
  --study examples/interpolation_test/particle_count_study.yml \
  --cluster config/schedulers/slurm_default.yml \
  --no-submit
```

Submit immediately:

```bash
./bin/picurv sweep \
  --study examples/interpolation_test/particle_count_study.yml \
  --cluster config/schedulers/slurm_default.yml
```

What to inspect before submission:

- `studies/<study_id>/cases/case_####/` contains one generated case per count
- each generated `interpolation_test.yml` has the intended
  `models.physics.particles.count`
- generated scheduler scripts under `studies/<study_id>/scheduler/`

What to inspect after completion:

- `studies/<study_id>/results/metrics_table.csv` (auto-collected by chained metrics job)
- per-case `logs/interpolation_error.csv`
- scheduler stdout/stderr files for any failed array task

If a case is killed by the walltime guard (common for large particle counts),
continue the study:

```bash
./bin/picurv sweep --continue --study-dir studies/<study_id>
```

To increase walltime for the remaining cases, pass a different cluster config:

```bash
./bin/picurv sweep --continue --study-dir studies/<study_id> \
  --cluster config/schedulers/slurm_longer.yml
```

Use this path when you want a reproducible study directory and cluster-managed
execution.

### Cluster Workflow Checklist

1. Push the repo changes from your workstation and pull them on the cluster.
2. Build or confirm `./bin/picurv`, `bin/simulator`, and `bin/postprocessor`.
3. Copy the interpolation example into a working case directory on the cluster.
4. Copy `particle_count_study.yml` into that directory.
5. Adjust `base_configs` paths in the study file if you changed the directory
   layout.
6. Create or edit a Slurm profile based on `config/schedulers/slurm_default.yml`
   with your account, partition, memory, walltime, and launcher settings.
7. Validate the base case files:

```bash
./bin/picurv validate \
  --case interpolation_test.yml \
  --solver Analytical-TGV.yml \
  --monitor Standard_Output.yml \
  --post interpolation_analysis.yml
```

8. Stage first with `--no-submit` and inspect generated study cases.
9. Submit the study after confirming the three particle counts are correct.
10. Review `results/metrics_table.csv` when all array tasks finish.

Operational note:

- The 5,000,000-particle case can require much more memory and runtime than the
  50k and 500k cases.
- Size `cluster.yml` for the most expensive case because the same scheduler
  profile is reused for every job in the study.
- Start with a conservative `execution.max_concurrent_array_tasks` setting.

### Which Workflow Should I Use?

- Use the local path for quick debugging and metric sanity checks.
- Use the job-based path for production studies, scheduler arrays, and
  `metrics_table.csv` aggregation.

---

## Verification Procedure

### Step 1: Load interpolation error output

```python
import numpy as np
import pandas as pd

df = pd.read_csv("<run_dir>/logs/interpolation_error.csv")
# Columns include: step, t, L2_u, L2_v, L2_w, Linf_u, Linf_v, Linf_w, ...
```

### Step 2: Check error norms

```python
print(f"Max L2 error (u):    {df['L2_u'].max():.6e}")
print(f"Max L2 error (v):    {df['L2_v'].max():.6e}")
print(f"Max Linf error (u):  {df['Linf_u'].max():.6e}")
print(f"Max Linf error (v):  {df['Linf_v'].max():.6e}")
print(f"Max L2 error (w):    {df['L2_w'].max():.6e}")  # Should be ~0 (w=0 analytically)
```

### Step 3: Pass/Fail Criteria

| Metric             | Pass Threshold      | Notes                                          |
|--------------------|---------------------|-------------------------------------------------|
| L2 error (u, v)    | < 0.05              | Trilinear on 32^3 for smooth TGV field          |
| Linf error (u, v)  | < 0.15              | Worst-case near cell corners                    |
| L2 error (w)       | ~ 0 (< 1e-10)       | w = 0 analytically; any nonzero is a bug        |
| Error stability     | Non-growing         | Error should not diverge over the run           |

---

## Interpreting Failures

### L2 error larger than expected
- Wrong interpolation method: verify `interpolation.method: "Trilinear"` in
  `Analytical-TGV.yml`. The legacy `CornerAveraged` method is first-order and
  will show larger error.
- Grid too coarse: at 16^3, errors roughly double due to second-order convergence.

### Linf error spikes
- Particles near domain walls may see boundary artefacts. TGV3D velocity is
  zero at the boundaries, so wall-adjacent particles should have near-zero
  interpolated velocity. Large Linf at walls suggests boundary handling issues.

### Nonzero w-component error
- The TGV3D field has w = 0 everywhere. Any nonzero w error indicates a bug in
  the analytical field setup or interpolation indexing.

### Error grows over time
- Particle positions drift due to interpolation error, placing particles in
  different cells over time. Mild growth is expected. Rapid divergence suggests
  an instability in the advection or interpolation pipeline.

---

## Grid Convergence Study

To verify second-order convergence of Trilinear interpolation, run at multiple
resolutions by editing `grid.programmatic_settings.im/jm/km` in
`interpolation_test.yml`:

| Resolution | im=jm=km | Expected L2 scaling |
|------------|----------|---------------------|
| Coarse     | 16       | ~4x baseline error  |
| Baseline   | 32       | 1x (reference)      |
| Fine       | 64       | ~0.25x baseline     |
| Very fine  | 128      | ~0.0625x baseline   |

The ratio of errors between successive doublings should approach 4 (second order).

---

## Comparing Interpolation Methods

Switch `interpolation.method` in `Analytical-TGV.yml` between `"Trilinear"`
and `"CornerAveraged"` to compare methods:

```yaml
interpolation:
  method: "CornerAveraged"   # Legacy first-order method
```

The CornerAveraged method should show roughly 2x the L2 error of Trilinear at
the same resolution, with first-order (rather than second-order) grid convergence.

---

## Reuse Notes

This example is intentionally modular:

- keep `interpolation_test.yml` when you want the same TGV3D particle setup,
- swap `Analytical-TGV.yml` to test other analytical fields or interpolation methods,
- reuse `Standard_Output.yml` or swap for a lighter/heavier monitor profile,
- adapt `interpolation_analysis.yml` to extract additional particle fields.

---

## Live Docs

- https://vishalkandala.me/picurv-docs/32_Analytical_Solutions.html
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
