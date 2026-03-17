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
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
