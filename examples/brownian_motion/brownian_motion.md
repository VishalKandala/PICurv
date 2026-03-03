# Brownian Motion Verification Test

## Purpose

Isolate and verify pure Brownian (stochastic) particle diffusion by comparing
simulated MSD to the Einstein relation. This is a component-level test — it
specifically exercises `UpdateAllParticlePositions` (Euler-Maruyama step in
`ParticleMotion.c`) and `ComputeEulerianDiffusivity` in isolation, with no
interference from mean advection, LES closures, or IEM micro-mixing.

---

## Analytical Reference

For a particle undergoing isotropic Brownian motion with diffusivity Γ,
the mean square displacement is:

    MSD(t) = 6 * Γ * t          [3D, isotropic]

With the test parameters:

    Re    = 1000
    Sc    = 1.0
    ν     = 1/Re = 1e-3             (non-dimensional kinematic viscosity)
    Γ_eff = ν / Sc = 1e-3           (non-dimensional, LES off so ν_t = 0)

    Predicted MSD slope = 6 * Γ_eff = 6.0e-3  [per unit non-dim time]

    At t* = 2.0 (step 2000):  MSD* = 0.012

The Euler-Maruyama discretization introduces O(dt) weak error.
With dt* = 0.001, the discretization bias is negligible relative to
the statistical noise at N = 50000.

---

## Run Command

```bash
# From project root
./bin/picurv run --solve --post-process -n 4 \
  --case  examples/brownian_motion/brownian_motion.yml \
  --solver examples/brownian_motion/Analytical-Zero.yml \
  --monitor examples/brownian_motion/Standard_Output.yml \
  --post examples/brownian_motion/brownian_analysis.yml
```

Output file: `<run_dir>/output/BrownianStats_msd.csv`

---

## Verification Procedure

### Step 1: Load MSD output

```python
import numpy as np
import pandas as pd

df = pd.read_csv("output/BrownianStats_msd.csv")
# Expected columns: step, time, msd  (or msd_x, msd_y, msd_z, msd_total)
t   = df["time"].values      # Non-dimensional time
msd = df["msd"].values       # Total MSD (or msd_x + msd_y + msd_z)
```

### Step 2: Linear fit (skip t=0 to avoid IC artefact)

```python
from scipy.stats import linregress

mask  = t > 0
slope, intercept, r, _, _ = linregress(t[mask], msd[mask])

print(f"Fitted slope:      {slope:.6f}")
print(f"Expected slope:    {6e-3:.6f}")
print(f"Relative error:    {abs(slope - 6e-3) / 6e-3 * 100:.2f}%")
print(f"R^2:               {r**2:.6f}")
print(f"Intercept:         {intercept:.6e}")
```

### Step 3: Pass/Fail Criteria

| Metric              | Pass Threshold         | Notes                                   |
|---------------------|------------------------|-----------------------------------------|
| Slope error         | < 2%                   | Monte Carlo noise ~ 0.45% for N=50000   |
| R²                  | > 0.999                | Brownian MSD is exactly linear in time  |
| Intercept           | < 1e-5 (abs)           | Point source → MSD(0) = 0 exactly       |
| MSD(t*=2.0)         | 0.012 ± 0.0003         | 3-sigma tolerance for N=50000           |

---

## Interpreting Failures

### Slope too high
- `nu_t > 0`: LES was not disabled. Check `-les 0` in solver.yml passthrough.
- Wrong Sc value passed to solver. Check `-schmidt_number` flag in control file.

### Slope too low
- `Gamma_eff` computation bug: diffusivity not computed or zeroed out.
- Particles not moving: check Euler-Maruyama noise term in `UpdateAllParticlePositions`.

### Non-linear MSD (R² < 0.999)
- Early transient: remove first 2-3 snapshots from fit (startup effects).
- Unexpected particle loss at walls: if parameters are changed significantly
  (larger Gamma, longer run, source off-center), the cloud can reach WALL
  boundaries and particles are removed. Removed particles are the farthest-
  displaced ones, which systematically underestimates MSD. At default
  parameters the wall is 8 sigma away so this cannot happen. Check particle
  count vs. time if the MSD curve flattens unexpectedly.

### Non-zero intercept
- Particles did not start at the point source. Check `init_mode: PointSource`
  and `point_source: {x: 0.5, y: 0.5, z: 0.5}` in case.yml.
- Initial displacement bias in walking search seeding. Inspect step-0 particle output.

### Anisotropic MSD (if per-component output available)
- MSD_x ≠ MSD_y ≠ MSD_z with uniform grid is a red flag for a directional
  bias in the noise generator (Box-Muller implementation) or interpolation.
- Each component should converge to 2 * Γ_eff * t independently.

---

## Secondary Check: Gaussian Profile

At each output time, the particle position PDF should be Gaussian:

    p(x, t) = (1 / (4πΓt)^(3/2)) * exp(-|x - x_0|² / (4Γt))

Check in ParaView by loading `Particle_*.vtp` and binning particle positions
in x, y, z at a late timestep. A Gaussian fit should match with width
σ = sqrt(2Γt*) = sqrt(2 * 1e-3 * 2.0) = 0.063.

---

## Known Limitations of This Test

1. MSD is computed from current positions, not unwrapped displacements.
   With the chosen parameters, particle spread is ~0.19 maximum, so
   periodic wrapping is rare. If parameters are changed to increase Γ or
   run time, verify spread stays well below L/2 = 0.5.

2. The drift correction term (∂Γ/∂x * dt) vanishes on a uniform Cartesian
   grid. This test does NOT verify the gradient interpolation pipeline.
   A separate stretched-grid test is needed for that.

3. This verifies the particle position update chain only. Scalar (Psi/IEM)
   correctness requires a separate micro-mixing test.

---

## Reuse Notes

This example is intentionally modular:

- keep `brownian_motion.yml` when you want the same Brownian particle setup,
- swap `Analytical-Zero.yml` for another analytical or solve-mode `solver.yml` only if you intentionally want a different carrier-flow model,
- reuse `brownian_analysis.yml` as a post recipe for other diffusion-focused particle tests.

For broader profile-composition examples, see:

- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
- https://vishalkandala.me/picurv-docs/32_Analytical_Solutions.html
