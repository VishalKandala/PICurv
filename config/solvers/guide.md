# PICurv Solver Profiles Library

This directory stores reusable `solver.yml` profiles that define numerical strategy. Solver profiles are where you tune stability, convergence, and pressure-correction behavior while keeping case geometry and postprocessing independent.

## 1. What Solver Profiles Control

- solver mode (`solve`, `load`, `analytical`),
- momentum strategy and tolerances,
- strategy-specific option blocks (for example `dual_time_picard_rk4`),
- pressure/Poisson/multigrid settings,
- raw PETSc passthrough options.

## 2. Included Profiles

- `Imp-MG-Standard.yml`: baseline implicit multigrid-oriented setup for general use.
- `Analytical-UniformFlow.yml`: reusable constant-velocity analytical profile for deterministic particle advection checks.
- `Analytical-Zero-Verification-LinearDiffusivity.yml`: reusable verification-oriented zero-flow profile with linear diffusivity override for `grad(Gamma)` drift checks.

For full schema coverage, see:
- `examples/master_template/master_solver.yml`

## 3. Usage Pattern

```bash
./bin/picurv run --solve -n 8 \
  --case my_study/case.yml \
  --solver config/solvers/Imp-MG-Standard.yml \
  --monitor config/monitors/Standard_Output.yml
```

## 4. Momentum Solver Selector Notes

Use `strategy.momentum_solver` with exact values:

- `Explicit RK4`
- `Dual Time Picard RK4`

Any newly introduced selector should only be exposed after parser normalization, runtime dispatch, and docs/tests are updated in one cohesive change.

Verification-only source overrides should stay under the structured `verification.*` namespace and be implemented in `verification_sources.*` rather than as one-off production flags.

## 5. CFD Tuning Order (Practical)

1. Stabilize timestep and pseudo-CFL behavior.
2. Tune momentum tolerances.
3. Tune Poisson tolerance and multigrid sweeps.
4. Revisit grid quality and BC consistency if instability persists.
