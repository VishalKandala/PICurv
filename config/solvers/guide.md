# PICurv Solver Profiles Library

This directory stores reusable `solver.yml` profiles that define numerical strategy. Solver profiles are where you tune stability, convergence, and pressure-correction behavior while keeping case geometry and postprocessing independent.

## 1. What Solver Profiles Control

- solver mode (`solve`, `load`, `analytical`),
- momentum strategy and tolerances,
- strategy-specific option blocks (`dual_time_picard_jameson_rk` and `newton_krylov`),
- scalar transport properties such as Schmidt numbers,
- pressure/Poisson/multigrid settings,
- raw PETSc passthrough options.

## 2. Included Profiles

- `Imp-MG-Standard.yml`: baseline implicit multigrid-oriented setup for general use.
- `Newton-Krylov-Standard.yml`: matrix-free Newton--Krylov setup for the validated periodic-flow scope.
- `Newton-Krylov-Frozen-Momentum-Point-Block.yml`: Newton--Krylov setup with the provisional frozen-momentum point-block preconditioner.
- `Analytical-Zero.yml`: zero-velocity carrier field for Brownian and stationary-particle checks.
- `Analytical-UniformFlow.yml`: reusable constant-velocity analytical profile for deterministic particle advection checks.
- `Analytical-TGV3D.yml`: three-dimensional Taylor--Green vortex for interpolation verification.
- `Analytical-Zero-Verification-LinearDiffusivity.yml`: reusable verification-oriented zero-flow profile with linear diffusivity override for `grad(Gamma)` drift checks.
- `Analytical-Zero-Verification-Scalar.yml`: zero-flow scalar-truth injection profile for scatter verification.

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
- `Dual Time Picard Jameson RK`
- `Newton Krylov`

These are two genuinely different implicit momentum approaches plus the explicit
verification path. `Dual Time Picard Jameson RK` is the established pseudo-time
fixed-point solver (tuned through pseudo-CFL and pseudo-iterations). `Newton
Krylov` is the newer matrix-free PETSc SNES/GMRES solver with a narrower validated
scope. See the momentum-solver overview and the dedicated Newton--Krylov page:
- https://vishalkandala.me/picurv-docs/31_Momentum_Solvers.html
- https://vishalkandala.me/picurv-docs/55_Newton_Krylov_Momentum_Solver.html

When `strategy.momentum_solver` is `Newton Krylov`, tune it under the structured
`momentum_solver.newton_krylov` block. Omitted fields keep the C/PETSc defaults
and select the matrix-free/no-preconditioner baseline:

```yaml
strategy:
  momentum_solver: "Newton Krylov"
momentum_solver:
  newton_krylov:
    jacobian:
      type: "finite_difference"
      finite_difference:
        mode: "matrix_free"
    preconditioner:
      model: "none"
    nonlinear_solver:
      method: "newtonls"
      absolute_tolerance: 1.0e-10
      relative_tolerance: 1.0e-8
      step_tolerance: 1.0e-12
      max_iterations: 25
      line_search:
        type: "bt"
    linear_solver:
      method: "gmres"
      absolute_tolerance: 1.0e-10
      relative_tolerance: 1.0e-6
      max_iterations: 400
      gmres:
        restart: 80
```

The other supported mathematical preconditioner is the provisional
`model: frozen_momentum_jacobian` with `structure.type: point_block`. It
approximates only same-cell frozen-momentum coupling; PETSc's block-Jacobi
implementation is selected internally, not in normal YAML. It does not broaden
the Newton--Krylov version-one physics or boundary-condition scope.

Any newly introduced selector should only be exposed after parser normalization, runtime dispatch, and docs/tests are updated in one cohesive change.

Verification-only source overrides should stay under the structured `verification.*` namespace and be implemented in `verification_sources.*` rather than as one-off production flags.

## 4b. The Multigrid Coarse Solve

`poisson_solver.multigrid.level_solvers.level_0` is the **coarse solve** at the
base of the V-cycle, not a smoother, and the evenly spaced `level_N` naming
hides that. `level_1..N` are smoothers and map to `-ps_mg_levels_N_*`; `level_0`
maps to `-ps_mg_coarse_*`.

Because multigrid is used here as a preconditioner, `level_0` must be a **fixed
linear operator**. A Krylov method there (`gmres`, `fgmres`, `cg`, `bcgs`, ...)
makes the whole preconditioner nonlinear, which decouples the outer FGMRES
tracked residual from the true residual `b - Ax` - the solver then reports
convergence on a number that no longer describes the incompressibility
constraint, and the failure is rank-dependent. Pinning `ksp_max_it` does not fix
it, because a fixed iteration count still builds a Krylov space out of its
input. PICurv logs a startup warning if you configure one.

Every shipped profile uses `{method: preonly, preconditioner: redundant}`, which
is appropriate up to roughly 1e4 coarse-grid unknowns; use `telescope` above
that. Size `levels` so the coarsest grid lands around 1e3-1e4 unknowns - a 5M
cell grid wants 5 levels, not 4.

Full discussion, the worked tracked-vs-true residual table, and the
coarsenability constraint:
https://vishalkandala.me/picurv-docs/25_Pressure_Poisson_GMRES_Multigrid.html

## 5. CFD Tuning Order (Practical)

1. Stabilize timestep and pseudo-CFL behavior.
2. Tune momentum tolerances.
3. Tune Poisson tolerance and multigrid sweeps.
4. Revisit grid quality and BC consistency if instability persists.

The standard implicit profile enables hardened momentum convergence with
`residual_relative_tol: 1.0e-3`. Omitting both residual tolerances, or setting
both non-positive, preserves the legacy update-only convergence rule. A looser
`1.0e-2` residual-relative tolerance can be useful for exploratory LES, but
should be treated as an explicit accuracy/runtime tradeoff.

For the Picard-Jameson controller, `max_iterations` limits total attempted
trials. The pseudo-CFL is global, rejected trials roll back globally, and the
controller-selected next CFL carries into the next physical timestep.

Detailed behavior and tuning guidance:
- https://vishalkandala.me/picurv-docs/24_Dual_Time_Picard_Jameson_RK.html
