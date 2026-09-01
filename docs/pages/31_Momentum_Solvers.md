@page 31_Momentum_Solvers Momentum Solver Implementations

@anchor _Momentum_Solvers

This page tracks momentum-solver options accepted by the current configuration and their runtime implementation status.

@tableofcontents

@section p31_selection_sec 1. Selection and Dispatch

Runtime selection is controlled by `-mom_solver_type`, produced from `solver.yml` (`strategy.momentum_solver`).
Dispatch currently happens in function @ref FlowSolver within the step orchestrator.

Accepted YAML values:

- `Explicit RK4` -> `EXPLICIT_RK`
- `Dual Time Picard Jameson RK` -> `DUALTIME_PICARD_JAMESON_RK`
- `Newton Krylov` -> `newton_krylov`

Only implemented values are exposed in the enum, parser, and dispatcher. New
solver values should be added only with a real implementation plus matching
parser, docs, and test updates.

For compatibility, the former `Dual Time Picard RK4` YAML display name and
`DUALTIME_PICARD_RK4` C CLI value still select the Jameson solver. New
configuration and code must use the Jameson names.

@section p31_approaches_sec 1a. The Two Momentum-Solution Approaches

Beyond the explicit RK4 verification path, PICurv provides **two distinct
implicit-in-physical-time momentum-solution approaches**. They are genuinely
different algorithms, not two names for the same method, and they have different
numerical behavior, controls, and maturity. Choose one with
`strategy.momentum_solver`.

**Dual-time Picard--Jameson** (`Dual Time Picard Jameson RK`) — the established,
comparatively robust default. It advances the implicit BDF2 update with a
fixed-point / pseudo-time iteration using staged Jameson RK smoothing. It is
controlled through pseudo-CFL and pseudo-iteration settings, may need
conservative pseudo-CFL values, and can converge slowly in demanding
high-Reynolds-number or near-inviscid regimes. It does **not** use SNES/GMRES
matrix-free Newton linearizations. Full details: @ref 24_Dual_Time_Picard_Jameson_RK.

**Newton--Krylov** (`Newton Krylov`) — a newer matrix-free nonlinear solver. It
solves the momentum residual with PETSc `SNES`, using matrix-free
Jacobian--vector products (finite-difference `Jv`), an inner GMRES Krylov solve,
and a backtracking line search. It exposes nonlinear, line-search, GMRES, and
preconditioner controls. The supported baseline is unpreconditioned
matrix-free differencing; the supported alternative is a provisional
frozen-momentum point-block preconditioner, whose PETSc block-Jacobi backend is
chosen internally. It requires a deterministic residual (its Cartesian boundary state is reconstructed
from the current trial vector before boundary conditions are applied). Its
convergence diagnostics and failure modes (`SNES`/`KSP` reasons) differ from the
Picard solver. It has a **narrower, explicitly validated scope** — see its
dedicated page: @ref 55_Newton_Krylov_Momentum_Solver.

Selection guidance (within the evidence available today):

- Prefer **Dual-time Picard--Jameson** for general production runs, complex
  geometries, and cases outside the Newton--Krylov version-one scope; it is the
  broadly exercised path.
- Consider **Newton--Krylov** on supported single-block cases when you want true
  Newton convergence behavior and `SNES`/`KSP`-style diagnostics, keeping its
  scope restrictions (Section 1 of @ref 55_Newton_Krylov_Momentum_Solver) in mind.
- Use **Explicit RK4** for verification or when the explicit stability limit is
  affordable.

@section p31_status_sec 2. Implementation Status Matrix

- `EXPLICIT_RK`: implemented by @ref MomentumSolver_Explicit_RungeKutta4
- `DUALTIME_PICARD_JAMESON_RK`: implemented by @ref MomentumSolver_DualTime_Picard_JamesonRK
- `newton_krylov`: implemented by @ref MomentumSolver_NewtonKrylov (matrix-free PETSc
  SNES/GMRES; `none` or frozen-momentum/point-block preconditioning; narrow
  validated version-one scope — see @ref 55_Newton_Krylov_Momentum_Solver)

@section p31_controls_sec 3. Numerical Controls In Use

Main controls consumed by implemented solvers:

- `-mom_max_pseudo_steps`
- `-mom_atol`
- `-mom_rtol`
- `-mom_resid_atol`, `-mom_resid_rtol`
- `-pseudo_cfl`, `-min_pseudo_cfl`, `-max_pseudo_cfl`
- `-pseudo_cfl_growth_factor`, `-pseudo_cfl_reduction_factor`
- `-mom_dt_jameson_residual_norm_noise_allowance_factor`
- `-mom_ratio_ema_alpha`

Defaults and final option ingestion are in function @ref CreateSimulationContext during startup parsing.

For the dual-time Jameson solver, `max_iterations` bounds **accepted** pseudo-iterations. A separate
hard cap of `3 × max_iterations` limits total attempts (accepted plus rejected) to prevent infinite
rejection loops. Convergence is decided by the residual: `residual_abs_pass` **OR**
(`residual_rel_pass` **AND** `update_pass`). The absolute residual test
(`|R| ≤ residual_absolute_tol · resid_ref`, dimensionless tolerance) is sufficient on its own;
the relative test is paired with the `relative_tol` update guard. `absolute_tol` takes no part
while a residual tolerance is set. Both residual tolerances default to enabled; setting both
non-positive selects the legacy update-only branch. See @ref p24_convergence_sec.

The dual-time controller uses one global pseudo-CFL and globally accepts or
rolls back a complete four-stage trial. The selected next pseudo-CFL is carried
directly into the next physical timestep. `step_tol`/`-imp_stol` remains
accepted only as a deprecated compatibility input and is unused by active
momentum solvers.

`pseudo_cfl.*` values are **dimensionless Courant numbers** (Phase 3+), not fractions of the physical
timestep `dt`. The solver computes `dtau = pseudo_cfl / lambda_max` where `lambda_max` is the global
maximum convective spectral radius of the current velocity field. This makes `pseudo_cfl` independent
of `dt`, grid size, and flow speed. The stable range for the 4-stage Jameson RK smoother is
`pseudo_cfl ≈ 0–2.83`; the shipped defaults are `initial: 0.5`, `maximum: 2.0`.

@section p31_testing_sec 4. Current test status

Current testing is uneven by solver path:

- dispatch and guardrails are directly covered through `FlowSolver`-side unit tests
- `MomentumSolver_DualTime_Picard_JamesonRK` is exercised mainly through smoke and runtime orchestration
- `MomentumSolver_Explicit_RungeKutta4` still needs a direct positive-path harness

That means the momentum stack is currently a stronger regression gate than bespoke debugging surface.

@section p31_rhs_cadence_sec 6. Call Cadence of the Shared RHS (Read Before Adding State)

Both momentum solvers are built on one residual implementation:

```
ComputeTotalResidual()            (src/momentumsolvers.c)
  └─ ComputeRHS()                 (src/rhs.c)
       └─ ComputeBodyForces()     → individual body forces
```

Picard reaches it once per Jameson RK stage; Newton--Krylov reaches it once per
residual evaluation, including every finite-difference probe of the matrix-free
Jacobian. **Neither calls it once per physical timestep**, and the ratio is not
fixed: it varies with pseudo-iteration count, line-search backtracks, and Krylov
iterations.

That makes the shared RHS a hazardous place to keep state. Anything advanced on
each call - a filter, a ramp, a moving average, an integral controller term, a
relaxation counter - becomes a function of *how many* evaluations happened
before it rather than of the solution. Two things break:

- `MomentumNewtonKrylov_FormResidual()` requires `F(X)` to be a deterministic
  function of the trial vector alone; otherwise the finite-difference Jacobian
  action `(F(X+hv) - F(X))/h` is inconsistent.
- The Picard shadow-Jacobian estimate treats body forces as constant forcing with
  zero velocity Jacobian.

The rule is to gate every such update on `simCtx->step` and reuse the resolved
value for the rest of the step. The same hazard applies to boundary handlers,
because `ApplyBoundaryConditions()` runs each handler's `PreStep` three times per
call and is itself called from inside both solvers' iteration loops.

**Worked example.** The driven periodic flow controller had this defect in two
places at once, and fixing only the first was not enough:

| State | Location | Symptom when ungated |
|---|---|---|
| `bulkVelocityCorrection` | `PreStep` in `src/BC_Handlers.c` | source recomputed from the trial vector mid-solve |
| smoothing EMA on the force | `ComputeDrivenChannelFlowSource()`, `src/BodyForces.c` | applied force walked 0.5, 0.75, 0.875 ... toward target within one step |

Measured on a 4-step run: 379 force evaluations produced **42 distinct force
values** before the fix and **3** after - one per timestep.
`tests/smoke/run_driven_periodic_regression.sh` asserts the force is piecewise
constant per step. See @ref p54_driven_cadence_sub and the contract in
`include/BodyForces.h`.

@section p31_extension_sec 7. Adding A New Momentum Solver

Required steps:

1. define solver implementation function in `src/momentumsolvers.c`,
2. ensure enum and parser mapping are present (`variables.h`, `setup.c`, `picurv_cli/core.py`),
3. add dispatch branch in function @ref FlowSolver for the new enum value,
4. expose and document solver-specific YAML options,
5. add smoke tests and docs updates.

For user-facing contract updates, also update:

- **@subpage 08_Solver_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 50_Modular_Selector_Extension_Guide**
