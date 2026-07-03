@page 55_Newton_Krylov_Momentum_Solver Newton--Krylov Momentum Solver

@anchor _Newton_Krylov_Momentum_Solver

This page documents PICurv's matrix-free Newton--Krylov momentum solver: what it
solves, how it is configured, how to read its convergence output, and the
residual-purity invariant that makes it work. It is one of two momentum-solution
approaches in PICurv; see @ref 31_Momentum_Solvers for how it compares to the
dual-time Picard--Jameson solver and how to select between them.

@tableofcontents

@section p55_scope_sec 1. Purpose and Current Scope

The Newton--Krylov solver advances the implicit momentum update by solving the
nonlinear momentum residual directly with PETSc's `SNES`, using a matrix-free
(Jacobian-free) Krylov linearization. It is implemented in
@ref MomentumSolver_NewtonKrylov (`src/momentum_newton_krylov.c`) and is selected
with `strategy.momentum_solver: "Newton Krylov"`.

Version one is deliberately narrow and validates its inputs up front
(@ref MomentumSolver_NewtonKrylov rejects anything outside this set):

- exactly **one block**, no immersed boundaries, no moving/rotating bodies, FSI,
  or reference frames;
- **no RANS, Clark, TwoD masking, or wall functions**;
- **fresh starts only** (`StartStep == 0`);
- no masked solid cells (`Nvert` must be fluid everywhere);
- boundary handlers limited to no-slip walls, the constant/parabolic/file inlets,
  the conservation outlet, and geometric periodic faces (with paired periodic
  axes).

Within that scope it is a drop-in alternative to the dual-time Picard--Jameson
solver and shares the same fractional-step projection, BDF time discretization,
boundary system, and pressure solve.

@section p55_equation_sec 2. What Is Being Solved

Each physical timestep advances the contravariant velocity `Ucont` by solving the
discrete momentum residual to zero:

\f[
F(\mathbf{U}) \;=\; -\,\mathrm{RHS}_{\text{spatial}}(\mathbf{U})
\;+\; \frac{a_0}{\Delta t}\,\mathbf{U}
\;-\; (\text{BDF history terms}) \;=\; 0,
\f]

where `RHS_spatial` is the convective + viscous + source assembly (@ref ComputeRHS)
and the time term is the BDF discretization added by @ref ComputeTotalResidual —

- **BDF1** (first physical step): \f$(\mathbf{U}^{n} - \mathbf{U}^{n-1})/\Delta t\f$, \f$a_0 = 1\f$;
- **BDF2** (subsequent steps): \f$(1.5\,\mathbf{U}^{n} - 2\,\mathbf{U}^{n-1} + 0.5\,\mathbf{U}^{n-2})/\Delta t\f$, \f$a_0 = 1.5\f$.

Order selection is centralized in `MomentumUsesBDF2()`/`MomentumBDFCoefficient()`,
shared with the momentum stability estimate. The Newton solver does **not** change
the residual arithmetic, BDF coefficients, conservation-outlet formulas, or the
number of boundary passes; it only changes how the resulting nonlinear system is
solved.

@section p55_snes_sec 3. SNES, Matrix-Free Jacobian, and Krylov Solve

The solver builds a per-step `SNES` (@ref MomentumSolver_NewtonKrylov):

- **`SNESNEWTONLS`** with a backtracking (`bt`) line search;
- a **matrix-free operator** from `MatCreateSNESMF`, whose action is the
  finite-difference directional derivative
  \f$J\mathbf{v} \approx [F(\mathbf{X}+h\mathbf{v}) - F(\mathbf{X})]/h\f$
  (`MatMFFDComputeJacobian`). No Jacobian is assembled;
- an inner **`KSPGMRES`** linear solve with `PCNONE` (unpreconditioned). `PCNONE`
  is required in version one and is enforced: option processing that selects any
  other preconditioner is rejected.

Because the Jacobian action is a finite difference of the residual, the residual
**must be a deterministic function of the trial vector** `X` (see
@ref p55_seed_sec). GMRES uses PETSc's default classical Gram--Schmidt
orthogonalization; modified Gram--Schmidt is **not** required and is not enabled.

The Newton loop is: evaluate `F(X)` → form Krylov solve of `J dX = -F` → line
search along `dX` → repeat until an `SNES` convergence test fires. On convergence
the solution is committed into `Ucont`; on failure the entry state is restored
(rollback) and the physical step is reported as not converged
(`simCtx->mom_last_converged`).

@section p55_residual_sec 4. Residual Evaluation and Boundary Application

One residual evaluation (@ref MomentumNewtonKrylov_FormResidual) performs, in order:

1. copy the trial `X` into the global `Ucont`;
2. synchronize periodic staggered `Ucont` and refresh local `lUcont` ghosts;
3. **deterministic Cartesian seed** (Section 5): reconstruct `Ucat` from the
   current `lUcont`, finalize periodic `Ucat`, and refresh `lUcat` ghosts;
4. @ref ApplyBoundaryConditions — the standard **three internal boundary passes**
   (inlet/farfield/wall/conservation-outlet handlers, then Cartesian
   reconstruction, dummy/periodic/corner updates) per call;
5. @ref ComputeTotalResidual — spatial RHS + BDF time terms + residual boundary
   enforcement;
6. `F = -Rhs`, then a **constrained-row** pass that replaces every
   non-independent row (fixed boundary-normal, homogeneous dummy/tangential, and
   periodic-duplicate rows) with an explicit algebraic equation so the matrix-free
   operator has no zero Jacobian rows.

The three internal boundary passes are unchanged and remain necessary: each pass
refreshes the Cartesian state after a boundary correction so the next pass sees a
consistent field.

@section p55_seed_sec 5. Deterministic Cartesian Seeding (Why It Is Required)

This is the invariant that makes the matrix-free solve correct, and it is easy to
break by "simplifying" the residual, so it is documented explicitly.

**Every evaluation of `F(X)` must start from velocity fields derived from that same
`X`.** The conservation-outlet handler reads the cell-centered Cartesian velocity
`lUcat` during its **first** boundary sweep (it measures the uncorrected outlet
flux and builds the outlet profile from it). If `lUcat` were left over from a
previous residual or matrix-free evaluation, `F(X)` would depend on that hidden
state — two evaluations at the same `X` could differ, and the finite-difference
Jacobian action would be inconsistent.

The velocity-state relationship is:

```
Ucont  : global staggered contravariant flux  (the SNES unknown, = X)
lUcont : local ghosted copy of Ucont
Ucat   : global cell-centered Cartesian velocity, reconstructed from lUcont
lUcat  : local ghosted copy of Ucat  (read by the conservation outlet)
```

So the residual seeds them in exactly this dependency order **before** the first
boundary pass:

```
X / Ucont
  -> synchronize periodic staggered Ucont, refresh lUcont
  -> Contra2Cart:                 reconstruct global Ucat from lUcont
  -> SynchronizePeriodicCellFields("Ucat"):  finalize periodic Ucat planes
  -> UpdateLocalGhosts("Ucat"):   refresh lUcat ghosts (outlet reads lUcat)
  -> ApplyBoundaryConditions (three internal passes)
```

Important subtleties, all captured in the source comment above the seed:

- `Contra2Cart()` alone is **not** sufficient: it rebuilds the interior of the
  global `Ucat` but does not refresh `lUcat` (nor `lUcont`), and the outlet reads
  the local ghosted `lUcat`.
- The reconstruction that already happens **inside** `ApplyBoundaryConditions()`
  runs *after* each handler sweep, so it prepares passes two and three — it cannot
  prepare the very first outlet read of pass one.
- `SynchronizePeriodicCellFields("Ucat")` must run before the ghost scatter so
  periodic duplicate planes are finalized consistently (it is a no-op when no
  direction is periodic).

Removing or shortening this sequence reintroduces a history-dependent residual and
invalidates the Newton directions. A permanent regression guards it (Section 10).

@section p55_config_sec 6. Configuration

Select the solver and (optionally) tune its PETSc controls. Omitted fields keep
the defaults established in `src/momentum_newton_krylov.c` and PETSc.

```yaml
strategy:
  momentum_solver: "Newton Krylov"   # -> -mom_solver_type newton_krylov

momentum_solver:
  newton_krylov:
    nonlinear_solver:
      method: "newtonls"             # -> -mom_nk_snes_type
      absolute_tolerance: 1.0e-10    # -> -mom_nk_snes_atol
      relative_tolerance: 1.0e-8     # -> -mom_nk_snes_rtol
      step_tolerance: 1.0e-12        # -> -mom_nk_snes_stol
      max_iterations: 25             # -> -mom_nk_snes_max_it
      line_search:
        type: "bt"                   # -> -mom_nk_snes_linesearch_type
    linear_solver:
      method: "gmres"                # -> -mom_nk_ksp_type
      absolute_tolerance: 1.0e-10    # -> -mom_nk_ksp_atol
      relative_tolerance: 1.0e-6     # -> -mom_nk_ksp_rtol
      max_iterations: 400            # -> -mom_nk_ksp_max_it
      gmres:
        restart: 80                  # -> -mom_nk_ksp_gmres_restart
      preconditioner:
        type: "none"                 # -> -mom_nk_pc_type (only "none" is supported)
```

Field-by-field mappings and validation rules (nonnegative tolerances, positive
iteration/restart counts, `preconditioner.type: none` only) are the authoritative
configuration reference in @ref 08_Solver_Reference "Solver Reference", section 4.
The complete annotated template is `examples/master_template/master_solver.yml`.

Three configuration layers interact, in increasing precedence:

1. **User-facing YAML** (`momentum_solver.newton_krylov.*`) — the supported surface.
2. **C/PETSc defaults** — used for any field you omit.
3. **`petsc_passthrough_options`** — raw PETSc options applied last; they can
   override structured values and are an advanced escape hatch.

The tolerances above are a reasonable starting point. Interpretation:

- `nonlinear_solver.absolute_tolerance` stops Newton when the nonlinear residual
  norm falls below it — the primary physical convergence gate.
- `nonlinear_solver.relative_tolerance` stops Newton relative to the initial
  residual norm.
- `linear_solver.relative_tolerance` controls how tightly each inner GMRES solve
  is converged; a loose `1e-6` inexact-Newton setting is typical and cheap.

@section p55_monitors_sec 7. Monitors and Log Output

Newton--Krylov monitors are enabled under `solver_monitoring.momentum` (see
@ref 09_Monitor_Reference):

- `newton_krylov_history` -> `-mom_nk_pic_monitor`: PICurv's own per-iteration
  nonlinear-norm history;
- `snes_monitor` -> `-mom_nk_snes_monitor`, `snes_converged_reason` ->
  `-mom_nk_snes_converged_reason`;
- `ksp_monitor` -> `-mom_nk_ksp_monitor`, `ksp_converged_reason` ->
  `-mom_nk_ksp_converged_reason`.

Independently of PETSc monitors, the solver writes structured rank-zero logs into
`log_dir`:

- `Momentum_Solver_Newton_Krylov_History_Block_<b>.log`: one row per Newton
  iteration (`step | block | newton | nonlinear_norm`);
- `Momentum_Solver_Newton_Krylov_Summary_Block_<b>.log`: one row per physical step
  (`SNES reason`, Newton iterations, residual evaluations, Krylov iterations,
  initial/final norm, and whether the result was committed or rolled back).

A healthy solve on the validated duct case shows the nonlinear norm dropping by
several orders of magnitude in about two Newton iterations with accepted line
search `lambda = 1`.

@section p55_reasons_sec 8. Convergence Reasons and Failure Modes

Do not treat all non-convergence the same — the `SNES`/`KSP` reason identifies the
failure class:

- **`CONVERGED_FNORM_ABS` / `CONVERGED_FNORM_RELATIVE`**: success (absolute or
  relative nonlinear tolerance met).
- **`DIVERGED_MAX_IT`**: hit `snes_max_it` without meeting a tolerance — usually
  under-resolved inner solves or too tight a nonlinear tolerance for the timestep;
  loosen `nonlinear_solver.relative_tolerance` or reduce `dt`.
- **`DIVERGED_LINEAR_SOLVE`**: an inner GMRES solve failed to converge — inspect
  `ksp_converged_reason`; raise `linear_solver.max_iterations` or `gmres.restart`,
  or loosen `linear_solver.relative_tolerance`.
- **`DIVERGED_LINE_SEARCH`**: the backtracking line search could not find a
  sufficient decrease — typically a poor Newton direction. In this solver that
  most often means the residual was **not deterministic** (a broken Cartesian
  seed, Section 5); it should not occur with the shipped residual.
- **GMRES breakdown / slow Krylov convergence**: a stagnating or breaking-down
  inner solve. With a pure residual the default classical Gram--Schmidt is
  sufficient; enabling modified Gram--Schmidt is **not** the intended remedy.
- **Nonlinear stagnation**: the nonlinear norm plateaus above tolerance — often a
  physically stiff step; reduce `dt`.

Troubleshooting workflow: enable `snes_monitor` + `snes_converged_reason` +
`ksp_converged_reason`, reproduce on a short run, and classify by reason before
changing tolerances. If you observe `DIVERGED_LINE_SEARCH` or non-repeatable
nonlinear norms, suspect residual determinism (Section 5) rather than the Krylov
settings.

@section p55_precond_sec 9. Preconditioning Status and Limitations

Version one supports and validates **only the unpreconditioned path**
(`preconditioner.type: none`, `PCNONE`), enforced at runtime. There is no
assembled operator to precondition against, so a matrix-free preconditioner would
be required; this is future work and is intentionally not exposed. Other current
limitations are the version-one scope restrictions in Section 1.

@section p55_validation_sec 10. Validation Coverage

The Newton--Krylov path is covered at two levels:

- **Default suite** (`unit-newton-krylov`, part of `make check`): constraint-row
  Jacobian structure, matrix-free vs direct differencing, small solve/rollback,
  and residual repeatability. The conservation-outlet conditioned-row derivative
  test doubles as a **seed-removal detector**: removing the deterministic
  Cartesian seed makes that row's self-derivative revert to the decoupled
  artifact and the test fails.
- **Opt-in integration regression** (`make unit-momentum-newton-boundary-fixedpoint`,
  one and four ranks): on the production-sized straight duct it advances a real
  physical step 1, then verifies (i) residual purity at the step-2 state
  (immediate and after real MFFD products), (ii) a complete step-2 solve with the
  default classical Gram--Schmidt, (iii) that the converged three-pass solution
  also zeros the 24-pass outlet residual, and (iv) clean pressure projection.

Validated behavior on that case: convergence in about two Newton iterations from
the true projected step-1 state, identical results with classical and modified
Gram--Schmidt, and divergence-free projection, on both one and four ranks.

@section p55_vs_jameson_sec 11. Differences from Dual-Time Picard--Jameson

| Aspect | Newton--Krylov | Dual-Time Picard--Jameson |
| --- | --- | --- |
| Linearization | true Newton (matrix-free `Jv`) | Picard fixed-point / pseudo-time smoothing |
| Inner solve | PETSc `SNES` + GMRES | staged Jameson RK pseudo-time |
| Main controls | SNES/KSP tolerances, GMRES restart | pseudo-CFL, pseudo-iterations |
| Maturity | newer, narrow validated scope (Section 1) | established, broadly exercised |
| Failure surface | SNES/KSP convergence reasons | pseudo-CFL rollback / rejection |

See @ref 24_Dual_Time_Picard_Jameson_RK for the Picard--Jameson solver and
@ref 31_Momentum_Solvers for selection guidance.

@section p55_refs_sec 12. Related Pages

- **@ref 31_Momentum_Solvers** — momentum-solver selection and status.
- **@ref 24_Dual_Time_Picard_Jameson_RK** — the alternative momentum solver.
- **@ref 08_Solver_Reference** — authoritative configuration reference.
- **@ref 09_Monitor_Reference** — monitor keys and log locations.
- **@ref 23_Fractional_Step_Method** — the projection workflow the solver plugs into.
