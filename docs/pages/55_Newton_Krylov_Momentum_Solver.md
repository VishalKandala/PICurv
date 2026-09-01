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
- **no RANS or TwoD masking**;
- fresh starts and Eulerian field restarts (including `--continue`); the first
  solved restart step uses BDF1 because only the checkpoint state is available;
- solid cells are permitted but the immersed-boundary method is not: a masked row
  carries no unknown and is constrained like any other such row, but the IBM velocity
  reconstruction does not run inside this solver's residual;
- boundary handlers limited to no-slip walls, the constant/parabolic/file inlets,
  the conservation outlet, and periodic faces (with paired periodic axes) using
  the `geometric`, `constant_flux`, or `initial_flux` handler.

@note The gradient (Clark) model and wall functions both reach this solver through the
shared residual it already evaluates. Two properties are worth knowing before using
them. The preconditioner does not model the Clark stress, so Krylov counts can rise
with it enabled. Wall functions obtain the friction velocity from an iterative
root-find, and a matrix-free Jacobian action amplifies that solve's tolerance by `1/h`,
so a loose inner tolerance shows up as a degraded Jacobian action rather than as an
error.

Within that scope it is a drop-in alternative to the dual-time Picard--Jameson
solver and shares the same fractional-step projection, BDF time discretization,
boundary system, and pressure solve.

@subsection p55_driven_ssec 1.1 Driven periodic faces

The two driven periodic handlers are inside the validated scope. What makes them
admissible is that their momentum source is **frozen for the whole timestep**.
That takes two independent gates, both keyed on `simCtx->step`, because two
pieces of state sit between the flux measurement and the applied force:

1. `bulkVelocityCorrection`, computed in the handler's `PreStep` from the field
   at the start of the step (`src/BC_Handlers.c`);
2. the 50/50 exponential moving average that smooths the force against the
   previous step, held in `SimCtx` and applied in
   `ComputeDrivenChannelFlowSource()` (`src/BodyForces.c`).

Both functions run once per residual evaluation, so gating only the first still
leaves the applied force walking toward its target across evaluations. With both
gated, every residual evaluation in the Newton solve sees the same body force.

This is worth stating explicitly because it is not automatic.
@ref ApplyBoundaryConditions runs the handler `PreStep` sweep three times per
call, and the residual callback calls it on **every** evaluation - so a handler
that recomputed its source in `PreStep` would let that source drift with the
trial vector, silently changing the operator the Newton solve is converging on
and invalidating the finite-difference Jacobian action. The per-timestep freeze
in `PreStep_PeriodicDrivenConstant` is what prevents that. The boundary trim
(`boundaryVelocityCorrection`) is deliberately *not* frozen, but it acts on
boundary fluxes rather than on the momentum source; see @ref p54_driven_cadence_sub.

Because the source is frozen at the start of a step, the flux controller lags the
target by one step - first order in `dt`. Halving `dt` halves the lag. That is the
same lag the Picard solver sees, so results from the two solvers are comparable.

The paired-face checks below still apply: both faces of the driven axis must be
`PERIODIC` and must carry the same handler.

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

@section p55_snes_sec 3. Residual, Jacobian, Preconditioning Matrix, and PC

The four layers have distinct roles. \f$F(U)=0\f$ is the unchanged nonlinear
equation. The Jacobian operator approximates \f$dF/dU\f$ and is the operator GMRES
uses for the Newton correction. The preconditioning matrix is a cheaper
mathematical approximation to that operator. Finally, the PETSc `PC` is an
internal algorithm that approximately applies the inverse of the supplied
preconditioning matrix.

The solver builds a per-step `SNES` (@ref MomentumSolver_NewtonKrylov):

- **`SNESNEWTONLS`** with a backtracking (`bt`) line search;
- a **matrix-free Jacobian operator** from `MatCreateSNESMF`, whose action is the
  finite-difference directional derivative
  \f$J\mathbf{v} \approx [F(\mathbf{X}+h\mathbf{v}) - F(\mathbf{X})]/h\f$
  (`MatMFFDComputeJacobian`); this is always authoritative;
- by default the preconditioning-matrix argument aliases the Jacobian operator
  and the inner **`KSPGMRES`** uses `PCNONE`;
- optionally, a separately assembled AIJ preconditioning matrix contains the
  provisional frozen-momentum point blocks; its validated internal PETSc
  backend is `PCPBJACOBI`.

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
   reconstruction, dummy/periodic/corner updates) per call. Note that this
   includes each handler's `PreStep`, so a handler carrying per-timestep
   controller state must guard it against being advanced here; see
   @ref p55_driven_ssec;
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
    jacobian:
      type: "finite_difference"
      finite_difference:
        mode: "matrix_free"
    preconditioner:
      model: "none"
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
```

The point-block alternative replaces the `preconditioner` block with:

```yaml
    preconditioner:
      model: "frozen_momentum_jacobian"
      structure:
        type: "point_block"
```

`jacobian.type: finite_difference` means that the complete deterministic
nonlinear residual \f$F(U)\f$ is differentiated numerically.
`jacobian.finite_difference.mode: matrix_free` means PETSc evaluates directional
products on demand and does not assemble the Jacobian. Finite difference is the
construction type; matrix free is one mode within that type.

The planned second finite-difference mode is `colored_sparse`, which would
assemble a sparse numerical Jacobian using coloring. It is future syntax and is
rejected today:

```yaml
    jacobian:
      type: finite_difference
      finite_difference:
        mode: colored_sparse
```

Planned frozen-momentum approximations are a different Jacobian type, not
storage modes of finite difference. Their future forms are:

```yaml
    jacobian:
      type: frozen_momentum_approximation
      frozen_momentum_approximation:
        structure: diagonal
```

or:

```yaml
    jacobian:
      type: frozen_momentum_approximation
      frozen_momentum_approximation:
        structure: full_sparse
```

Those forms are also rejected today.

The Jacobian fields map to `-mom_nk_jacobian_type finite_difference` and
`-mom_nk_jacobian_fd_mode matrix_free`. The preconditioner fields map to
application-owned `-mom_nk_preconditioner_*` selectors. They do not emit a
user-selected PETSc PC type. The released
`linear_solver.preconditioner.type: none` spelling is accepted as a deprecated
compatibility alias; it conflicts with a non-none new model. Field-by-field
mappings and validation rules (nonnegative tolerances and positive iteration/restart counts) are the authoritative
configuration reference in @ref 08_Solver_Reference "Solver Reference", section 4.
The complete annotated template is `examples/master_template/master_solver.yml`.

The planned operator-intent correspondence with legacy selectors is:

| Modern configuration | Legacy operator coverage |
|---|---|
| `finite_difference / matrix_free + none` | top-level `-imp 4` |
| `finite_difference / matrix_free + frozen point block` | `-imp 5 -imp_type 2` |
| `finite_difference / colored_sparse` | future `-imp 5 -imp_type 1` |
| `frozen momentum / diagonal` | future `-imp 5 -imp_type 3` |
| `frozen momentum / full_sparse` | future `-imp 5 -imp_type 4` |

This table maps numerical operator intent only. It does not reproduce legacy
defects, lifecycle, mutable-residual behavior, or historically unknown PETSc
runtime options.

Three configuration layers interact, in increasing precedence:

1. **User-facing YAML** (`momentum_solver.newton_krylov.*`) — the supported surface.
2. **C/PETSc defaults** — used for any field you omit.
3. **`petsc_passthrough_options`** — raw PETSc options applied last. A raw PC
   type must match the backend derived by the preconditioning engine or setup
   fails with an explicit incompatibility error.

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
  (the mathematical Jacobian and preconditioner selections, `SNES reason`,
  Newton iterations, residual evaluations, Krylov iterations, initial/final
  norm, and whether the result was committed or rolled back).

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

@section p55_precond_sec 9. Preconditioning Architecture and Status

The Jacobian interface owns creation, registration, update, naming, and cleanup
of the finite-difference/matrix-free operator. The preconditioner model
interface only describes a matrix structure and inserts interior physical
coefficients. The common engine owns matrix creation/preallocation, repeated
zeroing and assembly, constraint and periodic rows, PETSc backend selection,
alias/ownership tracking, and cleanup.

The optional frozen-momentum/point-block matrix is a separate AIJ matrix. For physical rows
it reproduces the audited same-cell 3x3 block from the reachable legacy mode-2 approximation,
with its sign reversed to match the modern residual convention; for
constraint rows it inserts the exact modern derivative (+1 identity for fixed
rows, or +1/-1 for periodic duplicates). Its viscous diagonal carries the same
effective viscosity the residual diffuses with, `nu + nu_t`, using the residual's own
face average of the eddy viscosity and its wall zeroing; omitting the eddy term left
the matrix modelling a viscous diagonal smaller than the operator's by the
eddy-to-molecular ratio, which on a developed large-eddy simulation is order one or
more. It still intentionally omits pressure, the eddy-viscosity *derivatives* with
respect to velocity, nonorthogonal viscous cross-couplings, the gradient (Clark)
stress, boundary-map and IBM derivatives, and body-force derivatives.

The legacy face inverse Jacobians are arithmetic averages of neighboring cell
inverse Jacobians, and its transverse metric terms average four squared metric-vector
norms. The modern directional fields `IAj/JAj/KAj` directly invert separately
constructed face determinants, while `ICsi/IEta/IZet`, `JCsi/JEta/JZet`, and
`KCsi/KEta/KZet` come from directional face-center constructions. They are not
algebraically equivalent to the legacy expressions on a general curvilinear grid:
in particular, an average of squared norms is not the squared norm of an average.
The explicit legacy formulas are therefore retained here. Directional fields remain
potential inputs to a future, separately specified modern preconditioner model.

The point-block coefficients have an independent test-only legacy transcription
covering nonuniform metrics and velocities, every block entry, both BDF coefficients,
constraint rows, and one- and multi-rank assembly. No performance claim is made.
`PCPBJACOBI` is only the current internal backend mapping; it is not a
user-facing numerical model and is not a historically proven legacy setting.

Future additions are localized as follows: add a Jacobian type/mode
beside `MomentumNewtonJacobian_Create/Update`; add a coefficient provider through
`MomentumPreconditionerModelOps`; add matrix metadata through
`MomentumPreconditionerDescription`; and add a validated structure-to-PETSc
backend mapping in `MomentumPreconditionerEngine_Create`. No placeholder modes
are exposed before their implementations exist.

@section p55_cap_pc_sec 9.1 Preconditioner Model Entries

@htmlinclude generated/capability_inventory_momentum_nk_preconditioner.html

@subsection p55_cap_pc_none_sub none

@anchor p55_cap_pc_none

**Identity.** `momentum_solver.newton_krylov.preconditioner.model: none` (the default) -> `-mom_nk_preconditioner_model none`.

**What it does.** Runs the Krylov solve unpreconditioned on the matrix-free Jacobian. No preconditioning matrix is created, assembled, or stored.

**When to choose it.** The default, and the right starting point: it has no assembly cost, no extra memory, and no approximation of its own to be wrong. Move to `frozen_momentum_jacobian` only once Krylov iteration counts are demonstrably the bottleneck on your case.

**Parameters it owns.** None. Setting `preconditioner.structure` alongside this model is a validation error - `none` accepts no matrix structure, because there is no matrix.

**Interactions.** Mutually exclusive with `frozen_momentum_jacobian`. The Krylov method and its tolerances are configured separately under `linear_solver`.

**Diagnostics.** The Krylov iteration count per Newton step is the signal. It is reported by the SNES/KSP monitors described in @ref p55_monitors_sec; a rising count across timesteps is what motivates preconditioning.

**Evidence.** Integration verified - `make unit-newton-krylov` exercises this path.

**Limitations.** Iteration counts grow with conditioning, so on stiff or highly stretched grids the unpreconditioned solve can dominate the timestep cost.

@subsection p55_cap_pc_frozen_momentum_jacobian_sub frozen_momentum_jacobian

@anchor p55_cap_pc_frozen_momentum_jacobian

**Identity.** `momentum_solver.newton_krylov.preconditioner.model: frozen_momentum_jacobian` with `preconditioner.structure.type: point_block` -> `-mom_nk_preconditioner_model frozen_momentum_jacobian` and `-mom_nk_preconditioner_structure point_block`.

**What it does.** Assembles a separate AIJ matrix holding the same-cell 3x3 momentum block, and uses it as the preconditioning operator for the matrix-free Jacobian. Interior rows reproduce the audited legacy same-cell approximation with its sign reversed to match the modern residual convention; constraint rows carry the exact modern derivative.

**When to choose it.** When Krylov iteration counts under `none` are the measured bottleneck. It buys roughly a 6% wall-clock improvement on the case it was characterised against, while iteration counts on that case rose from 26 to 137 as the problem stiffened - so the honest expectation is a modest gain, not a transformation. Measure before and after; do not enable it on the assumption that preconditioning always helps.

**Parameters it owns.** `preconditioner.structure.type`, which must be `point_block`. It is not an independent choice: the model determines it, and any other value is a validation error.

**Interactions.** Requires `structure.type: point_block` and is rejected without it. Its viscous diagonal carries the residual's own effective viscosity `nu + nu_t`, so a turbulence model is represented at the level the diagonal can represent it. The matrix still deliberately omits pressure, the eddy-viscosity derivatives with respect to velocity, nonorthogonal viscous cross-couplings, the gradient (Clark) stress, boundary-map and IBM derivatives, and body-force derivatives - so its quality degrades as those terms matter more.

**Diagnostics.** Krylov iteration counts before and after are the only meaningful diagnostic. `PCPBJACOBI` appears in PETSc output as the internal backend mapping; it is not a user-facing numerical model and should not be read as one.

**Evidence.** Integration verified - `make unit-newton-krylov`. The point-block coefficients additionally carry an independent test-only legacy transcription covering nonuniform metrics and velocities, every block entry, both BDF coefficients, constraint rows, and one- and multi-rank assembly.

**Limitations.** Experimental, and **no performance claim is made**. It costs an extra assembled matrix in memory and an assembly per update, and the omitted terms mean it is a same-cell approximation rather than an approximate Jacobian in any global sense.

@section p55_validation_sec 10. Validation Coverage

The Newton--Krylov path is covered at two levels:

- **Default suite** (`unit-newton-krylov`, part of `make check`): constraint-row
  Jacobian structure, matrix-free vs direct differencing, preconditioning-engine
  model/backend/ownership wiring, small solve/rollback, and residual repeatability.
  The conservation-outlet conditioned-row derivative
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
