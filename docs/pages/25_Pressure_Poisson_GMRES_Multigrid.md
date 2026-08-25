@page 25_Pressure_Poisson_GMRES_Multigrid Pressure-Poisson, GMRES, and Multigrid

@anchor _Pressure_Poisson_GMRES_Multigrid

This page describes the pressure-correction solve path used by projection in PICurv.

@tableofcontents

@section p25_equation_sec 1. Pressure-Correction Equation

The correction solve enforces incompressibility through:

\f[
\nabla^2 \phi = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*.
\f]

In code terms:

- RHS assembly: divergence source in @ref PoissonRHS
- LHS assembly: metric-aware operator in @ref PoissonLHSNew
- solve orchestration: @ref PoissonSolver_MG

Null-space handling is explicitly configured for Neumann-like pressure systems via function @ref PoissonNullSpaceFunction in the Poisson module.

@section p25_mg_sec 2. Multigrid/KSP Stack In Code

@ref PoissonSolver_MG currently:

1. assembles per-level operators,
2. configures `KSP` + `PCMG`,
3. sets restriction/interpolation operators (@ref MyRestriction and @ref MyInterpolation plus solid-aware variants),
4. applies level smoothers/coarse solve,
5. solves finest-level system for `Phi`.

After Poisson solve:

- pressure is updated by @ref UpdatePressure
- velocity is projected by @ref Projection

@section p25_config_sec 3. YAML Mapping and PETSc Options

From `solver.yml` via `picurv_cli/core.py`:

- `poisson_solver.method` -> `-ps_ksp_type`
- `poisson_solver.absolute_tolerance` -> `-ps_ksp_atol` and legacy `-poisson_tol`
- `poisson_solver.relative_tolerance` -> `-ps_ksp_rtol`
- `poisson_solver.max_iterations` -> `-ps_ksp_max_it`
- `poisson_solver.gmres.restart` -> `-ps_ksp_gmres_restart`
- `poisson_solver.preconditioner.type` -> `-ps_pc_type`; only `multigrid` is supported today
- `poisson_solver.multigrid.levels` -> `-mg_level`
- `poisson_solver.multigrid.pre_sweeps` -> `-mg_pre_it`
- `poisson_solver.multigrid.post_sweeps` -> `-mg_post_it`
- `poisson_solver.multigrid.semi_coarsening.{i,j,k}` -> `-mg_i_semi`, `-mg_j_semi`, `-mg_k_semi`
- `poisson_solver.multigrid.level_solvers.level_N.method` -> `-ps_mg_levels_N_ksp_type`
- `poisson_solver.multigrid.level_solvers.level_N.preconditioner` -> `-ps_mg_levels_N_pc_type`
- `poisson_solver.multigrid.cycle` and `.mode` are accepted structured keys; current supported values are `v` and `multiplicative`
- `pressure_solver` remains a legacy alias for `poisson_solver`
- `petsc_passthrough_options` -> advanced PETSc flags not exposed as structured YAML

Final option parsing happens in function @ref CreateSimulationContext during context creation.

MG level numbering follows PETSc/PICurv convention: `level_0` is the coarsest grid.
Larger level numbers are progressively finer. The default MG level preconditioner
is block Jacobi (`bjacobi`) when not specified. The current PETSc binding applies
one MG smoother count; if `pre_sweeps` and `post_sweeps` differ, PICurv uses the
larger value and logs a warning.

Common MG-level preconditioner notes:

- `jacobi` and `sor` are simple smoother PCs. SOR can use
  `-ps_mg_levels_N_pc_sor_omega <positive-real>` through passthrough when needed.
- `ilu` and `lu` are factor PCs. Useful passthrough knobs include
  `-ps_mg_levels_N_pc_factor_levels <nonnegative-integer>` and
  `-ps_mg_levels_N_pc_factor_shift_amount <nonnegative-real>`.
- `bjacobi` owns nested block solves; inspect exact nested PETSc prefixes with
  `-ps_ksp_view` before tuning sub-KSP/sub-PC options for a specific PETSc build.

@section p25_robustness_sec 4. Robustness Characteristics

Current implementation includes:

- periodic-boundary pressure synchronization,
- immersed-boundary-aware treatment paths (`Nvert`/solid checks),
- optional Poisson monitor logging via `monitor.yml -> solver_monitoring.poisson.pic_true_residual`.

If pressure solve quality degrades, check first:

1. BC consistency,
2. MG level/smoother settings,
3. grid metrics/orientation quality,
4. solver tolerances vs timestep.

## The V-cycle has two kinds of solver, and only one of them is a smoother

A V-cycle attacks the error at two different scales, with two different tools.

On every level above the coarsest, a **smoother** runs a fixed number of
relaxation sweeps. Relaxation is very good at removing error components whose
wavelength is comparable to the local mesh spacing, and almost useless against
error that is smooth on that mesh. That is the whole design: each level strips
out the high-frequency error it can see, restricts what is left to a coarser
mesh where the remaining error looks high-frequency again, and repeats.

At the base of the cycle sits the **coarse solve**. By construction it receives
the error that every smoother above it was blind to: the smooth, long-wavelength,
global component. Nothing below it will get another chance at that error, so the
coarse solve is expected to remove it essentially exactly, in one shot.

These are categorically different jobs, and PICurv's level naming actively hides
the distinction:

| YAML key | PETSc option prefix | Role |
|---|---|---|
| `level_solvers.level_0` | `-ps_mg_coarse_` | **Coarse solve** at the base of the V-cycle |
| `level_solvers.level_1` .. `level_N` | `-ps_mg_levels_N_` | **Smoothers**, coarse to fine |

`level_0` looks like just another entry in an evenly spaced list. It is not. It
is the one entry that is not a smoother, and configuring it as though it were is
the root of the defect described next.

## The coarse solve must be a fixed linear operator

Multigrid here is a **preconditioner**, not a standalone solver. The outer
Krylov method's convergence theory assumes the preconditioner is a fixed linear
operator: feeding it the same vector must always return the same vector.

A Krylov method violates that. Its iteration builds a subspace tailored to the
vector it was given, and its stopping test fires after however many iterations
that particular vector needs. Two different inputs get two different numbers of
inner iterations, so the map from input to output is not linear and not even
fixed. Put a Krylov method at `level_0` and the entire multigrid preconditioner
becomes a nonlinear operator.

**Pinning `ksp_max_it` does not fix this.** A fixed iteration count still builds
a Krylov space out of the input vector, so the operator still depends on its
input in a nonlinear way. This is exactly why the classical smoothers -
Richardson, Jacobi, Chebyshev, SOR - are the right tools for the smoothing
levels: each is a fixed linear operator, applied a fixed number of times.

The outer method PICurv uses, `fgmres`, is *flexible*: it tolerates a varying
preconditioner well enough to keep constructing a solution. What it cannot do is
keep its Arnoldi recurrence consistent with the true residual. FGMRES uses right
preconditioning, so under a fixed preconditioner its recurrence-tracked residual
equals the true residual `b - Ax`. Under a varying one, the two quantities
separate - and the convergence test reads the tracked one. The solver reports
convergence against a number that no longer describes the constraint it exists
to enforce.

### Choosing a coarse solver by coarse-grid size

| Coarse-grid unknowns | Recommended `level_0` | Why |
|---|---|---|
| up to ~1e4 | `{method: preonly, preconditioner: redundant}` | Every rank forms and factors the whole coarse operator with LU. Exact, fixed, linear, no iteration to vary. Cheap because the grid is tiny. |
| larger | `{method: preonly, preconditioner: telescope}` | Redundant LU on every rank stops being cheap. `PCTELESCOPE` moves the coarse solve onto a subset of ranks and does a direct solve there, still fixed and linear. |
| last resort | a Krylov method | Only when the coarse grid is genuinely too large to factor. Then set tolerances against the **true** residual, keep `pic_true_residual` on, and treat the reported convergence as unverified until you have checked the two norms agree. |

PICurv logs a startup warning when a Krylov `ksp_type` is configured at
`level_0`. It stays a warning rather than an error because the last-resort case
above is legitimate at large scale.

### How many levels

Two constraints bound `multigrid.levels` from opposite ends.

From below, aim for a coarsest grid of roughly **1e3 to 1e4 unknowns**. Smaller
wastes a level; larger makes the replicated LU factor expensive on every rank.
Because each level removes a factor of eight from a 3-D grid, this scales
logarithmically: a grid of ~5M cells wants 5 levels, not 4.

From above, coarsenability. Each level halves an axis as `IM -> (IM+1)/2`, so
`IM` must stay **odd at every level** for the coarsening to be exact. The chain
runs `IM_fine = 2 * IM_coarse - 1`, giving usable ladders such as

    5 -> 9 -> 17 -> 33 -> 65 -> 129 -> 257
    4 -> 7 -> 13 -> 25 -> 49 -> 97

An even count still runs, but logs

```
WARNING: Grid at level L, block B can't be consistently coarsened further.
```

and proceeds on a slightly misaligned coarse grid. The MPI bound described in
the next section applies on top of both of these.

### Scaling

`redundant` scales indefinitely **provided levels are added as the grid is
refined**. Keeping the level count fixed while refining grows the coarse grid
with the fine one, and the replicated LU factor eventually dominates. Adding a
level as the grid grows keeps the coarse grid, and therefore the coarse solve,
roughly constant.

`richardson` + `bjacobi` at `level_0` is a fixed linear operator and so is
*correct*, but it is not an exact coarse solve. As the coarse grid grows,
Richardson leaves more and more of the smooth error behind and the V-cycle stops
being mesh-independent: iteration counts creep up with resolution. It is a
reasonable stopgap, not a scalable answer.

## Diagnosing tracked-vs-true residual divergence

`DualKSPMonitor` (`src/logging.c`) prints both numbers to
`logs/Poisson_Solver_Convergence_History_Block_0.log` when
`monitor.yml -> solver_monitoring.poisson.pic_true_residual` is on:

- **`Unprecond Norm`** is PETSc's `rnorm`, carried by the Krylov recurrence.
  This is what the convergence test reads.
- **`True Norm`** is an explicit `KSPBuildResidual` recomputation of `b - Ax`.

Under a correct configuration the two agree to many digits at every iteration.
When they separate, the convergence test is passing on a fiction.

Worked example. A 32x32x192 curved, wall-clustered grid at Re = 40,000 on 8 MPI
ranks. Identical case and solver files; only `level_solvers.level_0` differs:

| `level_0` setting | Tracked residual | Recomputed `b-Ax` | Max divergence | Iterations |
|---|---|---|---|---|
| `{fgmres, bjacobi}` | 4.80e-10 | **1.52e-05** | **1.02e-08** | 14-16 |
| `{richardson, bjacobi}` | 1.67898e-09 | 1.67898e-09 | 5.18e-12 | 11 |
| `{preonly, redundant}` | 8.62828e-10 | 8.62848e-10 | 1.97e-12 | 10-11 |

Read the first row carefully. The tracked residual is the smallest of the three,
which is why the defect survived: by the number the solver reports, the Krylov
coarse solve looked *best*. The true residual is five orders of magnitude
larger, and the physical consequence follows directly - maximum divergence is
1e-08 rather than 1e-12. The incompressibility constraint the Poisson solve
exists to enforce was being violated by six orders of magnitude, silently. The
fixed-operator settings are simultaneously more accurate and cheaper, converging
in fewer iterations.

The failure was **rank-dependent**: 6 and 10 ranks were clean, 4 and 8 were
broken, with the same configuration. That follows from `bjacobi`, whose block
structure is inherited from the DMDA decomposition, so the coarse
preconditioner - and hence how nonlinear the coarse solve behaves - changes with
the rank layout. A configuration that looks fine on your development rank count
can be wrong on the production one.

Diagnostic sequence:

1. Turn on `pic_true_residual` and compare the two norms in the Poisson log.
2. If they disagree, look at `level_solvers.level_0` first.
3. Confirm the physical symptom in `logs/Continuity_Metrics.log`: a decoupled
   residual shows up as a max divergence orders of magnitude above the solver
   tolerance you asked for.
4. Re-run at a different rank count before concluding a configuration is clean.

`tests/smoke/run_driven_periodic_regression.sh` asserts this invariant at 4 and
10 ranks - two counts that previously disagreed.

## Multigrid depth is bounded by the MPI decomposition

`multigrid.levels` cannot be chosen independently of the rank layout. Every
level's `DMDA` must leave each rank at least `stencil_width` nodes along every
axis, and grid setup requests a stencil width of **3 whenever any axis is
periodic** (2 otherwise). Coarsening halves the node count per level, so a deep
hierarchy spread over many ranks eventually starves the coarsest grid and PETSc
aborts during DM creation, before the first timestep:

```
PETSC ERROR: Argument out of range
PETSC ERROR: Local x-width of domain x 2 is smaller than stencil width s 3
```

With `levels: L` and `N` cells along an axis, the coarsest grid holds
`N / 2^(L-1)` cells, hence `M = N / 2^(L-1) + 1` nodes. PETSc distributes those
over `P` ranks so the smallest rank holds `floor(M / P)`. The constraint is:

    floor( (N / 2^(L-1) + 1) / P ) >= stencil_width

Worked maxima for a triply periodic box (stencil width 3):

| Cells per axis | Ranks per axis | Total ranks | Max `levels` |
|---|---|---|---|
| 64 | 2 | 8 | 4 |
| 64 | 4 | 64 | 3 |
| 128 | 4 | 64 | 4 |
| 128 | 8 | 512 | 3 |
| 192 | 6 | 216 | 4 |

Two consequences worth planning around:

- Periodic cases are strictly more constrained than wall-bounded ones, because
  the stencil width is 3 rather than 2. A solver profile that runs for a channel
  can abort when reused for a fully periodic box at the same rank count.
- `picurv validate` cannot catch this. It does not see the runtime MPI
  decomposition, so the first evidence is the aborted job.

Reducing `levels` is the usual fix; also delete the now-unused
`level_solvers.level_N` entry for the level that no longer exists. Shallower
hierarchies converge more slowly but stay usable. Measure the cost in
`logs/Poisson_Solver_Convergence_History_Block_0.log`: if iteration counts grow
sharply, prefer lowering the rank count and restoring the deeper hierarchy.

@section p25_testing_sec 5. Current test status

Current direct tests are strongest for helper and invariant behavior:

- `PoissonLHSNew`
- `Projection`
- `PoissonNullSpaceFunction`
- RHS-related helpers used by `ComputeRHS`

The main remaining gap is `PoissonSolver_MG`: it is exercised in runtime smoke, but still lacks equivalent direct bespoke coverage for debugging. Periodic and immersed-boundary stencil branches also remain thinner than the core Cartesian helper surface.

@section p25_refs_sec 6. Related Pages

- **@subpage 23_Fractional_Step_Method**
- **@subpage 24_Dual_Time_Picard_Jameson_RK**
- **@subpage 31_Momentum_Solvers**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Pressure-Poisson, GMRES, and Multigrid** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
