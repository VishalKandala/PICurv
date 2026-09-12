- The frozen-momentum point-block preconditioner was transposed, and had been since it
  was written.
  - `FrozenMomentumJacobian_PointBlock` wrote its six off-diagonal entries the wrong way
    round. The assembler hands PETSc one row at a time and slices `&block[3*component]`,
    so `block` is row-major by construction - `block[3r+c]` is `dF_r/dU_c` - regardless
    of any `MAT_ROW_ORIENTED` setting. The coefficients had been transcribed from the
    legacy `FormJacobian_MF`, which emits the same 3x3 the other way round: it calls
    `MatSetValues(PJ, 3, row, 1, &col, val)` and so fills a column per call. Read back in
    row-major order the legacy block is exactly the transpose of what we stored. The `AJ`
    face factor belongs to the differentiated component - the column - not to the
    equation.
  - The diagonal is unaffected, which is why it survived. Two conditions have to hold
    together before the off-diagonals matter: more than one contravariant component
    significant in the same cell, and a mass term small enough not to swamp them. A run
    started from a uniform streamwise field fails the first, a small `dt` fails the
    second, and the preconditioner looks healthy in exactly the configurations that had
    been tested.
  - Measured on a two-rank bent duct at `Re=40000`, initialized with
    `ucont=(0.15,-0.1,1.0)` so all three components are live from step one. At `dt=0.01`
    (`dtc=150`) there is no signal at all - 30 Krylov iterations before, 29 after. At
    `dt=0.1` (`dtc=15`) the transposed block diverged at step 1 with 490 Krylov
    iterations, worse than running with no preconditioner at all, which reached step 6.
    The corrected block also reaches step 6. This reproduces the `DIVERGED_LINEAR_SOLVE`
    seen on a 193x193x513 cluster job. Rank count does not affect the counts, as the
    point blocks are cell-local.
  - Fixing the orientation restores the preconditioner from "worse than nothing" to "as
    good as its regime allows". It does not make point-block Jacobi strong: at `dt=0.1`
    the limit is the timestep, not the preconditioner, and no preconditioner arm gets
    past step 6. At `dt=0.01` the corrected block holds flat at 5 Krylov iterations per
    step over 400 steps with no divergences. Whether it holds at production grid
    resolution, where the near-wall viscous diagonal is far larger relative to `dtc`, is
    not established by that and needs a cluster run.
  - The engine test now pins the orientation: it sets the three contravariant components
    to distinct values and asserts `dF_i/dU_j` and `dF_j/dU_i` separately, so a transpose
    cannot pass. The previous oracle encoded the same transposition as the implementation
    and agreed with it, which is why the suite was green throughout.
- Newton-Krylov inner-solve observability, since the failure above is invisible in the
  per-step summary.
  - New `<run.runtime_logs>/Momentum_Solver_Newton_Krylov_Linear_History_Block_<n>.log`,
    one row per Krylov iteration, recording the requested relative tolerance beside
    PETSc's reported residual norm. That separates an Eisenstat-Walker tolerance change
    from genuine degradation of the linear solve, which the summary line cannot.
  - Eisenstat-Walker gains a structured `eisenstat_walker` configuration block under
    `momentum_solver.newton_krylov.nonlinear_solver`, mapping its seven fields onto
    PETSc's `-mom_nk_snes_ksp_ew*` options, instead of being reachable only as prefixed
    passthrough.
  - New `monitor.diagnostics.petsc.info`, with optional PETSc class filtering. It has to
    be set during `PetscInitialize`, so it cannot go through `petsc_passthrough_options`
    the way the other PETSc controls do.
- Page 55 said the point block used "the residual's own face average of the eddy
  viscosity and its wall zeroing". The face average is the residual's; the wall zeroing
  is not. `Viscous()` substitutes the wall-model eddy viscosity `lnu_wall` on a wall face
  and falls back to zero only when no wall model is configured. The page and the function
  now say what actually happens, including why the branch currently changes no assembled
  entry - the rows it can affect are exactly the ones `ClassifyMomentumRow()` reports as
  boundary-pinned, which never assemble a block - and why widening the stencil would make
  it reachable and wrong. Tracked as issue #8, to be fixed before any preconditioner with
  a wider stencil.
