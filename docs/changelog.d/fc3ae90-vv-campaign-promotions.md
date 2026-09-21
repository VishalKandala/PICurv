- A local verification campaign measured sixteen subsystems against known answers and
  promoted them to supported: the Picard and Explicit RK4 momentum solvers, the Poisson
  solver, periodic boundaries, domain definition, the Eulerian source, initial
  conditions, particle transport, the grid generator, monitoring, the post-processing
  pipeline, field statistics, spectra, the verification sources, case templates, and
  version management. Twenty-six measurement records carry the evidence. The campaign
  found and fixed the defects below; several change results.
  - Solution-monitoring means and norms counted both ghost layers, so `mean_speed` read
    1.7% and `mean_ke` 14% high on a periodic box. They now cover physical cells only.
  - Multi-rank `Volume` seeding gave each rank an equal count for an unequal subdomain,
    making particle density depend on the decomposition. Ranks now split the count by
    owned cells, and particle IDs are numbered contiguously across ranks.
  - The `Poiseuille` initial condition vanished at the first cell centres instead of on
    the walls and was applied across periodic axes; the `parabolic` inlet likewise
    vanished at the wall-adjacent cells and delivered 27% and 13% too little flux on 8
    and 16 cells. Both now vanish on the walls and are uniform along periodic axes.
  - With `dimensionalize: true`, grid coordinates were scaled once per processed step,
    particle positions and velocities were not scaled, and `Qcrit` carried `U_ref^2`
    rather than `(U_ref/L_ref)^2`. Each quantity is now scaled exactly once.
  - `poisson_solver.method` accepts only `fgmres` and `cg`: `gmres`, `lgmres` and `bcgs`
    converged their Krylov residual while the true residual stalled near 1e-3, and are
    now refused. Per-level `max_it`, `rtol` and `atol` reach the level solvers; they were
    emitted without PETSc's `ksp_` prefix and ignored. `poisson_solver.tolerance`, which
    nothing read, is now refused.
  - Explicit RK4 stops at the step that exceeds its stability limit, naming the limit,
    instead of failing later in the Poisson solve with a message about multigrid depth.
    `make smoke` now runs it on a stable step and on an unstable one.
  - The solver refuses a left-handed grid. It used to negate metric vectors instead,
    which gave a mirrored duct forty times its bulk velocity. The grid generator refuses a
    transform list that ends left-handed, and a new `reverse:axis=i|j|k` transform
    renumbers one logical axis to fix it. Unknown transform, wall and path keys are
    refused instead of silently ignored.
  - `--statistics-state carry` stages the restart bundle even for an analytical
    Eulerian source; it used to fail at start-up. `summarize --plot-spectrum` finds
    spectra that post recipes write into their own subdirectory.
  - The PETSc `-info` log survives the fresh-run log wipe; the reopen no longer leaks
    its filename. A profiling function name that matches no instrumented function
    draws a warning. `versions activate -- -j8` passes `-j8` to `make` instead of
    handing it to Git.
  - Documentation corrections: the convergence modes only log and never judge or stop a
    run; `CornerAveraged` interpolation is below second order even on uniform grids
    (order 1.68 in L2); `TGV3D` is a prescribed kinematic field, not a Navier-Stokes
    solution; the periodic wall-bounded pseudo-time stall no longer reproduces; the IEM
    scalar is inert because nothing seeds a non-zero `Psi`; and an external input
    reference's registration checksum is recorded but not compared, so a changed target
    is used as it now is, under a new asset identity.
  - A supported subsystem may now be demoted to experimental, with a recorded
    `demotion_reason` saying what the earlier claim did not establish.
