@page 18_Changelog Changelog

@anchor _Changelog

# Changelog

## Unreleased

- Wall-model observability, and a defect that made the wall model unusable.
  - `lFriction_Velocity` was created from `fda`, whose dof is 3, while every reader
    opened it through `da`, whose dof is 1. PETSc rejects that pairing outright, so
    enabling wall functions on a case with a `WALL` face aborted at the first boundary
    pass with `Vector local size N is not compatible with DMDA local sizes N/3`. The unit
    fixtures allocated it correctly and so never saw it. Storage is now created on `da`
    as a global/local pair, and `tests/c/test_boundaries.c` allocates the way production
    does and checks the allocation against the field catalog descriptor, so a DM or dof
    mismatch fails in the suite rather than only in a run.
  - The friction velocity was computed every step and discarded: no field, no column, no
    checkpoint. It is now the catalogued field `Utau`, written when a wall model is
    active. It is not read back on restart, because the first boundary pass of the
    restarted run recomputes it from the restored velocity.
  - New `<run.runtime_logs>/wall_model.csv`, one row per step: the friction velocity's
    mean, RMS and extrema over the corrected cells, the first-cell `y+` and wall
    distance, and the cell count. `y+` is formed inside the wall-model dispatch, where
    the wall distance is still in hand; nothing downstream can recover it. A step that
    corrects no cell warns and writes nothing rather than emitting zeros. Plottable as
    `picurv summarize --plot wall_model.y_plus_mean`.
  - The banner no longer prints a roughness height beside Werner-Wengle or Cabot, which
    discard it, and now reports the Yoshizawa coefficient with the diagnostics line.
- The wall model's stress now reaches the momentum equation. The correction sets a
  near-wall velocity, and the viscous operator reaches the wall through a viscosity times
  a gradient - so with the subgrid viscosity zeroed at wall faces, only the molecular
  fraction `u+/y+` of the modelled stress was delivered, about a fourteenth of it at
  `y+ = 267`. That zeroing is right for a wall-resolved run and backwards for a
  wall-modelled one. The wall pass now publishes an effective wall eddy viscosity,
  `nu_eff = tau_w y / u`, formed where the stress, the distance and the corrected speed
  are all still in hand, and the viscous flux uses it at that face. Reported as
  `nu_wall_over_nu_mean`; a 30-step channel shows the solutions with and without it
  separating monotonically.
- Three wall-model pairings are now rejected at validation, each naming why: a wall model
  with no turbulence model (there is no implicit-LES scheme here to stand in - the
  convection is QUICK, whose dissipation is not a subgrid model), and `cabot` or `werner`
  under RANS. A wall model on a laminar case is refused for the same reason from the other
  side. And because whether the first cell lands in the law's valid range depends on the
  mesh rather than the configuration, a `y+` excursion warns on every diagnostic sample
  and stops the run after ten consecutive ones.
- Logging is now a registered subsystem, `observability.logging`, at `supported`. It had
  been in no subsystem record, no capability family, no contract, no freshness surface,
  and in neither `src/guide.md` nor the module map of `13_Code_Architecture` - the shape
  that leaves a module unowned is precisely that everything reports through it and nothing
  depends on it for a result. A change to what a user sees therefore tripped no staleness
  suspicion anywhere. Both module maps now carry it, `09_Monitor_Reference` documents the
  append-and-continuation-marker lifecycle its runtime files follow across a restart, and
  a soft freshness surface watches `src/logging.c` and `include/logging.h`. The supported
  claim is about its behavioural contracts, which `tests/c/test_logging.c` covers; it is
  not a claim about any quantity another subsystem hands it to print.
- A toy channel exercising all three wall laws end to end found two defects that unit
  coverage could not reach, both now fixed.
  - Werner-Wengle reported one identical `u_tau` at every wall cell: `wall_function()`
    used its Newton initial guess as the answer. It now uses a closed-form
    `taw_Werner()`/`u_Werner_explicit()` pair - the profile integrated across the first
    cell, which is the Werner-Wengle model proper and inverts without any root-find, so
    the no-inner-iteration property the model is selected for is now real rather than
    claimed. The iterative `f_Werner`/`df_Werner`/`find_utau_Werner` helpers are retained
    but unwired, with the reasons recorded at their definitions: their branch threshold
    is written through the friction velocity being solved for, and their two branches
    invert different relations, so they reconstruct the pointwise `u_Werner()` only below
    `y+ = 11.81` - which is the only range their unit test samples.
  - The friction-velocity storage was allocated inside the `les || rans` block, so
    enabling a wall model without either of them left `ApplyWallFunction` dereferencing a
    null vector. `wall_function` is their sibling in the schema, not their child;
    allocation is now gated on the wall model alone.
  All three laws now vary cell to cell and differ at `y+ ~ 10-14`, and agree exactly at
  `y+ ~ 1.4`, where every wall model must reduce to `u+ = y+`.
- `les.clipping.max_cs` under a mode that imposes no ceiling is now an error rather than
  a warning. A ceiling that is not in force reads as one that is.
- The Newton-Krylov preconditioner's omission of the Clark stress is now documented as
  the deliberate choice it is - the term is higher order and not diagonally dominant, so
  representing it would change what kind of matrix the point block is - rather than as an
  unexplained gap. The Krylov cost of that choice has not been measured.
- Stale flag mappings in `07_Case_Reference` corrected: `-const_cs`, `-max_cs`,
  `-dynamic_freq` and `-testfilter_ik` no longer exist, and `les.max_cs` is not an
  accepted key. The list now matches what the CLI emits.
- `tests/c/test_les.c` was never listed in `tests/guide.md`, and the `2026-03-20` coverage
  snapshot for `src/les.c` described a file that has since been rewritten. Both corrected;
  the assembled LES model path has not been re-measured.

- Wall function boundary treatment: integration coverage, its own lifecycle record,
  and two corrected claims.
  - `ApplyWallFunction()` had no test. The velocity laws and their friction-velocity
    root-finds were already covered, but the wiring around them was not, and its
    plausible failures are silent: feeding the reference and boundary wall distances in
    the wrong order, or taking the reference velocity from the wrong cell, still
    produces a plausible corrected velocity. `tests/c/test_boundaries.c` now checks that
    the corrected first interior velocity is the modelled profile at that cell's wall
    distance for the friction velocity the integration stored, and separately that the
    stored friction velocity inverts the law at the *reference* cell's distance. Swapping
    the two distances in production is detected.
  - The startup banner now reports the wall-function model and roughness height.
    `07_Case_Reference` claimed it already did; `io.c` contained no such line. The claim
    is made true rather than removed, and the line is pinned by the user-facing
    reporting contract.
  - `turbulence.wall_function` has its own subsystem record. It was owned by the
    `turbulence.rans` record, which is the wrong home: a wall function is configured
    independently of both LES and RANS, and the RANS runtime update is incomplete while
    this treatment is separately usable.
  - The model selector now reaches the runtime. `-wallfunction` previously carried
    only an enable flag: the configured `model` was validated by the Python layer and
    then discarded, so every case ran the log law whatever it asked for, and
    `ApplyWallFunction()` dispatched to it unconditionally at all six faces. The flag
    now carries a `WallFunctionModel` code, one dispatcher replaces the six hardcoded
    call sites, and `werner` and `cabot` join `log_law` as selectable. Both were
    already implemented and unit-tested in `src/wallfunction.c` and had simply never
    been reachable. `WallFunctionModelToString()` reports the choice in the banner.
    Cabot's non-equilibrium pressure-gradient term is supplied as zero, reducing it to
    the equilibrium form; that is documented on its capability entry rather than
    implied.
  - Wall-overwritten cells are excluded from the dynamic procedure. The correction is
    an imposed analytic profile, so the Germano identity evaluated there measures that
    profile rather than resolved turbulence. No averaging mode escapes it: the biased
    set is the whole domain under `global`, the wall plane under `homogeneous` retaining
    the wall-normal direction, and that cell alone under `local`. A cell left with no
    data takes a zero coefficient, which is what a wall model already supplying the
    stress calls for. Carried through the averaging kernel's inclusion mask, which was
    already there for the statistics pipeline.
  - `PicurvSpatialRatioAverage()` now honours the target mask and the inclusion mask in
    its pointwise branch as well as its averaged one. It previously applied neither when
    no direction was selected, so an exclusion silently changed meaning depending on the
    averaging mode.
  - Cabot receives the resolved wall-parallel pressure gradient, taken from the
    near-wall pressure through `ComputeScalarFieldDerivatives()`. Supplying zero would
    have reduced it to its equilibrium form - not the model the selector names.
  - `roughness_height` is rejected under `werner` and `cabot`. Neither has a roughness
    formulation - Werner-Wengle has no such term and Cabot discards the argument - so
    accepting the key would silently drop a value set on purpose.
  - **The wall correction now reaches the LES closure.** `Contra2Cart` rebuilds the
    Cartesian field at the top of each timestep and discarded the previous solve's
    near-wall correction, so the eddy viscosity was computed from the uncorrected field
    while the momentum equation used the corrected one - the two halves of one timestep
    disagreeing about the velocity. The wall model is re-applied before the strain rates
    are formed. **This changes results for any case combining LES with a wall
    function.** The remaining caveat, that an imposed profile contaminates the dynamic
    coefficient at the wall-adjacent cell, is documented on the closure page.

- LES subgrid closure corrected and consolidated. **Results from earlier revisions of
  this tree are not reproducible under it, and neither LES model previously applied the
  dissipation it claimed to.**
  - The dynamic procedure now evaluates the Germano identity. `M_ij` previously built
    both of its terms from separately test-filtered quantities, so they shared an
    identical factor and the tensor collapsed to a fixed multiple of the filtered strain
    rate. The scale-similarity relation was therefore never evaluated: the model had no
    dynamic content, while still producing a positive, dimensionally correct,
    flow-responsive coefficient — which is why the defect survived. The second term is
    now the test filter of the product `|S| S_ij`, which is what the identity requires,
    and `M_ij` is projected trace-free so the Leonard stress's trace cannot enter the
    contraction through discrete divergence error.
  - The constant model no longer allocates a coefficient field. It previously filled a
    whole DMDA vector with one configured number, synchronized it, and checkpointed it —
    and on a fresh run never filled it at all, because its only call site hit the
    early-return that holds the coefficient at zero for the first steps. It now reads
    `constant_cs` directly each step, which both fixes the defect and removes a global
    vector, a local vector, a ghost exchange, and a periodic synchronization from that
    path.
  - Coefficient averaging is implemented, replacing a comment that said the logic
    "would go here". `averaging.mode` selects `local`, `homogeneous`, or `global`;
    homogeneous derives its directions from the case's periodic boundary pairs unless
    they are named, so a triply periodic box yields one coefficient and a channel yields
    a wall-normal profile with no extra configuration. Sums are volume-weighted, exclude
    solid cells and the periodic duplicate planes, and are reduced over the block's
    communicator, so the result does not depend on the decomposition.
  - Limiting is selectable. `clipping.mode: none` keeps the signed coefficient so
    backscatter survives, bounded instead by a total-viscosity floor that constrains the
    quantity which actually has to stay positive. The clipping modes remain available and
    remain biased: discarding only the negative tail raises the mean dissipation. The
    `max_cs` default drops from 0.5 to 0.3.
  - The grid filter width and the test-to-grid width ratio are exposed rather than
    hardcoded, and the dead `ComputeCellCharacteristicLengthScale` call that passed one
    output pointer twice is gone.
  - The solid mask is unified on the repository-wide `nvert < 0.1` fluid test; the
    dynamic pass previously used a different threshold from the eddy-viscosity pass.
  - `CS` now carries `FIELD_AVAILABILITY_LES_DYNAMIC` and holds the coefficient `C`
    that multiplies `Delta^2 |S|` — `Cs^2` in the classical notation, signed when
    backscatter is admitted. Restarting from a checkpoint written before this change is
    not supported for that field.
  - The Clark gradient term, live in the viscous flux but reachable only through PETSc
    passthrough, is now `les.gradient_model.enabled`. Combined with a Smagorinsky closure
    it is the mixed model.
  - Turbulence option ingress is namespaced under `-les_*` and grouped into one
    `LESConfig` block, with defaults defined once in `LESConfigSetDefaults()` so the
    control-file path and the test fixtures cannot drift. The dead `-mixed`,
    `-testfilter_1d`, and `-i/j/k_homo_filter` options are removed; the homogeneous-filter
    intent they carried is served by `averaging.directions`.
  - Per-step diagnostics are written to `<run.runtime_logs>/les_coefficient.csv`,
    including the pre-clipping backscattering and limited fractions, which no stored
    field preserves.
  - `les.c` is rebuilt from named kernels — symmetric-tensor algebra, strain rate, filter
    width, the two Germano tensors, averaging, limiting, and viscosity assembly — each
    unit-tested in the new `tests/c/test_les.c` (`make unit-les`), with a
    decomposition-independence case in `tests/c/test_mpi_kernels.c`.
  - The subsystem moves from `known-defective` to `experimental`, and both capability
    scope records are retired. It is not `supported`: the coefficient magnitude has not
    been validated against a reference flow. The check that would close that gap is
    `examples/decaying_isotropic_turbulence` with homogeneous averaging, where `Cs(t)`
    should settle near 0.16-0.17.
  - New page **72_LES_Turbulence_Closure** derives the formulation and documents the
    averaging, limiting, and diagnostic choices.
  - Follow-up consolidation, so the closure shares rather than restates:
    - the closure resolves its iteration domain through `SpatialTargetPlanCreate()`
      and its solid mask through `SpatialTargetPlanMaskAllows()`, retiring a
      local range helper and a duplicate fluid threshold;
    - the coefficient averaging moved out of `les.c` to
      `PicurvSpatialRatioAverage()` in `statistics_target.*`, generalised with an
      optional denominator, an optional inclusion mask, and a caller-chosen
      communicator. `PicurvWindowSpatialMean()` is now a thin caller of it; the
      statistics pipeline keeps count weighting and its existing communicator, so
      its output is unchanged. Volume weighting stays the LES caller's choice and
      is folded into the two contraction fields before the reduction;
    - per-axis periodicity comes from the resolved `SimCtx` flags rather than a
      second read of the boundary faces. The loader rejects a case whose opposite
      faces or whose blocks disagree, so the two cannot differ, and the flags are
      what both the DMDA and the spatial target were already built from;
    - the Clark gradient term in `rhs.c` uses `SymTensor` instead of nine
      hand-written components repeated in three places;
    - the six-component accumulator order in `statistics_accumulator.c` is tied to
      `SymTensor`'s member order by a compile-time assertion, so the two cannot
      drift apart unnoticed.
  - Newton-Krylov reach, so the closure is not confined to one momentum solver family:
    - the preconditioner's viscous diagonal now carries the residual's own effective
      viscosity `nu + nu_t`, using the same face average and wall zeroing `Viscous()`
      applies. Omitting the eddy term left the matrix modelling a viscous diagonal
      smaller than the operator's by the eddy-to-molecular ratio, which on a developed
      large-eddy simulation is order one or more; a plausible contributor to the
      recorded momentum-preconditioner bottleneck;
      The point-block oracle in `tests/c/test_momentum_newton_krylov.c` carries the
      term as an independent transcription from the residual's definition, with mutants
      for omitting it and for reading the cell value instead of the face average; both
      production mutations are detected;
    - the gradient (Clark) model and wall functions are no longer rejected. Both reach
      the solver through the shared residual it already evaluates, so neither needed
      solver-specific work. Two properties to know: the preconditioner does not model
      the Clark stress, so Krylov counts can rise with it enabled, and the wall
      function's iterative friction-velocity solve interacts with matrix-free
      differencing, whose `1/h` amplification surfaces a loose inner tolerance as a
      degraded Jacobian action;
    - solid cells no longer disqualify the solver. `ComputeRHS()` already zeroes their
      rows, which leaves a solver that marches on the residual nothing to do but leaves
      a solver that assembles a matrix with an undetermined unknown.
      `MomentumRowIsSolidMasked()` gives those rows the same identity treatment as any
      other row without an unknown. The immersed-boundary method stays unsupported,
      because its velocity reconstruction does not run inside this solver's residual.

  - Exposure sweep across the user-facing surfaces:
    - `les.diagnostics.yoshizawa_ci` is now a documented key. It was parsed by the C
      runtime but had no YAML spelling and no emitter, so it was reachable only
      through PETSc passthrough — the same hidden-control defect the Clark term had;
    - the startup banner reports the resolved closure: filter width, test filter and
      width ratio, dynamic cadence, averaging mode and directions, limiting mode and
      ceiling, viscosity floor, gradient term, and diagnostics state. Previously it
      printed the model name alone, so nothing confirmed a setting had landed. Three
      of the new lines are pinned by the user-facing reporting contract;
    - `picurv summarize` reads `les_coefficient.csv`, so `--list-plot-series` offers
      the `les.*` histories and `--plot les.cs_effective` renders the coefficient
      curve an LES run is judged on;
    - new LES recipe in the cookbook, troubleshooting section keyed to the diagnostic
      columns, and glossary entries for the closure vocabulary. The closure page is
      now reachable from both the user and developer index pages.

- Stopped the dual-time pseudo-CFL controller re-discovering the same stability limit, and
  stopped a rejected trial voting on the next decision.
  - The EMA-smoothed residual ratio is now committed only when its trial is accepted. A
    rejected trial is rolled back, so carrying its ratio forward was a state leak. Measured on
    the laminar channel: a diverging trial at cfl 2.0 (raw ratio 3.97) pushed the EMA to 1.61
    and the next trial at cfl 1.5 was rejected too, despite a raw ratio of 0.994 -- the
    residual had decreased. That false rejection discarded four ComputeRHS calls and drove the
    pseudo-CFL to 1.125, below the measured optimum.
  - Failed pseudo-CFLs are now remembered. Previously the controller climbed 1.36 -> 2.00,
    rejected, recovered and climbed straight back on every physical step, spending 28.8% of all
    trials above the measured stability limit. A dimensionless cap tightens to 0.60x a failed
    CFL and relaxes by 1.005 per accepted trial so a transient rejection is not permanent.
  - Rejections roughly halve and total pseudo-iterations fall 3-9% at the shipped ceiling:
    laminar 17x33x17 233 -> 226 trials (13.3% -> 4.9% rejected), plane channel 33^3 340 -> 308
    (20.6% -> 6.5%), flat channel 25x25x97 155 -> 151 (7.7% -> 5.3%). No case lost convergence.
  - The remaining cost is not rejection handling. Measured stability limits are 1.12 (33^3),
    1.25 (flat), 1.82 (laminar), 2.25 (LES retau180) and 2.50 (DIT) in pseudo-CFL, and the
    fastest-converging CFL sits at 0.67-0.80 of that limit -- the Jameson scheme is a smoother
    and smoothers damp worst near their limit. Operating at the optimum is worth 32-52%, which
    a stability-triggered cap cannot reach: the LES case rejects nothing before or after and is
    therefore unaffected.
  - Recorded for the record: promoting the shadow operator-scaled stability estimate
    (ComputeMomentumStabilityEstimate) was evaluated against these measurements and rejected.
    It is a *worse* normaliser of the measured stability limit than the estimate it would
    replace -- 3.20x spread across the five cases versus 2.23x -- so it stays in shadow.

- Made the momentum residual, not the velocity update, the criterion that decides
  dual-time convergence, and gave the absolute residual tolerance a portable scale.
  - `converged = residual_abs_pass OR (residual_rel_pass AND update_pass)`. The
    absolute residual test is sufficient on its own; the relative test keeps the
    `relative_tol` update guard. Requiring the guard alongside the absolute test made
    the latter unable to fire: on the laminar channel at step 793, `|R|` fell below
    the floor at pseudo-iteration 8 and the step still ran to 20.
  - `residual_absolute_tol` is now **dimensionless**, tested as
    `|R| <= tol * resid_ref` with `resid_ref = a0*|Ucont|_inf/dt`. A raw bound on `|R|`
    is not portable: across the shipped cases step-1 `|R|` spans 1.2e-2 (plane channels)
    to 2.0 (driven duct), a ~165x spread that normalising collapses to ~4x. Set to 1.0e-8
    in all shipped configs.
  - `absolute_tol` / `-mom_atol` no longer participates while a residual tolerance is
    set. `|dU| ~ dtau*|R|`, so bounding it absolutely is a disguised residual bound
    `|R| <= absolute_tol/dtau` that tightens as the controller grows `dtau`; it was in
    practice the criterion that gated convergence. It is retained, annotated, for the
    legacy update-only branch only. It is now **deprecated**: removed from all shipped
    configs (commented, with a note, in the master template), still accepted, and the CLI
    warns when it is set while it cannot take effect.
  - Residual-based convergence is now the **default** (`-mom_resid_atol` 1e-8,
    `-mom_resid_rtol` 1e-3, previously both 0.0). The update-only branch it replaces can
    converge falsely, since `|dU|` also goes small when `dtau` collapses. Setting both
    non-positive remains an explicit opt-out.
  - Verified on the 6000-step driven periodic laminar channel against its exact solution:
    identical `U_b` 0.999446, `u_max` 1.497795, `u_tau` 0.1729983 and profile error
    1.470e-03 before and after; fields agree to 5.06e-13 at step 6000 and the gap between
    the two runs *shrinks* from 2.2e-07 (step 500) to 3.4e-13, so error does not
    accumulate. Worst `|div u|` unchanged (7.457e-11 vs 7.443e-11). Total pseudo-iterations
    fell 28,247 -> 16,390 (-42%), concentrated in the settled regime (3% over steps 1-500,
    ~50% thereafter); a flow that never settles should expect the smaller figure.

- Restored future-architecture documentation and added the authoritative
  proposed field-statistics specification, including centered moment state,
  layout/mask rules, future spatial targets, restart behavior, postprocessing
  responsibilities, legacy removal, and phased implementation plan. Also
  restored the separate benchmark-gated function-identity specification.

- Field-statistics foundation (Phase 1 implementation):
  - moved physical-solution convergence configuration to
    `monitor.yml -> solution_monitoring.convergence`, preserving its existing
    flags, every-completed-step behavior, and log format while adding an
    effective enable switch;
  - kept the generated master `control` as the single C-ingress artifact and
    removed the premature observation sidecar, schema/version envelope,
    exposed cadence, and unsupported-window startup gate;
  - kept unfinished scientific field-statistics keys out of active user YAML
    until accumulation, checkpointing, and postprocessing work end to end;
  - removed `SimCtx::averaging`, legacy sum vectors, `su0`/`su1`/`su2`/`sp`
    read/write hooks, the dead averaging call, old templates, and the reserved
    averaged-field postprocessor passthrough; and
  - replaced flat restart files with immutable committed checkpoint bundles
    containing a text manifest, manifest SHA-256 commit marker, geometry/layout
    identity, catalog-selected PETSc-vector payloads, optional particle/model
    state, and the previous contravariant velocity needed for exact BDF2 restart;
  - routed solve output, continuation, restart staging, live post discovery, and
    final/graceful-shutdown output through the bundle contract while retaining
    the existing generic vector and swarm readers/writers;
  - kept solver/convergence/profiling and other persistent logs independent of
    checkpoint state; and
  - updated Python/C regression coverage and ingress/runtime/reference
    documentation. Scientific statistics accumulation remains Phase 2 work.

- Completed the field-statistics window stage with independent per-window PETSc
  accumulator state: `include/statistics_accumulator.h` and
  `src/statistics_accumulator.c` allocate one accumulator set per window from the
  vector factory, duplicating every vector from one the factory already built, and
  apply an accepted state pointwise through the centered kernels. A three-vector's
  second moment is the six symmetric co-moments between component pairs in fixed
  (xx, xy, xz, yy, yz, zz) order, not three per-component variances. Per-point
  occupancy is tracked separately from the moments because the fluid mask can move.
  Storage is released through the same teardown that releases the block's other
  vectors, and the runloop driver applies accepted weights to every block. Covered
  by a new `unit-statistics-accumulator` suite, including an integration case that
  drives the runloop entry point and asserts only scheduled states reach the fields.

- Added window lifecycle, scheduling, and weighting for the field-statistics
  pipeline: `include/statistics_window.h` and `src/statistics_window.c` decide
  whether a completed state is accepted and what weight it carries, applying
  right-rectangle weighting, final-interval clipping, and the rule that a
  zero-length interval is not a sample. Step and physical-time cadences are both
  supported, with time targets on an absolute grid so the schedule cannot drift,
  and a step overshooting several targets accepted exactly once. Duplicate offers
  of the same completed step are rejected. The runloop calls the window update
  immediately before physical-solution monitoring, after the Lagrangian block and
  before history rotation, and an optional console snapshot mirroring the particle
  console reports window-level progress through the logging sink. Covered by a new
  `unit-statistics-window` suite whose central case asserts that sample and
  physical-time weighting agree on a constant-timestep run.

- Routed the remaining out-of-factory allocations through the designated
  infrastructure. The coarsest-level immersed-boundary `KSKE` workspace was
  allocated and freed on every Poisson solve from inside the solver, with a
  defensive second free in teardown; it is now allocated once by
  `CreateAndInitializeAllVectors` and freed once in teardown, which is safe
  because `FullyBlocked` guards every read with its own flag array and so carries
  no state between solves. The corner-staging workspace was created lazily inside
  the interpolation routine and rebuilt whenever the block size changed; it is now
  two explicitly typed pairs, `CellScalarAtCorner` and `CellVectorAtCorner`,
  created by the same factory.
- Gave the corner-staging workspace typed field identities. Because the vectors
  are now fixed named members with global and local pairs, they satisfy the
  catalog's `offsetof` requirement, so the hand-rolled global-to-local scatter in
  the interpolation path is replaced by `UpdateLocalGhosts`. Corner anatomy
  logging takes the field identity from its caller rather than inferring the
  degree of freedom from a cached vector's block size.

- Added pointwise spatial target resolution for the field-statistics pipeline:
  `include/statistics_target.h` and `src/statistics_target.c` resolve, for one
  field on one block, an iteration domain that excludes PETSc halo storage and
  solver-layout boundary, dummy, and duplicate-periodic indices, which are
  distinct categories. Layout classification comes from the typed field catalog:
  cell-like dimensions span the shifted interior, node-like dimensions carry one
  more entry, and each face family is node-like only in its own direction.
  Periodic directions drop the wrapped duplicate plane, which leaves cell-like
  spans unchanged because their duplicates already sat outside. Component-
  staggered fields are rejected, since their components live on different face
  families and cannot share one pointwise domain. `SpatialTargetPlan` exists with
  only the pointwise identity mapping so later spatial bins extend rather than
  retrofit it. Covered by a new `unit-statistics-target` suite across cell, node,
  and I/J/K-face layouts in nonperiodic, mixed, and fully periodic domains, plus
  a multi-rank case asserting the resolved domain is decomposition independent.

- Added the weighted centered-moment kernels backing the field-statistics
  accumulator contract: `include/statistics_moments.h` and
  `src/statistics_moments.c` implement the stable weighted Welford update, the
  compatible co-moment update, and the weighted parallel-merge formula, plus
  weighted variance, covariance, and Kish effective sample size. The module is
  pure: no configuration, no PETSc vectors, and no window or scheduling
  knowledge. Non-positive sample weights are rejected rather than silently
  corrupting an accumulator. Covered by a new `unit-statistics` suite asserting
  bitwise-zero variance for constant signals, known scalar and two-sample
  results under equal and unequal weights, all six symmetric velocity
  components of a three-sample self-product, exact agreement between a
  self-paired co-moment and the scalar second moment, high-mean/low-fluctuation
  precision where a raw sum-of-squares would cancel, and merge-equals-sequential
  including empty-partition no-ops.

- Added the field-statistics Phase 2 implementation specification, settling the
  contracts page 58 deferred: right-rectangle endpoint quadrature for
  physical-time weighting (chosen because the difference from trapezoidal is an
  endpoint artifact decaying as 1/N, far below the statistical sampling floor,
  while trapezoidal would force a checkpointed field snapshot per accumulated
  field); final-interval clipping at a bounded window's end; and an interval
  convention under which a zero-length interval is not a sample, making
  initial-state handling identical under both weightings and removing the
  `include_initial` control. Also fixes the user contract, cross-field
  covariance spelling, window identity hash inputs, the indexed control-option
  scheme with its ingress-audit extension, the `statistics/` checkpoint
  namespace, the postprocessing contract, and the seven implementation stages.
  Amended page 58 accordingly.

- Field-statistics Phase 1 closeout:
  - added the same-step checkpoint acceptance test, covering both the idempotent
    no-op rewrite of an already-committed step and the rejection of a same-step
    write whose physical time disagrees with the committed bundle;
  - added a multi-rank smoke scenario that restarts one committed bundle at two
    different rank counts, asserting that load-and-store payloads survive the
    decomposition change byte identically and that the resumed run keeps
    marching;
  - documented that particle payloads are deliberately not block scoped in
    version 1 of the checkpoint format, since PICurv carries a single swarm; and
  - corrected stale status prose that still described the checkpoint contract as
    outstanding Phase 1 work; and
  - characterized restart fidelity and recorded it in the maintenance backlog:
    re-applying boundary conditions to loaded fields means a restart resumes from
    a boundary-consistent state that differs from the stored one at the outlet.
    The offset is outlet-local, decays, sits below solver stopping tolerances, and
    is invariant to them, so it cannot be tightened away. Softened the checkpoint
    contract's BDF2 wording accordingly: preserving `Ucont_rm1` preserves the order
    of the restart, not bitwise identity.

- Field identity infrastructure (phase 1):
  - added a 38-entry typed `FieldId`/`FieldDescriptor` catalog with canonical
    names, active aliases, DM family, degree of freedom, shifted/staggered
    layout, synchronization class, availability, capabilities, and existing
    `UserCtx` vector bindings;
  - added non-owning runtime `FieldView` resolution without changing PETSc
    vector allocation or teardown;
  - replaced string dispatch in ghost updates, periodic field synchronization,
    initial interior-field routing, and field-layout diagnostics with catalog
    metadata while preserving their numerical loops;
  - added a separate ten-entry `ParticleFieldId` catalog for DMSwarm component
    count, PETSc data type, registration ownership, initialization, model-update,
    and Eulerian-scatter metadata;
  - converted particle registration, initialization, scalar physics, analytical
    assignment, and Eulerian/particle interpolation/scatter orchestration to
    typed identities;
  - centralized fixed DMSwarm access names and made generic particle gather and
    restart/output casting use PETSc's registered field datatype;
  - retained runtime configuration, function logging/profiling, generic swarm
    I/O, postprocessing recipe dispatch, and all file formats as textual
    boundaries; and
  - expanded setup, periodic, runtime, solver, and logging regression coverage.

- Solver continuation now rejects `run_control.start_step: 0` instead of treating
  a fresh start as an in-place continuation and appending misleading step-zero
  separators to existing logs.

- Picard-Jameson spectral-radius pseudo-time step (semantic change):
  - changed the pseudo-time step from `dtau = pseudo_cfl × dt` to
    `dtau = pseudo_cfl / lambda_max`, where `lambda_max` is the global maximum
    convective spectral radius computed once per physical timestep. This makes
    `pseudo_cfl.*` true dimensionless Courant numbers, independent of `dt`, grid
    size, and flow speed.
  - changed shipped defaults to `pseudo_cfl.initial: 0.5` and
    `pseudo_cfl.maximum: 2.0` (stable 4-stage Jameson range ~0–2.83); existing
    configs with `0.1`/`1.0` still run but now produce physically smaller, safer
    `dtau` values.
  - added a BDF2 lower bound `lambda_max >= COEF_TIME_ACCURACY/dt` so zero-flow
    startup yields a finite `dtau`.
  - added the diagnostic field `mom_last_lambda_max` to the simulation context.
  - changed the per-trial momentum convergence log fields from `Pseudo-cfl`/
    `cfl_after` to `dtau`/`cfl_eff` and `dtau_after`/`cfl_eff_after`; `picurv
    summarize` parsing, JSON payload, console output, and plot series were
    updated to match.

- Initial-condition contract and workflow:
  - replaced user-facing legacy IC selectors with canonical `generated` and `file` modes.
  - added built-in zero, Cartesian constant, streamwise constant, and Poiseuille generators.
  - added single-block file-backed `Ucat` and `Ucont` startup through the existing field reader.
  - added `ic_gen` run/precompute orchestration, defaulting to `generators/ic.gen` with an optional compatible script override.
  - added the repository `generators/ic.gen` expression engine for grid-aware `Ucat` and staggered `Ucont` PETSc vectors.
  - made initial conditions subordinate to analytical, load, and restart Eulerian-state authority.
  - aligned grid, inlet-profile, and IC generators around repository defaults plus optional case-relative or absolute overrides.

- Picard-Jameson momentum-controller hardening:
  - changed pseudo-time trial acceptance and rollback to one global
    transactional decision across blocks and MPI ranks.
  - bounded runtime by counting rejected and accepted trials against
    `max_iterations`, and guaranteed nonfatal exits retain the last accepted
    finite state.
  - unified pseudo-CFL adaptation and carry-over, removed the extra hard-coded
    reduction multiplier and independent low-CFL rebound, and added controller
    bounds validation.
  - added optional `residual_absolute_tol` and `residual_relative_tol`;
    configurations without enabled residual tolerances retain legacy
    update-only convergence semantics.
  - deprecated unused `step_tol` while retaining compatibility ingestion.

- Search robustness observability:
  - added always-on aggregate search instrumentation for particle-enabled runs.
  - added `logs/search_metrics.csv` with timestep-level search, traversal, tie-break, boundary-clamp, bbox-guess, pass-depth, per-step loss, run-local cumulative loss, V2 population/outcome counters, and derived `search_failure_fraction`, `search_work_index`, and `re_search_fraction` signals.
  - added `LOG_SEARCH_METRICS` for compact DEBUG-gated console summaries when explicitly allow-listed.
  - added the `examples/search_robustness/` example family with Brownian Cartesian/curvilinear baselines plus deterministic Cartesian/curvilinear `UNIFORM_FLOW` migration-stress variants, a study starter, and a dedicated metrics-reference docs page.

- Profiling/logging cleanup:
  - removed the `LOG_PROFILE` log level from the C logging enum and all code paths.
  - per-step profiling summaries no longer use a special log level.
  - profiling is now configured explicitly via `profiling.timestep_output` (`off`, `selected`, `all`) and writes timestep rows to a dedicated profiling log file.
  - `profiling.final_summary.enabled` now controls whether the end-of-run `ProfilingSummary_*.log` file is written.
  - `LOG_LEVEL=PROFILE` is no longer a supported runtime setting.
  - removed the temporary `profiling.critical_functions` compatibility shorthand; monitors must now use `profiling.timestep_output`.

- Executable naming refresh:
  - renamed the conductor from `pic.flow` to `picurv`.
  - renamed the main solver executable from `picsolver` to `simulator`.
  - `init` now creates config-only case directories; binaries are resolved from `bin/` via PATH.
  - `init --pin-binaries` copies `simulator`/`postprocessor` into the case for version-pinning (protects running jobs from concurrent rebuilds).
  - `bin/picurv` is now a launcher for the `picurv_cli/picurv` source-tree entrypoint.
  - `sync-binaries` pins specific binary versions into a case directory (optional, equivalent to `--pin-binaries` after init).
  - removed the temporary `--copy-binaries` init flag.

- Grid contract correction:
  - `programmatic_c.im/jm/km` are now treated as cell counts as documented.
  - `picurv` now converts those values to node counts before emitting `-im/-jm/-km`.
  - `grid_gen` remains unchanged: `grid.gen` still accepts cell counts and writes node counts into `.picgrid`.
  - compatibility note: a historical `programmatic_c` case with `im=32` previously yielded 31 physical cells; it now yields the documented 32 physical cells.

- Interface contract cleanup:
  - removed deprecated numeric and alternate-string selector forms for field initialization, particle initialization, and analytical type selection; docs now use canonical YAML values only.
  - renamed the dual-time smoother from `Dual Time Picard RK4` to the more precise `Dual Time Picard Jameson RK`; deprecated RK4 selector, YAML-block, noise-control, C API, and C CLI spellings remain accepted as compatibility aliases.
  - removed the exposed placeholder Newton-Krylov momentum solver modes from parser and C runtime enums.
  - canonical `PICGRID` headers are required for file-based grids in C runtime ingestion.
  - added `grid.gen legacy1d` converter and optional `grid.legacy_conversion` wrapper in `picurv` for headerless 1D-axis legacy payload migration.

- Workflow launch policy update:
  - `run --num-procs` now applies to solver and field postprocessor stage sizing.
  - generated `post.sbatch` and sweep post arrays reuse the configured cluster resources.
  - manifests/dry-run plans now expose stage-specific MPI counts.

- Cluster orchestration and sweeps:
  - added `cluster.yml` Slurm contract support to `picurv run` (`--cluster`, `--scheduler`, `--no-submit`).
  - added `picurv submit` as the delayed-submit counterpart to `--no-submit` for existing run/study artifacts.
  - added `picurv cancel` so Slurm jobs can be stopped by `--run-dir` instead of manual job-id lookup; `picurv cancel --stage solve --graceful` requests solver final-output shutdown with `SIGUSR1`.
  - added scheduler artifact generation and submission metadata (`solver.sbatch`, `post.sbatch`, `submission.json`, `manifest.json`).
  - added `picurv sweep` for parameter studies using Slurm job arrays with post-stage dependency chaining.
  - added `picurv sweep --continue --study-dir <path>` for resuming partially-completed studies: detects per-case completion status, prepares checkpoint restarts via `resolve_restart_source`, and submits sparse solver arrays for incomplete cases only.
  - added `picurv sweep --reaggregate --study-dir <path>` for manual metrics re-aggregation on existing study outputs.
  - added automatic post-completion metrics aggregation via a chained Slurm job (`metrics_aggregate.sbatch`, `afterany` dependency on post array).
  - `detect_last_checkpoint_step` now falls back to particle checkpoint files (`position*.dat`) for analytical-mode cases with no eulerian output.
  - added study aggregation outputs (`metrics_table.csv`, `results/plots`, `summary.json`).
  - added templates: `master_cluster.yml`, `master_study.yml`.
  - added docs pages: cluster run guide and sweep/study guide.
  - documented signal-triggered final snapshot writes for impending walltime/termination handling (`SIGUSR1`, `SIGTERM`, `SIGINT`) with launcher-specific Slurm signal guidance.
  - added automatic runtime walltime guarding for generated Slurm solver jobs, with tunable `execution.walltime_guard` policy and exported batch metadata (`PICURV_JOB_START_EPOCH`, `PICURV_WALLTIME_LIMIT_SECONDS`).

- Documentation enforcement and CI:
  - enforced function-level Doxygen-compatible coverage across C product code, C tests, Python product scripts, and Python tests.
  - added `tests/tooling/audit_function_docs.py` as the repository-wide audit gate.
  - wired the audit into repository consistency tests and the GitHub Actions quality workflow as an explicit pre-pytest step.

- Run inspection and local log routing:
  - added `picurv summarize` for read-only per-step health summaries derived from existing run artifacts.
  - expanded `picurv summarize` with additive `--overview`, `--case`, `--solver`, and `--monitor` configuration views that work before timestep artifacts exist.
  - added `summarize --list-plot-series` and `--plot` time-history workflows backed by standalone `generators/plot.gen`, with append-order windows, automatic residual/norm log scaling, explicit saves, and headless fallback.
  - report missing `matplotlib` as an actionable `DEPENDENCY_MISSING` plotting error instead of a configuration-value error.
  - install `matplotlib` by default in the bootstrap-managed Python environment so plotting works after a standard install.
  - preserve the existing pip version during normal bootstrap runs; add explicit `--upgrade-pip` for environments that need it.
  - verify Python dependencies under the launcher's isolated environment and report the specific missing plotting dependency.
  - local wrapper stream logs for solver/post stages now write under `runs/<run_id>/scheduler/` instead of `runs/<run_id>/logs/`.
  - summary lookup can use continuity, particle metrics, momentum, Poisson, profiling, and sampled particle snapshot artifacts when present.
  - sampled particle snapshot summaries now include compact speed/bounds/rank/weight diagnostics plus sampled deltas against the previous snapshot when possible.

- Initial Doxygen docs scaffold, architecture, and developer guide pages.
- Main docs refresh:
  - updated tutorials/references for current YAML -> `picurv` -> C contract.
  - added method overview pages (CurvIB, fractional-step, dual-time Jameson RK, pressure Poisson/multigrid, walking search, interpolation/projection, IEM/averaging).
  - added non-dimensionalization page and linked it into nav/reference flow.
  - updated landing page visuals (`docs/assets/curv.gif`, `docs/assets/paraview_flat_channel.png`).
  - added maintenance backlog page for low-priority warning/refactor tracking.
- Repository organization updates:
  - moved build configs to `config/build/`.
  - moved shared `grid.gen` profile to `config/grids/coarse_square_tube_curved.cfg`.
  - converted `sandbox/` into explicit developer-sandbox documentation.
  - standardized internal directory guides to `guide.md` (non-root).
  - moved documentation media into `docs/assets/`.
  - removed non-root `README.md` files (single top-level `README.md` retained).
  - routed Doxygen warning output to `logs/doxygen.warnings`.
- Hardened YAML -> C config contract across case/solver/monitor/post:
  - scalar-only `da_processors_*` contract with validation.
  - structured mappings for analytical type and monitor subdirectories.
  - post input extension keys and statistics pipeline wiring.
- Added config contract documentation set:
  - `14_Config_Contract.md`
  - `15_Config_Ingestion_Map.md`
  - `16_Config_Extension_Playbook.md`
- Added ingress drift guard:
  - `tests/tooling/audit_ingress.py`
  - `tests/tooling/audit_ingress_manifest.json`
- Added documentation note on data-driven particle-closure expansion:
  - offline workflows supported now via solver/post artifacts.
  - tightly coupled inference requires runtime-selectable closure interface extension.
- Removed obsolete `stubs/` archive from repository after extracting useful documentation content.
