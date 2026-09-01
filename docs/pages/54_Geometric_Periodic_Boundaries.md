/**
@page p54_geometric_periodic Periodic Boundaries and Driven Flows

@section p54_scope_sec 1. Scope

This page covers two related things: the geometric contract that makes a
periodic axis legal at all (sections 2-4), and the driven-flow handlers that
sustain a bulk flow through such an axis (section 5). The page name kept its
original `geometric` anchor for reference stability; its subject is now periodic
boundaries generally.

PICurv supports geometric periodicity for Eulerian fields. Periodic axes are
derived exclusively from paired `PERIODIC` boundary conditions before the
DMDA is created. The resulting internal periodic flags select PETSc periodic
index wrapping; they are not user configuration options.

Particle positions are not currently wrapped across periodic boundaries.

@section p54_grid_sec 2. Grid Contract

For every active periodic axis, the two opposite nodal surfaces must match
pointwise under one nonzero constant Cartesian translation:

\f[
\mathbf{x}_{+}(p) - \mathbf{x}_{-}(p) = \mathbf{T}
\f]

for every corresponding surface point \f$p\f$. The translation may have
components in all Cartesian directions, but it may not vary over the surface.
Each periodic axis must contain at least four physical nodes. These checks are
performed independently, so mixed periodic/non-periodic domains are supported.

This admits translated curvilinear periods, including a repeating wavy
channel, while rejecting rotating, shearing, or otherwise nonmatching paired
surfaces.

@section p54_fields_sec 3. Field Synchronization

PETSc `DM_BOUNDARY_PERIODIC` wraps indices but does not translate coordinates
or resolve the ownership conventions of every staggered field. PICurv therefore
handles periodic Eulerian data by layout:

- cell-centered fields,
- persistent fields belonging to one face family,
- component-staggered fields such as `Ucont`,
- local-only staggered work fields.

The immutable field catalog records the conceptual layout and ghost-repair
class for each persistent field. The actual periodic axes still come from the
runtime DMDA and validated boundary pairs; field identity does not replace or
override that state. See @ref 56_Field_Identity_and_Layout_Catalog.

Face-center and cell-center coordinate ghosts are wrapped and shifted by the
validated translation. Higher-depth QUICK stencil preparation is handled
separately from ordinary local ghost refresh.

@section p54_config_sec 4. Configuration

Configure both faces of an axis as periodic:

```yaml
boundary_conditions:
  - face: "-Xi"
    type: PERIODIC
    handler: geometric
  - face: "+Xi"
    type: PERIODIC
    handler: geometric
```

Do not add `i_periodic`, `j_periodic`, or `k_periodic` under `models.domain`;
those legacy keys are rejected.


@section p54_driven_sec 5. Driven Periodic Flows

A streamwise-periodic domain has no inlet and no outlet, so nothing sustains a
bulk flow against wall drag: left alone, the flow simply decays. A driven
handler closes that loop. It measures the volumetric flux, compares it with a
target, and feeds the difference back as a body force in the momentum equation.

Two handlers implement this. They differ in exactly one respect - where the
target comes from - and are otherwise the same code.

| Handler | Target | `target_flux` parameter |
|---|---|---|
| `constant_flux` | Prescribed by the user | **Required** |
| `initial_flux`  | Measured from the field the run starts with, then held | **Rejected** - supplying one is a configuration error |

`initial_flux` is the right choice when you have seeded a field whose flux is
already what you want and you simply want it preserved - a restart from a
precursor, or a developed field carried over from another case. `constant_flux`
is the right choice when the flux is the thing you are prescribing, for instance
when matching a literature bulk Reynolds number.

@subsection p54_driven_control_sub 5.1 What the controller measures

`PreStep` takes two different flux measurements in one sweep of `lUcont`
(`MeasureDrivenFluxes` in `src/BC_Handlers.c`), because the controller needs a
fast signal and a stable one and they cannot be the same number.

- **Boundary flux** - the flux through the single periodic boundary plane. It
  responds immediately to anything happening at the seam where the domain wraps,
  and it is correspondingly noisy. It drives `boundaryVelocityCorrection`, a
  tactical trim applied directly to the boundary fluxes.
- **Plane-averaged flux** - the flux averaged over every cross-sectional plane
  in the driven direction. Local turbulent fluctuations average out, so it
  represents the bulk momentum of the fluid and changes only slowly. It drives
  `bulkVelocityCorrection`, which scales the momentum source. It is also the
  quantity `initial_flux` latches at startup.

Driving the body force from the noisy signal would make the forcing oscillate;
trimming the seam from the slow signal would leave the seam inconsistent. Hence
two sensors.

@subsection p54_driven_cadence_sub 5.2 Update cadence, and why it matters

These two corrections are refreshed on deliberately different cadences.

`ApplyBoundaryConditions()` runs the handler sweep three times per call, and it
is itself called once per Jameson RK stage under the Picard solver and once per
residual evaluation under Newton-Krylov. "Once per `PreStep`" is therefore very
far from "once per timestep".

- **`bulkVelocityCorrection` is computed once per physical timestep** and held
  for the rest of it. Both momentum solvers depend on this. The Picard
  shadow-Jacobian estimate treats the body force as a constant forcing with zero
  velocity Jacobian, and a matrix-free Newton solve needs a source that does not
  drift between residual evaluations. A source that moved with the trial vector
  would silently change the operator both solvers think they are solving.

  Freezing the correction is necessary but **not sufficient on its own**. The
  force actually applied is produced downstream by
  `ComputeDrivenChannelFlowSource()` (`src/BodyForces.c`), which smooths it
  against the previous step with a 50/50 exponential moving average held in
  `SimCtx`. That function runs from `ComputeRHS`, so it too executes once per RK
  stage or residual evaluation. Advancing the EMA on every call walked the
  applied force toward its target *within* a single timestep - 0.5, then 0.75,
  then 0.875 of the way - making the force depend on how many residual
  evaluations preceded it. Both the correction and the EMA are therefore gated
  on `simCtx->step`, and only the pair of them makes the source genuinely
  constant across a timestep.
- **`boundaryVelocityCorrection` is re-measured on every pass.** `Apply()`
  accumulates the trim into `Ucont`, and re-measuring is exactly what makes that
  accumulation self-limiting as the seam converges. Freezing it would make
  repeated `Apply()` calls add the same trim over and over.

Because the momentum source is frozen at the start of a step, the controller
holds the flux to first order in `dt`: the force applied during a step is the
one that would have corrected the flux measured at its beginning. Halving `dt`
halves the lag.

@subsection p54_driven_trim_sub 5.3 The seam trim (`enforce_seam_flux`)

The controller has two actuators, and this option controls the second one.

The **body force** is the gentle, global actuator: it pushes on every fluid cell
and is what actually sustains the flow. The **seam trim** is a local override.
It adds a correction straight into the face fluxes of the single boundary plane
where the domain wraps around:

```c
fluxTrim = boundaryVelocityCorrection * faceArea;
if (enforceSeamFlux) ucont[k_face][j][i].z += fluxTrim;    // edits the solution
uch[k_face][j][i].z = fluxTrim;                             // always recorded
```

"Trim" is aviation vocabulary - a small standing correction that holds a
quantity on target without touching the main control. Here the quantity is the
volumetric flux through the seam plane, and the main control is the body force.
The option was originally spelled `apply_trim`, which said that something was
trimmed but not what or why; it is now `enforce_seam_flux`, and the old spelling
is still accepted as a deprecated alias in both `case.yml` and `bcs.run`.

Note the asymmetry: `ucont` **accumulates** (`+=`) while `uch` **overwrites**
(`=`). The trim is re-measured on every boundary pass, and that re-measurement is
what makes the accumulation self-limiting - each pass sees a seam that is already
closer to target, so the next correction is smaller. This is why the trim
correction is deliberately *not* frozen per timestep the way the body force is
(section 5.2); freezing it would make repeated passes add the same correction
over and over.

**It is off by default, and that is the right default.** Reasons to leave it off:

- *It can fight the projection.* The pressure solve has just enforced a
  divergence-free field. The trim modifies face fluxes afterwards, so it can
  reintroduce a small divergence the projection removed. Where continuity is the
  accuracy limit, an imperfect seam flux is the cheaper error.
- *It is driven by the noisy sensor.* The trim uses the single-plane flux, which
  jitters step to step in turbulent flow. Applying it injects that jitter into
  the solution at one plane - a persistent, spatially localised perturbation
  exactly where periodicity is supposed to be seamless. For DNS or LES
  statistics gathered near the seam, that is contamination.
- *The body force is usually enough.* In a converged periodic flow the seam flux
  and the bulk flux agree, so the trim is near zero anyway. It earns its keep
  during startup transients, or when the seam is genuinely drifting.
- *It re-enters the residual path.* Because it is re-measured every pass, it is
  solution-dependent in a way the frozen body force is not. Leaving it off keeps
  the Newton-Krylov residual a cleaner function of the trial vector.

Reasons to turn it on: you have diagnosed an actual flux mismatch at the seam -
inspect `Bcs.Uch`, which records the trim whether or not it is applied - and the
body force alone is not closing it.

Practical workflow: run with it off, look at `Uch`. If the recorded trim is
negligible compared with the bulk flux, the seam is fine and enabling it would
only add noise. If it is persistently large, that is a real seam problem, and
enabling the trim treats the symptom while you look for the cause.

@subsection p54_driven_restart_sub 5.4 Restart behaviour

The two handlers restart differently, and the difference is deliberate.

`constant_flux` re-reads `target_flux` from the bcs file on every run, so
editing it between segments takes effect. Nothing is carried in the checkpoint.

`initial_flux` cannot do that: its target is a property of the field the run
*originally* started from. Re-measuring it after a restart would silently
retarget the controller at whatever the flux had drifted to. The latched value
is therefore written into the checkpoint manifest as

```
-checkpoint_driven_flux_latched true
-checkpoint_driven_flux_target  <value>
```

and restored by `RestoreDrivenFluxTarget()` in `src/io.c` before the first
`PreStep` of the resumed run. A checkpoint written before this metadata existed
has no entry; the controller then falls back to re-measuring and logs a warning
saying so.

One consequence worth knowing: the latch does **not** fire during setup.
Boundary handlers are initialized before `InitializeEulerianState()`, and the
setup-time boundary sweep runs before the initial condition has been scattered
into the ghosted `lUcont` the measurement reads. Latching there would record
zero. The handler therefore waits for the first `PreStep` of the first real
timestep, and does no work at all until then.

@subsection p54_driven_config_sub 5.5 Configuration

Both faces of the driven axis must carry the same handler, and the runtime
rejects any `INLET`, `OUTLET`, or `FARFIELD` face while a driven handler is
active - the flux controller and an inlet would be fighting over the same
degree of freedom.

```yaml
boundary_conditions:
  # Spanwise: plain periodic.
  - face: "-Xi"
    type: PERIODIC
    handler: geometric
  - face: "+Xi"
    type: PERIODIC
    handler: geometric
  # Wall-normal: no-slip.
  - face: "-Eta"
    type: WALL
    handler: noslip
  - face: "+Eta"
    type: WALL
    handler: noslip
  # Streamwise: periodic and driven.
  - face: "-Zeta"
    type: PERIODIC
    handler: constant_flux
    params:
      target_flux: 12.566370614359172   # U_b * cross-sectional area
  - face: "+Zeta"
    type: PERIODIC
    handler: constant_flux
    params:
      target_flux: 12.566370614359172
```

`target_flux` is a **volumetric** flux, `U_b` times the cross-sectional area,
not a velocity. For `initial_flux`, drop the `params` block entirely.

@subsection p54_driven_examples_sub 5.6 Worked examples

- `examples/periodic_test/driven_channel/laminar/` - plane channel with an exact
  closed-form answer. For half-height `h`, viscosity `nu` and the converged body
  force `f`:

      u(y)  = (f / 2nu) (h^2 - (y-h)^2)
      U_b   = f h^2 / 3nu      ->   f = 3 nu U_b / h^2
      u_tau = sqrt(f h)

  This verifies the source term and the flux controller with no turbulence model
  involved.
- `examples/periodic_test/driven_channel/turbulent_retau180/` and
  `.../turbulent_retau395/` - DNS-resolution plane channels against Moser, Kim &
  Mansour (1999), Phys. Fluids 11, 943.
- `examples/periodic_test/driven_channel/les_retau180/` - the coarse LES repeat.
  @note Both LES models are **experimental**: implemented and unit-tested, but not yet validated against a reference flow. The channel is periodic in xi and zeta, so `averaging.mode: homogeneous` derives those two directions from the boundary pairs and produces a wall-normal coefficient profile. See @ref p07_les_sec and @ref 72_LES_Turbulence_Closure.
- `examples/periodic_test/driven_duct/` - square duct at `Re_b = 4410`, which
  sustains turbulence-driven secondary flow of the second kind that a plane
  channel does not. References: Gavrilakis (1992), JFM 244, 101; Huser &
  Biringen (1993), JFM 257, 65.

Grids for all of these are generated from checked-in `.cfg` files under
`config/grids/`, not shipped as binary `.picgrid`. Wall spacing is a single
number per walled axis, expressed as a fraction of that axis's length, so `y+`
can be retargeted for a different `Re_tau` without regenerating anything by
hand.

@subsection p54_driven_limits_sub 5.7 Known limitations

**Initial-condition seeding.** There is no generator that produces a
divergence-free perturbed field respecting no-slip, so a wall-bounded turbulent
case has no good seed in-tree. `ic.gen`'s `spectral_random_velocity` is not
usable here: it explicitly requires `PERIODIC`/`geometric` on all six faces and
is rejected otherwise, which is correct - it cannot respect a no-slip wall. The
shipped turbulent cases seed with `streamwise_constant`, a laminar profile;
transition from it is slow and is not guaranteed without a finite-amplitude
perturbation. Until a seeded generator exists, develop a field at coarse
resolution (or from a precursor run) and carry it in, then use `initial_flux` to
hold its flux.

**Wall-normal profiles from statistics.** `field_statistics` is BC-agnostic - it
accumulates per-cell time moments and works unchanged under periodic boundaries.
What does not exist is a spatial reduction: nothing in the postprocessor averages
a statistics field over homogeneous directions to produce the `U+` vs `y+` and
RMS profiles a DNS comparison needs. That reduction currently has to be done
outside the postprocessor from the window payloads.

**Momentum convergence on periodic wall-bounded flow.**

@warning Status: previously observed; requires re-characterization at current
`HEAD`. The observation below was measured on 2026-08-24. Commits landed on
2026-08-25 moved momentum convergence onto a residual criterion, retired
`absolute_tol`, and changed the pseudo-CFL controller's wall-relearning behavior -
that is, they rewrote the exact criterion this was measured against. Do not treat
it as current solver behavior, and do not cite it, until the cases below are rerun.

As measured on 2026-08-24: under the Dual Time Picard Jameson RK solver these cases
did not reach the pseudo-time convergence tolerance. The pseudo-CFL controller drove
`dtau` to its floor and the residual ratio plateaued near 1, so every step reported
"reached N total attempts without convergence" and continued from the last accepted
state. The behavior was independent of the driven-flow machinery - an otherwise
identical case with plain `geometric` handlers behaved the same, while the shipped
wall-bounded inlet/outlet example converged cleanly - and independent of
`central_diff` and of the pseudo-CFL floor. The Poisson side was healthy throughout
(maximum divergence around 1e-14).

If the behavior reproduces at current `HEAD`, it must be resolved before the
turbulent campaigns above can be run. If it does not, replace this note with the
new measured result rather than deleting it silently.

@section p54_diag_sec 6. Diagnostics and Tests

At startup, the banner reports BC-derived periodic axes and each validated
translation. Use `make unit-periodic` for focused synchronization and geometry
coverage and `make smoke-periodic` for the real runtime periodic workflow.
Both are included in the standard `unit` and `check` gates.

For driven flows, `tests/smoke/run_driven_periodic_regression.sh` asserts that
the `initial_flux` target matches the analytic flux of the initial condition and
survives a restart unchanged, alongside the Poisson coarse-solve residual check.
The driven controller logs its target, both measured fluxes and both corrections
at `INFO` every step; the realized flux should track the target to the
first-order-in-`dt` lag described in section 5.2.

@section p54_refs_sec 7. Related Pages

- @subpage 44_Boundary_Conditions_Guide
- @subpage 07_Case_Reference
- @subpage 40_Testing_and_Quality_Guide
- @subpage 56_Field_Identity_and_Layout_Catalog
- @subpage 55_Newton_Krylov_Momentum_Solver
- @subpage 25_Pressure_Poisson_GMRES_Multigrid
*/
