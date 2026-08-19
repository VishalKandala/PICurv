@page 58_Field_Statistics Field Statistics

@anchor _Field_Statistics

PICurv accumulates weighted centered moments of Eulerian fields while the solver
runs. Turbulence statistics therefore depend on the whole simulated interval
rather than on however many instantaneous states happened to be written, and a
converged average costs one accumulator per field instead of a directory of
snapshots.

This page covers what the system computes, the semantics that change a result,
how state survives a restart, and where to extend it. The YAML keys themselves
are specified in @ref 09_Monitor_Reference and @ref 10_Post_Processing_Reference;
@ref 60_Field_Statistics_Planned_Extensions records what is not built yet.

@tableofcontents

@section p58_scope_sec 1. What It Computes

The solver accumulates only primitive fields and the centered moment state the
requested products need. Everything derived — Reynolds stresses, RMS, turbulent
kinetic energy, turbulent fluxes — is computed in post-processing from that
state. The split is deliberate: derived quantities are cheap to recompute and
expensive to store, and keeping them out of the solver means a recipe can be
changed without re-running.

Available for accumulation: `Ucat`, `P`, `Nvert`, `Phi`, `Psi`, `ParticleCount`,
`Nu_t`, and `CS`. The last four require their subsystem to be active — particles
for `Psi` and `ParticleCount`, a turbulence model for `Nu_t`, LES for `CS` — and
a field whose subsystem is off is rejected at validation rather than accumulated
as zeros.

Three boundaries are worth stating because they are easy to assume away.

**LES and RANS quantities stay model quantities.** `Nu_t` and `CS` can be
averaged like any other field, but modelled `k` is never relabelled as the
turbulent kinetic energy derived from resolved velocity fluctuations. The two
are different quantities and the pipeline does not conflate them.

**Full-field time histories are out of scope.** The accumulator holds moments,
not snapshots. A run that needs the states themselves selects ordinary
instantaneous field output.

**Particle ensemble statistics are a separate subsystem.** Mean-square
displacement and its relatives live under `post.yml -> statistics_pipeline` and
share no state with field windows. See @ref 28_IEM_and_Statistical_Averaging.

@section p58_window_sec 2. Windows

A window is an independently accumulating average over a named stretch of
physical time. Windows may overlap freely; each owns its own state, so two
windows covering the same interval with different cadences do not interact.

A window is `pending` before its start time, `active` while accumulating, and
`complete` once a bounded window passes its end. It records its requested and
effective bounds, sample count, total weight, represented physical time,
schedule position, last accepted event and time, resolved definition hash, and
restart lineage.

@subsection p58_cadence_sub 2.1 Cadence

Exactly one cadence is selected per window.

- `step_cadence: n` accepts a completed state every `n` steps, counted from the
  window's effective start.
- `time_cadence: dt_s` places targets on the absolute grid
  \f$t_{\text{start}} + k\,\mathit{dt_s}\f$ and accepts the **first completed
  state at or past** each target.

Targets sit on the absolute grid rather than relative to the last accepted state,
so the schedule cannot drift as `dt` varies. When one large step overshoots
several targets, that state is accepted exactly once and the schedule advances
past the current time. No time is lost or double counted, because the accepted
state's weight is the interval actually elapsed rather than the nominal cadence —
a direct consequence of the quadrature rule in @ref p58_weighting_sec.

@subsection p58_bounds_sub 2.2 Bounds, Clipping, and the First State

A window represents the half-open intervals \f$(t_{i-1}, t_i]\f$: each accepted
state carries the interval ending at it, measured from the previous accepted
state, or from the window's effective start for the first. **A state
representing a zero-length interval is not a sample.**

That single convention settles what would otherwise be a configuration key. A
state landing exactly on `start_time` anchors the interval origin without being
sampled; the next accepted state is a sample carrying the whole interval back to
the start. No sample is wasted, the represented interval is exactly
\f$[t_{\text{start}}, t_{\text{end}}]\f$, and the two weighting modes agree
exactly on a constant-timestep run:

\f[
\langle f\rangle_{\text{time}}
= \frac{\sum_{i=1}^{N} f_i \Delta t}{N\Delta t}
= \frac{f_1 + \cdots + f_N}{N}
= \langle f\rangle_{\text{sample}} .
\f]

Counting the anchoring state as a sample under `sample` weighting only would make
the two modes disagree by \f$O(1/N)\f$ on a fixed-timestep run for no physical
reason.

A bounded window's final state has its weight clipped to
\f$t_{\text{end}} - t_{\text{last}}\f$. Start and end clip by the identical rule;
neither boundary is a special case.

Two edge behaviours follow from the same convention and are stated so they are
not mistaken for defects:

- `effective_start` is the requested `start_time` unless the window is first
  observed later than that, in which case it is the first observed state time and
  the window records that it represents less than requested. A window never
  claims time it did not observe.
- With a coarse `step_cadence` the first sample legitimately represents up to `n`
  steps back to the effective start.

A window that reaches `complete` having accepted no sample is a hard error, not
an empty result for post-processing to divide by.

@section p58_weighting_sec 3. Weighting and Quadrature

`sample` weighting gives every accepted state equal weight. `physical_time`
weighting gives a state the weight of the interval it represents, which is what a
variable timestep requires for the average to mean anything.

Physical-time weighting uses the **right-rectangle rule**: an accepted state at
\f$t_i\f$ carries \f$w_i = t_i - t_{i-1}\f$.

The accuracy difference from the trapezoidal rule is an endpoint artefact rather
than an accumulating bias. For uniform \f$\Delta t\f$ the two rules differ by

\f[
\sum_{i=1}^{N} f_i \Delta t \;-\;
\Big[\tfrac{f_0\Delta t}{2} + \sum_{i=1}^{N-1} f_i \Delta t + \tfrac{f_N \Delta t}{2}\Big]
= \frac{(f_N - f_0)\,\Delta t}{2},
\f]

one endpoint term independent of \f$N\f$. Divided by the window length, the
difference in the mean decays as \f$1/N\f$ regardless of timestep size, which for
a converged turbulence average sits orders of magnitude below the sampling error
\f$\sigma/\sqrt{M}\f$.

The decisive argument is state cost. Right-rectangle weights depend only on the
past, so the only quadrature state crossing a restart is the last accepted time —
a scalar. Any trapezoidal or centered rule weights sample \f$i\f$ using
\f$t_{i+1}\f$, and because the weighted update is nonlinear in the weight,
retroactively weighting a past sample requires that sample's value. That would
force a full pending field snapshot per accumulated field to be held,
checkpointed, and restored exactly: roughly double the state and a new restart
failure mode, bought for a correction far below the noise floor.

@section p58_state_sec 4. Accumulated State and Component Order

@subsection p58_update_sub 4.1 Per-Point State and Update

Each point stores sample count, total weight \f$W\f$, squared-weight sum
\f$W_2\f$, mean \f$\mu\f$, and centered second-moment sum \f$M_2\f$. For a sample
\f$x\f$ with weight \f$w > 0\f$:

\f[
W' = W+w,\qquad \delta=x-\mu,\qquad
\mu'=\mu+\frac{w}{W'}\delta,\qquad
M_2'=M_2+w\,\delta(x-\mu').
\f]

The compatible co-moment update for a pair \f$(x, y)\f$ is

\f[
C' = C + w\,(x-\mu_x)\,(y-\mu_y'),
\f]

evaluated with the pre-update \f$\mu_x\f$ and the post-update \f$\mu_y'\f$, the
form that stays consistent with the \f$M_2\f$ update when \f$x = y\f$. Because
the update is centered and incremental, it stays stable for a high mean with
small fluctuations, where accumulating raw sums of squares would not.

\f$W_2\f$ supports effective-sample-size and uncertainty estimates and is
required because physical-time weighting produces unequal weights.

A stable weighted parallel-combination form is implemented alongside the
sequential update, so that offline window merging and spatial reductions inherit
a verified kernel:

\f[
\delta = \mu_B - \mu_A,\qquad W = W_A + W_B,\qquad
\mu = \mu_A + \frac{W_B}{W}\delta,\qquad
M_2 = M_{2,A} + M_{2,B} + \frac{W_A W_B}{W}\delta^2 .
\f]

@subsection p58_order_sub 4.2 Products and Component Order

Component order is explicit everywhere and never implicit.

- **Scalar `first`** stores \f$\mu\f$; **scalar `second`** adds \f$M_2\f$, one
  component.
- **Vector `first`** stores \f$\mu\f$ per component.
- **Vector `second`** stores the full symmetric self-product. For a
  three-component field that is **six** components in upper-triangular row-major
  order,

  \f[
  (xx,\; xy,\; xz,\; yy,\; yz,\; zz),
  \f]

  from index pairs \f$(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)\f$. These are the
  centered Reynolds-stress sums: `moments: [second]` on `Ucat` requests all six.
  **Reynolds stresses are not cross-field products and need no covariance entry.**
- **Cross-field covariance** between a three-component vector and a scalar stores
  three components in field component order
  \f$(x\!\cdot\!s,\; y\!\cdot\!s,\; z\!\cdot\!s)\f$ — the turbulent scalar flux.

Covariance requires both members to resolve to the same layout, which currently
means cell-centered with cell-centered. `Ucont` is component-staggered, so pairs
involving it are rejected; explicit face-flux products are a planned extension.
Vector-vector cross products are rejected, because they would need a full
nine-component tensor that nothing yet carries.

@section p58_target_sec 5. Spatial Target, Layout, and Mask

Field identity and layout are resolved through the typed catalog
([field_catalog.h](../../include/field_catalog.h)) and ownership through the
block's `DMDALocalInfo`, so no second field lookup, layout table, or range
convention exists. See @ref 56_Field_Identity_and_Layout_Catalog.

Iteration excludes two distinct categories: PETSc halo and ghost storage, and the
solver layout's boundary, dummy, and duplicate-periodic indices. On a periodic
axis the duplicate plane is dropped, so a triply periodic box counts every cell
exactly once.

The mask is `fluid`, defined by `Nvert < 0.1`. Because `Nvert` changes when an
immersed body moves, the mask is treated as moving: **every point keeps its own
valid count and valid weight**, and a `valid_fraction` diagnostic reports the
point's valid weight against the window's total. A point that is never valid keeps
zero weight and is reported, rather than producing a divide-by-zero downstream or
a zero that reads like a measurement.

@section p58_identity_sec 6. Window Identity and What May Change

At startup the resolved definition is hashed once with SHA-256
([checksum.h](../../include/checksum.h)). The hash is checkpoint metadata, not a
user option. Hashed, in fixed order: window name; `start_time`; weighting mode;
cadence kind and value; resolved field entries ordered by field identity, each
with its moment set; resolved covariance pairs, canonicalized within each pair and
ordered; mask identity; and target kind and layout semantics.

Each group is also hashed on its own, and the truncated group digests are stored
beside the full one. That is what lets a mismatch name the property that changed
rather than reporting only that two digests differ.

Explicitly excluded, so that changing them continues rather than invalidates a
saved window: `end_time` and `bounded`; `enabled`; the console cadence; and the
order in which fields and covariance pairs were listed.

**Extending a window forward is accepted** when the new end moves forward, does
not exceed the run horizon, and continuation resumes without an unsampled gap.
The gap is checked against the checkpoint's own time: reopening a window that
closed at \f$t_{\mathrm{end}}\f$ from a bundle written at \f$t_{\mathrm{ckpt}}\f$
is refused once \f$t_{\mathrm{ckpt}} - t_{\mathrm{end}}\f$ exceeds one timestep.
This follows from right-rectangle weighting: the first sample after the former end
carries the whole interval back to it, so a gap would be weighted as though it had
been observed. One step of slack is allowed because the closing state's own
clipping already leaves up to that much.

**Shortening a window is refused**, because an `end_time` earlier than the time
the saved state already represents would discard samples the metadata still
counts. Changing any hashed property requires a new window name.

@section p58_checkpoint_sec 7. Checkpoints, Restart, and Continuation

Window state rides in the existing committed checkpoint bundle. Payloads go
through the same vector writer Eulerian fields use, so rank-count independence is
inherited rather than reimplemented and no second binary serializer exists.

```text
output/checkpoints/step_000000001000/
  statistics/
    window_0000/
      block_0000/
        count.dat               # per-point accepted sample count
        weight.dat              # per-point valid weight W
        weight_sq.dat           # per-point squared-weight sum W2
        Ucat_mean.dat           # dof 3, the source field's own layout
        Ucat_m2.dat             # dof 6, order (xx,xy,xz,yy,yz,zz)
        Psi_mean.dat            # dof 1
        Psi_m2.dat              # dof 1
        Ucat_Psi_cm.dat         # dof 3, order (x*s, y*s, z*s)
```

Payloads are **block scoped**, so a multiblock run keeps one accumulator tree per
block under each window. Each product is **one payload carrying all its
components** rather than one payload per component: a symmetric tensor is a single
object, splitting it would cost six memory streams in the accumulation loop and
six collective gathers per checkpoint instead of one, and the split does not
generalize — it works for a three-vector only because six happens to equal three
plus three. Payload names come from the catalogued field name plus a fixed role
suffix, so they are stable across runs, rank counts, and configuration
reorderings. A field with no second moment requested contributes no product
payload.

Every statistics payload enters the manifest's payload inventory, qualified by
window as `<window>/<payload>`, and is validated for existence and size before the
bundle is accepted. Per-window scalar metadata is recorded in `checkpoint.meta`,
including the definition hash and its group digests, window state, sample count,
total weight, represented time, effective bounds, and restart lineage.

Schedule state is saved as well as reported state, and this is not cosmetic: the
activation step anchors the step stride, the next time target anchors the absolute
time grid, and the last event step is the duplicate-event guard. A restart that
restored only the reported quantities would resume on a silently shifted schedule
while still reporting a plausible sample count.

Restart rules:

- `--continue` resumes in place and continues a window when its hash matches, with
  a permitted `end_time` extension checked separately.
- `--restart-from` continues statistics only when the physical case identity and
  the definition hash both still match and continuation was explicitly requested.
  Otherwise the flow may restart while statistics reset under a new window name.
- Missing state for a requested continuation is fatal; it is never silently zeroed.
- A hash mismatch is fatal, naming the window and the first differing property.
- A different MPI rank count is permitted, because statistics payloads use the
  same natural ordering Eulerian payloads do.

Restoration happens with the flow state, from the same bundle and independently of
the Eulerian source, so an analytical restart still continues a window that was
accumulating before it.

@section p58_derived_sec 8. Derived Quantities and Post-Processing

A recipe at `post.yml -> field_statistics` names windows; it does not re-describe
them. The post-processor is launched with the run's own control file and resolves
the window list through the same parser the solver uses, so no second loader
exists to disagree with the first.

Derived quantities are normalized in one place from centered state:

\f[
R_{ij} = \frac{C_{ij}}{W},\qquad
\mathrm{RMS}_i = \sqrt{R_{ii}},\qquad
k = \tfrac{1}{2}\left(R_{xx}+R_{yy}+R_{zz}\right),\qquad
\overline{u_i'\psi'} = \frac{C_{i\psi}}{W}.
\f]

Small negative variances from floating-point cancellation are clamped only when
taking a square root, and only within a documented tolerance. Stored state is
never mutated.

The pipeline runs once per processed step, so a recipe spanning several steps
produces a convergence history rather than one picture repeated; pinning
`source_step` collapses it to a single bundle. Each window writes its own file,
which is forced rather than preferred: one window carrying every output already
produces fifteen fields against a per-file limit of twenty, so two windows cannot
share a file.

`vtk` writes the derived fields. `csv` appends one row per processed step with
sample count, total weight, represented time, the per-point valid-fraction range,
and the mean turbulent kinetic energy. The CSV is the artefact that answers
whether a window has run long enough, which no single snapshot can.

That mean is taken over **the fluid cells the window actually sampled**, not over
the stored vector. A derived field is zero outside the target domain and wherever
the mask excluded a point, and those zeros are absences rather than measurements:
averaging over them scales the result down by whatever fraction of the vector the
window never covered. The valid-fraction column in the same row reveals whether an
immersed body made the sampled set smaller than the targeted one; at one they are
the same set.

A window with no sample yet at a given step is skipped with a note rather than
treated as an error, because a recipe covering a whole run legitimately reaches
bundles from before that window began.

**Derived statistics are not dimensionalized**, even when
`global_operations.dimensionalize` is set. A Reynolds stress scales as velocity
squared and a co-moment as the product of two different scales; the per-field
scaling table expresses neither, so applying a velocity scale would be wrong
rather than merely incomplete.

@section p58_monitoring_sec 9. Runtime Monitoring and Logging

`io.statistics_console_output_frequency` emits a periodic console snapshot, built
on the same three-part gate as the particle snapshot it mirrors: the subsystem is
active, the cadence is positive, and the effective log level is at least
`LOG_INFO`. A zero cadence silences reporting without disabling accumulation.

The key sits under `io` rather than inside `field_statistics` because it controls
reporting rather than science. Changing it must not alter any accumulated result,
which is why it is excluded from the definition hash — and is verified by running
one field series at two cadences and requiring bit-identical accumulator state.

Each snapshot reports, per window: name, state, sample count, total weight,
represented time, progress toward a bounded end, and the range of per-point valid
fraction. These are window-level scalars; **the snapshot never dumps field data**.
The valid-fraction range is the mask-health indicator: a minimum of one means every
targeted point saw every state, below one means a moving body excluded some points
for part of the window, and zero means some point contributed nothing and its mean
is undefined. That last case cannot be seen by reading the mean field, because an
untouched point holds a zero that looks like a value, so it is also reported once
as a warning when the window completes. Computing the range is a collective
reduction, so it runs on the console cadence and at completion, never per step.

The startup banner always carries a statistics line — the configured cadence, or
`DISABLED` distinguishing a silenced console from a run with no window — so a log
records whether monitoring was active rather than leaving its absence ambiguous.

Lifecycle logging uses `LOG_ALLOW` at the levels the rest of the codebase uses:
`LOG_INFO` for once-per-run and once-per-window events, `LOG_DEBUG` for per-event
bookkeeping including states an active window did not sample, `LOG_WARNING` only
for conditions an operator must act on, and `LOG_TRACE` for per-point detail.
`LOG_ALLOW` carries the codebase's usual dual gate; see @ref p09_logging_sec for
the mechanism and @ref p11_logging_ssec for the workflow. The console snapshot is
the exception that uses plain `LOG`, because an accumulating window is an
operator-facing report rather than an instrument for debugging the code.

The functions to name in `logging.enabled_functions` when tracing this subsystem:

| Add to `logging.enabled_functions` | Reports |
| --- | --- |
| `LogResolvedWindow` | each window definition as resolved at startup |
| `PicurvWindowStorageCreate` | accumulator storage allocated per window |
| `FieldStatisticsUpdateWindows` | window activation and completion, and each accepted or skipped state |
| `WriteStatisticsFields` | window state written into a checkpoint bundle |
| `ReadStatisticsWindowState`, `RestoreFieldStatisticsState` | window state read back on continuation |
| `FieldStatisticsPipeline` | per-step derivation and output during post-processing |

The filter matches the function containing the call, so `LogResolvedWindow` is the
name to use for the resolved definitions rather than the configuration entry point
that calls it.

Accumulation sits on the timestep path, so debug-level paths perform no collective
reduction unless the level is already active.

@section p58_architecture_sec 10. Architecture and Extension Points

| Concern | Location |
| --- | --- |
| Centered moment and co-moment kernels | `statistics_moments.c` |
| Window lifecycle, scheduling, identity hash | `statistics_window.c` |
| Per-window storage, accumulation, derivation | `statistics_accumulator.c` |
| Spatial target plan, layout and mask resolution | `statistics_target.c` |
| Control resolution from the generated control file | `statistics_config.c` |
| Checkpoint write, validate, restore | `io.c` |
| Post-processing pipeline and derived output | `postprocessor.c`, `postprocessing_kernels.c` |

Two prefixes split the symbols by role rather than by file. `Picurv*` marks the
reusable library — the moment kernels, the window lifecycle, the accumulator —
which take plain arguments or a block context and hold no opinion about how a run
is configured. `FieldStatistics*` marks the run-integration surface: the predicate
the runloop, checkpoint writer, and console monitor share, the per-step driver, and
the post-processing pipeline. Configuration entry points keep the codebase's
verb-first form.

Three structural choices are recorded so a later change knows what was weighed.

**The accumulator module holds two halves**: the online accumulation that runs
every accepted step, and the offline derivation and spatial reduction that run only
in post-processing. They share exactly one thing, the symmetric component-pair
table that defines product order. The moment to split them is when a further
product family arrives and that table needs a home of its own.

**Several moment kernels have no production caller yet.** The reset functions are
their structs' constructors; the merge functions implement the weighted parallel
combination that makes spatial reduction exact after the fact, which is why
profiles, regions, and bins are post-processing operations rather than
accumulator kinds (see @ref p60_principle_sec); the scalar variance accessor is correct but
unused because every product is routed through the co-moment path so the diagonal
and off-diagonal share one update; and the effective-sample-size accessor needs a
window-level squared-weight sum surfaced before it can be reported. Each is covered
by the moment suite, so a future consumer inherits a tested kernel.

**Accumulator vectors are created in the one vector factory** that every other
vector a run owns is created in, and released through the same teardown. A
config-counted array is no exception; the window count is resolved from
configuration long before the factory runs. Component counts that neither the
scalar nor the vector DM provides are carried by a DM mirroring the block
decomposition at that degree of freedom, which copies the source DM's explicit
per-rank ownership ranges rather than letting PETSc re-decide the split. That is
what lets a pointwise loop read a source field and write an accumulator at the same
index.

**Accumulation requires no ghost exchange.** Each update reads a field value at an
owned point and writes the accumulator at that same point, with no neighbour
access. Where a later extension does need ghosts, it must go through
`UpdateLocalGhosts` rather than a hand-rolled exchange.

**The constraint to know before extending output.** `UpdateLocalGhosts` resolves
vectors through a view built on compile-time `UserCtx` offsets, and the nodal
averaging kernel maps fixed field names onto fixed members. Per-window accumulators
are config-counted and have no compile-time offset, so neither accepts them
directly. Derived results are therefore staged through two catalogued scratch
fields, written, converted, and copied out before the next reuses the buffer — the
same pattern the corner-staging buffers already establish. Lifting the constraint
properly means a view bindable to an explicitly supplied vector pair, after which
both nodal statistics output and accumulator field logging work unchanged. See
@ref 60_Field_Statistics_Planned_Extensions.

@section p58_validation_sec 11. How It Is Validated

An extension to this subsystem is expected to extend the suite that covers it, so
the map from behaviour to test is recorded here.

| Behaviour | Covered by |
| --- | --- |
| Moment and co-moment accuracy, zero variance on a constant field, high-mean precision, parallel merge | `tests/c/test_statistics_moments.c` |
| Cadence, stride, start and end clipping, overlapping windows, duplicate-event rejection, identity hashing | `tests/c/test_statistics_window.c` |
| Storage shape, accumulation, covariances, fixed and moving masks, derived quantities, spatial means | `tests/c/test_statistics_accumulator.c` |
| Cell, node, and face layouts; periodic duplicate-plane handling; fluid-mask threshold | `tests/c/test_statistics_target.c` |
| Control resolution, malformed configuration, console cadence independence | `tests/c/test_statistics_config.c` |
| Checkpoint round trip, continuation guards, payload validation, absence when disabled | `tests/c/test_io.c` |
| Console snapshot gate and content, silence when disabled, startup banner in every state | `tests/c/test_logging.c` |
| End-to-end accumulation and derived output, silence below `LOG_INFO`, equivalence across a changed MPI rank count | `tests/smoke/run_smoke.sh` |

Two properties are asserted in specific ways because a weaker test would pass
without meaning anything.

**Variable-step quadrature is asserted exactly.** The analytic variable-step tests
assert the right-rectangle weights exactly rather than within a loose tolerance,
because at small sample counts the difference between candidate rules is precisely
what a tolerance would hide.

**Rank-count equivalence uses an analytically prescribed field.** Accumulation is
pointwise, so the only things a decomposition change can disturb are the target
plan, the boundary and ghost exclusions, and the payload ordering. Comparing
committed payloads byte for byte tests exactly those, but only if the field itself
is identical across rank counts — which a solved field is not, because its Krylov
reductions differ in the last bits. The per-point count payload varies spatially,
so a rank boundary mistaken for a domain boundary changes it.

**Restart equivalence takes its tolerance from a measured floor.** Restart is not
bit-exact: re-applying boundary conditions on load perturbs the outlet, and the
offset is invariant to solver tolerance. See @ref 29_Maintenance_Backlog. That
floor was measured on a laminar case; in turbulent runs trajectories separate
exponentially, so equivalence must be tested on converged window means rather than
field by field.

Known gaps are recorded in @ref p29_stats_validation_sec rather than here, so that
this table describes only what is actually covered.

@section p58_related_sec 12. Related Pages

- @subpage 09_Monitor_Reference — the `field_statistics` and console-cadence keys
- @subpage 10_Post_Processing_Reference — the derived-statistics recipe
- @subpage 60_Field_Statistics_Planned_Extensions — spatial targets, further products, and histories
- @subpage 56_Field_Identity_and_Layout_Catalog — typed field identity and layout
- @subpage 28_IEM_and_Statistical_Averaging — particle mixing and particle statistics
- @subpage 29_Maintenance_Backlog — measured restart fidelity floor and open validation items
- @subpage 52_Run_Lifecycle_Guide — committed checkpoint bundles and restart modes
