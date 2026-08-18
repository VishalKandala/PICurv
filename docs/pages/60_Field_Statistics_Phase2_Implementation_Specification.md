@page 60_Field_Statistics_Phase2_Implementation_Specification Field Statistics Phase 2 Implementation Specification

@anchor _Field_Statistics_Phase2_Implementation_Specification

This page is the implementation specification for Phase 2 of the field-statistics
pipeline. @ref 58_Turbulence_Statistics_Pipeline_Specification remains the
authoritative scientific and architectural design; this page settles the items that
page deliberately deferred, fixes the exact contracts, and defines the
implementation and acceptance order.

Implementation status: Stages 1 through 6 are implemented. `monitor.yml ->
field_statistics` is validated, serialized into the generated `control`, resolved in
C, accumulated each accepted step, and carried across a restart; `post.yml ->
field_statistics` derives Reynolds stresses, RMS, turbulent kinetic energy, and
fluxes from a saved window and writes them as VTK. What remains is Stage 7, which
exposes both blocks in the shipped templates, adds a worked example, and updates the
reference documentation. Until then the keys work but are documented only on this
page, so they are usable without yet being a released contract.

@tableofcontents

@section p60_scope_sec 1. Scope and Relationship to Page 58

Page 58 §17 requires that Phase 2 receive its own detailed spec and approval before
implementation. Page 58 §4 additionally requires that Python validation, control
serialization, C resolution, online state, checkpoint continuation, and
postprocessing arrive as one usable contract, and that no inactive keys be accepted
early merely to reserve syntax. Both constraints shape this page: every contract is
settled here in full, and none of it is exposed to users until the whole path works.

Where this page and page 58 differ in detail, this page governs the implementation
and page 58 is amended to match. Two such amendments are required and are recorded
in @ref p60_amendments_sec.

@section p60_settled_sec 2. Decisions Settled by This Page

Page 58 §6 states that "the precise endpoint quadrature rule and closing-event
behavior must be fixed in the Phase 2 implementation spec and validated with
analytic variable-step signals before the YAML is released." That is settled here.

@subsection p60_quadrature_sub 2.1 Endpoint Quadrature: Right-Rectangle

Physical-time weighting uses the **right-rectangle rule**. An accepted state at
\f$t_i\f$ carries weight \f$w_i = t_i - t_{i-1}\f$, the interval **ending** at it.

The rationale is that the accuracy difference from the trapezoidal rule is an
endpoint artifact, not an accumulating bias. For uniform \f$\Delta t\f$ the two
rules differ by exactly

\f[
\sum_{i=1}^{N} f_i \Delta t \;-\;
\Big[\tfrac{f_0\Delta t}{2} + \sum_{i=1}^{N-1} f_i \Delta t + \tfrac{f_N \Delta t}{2}\Big]
= \frac{(f_N - f_0)\,\Delta t}{2},
\f]

a single endpoint term independent of \f$N\f$. Divided by the window length
\f$T = N\Delta t\f$, the difference in the mean is \f$(f_N - f_0)/(2N)\f$, which
decays as \f$1/N\f$ regardless of timestep size. For a converged turbulence average
this is several orders of magnitude below the statistical sampling error
\f$\sigma/\sqrt{M}\f$ over \f$M\f$ independent samples.

The decisive engineering argument is state cost. Right-rectangle weights depend only
on the past, so the only cross-restart state is the last accepted time, a scalar.
Any trapezoidal or centered rule assigns sample \f$i\f$ a weight depending on
\f$t_{i+1}\f$; because the weighted Welford update is nonlinear in the weight,
retroactively adding weight to a past sample requires that sample's value, forcing a
full pending field snapshot per accumulated field to be held, checkpointed, and
restored exactly. That roughly doubles statistics state and adds a restart failure
mode, in exchange for a correction far below the noise floor.

Right-rectangle is also the rule page 58 §6 already describes in prose: physical-time
weighting "uses the represented interval between accepted states."

@subsection p60_closing_sub 2.2 Closing-Event Behavior

A bounded window's final accepted state has its weight clipped to
\f$t_{\text{end}} - t_{\text{last}}\f$, where \f$t_{\text{end}}\f$ is the requested
`end_time` and \f$t_{\text{last}}\f$ is the previous accepted state. Start and end
clip by the identical rule; there is no special case at either boundary.

@subsection p60_initial_sub 2.3 Initial-State Inclusion

There is no `include_initial` user control. Its behavior is determined by the
interval convention and is identical under both weightings.

A window represents the half-open intervals \f$(t_{i-1}, t_i]\f$. Each accepted state
carries the interval ending at it, measured from the previous accepted state, or from
the window's effective start for the first. **A state representing a zero-length
interval is not a sample.**

Consequently a due state landing exactly on `start_time` anchors the interval origin
and is not a sample, under either weighting; and a due state after `start_time` is a
sample carrying \f$t_1 - t_{\text{start}}\f$, so no sample is wasted and the
represented interval is exactly \f$[t_{\text{start}}, t_{\text{end}}]\f$.

The property this protects is that **the two weightings agree on a constant-`dt`
run**:

\f[
\langle f\rangle_{\text{time}}
= \frac{\sum_{i=1}^{N} f_i \Delta t}{N\Delta t}
= \frac{f_1 + \cdots + f_N}{N}
= \langle f\rangle_{\text{sample}} .
\f]

Counting the initial state as a sample under `sample` weighting only would make the
two weightings disagree by \f$O(1/N)\f$ on a fixed-timestep run for no physical
reason. That is engineered out rather than documented around.

@section p60_yaml_sec 3. User Contract

The accepted location is `monitor.yml -> field_statistics`. This is the complete
Phase 2 surface:

```yaml
field_statistics:
  enabled: true
  windows:
    - name: production
      start_time: 50.0
      end_time: 250.0          # optional; window is open-ended if omitted
      weighting: physical_time # or: sample
      step_cadence: 1          # exactly one of step_cadence / time_cadence
      fields:
        - field: Ucat
          moments: [first, second]
        - field: Psi
          moments: [first, second]
      covariances:
        - [Ucat, Psi]
```

There is no `target:` block in Phase 2. Pointwise mapping and the `fluid` mask are
the only behaviors, so exposing single-valued `kind` and `mask` keys would reserve
syntax in violation of page 58 §4. Phase 3 introduces `target:` with `kind` and
`mask` as **optional** keys defaulting to Phase 2 behavior, so no configuration
written against this page breaks.

@subsection p60_validation_rules_sub 3.1 Python Validation Rules

Python rejects, with a specific message naming the window:

- duplicate window names within one `monitor.yml`;
- `start_time` greater than or equal to `end_time`;
- neither or both of `step_cadence` and `time_cadence`;
- non-positive `step_cadence` (must be a positive integer) or `time_cadence`
  (must be a positive real);
- `weighting` outside `{sample, physical_time}`;
- an empty `fields` list, or an empty `moments` list for any field;
- `moments` entries outside `{first, second}`;
- a field name that does not resolve in the Eulerian field catalog;
- a field whose subsystem is disabled for the run (for example `Psi` or
  `ParticleCount` when `particles.count` is zero, `Nu_t` when LES is off);
- a covariance pair whose members are not both present in `fields` with at least
  `first`, since the means are required to center against;
- a covariance pair that is not location-compatible (see @ref p60_products_sub); and
- a covariance pair naming the same field twice, which is the self-product already
  requested through `moments: [second]`.

@section p60_window_sec 4. Window Semantics, Scheduling, and Weighting

A window is `pending`, `active`, or `complete`, per page 58 §6. It records requested
and effective bounds, sample count, total weight, represented physical time, schedule
position, last accepted event key, last accepted time, resolved definition hash, and
restart lineage.

@subsection p60_cadence_sub 4.1 Cadence

Exactly one cadence is selected.

- `step_cadence: n` accepts a completed state every `n` steps, counted from the
  window's effective start.
- `time_cadence: dt_s` places nominal targets on the absolute grid
  \f$t_{\text{start}} + k\,\mathit{dt_s}\f$ and accepts the **first completed state at
  or past** each target.

Targets are placed on the nominal absolute grid, not relative to the last accepted
time, so the schedule cannot drift as `dt` varies.

When a single large step overshoots several nominal targets, that state is accepted
**exactly once** and the schedule position advances past the current time. No time is
lost or double counted, because under right-rectangle weighting the accepted state's
weight is the actual elapsed interval since the last accepted state rather than the
nominal cadence. This self-correction is a direct consequence of
@ref p60_quadrature_sub; under a nominal-cadence or centered weighting a skipped
target would silently mis-weight the average.

@subsection p60_bounds_sub 4.2 Effective Bounds and Edge Cases

- `effective_start` is normally the requested `start_time`. When a window is first
  observed later than that — for example a `--restart-from` resuming past
  `start_time` with fresh statistics — `effective_start` is the first observed state
  time, and the window records that its represented interval is shorter than
  requested. A window must never claim represented time it did not observe.
- `effective_end` is the requested `end_time` for a completed bounded window, the
  last accepted state time for an open window, or the run's final state time when a
  bounded window is still active at shutdown.
- A window that reaches `complete` with `sample_count == 0` is a **hard error**, not
  a silently empty result for postprocessing to divide by.
- With a coarse `step_cadence`, the first sample legitimately represents up to `n`
  steps back to `effective_start`. This is the same rule as every other sample and is
  stated only so it is not mistaken for a defect.
- An active off-schedule state changes no scientific state, per page 58 §6.

@section p60_accum_sec 5. Accumulator State and Products

@subsection p60_state_sub 5.1 Per-Point State

For each accumulated scalar quantity a point stores sample count, total weight
\f$W\f$, squared-weight sum \f$W_2\f$, mean \f$\mu\f$, and centered second-moment sum
\f$M_2\f$. The weighted update for sample \f$x\f$ with weight \f$w > 0\f$ is page 58
§7:

\f[
W' = W+w,\qquad \delta=x-\mu,\qquad
\mu'=\mu+\frac{w}{W'}\delta,\qquad
M_2'=M_2+w\,\delta(x-\mu').
\f]

The compatible co-moment update for a pair \f$(x, y)\f$ with means
\f$\mu_x, \mu_y\f$ is

\f[
C' = C + w\,(x-\mu_x)\,(y-\mu_y'),
\f]

evaluated with the pre-update \f$\mu_x\f$ and the post-update \f$\mu_y'\f$, which is
the form that stays consistent with the \f$M_2\f$ update when \f$x = y\f$.

\f$W_2\f$ is retained for effective-sample-size and uncertainty estimates and is
required because physical-time weighting produces unequal weights.

@subsection p60_products_sub 5.2 Products and Component Order

Page 58 §7 requires that component order be explicit and never implicit.

- **Scalar `first`** stores \f$\mu\f$. **Scalar `second`** additionally stores
  \f$M_2\f$, one component.
- **Vector `first`** stores \f$\mu\f$ per component, `dof` components.
- **Vector `second`** stores the full symmetric self-product. For a three-component
  field this is **six** components in upper-triangular row-major order:

  \f[
  (xx,\; xy,\; xz,\; yy,\; yz,\; zz)
  \f]

  corresponding to index pairs \f$(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)\f$. These
  are the centered Reynolds-stress sums; `moments: [second]` on `Ucat` requests all
  six. **Reynolds stresses are not cross-field products and require no covariance
  entry.**
- **Cross-field covariance** between a `dof = 3` vector and a `dof = 1` scalar stores
  **three** components in field component order \f$(x\!\cdot\!s,\; y\!\cdot\!s,\;
  z\!\cdot\!s)\f$. This is the turbulent scalar flux.

Location compatibility for covariance requires both members resolve to the same
`FieldDescriptor` layout. In Phase 2 that means cell-centered with cell-centered.
`Ucont` is component-staggered, so pairs involving it are rejected, consistent with
`Ucont` face-flux products being an explicit Phase 3 non-goal.

Only vector-scalar and scalar-scalar covariance pairs are accepted in Phase 2.
Vector-vector cross products are rejected.

@subsection p60_merge_sub 5.3 Parallel Merge

The stable weighted parallel-combination formula is implemented alongside the
sequential update, because page 58 §7 requires centered state be mergeable for offline
window merging and later spatial reductions. For partitions \f$A, B\f$:

\f[
\delta = \mu_B - \mu_A,\qquad W = W_A + W_B,\qquad
\mu = \mu_A + \frac{W_B}{W}\delta,\qquad
M_2 = M_{2,A} + M_{2,B} + \frac{W_A W_B}{W}\delta^2 .
\f]

Pointwise PETSc vector updates are local and do not themselves require the merge, but
it is unit tested in Stage 1 so later reductions inherit a verified kernel.

@section p60_target_sec 6. Spatial Target, Layout, and Mask

`SpatialTargetPlan` exists as a C abstraction with only the pointwise identity
mapping implemented, so Phase 3 extends rather than retrofits.

Field resolution reuses `FieldGetDescriptor` and `FieldGetView`
([field_catalog.h](../../include/field_catalog.h)) for identity and layout, and the
block's existing `DMDALocalInfo` for ownership. No new field lookup, layout table,
or range convention is introduced.

`GetOwnedCellRange` is deliberately not called. It reports cells by origin node
index and covers cell-centered layouts only, whereas the target plan works in
field-storage indices under the solver's shifted convention and must also resolve
node and face layouts. The two describe the same cells in different index spaces:
the cell with origin `j` is stored at index `j + 1`. That relationship is asserted
by a unit test rather than left implicit, so the plan and the established helper
cannot drift apart.

Iteration must exclude **both** categories, which are distinct:

- PETSc halo/ghost storage, and
- solver-layout boundary, dummy, and duplicate-periodic indices.

The default and only Phase 2 mask is `fluid`, defined by `Nvert < 0.1`. Because
`Nvert` can change when immersed bodies move, the mask is treated as potentially
moving: **every point keeps its own valid count and valid weight**, and a
`valid_fraction` diagnostic is derived as the point's valid weight divided by the
window's total weight. A point that is never valid retains zero weight and is
reported rather than silently producing a divide-by-zero in postprocessing.

@section p60_hash_sec 7. Window Identity Hash

At startup C resolves the window definition once and computes a SHA-256 over its
canonical serialization using [checksum.h](../../include/checksum.h). The hash is
internal checkpoint metadata, not a user option.

Hashed, in this fixed order:

1. window name;
2. `start_time`;
3. weighting mode;
4. cadence kind and value;
5. resolved field entries, ordered by `FieldId`, each as canonical field name plus
   its requested moment set in fixed `first,second` order;
6. resolved covariance pairs, each ordered by `FieldId` within the pair and the list
   ordered by first then second member;
7. mask identity; and
8. target kind and layout semantics.

Each of those eight groups is additionally hashed on its own and the truncated group
digests are stored beside the full one, which is what lets a mismatch name the
property that changed rather than only reporting that two digests differ. Group 7,
the mask, hashes a constant in Phase 2 because exactly one mask is resolved; it
exists now so a Phase 3 mask key extends that group instead of renumbering the groups
after it and invalidating every saved window.

Field and covariance entries are serialized in catalog order, and each covariance
pair is canonicalized within itself, so a reordered but otherwise identical
configuration continues rather than being rejected. Covariance is symmetric, so
listing its members either way describes the same quantity.

Explicitly **excluded**: `end_time`, so a bounded window may be extended forward per
page 58 §10; and `enabled`, so switching statistics off and on does not invalidate
saved state; and the order in which fields and covariance pairs were listed.

Extension is accepted only when the new end moves forward, does not exceed the
configured simulation horizon, and continuation occurs without an unsampled gap after
the former end. Shortening a window, or changing any hashed property, requires a new
window name.

The gap condition is enforced against the checkpoint's own physical time: reopening a
window that closed at \f$t_{\mathrm{end}}\f$ from a bundle written at
\f$t_{\mathrm{ckpt}}\f$ is refused once \f$t_{\mathrm{ckpt}} - t_{\mathrm{end}}\f$
exceeds one timestep. This is a direct consequence of right-rectangle weighting: the
first sample after the former end carries the whole interval back to it, so a gap
would be weighted as though the window had observed it. One step of slack is allowed
because the closing state's own clipping already leaves up to that much between the
window's end and the state that closed it.

`bounded` is excluded from the hash alongside `end_time`, since it only says whether
an end exists. The dangerous direction — turning an open window into one that ends
before the time it already represents — is caught by the shortening check rather than
by the hash.

@section p60_ingress_sec 8. Control Serialization and C Ingress

The generated master `control` remains the only C-ingress artifact. No sidecar file
is introduced, per page 58 §2 and §4.

@subsection p60_control_sub 8.1 Option Families

Because a window list is variable-arity, option names are constructed rather than
literal. The scheme mirrors the checkpoint manifest's indexed style:

```
-field_statistics_enabled                        true
-field_statistics_window_count                   2
-field_statistics_window_0_name                  production
-field_statistics_window_0_start_time            50.0
-field_statistics_window_0_end_time              250.0
-field_statistics_window_0_weighting             physical_time
-field_statistics_window_0_step_cadence          1
-field_statistics_window_0_field_count           2
-field_statistics_window_0_field_0_name          Ucat
-field_statistics_window_0_field_0_moments       first,second
-field_statistics_window_0_covariance_count      1
-field_statistics_window_0_covariance_0          Ucat,Psi
```

Python emits these alongside the existing monitor flags, through a
`normalize_field_statistics_config` and `resolve_field_statistics_flags` pair modeled
on the solution-monitoring pair.

@subsection p60_audit_sub 8.2 Ingress Audit Extension

[audit_ingress.py](../../tests/tooling/audit_ingress.py) currently matches only
`PetscOptionsGet*(NULL, NULL, "-literal", ...)`. Every control read in the codebase
today is such a literal. Constructed names would be invisible to the audit **in any
file**, so relocating the parse site does not address this; the audit itself must be
extended.

The extension relies on the fact that a constructed name still has a literal format
specifier in the source:

```c
PetscCall(PetscSNPrintf(key, sizeof(key), "-field_statistics_window_%d_field_%d_moments", w, f));
PetscCall(PetscOptionsGetString(NULL, NULL, key, buf, sizeof(buf), &found));
```

The audit therefore gains three behaviors:

1. extract `PetscSNPrintf` format strings beginning with `-` as **family patterns**;
2. compare them against a new `known_petsc_option_families` manifest key using the
   same exact-set-equality rule already applied to `known_petsc_options`; and
3. **fail** on any `PetscOptionsGet*(NULL, NULL, <non-literal>, ...)` whose name
   variable is not traceable to a declared family.

Behavior 3 is what preserves the guarantee rather than widening a loophole. After the
extension the invariant is that every option name reaching C is either a declared
literal or a declared family.

Note that the checkpoint manifest's indexed reads pass a private `PetscOptions`
object rather than `NULL` and are correctly outside this audit's scope: they are
checkpoint state, not user ingress.

@subsection p60_site_sub 8.3 Parse Site

Statistics option parsing lives in a new `src/statistics_config.c`, not in
`src/setup.c`, which already exceeds 3,700 lines. The file is added to the audit's
scan paths and to the manifest's `sources` key so the code and the audit cannot
drift apart.

@section p60_checkpoint_sec 9. Checkpoint Namespace

Window state is written into the existing committed bundle by extending
`WriteCheckpointBundle`, `WriteCheckpointManifest`, and `ValidateCheckpointBundle`
rather than adding a parallel coordinator. Payloads go through the generic vector
writer; no second binary serializer is introduced, per page 58 §12.

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

Two properties of that tree are load bearing.

Payloads are **block scoped**, exactly as Eulerian payloads are, so a multiblock run
keeps one accumulator tree per block under each window.

Each product and co-moment is **one payload carrying all of its components**, not one
payload per component. A symmetric second-order tensor is a single physical object,
and splitting it would cost six memory streams in the per-step accumulation loop
instead of one cache line, and six collective gathers per checkpoint instead of one.
Splitting also fails to generalize: it works for a three-vector only because six
happens to equal three plus three, and no comparable split exists for a dof-2 field's
three-component product or a third moment's ten.

Component counts that neither `da` nor `fda` provides are carried by a DM mirroring
the block decomposition at that degree of freedom, created through
`CreateCompatibleBlockDM` ([grid.c](../../src/grid.c)). This is the pattern the
`fda2` slot on `UserCtx` already establishes for dof-2 RANS fields; the helper is
what that slot needs as well. Mirroring copies the source DM's explicit per-rank
ownership ranges rather than letting PETSc re-decide the split, so the decomposition
is identical and a pointwise loop can read a source field and write an accumulator at
the same index.

Payload names come from the catalogued field name plus a fixed role suffix, so they
are stable across runs, rank counts, and configuration reorderings. A field for which
no second moment was requested contributes no product payload at all.

Scalar per-window metadata is recorded in `checkpoint.meta`:

```
-checkpoint_statistics_window_count            1
-checkpoint_statistics_window_0_name           production
-checkpoint_statistics_window_0_hash           <sha256>
-checkpoint_statistics_window_0_hash_groups    <8 truncated group digests>
-checkpoint_statistics_window_0_state          active
-checkpoint_statistics_window_0_sample_count   1500
-checkpoint_statistics_window_0_total_weight   1.5
-checkpoint_statistics_window_0_represented_time 1.5
-checkpoint_statistics_window_0_last_accepted_time 51.5
-checkpoint_statistics_window_0_effective_start  50.0
-checkpoint_statistics_window_0_effective_end    51.5
-checkpoint_statistics_window_0_activation_step  1000
-checkpoint_statistics_window_0_last_event_step  1500
-checkpoint_statistics_window_0_next_time_target 30
-checkpoint_statistics_window_0_restart_count    2
```

The window count is written even when no window is configured, so a restart can tell
an absent window list from an unreadable bundle.

Everything after `effective_start` is **schedule state**, and omitting it is not a
cosmetic loss. `activation_step` anchors the step stride, `next_time_target` anchors
the absolute time grid, and `last_event_step` is the duplicate-event guard. A restart
that restored only the reported quantities would resume on a silently shifted
schedule while continuing to report a correct-looking sample count.

`restart_count` is the concrete form of the restart lineage page 58 §6 requires: the
number of restart segments the state descends from. Physical case identity is carried
by the bundle's existing `-checkpoint_geometry_sha256`, so no second identity scheme
is introduced.

`hash_groups` holds one truncated digest per hashed property group, comma separated in
the order of @ref p60_hash_sec. The saved definition itself is never stored, so
without these a mismatch could only report that two opaque digests differ; with them
the restart names the property that changed.

Every statistics payload appears in the existing payload inventory with its field
name, block, layout, component count, logical type, global size, encoding, relative
path, and byte size, and is validated for existence and size before the bundle is
accepted. That validation needs no new code: entering the inventory is what makes a
payload subject to the loop the validator already runs. Because that is easy to lose
silently, acceptance damages a statistics payload and requires the bundle to be
rejected, rather than asserting the property by inspection. The inventory field name is qualified by window, as `<window>/<payload>`,
so two windows accumulating the same field stay distinguishable to anything reading
the manifest alone.

`last_accepted_time` is the only quadrature state required across a restart. That it
is a single scalar rather than a field snapshot is the direct consequence of choosing
right-rectangle in @ref p60_quadrature_sub.

@section p60_restart_sec 10. Restart and Continuation

- `--continue` resumes in the same run directory and continues a window when the
  resolved definition hash matches. An allowed `end_time` extension is checked
  separately per @ref p60_hash_sec.
- `--restart-from` continues statistics **only** when the physical case identity and
  resolved definition hash still match and continuation was explicitly requested.
  Otherwise the flow may restart while statistics reset under a new window name.
- Missing required statistics state for a requested continuation is **fatal**. It is
  never silently zeroed.
- A hash mismatch is fatal with a message naming the window and the first differing
  resolved property.
- A different MPI rank count is permitted. Statistics payloads use the same natural
  ordering as Eulerian payloads, so they are decomposition independent.
- Shortening a window is refused: a requested `end_time` earlier than the time the
  saved state already represents would discard samples the metadata still counts.
- Restoration happens in `InitializeEulerianState`
  ([initialcondition.c](../../src/initialcondition.c)), from the same bundle the flow
  state came from, and independently of the Eulerian source, so an analytical restart
  still continues a window that was accumulating before it.

@section p60_post_sec 11. Postprocessing Contract

Recipes are spelled at `post.yml -> field_statistics`, kept separate from the existing
particle `statistics_pipeline`.

Two structural changes are required in [postprocessor.c](../../src/postprocessor.c):

1. **Decouple statistics from the particle stage.** `needs_particle_stage` is
   currently true whenever any statistics pipeline is non-empty, and the run
   hard-errors when `np == 0`. Task dispatch becomes capability based: a
   `field_statistics` recipe declares `needs_statistics_checkpoint` instead.
2. **Add a pipeline beside the particle one.** `FieldStatisticsPipeline(user, pps, ti)`
   mirrors `GlobalStatisticsPipeline`, which is the established shape for a
   per-timestep stage that writes its own output family. Statistics fields are a
   fourth output family alongside Eulerian `.vts`, particle `.vtp`, and particle
   statistics `.csv`, so they share the VTK assembly helpers but not the Eulerian
   writer's fixed field list.

**Where the pipeline runs.** `FieldStatisticsPipeline` is called once per processed
step from the post loop, beside the particle `GlobalStatisticsPipeline` it mirrors. A
recipe spanning several steps therefore produces a convergence history rather than
one picture repeated; pinning `source_step` collapses it to a single bundle. This
supersedes the earlier requirement for a stage outside the loop: the loop's own
cadence is what makes the history possible, and running inside it costs nothing when
a recipe covers one step.

Each requested window is written to its own file. That is forced rather than
preferred: `MAX_POINT_DATA_FIELDS` is 20 and one window carrying every output already
produces 15 fields, so two windows cannot share a file.

Outputs come in two formats. `vtk` writes the derived fields; `csv` appends one row
per processed step recording sample count, total weight, represented time, the
per-point valid-fraction range, and the mean turbulent kinetic energy. The CSV is the
artifact that answers whether a window has run long enough, which no single field
snapshot can.

That mean is taken over the **fluid cells the window actually sampled**, not over the
stored vector. A derived field is zero outside the target domain and at any point the
mask excluded, and those zeros are absences rather than measurements: averaging over
them scales the result down by whatever fraction of the vector the window never
covered. For a cell-centred field the targeted span is exactly the real cells, since
the first and last index of each dimension are dummy layers, so no cell count is
computed from grid dimensions anywhere. The `valid_fraction_min` column in the same
row is what reveals whether an immersed body made the sampled set smaller than the
targeted one; at one they are the same set.

A window that has not yet accumulated a sample at a given step is skipped with a note
rather than treated as an error, because a recipe covering a whole run legitimately
reaches bundles from before that window began.

**Derived statistics are not dimensionalized**, even when
`global_operations.dimensionalize` is set. A Reynolds stress scales as velocity
squared and a co-moment as the product of two different scales; the existing
per-field scaling table expresses neither, so applying a velocity scale would be
wrong rather than merely incomplete. Dimensional derived statistics require scaling
rules of their own.

**How the postprocessor learns the windows.** It does not re-describe them. The
postprocessor is launched with the run's own `-control_file`, and
`CreateSimulationContext` resolves it through the same `ParseFieldStatisticsConfig`
the solver uses, so the window list, the accumulator storage, and
`RestoreFieldStatisticsState` all apply unchanged. A recipe therefore names a window
and nothing more, and no second loader exists to disagree with the first.

**Resolution of the nodal-output constraint.** Derived results are staged through two
catalogued fields, `PostScalar` and `PostVector`, created beside the corner-staging
workspace. A derived quantity is written there, reached by catalogued name through
`UpdateLocalGhosts` and `ComputeNodalAverage`, and copied out to VTK before the next
one reuses the buffer. This works *with* the compile-time offset constraint rather
than against it: the alternative, a `FieldView` bound to explicit vectors, would still
need the same local and nodal scratch vectors and would additionally have to move
`dof`, layout, and sync class into the view for every consumer. The precedent is
already in the catalog, where `CellScalarAtCorner` and `CellVectorAtCorner` are
staging buffers rather than simulation state.

Accumulator field logging is reachable the same way, by staging a window's vector
before logging it, so no separate diagnostic path is required.

**Known constraint on nodal output.** `UpdateLocalGhosts` resolves its vectors
through `FieldGetView`, which reads compile-time `UserCtx` offsets stored in the
catalog, and `ComputeNodalAverage`
([postprocessing_kernels.c](../../src/postprocessing_kernels.c)) maps field names
onto fixed nodal members. Per-window accumulators are config-counted and therefore
have no compile-time offset, so neither routine accepts them as they stand. This
must be resolved before nodal statistics output, and the resolution must extend
the existing mechanism rather than duplicate it: the preferred route is a
`FieldView` that can be bound to an explicitly supplied vector pair, after which
the whole ghost path applies unchanged. Emitting cell-centred statistics only is
the fallback if that extension is judged out of scope.

Derived quantities are normalized in one place from centered state:

\f[
R_{ij} = \frac{C_{ij}}{W},\qquad
\mathrm{RMS}_i = \sqrt{R_{ii}},\qquad
k = \tfrac{1}{2}\left(R_{xx}+R_{yy}+R_{zz}\right),\qquad
\overline{u_i'\psi'} = \frac{C_{i\psi}}{W}.
\f]

Small negative variances arising from floating-point cancellation are clamped **only**
when taking a square root and only within a documented tolerance. Stored state is
never mutated. Existing dimensionalization, nodal-conversion, VTK, and CSV kernels are
reused.

@section p60_logging_sec 12. Runtime Logging and Console Monitoring

The pipeline reports through the existing logging infrastructure rather than
printing directly, and follows the conventions already established for particle
monitoring. Every stage carries its logging obligation; it is not deferred to a
cleanup pass.

@subsection p60_console_sub 12.1 Console Monitoring

Field statistics gains an optional periodic console snapshot, modelled on the
particle console snapshot so that operators reading a running job see a familiar
shape. Each element mirrors an existing one:

| Concern | Particle (existing) | Statistics (new) |
| --- | --- | --- |
| YAML key | `io.particle_console_output_frequency` | `io.statistics_console_output_frequency` |
| Control flag | `-particle_console_output_freq` | `-statistics_console_output_freq` |
| `SimCtx` field | `particleConsoleOutputFreq` | `statisticsConsoleOutputFreq` |
| Enabled predicate | `IsParticleConsoleSnapshotEnabled` | `IsStatisticsConsoleSnapshotEnabled` |
| Cadence predicate | `ShouldEmitPeriodicParticleConsoleSnapshot` | `ShouldEmitPeriodicStatisticsConsoleSnapshot` |
| Emitter | `EmitParticleConsoleSnapshot` | `EmitStatisticsConsoleSnapshot` |

The key sits under `io` rather than inside `field_statistics` because that is
where the existing output cadences live, and because it controls reporting rather
than science: changing it must not alter any accumulated result, and it is
therefore excluded from the window definition hash.

The enabled predicate uses the same three-part gate as the particle one: the
subsystem is active, the configured cadence is positive, and the effective log
level is at least `LOG_INFO`. A zero cadence disables the snapshot without
disabling accumulation, exactly as a zero particle cadence does.

Each snapshot reports, per active window: name, state, accumulated sample count,
total weight, represented physical time, progress toward a bounded end, and the
range of per-point valid fraction as a mask-health indicator. These are all
window-level scalars. **The console snapshot never dumps field data**; that is
what the debug-level field loggers are for.

The mask-health indicator is the range of per-point valid fraction, where a point's
fraction is its own accumulated count divided by the window's accepted-sample count.
A minimum of one means every targeted point saw every state; a minimum below one means
a moving body excluded some points for part of the window; a minimum of zero means
some point contributed nothing and its mean is undefined. That last case cannot be
seen by reading the mean field, because an untouched point holds a zero that looks
like a legitimate value, so it is also reported once as a warning when the window
completes.

Computing the range is a collective reduction, so it runs only on the console cadence
and at window completion, never per step.

The startup banner gains a line in the shape the existing `Field/Restart Cadence`
entry uses, printing the configured cadence or `DISABLED`, so a run's log records
whether monitoring was active.

@subsection p60_levels_sub 12.2 Log Levels

Statistics code uses `LOG_ALLOW` at the levels the rest of the codebase uses:

- `LOG_INFO` for once-per-run and once-per-window lifecycle events: resolved
  window definitions at startup, a window becoming active, a window completing
  with its final sample count and represented time, and checkpoint save or
  restore of window state.
- `LOG_DEBUG` for per-event bookkeeping: a state accepted with its weight and the
  window's running totals, and a state an active window did not sample, with the
  last accepted time that explains it. Both are expected outcomes, which is why
  neither is a warning.
- `LOG_WARNING` reserved for conditions an operator must act on: a window completing
  without accepting any sample, and a window completing with points whose valid
  fraction is zero. Expected behaviour is never logged as a warning.
- `LOG_TRACE` for per-point or per-block detail that would otherwise flood a run.

Accumulation is on the timestep path, so `LOG_DEBUG` and below must not perform
collective reductions unless the level is already active.

`LOG_ALLOW` is gated by the per-function allow list as everywhere else in the
codebase, whose default admits only `main` and `CreateSimulationContext`. Statistics
lifecycle lines therefore appear only when their functions are named in
`logging.enabled_functions`, exactly as any other module's do. This is deliberate and
is why the console snapshot uses `LOG` instead: monitoring an accumulating window is
an operator-facing report that must not depend on naming an internal function, while
the lifecycle lines are diagnostics.

@subsection p60_fieldid_sub 12.3 Field Identity in Logging

`LOG_FIELD_MIN_MAX` and `LOG_FIELD_ANATOMY` are `FieldId`-routed, resolving the
DM, degree of freedom, and layout through `FieldGetView`. **Source fields being
accumulated are catalogued, so they use these directly.** Accumulator vectors are
config-counted, have no compile-time `UserCtx` offset, and therefore cannot be
addressed by them, which is the same constraint recorded against nodal output in
@ref p60_post_sec. Until a `FieldView` bindable to an explicitly supplied vector
pair exists, accumulator diagnostics report window-level scalars rather than field
anatomy. That extension would make both surfaces work at once, which is a further
argument for it over a statistics-specific logging path.

@section p60_stages_sec 13. Implementation Stages

Each stage is independently testable. No user-facing YAML is accepted until Stage 7.

1. **Centered moment kernels.** `include/statistics_moments.h`,
   `src/statistics_moments.c`. Weighted Welford, co-moment, and parallel merge. Pure:
   no configuration, no PETSc state.
2. **Field resolution, masks, and the pointwise target plan.** `SpatialTargetPlan`
   with identity mapping only; reuse of catalog and range helpers.
3. **Windows, scheduling, and the runloop hook.** Independent PETSc state per named
   window. The completed-state event is inserted in
   [runloop.c](../../src/runloop.c) immediately before `LOG_SOLUTION_CONVERGENCE`,
   after the Lagrangian block and before history rotation and checkpoint output.

   Accumulator vectors are created inside `CreateAndInitializeAllVectors`
   ([setup.c](../../src/setup.c)) and released through the teardown that
   `DestroyUserVectors` participates in. Every vector a run owns is created in that
   one factory; a config-counted array is no exception, since the window count is
   resolved from configuration long before the factory runs.

   The mechanism is the one the convergence state already uses: `PetscCalloc1` for
   the array, then a global vector from the DM that carries the accumulator's
   component count. That inherits the block decomposition, so no accumulator can
   misalign with the field it accumulates. When a future mask makes a window smaller
   than the full field, only the source DM changes; the placement does not.

   Allocation happens once at setup and release once at teardown. The runloop hook
   only updates state that already exists and allocates nothing.

   Window lifecycle logging lands here: resolved definitions at startup, activation,
   and completion, at the levels in @ref p60_levels_sub.

   Accumulation requires **no ghost exchange**. Each update reads a field value at
   an owned point and writes the accumulator at that same point, with no neighbour
   access, so there is nothing to halo-exchange. No `DMGlobalToLocal` call belongs
   anywhere in this stage; where a later stage does need ghosts, it must go through
   `UpdateLocalGhosts` rather than a hand-rolled exchange.
4. **Checkpoint continuation.** Extend the existing coordinator, manifest writer, and
   validator as specified in @ref p60_checkpoint_sec. Payloads go through the same
   `WriteFieldData`/`ReadFieldData` pair Eulerian fields use, so no second binary
   serializer appears and rank-count independence is inherited rather than
   reimplemented.

   Persistence is driven by one enumeration of a window's storage, walked identically
   by the manifest inventory, the writer, and the reader. A payload therefore cannot
   be written under one name and looked for under another.

   The definition hash lands here rather than with ingress, because it is checkpoint
   metadata: it is what a restart compares against, and it operates on a resolved
   definition regardless of how that definition was produced.

   The continuation switch defaults to the run's continue mode. Stage 5 exposes the
   explicit request `--restart-from` requires.
5. **Ingress.** Python schema, normalize/resolve pair, control emission, C resolution
   and hashing in `src/statistics_config.c`, and the audit extension. This stage also
   carries `io.statistics_console_output_frequency` through the same chain, and adds
   the console snapshot predicates and emitter in `logging.c` plus the startup banner
   line, per @ref p60_console_sub.
6. **Postprocessing.** Capability dispatch, the out-of-loop stage, and derived
   normalization.
7. **Release the contract.** Expose the YAML in templates, add a worked example, and
   update the reference documentation. Each gate document currently carries an
   explicit "not yet active" statement that must be **replaced**, not merely deleted.

@section p60_acceptance_sec 14. Validation and Acceptance

Page 58 §16 applies in full. Its items are carried by these suites:

- constant scalar and vector fields yielding exactly zero covariance, analytic signals
  with known sample- and time-weighted moments, and high-mean, low-fluctuation
  precision — [test_statistics_moments.c](../../tests/c/test_statistics_moments.c) and
  [test_statistics_accumulator.c](../../tests/c/test_statistics_accumulator.c);
- variable `dt`, cadence stride, start and end clipping, overlapping windows, and
  duplicate-event rejection — [test_statistics_window.c](../../tests/c/test_statistics_window.c);
- fixed and moving masks with per-point valid weights — the valid-fraction tests in
  [test_statistics_accumulator.c](../../tests/c/test_statistics_accumulator.c), which
  blank a point partway through a window;
- nonperiodic, mixed, and fully periodic domains with no duplicate planes, and cell,
  node, and I/J/K-face layouts — [test_statistics_target.c](../../tests/c/test_statistics_target.c);
- uninterrupted versus restart or continue equivalence — the checkpoint round trip and
  continuation guards in [test_io.c](../../tests/c/test_io.c);
- serial and rank-count-change equivalence — the statistics smoke scenarios in
  [run_smoke.sh](../../tests/smoke/run_smoke.sh), which run one analytically prescribed
  case at one rank and at the harness rank count and require the committed accumulator
  payloads to be byte-identical. An analytic field is used deliberately: it is
  evaluated pointwise, so a solved field's Krylov reductions cannot mask a real
  decomposition dependence in the last bits; and
- unchanged flow, monitoring, logging, profiling, particle, and model behavior when
  field statistics are disabled — the disabled-run smoke variant, together with the
  bundle and console gates in [test_io.c](../../tests/c/test_io.c) and
  [test_logging.c](../../tests/c/test_logging.c).

Logging and monitoring carry their own acceptance items, covered by the statistics
console tests in [test_logging.c](../../tests/c/test_logging.c) and the three smoke
variants:

- the console snapshot is emitted at the configured cadence and **not at all** when
  the cadence is zero or the log level is below `LOG_INFO`, mirroring the existing
  particle-console tests. The below-`LOG_INFO` half is asserted from a solver run
  rather than a unit test, because the effective level is resolved once per process;
- the startup banner reports the configured cadence, or `DISABLED`, in both cases;
- changing the console cadence does not change any accumulated result, which is the
  observable form of its exclusion from the definition hash. The comparison runs where
  the log level admits `LOG_INFO`, so the reporting run genuinely reports; and
- a run with statistics disabled emits no statistics logging on any path.

Two items remain open, and neither is Stage 7 work:

- **Multiblock equivalence** has no coverage, because no multiblock runtime harness
  exists at any level. Establishing one is the prerequisite, and it is not specific to
  statistics.
- **Graceful-shutdown timing** is covered only through the writer it shares with the
  ordinary cadence: the committed bundle at a given step carries that step's
  contribution, which the smoke scenarios assert. That the run loop offers the
  completed state to the windows *before* it reaches the shutdown writer is currently
  guaranteed by the ordering in [runloop.c](../../src/runloop.c) rather than by a test.

Two acceptance constraints are specific to this design:

- **Quadrature tests assert exact weights.** The analytic variable-step tests
  required by page 58 §6 must assert the right-rectangle weights exactly rather than
  within a loose tolerance, because at small sample counts the difference between
  candidate rules is precisely what would otherwise be masked.
- **Restart-equivalence tolerance comes from a measured floor.** Restart is not
  bit-exact: re-applying boundary conditions on load perturbs the outlet, and the
  offset is invariant to solver tolerance. The measured floor is documented in
  @ref 29_Maintenance_Backlog §8 and restart-equivalence tolerances must be derived
  from it, not from machine epsilon. That floor was measured on a laminar case; for
  turbulent runs, trajectories separate exponentially and equivalence must be tested
  on converged window means rather than field by field.

@section p60_amendments_sec 15. Required Amendments to Page 58

1. §4's illustrative YAML shows `include_initial: false` and lists "initial-state
   inclusion is explicit and defaults to false" among settled semantics. Both are
   removed; @ref p60_initial_sub governs. Leaving the key would itself violate §4's
   prohibition on reserving inactive syntax.
2. §4's illustrative YAML shows a `target:` block with `kind` and `mask`. Both are
   deferred to Phase 3 as optional keys; see @ref p60_yaml_sec.

@section p60_nongoals_sec 16. Non-Goals

Deferred to Phase 3 and beyond, each opt-in and built on the same target and product
abstractions: spatial bins, profile layers, named regions, grid and immersed
surfaces, and point probes; explicit `Ucont` face-flux and single-component products;
reduced histories, PDFs, and higher moments, which cannot be recovered from centered
state after the fact; and replacing `ComputeCurrentFlowObservables` with the shared
reducer, which is permitted only after tests prove identical results.

Explicitly not planned at any phase: a legacy `su0`/`su1`/`su2`/`sp` importer, a
dual-write period, or any compatibility workflow.

@subsection p60_naming_sub 15.1 Naming Convention

The pipeline uses two prefixes, split by what a symbol is for rather than where it
lives, so a name can be predicted from its role.

`Picurv*` marks the reusable statistics library: the centered-moment kernels, the
window lifecycle, and the accumulator. These take plain arguments or a block context
and hold no opinion about how a run is configured, and their bare nouns — `Window`,
`MomentState`, `ProductComponentCount` — are generic enough to collide in a single
namespace. This follows `PicurvSHA256*` in [checksum.h](../../include/checksum.h),
which carries the prefix for the same reason.

`FieldStatistics*` marks the run-integration surface: the predicate the runloop,
checkpoint writer, and console monitor share, the per-step driver, and the
post-processing pipeline. These take a `SimCtx` and exist to be called from the
solver and post-processor, so they follow the `Field*` family of the typed catalog
rather than the library prefix. Configuration entry points keep the codebase's
verb-first form, matching `ParsePostProcessingSettings`.

`SpatialTargetPlan*` carries neither, because the noun is already specific enough not
to collide, in the same way `FieldView` and `BoundingBox` need no qualifier.

Loop variables follow the surrounding code: `i`, `j`, `k` for grid traversal, and a
descriptive name for anything indexing a configured list, as the checkpoint writer
already does with `block` and `payload`.

@subsection p60_carried_sub 16.1 Carried Deliberately

Three structural choices are left as they are, recorded here so a later change knows
what was weighed rather than rediscovering it.

**The accumulator module holds two halves.**
[statistics_accumulator.c](../../src/statistics_accumulator.c) carries both the
online accumulation that runs every accepted step and the offline derivation and
spatial reduction that run only in post-processing. They are close to equal in size
and share exactly one thing: the symmetric component-pair table that defines the
order of a three-vector's self-product. Splitting the derivation into its own module
would separate a hot path from a cold one, at the cost of exporting that table across
a file boundary. The moment to do it is when a further product family arrives and the
table needs a home of its own; until then the shared order is a stronger argument for
one file than the size is against it.

**Several centered-moment kernels have no production caller.** `PicurvMomentStateReset`
and `PicurvCoMomentStateReset` are the constructors of their structs.
`PicurvMomentStateMerge` and `PicurvCoMomentStateMerge` implement the stable weighted
parallel combination that page 58 §7 requires centered state to support; pointwise
accumulation never merges, so the first caller will be a spatial reduction in Phase 3.
`PicurvMomentStateVariance` is the only reader of the second moment that
`PicurvMomentStateUpdate` already maintains — the pipeline routes every product
through the co-moment path so the diagonal and off-diagonal share one update, which
leaves the scalar accessor correct but unused. `PicurvMomentStateEffectiveCount`
returns the Kish effective sample size, which is what a physical-time-weighted window
should report instead of a raw sample count when the timestep varies; surfacing it
needs a window-level sum of squared weights, and therefore a new checkpoint field and
restore path. Each is verified by the moment suite, so a future consumer inherits a
tested kernel rather than writing one.

**Comma-separated lists are walked in several places rather than through one helper.**
The configuration parser splits into a token array, the derivation matches tokens
against a keyword table, and the post-processing pipeline iterates and acts directly.
The shapes differ enough that a shared splitter would take a callback or an output
array that most callers discard, and
[postprocessor.c](../../src/postprocessor.c) uses `strtok` throughout, so a bespoke
helper there would read as the exception rather than the rule.

@section p60_related_sec 17. Related Pages

- @ref 58_Turbulence_Statistics_Pipeline_Specification — authoritative design
- @ref 56_Field_Identity_and_Layout_Catalog — typed field identity and layout
- @ref 29_Maintenance_Backlog — measured restart fidelity floor
- @ref 09_Monitor_Reference — monitor contract, updated at Stage 7
- @ref 10_Post_Processing_Reference — postprocessing contract, updated at Stage 7
