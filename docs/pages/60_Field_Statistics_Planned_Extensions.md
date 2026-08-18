@page 60_Field_Statistics_Planned_Extensions Field Statistics Planned Extensions

@anchor _Field_Statistics_Planned_Extensions

@ref 58_Field_Statistics describes what the statistics pipeline does today:
pointwise weighted centered first and second moments, accumulated per named
window, carried across a restart, and derived into Reynolds stresses, RMS,
turbulent kinetic energy, and fluxes.

This page records what it deliberately does not do, why each extension was left
out, and what it depends on. Each is opt-in and builds on the target and product
abstractions that already exist, so none of them requires reworking what is
there. Configuration written against the current contract keeps working through
all of them.

@tableofcontents

@section p60_order_sec 1. Dependency Order

| Extension | Depends on | Unlocks |
| --- | --- | --- |
| @ref p60_view_sec | nothing | nodal statistics output, accumulator field logging |
| @ref p60_targets_sec | nothing | bins, profiles, regions, surfaces, probes |
| @ref p60_products_sec | @ref p60_targets_sec for face targets | face-flux and single-component products |
| @ref p60_ess_sec | nothing | uncertainty reporting under unequal weights |
| @ref p60_history_sec | @ref p60_targets_sec | histories, PDFs, higher moments |
| @ref p60_reducer_sec | @ref p60_targets_sec | one reducer behind flow observables |
| @ref p60_dimensional_sec | nothing | dimensional derived statistics |

The view extension is listed first because it is the only one that two separate
surfaces are already waiting on, and because it is self-contained.

@section p60_view_sec 2. A View Bindable To An Explicit Vector Pair

**What is missing.** `UpdateLocalGhosts` resolves its vectors through a field view
built on compile-time `UserCtx` offsets recorded in the catalog, and the nodal
averaging kernel maps fixed field names onto fixed members. Per-window
accumulators are counted from configuration and therefore have no compile-time
offset, so neither routine accepts one.

**Current workaround.** Derived results are staged through two catalogued scratch
fields: a derived quantity is written there, reached by catalogued name, converted
to nodal values, and copied out before the next quantity reuses the buffer. This
works with the offset constraint rather than against it, and follows the precedent
of the corner-staging buffers, which are also workspace rather than simulation
state.

**What the extension is.** A `FieldView` that can be bound to an explicitly
supplied global and local vector pair, carrying its own degree of freedom, layout,
and synchronization class rather than reading them from a compile-time table.
After that the whole ghost path applies to an accumulator unchanged.

**Why it is worth doing.** Two surfaces depend on it at once. Nodal statistics
output stops needing a staging round trip, and accumulator vectors become
addressable by the field-anatomy and min/max loggers, which today can only report
window-level scalars. A statistics-specific logging path would solve half of that
and add a second mechanism doing what the first already does.

@section p60_targets_sec 3. Spatial Targets Beyond Pointwise

`SpatialTargetPlan` exists as an abstraction with only the pointwise identity
mapping implemented, so these extend it rather than retrofit it.

Planned target kinds:

- **Bins** over a coordinate direction or a derived coordinate.
- **Profile layers**, the wall-normal case that channel and boundary-layer work
  needs.
- **Named regions**, an axis-aligned or geometric subset accumulated as one
  reduced quantity.
- **Grid and immersed surfaces**, accumulating over a boundary rather than a
  volume.
- **Point probes**, a small named set of locations sampled at every accepted
  state.

Each is a different mapping from field storage indices to accumulator slots, which
is exactly what the plan abstraction is for. The reduction they need is the
weighted parallel-combination form that already exists and is tested but has no
production caller yet.

**The YAML is already reserved for this compatibly.** A `target:` block with
`kind` and `mask` arrives as **optional** keys defaulting to current behaviour, so
no configuration written today breaks. They are not exposed now because a
single-valued key reserves syntax without earning it.

**The identity hash is already prepared.** The mask occupies its own hash group
even though exactly one mask is resolved today, so adding a mask key extends that
group rather than renumbering the groups after it and invalidating every saved
window.

@section p60_products_sec 4. Further Products

**Explicit `Ucont` face-flux products.** `Ucont` is component-staggered, so pairs
involving it are rejected on layout compatibility today. Face-flux statistics need
a face target and a product definition that respects which component lives on
which face, not a relaxation of the compatibility check.

**Single-component products.** Requesting one component of a vector's self-product
rather than all six, for a run that needs only a shear stress and cannot afford
the rest.

**Vector-vector cross products.** These need a full nine-component carrier, which
nothing currently allocates. The symmetric six-component form exists precisely
because a self-product is symmetric; a cross product between two different vectors
is not.

@section p60_ess_sec 5. Effective Sample Size Reporting

Under physical-time weighting the weights are unequal, so a raw sample count
overstates how much independent information a window holds. The Kish effective
sample size \f$W^2/W_2\f$ is the right quantity to report instead, and the kernel
that computes it already exists and is tested.

What is missing is plumbing rather than mathematics: a window-level sum of squared
weights, which means a new checkpoint metadata field and a restore path for it.
Per-point squared weight is already accumulated and checkpointed; only the
window-level scalar is absent.

@section p60_history_sec 6. Histories, PDFs, and Higher Moments

Reduced time histories, probability density functions, and moments beyond the
second **cannot be recovered from centered state after the fact**. They require
their own accumulation, and are therefore additions rather than derivations.

They are constrained by a boundary that must hold: whatever retains history must
stay bounded, for reduced probes or regions only. Turning the statistics manager
into a second full-field output system is not on this list at any point. A run that
needs the states themselves selects ordinary instantaneous field output.

@section p60_reducer_sec 7. One Reducer Behind Flow Observables

`ComputeCurrentFlowObservables` computes its own domain reductions for the runtime
convergence monitor. The spatial reducer these extensions introduce would be able
to serve it, removing a second implementation of the same traversal.

This is permitted **only after tests demonstrate identical results**, and the
replacement must correctly follow field layout, periodicity, masks, blocks, and
MPI ownership. It is listed last among the functional extensions because the
existing implementation is correct and the gain is structural.

Function-level logging, logging allow-lists, PETSc monitors, and performance
profiling are not subsumed into scientific statistics at any point.

@section p60_dimensional_sec 8. Dimensional Derived Statistics

Derived statistics are written in nondimensional form even when
`global_operations.dimensionalize` is set. A Reynolds stress scales as velocity
squared, and a co-moment as the product of two different scales; the existing
per-field scaling table expresses neither.

The extension is a scaling rule per derived quantity rather than per field, so a
product knows the scales of both its factors. Applying the velocity scale as things
stand would be wrong rather than merely incomplete, which is why the current
behaviour is to leave derived output nondimensional and say so.

@section p60_notplanned_sec 9. Not Planned

These are recorded so they are not proposed again as gaps.

- **A legacy `su0`/`su1`/`su2`/`sp` importer.** The old averaging state was removed
  rather than retained as a second workflow, and no reader, writer, or converter
  for it is planned.
- **A dual-write period or compatibility workflow** between the old averaging
  surface and field statistics.
- **Relabelling modelled turbulence quantities.** `Nu_t` and `CS` can be averaged
  like any other field, but modelled `k` is never presented as turbulent kinetic
  energy derived from resolved fluctuations.
- **Full-field time histories**, as above.

@section p60_related_sec 10. Related Pages

- @subpage 58_Field_Statistics — what the pipeline does today
- @subpage 56_Field_Identity_and_Layout_Catalog — the catalog these extensions build on
- @subpage 57_Future_Architecture_Specifications — status and sequencing of proposed architecture
- @subpage 29_Maintenance_Backlog — open validation items
