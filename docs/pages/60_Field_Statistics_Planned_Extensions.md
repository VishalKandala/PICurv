@page 60_Field_Statistics_Planned_Extensions Field Statistics Planned Extensions

@anchor _Field_Statistics_Planned_Extensions

@ref 58_Field_Statistics describes what the statistics pipeline does today:
pointwise weighted centered first and second moments, accumulated per named
window, carried across a restart, and derived into Reynolds stresses, RMS,
turbulent kinetic energy, and fluxes.

This page records what it deliberately does not do, why each extension was left
out, and what it depends on. Each is opt-in and builds on the abstractions that
already exist, so none of them requires reworking what is there. Configuration
written against the current contract keeps working through all of them.

@tableofcontents

@section p60_principle_sec 1. What Must Be Online, And What Must Not

One rule decides where each extension belongs, and it is worth stating before
the list because it removes more work than it adds.

**Averaging commutes with linear operations and not with anything else.**

A spatial reduction over accumulated pointwise state is linear, so it is exact
after the fact. Reducing per-point accumulators over a region reproduces the
region-wide centered moment about the region mean, including the spatial
variance of the local means, through the same weighted combination
`PicurvMomentStateMerge` already implements:

\f[
M_{2,AB} = M_{2,A} + M_{2,B} + \frac{W_A W_B}{W_A + W_B}\,(\mu_B - \mu_A)^2 .
\f]

That identity extends to higher orders, so whatever order is accumulated
pointwise can be reduced to any spatial subset at that order or below without
approximation. **Bins, profiles, regions, and clipping are therefore
post-processing operations, not accumulator kinds.** They need no online
machinery, no additional checkpoint state, and no monitor-side configuration.

What genuinely must happen online is everything the identity does not reach:

- **The moment order.** Order cannot be raised after the fact; the second moment
  is what is accumulated today.
- **Which products.** A cross-moment is not derivable from its marginals, so a
  pair not requested at solve time is unavailable afterwards.
- **Statistics of nonlinear or interpolated quantities.** For a probe or surface
  value built from a stencil, \f$P = \sum_i w_i u_i\f$, the mean interpolates
  but the variance does not: \f$\mathrm{Var}(P)\f$ needs the covariance
  *between* neighbouring points, which pointwise state does not hold. The same
  applies to the statistics of any derived nonlinear field.
- **Conditional sampling.** Unconditional weights cannot be split afterwards.

Memory is the one remaining argument for accumulating a reduced quantity online
rather than reducing later, and it is stratified: negligible for moments,
decisive for per-point histograms, prohibitive for retained histories.

A physical caveat applies to spatial reduction however it is computed: merging
across points is only *meaningful* where those points are statistically
equivalent. Over an inhomogeneous region the result is a well-defined total
variance that mixes in the spatial variation of the mean, which is rarely what
was wanted. Homogeneous-direction averaging is the case this serves.

@section p60_order_sec 2. Dependency Order

| Extension | Depends on | Unlocks |
| --- | --- | --- |
| @ref p60_reduction_sec | nothing | profiles, regions, bins, clipping |
| @ref p60_view_sec | nothing | nodal statistics output, accumulator field logging |
| @ref p60_ess_sec | nothing | uncertainty reporting under unequal weights |
| @ref p60_order_extension_sec | nothing | skewness, flatness, and their spatial reductions |
| @ref p60_products_sec | online face targets, for the face case | face-flux and single-component products |
| @ref p60_targets_sec | @ref p60_view_sec for surface output | surfaces, probes, conditional sampling |
| @ref p60_reducer_sec | @ref p60_reduction_sec | one reducer behind flow observables |
| @ref p60_dimensional_sec | nothing | dimensional derived statistics |
| @ref p60_spectra_partial_sec | nothing | spectra for channels, ducts, and boundary layers |
| @ref p60_spectra_temporal_sec | bounded probe history | frequency spectra where nothing is homogeneous |
| @ref p60_spectra_online_sec | @ref p60_spectra_temporal_sec for the transpose | spatial spectra without a snapshot series |

Spatial reduction is listed first because it is the largest capability for the
least new code: the kernel it needs is already written and tested, and it has no
online counterpart to design.

@section p60_view_sec 3. A View Bindable To An Explicit Vector Pair

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

@section p60_reduction_sec 4. Spatial Reduction In Post-Processing

Profiles, regions, bins, and clipping are derived from the pointwise state a
window already holds, by the reduction described in @ref p60_principle_sec. This
is where the majority of the remaining scientific capability lives, and it needs
no new online state.

What it requires:

- **A region selector in the post recipe.** Coordinate ranges, index ranges, a
  homogeneous-direction collapse such as averaging over `i` and `k` at fixed `j`,
  or a named region. This sits in `post.yml` because the choice is made after
  the run, not committed at solve time.
- **A reduction driver** that walks the selected points and combines their
  accumulators through the existing weighted merge, in a pairwise or tree order
  rather than sequentially, so the combination stays well conditioned across a
  large point count. `PicurvWindowSpatialMean` is the shape to follow; it already
  performs the domain-wide case of exactly this reduction, and it already
  excludes points the mask never sampled.
- **An output form for reduced results.** A profile is a one-dimensional table
  rather than a field, so it belongs in a CSV alongside the convergence history
  rather than in the `.vts`.

The practical gain beyond cost is that regions are chosen *after* the solve.
Because the pointwise state is checkpointed, a wake window or a wall-normal
profile that nobody anticipated can be extracted from a completed run without
re-solving, and revised as often as the analysis demands.

@section p60_targets_sec 5. Online Targets For The Non-Commuting Cases

A target must be resolved online only where @ref p60_principle_sec says the
reduction cannot follow:

- **Immersed and grid surfaces**, where the quantity is interpolated onto a
  surface that cuts through cells. The surface mean can be interpolated from
  pointwise means, but its variance cannot, because that needs the covariance
  between the stencil points.
- **Point probes at off-grid locations**, for the same reason.
- **Conditional sampling**, where a state is accumulated only when a condition
  holds, since unconditional weights cannot be split afterwards.

`SpatialTargetPlan` exists as an abstraction with only the pointwise identity
mapping implemented, so these extend it rather than retrofit it, and the
identity hash already reserves its own group for the mask so a later mask key
extends that group instead of renumbering the groups after it and invalidating
every saved window.

The user-facing surface for these is a `target:` block with `kind` and `mask` as
**optional** keys defaulting to current behaviour, so no configuration written
today breaks. They are not exposed now because a single-valued key reserves
syntax without earning it.

@section p60_order_extension_sec 6. Raising The Accumulated Moment Order

The second moment is what is accumulated today, which fixes the ceiling on what
any later reduction can produce: order cannot be raised after the fact. Raising
it to the third and fourth moments would make skewness and flatness available,
and — because the merge identity extends to higher orders — would make them
available for every spatial reduction at the same time, not only for the whole
domain.

What this needs: additional per-point state, the higher-order terms of the
weighted combination, the corresponding checkpoint payloads, and moment entries
in the window definition so the requested order is part of the identity hash.
The per-point cost is one additional field per order per accumulated quantity,
which is modest next to what a per-point histogram would cost.

This is the extension that most enlarges what post-processing can do, because it
raises the ceiling rather than adding a mechanism.

@section p60_pdf_sec 7. Distributions, Deferred

Per-point probability density functions and histograms are deferred rather than
designed. They would reduce spatially without difficulty — histograms over shared
bin edges are additive, so a region histogram is the sum of its points' — but the
per-point cost is one field per bin per quantity, which is one to two orders of
magnitude above the moment path. That cost, not the mathematics, is what defers
them.

Should they arrive, bin edges must be fixed and global rather than per-point
adaptive, since additivity is what makes the spatial reduction exact.

Retained time histories remain outside this pipeline at any point. Whatever
retains history must stay bounded, for reduced probes or regions only; turning
the statistics manager into a second full-field output system is not on this
list. A run that needs the states themselves selects ordinary instantaneous
field output.

@section p60_products_sec 8. Further Products

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

@section p60_ess_sec 9. Effective Sample Size Reporting

Under physical-time weighting the weights are unequal, so a raw sample count
overstates how much independent information a window holds. The Kish effective
sample size \f$W^2/W_2\f$ is the right quantity to report instead, and the kernel
that computes it already exists and is tested.

What is missing is plumbing rather than mathematics: a window-level sum of squared
weights, which means a new checkpoint metadata field and a restore path for it.
Per-point squared weight is already accumulated and checkpointed; only the
window-level scalar is absent.

@section p60_reducer_sec 10. One Reducer Behind Flow Observables

`ComputeCurrentFlowObservables` computes its own domain reductions for the runtime
convergence monitor. The reduction driver in @ref p60_reduction_sec performs the
same traversal, so it could serve that monitor too and remove one of the two
implementations.

This is permitted **only after tests demonstrate identical results**, and the
replacement must correctly follow field layout, periodicity, masks, blocks, and
MPI ownership. It is listed last among the functional extensions because the
existing implementation is correct and the gain is structural.

Function-level logging, logging allow-lists, PETSc monitors, and performance
profiling are not subsumed into scientific statistics at any point.

@section p60_dimensional_sec 11. Dimensional Derived Statistics

Derived statistics are written in nondimensional form even when
`global_operations.dimensionalize` is set. A Reynolds stress scales as velocity
squared, and a co-moment as the product of two different scales; the existing
per-field scaling table expresses neither.

The extension is a scaling rule per derived quantity rather than per field, so a
product knows the scales of both its factors. Applying the velocity scale as things
stand would be wrong rather than merely incomplete, which is why the current
behaviour is to leave derived output nondimensional and say so.

@section p60_spectra_partial_sec 12. Offline Line And Plane Spectra

**What exists today.** `post.yml -> spectra` implements one task, `shell_spectrum`,
which requires all three directions to be periodic and homogeneous. See
@ref p10_spectra_sec.

**What is missing.** Everything between "homogeneous in three directions" and
"homogeneous in none": a channel homogeneous in two directions and bounded in the
third, a straight duct homogeneous in one, a spatially developing boundary layer
homogeneous in the spanwise direction alone. This is the most common wall-bounded
turbulence configuration there is, and it is presently unserved — a channel has a
perfectly good spatial spectrum that the pipeline cannot produce.

Note this is *not* what @ref p60_spectra_temporal_sec covers. That section is for
geometries with no homogeneous direction at all, where only a frequency spectrum is
defined. A channel has homogeneous directions, and what is wanted there is a
wavenumber spectrum resolving scales, not a frequency spectrum at a probe.

**Why it is first among the three spectral extensions.** It needs no new online state,
no retained history, no FFT in the solver, and no MPI transpose. It reuses the offline
stage exactly as it stands: read a committed checkpoint, transform, bin, write. Of the
three, it is the only one that is purely additive to machinery already running.

@subsection p60_spectra_partial_shape_sub The reduction is the choice, not the topology

The obvious design — infer the task from how many directions are periodic — does not
survive contact with what people actually want. Two homogeneous directions admit at
least four distinct deliverables, and they answer different questions:

- \f$E(k_x)\f$ at each wall-normal station,
- \f$E(k_z)\f$ at each wall-normal station,
- the full two-dimensional \f$E(k_x, k_z)\f$ at one chosen station,
- premultiplied \f$k E(k)\f$, which is what log-layer scaling is read from.

The most common channel deliverable is not "the plane spectrum" at all: it is a plane
transform *reduced two different ways*. So the map from topology to task is
one-to-many, and the topology cannot pick for the user.

**Periodicity is also necessary rather than sufficient.** The pairwise consistency
check in `validate_and_prepare_boundary_conditions` makes "this axis is periodic" a
reliable, already-validated property, but homogeneity is the real requirement. An
immersed body inside a periodic box breaks it while every boundary condition still
reads `PERIODIC`, and the `Nvert` mask that reveals it is a runtime property that moves
with the body, so no configuration-time inference can see it.

**What may be inferred.** The *direction* may, where it is unambiguous: a
`line_spectrum` on a case with exactly one periodic axis has only one sensible
`direction`, and defaulting it removes real tedium without removing a scientific
decision. The *task* may not, for the reasons above.

That distinction matters beyond convenience. The value of the precondition gate is that
it refuses: ask for something undefined and the reason is named before any compute is
spent. Inferring the task inverts that — when the inference is right a configuration
line is saved, and when it is wrong there is no error at all, only a plausible curve.
A pipeline built on @ref p60_principle_sec should not trade a loud failure for a silent
one to save a line of YAML.

What is worth doing instead is making the refusal constructive: a case whose spanwise
axis is periodic could be told that `line_spectrum` along that axis would be valid,
rather than only being told what is wrong.

@subsection p60_spectra_partial_cost_sub What it costs

The transform is the easy part. What genuinely changes:

- **Output stops being one-dimensional.** `shell_spectrum` yields \f$E(k)\f$ per
  step. A line spectrum yields \f$E(k_1)\f$ per retained index, a plane spectrum
  \f$E(k_x,k_z)\f$ per station. The long-format CSV generalizes with an index column;
  the plot does not, because "a spectrum at every wall-normal station" is a family or a
  contour rather than a curve. `--plot-spectrum` needs a station selector.
- **Masking becomes per-line.** `shell_spectrum` sidesteps immersed boundaries by
  refusing them outright. A line spectrum must test each transform line for blanked
  points and drop or refuse that line individually.
- **Serial memory becomes the ceiling.** `generators/spectra.gen` reads whole fields in
  Python. A 64-cubed box is under seven megabytes; a 512-cubed channel is a few
  gigabytes per snapshot. For channel-scale cases this, rather than sample count, is
  what eventually forces @ref p60_spectra_online_sec.

**Test case.** `examples/flat_channel`, which is homogeneous in the streamwise and
spanwise directions. Not `examples/bent_channel`: a 90-degree bend develops
streamwise and is bounded on all four sides, so it has no homogeneous direction and
belongs to @ref p60_spectra_temporal_sec.

@section p60_spectra_temporal_sec 13. Temporal Spectra From Bounded Probe Histories

**What exists today.** `post.yml -> spectra` measures spatial spectra offline from the
instantaneous velocity in each committed checkpoint, documented in
@ref p10_spectra_sec. The one implemented task, `shell_spectrum`, requires a triply
periodic uniform Cartesian box.

**What is missing, and why it is not a gap in the above.** A spatial spectrum needs a
direction that is uniformly spaced, periodic, and statistically homogeneous. A curved
duct, a wall-bounded channel section, or any immersed geometry supplies none of the
three: the streamwise direction develops, and the cross-stream directions are bounded.
For these the meaningful quantity is a **frequency** spectrum \f$E(f)\f$ at a point,
which needs no homogeneity at all and is what resolves shedding and Strouhal content.

A frequency spectrum needs the signal itself, retained in time. That is the one thing
@ref p60_pdf_sec rules out at full-field scale and explicitly leaves open at reduced
scale: *whatever retains history must stay bounded, for reduced probes or regions
only*. A handful of probes sampled every step is kilobytes; the same information as
full-field snapshots is terabytes.

**Why the snapshot path cannot substitute.** A temporal transform of the existing
checkpoint series is limited to
\f$f_{\max} = 1/(2\,\Delta t\,\texttt{data\_output\_frequency})\f$. Reaching a
useful band means writing full three-dimensional fields every step or two, which is
precisely the cost this extension avoids. The snapshot series can reach a
low-frequency peak; it cannot reach an inertial range.

What it requires:

- **A probe set in the monitor configuration**, since retaining history is a solve-time
  decision rather than a post-processing one.
- **Bounded per-probe history**, appended as the run proceeds, with the sampling cadence
  fixed at configuration time so the transform has a uniform time base.
- **Interpolation for off-grid probes.** Note that the caveat in
  @ref p60_targets_sec does not apply here: that section concerns the *variance* of an
  interpolated quantity, which needs the covariance between stencil points. A retained
  time series holds the actual interpolated values, so its spectrum is exact.
- **A `temporal_spectrum` task** in the spectra recipe, with Welch averaging over
  segments and a window function, since a finite record is not periodic in time.

This is the extension that unlocks spectra for every geometry the spatial tasks refuse,
and it is the smaller of the two remaining spectral pieces.

@section p60_spectra_online_sec 14. Online Spatial Spectral Accumulation

**What is missing.** Spatial spectra are measured offline from written snapshots, so
the number of independent samples is bounded by how many full fields a run can afford
to write, and the snapshot cadence has to be chosen before the run.

**What the extension is.** Accumulate \f$\langle |\hat{u}(k)|^2 \rangle\f$ online:
transform along the homogeneous direction each accepted state and add the modulus
squared into a per-wavenumber accumulator.

**Why it must be online rather than derived later.** The transform is linear, so a
time-averaged \f$\hat{u}\f$ is recoverable from a time-averaged field. The modulus
squared is not: \f$\langle |\hat{u}|^2 \rangle \neq |\langle \hat{u} \rangle|^2\f$.
This is the same rule @ref p60_principle_sec states, applied in Fourier space, and it
is also why no reduction of the existing pointwise accumulators can yield a spectrum —
a spectrum is the transform of a two-point correlation, and pointwise state holds no
covariance between distinct points.

**Why it is nonetheless second in priority.** It is an efficiency change, not a
capability change: it produces spectra the snapshot path already produces, for runs
where writing the snapshots is the binding constraint. @ref p60_spectra_temporal_sec
produces spectra that are otherwise unavailable at all.

What it requires:

- **A spectral target kind and product**, slotting into the existing window
  machinery. Scheduling, weighting, restart merging, and checkpointing all apply
  unchanged; the accumulator is a new product, not a new lifecycle.
- **An FFT in the solver.** PETSc's `MatCreateFFT` needs an FFTW-enabled build, which
  the current configuration does not have. A single-direction transform is cheap
  enough that gathering the pencil and transforming directly is a viable alternative.
- **A transpose or line gather** when the homogeneous direction is distributed across
  ranks.

Storage is negligible — wavenumber by wall-normal index by component — which is orders
of magnitude below the per-point histograms @ref p60_pdf_sec defers.

@section p60_notplanned_sec 15. Not Planned

These are recorded so they are not proposed again as gaps.

- **A legacy `su0`/`su1`/`su2`/`sp` importer.** The old averaging state was removed
  rather than retained as a second workflow, and no reader, writer, or converter
  for it is planned.
- **A dual-write period or compatibility workflow** between the old averaging
  surface and field statistics.
- **Relabelling modelled turbulence quantities.** `Nu_t` and `CS` can be averaged
  like any other field, but modelled `k` is never presented as turbulent kinetic
  energy derived from resolved fluctuations.
- **Full-field time histories**, as above. Per-point distributions are a
  different matter: they are deferred on cost rather than ruled out, and
  @ref p60_pdf_sec records the terms under which they could arrive.

@section p60_related_sec 16. Related Pages

- @subpage 58_Field_Statistics — what the pipeline does today
- @subpage 10_Post_Processing_Reference — the spectra recipe these two extend
- @subpage 56_Field_Identity_and_Layout_Catalog — the catalog these extensions build on
- @subpage 57_Future_Architecture_Specifications — status and sequencing of proposed architecture
- @subpage 29_Maintenance_Backlog — open validation items
