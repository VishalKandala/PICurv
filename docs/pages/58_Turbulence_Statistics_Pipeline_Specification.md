@page 58_Turbulence_Statistics_Pipeline_Specification Field Statistics Pipeline Specification

@anchor _Turbulence_Statistics_Pipeline_Specification

This is the authoritative scientific and architectural specification for the
planned Eulerian field-statistics pipeline. It records settled decisions while
keeping unfinished features out of the active user contract.

Implementation status: Phase 1 checkpoint infrastructure is implemented on the
working branch. The typed field catalogs drive persistence, the legacy
averaging system is retired, existing physical-solution convergence settings
live at `monitor.yml -> solution_monitoring.convergence`, and committed bundles
are shared by solve, restart, continuation, and postprocessing. Scientific
statistics YAML, accumulators, and derived-statistics postprocessing are not
active yet; they remain Phase 2 work after Phase 1 review and merge.

@tableofcontents

@section p58_goals_sec 1. Goals and Boundaries

The completed system must provide:

- independent and overlapping named statistics windows;
- sample-weighted and physical-time-weighted first and second moments;
- numerically stable online accumulation without retaining every snapshot;
- field-catalog-aware handling of cell, node, face, and component-staggered
  layouts;
- physical-domain, periodicity, mask, MPI, multiblock, and moving-geometry
  correctness;
- restart-safe state with no duplicate or missing samples;
- a pointwise first implementation that can later support common spatial
  targets through the same reducer abstraction; and
- one postprocessor path for normalization and derived resolved-turbulence
  quantities.

The online solver collects only registered primitive fields and the centered
moment state required by the selected products. Reynolds stress, RMS, resolved
TKE, turbulent fluxes, dimensional normalization, profiles, and visualization
outputs are postprocessor responsibilities.

LES and RANS fields remain model quantities. They can later be sampled
explicitly, but modeled `k` is never relabeled as TKE derived from resolved
velocity fluctuations. Particle ensemble tasks such as MSD remain in the
particle-statistics subsystem.

Full-field time histories are explicitly outside this accumulator. Users who
need them select ordinary instantaneous field output. A later extension may
retain bounded histories for reduced probes or regions, but it must not turn
the statistics manager into a second full-field output system.

@section p58_decisions_sec 2. Settled Design Decisions

| Question | Decision |
| --- | --- |
| Existing physical monitoring YAML | `monitor.yml -> solution_monitoring.convergence` |
| Future scientific statistics YAML | `monitor.yml -> field_statistics` |
| Future derived-statistics recipes | `post.yml -> field_statistics` |
| Generated C ingress | the master `control` file only |
| Online first/second-moment representation | weighted centered state |
| First spatial implementation | pointwise physical locations |
| State ownership | independent PETSc state for every named window |
| Default velocity source | cell-centered Cartesian `Ucat` |
| `Ucont` statistics | not collected by default; only future explicit face/flux products |
| Derived Reynolds statistics | postprocessor only |
| Runtime monitoring | preserve existing writer and every-completed-step behavior; do not checkpoint persistent logs |
| Scientific compatibility identity | C-computed SHA-256 of the resolved immutable window definition |
| Legacy statistics | removed; no compatibility reader, writer, or importer required |
| Function IDs | separate, benchmark-gated future specification |

The word `observation` is not a user-facing container. Runtime solution
monitoring and scientific field statistics are related consumers but have
separate, direct vocabulary.

`control` is generated afresh by `picurv` whenever a solve is staged and is not
a user configuration surface. No `observations.run`, public plan schema, or
external control-file checksum is required. The YAML inputs and run manifest
remain the user-facing provenance. C is authoritative for resolving fields and
computing scientific compatibility identity.

@section p58_active_monitor_sec 3. Active Monitoring Contract

The currently active configuration is:

```yaml
solution_monitoring:
  convergence:
    enabled: true
    mode: statistical_steady
    statistical_steady:
      window_steps: 100
```

It maps directly to the existing `-solution_convergence_*` options in the
generated master control. The current `logs/solution_convergence.log` writer,
columns, mode behavior, and observable implementation remain in place.
Monitoring runs after every completed timestep, so no cadence option is
exposed.

`ComputeCurrentFlowObservables` is retained for now. A future shared reducer may
replace or sit behind it only after tests demonstrate identical results and the
new implementation correctly follows field layout, periodicity, masks, blocks,
and MPI ownership. Function logging, logging allowlists, PETSc monitors, and
performance profiling are not subsumed into scientific statistics.

@section p58_future_yaml_sec 4. Future Scientific User Contract

The intended location, exposed only once the feature works end to end, is:

```yaml
field_statistics:
  windows:
    - name: production
      start_time: 50.0
      end_time: 250.0
      include_initial: false
      schedule:
        every_steps: 1
      weighting: physical_time
      target:
        kind: pointwise
        mask: fluid
      fields:
        - field: Ucat
          moments: [first, second]
        - field: P
          moments: [first]
```

The exact Phase 2 spelling will be reviewed before implementation, but these
semantics are settled:

- a window name is stable and unique within its lineage;
- `start_time` is required and `end_time` may be open;
- one positive step- or physical-time schedule is selected;
- weighting is `sample` or `physical_time`;
- initial-state inclusion is explicit and defaults to false;
- fields resolve once through the typed field catalog;
- requested products must be scientifically compatible with field locations;
- changing products, layouts, masks, weighting, or schedule starts a new named
  window; and
- full-field pointwise history requests are rejected in favor of normal output.

No inactive keys are accepted early merely to reserve syntax. Python validation,
control serialization, C resolution, online state, checkpoint continuation,
and postprocessing must arrive as one usable contract.

@section p58_lifecycle_sec 5. Runtime Lifecycle

The high-level solver lifecycle remains unchanged:

```text
initialize PETSc/MPI and the simulation context
  -> create the existing DMs and solver-owned vectors
  -> initialize boundaries, models, and optional particles
  -> resolve requested statistics fields from the catalog
  -> create or restore only the requested window state
  -> enter the existing time loop
```

The statistics system does not create or own source flow fields. It obtains
non-owning `FieldView` objects after the existing vectors exist and owns only
its window schedules, accumulator vectors, and persistence metadata.

After a timestep has completed, the ordering is:

```text
final Eulerian and model state
  -> final particle motion and particle-to-grid state, when enabled
  -> one completed-state event
  -> update each due statistics window exactly once
  -> update physical-solution monitoring
  -> existing log writers
  -> solver history
  -> checkpoint/output when due
```

This ensures every sampled source belongs to the same physical state. A saved
last-event key prevents duplicate accumulation across ordinary continuation,
restart-from, graceful shutdown, and checkpoint replay.

Each window owns separate vectors. The implementation does not reuse one vector
and fill off-window timesteps with zeros. A pending or off-schedule window makes
no accumulator update. Storage may be allocated during setup or lazily at first
activation, provided restart and memory behavior are identical and testable.

@section p58_window_sec 6. Window and Weighting Semantics

A window is `pending`, `active`, or `complete`. It records its requested and
effective bounds, sample count, total weight, represented physical time,
schedule position, last event key, resolved definition hash, and restart
lineage.

- Before the requested start, it receives no samples.
- The first eligible due state records the effective start.
- An active off-schedule state changes no scientific state.
- A due completed state is accepted exactly once.
- Overlapping windows read a source state once where practical but update
  independent accumulator state.
- A bounded window stops at its configured end; physical-time weighting clips
  the final represented interval to that bound.

Sample weighting assigns equal weight to accepted states. Physical-time
weighting uses the represented interval between accepted states, including
variable `dt` and a clipped final interval. `sample_count` and represented time
are never inferred from the global timestep number.

The precise endpoint quadrature rule and closing-event behavior must be fixed
in the Phase 2 implementation spec and validated with analytic variable-step
signals before the YAML is released.

@section p58_centered_sec 7. Centered Accumulator Contract

For a scalar sample `x` with positive weight `w`, a point stores sample count,
total weight `W`, optional squared-weight sum `W2`, mean `mu`, and centered
second-moment sum `M2`. The stable weighted update is:

\f[
W' = W+w,\qquad \delta=x-\mu,
\f]

\f[
\mu'=\mu+\frac{w}{W'}\delta,\qquad
M_2'=M_2+w\,\delta(x-\mu').
\f]

Compatible covariance uses the corresponding stable co-moment update. A
three-component `Ucat` self-product stores the six symmetric centered
components in fixed catalog metadata order. Every moving-mask point keeps its
own valid count and weight.

This representation retains all information needed for first and second
moments. Raw sums can be reconstructed, up to normal floating-point effects:

\f[
S_x=W\mu,\qquad Q_{xx}=M_2+W\mu^2.
\f]

Therefore storing raw sums as a second authoritative copy would add memory,
checkpoint volume, synchronization, and consistency risk without adding
first/second-moment information. Centered state is also mergeable using the
standard weighted parallel-combination formula, which is required for offline
window merging and later spatial reductions even though pointwise PETSc vector
updates themselves are local.

What centered state does not retain is temporal ordering, a PDF, arbitrary
higher moments, or the ability to redefine a smaller time subwindow after the
run. Those require explicitly requested future products or saved instantaneous
snapshots. A saved snapshot can be converted into a fluctuating field by
subtracting the stored window mean, but an unsaved instantaneous field cannot
be reconstructed.

@section p58_layout_sec 8. Field, Grid, Boundary, and Mask Contract

Every request resolves through @ref FieldDescriptor and existing range/layout
helpers. Iteration excludes both:

1. PETSc/MPI halo storage; and
2. solver-layout boundary, dummy, unused, or duplicate periodic indices.

This is not a generic "skip ghost cells" rule. Cell-centered fields, node
fields, I/J/K-face fields, and packed component-staggered fields have different
physical index domains. Periodic boundary selection changes DM topology and
representative ranges and must use the already-established periodic logic.

The default `fluid` mask accepts a physical cell only when `Nvert < 0.1`.
`all_physical` still excludes nonphysical storage. Existing mask, geometry,
surface, and range helpers are reused rather than reimplemented. With moving
immersed geometry, every point retains independent valid weight so a cell that
was fluid for only part of a window is not normalized as if it were always
available.

`Ucat` is the default velocity primitive because its Cartesian components are
collocated at cell centers and directly support the six resolved covariance
components. `Ucont` is not accumulated by default. It is component-staggered,
so a generic three-component covariance is scientifically invalid. A future
explicit face-flux or single-component request may be added if a concrete use
case needs it and the catalog advertises the necessary capability.

@section p58_spatial_sec 9. Spatial Targets and Future Bins

The first implementation is pointwise because it preserves the greatest
spatial information: one centered state at every valid source location.
Postprocessing can merge those states with volume or area measures to produce
many spatial summaries without rerunning the solve.

The backend must nevertheless use a general `SpatialTargetPlan` abstraction.
Its initial pointwise mapping is the identity mapping from a physical location
to its accumulator. Later implementations may map locations into:

- one global domain bin;
- one bin per computational or physical-coordinate profile layer;
- a named volume region;
- a named grid or immersed surface; or
- a point probe represented as a one-location bin.

A bin is simply a named group of physical locations plus its selection and
spatial weighting rule. Front-end vocabulary can provide convenient profiles,
probes, regions, and surfaces while reusing that backend. Existing Jacobian,
centroid, face-area, immersed-surface, and MPI/block reduction code must be used
where applicable.

Pointwise postprocessing can form a spatial distribution of pointwise temporal
statistics. It cannot reconstruct the time history of a spatially reduced
observable. If a later user needs a PDF or time history of a plane mean, probe,
or region integral, that reduced product must be explicitly collected online.

@section p58_hash_sec 10. Window Identity and Allowed Changes

At startup, C resolves canonical field descriptors, layouts, products, mask,
weighting, schedule, target semantics, and other scientifically immutable
properties. It serializes that small resolved definition canonically and
computes its SHA-256 once. It does not hash PETSc vector contents or recompute
the hash each timestep, so the cost is negligible.

The hash is internal checkpoint metadata, not a user option and not a checksum
for the generated control file. It answers one question on restore: "Is this
saved accumulator scientifically the same window state requested now?"

The requested `end_time` is excluded from the immutable definition hash so an
active window may be extended. Extension is accepted separately only when:

- the new end moves forward;
- it does not exceed the configured simulation horizon; and
- continuation occurs without an unsampled gap after the former end.

Shortening an accumulated window or changing its start, schedule, weighting,
mask, fields, products, layout semantics, or target definition requires a new
named window.

@section p58_restart_sec 11. Continue and Restart-From Semantics

`--continue` resumes in the same run directory and requires the same physical
case and compatible window definition. Solver parameters, monitoring,
postprocessing, resource settings, and continuation bookkeeping may change.
An allowed window-end extension is checked separately as described above.

`--restart-from` creates a new run directory and may branch the flow case while
retaining the required geometry/layout compatibility of the restart state.
Statistics continue only when the physical case identity and resolved window
definition still match and continuation was explicitly requested. If physics
changed, or a different statistic is desired, the flow may restart but the
statistics must reset under a new window name.

Both paths load flow, optional model/particle state, window accumulators,
schedules, weights, and last-event identity from one matching committed
checkpoint. Missing required statistics state is fatal for a requested
continuation; it is never silently zeroed. A different MPI rank count remains
permitted when the underlying PETSc ordering and geometry/layout are compatible.

@section p58_checkpoint_sec 12. Committed Checkpoint Contract

One checkpoint is an immutable directory:

```text
output/checkpoints/step_000000001000/
  checkpoint.meta
  COMMITTED
  eulerian/
    block_0000/
      Ucat.dat
      Ucont.dat
      Ucont_rm1.dat
      P.dat
      Nvert.dat
      ...enabled model or particle-derived Eulerian fields...
  particles/                 # only when particle state is enabled
    position.dat
    velocity.dat
    ...catalogued restart fields...
```

The format identity in metadata is `picurv-checkpoint`, version 1. `picurv` is
the project name in the format identifier; it is not a new binary encoding.
Every `.dat` payload remains an ordinary PETSc binary vector and is written and
read through the existing generic field/swarm I/O functions. Eulerian vectors
use natural ordering, so compatible grid data are independent of the MPI
decomposition used to write them.

Eulerian payloads are block scoped: each block writes its own
`eulerian/block_%04d` subdirectory and each inventory entry records its block
index. Particle payloads are deliberately **not** block scoped. PICurv carries a
single swarm for the whole simulation rather than one swarm per block, so
`particles/` sits directly under the bundle root and its inventory entries
record block `-1`, meaning "not owned by any single block". This is a stated
limitation of version 1 of the format rather than an oversight. Introducing
per-block swarms would require an additional directory level and a format
version bump, exactly as the Eulerian side already carries.

`checkpoint.meta` is a small PETSc-options text manifest. It records at least:

- format version, step, authoritative physical time, `dt`, and write reason;
- block dimensions, periodicity, and a C-computed SHA-256 of natural-order grid
  coordinates plus layout-defining metadata;
- enabled LES, RANS, and particle-state flags and the particle count; and
- an indexed payload inventory with field name, block, layout, component count,
  logical type, global vector size, encoding, relative path, and byte size.

`COMMITTED` contains the SHA-256 of `checkpoint.meta`. The checksum protects the
small identity/inventory record cheaply; the implementation does not rehash
all full-field payload bytes. Every inventoried payload must exist with its
recorded size before the bundle is accepted. C then validates the actual grid
geometry before loading fields.

Writes go to a sibling
`.step_<12 digits>.incomplete.<writer-pid>` directory. After every existing
payload writer closes, rank zero writes and flushes the manifest, writes and
flushes `COMMITTED`, and atomically renames the directory to its final name.
Discovery considers only final directories whose marker, manifest hash, and
complete payload inventory validate. A crash may lose work after the previous
commit, but cannot expose a mixed flow/statistics state as restartable.

The coordinator is called for the initialized state, ordinary
`io.data_output_frequency` cadence, graceful signal/walltime shutdown, and the
final completed state. The final state is committed even when it is off the
ordinary cadence. A duplicate request for an already committed same-step,
same-time bundle is idempotent.

Core Eulerian restart state includes `Ucat`, `Ucont`, pressure, `Nvert`, and
`Ucont_rm1`; preserving `Ucont_rm1` allows BDF2 to continue without an artificial
first-order restart step. This preserves the *order* of the restart, not bitwise
identity: a restart re-applies boundary conditions to the loaded fields, so
`Ucont` and `Ucat` resume from a boundary-consistent state that differs from the
stored one at the outlet. The offset is bounded, decays, sits below solver
tolerances, and is invariant to them. See @ref 29_Maintenance_Backlog section 8
for the measurements and for the tolerance floor that statistics
restart-equivalence tests must be written against.

Catalogued optional fields are present only when the corresponding subsystem is
enabled. On `--restart-from`, a newly enabled
subsystem may initialize normally (for example particles with
`restart_mode: init`) and becomes part of all new bundles. Explicit `load`
requires its saved payloads.

`--continue` reads the requested committed bundle directly from the same run's
output and is restricted to the same physical `case.yml` identity, excluding
run-control timing and particle load/init policy. `--restart-from` copies one
validated bundle transactionally into the new run's restart staging directory;
multi-step `eulerian_field_source: load` reads committed bundles directly from
the source output. Both Python discovery and C loading reject partial bundles.

Solver iteration logs, solution-convergence logs, particle/search/scatter
metrics, runtime-memory logs, PETSc logs, and profiling files remain independent
persistent logs. They are neither payloads nor restored checkpoint state. The
future statistics subsystem adds its own accumulator vectors and window
metadata to this same bundle coordinator; it must reuse the generic vector
read/write path rather than introduce another binary serializer.

There is no reader, dual writer, or compatibility mode for the former flat
`ufieldNNNNN_0.dat` restart layout.

@section p58_post_sec 13. Postprocessing Contract

The future postprocessor surface is separate from the existing particle
`statistics_pipeline` and is intended to look like:

```yaml
field_statistics:
  - window: production
    checkpoint: latest_committed
    derive:
      - mean: Ucat
      - mean: P
      - reynolds_stress: Ucat
      - rms: Ucat
      - tke: Ucat
    outputs:
      pointwise_vtk: true
      profile_csv: true
```

The exact recipe spelling is released with the end-to-end feature. The
postprocessor loads centered state, normalizes it in one implementation, and
derives resolved covariance, RMS, and TKE. Small negative variances may be
clamped only when taking a square root and only within a documented roundoff
tolerance; stored state is not modified.

Profiles and region summaries combine pointwise weighted centered states using
the stable parallel-merge formula and existing volume/area measures. Offline
replay may feed saved instantaneous checkpoints through the same accumulation
kernels, but cannot recover omitted timesteps or infer timestamps when variable
time stepping was used.

@section p58_legacy_sec 14. Legacy Retirement

The retired surface includes:

- `Do_averaging` and its active/dead calls;
- `SimCtx::averaging` and `-averaging`;
- `Ucat_sum`, `P_sum`, `Ucat_square_sum`, and `Ucat_cross_sum` ownership;
- `ReadStatisticalFields`/`WriteStatisticalFields` and
  `su0`/`su1`/`su2`/`sp` payloads;
- `case.yml -> models.statistics.time_averaging`; and
- reserved averaged-field postprocessor output that implied the old format.

There is no legacy importer, dual-write period, or compatibility workflow.
Current users either postprocess saved instantaneous states or wait for the new
end-to-end field-statistics contract.

@section p58_reuse_sec 15. Required Reuse and Non-Redundancy

Before adding a helper, implementation must audit and reuse, where applicable:

- @ref FieldGetDescriptor, @ref FieldGetView, and catalog capabilities;
- @ref GetOwnedCellRange and existing face/component/periodic range contracts;
- `Nvert`, Jacobians, centroids, face geometry, and immersed selectors;
- existing block and MPI reduction patterns;
- current flow, model, particle, forced-output, and restart serializers;
- current convergence log writer and rolling-history behavior; and
- current dimensionalization, nodal conversion, VTK, CSV, and postprocessing
  kernels.

New code is justified only for behavior that these surfaces do not supply:
completed-state scheduling, generic centered-moment accumulation, spatial-target
mapping, scientific metadata, transactional checkpoint coordination, and
derived-statistics normalization.

@section p58_validation_sec 16. Validation and Acceptance

Required statistics tests include:

- constant scalar/vector fields with exactly zero covariance;
- analytic periodic signals with known sample- and time-weighted moments;
- high-mean/low-fluctuation precision;
- variable `dt`, schedule stride, start/end clipping, and overlapping windows;
- uninterrupted versus restart/continue equivalence;
- duplicate-event rejection and graceful shutdown timing;
- fixed and moving masks with per-point valid weights;
- nonperiodic and mixed/fully periodic domains without duplicate planes;
- cell, node, I/J/K-face, and explicitly supported staggered requests;
- serial/MPI, rank-count change, and multiblock equivalence; and
- unchanged flow, monitoring, logging, profiling, particles, and model behavior
  when field statistics are disabled.

Checkpoint acceptance additionally requires partial-write rejection,
same-step state consistency, forced/final output, continue, restart-from,
different-run-directory, and solve-to-post tests. Every implementation branch
runs its focused tests plus the relevant repository, ingress, MPI, periodic,
restart, postprocessing, and documentation gates.

@section p58_plan_sec 17. Phase Boundary

Phase 1 is the common infrastructure foundation:

1. typed Eulerian and particle field identity/layout catalog (merged);
2. legacy statistics retirement and direct `solution_monitoring` control ingress
   while preserving current logger behavior (implemented on the Phase 1 branch);
3. committed-checkpoint format and coordinator; and
4. solve/restart/post discovery and regression coverage for that coordinator.

After Phase 1 is reviewed and merged, Phase 2 receives its own detailed spec and
approval before implementation. Its settled direction is pointwise named
windows with weighted centered first/second moments, C-resolved scientific
identity, exact checkpoint continuation, and postprocessor-derived turbulence
quantities. Spatial bins, reduced histories, PDFs, and higher moments remain
later opt-in extensions built on the same target and product abstractions.

Function logging/profiling identity remains outside these phases and is governed
by @ref 59_Function_Identity_and_Observability_Specification.

@section p58_related_sec 18. Related Pages

- **@subpage 56_Field_Identity_and_Layout_Catalog**
- **@subpage 57_Future_Architecture_Specifications**
- **@subpage 09_Monitor_Reference**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@subpage 46_C_Runtime_Execution_Map**
- **@subpage 52_Run_Lifecycle_Guide**
