@page 10_Post_Processing_Reference Configuration Reference: Postprocessor YAML

@anchor _Post_Processing_Reference

For the full commented template, see:

@verbinclude master_template/master_postprocessor.yml

`post.yml` defines postprocessing input range, processing pipelines, statistics tasks, and VTK output selection.

It carries two separate statistics contracts, which are deliberately not
combined. `statistics_pipeline` reduces particle trajectories; `field_statistics`
derives the Eulerian field windows the solver accumulated. Particle MSD and
resolved turbulence statistics answer different questions and share no state.
A third block, `spectra`, measures energy spectra from the instantaneous velocity;
unlike every other block here it runs in the conductor rather than in the C
post-processor.

@tableofcontents

@section p10_structure_sec 1. File Structure

```yaml
run_control:
  start_step: 100
  end_step: 1000
  step_interval: 100

source_data:
  directory: "<solver_output_dir>"

global_operations:
  dimensionalize: true

eulerian_pipeline:
  - task: nodal_average
    input_field: Ucat
    output_field: Ucat_nodal
  - task: q_criterion

lagrangian_pipeline:
  - task: specific_ke
    input_field: velocity
    output_field: SpecificKE

statistics_pipeline:
  output_prefix: "Stats"
  tasks:
    - task: msd

field_statistics:
  windows: [production]
  outputs: [mean, reynolds_stress, rms, tke, flux]
  formats: [vtk, csv]

spectra:
  output_prefix: "Spectrum"
  tasks:
    - task: shell_spectrum
      field: Ucat

io:
  output_directory: "viz"
  output_filename_prefix: "Field"
  particle_filename_prefix: "Particle"
  output_particles: true
  eulerian_fields: [Ucat_nodal, Qcrit]
  particle_fields: [velocity, SpecificKE]
```

@section p10_run_control_sec 2. run_control

Mappings in generated `post.run`:
- `start_step` -> `startTime`
- `end_step` -> `endTime`
- `step_interval` -> `timeStep`

**Requested steps must be steps the solver committed.** The solver writes a
checkpoint every `monitor.yml -> io.data_output_frequency` completed steps, plus
the initial and final states. `step_interval` must therefore be a multiple of that
cadence; a finer stride names steps that were never written, and post-processing
stops at the first one it cannot find rather than skipping it. `picurv validate`
rejects a mismatch, and warns when `start_step` is off cadence — that one is only
valid if it is the run's own starting step, which is committed off cadence.

Operational semantics when launched through `picurv`:
- keep `start_step` and `end_step` as the full logical analysis window you want the recipe to represent.
- `picurv run --post-process --continue --run-dir ... --post post.yml` computes an internal effective start step for the same recipe lineage, so you do not need to keep editing `start_step` during batch catch-up.
- if the recipe changes in a way that affects generated outputs, such as `step_interval`, pipeline tasks, output prefixes, or selected fields, PICurv starts from the configured `start_step` instead of inheriting completion from the previous recipe.
- on live solver runs, `end_step: -1` still means "up to the last available step", but PICurv now caps each invocation to the highest fully available contiguous source prefix before generating `post.run`.

@section p10_source_sec 3. source_data

- `source_data.directory` -> `source_directory`
- `<solver_output_dir>` is a supported placeholder resolved by `picurv`.
- for live post-processing while the solver is still running, PICurv treats a
  timestep as source-available only when its checkpoint bundle has a valid
  `checkpoint.meta`/`COMMITTED` pair and complete payload inventory. A recipe
  that needs particles also requires committed particle state.

@section p10_pipelines_sec 4. Processing Pipelines

Eulerian tasks (`eulerian_pipeline`):
- `q_criterion` -> `ComputeQCriterion`
- `nodal_average` -> `CellToNodeAverage:<in>><out>`
- `normalize_field` -> `NormalizeRelativeField:<field>`

Global operation:
- `global_operations.dimensionalize: true` prepends `DimensionalizeAllLoadedFields`

Lagrangian tasks (`lagrangian_pipeline`):
- `specific_ke` -> `ComputeSpecificKE:<in>><out>`

@section p10_stats_sec 5. Particle Statistics Pipeline

This is the particle-reduction path. For Eulerian turbulence statistics see
@ref p10_field_stats_sec.

`statistics_pipeline` supports either:
- list form, or
- mapping with `tasks` and optional `output_prefix`

Currently supported statistics task is `msd`, which `picurv` serializes as the
`ComputeMSD` pipeline token consumed by the C dispatcher before it calls
`ComputeParticleMSD`.

Mappings:
- tasks -> `statistics_pipeline`
- `output_prefix` -> `statistics_output_prefix`

@section p10_field_stats_sec 6. field_statistics

Derives the accumulated field windows configured at
`monitor.yml -> field_statistics` into scientific turbulence quantities. Omit the
block to derive nothing.

```yaml
field_statistics:
  windows: [production]
  source_step: 4000        # optional; defaults to the step being processed
  outputs: [mean, reynolds_stress, rms, tke, flux]
  formats: [vtk, csv]
```

Mappings:
- `windows` -> `field_statistics_windows`
- `outputs` -> `field_statistics_outputs`
- `formats` -> `field_statistics_formats`
- `source_step` -> `field_statistics_source_step`

The recipe **names** windows rather than redescribing them. The post-processor is
launched with the run's own control file and resolves the window definitions
through the same parser the solver used, so there is no second description to fall
out of step with the first.

Available outputs, each resolved against what the window actually accumulated:

| Output | Produces | Requires |
| --- | --- | --- |
| `mean` | the window mean of every accumulated field | — |
| `reynolds_stress` | centered second moments; six components for a vector | `moments: [second]` |
| `rms` | per-component root mean square of the fluctuation | `moments: [second]` |
| `tke` | turbulent kinetic energy | `moments: [second]` on a vector |
| `flux` | the cross-field co-moments | a `covariances` entry |

A recipe that would produce an empty file — asking for `rms` from a window that
kept only first moments, for instance — is rejected before the run rather than
after it.

`vtk` writes one `.vts` per window per processed step. `csv` appends one row per
processed step to a per-window convergence history recording sample count,
accumulated weight, represented time, mask coverage, and the domain-mean TKE. The
CSV is what answers whether a window has run long enough; no single snapshot can.

Because the pipeline runs once per processed step, a recipe spanning a range
produces a history. Pinning `source_step` reads every processed step's output from
that one bundle instead, collapsing the history to a single answer.

Each window writes its own files, so a window name may not be listed twice.

**Derived statistics are not dimensionalized**, even under
`global_operations.dimensionalize`. A Reynolds stress scales as velocity squared
and a co-moment as the product of two different scales, and the per-field scaling
table expresses neither. See @ref p58_derived_sec.

@section p10_spectra_sec 7. spectra

Measures turbulent energy spectra from the instantaneous velocity in each committed
checkpoint. Omit the block to measure nothing.

```yaml
spectra:
  output_prefix: "Spectrum"
  tasks:
    - task: shell_spectrum
      field: Ucat
      block: 0
      symbol: continuum          # continuum | discrete
      subtract_mean: none        # none | domain | window:<name>
      mean_source_step: 4000     # optional; only with window:<name>
```

Unlike every other block on this page, spectra are **not** computed by the C
post-processor. They are measured by `generators/spectra.gen`, which the conductor
runs against the raw checkpoint payloads, so the stage needs no field output and does
not depend on `eulerian_pipeline`. The only key that reaches the C recipe is
`spectra_signature`, a digest the post-processor accepts and ignores; it exists so a
changed spectra recipe reaches the recipe fingerprint that `--continue` compares.

@subsection p10_spectra_only_sub Running spectra alone

```bash
picurv run --post-process --only spectra --run-dir runs/my_run --post post.yml
```

`--only` selects among the post stages, `fields` (the field post-processor) and
`spectra`, and defaults to both. Re-measuring a spectrum is cheap; rebuilding the
field output is not, so the analysis loop belongs behind `--only spectra`.

@subsection p10_spectra_preconditions_sub Preconditions

A spectrum is only defined where the transform direction is uniformly spaced,
periodic, and statistically homogeneous. Each task declares what it needs, and
`picurv validate` checks it against the case **before any field is read**:

| Task | Status | Requires |
| --- | --- | --- |
| `shell_spectrum` | implemented | every face `PERIODIC`, a single block, and a uniform axis-aligned Cartesian grid |
| `plane_spectrum` | planned | two periodic axes; see @ref p60_spectra_partial_sec |
| `line_spectrum` | planned | one periodic axis; see @ref p60_spectra_partial_sec |
| `temporal_spectrum` | planned | no homogeneous direction needed; see @ref p60_spectra_temporal_sec |

Only `shell_spectrum` exists today, so a case homogeneous in one or two directions —
a channel, a straight duct, a boundary layer — has a spatial spectrum the pipeline
cannot yet produce. The task table is built to hold the others; the reservation is
recorded rather than implied.

Boundary conditions and block count are checked from `case.yml`; grid uniformity is
checked by the generator against the staged PICGRID, which is the only place the node
coordinates are known. A wall-bounded or curvilinear case is refused by name rather
than yielding a curve with no meaning:

```
spectra task 'shell_spectrum' requires every face of block 0 to be PERIODIC,
because a shell-averaged spectrum is only defined for a triply periodic
homogeneous box. Non-periodic faces: ['+Eta', '+Xi', ...].
```

This is why a 90-degree bend cannot produce a spatial spectrum at all: it develops
streamwise and is bounded on all four sides, so it has no homogeneous direction and
needs @ref p60_spectra_temporal_sec. A *straight* duct or channel is a different case —
it does have homogeneous directions, and awaits @ref p60_spectra_partial_sec rather
than the temporal path.

@subsection p10_spectra_keys_sub Keys

- `symbol` selects the binning abscissa. `continuum` bins by \f$|k|\f$;
  `discrete` bins by the centered-difference symbol \f$\sin(k\,\Delta x)/\Delta x\f$,
  which is what the solver's own operator resolves. Comparing a measured rolloff
  against the continuum abscissa alone will attribute the scheme's damping to physics.
- `subtract_mean` chooses the fluctuation. `none` transforms the field as stored;
  `domain` removes the volume mean of that snapshot; `window:<name>` subtracts the
  accumulated mean of a field-statistics window, which is the better estimate wherever
  the flow is not homogeneous. A window mean is only defined where the window actually
  sampled, so a window that has not accumulated over the domain is refused rather than
  subtracted as zero.
- `mean_source_step` pins the bundle the window mean is read from, independently of the
  step being transformed. It applies only to `window:<name>`. Left unset, each step
  reads its own bundle, whose mean is the *running* mean at that step — undefined at the
  moment the window activates, and built from few samples shortly after. Pin it to a
  step at or after the window closes to subtract the converged mean, which is almost
  always what is wanted. This mirrors `field_statistics.source_step`.
- `block` selects the block; `field` is `Ucat`. `Ucont` is component-staggered and has
  no single cell-centered spectrum, per @ref p60_products_sec.

Two tasks differing in task, field, block, or symbol write different files and may
coexist; two identical tasks are refused, since the second would overwrite the first.

@subsection p10_spectra_outputs_sub Outputs

Each task writes two CSVs under `<monitor output>/spectra/`:

| File | Columns | Purpose |
| --- | --- | --- |
| `<prefix>_<task>_<field>_block<nnnn>_<symbol>.csv` | `step,time,k,energy` | the spectra themselves, long format |
| `..._history.csv` | `step,time,` then the scalars below | one row per processed step |

The spectrum file is long format so a family of curves stays one file, which is also
what makes comparison across a sweep a glob rather than a special case.

Recorded per step: `resolved_kinetic_energy`, `spectrum_total_energy`,
`parseval_residual`, `spectrum_peak_k`, `zero_mode_energy`, `integral_length_scale`,
`taylor_microscale`, and `dissipation_over_viscosity`.

`parseval_residual` is a self-check rather than an input: summed shell energy must
equal the resolved kinetic energy of the transformed field, and a residual above
round-off means the measurement is wrong.

Dissipation is reported **divided by the kinematic viscosity** —
\f$\varepsilon = \nu \times\f$ `dissipation_over_viscosity` — because the
measurement needs no fluid properties. Both length scales need none either;
\f$\nu\f$ cancels out of their isotropic definitions.

@subsection p10_spectra_crosscheck_sub Cross-checking against field statistics

`parseval_residual` verifies the transform against itself. An independent check is
available against @ref p10_field_stats_sec, which computes the same turbulent kinetic
energy through an entirely separate path — a C accumulator over every accepted state
rather than a Python transform of written snapshots:

```yaml
field_statistics:
  windows: [early_decay]
  outputs: [tke]
  formats: [csv]
  source_step: 400
spectra:
  tasks:
    - task: shell_spectrum
      field: Ucat
      subtract_mean: window:early_decay
      mean_source_step: 400
```

Time-average the spectra history's `resolved_kinetic_energy` across the window and
compare with `mean_tke` in the statistics CSV. Weight by the represented interval, not
by sample: the window itself weights by physical time, and a plain mean over a curved
history is biased by its endpoints.

Exact agreement is not expected and is not the point. The window accepts states at its
own `step_cadence`, while spectra can only read the steps a checkpoint exists at, so
the two sides average different sample sets. Agreement to a few percent confirms the
paths measure the same quantity; a large gap does not.

@subsection p10_spectra_plot_sub Viewing them

```bash
picurv summarize --run-dir runs/my_run --plot-spectrum
picurv summarize --run-dir runs/my_run --plot spectra.resolved_kinetic_energy
```

`--plot-spectrum` draws one curve per processed step on log-log axes, overlaid on
`diagnostics/initial_condition_spectrum.csv` when the run has one. Pass a task-name
substring when a recipe measured more than one spectrum. The per-step scalars are
ordinary series, so `--list-plot-series` finds them and `--plot` draws them.

**A decaying flow is a different statistical state at every step**, so the family is
plotted rather than averaged, and the stage never averages across steps. Smoothing a
spectrum comes from ensembling over seeds at matched times, not from pooling times.

@section p10_io_sec 8. io

Mappings:
- `output_directory` + `output_filename_prefix` -> `output_prefix`
- `output_directory` + `particle_filename_prefix` -> `particle_output_prefix`
- `output_particles` -> `output_particles`
- `particle_subsampling_frequency` -> `particle_output_freq`
- `eulerian_fields` -> `output_fields_instantaneous`
- `particle_fields` -> `particle_fields_instantaneous` — **must be non-empty when `output_particles: true`**; an empty list produces no particle VTP output even if particle output is enabled. Standard swarm field names: `position`, `velocity`, `pid`, `CellID`, `weight`.
- `input_extensions.eulerian/particle`, when present, must be `dat`; committed
  bundle payload names are fixed

The postprocessor reads the catalogued `.dat` PETSc vectors inventoried by each
committed checkpoint; input extensions are not runtime-selectable.
`statistics_pipeline.output_prefix` is independent of `io.output_directory`; bare basenames
default under `<monitor output>/statistics/`, while explicit relative or absolute paths are preserved.
When the same timestep is post-processed again, PICurv now rewrites same-step VTK/VTP outputs and rewrites same-step statistics rows so the final CSV still contains one row per step.

@section p10_next_steps_sec 9. Next Steps

Proceed to **@subpage 11_User_How_To_Guides** for goal-oriented recipes.

For mapping and extension details:
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 50_Modular_Selector_Extension_Guide**
- **@subpage 58_Field_Statistics** (window semantics and derived-quantity definitions)

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Configuration Reference: Postprocessor YAML** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
