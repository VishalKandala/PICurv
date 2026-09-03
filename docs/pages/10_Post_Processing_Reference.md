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
  and reports every derived product in physical units. It reaches three producers:
  the field pipeline scales loaded fields, the accumulator scales derived statistics,
  and the spectra generator scales its own outputs.

Lagrangian tasks (`lagrangian_pipeline`):
- `specific_ke` -> `ComputeSpecificKE:<in>><out>`

@section p10_cap_eul_sec 4.1 Eulerian Post-Processing Kernel Entries

@htmlinclude generated/capability_inventory_post_eulerian_task.html

@subsection p10_cap_eul_q_criterion_sub q_criterion

@anchor p10_cap_eul_q_criterion

**Identity.** `eulerian_pipeline: [{task: q_criterion}]` -> `ComputeQCriterion` in the `-process_pipeline` chain.

**What it does.** Computes the Q-criterion, the second invariant of the velocity-gradient tensor, as a new cell field. Positive Q marks regions where rotation dominates strain.

**When to choose it.** For vortex visualisation - isosurfaces of Q are the standard way to see coherent structures in a turbulent field. It answers 'where are the vortices', not 'how strong is the turbulence'; use field statistics for the latter.

**Parameters it owns.** None. It reads the loaded velocity field and writes one derived field.

**Interactions.** Requires `Ucat` to be loaded, which the standard source data provides. Order within `eulerian_pipeline` matters only against tasks that rewrite velocity; it does not consume the output of `nodal_average`.

**Diagnostics.** The derived field appears in the written `.vts` and in the post-processor's field listing. Its absence there means the task never ran.

**Evidence.** Unit verified - `make unit-post`.

**Limitations.** Q is a diagnostic, not a threshold: the isovalue that reveals structure is flow-dependent and this task chooses none for you.

@subsection p10_cap_eul_normalize_field_sub normalize_field

@anchor p10_cap_eul_normalize_field

**Identity.** `eulerian_pipeline: [{task: normalize_field, field: P, reference_point: [i,j,k]}]` -> `NormalizeRelativeField:<field>`.

**What it does.** Subtracts the value at a chosen reference point from a field, so what is written is the field relative to that point rather than its absolute level.

**When to choose it.** For pressure, always, in any flow without a pressure boundary condition. The incompressible pressure is determined only up to a constant, so absolute values are not comparable between runs, timesteps, or solvers - and a plot of them is a plot of an arbitrary offset.

**Parameters it owns.** `field` (currently `P` only) and `reference_point`, a three-item cell index that becomes `reference_ip/jp/kp`.

**Interactions.** Applies to the loaded field before any writer sees it. Choosing a reference point inside a recirculation or a boundary layer makes the normalized field harder to compare, not easier - prefer a quiescent interior cell.

**Diagnostics.** The written field carries the same name; the difference is visible as a shifted range. A reference point outside the block is a validation error naming the index.

**Evidence.** Unit verified - `make unit-post`.

**Limitations.** Only `P` is accepted today, and the reference point is a fixed index rather than a coordinate, so it does not follow a moving feature or survive a grid change.

@subsection p10_cap_eul_nodal_average_sub nodal_average

@anchor p10_cap_eul_nodal_average

**Identity.** `eulerian_pipeline: [{task: nodal_average, input_field: X, output_field: Y}]` -> `CellToNodeAverage:X>Y`.

**What it does.** Averages a cell-centred field onto grid nodes, writing the result as a separate named field.

**When to choose it.** When a downstream consumer expects nodal data - some visualisation filters and line-extraction tools interpolate badly from cell data - or when comparing against a reference that is defined at nodes. Leave cell fields alone otherwise: averaging is a smoothing operation and loses information.

**Parameters it owns.** `input_field` and `output_field`, which must differ. The audit rejects an in-place average because it would destroy the source the rest of the pipeline reads.

**Interactions.** Runs after any task that produces the input field, so a `q_criterion` output can be averaged by listing `nodal_average` later in the pipeline.

**Diagnostics.** Both fields appear in the written `.vts`. A missing output field means the input name did not match anything loaded.

**Evidence.** Unit verified - `make unit-post`.

**Limitations.** Averaging is a low-pass filter: peak values move toward their neighbourhood mean, so nodal fields understate extrema and should not be used for max-value claims.

@section p10_cap_lag_sec 4.2 Lagrangian Post-Processing Kernel Entries

@htmlinclude generated/capability_inventory_post_lagrangian_task.html

@subsection p10_cap_lag_specific_ke_sub specific_ke

@anchor p10_cap_lag_specific_ke

**Identity.** `lagrangian_pipeline: [{task: specific_ke, input_field: X, output_field: Y}]` -> `ComputeSpecificKE:X>Y`.

**What it does.** Computes specific kinetic energy - half the squared magnitude of a per-particle vector field - and stores it as a new scalar swarm field.

**When to choose it.** When particle energetics are the question: comparing swarm energy against the carrier flow, or watching a swarm equilibrate after seeding. For dispersion questions use the particle statistics pipeline instead.

**Parameters it owns.** `input_field` (a particle vector field, typically the particle velocity) and `output_field` (the scalar to write).

**Interactions.** Reads a swarm field, so it requires particle data in the source bundle. It does not interact with the Eulerian pipeline, which runs over a different data structure.

**Diagnostics.** The derived scalar appears in the written particle output. An input name that matches no swarm field leaves the output absent rather than zero.

**Evidence.** Unit verified - `make unit-post`.

**Limitations.** Specific energy per particle, not per unit mass of a distribution: it carries no particle mass or number weighting, so summing it across a swarm is not a physical total unless every particle represents the same mass.

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

**Derived statistics follow `global_operations.dimensionalize`.** Each carries the
source field's reference scale raised to its own power: a mean and an RMS are linear
in it, a Reynolds stress, a turbulent kinetic energy and a co-moment flux quadratic.
A field with no declared reference scale stays non-dimensional and says so in the log.
See @ref p58_derived_sec.

@section p10_cap_fso_sec 6.1 Field Statistics Output Entries

@htmlinclude generated/capability_inventory_post_field_statistics_output.html

@note **`formats` is a parameter, not a choice between behaviours.** `vtk` writes the
derived fields listed below into the window's bundle; `csv` appends one row per
processed step carrying sample count, total weight, represented time, the per-point
valid-fraction range, and the mean turbulent kinetic energy. Both are written from the
same accumulated state, so requesting both costs output volume rather than computation.
The `csv` is the artefact that answers whether a window converged; the `vtk` is the one
that shows what it converged to. `tests/tooling/family_census.json` records why this is
classified as a parameter of these entries rather than as a family of its own.

@subsection p10_cap_fso_mean_sub mean

@anchor p10_cap_fso_mean

**Identity.** `field_statistics.outputs: [mean]`.

**What it does.** Writes the accumulated first moment of each tracked field, divided by the window's total weight - the time-averaged field over the window.

**When to choose it.** Always, effectively: every other output is defined against the mean, and a mean field is what most comparisons against a reference profile need. It is also the cheapest output, since the first moment is accumulated regardless.

**Parameters it owns.** None. The fields averaged, the window bounds, and the weighting are properties of the window, not of this output.

**Interactions.** Requires an accumulated window. Reads the same state `reynolds_stress`, `rms`, and `tke` are derived from, so requesting several outputs costs one accumulation.

**Diagnostics.** Written into the window's `vtk` bundle; the `csv` row reports the sample count and total weight behind it.

**Evidence.** Unit verified - `make unit-statistics`.

**Limitations.** A mean over a window that has not converged is still a mean: the output carries no statement about whether the averaging interval was long enough. The `csv` convergence history is what answers that.

@subsection p10_cap_fso_reynolds_stress_sub reynolds_stress

@anchor p10_cap_fso_reynolds_stress

**Identity.** `field_statistics.outputs: [reynolds_stress]`.

**What it does.** Writes the full Reynolds stress tensor R_ij = C_ij / W, the centred second moment divided by total weight - six independent components.

**When to choose it.** For anything that needs the anisotropy of the turbulence rather than its magnitude: budget terms, near-wall structure, comparison against DNS stress profiles. Choose `rms` instead if only the diagonal matters, and `tke` if only the trace does.

**Parameters it owns.** None; it is derived from the accumulated second moment.

**Interactions.** Requires the window to accumulate second moments. `rms` and `tke` are functions of the same tensor, so requesting all three adds output volume, not computation.

**Diagnostics.** Six fields per window in the `vtk` bundle. One window carrying every output already writes fifteen fields against a per-file limit of twenty, which is why windows cannot share a file.

**Evidence.** Unit verified - `make unit-statistics`.

**Limitations.** Off-diagonal components converge more slowly than the diagonal, so a window long enough for `rms` is not necessarily long enough for the shear stress.

@subsection p10_cap_fso_rms_sub rms

@anchor p10_cap_fso_rms

**Identity.** `field_statistics.outputs: [rms]`.

**What it does.** Writes the root-mean-square fluctuation of each component, sqrt(R_ii) - the diagonal of the Reynolds stress tensor under a square root.

**When to choose it.** When the question is fluctuation magnitude per direction: turbulence intensity profiles, and the usual first comparison against a reference. It is `reynolds_stress` with the cross terms dropped and a square root applied.

**Parameters it owns.** None.

**Interactions.** Derived from the same second moment as `reynolds_stress` and `tke`. Small negative variances from floating-point cancellation are clamped only when taking the square root, and only within a documented tolerance - stored state is never mutated.

**Diagnostics.** Three fields per window. The `csv` row's valid-fraction range is the check on whether every point had enough samples.

**Evidence.** Unit verified - `make unit-statistics`.

**Limitations.** The clamp means an RMS of exactly zero can mean 'no fluctuation' or 'variance below the cancellation floor'; the sample count distinguishes them.

@subsection p10_cap_fso_tke_sub tke

@anchor p10_cap_fso_tke

**Identity.** `field_statistics.outputs: [tke]`.

**What it does.** Writes turbulent kinetic energy, k = (R_xx + R_yy + R_zz) / 2 - half the trace of the Reynolds stress tensor.

**When to choose it.** For a single scalar summarising turbulence level, and for anything that plots decay or growth of turbulence over time. It is the natural quantity for isotropic decay studies, where the anisotropy `reynolds_stress` reports is not the point.

**Parameters it owns.** None.

**Interactions.** A scalar contraction of the same tensor `reynolds_stress` and `rms` come from. The `csv` output reports the mean value of this field per processed step, which is the convergence signal for the window as a whole.

**Diagnostics.** One field per window, plus the per-step mean in the `csv`.

**Evidence.** Unit verified - `make unit-statistics`.

**Limitations.** Being a trace, it is blind to anisotropy: two flows with very different stress structure can report the same k.

@subsection p10_cap_fso_flux_sub flux

@anchor p10_cap_fso_flux

**Identity.** `field_statistics.outputs: [flux]`.

**What it does.** Writes the centred cross-moment between velocity and a scalar, u_i'psi' = C_ipsi / W - the turbulent flux of that scalar.

**When to choose it.** When a scalar is being transported and the question is how turbulence moves it: heat flux, species flux, any gradient-transport closure being assessed. It requires a scalar to exist, so it is inert in a pure momentum run.

**Parameters it owns.** None of its own; which scalar is tracked is a property of the window.

**Interactions.** Requires scalar transport to be active and the scalar included in the window's tracked fields. Otherwise the cross-moment has nothing to accumulate.

**Diagnostics.** Three fields per window, one per velocity component.

**Evidence.** Unit verified - `make unit-statistics`.

**Limitations.** A cross-moment converges more slowly than either factor's own variance, so flux is the output most likely to be under-converged in a given window.

@section p10_stat_entries_sec 7. Particle Statistics Task Entries

@htmlinclude generated/capability_inventory_statistics_task.html

@subsection p10_cap_stat_msd_sub msd

@anchor p10_cap_stat_msd

**Identity.** Statistics pipeline task `msd` -> `ComputeMSD`.

**What it does.** Computes the mean squared displacement of the particle swarm over time
and writes it as a CSV time series.

**When to choose it.** Characterizing dispersion. MSD against time is the standard way to
read a diffusion coefficient out of a particle simulation, which is what the Brownian and
drift verification cases do.

**Parameters it owns.** The statistics pipeline block that selects it; the output file
name follows the configured statistics prefix.

**Interactions.** Requires particles. It reads particle positions only, so it is unaffected
by the Eulerian field source.

**Diagnostics.** The emitted CSV is the diagnostic: a linear MSD in time indicates
diffusive behaviour, a quadratic one indicates ballistic transport.

**Evidence.** Unit verified - `make unit-statistics` covers the MSD kernel including its
empty-swarm behaviour.

**Limitations.** The only statistics task currently exposed. It is a whole-swarm measure
with no spatial conditioning.

@section p10_spectra_sec 8. spectra

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
`<run.analysis>/spectra/initial_condition_spectrum.csv` when the run has one. Pass a task-name
substring when a recipe measured more than one spectrum. The per-step scalars are
ordinary series, so `--list-plot-series` finds them and `--plot` draws them.

**A decaying flow is a different statistical state at every step**, so the family is
plotted rather than averaged, and the stage never averages across steps. Smoothing a
spectrum comes from ensembling over seeds at matched times, not from pooling times.

@section p10_cap_spec_sec 8.1 Spectra Task Entries

@htmlinclude generated/capability_inventory_post_spectra_task.html

@subsection p10_cap_spec_shell_spectrum_sub shell_spectrum

@anchor p10_cap_spec_shell_spectrum

**Identity.** `spectra.tasks: [{task: shell_spectrum, field: Ucat, block: 0}]`, measured by `generators/spectra.gen` rather than by the C post-processor.

**What it does.** Computes a three-dimensional energy spectrum by transforming the velocity field and binning the result into spherical shells in wavenumber space, producing E(k).

**When to choose it.** For homogeneous turbulence where the question is how energy is distributed across scales - resolving whether an inertial range exists, or where the dissipation cutoff sits. It is the wrong tool for an inhomogeneous flow, where a shell average mixes physically distinct directions.

**Parameters it owns.** `field` (the vector field to transform), `block`, `symbol` (`continuum` or `discrete`, the binning abscissa), `subtract_mean` (`none`, `domain`, or `window:<name>`), and `mean_source_step` when a window mean is used.

**Interactions.** Requires a uniformly spaced, geometrically periodic, single-block domain, and `picurv validate` checks all three against the case **before any field is read**. Independent of `eulerian_pipeline`: it reads raw checkpoint payloads, so no field output stage is needed. `--only spectra` runs it alone.

**Diagnostics.** Writes one file per task per checkpoint under the configured `output_prefix`. A precondition failure is a validation error naming the violated requirement, not a silent empty spectrum.

**Evidence.** Implemented only. No shipped case gates a numerical acceptance on a measured spectrum, so nothing establishes the binning against a reference result.

**Limitations.** Experimental. Shell averaging assumes isotropy that a wall-bounded or bent geometry does not have, and `subtract_mean: window:<name>` depends on an accumulated window existing at the named step. No convergence or windowing guidance is established.

@section p10_io_sec 9. io

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

@section p10_next_steps_sec 10. Next Steps

Proceed to **@subpage 11_User_How_To_Guides** for goal-oriented recipes.

For mapping and extension details:
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 50_Modular_Selector_Extension_Guide**
- **@subpage 58_Field_Statistics** (window semantics and derived-quantity definitions)
