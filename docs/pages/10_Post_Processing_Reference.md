@page 10_Post_Processing_Reference Configuration Reference: Postprocessor YAML

@anchor _Post_Processing_Reference

For the full commented template, see:

@verbinclude master_template/master_postprocessor.yml

`post.yml` defines postprocessing input range, processing pipelines, statistics tasks, and VTK output selection.

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

Operational semantics when launched through `picurv`:
- keep `start_step` and `end_step` as the full logical analysis window you want the recipe to represent.
- `picurv run --post-process --continue --run-dir ... --post post.yml` computes an internal effective start step for the same recipe lineage, so you do not need to keep editing `start_step` during batch catch-up.
- if the recipe changes in a way that affects generated outputs, such as `step_interval`, pipeline tasks, output prefixes, or selected fields, PICurv starts from the configured `start_step` instead of inheriting completion from the previous recipe.
- on live solver runs, `end_step: -1` still means "up to the last available step", but PICurv now caps each invocation to the highest fully available contiguous source prefix before generating `post.run`.

@section p10_source_sec 3. source_data

- `source_data.directory` -> `source_directory`
- `<solver_output_dir>` is a supported placeholder resolved by `picurv`.
- for live post-processing while the solver is still running, PICurv treats a timestep as source-available only when the full required source set for the current recipe exists. For Eulerian-only recipes that means the mandatory field files; for particle/statistics recipes it also requires the particle position file for that timestep.

@section p10_pipelines_sec 4. Processing Pipelines

Eulerian tasks (`eulerian_pipeline`):
- `q_criterion` -> `ComputeQCriterion`
- `nodal_average` -> `CellToNodeAverage:<in>><out>`
- `normalize_field` -> `NormalizeRelativeField:<field>`

Global operation:
- `global_operations.dimensionalize: true` prepends `DimensionalizeAllLoadedFields`

Lagrangian tasks (`lagrangian_pipeline`):
- `specific_ke` -> `ComputeSpecificKE:<in>><out>`

@section p10_stats_sec 5. Statistics Pipeline

`statistics_pipeline` supports either:
- list form, or
- mapping with `tasks` and optional `output_prefix`

Currently supported statistics task is `msd`, which `picurv` serializes as the
`ComputeMSD` pipeline token consumed by the C dispatcher before it calls
`ComputeParticleMSD`.

Mappings:
- tasks -> `statistics_pipeline`
- `output_prefix` -> `statistics_output_prefix`

@section p10_io_sec 6. io

Mappings:
- `output_directory` + `output_filename_prefix` -> `output_prefix`
- `output_directory` + `particle_filename_prefix` -> `particle_output_prefix`
- `output_particles` -> `output_particles`
- `particle_subsampling_frequency` -> `particle_output_freq`
- `eulerian_fields` -> `output_fields_instantaneous`
- `eulerian_fields_averaged` -> `output_fields_averaged` (reserved)
- `particle_fields` -> `particle_fields_instantaneous`
- `input_extensions.eulerian` -> `eulerianExt`
- `input_extensions.particle` -> `particleExt`

Default post input extension remains `dat` unless overridden.
`statistics_pipeline.output_prefix` is independent of `io.output_directory`; bare basenames
default under `<monitor output>/statistics/`, while explicit relative or absolute paths are preserved.
When the same timestep is post-processed again, PICurv now rewrites same-step VTK/VTP outputs and rewrites same-step statistics rows so the final CSV still contains one row per step.

@section p10_next_steps_sec 7. Next Steps

Proceed to **@subpage 11_User_How_To_Guides** for goal-oriented recipes.

For mapping and extension details:
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 50_Modular_Selector_Extension_Guide**

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
