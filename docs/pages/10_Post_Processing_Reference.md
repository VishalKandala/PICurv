@page 10_Post_Processing_Reference Configuration Reference: post.yml

For the full commented template, see:

@verbinclude master_template/master_postprocessor.yml

`post.yml` defines postprocessing input range, processing pipelines, statistics tasks, and VTK output selection.

@tableofcontents

@section structure_sec 1. File Structure

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

@section run_control_sec 2. run_control

Mappings in generated `post.run`:
- `start_step` -> `startTime`
- `end_step` -> `endTime`
- `step_interval` -> `timeStep`

@section source_sec 3. source_data

- `source_data.directory` -> `source_directory`
- `<solver_output_dir>` is a supported placeholder resolved by `pic.flow`.

@section pipelines_sec 4. Processing Pipelines

Eulerian tasks (`eulerian_pipeline`):
- `q_criterion` -> `ComputeQCriterion`
- `nodal_average` -> `CellToNodeAverage:<in>><out>`
- `normalize_field` -> `NormalizeRelativeField:<field>`

Global operation:
- `global_operations.dimensionalize: true` prepends `DimensionalizeAllLoadedFields`

Lagrangian tasks (`lagrangian_pipeline`):
- `specific_ke` -> `ComputeSpecificKE:<in>><out>`

@section stats_sec 5. Statistics Pipeline

`statistics_pipeline` supports either:
- list form, or
- mapping with `tasks` and optional `output_prefix`

Currently supported task alias set maps to C kernel `ComputeMSD`.

Mappings:
- tasks -> `statistics_pipeline`
- `output_prefix` -> `statistics_output_prefix`

@section io_sec 6. io

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

@section next_steps_sec 7. Next Steps

Proceed to **@subpage 11_User_How_To_Guides** for goal-oriented recipes.

For mapping and extension details:
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
