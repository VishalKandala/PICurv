@page 14_Config_Contract Configuration Contract (YAML -> Generated Artifacts -> Runtime)

This page is the user-facing source of truth for the configuration contract implemented by `pic.flow`.

@tableofcontents

@section inputs_sec 1. Required Input Roles

PIC-flow composes a run from five logical inputs:

1. `case.yml`: physics, grid, BC definitions, run control.
2. `solver.yml`: numerical strategy and solver parameters.
3. `monitor.yml`: I/O and logging/profiling controls.
4. `post.yml`: post-processing recipe.
5. MPI launch settings (`-n`, executable stage selection).

You can name files however you want. File names are not hardcoded on the C side; `pic.flow` resolves paths and emits generated artifacts.

@section artifacts_sec 2. Generated Artifacts

For each run, `pic.flow` generates:

- `<run_id>.control`: master PETSc/control flags for solver/post setup.
- `bcs.run` or `bcs_block*.run`: boundary condition definitions.
- `whitelist.run`: logging function allow-list.
- `profile.run`: profiling critical-function list.
- `post.run`: key=value post-processing recipe consumed by C post parser.

@section case_sec 3. Case Contract Highlights

- `grid.mode` supports: `file`, `programmatic_c`, `grid_gen`.
- For `programmatic_c`, per-block arrays are supported for geometry (`im/jm/km`, bounds, stretching).
- `da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block MPI decomposition is not currently supported.
- `boundary_conditions` supports single-block list or multi-block list-of-lists.
- `solver_parameters` is an advanced passthrough map for raw flags not yet modeled in schema.

@section solver_sec 4. Solver Contract Highlights

- `operation_mode.eulerian_field_source` -> `-euler_field_source`
- `operation_mode.analytical_type` -> `-analytical_type`
- `strategy.momentum_solver` -> `-mom_solver_type` via normalized names.
- Solver-specific block support currently includes `momentum_solver.dual_time_picard_rk4`.
- `petsc_passthrough_options` remains the escape hatch for advanced PETSc/C flags.

@section monitor_sec 5. Monitor Contract Highlights

- `io.data_output_frequency` -> `-tio`
- `io.particle_log_interval` -> `-logfreq`
- `io.directories.output/restart/log` -> `-output_dir/-restart_dir/-log_dir`
- `io.directories.eulerian_subdir/particle_subdir` -> `-euler_subdir/-particle_subdir`
- `solver_monitoring` maps raw flags directly into control output.

@section post_sec 6. Post Contract Highlights

- Pipelines are serialized into semicolon-delimited C pipeline strings.
- `io.eulerian_fields` -> `output_fields_instantaneous`
- `io.eulerian_fields_averaged` -> `output_fields_averaged` (reserved/no-op in current writer path)
- `io.particle_fields` -> `particle_fields_instantaneous`
- `io.input_extensions.eulerian/particle` -> `eulerianExt/particleExt` for post input readers
- `source_data.directory` -> `source_directory`

@section passthrough_sec 7. Escape Hatches and Defaults

- Escape hatches stay supported:
  - `case.solver_parameters`
  - `solver.petsc_passthrough_options`
  - `monitor.solver_monitoring`
- Default post input/output extension is `dat` unless overridden.
- Missing optional keys preserve existing C defaults.
