@page 07_Case_Reference Configuration Reference: Case YAML

@anchor _Case_Reference

For the full commented template, see:

@verbinclude master_template/master_case.yml

`case.yml` defines physical setup, grid source, domain topology, and boundary conditions.
It is intentionally modular: the same `case.yml` can be paired with different `solver.yml`, `monitor.yml`, and `post.yml` profiles when the combination remains contract-compatible.

@tableofcontents

@section p07_properties_sec 1. properties

```yaml
properties:
  scaling:
    length_ref: 0.1
    velocity_ref: 1.5
  fluid:
    density: 1000.0
    viscosity: 0.001
  initial_conditions:
    mode: "Constant"
    u_physical: 1.5
    v_physical: 0.0
    w_physical: 0.0
```

Key mappings:
- `scaling.length_ref` -> `-scaling_L_ref`
- `scaling.velocity_ref` -> `-scaling_U_ref`
- `fluid.density` and `fluid.viscosity` are used by `picurv` to compute Reynolds number -> `-ren`
- `initial_conditions.mode` -> `-finit` (`Zero`, `Constant`, `Poiseuille`)
- `u_physical/v_physical/w_physical` -> `-ucont_x/-ucont_y/-ucont_z`
- `peak_velocity_physical` (Poiseuille only) -> mapped by `picurv` to the inlet-aligned `-ucont_*` component

For the scaling model and conversion logic, see **@subpage 19_Nondimensionalization**.
For detailed startup behavior of field initialization modes, see **@subpage 33_Initial_Conditions**.

Practical contract notes:

- `initial_conditions.mode` should be set explicitly. The launcher now requires the choice instead of silently inferring it.
- `mode: "Zero"` may omit velocity components entirely.
- `mode: "Constant"` requires explicit `u_physical/v_physical/w_physical`.
- `mode: "Poiseuille"` may use either:
  - `peak_velocity_physical` for the scalar centerline-speed input, or
  - `u_physical/v_physical/w_physical` for an explicit component override.

@section p07_run_control_sec 2. run_control

```yaml
run_control:
  dt_physical: 0.0001
  start_step: 0
  total_steps: 2000
```

Mappings:
- `dt_physical` -> `-dt` (non-dimensionalized)
- `start_step` -> `-start_step`
- `total_steps` -> `-totalsteps`

Restart is handled via CLI flags rather than `case.yml` keys:
- `--restart-from <previous_run_dir>` -> `picurv` resolves the prior run's actual restart source directory and emits `-restart_dir <absolute_previous_source>`
- `--continue` -> shorthand for resuming from the most recent run of the same case

@section p07_grid_sec 3. grid

Supported modes:
- `programmatic_c`
- `file`
- `grid_gen`

Mode compatibility note:

- for normal `solve` and `load` workflows, all three grid modes are supported.
- for `solver.yml -> operation_mode.eulerian_field_source: analytical`, the current contract requires `grid.mode: programmatic_c`.

Optional global DMDA layout hints apply to all grid modes:

```yaml
grid:
  da_processors_x: 4
  da_processors_y: 2
  da_processors_z: 2
```

These are scalar global values, not per-block vectors. Legacy placement under
`grid.programmatic_settings.da_processors_*` is still accepted for compatibility,
but the shared top-level `grid.da_processors_*` form is preferred.

@subsection p07_grid_prog_ssec 3.1 mode: programmatic_c

`programmatic_settings` supports per-block lists for geometry arrays:
- `im/jm/km`
- `xMins/xMaxs`, `yMins/yMaxs`, `zMins/zMaxs`
- `rxs/rys/rzs`

Dimension contract:
- `im/jm/km` in YAML are cell counts.
- `picurv` converts them to node counts before emitting `-im/-jm/-km` for the C runtime.

Important constraint:
- `grid.da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block processor decomposition is not implemented.

@subsection p07_grid_file_ssec 3.2 mode: file

```yaml
grid:
  mode: file
  source_file: my_grid.picgrid
```

`picurv` validates existence and prepares normalized grid data for C-side ingestion.

Optional legacy conversion path (headerless 1D-axis payloads):

```yaml
grid:
  mode: file
  source_file: legacy_flat.grid
  legacy_conversion:
    enabled: true
    format: legacy1d           # aliases: legacy_1d, les_flat_1d, les-flat-1d
    axis_columns: [0, 1, 2]    # preferred source columns for X/Y/Z axis rows
    strict_trailing: true
```

When enabled, `picurv` first calls `scripts/grid.gen legacy1d` to create a canonical PICGRID file in the run config, then runs the normal validation/non-dimensionalization staging path.

@subsection p07_grid_gen_ssec 3.3 mode: grid_gen

```yaml
grid:
  mode: grid_gen
  generator:
    config_file: config/grids/coarse_square_tube_curved.cfg
    grid_type: cpipe
    cli_args: ["--ncells-i", "96", "--ncells-j", "96"]
```

This runs `scripts/grid.gen` before solver launch and stages generated grid artifacts into the run config.
`grid.generator.config_file` is required today; `picurv` does not synthesize a temporary grid config.
`grid.gen` accepts cell-count inputs (`ncells_*` / `--ncells-*`) and writes node counts into the generated `.picgrid` header.

For direct `grid.gen` usage, generator types, and config-file structure, see **@subpage 48_Grid_Generator_Guide**.

@section p07_models_sec 4. models

```yaml
models:
  domain:
    blocks: 1
    i_periodic: false
    j_periodic: false
    k_periodic: false
  physics:
    dimensionality: "3D"
    turbulence:
      les: false
    particles:
      count: 0
      init_mode: "Surface"
      restart_mode: "init"
```

Common mappings:
- `domain.blocks` -> `-nblk`
- periodic flags -> `-i_periodic/-j_periodic/-k_periodic`
- `physics.dimensionality: "2D"` -> `-TwoD 1`
- `physics.turbulence.les` -> `-les`
- `physics.particles.count` -> `-numParticles`
- `physics.particles.init_mode` -> `-pinit` (`Surface`, `Volume`, `PointSource`, `SurfaceEdges`)
- `physics.particles.restart_mode` -> `-particle_restart_mode`
- point source coordinates -> `-psrc_x/-psrc_y/-psrc_z`

Restart note:

- if `run_control.start_step > 0`, particles are enabled, and `restart_mode` is omitted, `picurv` warns that C will default to `load`.

For mode-specific particle behavior and restart flow, see **@subpage 45_Particle_Initialization_and_Restart**.

@section p07_bc_sec 5. boundary_conditions

Single-block syntax: list of 6 face entries.
Multi-block syntax: list-of-lists, one 6-face list per block.

Supported face names:
- `-Xi`, `+Xi`, `-Eta`, `+Eta`, `-Zeta`, `+Zeta`

Supported type/handler combinations:
- `INLET` + `constant_velocity` (`vx/vy/vz`)
- `INLET` + `parabolic` (`v_max`)
- `OUTLET` + `conservation`
- `WALL` + `noslip`
- `PERIODIC` + `geometric`
- `PERIODIC` + `constant_flux` (`target_flux`, optional `apply_trim`)

All six faces must be explicitly provided for each block.
For detailed handler semantics, validation constraints, and C dispatch path, see **@subpage 44_Boundary_Conditions_Guide**.

@section p07_passthrough_sec 6. solver_parameters (Advanced)

Optional escape hatch for flags not yet exposed in structured schema:

```yaml
solver_parameters:
  -read_fields: true
  -some_new_flag: 123
```

Use sparingly and prefer structured keys when available.

@section p07_modular_sec 7. Mixing With Other Profiles

`case.yml` is designed to be combined with reusable profiles for the other config roles.

Common patterns:

- one `case.yml` + multiple `monitor.yml` files (debug vs production output),
- one `case.yml` + multiple `post.yml` recipes (quick scalar check vs heavy VTK/statistics),
- one `solver.yml` reused across many `case.yml` files,
- one `cluster.yml` reused across many runs and sweeps.

For worked combinations, see **@subpage 49_Workflow_Recipes_and_Config_Cookbook**.

@section p07_next_steps_sec 8. Next Steps

Proceed to **@subpage 08_Solver_Reference**.

Cross-file contract/mapping:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 48_Grid_Generator_Guide**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 32_Analytical_Solutions**
- **@subpage 33_Initial_Conditions**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Configuration Reference: Case YAML** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
