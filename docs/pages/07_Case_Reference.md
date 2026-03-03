@page 07_Case_Reference Configuration Reference: Case YAML

For the full commented template, see:

@verbinclude master_template/master_case.yml

`case.yml` defines physical setup, grid source, domain topology, and boundary conditions.
It is intentionally modular: the same `case.yml` can be paired with different `solver.yml`, `monitor.yml`, and `post.yml` profiles when the combination remains contract-compatible.

@tableofcontents

@section properties_sec 1. properties

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

@section run_control_sec 2. run_control

```yaml
run_control:
  dt_physical: 0.0001
  start_step: 0
  total_steps: 2000
  restart_from_run_dir: "../runs/flat_channel_20260303-120000"  # optional
```

Mappings:
- `dt_physical` -> `-dt` (non-dimensionalized)
- `start_step` -> `-start_step`
- `total_steps` -> `-totalsteps`
- `restart_from_run_dir` -> `picurv` resolves the prior run's actual restart source directory and emits `-restart_dir <absolute_previous_source>`

@section grid_sec 3. grid

Supported modes:
- `programmatic_c`
- `file`
- `grid_gen`

Mode compatibility note:

- for normal `solve` and `load` workflows, all three grid modes are supported.
- for `solver.yml -> operation_mode.eulerian_field_source: analytical`, the current contract requires `grid.mode: programmatic_c`.

@subsection grid_prog_ssec 3.1 mode: programmatic_c

`programmatic_settings` supports per-block lists for geometry arrays:
- `im/jm/km`
- `xMins/xMaxs`, `yMins/yMaxs`, `zMins/zMaxs`
- `rxs/rys/rzs`

Dimension contract:
- `im/jm/km` in YAML are cell counts.
- `picurv` converts them to node counts before emitting `-im/-jm/-km` for the C runtime.

Important constraint:
- `da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block processor decomposition is not implemented.

@subsection grid_file_ssec 3.2 mode: file

```yaml
grid:
  mode: file
  source_file: my_grid.picgrid
```

`picurv` validates existence and prepares normalized grid data for C-side ingestion.

@subsection grid_gen_ssec 3.3 mode: grid_gen

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

@section models_sec 4. models

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

@section bc_sec 5. boundary_conditions

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

@section passthrough_sec 6. solver_parameters (Advanced)

Optional escape hatch for flags not yet exposed in structured schema:

```yaml
solver_parameters:
  -read_fields: true
  -some_new_flag: 123
```

Use sparingly and prefer structured keys when available.

@section modular_sec 7. Mixing With Other Profiles

`case.yml` is designed to be combined with reusable profiles for the other config roles.

Common patterns:

- one `case.yml` + multiple `monitor.yml` files (debug vs production output),
- one `case.yml` + multiple `post.yml` recipes (quick scalar check vs heavy VTK/statistics),
- one `solver.yml` reused across many `case.yml` files,
- one `cluster.yml` reused across many runs and sweeps.

For worked combinations, see **@subpage 49_Workflow_Recipes_and_Config_Cookbook**.

@section next_steps_sec 8. Next Steps

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
