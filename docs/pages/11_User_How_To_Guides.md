@page 11_User_How_To_Guides User How-To Guides

Quick recipes for common PICurv tasks.

@tableofcontents

@section setup_sec 1. Setup and Physics

@subsection reynolds_ssec 1.1 Change Reynolds Number

Edit `case.yml -> properties` values used by `Re = rho * U * L / mu`:

```yaml
properties:
  scaling:
    length_ref: 0.1
    velocity_ref: 1.5
  fluid:
    density: 1000.0
    viscosity: 0.001
```

@subsection twod_ssec 1.2 Run in 2D

```yaml
models:
  physics:
    dimensionality: "2D"

grid:
  mode: programmatic_c
  programmatic_settings:
    km: [3]
```

@subsection gridres_ssec 1.3 Increase Grid Resolution

For `programmatic_c`, adjust `im/jm/km`.
For `file`, provide a finer external grid and update `source_file`.
For `grid_gen`, increase generator resolution controls in `grid.generator.cli_args` or config.

@section bc_sec 2. Boundary Conditions

@subsection bc_simple_ssec 2.1 Set a Basic Inlet + Wall

```yaml
boundary_conditions:
  - face: "-Zeta"
    type: "INLET"
    handler: "constant_velocity"
    params: {vx: 0.0, vy: 0.0, vz: 1.5}
  - face: "-Eta"
    type: "WALL"
    handler: "noslip"
  # ... include all remaining faces ...
```

@subsection bc_periodic_ssec 2.2 Set a Periodic Direction

```yaml
models:
  domain:
    i_periodic: true
```

Use `PERIODIC` type on both paired faces (`-Xi` and `+Xi`) with supported handlers:
- `geometric`
- `constant_flux` (requires `target_flux`)

@section run_sec 3. Run and Monitoring

@subsection mpi_ssec 3.1 Run in Parallel and Set DMDA Layout

```bash
./scripts/pic.flow run -n 16 --solve ...
```

Optional global layout hint:

```yaml
grid:
  mode: programmatic_c
  programmatic_settings:
    da_processors_x: 4
    da_processors_y: 2
    da_processors_z: 2
```

`da_processors_*` are scalar globals, not per-block arrays.

@subsection cluster_run_ssec 3.2 Run on Slurm (Generate + Submit)

```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml
```

To only generate scripts/configs:

```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml \
  --no-submit
```

@subsection restart_ssec 3.3 Restart from a Saved Step

```yaml
run_control:
  start_step: 1000
  total_steps: 2000
```

This advances from step 1000 to 3000.

@subsection logging_ssec 3.4 Targeted Debug Logging

```yaml
logging:
  verbosity: "DEBUG"
  enabled_functions:
    - Projection
    - UpdatePressure
```

@section post_sec 4. Post-Processing

@subsection post_existing_ssec 4.1 Postprocess an Existing Run

```bash
./scripts/pic.flow run --post-process \
  --run-dir runs/flat_channel_20240401-153000 \
  --post my_study/standard_analysis.yml
```

@subsection qcrit_ssec 4.2 Add Q-Criterion

```yaml
eulerian_pipeline:
  - task: nodal_average
    input_field: Ucat
    output_field: Ucat_nodal
  - task: q_criterion

io:
  eulerian_fields:
    - Ucat_nodal
    - Qcrit
```

@subsection stats_ssec 4.3 Enable Statistics Output (MSD)

```yaml
statistics_pipeline:
  output_prefix: "Stats"
  tasks:
    - task: msd
```

@section sweep_sec 5. Sweeps and Study Arrays

```bash
python3 scripts/pic.flow sweep \
  --study my_study/study.yml \
  --cluster my_study/cluster.yml
```

This creates `studies/<study_id>/...`, submits solver/post array jobs with dependency chaining, and writes aggregated metrics/plots.

@section next_steps_sec 6. Next Steps

Proceed to **@subpage 12_Capabilities_Summary**.

For full option reference, also see:
- **@subpage 07_Case_Reference**
- **@subpage 08_Solver_Reference**
- **@subpage 09_Monitor_Reference**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 33_Initial_Conditions**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 38_Start_Here_10_Minutes**
- **@subpage 39_Common_Fatal_Errors**
