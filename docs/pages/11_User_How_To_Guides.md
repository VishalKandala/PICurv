@page 11_User_How_To_Guides User How-To Guides

This page provides operational recipes for common PICurv tasks.
Each recipe includes what to change, why it matters, and a quick verification action.

@tableofcontents

@section setup_sec 1. Setup and Physics

@subsection reynolds_ssec 1.1 Change Reynolds Number

What to change (in `case.yml`):

```yaml
properties:
  scaling:
    length_ref: 0.1
    velocity_ref: 1.5
  fluid:
    density: 1000.0
    viscosity: 0.001
```

Why:

\f[
Re = \frac{\rho U L}{\mu}
\f]

These values set the non-dimensional operating point consumed by solver controls.

Quick check:

- rerun `pic.flow validate ...`,
- inspect generated `.control` file for updated values.

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

Why:

- 2D mode still uses a thin third dimension for structured-grid machinery,
- small `km` is typically enough for planar scenarios.

Quick check:

- confirm generated run uses expected Z resolution.

@subsection gridres_ssec 1.3 Increase Grid Resolution

- `programmatic_c`: increase `im/jm/km` arrays,
- `file`: use finer `.picgrid`,
- `grid_gen`: increase generator resolution args.

Verification:

- compare runtime memory/cost and key output metrics across resolutions.

@section bc_sec 2. Boundary Conditions

@subsection bc_simple_ssec 2.1 Set a Constant-Velocity Inlet and Walls

```yaml
boundary_conditions:
  - face: "-Zeta"
    type: "INLET"
    handler: "constant_velocity"
    params: {vx: 0.0, vy: 0.0, vz: 1.5}
  - face: "-Eta"
    type: "WALL"
    handler: "noslip"
  # define all remaining faces explicitly
```

Why:

- handler/type compatibility is validated,
- all faces must be covered for each block.

Verification:

- use `validate` first and check BC generation files under `runs/<run_id>/config/`.

@subsection bc_periodic_ssec 2.2 Enable Periodicity in One Direction

```yaml
models:
  domain:
    i_periodic: true
```

Use `PERIODIC` BC type on both paired faces with supported handlers:

- `geometric`
- `constant_flux` (requires `target_flux`)

Verification:

- confirm paired-face consistency in validation output.

@section run_sec 3. Running and Monitoring

@subsection mpi_ssec 3.1 Run in Parallel and Control DMDA Layout

```bash
./bin/pic.flow run -n 16 --solve ...
```

Optional partition hints:

```yaml
grid:
  mode: programmatic_c
  programmatic_settings:
    da_processors_x: 4
    da_processors_y: 2
    da_processors_z: 2
```

Note: `da_processors_*` are scalar globals, not per-block vectors.

@subsection cluster_run_ssec 3.2 Run on Slurm (Generate and Submit)

```bash
./bin/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml
```

Generate-only mode:

```bash
./bin/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml \
  --no-submit
```

Verification:

- inspect `scheduler/*.sbatch` and `submission.json` in run directory.

@subsection restart_ssec 3.3 Restart from a Saved Step

```yaml
run_control:
  start_step: 1000
  total_steps: 2000
```

Meaning:

- run resumes from step 1000 and advances 2000 additional steps.

Verification:

- confirm restart directory and step indices in run logs.

@subsection logging_ssec 3.4 Enable Targeted Debug Logging

```yaml
logging:
  verbosity: "DEBUG"
  enabled_functions:
    - Projection
    - UpdatePressure
```

Use this for local diagnosis of instability or boundary anomalies.
Prefer narrow function lists to keep logs manageable.

@section post_sec 4. Post-Processing Recipes

@subsection post_existing_ssec 4.1 Postprocess an Existing Run

```bash
./bin/pic.flow run --post-process \
  --run-dir runs/flat_channel_20240401-153000 \
  --post my_study/standard_analysis.yml
```

Use when solver outputs already exist and you are iterating only on analysis pipeline.

@subsection qcrit_ssec 4.2 Add Q-Criterion to Eulerian Pipeline

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

Verification:

- open VTK output and confirm `Qcrit` field is present.

@subsection stats_ssec 4.3 Enable Statistics Output (MSD)

```yaml
statistics_pipeline:
  output_prefix: "Stats"
  tasks:
    - task: msd
```

Verification:

- check `Stats_msd.csv` in configured output location.

@section sweep_sec 5. Sweep Studies

```bash
./bin/pic.flow sweep \
  --study my_study/study.yml \
  --cluster my_study/cluster.yml
```

What you get:

- expanded case matrix,
- scheduler array scripts,
- aggregated metrics table,
- optional plots.

See **@subpage 37_Sweep_Studies_Guide** for full contract details.

@section next_steps_sec 6. Next Steps

- Full capability map: **@subpage 12_Capabilities_Summary**
- Config references: **@subpage 07_Case_Reference**, **@subpage 08_Solver_Reference**, **@subpage 09_Monitor_Reference**, **@subpage 10_Post_Processing_Reference**
- Troubleshooting and quality: **@subpage 39_Common_Fatal_Errors**, **@subpage 40_Testing_and_Quality_Guide**
