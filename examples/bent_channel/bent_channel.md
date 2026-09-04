# Case Template: Laminar Flow in a Bent Channel

## 1. Description

This template demonstrates flow through a curved channel using a file-based curvilinear grid.

Default setup:
- Flow type: laminar
- Reynolds number: 200
- Grid mode: `file`

## 2. Files in this Template

- `bent_channel.yml`: case definition (`grid.mode: file`)
- `Imp-MG-Standard.yml`: solver profile
- `Standard_Output.yml`: monitor profile
- `standard_analysis.yml`: postprocessor recipe
- `slurm_cluster.yml`: sample Slurm scheduler profile
- `execution.example.yml`: optional shared site launcher example
- `timestep_sensitivity_study.yml`: sample timestep sensitivity sweep

The grid is not a bundled file. `bent_channel.yml` selects `grid.mode: grid_gen` and points
at `config/grids/bent_channel_coarse.cfg`, which describes the geometry as a square
cross-section swept through a quarter turn. `picurv` builds the mesh, its `.vts` preview and
its `.info` quality report into the run's own asset store on every run, so they are current
rather than committed.

## 3. Typical Workflow

Initialize study:

```bash
./bin/picurv init bent_channel --dest my_bent_case
```

If cluster login-node or batch launches need site-specific MPI commands, edit `my_bent_case/.picurv-execution.yml` after `init`. That file is created with inert defaults, so it is safe to leave unchanged on ordinary local machines. Keep `slurm_cluster.yml` for scheduler settings and any batch-only override.

Run solve + post:

```bash
./bin/picurv run \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml \
  -n 4 --solve --post-process
```
`-n 4` runs both the solver and field postprocessor with four MPI ranks.

## 4. Changing Geometry

To use another grid:
1. Put new `.picgrid` in study directory.
2. Update `grid.source_file` in `bent_channel.yml`.
3. Re-run.

For generated grids, you can also switch case to `grid.mode: grid_gen` and configure `grid.generator`.
This makes `bent_channel.yml` a useful bridge case when comparing:

- fixed file-based geometry (`grid.mode: file`),
- generated file staging (`grid.mode: grid_gen`),
- shared solver/monitor/post profiles across both approaches.

## 5. Output Check

Open `<run.visualization>/standard_analysis/eulerian_data_*.vts` in ParaView and
inspect streamlines colored by `Ucat_nodal` magnitude.

## 6. Cluster and Sweep Examples

Single run on Slurm:

```bash
./bin/picurv run --solve --post-process \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml \
  --cluster my_bent_case/slurm_cluster.yml
```

Timestep sensitivity sweep:

```bash
./bin/picurv sweep \
  --study my_bent_case/timestep_sensitivity_study.yml \
  --cluster my_bent_case/slurm_cluster.yml
```

## 7. Live Docs

- https://vishalkandala.me/docs/picurv/03_Tutorial_File-Based_Grid.html
- https://vishalkandala.me/docs/picurv/36_Cluster_Run_Guide.html
- https://vishalkandala.me/docs/picurv/37_Sweep_Studies_Guide.html
- https://vishalkandala.me/docs/picurv/48_Grid_Generator_Guide.html
