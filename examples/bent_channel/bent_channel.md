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
- `timestep_sensitivity_study.yml`: sample timestep sensitivity sweep
- `bent_channel_coarse.picgrid`: solver input grid
- `bent_channel_coarse.vts`: grid visualization helper
- `bent_channel_coarse.info`: grid metrics summary

## 3. Typical Workflow

Initialize study:

```bash
./scripts/pic.flow init bent_channel --dest my_bent_case
```

Run solve + post:

```bash
./scripts/pic.flow run \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml \
  -n 4 --solve --post-process
```
`-n 4` sizes the solver stage. Post-processing defaults to one rank.

## 4. Changing Geometry

To use another grid:
1. Put new `.picgrid` in study directory.
2. Update `grid.source_file` in `bent_channel.yml`.
3. Re-run.

For generated grids, you can also switch case to `grid.mode: grid_gen` and configure `grid.generator`.

## 5. Output Check

Open `runs/<run_id>/viz/Field_*.vts` in ParaView and inspect streamlines colored by `Ucat_nodal` magnitude.

## 6. Cluster and Sweep Examples

Single run on Slurm:

```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml \
  --cluster my_bent_case/slurm_cluster.yml
```

Timestep sensitivity sweep:

```bash
python3 scripts/pic.flow sweep \
  --study my_bent_case/timestep_sensitivity_study.yml \
  --cluster my_bent_case/slurm_cluster.yml
```
