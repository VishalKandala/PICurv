# Case Template: Laminar Flow in a Flat Channel

## 1. Description

This template runs laminar incompressible channel flow on a programmatically generated grid.

Default setup:
- Flow type: laminar
- Reynolds number: 200 (from physical scaling inputs)
- Grid mode: `programmatic_c`

## 2. Files in this Template

- `flat_channel.yml`: case definition (physics, grid, BCs, run control)
- `Imp-MG-Standard.yml`: solver profile
- `Standard_Output.yml`: monitor profile
- `standard_analysis.yml`: postprocessor recipe
- `slurm_cluster.yml`: sample Slurm scheduler profile
- `grid_independence_study.yml`: sample sweep definition

## 3. Typical Workflow

Initialize study:

```bash
./scripts/pic.flow init flat_channel --dest my_channel_case
```

Run solve + post:

```bash
./scripts/pic.flow run \
  --case my_channel_case/flat_channel.yml \
  --solver my_channel_case/Imp-MG-Standard.yml \
  --monitor my_channel_case/Standard_Output.yml \
  --post my_channel_case/standard_analysis.yml \
  -n 4 --solve --post-process
```
`-n 4` sizes the solver stage. Post-processing defaults to one rank.

## 4. Common Edits

- Reynolds number: `properties.scaling` and `properties.fluid`
- Runtime: `run_control`
- Grid resolution: `grid.programmatic_settings.im/jm/km`
- Momentum strategy: `solver.strategy.momentum_solver`

## 5. Cluster and Sweep Examples

Single run on Slurm:

```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_channel_case/flat_channel.yml \
  --solver my_channel_case/Imp-MG-Standard.yml \
  --monitor my_channel_case/Standard_Output.yml \
  --post my_channel_case/standard_analysis.yml \
  --cluster my_channel_case/slurm_cluster.yml
```

Grid-independence sweep:

```bash
python3 scripts/pic.flow sweep \
  --study my_channel_case/grid_independence_study.yml \
  --cluster my_channel_case/slurm_cluster.yml
```

## 6. Output Check

Open `runs/<run_id>/viz/Field_*.vts` in ParaView and color by `Ucat_nodal` to inspect channel profile development.
