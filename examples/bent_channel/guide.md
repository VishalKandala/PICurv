# Bent Channel Example Guide

This example demonstrates a generated curvilinear grid (`grid.mode: grid_gen`, `sweep`
geometry): a square cross-section swept through a quarter turn. For the file-based
curvilinear-grid workflow (`grid.mode: file`) instead, see `examples/search_robustness`,
which carries the one `.picgrid` the examples still ship - including this same bent-duct
geometry, for direct comparison.

Compared with programmatic-grid examples, this case is useful for understanding generated
curvilinear geometry, metadata consistency, and curved-geometry boundary behavior in a
realistic setup.

## Included Files

- `bent_channel.yml`, `Imp-MG-Standard.yml`, `Standard_Output.yml`, `standard_analysis.yml`
- the grid recipe `config/grids/bent_channel_coarse.cfg`, built at run time
- `slurm_cluster.yml`, `execution.example.yml`, and `timestep_sensitivity_study.yml`

## Quick Start

```bash
./bin/picurv init bent_channel --dest my_bent_case
./bin/picurv validate \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml
./bin/picurv run --solve --post-process -n 4 \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml
```

Analytical uniform-flow variant: `config/solvers/Analytical-UniformFlow.yml` needs
`grid.mode: file` or `programmatic_c` - `grid_gen` is refused for any analytical carrier
other than `TGV3D` (see `docs/pages/14_Config_Contract.md`). To pair it with this case's bent
geometry, override `grid` in `bent_channel.yml` to `mode: file` pointing at a copy of
`examples/search_robustness/bent_channel_coarse.picgrid` placed in the case directory
(`source_file` resolves relative to the case root). `tests/smoke/run_smoke.sh`'s
`bent_analytical_uniform` workflow does exactly this and is the reference for the exact
steps.

If you run this example on a cluster and need site-specific MPI launcher tokens, edit the `.picurv-execution.yml` that `init` created in the case directory. Existing cases or repo-root site configs can still start from `execution.example.yml`. That same file can drive both login-node runs and generated batch jobs. Keep `slurm_cluster.yml` reserved for scheduler policy and batch-only overrides. Generated Slurm solver jobs now enable the runtime walltime guard by default; override it in `slurm_cluster.yml -> execution.walltime_guard` only when needed.

## Recommended Uses

- validating the `sweep` generator geometry and its config grammar,
- checking BC behavior on curved geometry,
- testing timestep sensitivity studies,
- a starting point for other bent-duct or multi-segment centreline studies.

Run naming note:

- `picurv` will create `runs/bent_channel_<timestamp>/` automatically for this example.
- `slurm_cluster.yml` does not provide a separate run-name key.

## Common Pitfalls

- Editing `config/grids/bent_channel_coarse.cfg`'s cell counts without checking the
  multigrid ladder (`Imp-MG-Standard.yml`'s levels need node counts of the form
  `2^(L-1)*k + 1`); `mg_levels` in the config makes the generator refuse a bad count
  instead of the solver silently coarsening less than requested.
- Shortening the bend's `arc` segment without checking `Max_Non_Orthogonality_deg` in the
  generated `.info` report - a tighter turn on the same cell count costs orthogonality.
- Interpreting pressure/projection instability as purely solver-tolerance issues when grid
  quality is the root cause.

## Related Docs

- `bent_channel.md`
- https://vishalkandala.me/docs/picurv/48_Grid_Generator_Guide.html
- https://vishalkandala.me/docs/picurv/03_Tutorial_File-Based_Grid.html (`grid.mode: file`,
  demonstrated by `examples/search_robustness`)
