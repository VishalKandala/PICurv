# Bent Channel Example Guide

This example demonstrates the file-based curvilinear-grid workflow (`grid.mode: file`) and is the baseline template for users who need to run on externally generated/curated grids.

Compared with programmatic-grid examples, this case is useful for understanding grid ingestion, metadata consistency, and curved-geometry boundary behavior in a realistic setup.

## Included Files

- `bent_channel.yml`, `Imp-MG-Standard.yml`, `Standard_Output.yml`, `standard_analysis.yml`
- `bent_channel_coarse.picgrid` and companion metadata files
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

This case expects `grid.source_file: bent_channel_coarse.picgrid`, which is bundled with the template.

If you run this example on a cluster and need site-specific MPI launcher tokens, copy `execution.example.yml` to `.picurv-execution.yml` in the repo root or case directory. That same file can drive both login-node runs and generated batch jobs. Keep `slurm_cluster.yml` reserved for scheduler policy and batch-only overrides.

## Recommended Uses

- validating the file-grid ingestion path,
- checking BC behavior on curved geometry,
- testing timestep sensitivity studies,
- establishing a reference for imported production grids.

Run naming note:

- `picurv` will create `runs/bent_channel_<timestamp>/` automatically for this example.
- `slurm_cluster.yml` does not provide a separate run-name key.

## Common Pitfalls

- Moving YAML files without moving associated `.picgrid`/metadata files.
- Changing decomposition or dimensions inconsistently with grid metadata.
- Interpreting pressure/projection instability as purely solver-tolerance issues when grid quality is the root cause.

## Related Docs

- `bent_channel.md`
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
