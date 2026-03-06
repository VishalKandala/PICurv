# Bent Channel Example Guide

This example demonstrates the file-based curvilinear-grid workflow (`grid.mode: file`) and is the baseline template for users who need to run on externally generated/curated grids.

Compared with programmatic-grid examples, this case is useful for understanding grid ingestion, metadata consistency, and curved-geometry boundary behavior in a realistic setup.

## Included Files

- `bent_channel.yml`, `Imp-MG-Standard.yml`, `Standard_Output.yml`, `standard_analysis.yml`
- `bent_channel_coarse.picgrid` and companion metadata files
- `slurm_cluster.yml` and `timestep_sensitivity_study.yml`

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

## Recommended Uses

- validating the file-grid ingestion path,
- checking BC behavior on curved geometry,
- testing timestep sensitivity studies,
- establishing a reference for imported production grids.

## Common Pitfalls

- Moving YAML files without moving associated `.picgrid`/metadata files.
- Changing decomposition or dimensions inconsistently with grid metadata.
- Interpreting pressure/projection instability as purely solver-tolerance issues when grid quality is the root cause.

## Related Docs

- `bent_channel.md`
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
