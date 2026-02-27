# Bent Channel Example Guide

This example demonstrates file-based curvilinear-grid workflow (`grid.mode: file`).

## Included Files

- `bent_channel.yml`, `Imp-MG-Standard.yml`, `Standard_Output.yml`, `standard_analysis.yml`
- `bent_channel_coarse.picgrid` and companion metadata files
- `slurm_cluster.yml` and `timestep_sensitivity_study.yml`

## Quick Start

```bash
./scripts/pic.flow init bent_channel --dest my_bent_case
./scripts/pic.flow validate \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml
./scripts/pic.flow run --solve --post-process -n 4 \
  --case my_bent_case/bent_channel.yml \
  --solver my_bent_case/Imp-MG-Standard.yml \
  --monitor my_bent_case/Standard_Output.yml \
  --post my_bent_case/standard_analysis.yml
```

This case expects `grid.source_file: bent_channel_coarse.picgrid`, which is bundled in the template.

## Recommended Uses

- validating file-grid ingestion path,
- checking BC behavior on curved geometry,
- testing timestep sensitivity studies.

## Related Docs

- `bent_channel.md`
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
