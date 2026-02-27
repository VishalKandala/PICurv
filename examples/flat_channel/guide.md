# Flat Channel Example Guide

This example is the baseline first-run case for PICurv using programmatic grid generation.

## Included Files

- `flat_channel.yml`: case configuration.
- `Imp-MG-Standard.yml`: solver configuration.
- `Standard_Output.yml`: monitor/logging settings.
- `standard_analysis.yml`: postprocessing recipe.
- `slurm_cluster.yml`: sample cluster scheduler config.
- `grid_independence_study.yml`: sample sweep study definition.

## Quick Start

```bash
./scripts/pic.flow init flat_channel --dest my_case
./scripts/pic.flow validate \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
./scripts/pic.flow run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

If project binaries are already built (`pic.flow build`), `init` also links/copies `picsolver` and `postprocessor` into the new case directory.

## Recommended Uses

- first local smoke run,
- solver/post pipeline validation,
- baseline for parameter perturbation studies.

## Related Docs

- `flat_channel.md`
- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
