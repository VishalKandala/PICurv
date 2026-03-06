# Flat Channel Example Guide

This example is the default first-run baseline for PICurv using programmatic grid generation. It is the fastest path to validate a local install, understand config-role composition, and run the full solver-post workflow.

Because it avoids external grid files, this template is also useful for controlled solver tuning experiments where geometry complexity is not the primary variable.

## Included Files

- `flat_channel.yml`: case configuration.
- `Imp-MG-Standard.yml`: solver configuration.
- `Standard_Output.yml`: monitor/logging settings.
- `standard_analysis.yml`: postprocessing recipe.
- `slurm_cluster.yml`: sample cluster scheduler config.
- `grid_independence_study.yml`: sample sweep study definition.

## Quick Start

```bash
./bin/picurv init flat_channel --dest my_case
./bin/picurv validate \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
./bin/picurv run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

If project binaries are already built, `init` copies executables (`picurv`, `simulator`, `postprocessor`) into the case directory and records source linkage in `.picurv-origin.json` for later maintenance commands.

## Recommended Uses

- first local smoke run,
- solver/post pipeline validation,
- baseline for parameter perturbation studies,
- reference behavior before introducing geometry complexity.

## Practical Debugging Tips

- Start with default tolerances and adjust one numerical control at a time.
- Keep monitor cadence high during early debugging, then reduce for long runs.
- Use this case as a control when comparing against failures in more complex geometries.

## Related Docs

- `flat_channel.md`
- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
