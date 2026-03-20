# Flat Channel Example Guide

This example is the default first-run baseline for PICurv using programmatic grid generation. It is the fastest path to validate a local install, understand config-role composition, and run the full solver-post workflow.

Because it avoids external grid files, this template is also useful for controlled solver tuning experiments where geometry complexity is not the primary variable.

## Included Files

- `flat_channel.yml`: case configuration.
- `Imp-MG-Standard.yml`: solver configuration.
- `Standard_Output.yml`: monitor/logging settings.
- `standard_analysis.yml`: postprocessing recipe.
- `slurm_cluster.yml`: sample cluster scheduler config.
- `execution.example.yml`: optional shared local/login-node + batch launcher defaults.
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

`init` creates the case directory with config files — no binaries are copied. Runtime executables are resolved from the project `bin/` directory via PATH. `init` records source linkage in `.picurv-origin.json` and writes an inert `.picurv-execution.yml` for later site-specific launcher overrides. To pin specific binary versions, run `picurv sync-binaries`.

If you run this example on a cluster and need site-specific MPI launcher tokens, edit the `.picurv-execution.yml` that `init` created in the case directory. Existing cases or repo-root site configs can still start from `execution.example.yml`. That same file can cover both login-node runs and generated batch jobs. Keep `slurm_cluster.yml` for scheduler policy and any batch-only override. Generated Slurm solver jobs now enable the runtime walltime guard by default; override it in `slurm_cluster.yml -> execution.walltime_guard` only when the default warmup/headroom policy is not a good fit.

## Recommended Uses

- first local smoke run,
- solver/post pipeline validation,
- baseline for parameter perturbation studies,
- reference behavior before introducing geometry complexity.

Run naming note:

- `picurv` will create `runs/flat_channel_<timestamp>/` automatically for this example.
- `slurm_cluster.yml` does not have a separate run-name field.

## Practical Debugging Tips

- Start with default tolerances and adjust one numerical control at a time.
- Keep monitor cadence high during early debugging, then reduce for long runs.
- Use this case as a control when comparing against failures in more complex geometries.

## Related Docs

- `flat_channel.md`
- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
