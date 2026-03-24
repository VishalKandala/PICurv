# PICurv Post-Processing Profiles Library

This directory stores reusable `post.yml` analysis recipes for `./bin/picurv run --post-process`. A post profile defines what data is read, what transforms/statistics are computed, and what artifacts are written for inspection or downstream analysis.

## 1. What A Post Profile Controls

- timestep window selection (`run_control`),
- input-source location and extension mapping (`source_data`, `io.input_extensions`),
- Eulerian field operations,
- Lagrangian particle operations,
- statistics pipeline tasks,
- output naming and directory policy.

## 2. Included Profiles

- `standard_analysis.yml`: baseline Eulerian analysis (dimensionalization, nodal fields, Q-criterion).

For full schema coverage, use:
- `examples/master_template/master_postprocessor.yml`

## 3. How To Use

Postprocess an existing run:

```bash
./bin/picurv run --post-process \
  --run-dir runs/<run_id> \
  --post config/postprocessors/standard_analysis.yml
```

Catch up the same post profile on an existing run without editing `run_control.start_step`:

```bash
./bin/picurv run --post-process --continue \
  --run-dir runs/<run_id> \
  --post config/postprocessors/standard_analysis.yml
```

Solve and postprocess in one command:

```bash
./bin/picurv run --solve --post-process -n 8 \
  --case my_study/case.yml \
  --solver my_study/solver.yml \
  --monitor my_study/monitor.yml \
  --post config/postprocessors/standard_analysis.yml
```

`-n/--num-procs` applies to solver execution. Postprocessing defaults to single-rank execution.

## 4. Notes on Newer Options

- `statistics_pipeline` is supported (canonical task: `msd`).
- `io.input_extensions.eulerian/particle` map reader expectations to generated filenames.
- `io.eulerian_fields_averaged` is accepted as reserved passthrough.
- keep the full desired timestep window in `run_control` for a reusable post profile. With `--continue`, PICurv computes the effective restart step for the same recipe lineage instead of requiring you to edit `start_step` by hand.
- if the solver is still writing outputs, PICurv processes only the highest fully available contiguous source prefix for the current recipe.
- PICurv enforces a single post writer per run directory so two post jobs cannot race on the same `viz/` or `statistics/` outputs.

## 5. Validation Tips

- Start with a small timestep window when adding a new post operation.
- Validate output file count and naming before long batch processing.
- Keep post profiles under version control beside study definitions for reproducibility.
