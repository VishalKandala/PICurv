# PICurv Post-Processing Profiles Library

## 1. Purpose

This directory stores reusable `post.yml` analysis recipes for `./bin/picurv run --post-process`.

A post profile controls:
- timestep range (`run_control`),
- source location (`source_data`),
- Eulerian/Lagrangian/statistics pipelines,
- output field selection and file naming (`io`).

## 2. Included Profiles

- `standard_analysis.yml`: baseline Eulerian analysis (dimensionalize, nodal fields, Q-criterion).

For full schema coverage, use the master reference:
- `examples/master_template/master_postprocessor.yml`

## 3. How to Use

Postprocess an existing run:

```bash
./bin/picurv run --post-process \
  --run-dir runs/<run_id> \
  --post config/postprocessors/standard_analysis.yml
```

Solve + post in one command:

```bash
./bin/picurv run --solve --post-process -n 8 \
  --case my_study/case.yml \
  --solver my_study/solver.yml \
  --monitor my_study/monitor.yml \
  --post config/postprocessors/standard_analysis.yml
```
`-n/--num-procs` applies to the solver stage. Post-processing defaults to single-rank execution.

## 4. Notes on Newer Options

- `statistics_pipeline` is supported (currently the canonical `msd` task).
- `io.input_extensions.eulerian/particle` map to post reader extensions.
- `io.eulerian_fields_averaged` is accepted as reserved passthrough.
