@page 37_Sweep_Studies_Guide Sweep and Study Guide

This page documents Slurm array studies using `pic.flow sweep`.

@tableofcontents

@section sweep_inputs_sec 1. Inputs

Required:
- `study.yml`
- `cluster.yml`

Reference templates:
- `examples/master_template/master_study.yml`
- `examples/master_template/master_cluster.yml`

@section sweep_cmd_sec 2. Command

```bash
python3 scripts/pic.flow sweep \
  --study study.yml \
  --cluster cluster.yml
```

Generate only:

```bash
python3 scripts/pic.flow sweep \
  --study study.yml \
  --cluster cluster.yml \
  --no-submit
```

@section sweep_contract_sec 3. Study Contract

- `base_configs` defines base case/solver/monitor/post files.
- `parameters` expands as a cartesian matrix.
- keys use `<target>.<yaml.path>` where target is `case|solver|monitor|post`.
- `metrics` defines file-based extractors (CSV/log).
- `plotting` controls plot generation.
- `execution.max_concurrent_array_tasks` maps to Slurm `%N` array throttling.

@section sweep_outputs_sec 4. Outputs

`studies/<study_id>/` contains:
- `cases/<case_i>/...` materialized run trees
- `scheduler/case_index.tsv`
- `scheduler/solver_array.sbatch`
- `scheduler/post_array.sbatch`
- `scheduler/submission.json`
- `results/metrics_table.csv`
- `results/plots/*` (when plottable metrics are available)
- `study_manifest.json`

@section sweep_notes_sec 5. Notes

- Post array is submitted with `afterok:<solver_array_jobid>`.
- Metrics/plots rely on existing post/stat outputs (no C-side changes required in v1).
- If files are not available yet at generation time, aggregation may produce sparse metrics.

