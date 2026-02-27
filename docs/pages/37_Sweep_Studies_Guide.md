@page 37_Sweep_Studies_Guide Sweep and Study Guide

`pic.flow sweep` orchestrates parameter studies with generated run variants, scheduler arrays, and aggregate metrics.

@tableofcontents

@section inputs_sec 1. Inputs and Templates

A sweep/study commonly uses:

- base case/solver/monitor/post templates,
- `study.yml` defining parameter combinations and metrics,
- optional cluster scheduler settings for array submission.

Starter templates are available under `examples/*/*study*.yml` and `examples/master_template/`.

@section command_sec 2. Core Sweep Command

```bash
./scripts/pic.flow sweep --study <study.yml> --scheduler slurm
```

Typical optional flags include dry-run or no-submit behaviors depending on workflow stage.

@section contract_sec 3. Study Contract Essentials

A study definition usually specifies:

- parameter space (explicit lists/ranges),
- mapping from parameter key -> YAML target path,
- execution controls,
- metric extraction definitions.

Each combination yields a generated run with fully materialized config set.

@section outputs_sec 4. Outputs and Aggregates

Expected study outputs include:

- per-combination run directories,
- aggregate table (for example `metrics_table.csv`),
- summary metadata (for example `summary.json`),
- optional plot outputs in results/plots.

This keeps raw run data and comparative study diagnostics in one reproducible structure.

@section operations_sec 5. Operational Workflow

Recommended workflow:

1. run a tiny subset locally or with `--no-submit`,
2. verify parameter substitution and metric extraction,
3. launch full array,
4. inspect aggregate outputs,
5. archive the exact study file with results for reproducibility.

For fragile metrics, add smoke tests or fixture-based validation before large queue submissions.

@section refs_sec 6. Related Pages

- **@subpage 36_Cluster_Run_Guide**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 40_Testing_and_Quality_Guide**
