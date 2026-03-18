@page 37_Sweep_Studies_Guide Sweep and Study Guide

@anchor _Sweep_Studies_Guide

`picurv sweep` orchestrates parameter studies with generated run variants, scheduler arrays, and aggregate metrics.

@tableofcontents

@section p37_inputs_sec 1. Inputs and Templates

A sweep/study commonly uses:

- base case/solver/monitor/post templates,
- `study.yml` defining parameter combinations and metrics,
- optional cluster scheduler settings for array submission.

Starter templates are available under `examples/*/*study*.yml` and `examples/master_template/`.

@section p37_command_sec 2. Core Sweep Command

```bash
./bin/picurv sweep --study <study.yml> --cluster <cluster.yml>
```

Optional generation-only mode:

```bash
./bin/picurv sweep --study <study.yml> --cluster <cluster.yml> --no-submit
```

Delayed submit from existing staged study artifacts:

```bash
./bin/picurv submit --study-dir studies/<study_id>
```

There is no dedicated `--dry-run` flag on `sweep`; use `--no-submit` for non-submitting artifact generation.

@section p37_contract_sec 3. Study Contract Essentials

A study definition usually specifies:

- `base_configs`:
  - `case`, `solver`, `monitor`, `post` paths (all required)
- `study_type`:
  - one of `grid_independence`, `timestep_independence`, `sensitivity`
- `parameters`:
  - non-empty mapping of `<target>.<yaml.path>` -> non-empty list of values
  - `<target>` must be one of `case`, `solver`, `monitor`, `post`
- `metrics` (optional):
  - list of metric specs or metric names for aggregation
- `plotting` (optional):
  - output controls (`enabled`, `output_format`)
- `execution` (optional):
  - controls like `max_concurrent_array_tasks` for Slurm array throttling

Each combination yields a generated run with fully materialized config set.

Parameter keys can target nested case/solver/monitor/post values such as:

- `case.models.physics.particles.count`
- `case.run_control.dt_physical`

Not every study should use the default `msd_final` metric shorthand. Cases that
write other scalar diagnostics, such as `logs/interpolation_error.csv`, should
define explicit CSV metrics instead.

@section p37_outputs_sec 4. Outputs and Aggregates

Expected study outputs include:

- `studies/<study_id>/cases/case_####/` per-combination run directories
- `studies/<study_id>/scheduler/case_index.tsv`
- `studies/<study_id>/scheduler/solver_array.sbatch`
- `studies/<study_id>/scheduler/post_array.sbatch`
- `studies/<study_id>/scheduler/metrics_aggregate.sbatch`
- `studies/<study_id>/scheduler/solver_<array_jobid>_<taskid>.out/.err` after submission
- `studies/<study_id>/scheduler/post_<array_jobid>_<taskid>.out/.err` after submission
- `studies/<study_id>/scheduler/submission.json` (when jobs are submitted)
- `studies/<study_id>/results/metrics_table.csv`
- `studies/<study_id>/results/plots/*` (when plotting is enabled and matplotlib is available)
- `studies/<study_id>/study_manifest.json`

This keeps raw run data and comparative study diagnostics in one reproducible structure.

Metrics aggregation runs automatically as a Slurm job chained after the post-processing
array (`afterany` dependency). If the automatic metrics job fails (e.g. Python unavailable
on compute nodes), use `--reaggregate` manually.

@section p37_operations_sec 5. Operational Workflow

Recommended workflow:

1. run a tiny subset locally or with `--no-submit`,
2. verify parameter substitution and metric extraction,
3. launch full array, either directly with `picurv sweep ...` or later with `picurv submit --study-dir ...`,
4. inspect aggregate outputs (auto-collected by the metrics Slurm job, or via `--reaggregate`),
5. archive the exact study file with results for reproducibility.

`picurv sweep` is the scheduler-backed study path. For local parameter studies,
repeat `picurv run` manually across a small set of edited case variants and
compare the resulting run directories.

For fragile metrics, add smoke tests or fixture-based validation before large queue submissions.

Implementation details worth knowing:

- case expansion uses cartesian product over all `parameters.*` lists.
- generated case configs are revalidated through the same solver/post validators used by `picurv run`.
- submission chain: solver array → post array (`afterok`) → metrics job (`afterany`).
- `scheduler/submission.json` is the study-directory contract consumed by `picurv submit --study-dir ...`.
- generator/file grid external paths are rewritten to absolute paths during case materialization so they remain valid in `studies/<study_id>/cases/...`.
- generated `solver_array.sbatch` exports walltime metadata for the runtime walltime guard, while `post_array.sbatch` remains a plain post-processing launcher.

@section p37_continue_sec 6. Continuing a Partially-Completed Study

If any solver case is killed (e.g. by the walltime guard or Slurm time limit),
the entire post array is cancelled (`afterok` dependency). Use `--continue` to
resume the study:

```bash
./bin/picurv sweep --continue --study-dir studies/<study_id>
```

To override cluster resources (e.g. increase walltime):

```bash
./bin/picurv sweep --continue --study-dir studies/<study_id> \
  --cluster cluster_more_time.yml
```

What `--continue` does:

1. Reads the original `study.yml` and `case_index.tsv` from the study directory.
2. Classifies each case as **complete**, **partial**, or **empty** by scanning checkpoints.
3. If all cases are complete, auto-aggregates metrics and exits (no jobs submitted).
4. For partial cases: updates `case.yml` (`start_step`, `total_steps`), sets particle
   `restart_mode` to `load` when a checkpoint exists, and delegates to `resolve_restart_source`
   for the full restart scenario matrix.
5. For empty cases (no checkpoint): re-runs from scratch with unmodified control files.
6. Submits a sparse solver array (incomplete cases only) → full post array → metrics aggregation.

Repeated continuation is safe: the target step count is always computed from the
original `study.yml`, not from the (potentially modified) per-case `case.yml`.

@section p37_reaggregate_sec 7. Manual Metrics Re-Aggregation

If the automatic metrics Slurm job fails or you want to re-collect metrics after
manual intervention:

```bash
./bin/picurv sweep --reaggregate --study-dir studies/<study_id>
```

This reads all case outputs, writes `results/metrics_table.csv`, and generates
plots (if enabled in `study.yml`).

@section p37_refs_sec 8. Related Pages

- **@subpage 36_Cluster_Run_Guide**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 40_Testing_and_Quality_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Sweep and Study Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
