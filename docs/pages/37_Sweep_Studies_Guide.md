@page 37_Sweep_Studies_Guide Sweep and Study Guide

@anchor _Sweep_Studies_Guide

`picurv sweep` orchestrates parameter studies with generated run variants, scheduler arrays, and aggregate metrics.

For protecting or offloading a whole study or selected `case_####` members, see **@subpage 61_Storage_Management_Guide**.

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
- `parameter_sets`:
  - optional alternative to `parameters` for coupled overrides that should move together
  - non-empty list of explicit `<target>.<yaml.path>` -> scalar bundles
  - provide exactly one of `parameters` or `parameter_sets`
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
write other scalar diagnostics, such as `<run.runtime_logs>/interpolation_error.csv`, should
define explicit CSV metrics instead. Search and migration characterization
studies can aggregate `<run.runtime_logs>/search_metrics.csv` columns such as
`search_failure_fraction`, `search_work_index`, `re_search_fraction`, or
normalized run-level signals derived from `lost_cumulative`.

CSV metric specs also support `reduction: p95`, per-row ratios via
`numerator_column` plus `denominator_column`, and scalar normalization through
`normalize_by_parameter` for observables such as run loss fraction.

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

- case expansion uses the full cross-product over all `parameters.*` lists, or explicit paired bundles from `parameter_sets` when coupled overrides must move together.
- within `parameters`, mapping a parent key to a list of dicts sweeps multiple sub-fields together as a unit — they are NOT cartesian-producted with each other, but ARE cross-producted against other parameters. Example: `case.run_control: [{total_steps: 400, dt_physical: 0.005}, {total_steps: 800, dt_physical: 0.0025}]` sweeps the whole `run_control` block as two groups.
- generated case configs are revalidated through the same solver/post validators used by `picurv run`.
- submission chain: solver array → post array (`afterok`) → metrics job (`afterany`).
- `<run.scheduler>/submission.json` is the study-directory contract consumed by `picurv submit --study-dir ...`.
- generator/file grid external paths are rewritten to absolute paths during case materialization so they remain valid in `studies/<study_id>/cases/...`.
- generated `solver_array.sbatch` exports walltime metadata for the runtime walltime guard, while `post_array.sbatch` remains a plain post-processing launcher.
- `post_array.sbatch` uses the same `nodes * ntasks_per_node` allocation as the solver array; conflicting launcher `-n`/`-np` flags are rewritten to that task count.

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

@section p37_cap_study_sec 7.1 Study Type Entries

@htmlinclude generated/capability_inventory_study_type.html

@subsection p37_cap_study_grid_independence_sub grid_independence

@anchor p37_cap_study_grid_independence

**Identity.** `study.yml -> study_type: grid_independence`.

**What it does.** Expands a cross-product over grid parameters and aggregates the declared metrics against resolution, so the study answers whether the result has stopped changing with the mesh.

**When to choose it.** Before trusting any quantitative result from a new geometry. It is the first study to run and the one most often skipped; a converged solver on an unconverged grid is still the wrong answer.

**Parameters it owns.** `parameters` naming the grid keys to vary - typically `case.grid.programmatic_settings.im/jm/km` - and `metrics` naming what to compare.

**Interactions.** The sweep recognises the grid keys and orders the aggregated table by resolution. Combines naturally with `grid.mode: grid_gen`, which regenerates the mesh per case; with `mode: file` every resolution needs its own staged mesh.

**Diagnostics.** `results/metrics_table.csv` carries one row per case. A metric that does not settle across the ladder is the finding, not a failure of the study.

**Evidence.** Production exercised - `examples/flat_channel/grid_independence_study.yml`.

**Limitations.** It varies what you list and nothing else: refining only one direction produces a table that looks converged while the unrefined direction still controls the error.

@subsection p37_cap_study_timestep_independence_sub timestep_independence

@anchor p37_cap_study_timestep_independence

**Identity.** `study.yml -> study_type: timestep_independence`.

**What it does.** Expands a cross-product over timestep parameters and aggregates metrics against dt, separating temporal discretisation error from spatial.

**When to choose it.** After grid independence, and before any claim about a time-dependent quantity. Run it separately from the grid ladder: varying both at once cannot attribute the change to either.

**Parameters it owns.** `parameters` naming the timestep keys, and `metrics`.

**Interactions.** Independent of the momentum solver chosen, but the useful dt range is not: an explicit path is stability-limited where the implicit paths are not, so the ladder that is informative for one may not run at all for the other.

**Diagnostics.** Same `results/metrics_table.csv` shape as the grid ladder, ordered by dt.

**Evidence.** Implemented only. No shipped study selects it.

**Limitations.** Nothing here detects that the ladder crossed a stability boundary rather than a convergence one; a diverged case appears as a missing or absurd metric row.

@subsection p37_cap_study_sensitivity_sub sensitivity

@anchor p37_cap_study_sensitivity

**Identity.** `study.yml -> study_type: sensitivity`.

**What it does.** Expands either a cross-product or explicit `parameter_sets` over arbitrary configuration keys and aggregates the declared metrics, without assuming the varied parameter is a discretisation control.

**When to choose it.** For everything that is not a convergence ladder: Reynolds number, model coefficients, boundary treatments, solver tolerances. Use `parameter_sets` rather than `parameters` when the combinations are coupled and a full cross-product would include cases that make no sense.

**Parameters it owns.** `parameters` (cross-product) or `parameter_sets` (explicit coupled sets), and `metrics`.

**Interactions.** The only study type that does not order its table by a convergence quantity, so the aggregation is a plain comparison. `execution.max_concurrent_array_tasks` matters most here, because a sensitivity sweep is usually the widest.

**Diagnostics.** `results/metrics_table.csv` with one row per case and the varied keys as columns.

**Evidence.** Implemented only. No shipped study selects it.

**Limitations.** It compares what you asked for and infers nothing: there is no sensitivity index, no ranking, and no detection that two varied parameters interact.

@section p37_refs_sec 8. Related Pages

- **@subpage 36_Cluster_Run_Guide**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 40_Testing_and_Quality_Guide**
