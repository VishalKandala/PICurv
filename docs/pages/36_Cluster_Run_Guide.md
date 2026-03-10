@page 36_Cluster_Run_Guide Cluster Run Guide (Slurm)

@anchor _Cluster_Run_Guide

This guide documents how `picurv` converts case configs into scheduler artifacts for cluster execution.

@tableofcontents

@section p36_inputs_sec 1. Required Inputs

Typical cluster-enabled run uses:

- `case.yml`
- `solver.yml`
- `monitor.yml`
- `post.yml`
- `cluster.yml` (scheduler contract)

Initialize templates from examples, then customize per cluster/account policy.

@section p36_command_sec 2. Core Command Patterns

Generate and submit:

```bash
./bin/picurv run \
  --solve --post-process \
  --case <case.yml> \
  --solver <solver.yml> \
  --monitor <monitor.yml> \
  --post <post.yml> \
  --cluster <cluster.yml>
```

Generate only (no submission):

```bash
./bin/picurv run ... --cluster <cluster.yml> --no-submit
```

Explicit scheduler cross-check (optional):

```bash
./bin/picurv run ... --cluster <cluster.yml> --scheduler slurm
```

Notes:

- `--scheduler` is optional and only valid with `--cluster`.
- if provided, it must match `cluster.yml:scheduler.type`.
- in cluster mode, solver stage rank count is resolved from `cluster.yml` (`nodes * ntasks_per_node`).
- post stage stays single-task by default (`nodes=1`, `ntasks_per_node=1`) even when solver stage is multi-rank.
- `cluster.yml` does not currently name the run directory. `picurv` generates `run_id` as `<case_basename>_<timestamp>`, then derives Slurm job names like `<run_id>_solve` and `<run_id>_post`.
- batch launcher precedence is: `cluster.yml.execution` first, then nearest `.picurv-execution.yml` (`cluster_execution`, then `default_execution`), then built-in `srun`.

@section p36_artifacts_sec 3. Generated Scheduler Artifacts

In run directory, scheduler generation typically produces:

- `runs/<run_id>/scheduler/solver.sbatch`
- `runs/<run_id>/scheduler/post.sbatch`
- `runs/<run_id>/manifest.json`
- `runs/<run_id>/scheduler/submission.json` (always in cluster mode when run artifacts are written; contains launch metadata and submission IDs when present)

These coexist with standard runtime control artifacts used by solver/postprocessor binaries.

@section p36_flow_sec 4. Submission Flow

1. YAML validation and contract checks,
2. run directory + control artifact generation,
3. sbatch script rendering from scheduler settings,
4. optional solver submission (`solver.sbatch`),
5. optional post submission (`post.sbatch`), with `afterok:<solver_jobid>` dependency when solve+post are both requested in one command.

This allows consistent local dry-run and cluster production flow from the same inputs.

@section p36_notes_sec 5. Operational Notes

- Prefer `--no-submit` first when validating new scheduler settings.
- Keep cluster defaults in reusable templates (`examples/master_template/master_cluster.yml`).
- For personal cluster copies, prefer local operational files such as `short_job.local.yml` / `long_job.local.yml` instead of committing account/module-specific scheduler profiles.
- If queue policies differ by partition/account, encode them in `cluster.yml` instead of editing generated scripts manually.
- Solver stage uses `cluster.yml` resources directly.
- Post stage defaults to single-task scheduling (`nodes=1`, `ntasks_per_node=1`) in generated `post.sbatch`.
- `picurv init` now creates `.picurv-execution.yml` in each new case with inert defaults.
- If your cluster needs the same MPI launcher tokens for login-node and batch runs, edit that file and let `cluster.yml` override only when batch jobs differ.
- For one-off interactive multi-rank runs on cluster login nodes, `PICURV_MPI_LAUNCHER` still overrides everything.
- `--num-procs` in cluster mode is a consistency guard, not an independent rank selector:
  - allowed values are `1` (auto) or exactly `nodes * ntasks_per_node`.
  - any other value is rejected as an inconsistent configuration combo.
- use `picurv run --dry-run --format json` before submission when integrating with external launch wrappers.

See also:

- **@subpage 37_Sweep_Studies_Guide**
- **@subpage 05_The_Conductor_Script**
- **@subpage 52_Run_Lifecycle_Guide**
- **@subpage 39_Common_Fatal_Errors**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Cluster Run Guide (Slurm)** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
