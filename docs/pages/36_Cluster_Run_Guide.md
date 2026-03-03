@page 36_Cluster_Run_Guide Cluster Run Guide (Slurm)

This guide documents how `pic.flow` converts case configs into scheduler artifacts for cluster execution.

@tableofcontents

@section inputs_sec 1. Required Inputs

Typical cluster-enabled run uses:

- `case.yml`
- `solver.yml`
- `monitor.yml`
- `post.yml`
- `cluster.yml` (scheduler contract)

Initialize templates from examples, then customize per cluster/account policy.

@section command_sec 2. Core Command Patterns

Generate and submit:

```bash
./bin/pic.flow run \
  --case <case.yml> \
  --solver <solver.yml> \
  --monitor <monitor.yml> \
  --post <post.yml> \
  --cluster <cluster.yml> \
  --scheduler slurm
```

Generate only (no submission):

```bash
./bin/pic.flow run ... --scheduler slurm --no-submit
```

@section artifacts_sec 3. Generated Scheduler Artifacts

In run directory, scheduler generation typically produces:

- `solver.sbatch`
- `post.sbatch`
- `manifest.json`
- `submission.json` (when submitted)

These coexist with standard runtime control artifacts used by solver/postprocessor binaries.

@section flow_sec 4. Submission Flow

1. YAML validation and contract checks,
2. run directory + control artifact generation,
3. sbatch script rendering from scheduler settings,
4. optional solver submission,
5. optional post job with dependency linkage.

This allows consistent local dry-run and cluster production flow from the same inputs.

@section notes_sec 5. Operational Notes

- Prefer `--no-submit` first when validating new scheduler settings.
- Keep cluster defaults in reusable templates (`examples/master_template/master_cluster.yml`).
- If queue policies differ by partition/account, encode them in `cluster.yml` instead of editing generated scripts manually.
- Solver stage uses `cluster.yml` resources directly.
- Post stage defaults to single-task scheduling (`nodes=1`, `ntasks_per_node=1`) in generated `post.sbatch`.

See also:

- **@subpage 37_Sweep_Studies_Guide**
- **@subpage 05_The_Conductor_Script**
- **@subpage 39_Common_Fatal_Errors**
