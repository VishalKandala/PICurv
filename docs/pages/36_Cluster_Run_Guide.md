@page 36_Cluster_Run_Guide Cluster Run Guide (Slurm)

This page documents single-case cluster execution using `pic.flow run`.

@tableofcontents

@section cluster_inputs_sec 1. Inputs

Required files:
- `case.yml`
- `solver.yml`
- `monitor.yml`
- `post.yml` (if `--post-process`)
- `cluster.yml`

Reference template:
- `examples/master_template/master_cluster.yml`

@section cluster_run_sec 2. Commands

Generate + submit:
```bash
python3 scripts/pic.flow run --solve --post-process \
  --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml \
  --cluster cluster.yml
```

Generate only:
```bash
python3 scripts/pic.flow run --solve --post-process \
  --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml \
  --cluster cluster.yml --no-submit
```

@section cluster_outputs_sec 3. Outputs

In `runs/<run_id>/scheduler/`:
- `solver.sbatch`
- `post.sbatch`
- `submission.json`

In `runs/<run_id>/`:
- `manifest.json`

`post.sbatch` is submitted with `afterok:<solver_jobid>` dependency when both stages are requested and submission is enabled.

@section cluster_notes_sec 4. Notes

- Slurm is the only scheduler supported in v1.
- Cluster mode uses `nodes * ntasks_per_node` from `cluster.yml` as effective MPI task count.
- `--no-submit` is useful for site review and manual `sbatch` execution.

