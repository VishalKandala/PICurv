@page 05_The_Conductor_Script The Conductor Script: pic.flow

`pic.flow` is the workflow orchestrator for PICurv. It validates YAML inputs, generates C runtime artifacts, and runs or schedules solver/postprocessing stages.

@tableofcontents

@section usage_sec 1. General Usage

```bash
python3 scripts/pic.flow [COMMAND] [ARGS...]
```

Primary commands:
- `init`
- `build`
- `run`
- `sweep`
- `validate`

Help:
```bash
python3 scripts/pic.flow --help
python3 scripts/pic.flow run --help
python3 scripts/pic.flow sweep --help
python3 scripts/pic.flow validate --help
```

@section run_sec 2. run: Single-Case Workflow

```bash
python3 scripts/pic.flow run [STAGES] [INPUTS] [OPTIONS]
```

Stages:
- `--solve`
- `--post-process`

Inputs for `--solve`:
- `--case <case.yml>`
- `--solver <solver.yml>`
- `--monitor <monitor.yml>`

Inputs for `--post-process`:
- `--post <post.yml>`
- either same invocation with `--solve`, or `--run-dir <existing_run_dir>`

MPI/local options:
- `-n, --num-procs`
  - applies to solver stage launch sizing.
  - post-processing defaults to one rank/task.

Cluster/Slurm options:
- `--cluster <cluster.yml>`
- `--scheduler slurm` (optional explicit selector)
- `--no-submit` (generate scripts/manifests only)

Preflight options:
- `--dry-run` (resolve and print launch/artifact plan only)
- `--format json` (machine-readable output for `--dry-run`)

Local example:
```bash
python3 scripts/pic.flow run --solve --post-process -n 8 \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml
```
In this command, solver runs with 8 ranks; post-processing still defaults to 1 rank.

Slurm example (generate + submit):
```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml
```

Slurm example (generate only):
```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml \
  --no-submit
```

Dry-run example (no file writes):
```bash
python3 scripts/pic.flow run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --dry-run --format json
```

@section sweep_sec 3. sweep: Parameter Study via Slurm Arrays

```bash
python3 scripts/pic.flow sweep \
  --study my_study/study.yml \
  --cluster my_study/cluster.yml [--no-submit]
```

Behavior:
- expands parameter matrix from `study.yml`
- materializes case directories under `studies/<study_id>/cases/`
- generates `solver_array.sbatch` and `post_array.sbatch`
- submits post array with `afterok:<solver_jobid>` dependency (unless `--no-submit`)
- aggregates metrics and emits plots in `studies/<study_id>/results/`

@section validate_sec 4. validate: Config-Only Checks

```bash
python3 scripts/pic.flow validate \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml
```

`validate` does not launch solver/post and does not create run/study artifacts.

@section artifacts_sec 5. Generated Runtime Artifacts

Single run (`run`):
- `runs/<run_id>/config/*.control`, `bcs*.run`, `whitelist.run`, `profile.run`, `post.run`
- `runs/<run_id>/scheduler/solver.sbatch`, `post.sbatch` (cluster mode)
- `runs/<run_id>/scheduler/submission.json` (cluster mode)
- `runs/<run_id>/manifest.json`

Sweep (`sweep`):
- `studies/<study_id>/cases/<case_i>/...`
- `studies/<study_id>/scheduler/case_index.tsv`
- `studies/<study_id>/scheduler/solver_array.sbatch`
- `studies/<study_id>/scheduler/post_array.sbatch`
- `studies/<study_id>/scheduler/submission.json`
- `studies/<study_id>/results/metrics_table.csv`
- `studies/<study_id>/results/plots/*.png` (if plotting enabled and matplotlib available)
- `studies/<study_id>/study_manifest.json`

@section next_steps_sec 6. Next Steps

- Config contract: **@subpage 14_Config_Contract**
- User workflows: **@subpage 11_User_How_To_Guides**
- Extensibility: **@subpage 17_Workflow_Extensibility**
- First-run path: **@subpage 02_Tutorial_Programmatic_Grid**
- Troubleshooting map: **@subpage 39_Common_Fatal_Errors**
- Smoke/CI verification details: **@subpage 40_Testing_and_Quality_Guide**
