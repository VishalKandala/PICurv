# PICurv

A parallel Eulerian-Lagrangian solver for incompressible flow and particle transport on curvilinear structured grids.

![PICurv Preview](docs/assets/curv.gif)

## What It Provides

- Incompressible flow solve (fractional-step projection) on curvilinear grids
- Particle tracking with PETSc `DMSwarm`
- Grid-particle interpolation and particle-grid projection
- Analytical flow modes for verification (`TGV3D`, `ZERO_FLOW`)
- YAML-driven orchestration through `./scripts/pic.flow`
- Slurm job generation/submission from YAML (`cluster.yml`)
- Parameter sweep orchestration with Slurm arrays (`study.yml`)
- Solver and postprocessor executables from one build system

## Requirements

- PETSc build available via `PETSC_DIR` (and `PETSC_ARCH` when required by your install)
- C toolchain + MPI (`gcc/clang`, `mpicc`, GNU Make)
- Python 3 with `pyyaml` and `numpy`
- Optional for plotting in sweep post-processing: `matplotlib`

## Quick Start (Local)

1. Set PETSc environment variables:
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
```

2. Build:
```bash
./scripts/pic.flow build
```

3. Initialize an example:
```bash
./scripts/pic.flow init flat_channel --dest my_case
```

`pic.flow` treats `case.yml`, `solver.yml`, `monitor.yml`, and `post.yml` as modular profiles.
You can reuse and recombine them instead of rewriting a monolithic config for every run.

4. Validate configs (no run yet):
```bash
./scripts/pic.flow validate \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

5. Preview planned launch/artifacts:
```bash
./scripts/pic.flow run --solve --post-process \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml \
  --dry-run
```

6. Run solver + post:
```bash
./scripts/pic.flow run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```
`-n/--num-procs` applies to the solver stage. Post-processing defaults to single-rank execution.

## Cluster and Sweep Workflow

Run on a cluster (Slurm):
```bash
./scripts/pic.flow run --solve --post-process \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml \
  --cluster my_case/slurm_cluster.yml
```

Launch a sweep study:
```bash
./scripts/pic.flow sweep \
  --study my_case/grid_independence_study.yml \
  --cluster my_case/slurm_cluster.yml
```

## Generated Artifacts

- Single run:
  - `runs/<run_id>/config/` generated C-facing runtime artifacts (`*.control`, `bcs*.run`, `post.run`, etc.)
  - `runs/<run_id>/scheduler/` scheduler scripts/manifests when cluster mode is used
  - `runs/<run_id>/visualization/` and/or post outputs
- Parameter sweep:
  - `studies/<study_id>/cases/` materialized case variants
  - `studies/<study_id>/scheduler/` array scripts and submission metadata
  - `studies/<study_id>/results/metrics_table.csv` and plots

## Smoke and Quality Checks

Manual smoke commands:

```bash
./scripts/pic.flow --help
./scripts/pic.flow run --help
./scripts/pic.flow validate --help

./scripts/pic.flow validate \
  --case tests/fixtures/valid/case.yml \
  --solver tests/fixtures/valid/solver.yml \
  --monitor tests/fixtures/valid/monitor.yml \
  --post tests/fixtures/valid/post.yml

./scripts/pic.flow run --solve --post-process \
  --case tests/fixtures/valid/case.yml \
  --solver tests/fixtures/valid/solver.yml \
  --monitor tests/fixtures/valid/monitor.yml \
  --post tests/fixtures/valid/post.yml \
  --dry-run --format json

python3 scripts/check_markdown_links.py
```

Automated smoke checks:

- Local (when pytest is available):
  - `pytest -q`
- CI:
  - `.github/workflows/quality.yml`

Detailed guide:
- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html

## Documentation (Live)

- Docs home: https://vishalkandala.me/picurv-docs/
- Getting started index: https://vishalkandala.me/picurv-docs/41_Getting_Started_Index.html
- 10-minute start: https://vishalkandala.me/picurv-docs/38_Start_Here_10_Minutes.html
- Installation guide: https://vishalkandala.me/picurv-docs/01_Installation.html
- Conductor CLI: https://vishalkandala.me/picurv-docs/05_The_Conductor_Script.html
- Workflow recipes and config cookbook: https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
- Config contract: https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- Config ingestion map: https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- Grid generator guide: https://vishalkandala.me/picurv-docs/48_Grid_Generator_Guide.html
- Extension playbook: https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
- Cluster guide: https://vishalkandala.me/picurv-docs/36_Cluster_Run_Guide.html
- Sweep guide: https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
- Testing and quality: https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- Common fatal errors: https://vishalkandala.me/picurv-docs/39_Common_Fatal_Errors.html
- Repository navigation: https://vishalkandala.me/picurv-docs/30_Repository_Navigation.html

## Repository Navigation

Top-level guides:
- `config/guide.md`
- `docs/guide.md`
- `examples/guide.md`
- `include/guide.md`
- `src/guide.md`
- `scripts/guide.md`
- `sandbox/guide.md`
- `logs/guide.md`

Build config location:
- `config/build/` (`config.local.mk`, `config.cluster.mk`, `config.petsc.mk`)

Grid profile library:
- `config/grids/` (shared `grid.gen` configs)

Scheduler/study profiles:
- `config/schedulers/`
- `config/studies/`

Build docs locally (if dependencies are available):
```bash
make build-docs
```

Docs output:
- `docs_build/html/index.html`

Doxygen warnings log:
- `logs/doxygen.warnings`

## Notes

- Only one `README.md` is maintained at repository root.
- Internal per-directory usage notes are standardized as `guide.md`.
