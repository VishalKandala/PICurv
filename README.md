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

## Quick Start

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

4. Run solver + post:
```bash
./scripts/pic.flow run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

5. Run on a cluster (Slurm):
```bash
./scripts/pic.flow run --solve --post-process \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml \
  --cluster my_case/slurm_cluster.yml
```

6. Launch a sweep study:
```bash
./scripts/pic.flow sweep \
  --study my_case/grid_independence_study.yml \
  --cluster my_case/slurm_cluster.yml
```

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

## Documentation

Authoritative docs entry points:
- `docs/mainpage.md`
- `docs/pages/14_Config_Contract.md`
- `docs/pages/15_Config_Ingestion_Map.md`
- `docs/pages/16_Config_Extension_Playbook.md`
- `docs/pages/30_Repository_Navigation.md`

Build docs:
```bash
make build-docs
```

Docs output:
- `docs/docs_build/html/index.html`

Doxygen warnings log:
- `logs/doxygen.warnings`

## Notes

- Only one `README.md` is maintained at repository root.
- Internal per-directory usage notes are standardized as `guide.md`.
