@mainpage PICurv Solver Documentation

@section intro_sec Introduction

# PICurv: A Hybrid Eulerian-Lagrangian Framework for Scalar Transport

PICurv is a parallel CFD and particle framework for incompressible flow and scalar transport. It combines:

- Eulerian flow solves on curvilinear grids (PETSc `DMDA` based)
- Lagrangian particle transport (`DMSwarm` based)
- A YAML-driven run pipeline through `./scripts/pic.flow`

The user workflow is configuration-first: YAML inputs are validated and translated into generated runtime artifacts consumed by the C solver and postprocessor.

@section get_started_main To Begin, Choose Your Path

- **@subpage 38_Start_Here_10_Minutes**: one-page first run path (validate -> dry-run -> execute).
- **@subpage 01_Installation**: install dependencies, build binaries, run first case.
- **@subpage 05_The_Conductor_Script**: YAML-first workflow and command patterns.
- **@subpage 39_Common_Fatal_Errors**: map frequent failures to fix commands.
- **@subpage 13_Code_Architecture**: C-side architecture and extension touchpoints.
- **@subpage 30_Repository_Navigation**: repository map and per-directory `guide.md` pointers.

@section preview_sec Quick Preview

@image html assets/curv.gif "PICurv simulation preview" width=900px

@image html assets/paraview_flat_channel.png "Example visualization output (ParaView)" width=900px

@section methodology_sec Core Methodology

PICurv uses a two-way coupled Eulerian-Lagrangian strategy:

1. Eulerian phase: velocity/pressure evolution on a structured curvilinear grid.
2. Lagrangian phase: particle motion and particle-carried quantities.
3. Coupling:
   - Grid -> particle interpolation for advection/source evaluation.
   - Particle -> grid projection for reconstructed Eulerian particle-derived fields.

@section features_sec Key Features

- YAML-driven pipeline (`case.yml`, `solver.yml`, `monitor.yml`, `post.yml`)
- Cluster scheduler integration (`cluster.yml`) for Slurm job generation/submission
- Study/sweep orchestration (`study.yml`) with Slurm arrays, dependency chaining, and metric plots
- Multiple grid ingestion modes: `programmatic_c`, `file`, and `grid_gen`
- Named momentum solver strategy with solver-specific option blocks
- Structured post-processing pipelines (Eulerian, Lagrangian, statistics)
- Dedicated config contract and ingestion mapping docs for maintenance and extension

For release notes and recent contract changes, see **@subpage 18_Changelog**.
For solver-method overviews, see **@subpage 21_Methods_Overview**.
