@page 12_Capabilities_Summary Capabilities Summary

This page summarizes what is currently available through YAML + `pic.flow` without editing C source.

@tableofcontents

@section grid_sec 1. Grid and Domain

- Grid modes:
  - `programmatic_c` (C-side structured generation)
  - `file` (external `.picgrid` ingestion with validation and scaling)
  - `grid_gen` (pre-run Python grid generation orchestration)
- Domain topology:
  - single- and multi-block support
  - per-direction periodic flags
- Programmatic grid controls:
  - per-block geometry arrays (`im/jm/km`, bounds, stretching)
  - global DMDA partition hints (`da_processors_x/y/z`, scalar only)

@section flow_sec 2. Flow and Particle Physics

- Incompressible Navier-Stokes flow
- 2D/3D selection
- LES toggle
- Eulerian source modes:
  - `solve`
  - `load`
  - `analytical` + `analytical_type`
- Particle controls:
  - particle count
  - init modes (`Surface`, `Volume`, `PointSource`, `SurfaceEdges`)
  - restart mode (`init`, `load`)

@section numerics_sec 3. Numerics and Solvers

- Named momentum solver selection (extensible string selector)
- Solver-specific option blocks (`momentum_solver.dual_time_picard_rk4`)
- Pseudo-CFL and convergence tuning for dual-time method
- Pressure Poisson control:
  - multigrid levels/sweeps/semi-coarsening
  - PETSc passthrough options and monitors

@section bc_sec 4. Boundary Conditions

Supported combinations include:
- `INLET`: `constant_velocity`, `parabolic`
- `OUTLET`: `conservation`
- `WALL`: `noslip`
- `PERIODIC`: `geometric`, `constant_flux`

`pic.flow` validates face completeness, type-handler compatibility, and periodic pair consistency.

@section io_sec 5. Runtime I/O and Monitoring

- Output/restart/log root directory controls
- Eulerian/particle subdirectory controls
- Function whitelist logging and critical-function profiling files
- Raw PETSc monitor passthrough from monitor profile

@section post_sec 6. Postprocessing and Statistics

- Eulerian pipeline tasks:
  - dimensionalize
  - nodal averaging
  - Q-criterion
  - relative field normalization
- Lagrangian pipeline tasks:
  - specific kinetic energy
- Statistics pipeline:
  - MSD (`ComputeMSD`) with configurable output prefix
- Configurable post input extensions (`eulerianExt`, `particleExt`)

@section data_driven_sec 7. Data-Driven Closure Readiness

- Offline ML workflows are supported today via stable solver/post outputs.
- Tightly coupled in-solver inference requires extension work (documented in playbook).

@section orchestration_sec 8. Cluster and Study Orchestration

- Slurm single-run orchestration from `pic.flow run --cluster ...`
  - solver/post script generation
  - optional direct submission
  - dependency chaining (`post` after `solve`)
  - per-run submission and manifest metadata
- Slurm study sweeps from `pic.flow sweep`
  - parameter matrix expansion
  - solver/post array script generation
  - metrics aggregation from CSV/log outputs
  - convergence/sensitivity plot generation (matplotlib-enabled environments)

@section next_steps_sec 9. Next Steps

Proceed to **@subpage 13_Code_Architecture**.

For extension and maintenance details:
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 17_Workflow_Extensibility**
- **@subpage 21_Methods_Overview**
- **@subpage 31_Momentum_Solvers**
- **@subpage 32_Analytical_Solutions**
- **@subpage 33_Initial_Conditions**
- **@subpage 34_Particle_Model_Overview**
