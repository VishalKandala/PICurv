@page 12_Capabilities_Summary Capabilities Summary

This page summarizes current capabilities from YAML + `pic.flow` without editing C source.
It is organized by workflow stage rather than just a feature bullet list.

@tableofcontents

@section ingest_sec 1. Input and Grid Capabilities

PICurv currently supports three grid ingestion modes:

- `programmatic_c`: C-side structured grid generation,
- `file`: external `.picgrid` read path with scaling/validation,
- `grid_gen`: pre-run Python generator orchestration.

Domain controls include:

- single- and multi-block support,
- per-direction periodicity,
- optional DMDA partition hints (`da_processors_x/y/z`).

@section physics_sec 2. Physics and Model Selection

Supported high-level operation modes:

- solve from numerically evolved Eulerian fields,
- load/restart from prior field outputs,
- analytical Eulerian field modes (`TGV3D`, `ZERO_FLOW`).

Particle controls include:

- particle count,
- initialization modes (`Surface`, `Volume`, `PointSource`, `SurfaceEdges`),
- restart modes (`init`, `load`),
- scalar micromixing update path (IEM-style `Psi` model).

@section solver_sec 3. Numerical Solver Stack

Momentum:

- named momentum strategy selection,
- active implementations: explicit RK and dual-time Picard RK4,
- tunable tolerances and pseudo-CFL controls.

Pressure:

- multigrid Poisson workflow,
- level/sweep/semi-coarsening controls,
- PETSc passthrough flags for advanced tuning.

See method details in **@subpage 21_Methods_Overview**.

@section bc_sec 4. Boundary and Runtime Controls

Boundary capabilities include validated type-handler pairings across inlet/outlet/wall/periodic classes.
Runtime controls include:

- output/restart/log directory selection,
- function-level logging allowlists,
- profiling critical function lists,
- monitor verbosity and cadence controls.

@section post_sec 5. Post-Processing and Statistics

Pipeline capabilities include:

- Eulerian transforms (dimensionalization, nodal averaging, Q-criterion, normalization),
- Lagrangian particle tasks,
- statistics reduction pipeline (currently MSD family),
- configurable input extensions and output field selection.

@section orchestration_sec 6. Cluster and Study Orchestration

Single-run cluster flow (`run --cluster ...`):

- scheduler script generation,
- optional submission,
- solver/post dependency chaining,
- run manifests.

Study flow (`sweep`):

- parameter matrix expansion,
- array script generation,
- metric aggregation and optional plots,
- study manifest and reproducible directory structure.

@section extension_sec 7. Extensibility Status

Current extension pathways are documented and active for:

- YAML contract extension,
- ingestion mapping updates,
- workflow orchestration growth,
- method-level and model-level solver extension.

Reference pages:

- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 17_Workflow_Extensibility**

@section next_steps_sec 8. Suggested Reading Order

1. **@subpage 13_Code_Architecture**
2. **@subpage 21_Methods_Overview**
3. **@subpage 31_Momentum_Solvers**
4. **@subpage 34_Particle_Model_Overview**
