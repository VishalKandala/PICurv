@page 12_Capabilities_Summary Capabilities Summary

@anchor _Capabilities_Summary

This page summarizes current capabilities from YAML + `picurv` without editing C source.
It is organized by workflow stage rather than just a feature bullet list.

@tableofcontents

@section p12_ingest_sec 1. Input and Grid Capabilities

PICurv currently supports three grid ingestion modes:

- `programmatic_c`: C-side structured grid generation,
- `file`: external `.picgrid` read path with scaling/validation,
- `grid_gen`: pre-run Python generator orchestration.

Domain controls include:

- single- and multi-block support,
- per-direction periodicity,
- optional DMDA partition hints (`da_processors_x/y/z`).

@section p12_physics_sec 2. Physics and Model Selection

Supported high-level operation modes:

- solve from numerically evolved Eulerian fields,
- load/restart from prior field outputs,
- analytical Eulerian field modes (`TGV3D`, `ZERO_FLOW`, `UNIFORM_FLOW`).
- file-grid analytical support for the non-custom analytical modes (`ZERO_FLOW`, `UNIFORM_FLOW`).

Particle controls include:

- particle count,
- initialization modes (`Surface`, `Volume`, `PointSource`, `SurfaceEdges`),
- restart modes (`init`, `load`),
- grid-to-particle interpolation method (`Trilinear` direct cell-center or `CornerAveraged` legacy),
- scalar micromixing update path (IEM-style `Psi` model).

@section p12_solver_sec 3. Numerical Solver Stack

Momentum:

- named momentum strategy selection,
- active implementations: explicit RK and dual-time Picard RK4,
- tunable tolerances and pseudo-CFL controls.

Pressure:

- multigrid Poisson workflow,
- level/sweep/semi-coarsening controls,
- PETSc passthrough flags for advanced tuning.

See method details in **@subpage 21_Methods_Overview**.

@section p12_bc_sec 4. Boundary and Runtime Controls

Boundary capabilities include validated type-handler pairings across inlet/outlet/wall/periodic classes.
Runtime controls include:

- output/restart/log directory selection,
- function-level logging allowlists,
- profiling critical function lists,
- monitor verbosity and cadence controls.

@section p12_post_sec 5. Post-Processing and Statistics

Pipeline capabilities include:

- Eulerian transforms (dimensionalization, nodal averaging, Q-criterion, normalization),
- Lagrangian particle tasks,
- statistics reduction pipeline (currently MSD family),
- configurable input extensions and output field selection.

@section p12_orchestration_sec 6. Cluster and Study Orchestration

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

@section p12_extension_sec 7. Extensibility Status

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

@section p12_next_steps_sec 8. Suggested Reading Order

1. **@subpage 13_Code_Architecture**
2. **@subpage 21_Methods_Overview**
3. **@subpage 31_Momentum_Solvers**
4. **@subpage 34_Particle_Model_Overview**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Capabilities Summary** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

