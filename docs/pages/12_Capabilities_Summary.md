@page 12_Capabilities_Summary Capabilities Summary

@anchor _Capabilities_Summary

@pagemeta{Reference, All readers, Capability statuses enforced by make audit-capability}

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
- per-direction geometric periodicity for Eulerian fields, derived from paired
  BCs and requiring matching surfaces under a constant translation,
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

Particle positions are not currently wrapped across periodic boundaries.

@subsection p12_turbulence_ssec 2.1 Turbulence Models

Two LES closures are selectable. `constant_smagorinsky` applies a prescribed
coefficient and allocates no coefficient field. `dynamic_smagorinsky` measures the
coefficient each update through the Germano identity and Lilly's least-squares
contraction, with selectable grid filter width, test-filter kernel and width ratio,
coefficient averaging set, and limiting policy. Entries and full detail at
@ref p07_les_sec; the formulation is derived in @ref 72_LES_Turbulence_Closure.

@note Both models are **experimental**: they are implemented and unit-tested, but
neither has been validated against a reference flow. The check that would settle the
dynamic model is decaying isotropic turbulence with homogeneous averaging, where
`Cs(t)` should settle near 0.16-0.17. Until such a run is recorded, treat coefficient
magnitudes as uncharacterized.

RANS (`k_omega`) is accepted by the configuration layer but its runtime update is
incomplete. Wall functions are configured separately from both, and offer three
laws - `log_law`, `werner`, and `cabot`. The correction is applied inside the
momentum solve and again before the LES strain rates are formed, so a
wall-modelled large-eddy simulation is coupled in both directions.

@section p12_solver_sec 3. Numerical Solver Stack

Momentum:

- named momentum strategy selection,
- active implementations: `Explicit RK4`, `Dual Time Picard Jameson RK`, and
  `Newton Krylov` (see @ref p08_entries_sec for the comparison and status of each),
- tunable tolerances and pseudo-CFL controls (Picard-Jameson only).

@note `Explicit RK4` is `experimental`: it has no positive-path verification harness.
`Dual Time Picard Jameson RK` is the production default.

Pressure:

- multigrid Poisson workflow,
- level/sweep/semi-coarsening controls,
- PETSc passthrough flags for advanced tuning.

See method details in **@subpage 21_Methods_Overview**.

@section p12_bc_sec 4. Boundary and Runtime Controls

Boundary capabilities include validated type-handler pairings across inlet/outlet/wall/periodic classes.
Runtime controls include:

- <run.solver_output>/restart/log directory selection,
- function-level logging allowlists,
- profiling critical function lists,
- monitor verbosity and cadence controls,
- online field-statistics windows accumulated during the solve, checkpointed with
  the flow state and resumed on continuation.

@section p12_post_sec 5. Post-Processing and Statistics

Pipeline capabilities include:

- Eulerian transforms (dimensionalization, nodal averaging, Q-criterion, normalization),
- Lagrangian particle tasks,
- particle statistics reduction pipeline (currently MSD family),
- derived Eulerian field statistics: Reynolds stresses, RMS, turbulent kinetic
  energy, and turbulent fluxes from windows the solver accumulated online,
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

@section p12_inspection_sec 7. Run Inspection and Plotting

`picurv summarize` supports:

- curated run, case, solver, and monitor configuration dashboards,
- selected-step health summaries from continuity, particles, momentum, Poisson, profiling, memory, and convergence logs,
- discovery of available plottable numeric histories with `--list-plot-series`,
- full or last-N append-order time-history plots with `--plot` and `--last`,
- per-block and per-function lines, automatic positive residual/norm log scaling, explicit saves, and headless fallback.

PICurv owns run/log interpretation; standalone `generators/plot.gen` renders
versioned normalized JSON requests without knowing run-directory layouts.

@section p12_extension_sec 8. Extensibility Status

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

@section p12_next_steps_sec 9. Suggested Reading Order

1. **@subpage 13_Code_Architecture**
2. **@subpage 21_Methods_Overview**
3. **@subpage 31_Momentum_Solvers**
4. **@subpage 34_Particle_Model_Overview**
