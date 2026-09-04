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

`grid_gen` offers two composed geometries rather than a list of named shapes: `box`, a
Cartesian block whose bounding walls are piecewise height fields, and `sweep`, a
cross-section carried along a piecewise centreline. Both accept an ordered placement and
similarity transform list. Both are experimental. See **@subpage 48_Grid_Generator_Guide**.

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
wall-modelled large-eddy simulation is coupled in both directions, and the modelled
stress reaches the momentum equation through an effective eddy viscosity installed at
the wall face.

@note A wall model is not independent of the turbulence model in the way its
configuration placement suggests. It supplies the stress of a layer the mesh does not
resolve, which is only meaningful if the unresolved motions are modelled somewhere, so
a wall model with LES and RANS both disabled is rejected - there is no implicit-LES
scheme here to stand in. `cabot` and `werner` are rejected under RANS, being large-eddy
constructs, and a wall model on a laminar case is rejected outright. Whether the first
cell falls in the selected law's valid range depends on the mesh and is checked at
runtime instead.

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

- fixed run-owned output, restart-input, log, and analysis homes,
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

@section p12_lifecycle_sec 7. Workspace, Assets, and Data Lifecycle

`picurv init` creates a workspace whose editable configuration, imported inputs,
reusable assets, runs, and studies each have a fixed home. Run-owned directory names
are part of that contract rather than a monitor setting; the topology and the guards
that hold it are in @ref p71_topology_sec.

External files enter a workspace explicitly. `picurv inputs import <kind> <source>`
records a checksum and an ownership mode, and `reference` mode registers a path
without copying it, so an unavailable external target fails loudly instead of
silently. Entries at @ref p05_cap_input_mode_sec.

`picurv precompute` resolves the grid, initial-condition, and inlet-profile provider
graph and publishes each result as an immutable content-addressed object in the
workspace asset store. Identity covers the normalized provider settings, the checksums
of referenced files, and the PICurv build, so an unchanged input reuses its object and
any change selects a new one. Runs materialize what they need by reflink, hardlink, or
copy and record the mapping. `--require-precomputed` refuses to build a missing object
instead of quietly rebuilding it; `--fetch-missing` looks in configured storage first.
A provider that only exists inside the simulator is reported rather than imitated.

Branching a run with `--restart-from` carries the checkpoint's physical fields but
starts field-statistics accumulators empty unless `--statistics-state carry` asks for
the saved window state. Entries at @ref p52_cap_restart_stats_sec.

`picurv storage` moves finished runs, studies, and study members through an `rclone`
remote. `protect` uploads and verifies while keeping every local file; `offload` does
the same and then prunes the verified payload according to an offload policy that
chooses what stays local, with `--retain`/`--drop` adjusting that preset one component
at a time. Compression policy entries are at @ref p61_cap_comp_sec, offload policy
entries at @ref p61_cap_policy_sec, and retention component entries at
@ref p61_cap_retain_sec. Cold data is marked, so continuation,
post-processing, submission, and study reaggregation refuse it with the restore command
rather than mistaking pruned output for output that was never produced.

`picurv version`, `picurv versions`, and `picurv source` report and select the release
and build identity that runs, checkpoints, and archives are stamped with.
`picurv version status` additionally validates that the conductor, both executables, and
any workspace version requirement agree, and exits non-zero when they do not.

@section p12_inspection_sec 8. Run Inspection and Plotting

`picurv summarize` supports:

- curated run, case, solver, and monitor configuration dashboards,
- selected-step health summaries from continuity, particles, momentum, Poisson, profiling, memory, and convergence logs,
- discovery of available plottable numeric histories with `--list-plot-series`,
- full or last-N append-order time-history plots with `--plot` and `--last`,
- per-block and per-function lines, automatic positive residual/norm log scaling, explicit saves, and headless fallback.

PICurv owns run/log interpretation; standalone `generators/plot.gen` renders
versioned normalized JSON requests without knowing run-directory layouts.

@section p12_extension_sec 9. Extensibility Status

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

@section p12_next_steps_sec 10. Suggested Reading Order

1. **@subpage 13_Code_Architecture**
2. **@subpage 21_Methods_Overview**
3. **@subpage 31_Momentum_Solvers**
4. **@subpage 34_Particle_Model_Overview**
