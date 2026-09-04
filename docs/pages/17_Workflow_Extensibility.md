@page 17_Workflow_Extensibility Workflow Extensibility Guide

@anchor _Workflow_Extensibility

This page captures practical extension directions that are already compatible with the current architecture.

@tableofcontents

@section p17_goals_sec 1. Design Goal

Keep the pipeline easy to extend while preserving a strict, debuggable contract:

`YAML -> generated artifacts -> C ingestion -> runtime consumers`

This contract is already explicit in pages 14/15/16 and the ingress audit tooling.

@section p17_grid_sec 2. Grid Workflow Extensions

Supported today:
- External grid ingestion (`grid.mode: file`) with validation and non-dimensionalization in `picurv`
- Pre-run generation (`grid.mode: grid_gen`) via `generators/grid.gen`, composed from `box` (a
  Cartesian block with piecewise-shaped walls) and `sweep` (a cross-section swept along a
  centreline); see **@subpage 48_Grid_Generator_Guide**
- Structural pre-run refusals in `grid.gen` itself: a multigrid-illegal node count, a
  periodic axis with an asymmetric seam or a rotation applied to it, and a wall or path
  segment too short to represent

Extension-friendly next steps:
- Add additional `grid_type` geometries beyond `box`/`sweep`
- Add optional strict checks that block a run when the *computed* quality metrics
  (`Max_Non_Orthogonality_deg`, `Max_Aspect_Ratio`, ...) exceed a configured threshold,
  rather than only reporting them in `.info` - the structural refusals above catch what
  cannot work; nothing yet catches what merely works badly

@section p17_orchestration_sec 3. Multi-Run Orchestration (Sweeps/Studies)

Implemented in current `picurv`:
- `run --cluster <cluster.yml>` for Slurm-backed single-run scheduling
- `sweep --study <study.yml> --cluster <cluster.yml>` for matrix expansion + array submission
- post-stage dependency chaining (`afterok`)
- run/study manifests and submission metadata for reproducibility
- metrics aggregation and basic plot generation from existing outputs

Extension-ready next steps:
1. add more metric extractors (custom CSV schemas, multi-file reducers),
2. add sweep resume/retry policies,
3. add scheduler backends beyond Slurm without changing C contracts.

@section p17_completion_sec 4. Runtime Completion Extensions

The current runtime already supports safe early shutdown for walltime and signal
events. Physics-based completion is a separate future extension and should keep
that same checkpoint discipline.

Extension-ready next steps:
1. add optional convergence-stop criteria under
   `monitor.yml -> solution_monitoring.convergence`,
2. limit deterministic fixed-state stopping to `mode: steady_deterministic`,
3. use explicit tolerances, minimum samples, and dwell windows before declaring completion,
4. keep `mode: transient` diagnostic-only because nonzero drift can be the expected physics,
5. route any automatic completion through the existing final-write/safe-checkpoint path.

For statistically steady turbulence, prefer postprocessing-heavy turbulence
analysis after runtime window metrics indicate stationarity; runtime completion
criteria should remain lightweight and conservative.

@section p17_ml_sec 5. Data-Driven Particle Closure Integration

@subsection p17_ml_offline_ssec 5.1 Offline (Recommended First)

- Use solver/post outputs as training/inference datasets.
- Keep ML scripts external (Python pipeline).
- Reinject inferred coefficients/flags through YAML for later runs.

This requires minimal C changes and is ideal for rapid iteration.

@subsection p17_ml_coupled_ssec 5.2 Tightly Coupled Inference (Runtime)

For in-solver inference:
1. add runtime-selectable closure model interface in `ParticlePhysics.c`,
2. expose model type and parameters in schema,
3. map/validate in `picurv`,
4. ingest in `setup.c`,
5. wire in particle update path.

Keep deterministic fallback as default for robustness.

@section p17_guardrails_sec 6. Guardrails for Safe Growth

- Keep new user options structured first; reserve passthrough for temporary gaps.
- Keep ingestion mostly in `setup.c`/`io.c`.
- Update mapping docs and ingress manifest together with code changes.
- Require at least one template/example update per new user-visible feature.

@section p17_related_sec 7. Related Docs

- User contract: **@subpage 14_Config_Contract**
- Ingestion map: **@subpage 15_Config_Ingestion_Map**
- Extension playbook: **@subpage 16_Config_Extension_Playbook**
