@page 17_Workflow_Extensibility Workflow Extensibility Guide

This page captures practical extension directions that are already compatible with the current architecture.

@tableofcontents

@section goals_sec 1. Design Goal

Keep the pipeline easy to extend while preserving a strict, debuggable contract:

`YAML -> generated artifacts -> C ingestion -> runtime consumers`

This contract is already explicit in pages 14/15/16 and the ingress audit tooling.

@section grid_sec 2. Grid Workflow Extensions

Supported today:
- External grid ingestion (`grid.mode: file`) with validation and non-dimensionalization in `pic.flow`
- Pre-run generation (`grid.mode: grid_gen`) via `scripts/grid.gen`

Extension-friendly next steps:
- Add additional generator wrappers under `grid.generator` (new `grid_type` aliases)
- Add stronger pre-run grid quality checks in Python before solver launch
- Add optional strict checks that block run if metrics exceed thresholds

@section orchestration_sec 3. Multi-Run Orchestration (Sweeps/Studies)

Implemented in current `pic.flow`:
- `run --cluster <cluster.yml>` for Slurm-backed single-run scheduling
- `sweep --study <study.yml> --cluster <cluster.yml>` for matrix expansion + array submission
- post-stage dependency chaining (`afterok`)
- run/study manifests and submission metadata for reproducibility
- metrics aggregation and basic plot generation from existing outputs

Extension-ready next steps:
1. add more metric extractors (custom CSV schemas, multi-file reducers),
2. add sweep resume/retry policies,
3. add scheduler backends beyond Slurm without changing C contracts.

@section ml_sec 4. Data-Driven Particle Closure Integration

@subsection ml_offline_ssec 4.1 Offline (Recommended First)

- Use solver/post outputs as training/inference datasets.
- Keep ML scripts external (Python pipeline).
- Reinject inferred coefficients/flags through YAML for later runs.

This requires minimal C changes and is ideal for rapid iteration.

@subsection ml_coupled_ssec 4.2 Tightly Coupled Inference (Runtime)

For in-solver inference:
1. add runtime-selectable closure model interface in `ParticlePhysics.c`,
2. expose model type and parameters in schema,
3. map/validate in `pic.flow`,
4. ingest in `setup.c`,
5. wire in particle update path.

Keep deterministic fallback as default for robustness.

@section guardrails_sec 5. Guardrails for Safe Growth

- Keep new user options structured first; reserve passthrough for temporary gaps.
- Keep ingestion mostly in `setup.c`/`io.c`.
- Update mapping docs and ingress manifest together with code changes.
- Require at least one template/example update per new user-visible feature.

@section related_sec 6. Related Docs

- User contract: **@subpage 14_Config_Contract**
- Ingestion map: **@subpage 15_Config_Ingestion_Map**
- Extension playbook: **@subpage 16_Config_Extension_Playbook**
