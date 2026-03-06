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
- Pre-run generation (`grid.mode: grid_gen`) via `scripts/grid.gen`

Extension-friendly next steps:
- Add additional generator wrappers under `grid.generator` (new `grid_type` aliases)
- Add stronger pre-run grid quality checks in Python before solver launch
- Add optional strict checks that block run if metrics exceed thresholds

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

@section p17_ml_sec 4. Data-Driven Particle Closure Integration

@subsection p17_ml_offline_ssec 4.1 Offline (Recommended First)

- Use solver/post outputs as training/inference datasets.
- Keep ML scripts external (Python pipeline).
- Reinject inferred coefficients/flags through YAML for later runs.

This requires minimal C changes and is ideal for rapid iteration.

@subsection p17_ml_coupled_ssec 4.2 Tightly Coupled Inference (Runtime)

For in-solver inference:
1. add runtime-selectable closure model interface in `ParticlePhysics.c`,
2. expose model type and parameters in schema,
3. map/validate in `picurv`,
4. ingest in `setup.c`,
5. wire in particle update path.

Keep deterministic fallback as default for robustness.

@section p17_guardrails_sec 5. Guardrails for Safe Growth

- Keep new user options structured first; reserve passthrough for temporary gaps.
- Keep ingestion mostly in `setup.c`/`io.c`.
- Update mapping docs and ingress manifest together with code changes.
- Require at least one template/example update per new user-visible feature.

@section p17_related_sec 6. Related Docs

- User contract: **@subpage 14_Config_Contract**
- Ingestion map: **@subpage 15_Config_Ingestion_Map**
- Extension playbook: **@subpage 16_Config_Extension_Playbook**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Workflow Extensibility Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

