@page 53_Search_Robustness_Metrics_Reference Search Robustness Metrics Reference

@anchor _Search_Robustness_Metrics_Reference

PICurv writes `logs/search_metrics.csv` for all particle-enabled runs. This page
is the canonical reference for the meaning, scope, and intended use of those
metrics.

@tableofcontents

@section p53_purpose_sec 1. Purpose

The search metrics serve four roles:

- runtime monitoring during active runs,
- debugging and branch diagnosis,
- study reduction for sweeps, and
- research-grade validation of the walking-search and migration pipeline.

The two official public signals are:

- `search_failure_fraction`
- `search_work_index`

The optional third supporting signal is:

- `re_search_fraction`

@section p53_taxonomy_sec 2. Metric Taxonomy

Raw per-timestep counters written directly by the runtime:

- `search_attempts`
- `search_population`
- `search_located_count`
- `search_lost_count`
- `traversal_steps_sum`
- `re_search_count`
- `max_traversal_steps`
- `max_traversal_fail_count`
- `tie_break_count`
- `boundary_clamp_count`
- `bbox_guess_success_count`
- `bbox_guess_fallback_count`
- `max_particle_pass_depth`
- `migrated`
- `migration_passes`
- `lost`
- `lost_cumulative`
- `load_imbalance`

Derived per-timestep signals emitted by `logging.c`:

- `mean_traversal_steps = traversal_steps_sum / max(1, search_attempts)`
- `search_failure_fraction = search_lost_count / max(1, search_population)`
- `search_work_index = traversal_steps_sum / max(1, search_population)`
- `re_search_fraction = re_search_count / max(1, search_population)`

Run-level reductions typically used in sweeps or papers:

- `run_loss_fraction = last(lost_cumulative) / initial_particle_count`
- `run_swi_p95 = p95_t(search_work_index)`
- `run_research_p95 = p95_t(re_search_fraction)`
- `mean_migration_fraction = mean_t(migrated / max(1, search_population))`

@section p53_semantics_sec 3. Semantics and Invariants

Important distinctions:

- `lost` is per-step cleanup loss.
- `lost_cumulative` is the run-local cumulative loss over the current run only; it is intentionally not restored from checkpoint state.
- `migrated` measures handoff exposure and stress. It is not a failure metric by itself.
- `search_attempts` counts robust walking searches, not every search-related operation.
- `mean_traversal_steps` measures average cost per robust walk.
- `search_work_index` measures total robust-search effort per original particle entering settlement.

Nominal invariant after settlement convergence:

- `search_population = search_located_count + search_lost_count`

@section p53_definitions_sec 4. Definitions

- `search_population`: number of particles entering settlement on pass 1 before any removal.
- `search_located_count`: number of particles from `search_population` whose final settlement status is `ACTIVE_AND_LOCATED`.
- `search_lost_count`: number of particles from `search_population` whose final settlement status is `LOST` before cleanup.
- `traversal_steps_sum`: global sum of robust-walk traversal steps over the timestep.
- `re_search_count`: number of robust walks executed on settlement passes greater than 1.
- `max_traversal_fail_count`: number of searches that fail by exhausting `MAX_TRAVERSAL`.
- `tie_break_count`: number of stuck-boundary tie-break acceptances.
- `boundary_clamp_count`: number of global-boundary clamp events during walking search.
- `bbox_guess_success_count`: number of invalid-cell particles whose bbox guess immediately resolved a remote owner.
- `bbox_guess_fallback_count`: number of invalid-cell particles that still required the robust walk after the bbox stage.
- `max_particle_pass_depth`: maximum settlement-pass depth reached by any particle during the timestep.
- `load_imbalance`: max-rank particle load divided by mean-rank particle load.

@section p53_interpretation_sec 5. Interpretation

Healthy baseline behavior usually looks like:

- `search_failure_fraction` near zero,
- `search_work_index` low and stable,
- `re_search_fraction` low for gentle cases and nonzero for migration-heavy cases,
- `max_traversal_steps` small,
- `tie_break_count` and `boundary_clamp_count` near zero in nominal cases.

Stressed-but-healthy behavior usually looks like:

- `search_failure_fraction` still near zero,
- `search_work_index` rising,
- `re_search_fraction` rising,
- `migrated` increasing,
- `migration_passes` or `max_particle_pass_depth` greater than 1.

Failing behavior usually looks like:

- `search_failure_fraction` nonzero,
- `lost` and `lost_cumulative` increasing,
- `search_work_index` large,
- `max_traversal_fail_count` nonzero,
- strong `boundary_clamp_count` or repeated fallback behavior.

@section p53_examples_sec 6. Worked Example Classes

Healthy Cartesian Brownian baseline:

- zero loss,
- early migration,
- `search_work_index` low,
- `re_search_fraction` nonzero once migration begins.

Healthy curvilinear Brownian baseline:

- zero loss,
- delayed migration,
- `search_work_index` low,
- `re_search_fraction` near zero early and growing later.

Pathological bad-seed failure case:

- high `search_failure_fraction`,
- high `search_work_index`,
- large `boundary_clamp_count`,
- `lost_cumulative` rises immediately.

@section p53_usage_sec 7. Recommended Usage

For live monitoring, watch:

- `search_failure_fraction`
- `search_work_index`
- `re_search_fraction`

For studies and papers, report at minimum:

- `run_loss_fraction`
- `run_swi_p95`

Supporting plots can add:

- `run_research_p95`
- `mean_migration_fraction`

@section p53_related_sec 8. Related Pages

- **@subpage 26_Walking_Search_Method**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 37_Sweep_Studies_Guide**
- `examples/search_robustness/search_robustness.md`
