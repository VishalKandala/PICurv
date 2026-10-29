@page 53_Search_Robustness_Metrics_Reference Search Robustness Metrics Reference

@anchor _Search_Robustness_Metrics_Reference

PICurv writes `<run.analysis.metrics>/search_metrics.csv` for all particle-enabled runs. This page
is the canonical reference for the meaning, scope, formulas, and intended use
of those metrics.

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

These three signals are intentionally compact. They are meant to summarize
correctness, effort, and churn without forcing users to interpret every raw
counter during a run.

@section p53_taxonomy_sec 2. Metric Taxonomy

The metrics fall into three layers.

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
- `search_exited_count`
- `exited`
- `exited_cumulative`
- `blocked_walk_retry_count`
- `reflections`
- `reflection_fallbacks`
- `wraps`
- `occupied_cells`
- `particles_per_cell_min`
- `particles_per_cell_mean`
- `particles_per_cell_max`

Each row is indexed by `step`, `time` (solver time), and `physical_time` (seconds).
`lost_psi_sum` and `exited_psi_sum` are the sums of `Psi` over the lost and the exited
particles removed that step: with a 0/1 label and `iem_constant: 0`, `exited_psi_sum` is the
number of labelled particles that left the domain.

These are cheap, always-on aggregate counters. They are the authoritative source
for later derived signals and study reductions.

Derived per-timestep signals emitted by `logging.c`:

- `mean_traversal_steps = traversal_steps_sum / max(1, search_attempts)`
- `search_failure_fraction = search_lost_count / max(1, search_population)`
- `search_work_index = traversal_steps_sum / max(1, search_population)`
- `re_search_fraction = re_search_count / max(1, search_population)`

These are the metrics that should usually be inspected first in live monitoring
and figure-making because they normalize the raw counters into comparable
signals.

Run-level reductions typically used in sweeps or papers:

- `run_loss_fraction = last(lost_cumulative) / initial_particle_count`
- `run_swi_p95 = p95_t(search_work_index)`
- `run_research_p95 = p95_t(re_search_fraction)`
- `mean_migration_fraction = mean_t(migrated / max(1, search_population))`

These are not runtime counters. They are computed from the emitted CSV after the
run and are intended for comparison across parameter sweeps or example classes.

@section p53_semantics_sec 3. Semantics and Invariants

Important distinctions:

- `lost` counts search failures removed in cleanup: particles believed inside the domain that the search could not place. It is a failure signal; a nonzero value should be investigated.
- `exited` counts particles that left through a face whose `particle_policy` is `remove` (@ref p44_ppol_sec). Leaving is physics, not failure.
- `lost_cumulative` and `exited_cumulative` are run-local cumulative counts over the current run only; they are intentionally not restored from checkpoint state.
- `migrated` measures handoff exposure and stress. It is not a failure metric by itself.
- `search_attempts` counts robust walking searches, not every search-related operation.
- `mean_traversal_steps` measures average cost per robust walk.
- `search_work_index` measures total robust-search effort per original particle entering settlement.

Nominal invariant after settlement convergence:

- `search_population = search_located_count + search_lost_count + search_exited_count`

That invariant is important because it separates particles that entered search on
pass 1 from later re-search activity. Migrating particles may trigger extra
search work, but they should still end the timestep located, lost or exited. A lost or
exited particle is final within one settlement: only pass 1 searches it, so it is never
re-walked or re-migrated on later passes.

@section p53_definitions_sec 4. Detailed Metric Definitions

@subsection p53_definitions_population 4.1 Population and Outcome Counters

- `search_population`: number of particles entering settlement on pass 1 before any removal.
- `search_located_count`: number of particles from `search_population` whose final settlement status is `ACTIVE_AND_LOCATED`.
- `search_lost_count`: number of particles from `search_population` whose final settlement status is `LOST` before cleanup.

These three counters describe the top-level settlement outcome. `search_population`
is the most important denominator in the V2 design because it counts the original
problem size for the timestep. `search_located_count` and `search_lost_count`
then describe how that population resolved.

`search_population` is intentionally different from `search_attempts`. A single
particle can trigger more than one robust walk if migration or later settlement
passes require re-search. That is why `search_population` is the right
denominator for robustness and total-effort observables.

@subsection p53_definitions_work 4.2 Work and Churn Counters

- `search_attempts`: number of robust walking searches executed in the timestep.
- `traversal_steps_sum`: global sum of robust-walk traversal steps over the timestep.
- `re_search_count`: number of robust walks executed on settlement passes greater than 1.
- `max_traversal_steps`: largest traversal depth seen in any single robust walk during the timestep.
- `max_traversal_fail_count`: number of searches that fail by exhausting `MAX_TRAVERSAL`.

These counters describe how hard the search algorithm had to work.

`search_attempts` answers “how many walks were executed?” `traversal_steps_sum`
answers “how much total walk work was done?” `re_search_count` answers “how much
of that work came from churn after migration rather than from the original
population?”

`max_traversal_steps` is useful as a safety signal, but it should not be used as
a primary figure-of-merit because a single worst particle can dominate it.
`max_traversal_fail_count` is more serious: it indicates that some searches did
not merely get expensive, but exhausted the configured traversal budget.

@subsection p53_definitions_supporting 4.3 Supporting Diagnostic Counters

- `tie_break_count`: number of stuck-boundary tie-break acceptances.
- `boundary_clamp_count`: number of global-boundary clamp events during walking search.
- `bbox_guess_success_count`: number of invalid-cell particles whose bbox guess immediately resolved a remote owner.
- `bbox_guess_fallback_count`: number of invalid-cell particles that still required the robust walk after the bbox stage.
- `max_particle_pass_depth`: maximum settlement-pass depth reached by any particle during the timestep.
- `migrated`: number of particles migrated across rank ownership during the timestep.
- `migration_passes`: total settlement passes needed for the timestep.
- `load_imbalance`: max-rank particle load divided by mean-rank particle load.
- `occupied_cells`: physical cells holding at least one particle after the step.
- `particles_per_cell_min`, `particles_per_cell_max`: the fewest and most particles in any
  physical cell after the step; dummy cells are excluded.
- `particles_per_cell_mean`: particles counted in cells divided by the number of physical
  cells. A low minimum marks cells whose per-cell averages, such as the IEM mixing mean,
  rest on few particles.

These counters are not the main proposal-level observables, but they provide the
context needed to explain why search work changed.

`migrated` is especially important to interpret correctly. A large migrated count
can simply mean that the case is strongly exercising MPI handoff. It should be
read as a stress or exposure metric, not as a failure metric. Similarly,
`migration_passes` and `max_particle_pass_depth` indicate how much settlement
iteration was required, not whether the method was wrong.

The bbox-guess counters are useful for diagnosis of invalid-prior-cell paths.
`bbox_guess_success_count` indicates that owner resolution succeeded cheaply at
the bbox stage. `bbox_guess_fallback_count` indicates that the bbox stage did
not finish the job and the robust walk still had to be used.

@subsection p53_definitions_legacy 4.4 Loss and Exit Counters

- `search_exited_count`: number of particles from `search_population` whose final settlement status is `EXITED`.
- `lost`: per-step removals of `LOST` particles (search failures) during cleanup.
- `lost_cumulative`: run-local cumulative count of `lost`.
- `exited`: per-step removals of `EXITED` particles, which left through a removing face.
- `exited_cumulative`: run-local cumulative count of `exited`.

`lost` is the instantaneous failure signal. It is useful for spotting when a run first
begins to fail; in a closed or periodic domain it should be zero. `exited` is outflow:
in an open domain it should match the particle flux through the outlet faces.

In nominal paths, `lost` matches `search_lost_count` and `exited` matches
`search_exited_count`, because the particles that finish settlement `LOST` or `EXITED`
are the ones removed in cleanup.

@subsection p53_definitions_boundary 4.4a Boundary Event Counters

- `reflections`: wall reflections performed during the walk. Each is one mirror of the particle about the plane of a boundary face whose policy is `reflect`; an edge or corner crossing reflects about each crossed face in turn.
- `reflection_fallbacks`: particles that reached the reflection cap in one search and were placed just inside their boundary cell. Expected to be zero; a nonzero value also prints a WARNING line.
- `blocked_walk_retry_count`: walks that ended far beyond a boundary face (more than the blocked-walk guard of three cell sizes), treated as blocked by the domain's shape rather than as a real crossing, and searched again from scratch in the same pass.
- `wraps`: periodic wraps. Each crossed periodic face counts once, so an edge or corner crossing of two seams counts two.

Reflections and wraps measure boundary traffic, not trouble. In a statistically steady
flow, `reflections` should track the particle flux into the walls and `wraps` the flux
through the periodic seams. A per-step INFO line ("boundary events: ... wall reflections,
... periodic wraps, ... blocked walks searched again") repeats them when
`LOG_SEARCH_METRICS` is in `logging.enabled_functions`.

@subsection p53_definitions_derived 4.5 Derived Per-Timestep Signals

- `mean_traversal_steps = traversal_steps_sum / max(1, search_attempts)`
- `search_failure_fraction = search_lost_count / max(1, search_population)`
- `search_work_index = traversal_steps_sum / max(1, search_population)`
- `re_search_fraction = re_search_count / max(1, search_population)`

`mean_traversal_steps` is the average cost per robust walk. It is useful for
understanding how expensive an individual walk was on average, but it does not
capture extra work caused by repeated searches.

`search_failure_fraction` is the primary correctness observable. It answers: of
the particles that entered search this timestep, what fraction ended up lost?

`search_work_index` is the primary effort observable. It answers: how many total
robust-walk traversal steps were spent per original particle entering search?
This is why it is different from `mean_traversal_steps`. If migration churn or
re-search increases, `mean_traversal_steps` can stay modest while
`search_work_index` rises because the algorithm is doing more total work for the
same incoming population.

`re_search_fraction` is the primary churn observable. It answers: what fraction
of the incoming timestep population triggered later-pass robust walks? This is
why it is a good supporting metric for MPI handoff and settlement churn.

@subsection p53_definitions_reductions 4.6 Run-Level Reductions and Percentiles

- `run_loss_fraction = last(lost_cumulative) / initial_particle_count`
- `run_swi_p95 = p95_t(search_work_index)`
- `run_research_p95 = p95_t(re_search_fraction)`
- `mean_migration_fraction = mean_t(migrated / max(1, search_population))`

`p95_t(...)` means the 95th percentile over timesteps in one run. In plain
language, `run_swi_p95` is the value that 95 percent of timestep-level
`search_work_index` samples stay below. It is a “high but still typical” search
work level.

This is often a better proposal- or paper-level summary than a simple mean or a
max. A mean can hide short stressed periods. A max can be dominated by a single
outlier timestep. The 95th percentile captures the upper stressed behavior
without letting one pathological instant define the whole run.

`run_research_p95` is interpreted the same way, but for churn rather than work.
`mean_migration_fraction` is the average handoff exposure level across the run.

@section p53_interpretation_sec 5. Interpretation

Healthy baseline behavior usually looks like:

- `search_failure_fraction` near zero,
- `search_work_index` low and stable,
- `re_search_fraction` low for gentle cases and nonzero for migration-heavy cases,
- `max_traversal_steps` small,
- `tie_break_count` and `boundary_clamp_count` near zero in nominal cases,
- `reflection_fallbacks` and `blocked_walk_retry_count` zero; `exited` zero in a closed or
  periodic domain.

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
- strong `boundary_clamp_count` or repeated fallback behavior,
- `reflection_fallbacks` nonzero, or `blocked_walk_retry_count` large outside a strongly
  non-convex domain.

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
