# Search Robustness Characterization

## Purpose

This example family instruments the particle walking-search and migration
pipeline end to end. It is designed to answer:

- are particles being lost unexpectedly?
- how hard is the search working to settle them?
- how much migration and handoff churn is occurring?

The runtime writes `logs/search_metrics.csv` for all particle-enabled runs.
That CSV is the authoritative quantitative artifact for this family.

## Case Families

Brownian baselines:

- `search_robustness_cartesian.yml`: Cartesian `ZERO_FLOW` baseline with early migration on the documented `2 x 1 x 1` split
- `search_robustness_curvilinear.yml`: curvilinear `ZERO_FLOW` baseline with healthy but delayed migration on the bundled file grid

Deterministic migration stress cases:

- `search_robustness_cartesian_uniform.yml`: Cartesian analytical `UNIFORM_FLOW` transport/search stress case with repeatable advection-driven search work
- `search_robustness_curvilinear_uniform.yml`: curvilinear analytical `UNIFORM_FLOW` transport/search stress case on the bundled file grid

Shared assets:

- `Analytical-Zero.yml`: shared analytical zero-flow solver profile
- `Analytical-UniformFlow.yml`: shared analytical uniform-flow solver profile
- `Imp-MG-Standard.yml`: solve-mode curvilinear solver profile retained for the Brownian file-grid baseline
- `Search_Robustness_Output.yml`: shared monitor profile
- `search_robustness_analysis.yml`: optional particle-visualization/MSD post recipe
- `search_velocity_study.yml`: sweep-study starter for the deterministic curvilinear migration family
- `search_robustness_cartesian_dt_study.yml`: paired `dt`/`total_steps` sweep for the Cartesian Brownian baseline
- `search_robustness_curvilinear_dt_study.yml`: paired `dt`/`total_steps` sweep for the curvilinear Brownian baseline

## Why `ZERO_FLOW` Still Moves Particles

These cases use analytical `ZERO_FLOW`, which sets the carrier velocity to
zero. Particles still move because the particle update includes Brownian
diffusion when the effective diffusivity is nonzero. That makes these examples
good for isolating search robustness from carrier-flow solver and interpolation
effects.

## Metrics Contract

`logs/search_metrics.csv` is written once per timestep after particle settlement.

V1 compatibility columns are preserved and V2 adds:

- `search_population`
- `search_located_count`
- `search_lost_count`
- `traversal_steps_sum`
- `re_search_count`
- `max_traversal_fail_count`
- `search_failure_fraction`
- `search_work_index`
- `re_search_fraction`

Use **@subpage 53_Search_Robustness_Metrics_Reference** for the exact
formulas, prose definitions, and interpretation guidance for every V2 metric
and run-level reduction.

## Recommended Runs

Cartesian Brownian baseline:

```bash
./bin/picurv run --solve --post-process -n 2   --case  examples/search_robustness/search_robustness_cartesian.yml   --solver examples/search_robustness/Analytical-Zero.yml   --monitor examples/search_robustness/Search_Robustness_Output.yml   --post examples/search_robustness/search_robustness_analysis.yml
```

Curvilinear Brownian baseline:

```bash
./bin/picurv run --solve --post-process -n 2   --case  examples/search_robustness/search_robustness_curvilinear.yml   --solver examples/search_robustness/Imp-MG-Standard.yml   --monitor examples/search_robustness/Search_Robustness_Output.yml   --post examples/search_robustness/search_robustness_analysis.yml
```

Cartesian deterministic migration stress:

```bash
./bin/picurv run --solve --post-process -n 2   --case  examples/search_robustness/search_robustness_cartesian_uniform.yml   --solver examples/search_robustness/Analytical-UniformFlow.yml   --monitor examples/search_robustness/Search_Robustness_Output.yml   --post examples/search_robustness/search_robustness_analysis.yml
```

Curvilinear deterministic migration stress:

```bash
./bin/picurv run --solve --post-process -n 2   --case  examples/search_robustness/search_robustness_curvilinear_uniform.yml   --solver examples/search_robustness/Analytical-UniformFlow.yml   --monitor examples/search_robustness/Search_Robustness_Output.yml   --post examples/search_robustness/search_robustness_analysis.yml
```

## Interpretation Notes

- The Cartesian Brownian case is the fast baseline for healthy search plus early migration.
- The curvilinear Brownian case is a healthy delayed-migration baseline, not an immediate migration stress test.
- The two `UNIFORM_FLOW` cases are deterministic transport/search stress references; exact handoff onset depends on the grid geometry, seed placement, and imposed drift magnitude.
- `migrated` is a stress/exposure metric, not a failure metric.
- The official paper-grade signals for this family are `search_failure_fraction` and `search_work_index`.
- `re_search_fraction` is the optional third live/supporting signal.

## Console Behavior

The CSV is always written for particle-enabled runs. Compact console summaries
from `LOG_SEARCH_METRICS` appear only when:

- monitor verbosity is raised to `DEBUG` or higher, and
- `LOG_SEARCH_METRICS` remains in the function allow-list

That keeps production runs quiet while still preserving machine-readable search
observability.
