# Search Robustness Characterization

## Purpose

This example family instruments the particle walking-search and migration
pipeline end to end. It is designed to answer:

- are particles being lost unexpectedly?
- how hard is the search working to settle them?
- how much migration and handoff churn is occurring?

The runtime writes `logs/search_metrics.csv` for all particle-enabled runs.
That CSV is the authoritative quantitative artifact for this family.

## Why `ZERO_FLOW` Still Moves Particles

These cases use analytical `ZERO_FLOW`, which sets the carrier velocity to
zero. Particles still move because the particle update includes Brownian
diffusion when the effective diffusivity is nonzero. That makes these examples
good for isolating search robustness from carrier-flow solver and interpolation
effects.

## Files In This Family

- `search_robustness_cartesian.yml`: localized Cartesian stress case
- `search_robustness_curvilinear.yml`: localized curvilinear stress case
- `bent_channel_coarse.picgrid`: bundled curvilinear grid for the file-based variant
- `Analytical-Zero.yml`: shared analytical zero-flow solver profile
- `Imp-MG-Standard.yml`: solve-mode curvilinear solver profile for the file-grid variant
- `Search_Robustness_Output.yml`: shared monitor profile
- `search_robustness_analysis.yml`: optional particle-visualization/MSD post recipe

## Search Metrics CSV Contract

`logs/search_metrics.csv` is written once per timestep after particle settlement.

Columns:

- `step`
- `time`
- `total_particles`
- `lost`
- `lost_cumulative`
  cumulative lost-particle count over the current run only; this is intentionally not restored from restart state
- `migrated`
- `migration_passes`
- `search_attempts`
- `mean_traversal_steps`
- `max_traversal_steps`
- `tie_break_count`
- `boundary_clamp_count`
- `bbox_guess_success_count`
- `bbox_guess_fallback_count`
- `max_particle_pass_depth`
- `load_imbalance`

Interpretation:

- failure metrics:
  - `lost`
- effort and stability metrics:
  - `mean_traversal_steps`
  - `max_traversal_steps`
  - `tie_break_count`
  - `boundary_clamp_count`
  - `migration_passes`
  - `max_particle_pass_depth`
- exposure and context metrics:
  - `migrated`
  - `bbox_guess_success_count`
  - `bbox_guess_fallback_count`
  - `load_imbalance`

`migrated` is not a failure metric by itself. It tells you how much stress the
case is placing on the ownership-handoff logic.

## Recommended Runs

Cartesian localized stress:

```bash
./bin/picurv run --solve --post-process -n 2 \
  --case  examples/search_robustness/search_robustness_cartesian.yml \
  --solver examples/search_robustness/Analytical-Zero.yml \
  --monitor examples/search_robustness/Search_Robustness_Output.yml \
  --post examples/search_robustness/search_robustness_analysis.yml
```

Curvilinear localized stress:

```bash
./bin/picurv run --solve --post-process -n 2 \
  --case  examples/search_robustness/search_robustness_curvilinear.yml \
  --solver examples/search_robustness/Imp-MG-Standard.yml \
  --monitor examples/search_robustness/Search_Robustness_Output.yml \
  --post examples/search_robustness/search_robustness_analysis.yml
```

The Cartesian case fixes the seed on the documented x-rank interface for the
`2 x 1 x 1` decomposition in the case file. The curvilinear case uses the
bundled bent-channel grid in this directory with a mid-domain localized seed to
stress geometry and migration together. This curvilinear variant intentionally
keeps a solve-mode profile with zero initial conditions so search robustness is
characterized against the standard file-grid solve path, even though analytical
file-grid runs are now supported for `ZERO_FLOW` and `UNIFORM_FLOW`.

## Console Behavior

The CSV is always written for particle-enabled runs. Compact console summaries
from `LOG_SEARCH_METRICS` appear only when:

- monitor verbosity is raised to `DEBUG` or higher, and
- `LOG_SEARCH_METRICS` remains in the function allow-list

That keeps production runs quiet while still preserving machine-readable search
observability.
