# Study Config Guide

This directory stores reusable `study.yml` definitions for `picurv sweep` workflows. A study file turns one validated base workflow into a controlled parameter experiment with reproducible case generation and automated metric collection.

## Typical Study Contents

- parameter-space definitions (what is varied and over what range),
- optional explicit `parameter_sets` bundles when multiple keys must move together,
- metric extraction specs (what performance/accuracy signal to compare),
- optional plotting directives,
- execution controls for scheduler array behavior.

Parameter keys may target nested config values such as
`case.models.physics.particles.count`. Use `parameters` for cross-product sweeps and
`parameter_sets` when you need explicit coupled bundles instead.

Use explicit metric specs whenever the study signal comes from a non-default
CSV or log artifact. For example, interpolation studies can aggregate values
from `logs/interpolation_error.csv` instead of using `msd_final`. Search and
migration characterization studies can likewise aggregate columns from
`logs/search_metrics.csv`, such as `search_failure_fraction`,
`search_work_index`, `re_search_fraction`, or normalized run-level signals
derived from `lost_cumulative`.

Study CSV metrics also support `reduction: p95`, per-row ratios via
`numerator_column` plus `denominator_column`, and normalization against the
study parameter space with `normalize_by_parameter`.

## How To Use

```bash
./bin/picurv sweep \
  --study config/studies/grid_independence_example.yml \
  --cluster <cluster.yml>
```

## Recommended Workflow

1. Validate one base case end-to-end before creating sweeps.
2. Keep study parameters physically meaningful and isolated where possible.
3. Start with a small sweep subset to verify orchestration and metric extraction.
4. Scale up parameter space only after logs/results schema look correct.

If you are not using a scheduler, emulate a study locally by repeating
`picurv run` with a few controlled config variants instead of using
`picurv sweep`.

## Reference Files

- `examples/master_template/master_study.yml`
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
