# Study Config Guide

This directory stores reusable `study.yml` definitions for `picurv sweep` workflows. A study file turns one validated base workflow into a controlled parameter experiment with reproducible case generation and automated metric collection.

## Typical Study Contents

- parameter-space definitions (what is varied and over what range),
- metric extraction specs (what performance/accuracy signal to compare),
- optional plotting directives,
- execution controls for scheduler array behavior.

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

## Reference Files

- `examples/master_template/master_study.yml`
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
