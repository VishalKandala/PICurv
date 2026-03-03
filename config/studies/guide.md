# Study Config Guide

This directory stores reusable `study.yml` definitions for `picurv sweep` workflows.

## Typical Contents

- parameter-space definitions,
- metric extraction specs,
- optional plotting controls,
- execution controls for scheduler arrays.

## How To Use

```bash
./bin/picurv sweep \
  --study config/studies/grid_independence_example.yml \
  --cluster <cluster.yml>
```

## Reference Files

- `examples/master_template/master_study.yml`
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
