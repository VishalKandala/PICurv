# Profile Generation Guide

This directory documents reusable inlet-profile patterns for `prescribed_flow`
boundary conditions. The actual profile choice lives in `case.yml` because the
BC face and grid dimensions determine the required `PICSLICE` shape.

## Supported Patterns

File-backed profile:

```yaml
boundary_conditions:
  - face: "-Zeta"
    type: INLET
    handler: prescribed_flow
    params:
      source:
        type: file
        path: profiles/inlet.picslice
```

Generated square-duct Poiseuille profile:

```yaml
boundary_conditions:
  - face: "-Zeta"
    type: INLET
    handler: prescribed_flow
    params:
      source:
        type: generated
        generator: square_duct_poiseuille
        params:
          bulk_velocity: 1.0
          n_terms: 101
```

`picurv run --solve` automatically materializes generated profiles under
`runs/<run_id>/config/`, writes `profile.info`, stages a nondimensional solver
copy, and passes the staged `source_file` to the C runtime.

Use `picurv precompute --case case.yml --output-dir precomputed/<name>` when you
want to generate and inspect profile artifacts before launching the solver.

## Standalone Engine

`scripts/profile.gen` is the low-level generator:

```bash
python3 scripts/profile.gen square_duct_poiseuille \
  --output inlet.picslice \
  --dims 63 63 \
  --bulk-velocity 1.0 \
  --n-terms 101
```

Prefer the `picurv` conductor for normal cases, because it derives dimensions
from `case.yml`, `grid.mode`, block index, and BC face.
