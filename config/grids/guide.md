# Grid Generator Config Library

This directory stores reusable `grid.gen` configuration files for programmatic structured-grid generation. These configs are useful when you want deterministic geometry generation controlled by parameters rather than manually curated grid files.

## Current Profile

- `coarse_square_tube_curved.cfg`: baseline curved square-tube style setup used for template-level workflows.

## How To Use

Generate a grid directly:

```bash
python3 scripts/grid.gen --config config/grids/coarse_square_tube_curved.cfg
```

Use generated grids through `case.yml` with `grid.mode: grid_gen` and `grid.generator.config_file`.

## Reproducibility Strategy

You can point `grid.generator.config_file` to either:

- a study-local file (recommended for frozen, reproducible case snapshots), or
- a shared library file in `config/grids/` (recommended while iterating).

Recommended workflow:

1. Develop and tune against a shared library config in `config/grids/`.
2. Once a study stabilizes, copy the config beside the case.
3. Reference that local copy from `case.yml` for long-term archival reproducibility.

## Common Failure Signals

- Cell orientation warnings after ingestion suggest generator parameter inconsistency.
- Unexpected domain extents often indicate stale config values or wrong config path.
- Projection/Poisson instability can be downstream symptoms of poor generated grid quality.

## Related Docs

- https://vishalkandala.me/picurv-docs/48_Grid_Generator_Guide.html
- https://vishalkandala.me/picurv-docs/07_Case_Reference.html
