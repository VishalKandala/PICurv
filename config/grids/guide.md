# Grid Generator Config Library

This directory stores reusable `grid.gen` configuration files.

Current profile:
- `coarse_square_tube_curved.cfg`: baseline curved square-tube style setup.

Use with:
```bash
python3 scripts/grid.gen --config config/grids/coarse_square_tube_curved.cfg
```

In `case.yml` (`grid.mode: grid_gen`), `grid.generator.config_file` can point to:
- a study-local file (recommended for reproducible case snapshots), or
- a shared library file in `config/grids/`.

Recommended workflow:

1. develop and tune against a shared library config in `config/grids/`,
2. once a study stabilizes, copy the config beside the case for a self-contained snapshot,
3. reference that local copy from `case.yml` if long-term reproducibility matters more than central reuse.

Related docs:

- https://vishalkandala.me/picurv-docs/48_Grid_Generator_Guide.html
- https://vishalkandala.me/picurv-docs/07_Case_Reference.html
