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
