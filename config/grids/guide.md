# Grid Generator Config Library

This directory stores reusable `grid.gen` configuration files for programmatic structured-grid generation. These configs are useful when you want deterministic geometry generation controlled by parameters rather than manually curated grid files.

## Current Profiles

- `coarse_square_tube_curved.cfg`: baseline curved square-tube style setup used for template-level workflows.

Driven-periodic validation grids. All use grid type `warp` with
`amp_A = amp_B = amp_C = 0`, which reduces the warping map to the identity and
leaves an exactly Cartesian box, while still giving access to the shared
stretching controls on the parent parser:

| Config | Case | Cells |
|---|---|---|
| `plane_channel_laminar.cfg` | laminar exact verification | 17 x 33 x 17 |
| `plane_channel_les_retau180.cfg` | constant-Smagorinsky LES, `Re_tau = 180` | 33 x 65 x 33 |
| `plane_channel_retau180.cfg` | channel DNS, `Re_tau = 180` | 129^3 |
| `plane_channel_retau395.cfg` | channel DNS, `Re_tau = 395` | 129 x 257 x 257 |
| `square_duct_reb4410.cfg` | square duct DNS, `Re_b = 4410` | 129 x 129 x 257 |
| `plane_channel_regression.cfg` | rank-robust fixture for `make smoke-driven-periodic` | 33^3 |

## Two Things That Will Bite You

**`first_cell_*` is a fraction of the axis length, not a physical length.** The
stretching is computed on a normalized `[0, 1]` axis and then scaled by the
axis extent, so for a wall-normal axis of length `Ly`:

    first_cell_j = (target physical spacing) / Ly

For a channel at a given `Re_tau` that means
`first_cell_j = (y_plus_target / Re_tau) / Ly`. This is the single number to
edit when retargeting `y+`; `stretch_j` is only the initial guess for the tanh
fit, and the generator reports the factor it actually solved for.

**Cell counts must stay odd at every multigrid level.** Each level coarsens as
`IM -> (IM+1)/2`, so the usable ladder is
`5 -> 9 -> 17 -> 33 -> 65 -> 129 -> 257`. An even count still runs but logs
`can't be consistently coarsened further` and proceeds on a slightly misaligned
coarse grid.

## How To Use

Generate a grid directly:

```bash
python3 generators/grid.gen --config config/grids/coarse_square_tube_curved.cfg
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
- A `config setting '<key>' is not a '<type>' option and was ignored` warning means a
  setting in your `.cfg` never reached the generator and its built-in default stayed in
  force. Fix the key rather than assuming the value took effect.

## Related Docs

- https://vishalkandala.me/picurv-docs/48_Grid_Generator_Guide.html
- https://vishalkandala.me/picurv-docs/07_Case_Reference.html
