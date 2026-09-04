# Grid Generator Config Library

Reusable `generators/grid.gen` configuration files. Use these when you want geometry
generated deterministically from parameters rather than from a manually curated grid file.

## Shipped Profiles

| Config | Geometry | Case | Cells |
|---|---|---|---|
| `plane_channel_laminar.cfg` | `box` | laminar exact verification | 17 x 33 x 17 |
| `plane_channel_regression.cfg` | `box` | rank-robust fixture for `make smoke-driven-periodic` | 33^3 |
| `plane_channel_les_retau180.cfg` | `box` | constant-Smagorinsky LES, `Re_tau = 180` | 33 x 65 x 33 |
| `plane_channel_retau180.cfg` | `box` | channel DNS, `Re_tau = 180` | 129^3 |
| `plane_channel_retau395.cfg` | `box` | channel DNS, `Re_tau = 395` | 129 x 257 x 257 |
| `square_duct_reb4410.cfg` | `box` | square duct DNS, `Re_b = 4410` | 129 x 129 x 257 |
| `isotropic_box_64.cfg` | `box` | decaying isotropic turbulence | 64^3 |
| `coarse_square_tube_curved.cfg` | `sweep` | bent square duct | 64 x 64 x 280 |

A `box` with no wall lists is an exactly Cartesian block, which is what every channel and
duct profile above is. Shaping a wall is what turns the same geometry into a step, a hill
or a rib array; see the guide linked below.

## Three Things That Will Bite You

**`first_cell_*` is a fraction of the axis length, not a physical length.** For a
wall-normal axis of length `Ly`:

    first_cell_j = (target physical spacing) / Ly

which for a channel at a given `Re_tau` means
`first_cell_j = (y_plus_target / Re_tau) / Ly`.

You no longer have to trust that arithmetic. Give the generator the reference scales and
it reports the `y+` the spacing actually realizes:

```bash
python3 generators/grid.gen --config config/grids/plane_channel_retau180.cfg box \
        --re-tau 180 --length-ref 1.0 --nu 3.5714286e-04 --stats-file grid.info
```

The `.info` then carries `First_Cell_j_Start_Plus`, the per-axis maximum spacing in wall
units, and each axis extent in wall units. When `picurv` drives the generator,
`length_ref`, `nu` and `velocity_ref` come from the case automatically; `re_tau` is a
design target and belongs in the config. Wall units are labelled nominal because `Re_tau`
is a flow outcome, not a measurement.

**Cell counts must stay odd at every multigrid level.** Each level coarsens as
`IM -> (IM+1)/2`, so `L` levels need node counts of the form `2^(L-1)*k + 1`. Declare
`mg_levels` and the generator refuses a count that cannot coarsen that far, instead of the
solver logging `can't be consistently coarsened further` much later and proceeding on a
misaligned coarse grid.

**Segment lengths must sum to the axis extent.** A wall or path list that stops short
leaves the remainder undefined, and silently rescaling it would move geometry you did not
ask to move. The generator says so rather than guessing.

## How To Use

```bash
python3 generators/grid.gen --config config/grids/plane_channel_retau180.cfg box
```

Through `case.yml`, with `grid.mode: grid_gen` and `grid.generator.config_file`. Output
destinations are not configurable: `picurv` writes the grid, its `.vts` preview and its
`.info` report into the run's own asset store.

## Reproducibility Strategy

Point `grid.generator.config_file` at either a study-local file (for a frozen, reproducible
case snapshot) or a shared library file here (while iterating). Develop against the shared
file, then copy it beside the case once the study stabilizes.

## Common Failure Signals

- `Max_Non_Orthogonality_deg` high in the `.info`: a corner or bend is too sharp for the
  cell count crossing it. Lengthen the segment or refine the axis.
- `Right_Handed = False`: a mirror or an odd axis permutation inverted the grid.
  `Jacobian_Sign_Consistent` cannot see this, which is why the two are separate.
- A segment reported as spanning under four cells: representable but coarse.
- Unexpected domain extents usually mean a stale config value or the wrong config path.
- A rejected config key means the setting never reached the generator, so its built-in
  default would otherwise have stayed in force.

## Related Docs

- https://vishalkandala.me/docs/picurv/48_Grid_Generator_Guide.html
- https://vishalkandala.me/docs/picurv/07_Case_Reference.html
