@page 48_Grid_Generator_Guide Grid Generator Guide: generators/grid.gen

@anchor _Grid_Generator_Guide

`generators/grid.gen` is the standalone structured-grid utility used by the `grid.mode: grid_gen` workflow.
This page documents what it can generate, how its config files work, and how `picurv` wraps it.

@tableofcontents

@section p48_role_sec 1. What `grid.gen` Is For

Use `grid.gen` when you want:

- a reproducible Python-side grid-generation step before solver launch,
- generated `.picgrid` geometry instead of C-side `programmatic_c` geometry,
- reusable generator configs that can be versioned and shared across cases,
- optional quick mesh diagnostics (`.info`) and visualization output (`.vts`).

This is distinct from:

- `grid.mode: programmatic_c`: geometry is generated inside the C runtime,
- `grid.mode: file`: an existing `.picgrid` file is staged and used directly.

@section p48_invocation_sec 2. Direct CLI Usage

General form:

```bash
python3 generators/grid.gen [--config <cfgfile>] <geometry> [options...]
```

```bash
python3 generators/grid.gen --config config/grids/plane_channel_retau180.cfg box
python3 generators/grid.gen --config config/grids/coarse_square_tube_curved.cfg sweep --ncells-k 200
python3 generators/grid.gen box --bounds-y 0 2 --stretch-j 2.0 --first-cell-j-start 1.4e-3
python3 generators/grid.gen sweep --cross-section circle --side-lengths 1 1 \
        --path straight:len=5 arc:radius=4,deg=180 straight:len=5
```

`grid.gen` accepts both a config file (`--config`) for reusable defaults and direct CLI
flags for overrides. CLI values override config-file values. A config key matching no
option is an error, not a warning: it would otherwise leave the built-in default in force,
so the file would say one thing and the grid be another.

@subsection p48_option_matrix_ssec 2.1 Option Families

Common to every geometry:

- `--config <file>`, `--output <picgrid_path>`, `--vts <vtk_path>`, `--stats-file <info_path>`
- `--origin X Y Z` - places the bounding-box minimum corner
- `--transforms OP...` - see @ref p48_cap_xform_sec
- `--show-stats` / `--no-show-stats`, `--write-vtk` / `--no-write-vtk`

Stretching, on every axis:

- `--stretch-i`, `--stretch-j`, `--stretch-k`
- `--first-cell-{i,j,k}-{start,end}` - **a fraction of the axis length, not a physical
  length**. See @ref p48_units_sec.

Validation, on every geometry:

- `--mg-levels N` - refuse node counts the multigrid hierarchy cannot coarsen N times
- `--periodic AXIS...` - refuse a seam joining unequal cells, and refuse a rotation

Reference scales, reporting only (see @ref p48_units_sec):

- `--length-ref`, `--nu`, `--velocity-ref`, `--re-tau`, `--u-tau`

Geometry-specific:

- `box`: `--ncells-i`, `--ncells-j`, `--ncells-k`, `--bounds-x`, `--bounds-y`, `--bounds-z`,
  `--wall-j-lo`, `--wall-j-hi`, `--wall-j-lo-span`, `--wall-j-hi-span`,
  `--amp-A`, `--amp-B`, `--amp-C`
- `sweep`: `--ncells-i`, `--ncells-j`, `--ncells-k`, `--cross-section`, `--side-lengths`,
  `--path`, `--cross-section-scale`

Notes for CFD users:

- `ncells-*` are cell counts at input.
- PICGRID header dimensions are node counts (`cells + 1`) after conversion.
- this conversion behavior is regression-tested in `tests/test_cli_smoke.py`.

@section p48_types_sec 3. Supported Geometries

- `box`: a Cartesian block whose two bounding walls may each be shaped by a piecewise
  height field. See @ref p48_cap_geom_box.
- `sweep`: a cross-section carried along a piecewise centreline. See @ref p48_cap_geom_sweep.

Geometry comes from composition rather than from a long list of named shapes. `box` with a
step wall is a backward-facing step; the same wall repeated is a rib-roughened channel;
`sweep` with a 180-degree arc is a U-bend, and with a scale field a nozzle. What that
reaches, and what it does not, is in @ref p48_cap_geom_box and @ref p48_cap_geom_sweep.

@subsection p48_grammar_ssec 3.1 The Ordered-List Grammar

Wall fields, centreline paths and transforms are all ordered lists written the same way,
in the config file and on the command line alike:

```
kind:key=value,key=value
```

whitespace-separated. A config file may spread one across several indented lines:

```ini
wall_j_lo =
    flat:len=4.0,y=1.0
    step:len=0.2,dy=-1.0
    flat:len=15.8
```

The kinds each list accepts are enumerated in @ref p48_cap_wall_sec, @ref p48_cap_path_sec
and @ref p48_cap_xform_sec.

@section p48_config_sec 4. Config File Model

Reusable config files live in `config/grids/`, or beside a study for a reproducibility
snapshot. The file is INI-style:

- the section name is the geometry (`[box]`, `[sweep]`)
- values in that section provide defaults
- CLI arguments override them

Shipped profiles:

| Config | Geometry |
|---|---|
| `plane_channel_laminar.cfg`, `plane_channel_retau180.cfg`, `plane_channel_retau395.cfg`, `plane_channel_les_retau180.cfg`, `plane_channel_regression.cfg` | `box`, plane channel |
| `square_duct_reb4410.cfg` | `box`, square duct |
| `isotropic_box_64.cfg` | `box`, periodic cube |
| `coarse_square_tube_curved.cfg` | `sweep`, bent square duct, 90 deg |
| `bent_channel_coarse.cfg` | `sweep`, `examples/bent_channel`'s bent duct |

@section p48_outputs_sec 5. Outputs

- `.picgrid`: solver-ingestible structured grid file
- `.info`: grid quality and units report
- `.vts`: visualization helper for ParaView

When `picurv` drives the generator it writes all three unconditionally, into the run's own
asset store. Their destinations are not configurable: `grid.generator.output_file`,
`stats_file` and `vts_file` are rejected.

@subsection p48_quality_ssec 5.1 Reading the Quality Report

Two entries decide whether a shaped mesh is usable:

- `Max_Non_Orthogonality_deg` is the cost of every corner in the geometry. A step corner
  spanning one cell reports about 84 degrees, four cells 74, ten 56, twenty 37.
- `Right_Handed` is reported separately from `Jacobian_Sign_Consistent`, because
  consistency only says every cell agrees on a sign, not which sign; a uniformly inverted
  mesh satisfies it. Mirroring or an odd axis permutation is how one gets made.

The generator refuses a wall segment spanning less than one cell, because a shorter
segment then produces the identical single-cell jump: the refinement asked for would not
exist. Between one and four cells it reports and proceeds.

@section p48_units_sec 5.2 Units In The Report

`grid.gen` writes a **dimensional** grid. The launcher divides coordinates by
`length_ref` once, in `validate_and_nondimensionalize_picgrid`, for every grid whatever
its origin - the generator never scales coordinates itself.

Given reference scales it reports the grid in three sets of units:

- as written, dimensional;
- in **solver units**, if `--length-ref` is supplied - the grid the solver receives;
- in **nominal wall units**, if `--re-tau` (with `--length-ref`) or `--u-tau` (with
  `--nu`) is supplied.

Wall units are nominal by construction: `Re_tau` is a flow outcome, not an input, so what
is reported is what the supplied target implies, not a measurement. When `picurv` drives
the generator, `length_ref`, `nu` and `velocity_ref` are taken from the case
automatically; `re_tau` is a design target and belongs in the generator config.

This is also the way to stop computing `y+` by hand. `first_cell_*` is a fraction of the
axis length, so a wall-normal axis of length `Ly` needs
`first_cell_j = (y_plus_target / Re_tau) / Ly`. Supply the reference scales and the report
states the `y+` that spacing actually realizes, rather than leaving the arithmetic in a
comment that nothing checks.

@section p48_case_sec 6. Using `grid.gen` Through `case.yml`

```yaml
grid:
  mode: grid_gen
  generator:
    config_file: config/grids/plane_channel_retau180.cfg
    grid_type: box
    cli_args: ["--ncells-i", "96", "--ncells-j", "96"]
```

Contract notes:

- `grid.generator.config_file` is required.
- `picurv` does not generate a temporary `grid.cfg` for `grid.gen`.
- `grid.gen` input dimensions (`ncells_*`, `--ncells-*`) are cell counts, converted to
  node counts before the `.picgrid` header is written.
- `cli_args` is passed through untouched, but closed-choice values inside it are checked
  at validation time, so a misspelled selector fails before the run rather than partway
  into it.

@section p48_choose_sec 7. When To Use Which Grid Path

Use `programmatic_c` when:

- you want the simplest structured box-style runtime setup,
- you do not need an external generated mesh artifact,
- you want direct C-side control of the structured grid box inputs.

Use `file` when:

- you already have a stable `.picgrid`,
- geometry should be fixed and explicit,
- you want no generator step during launch.

Use `grid_gen` when:

- you want generated geometry but still want file-based staging,
- you want a reusable mesh recipe with controlled overrides,
- you want generated `.info`/`.vts` outputs alongside the mesh.

For a generated grid used with geometric-periodic BCs, each opposite periodic
surface must be a pointwise translated copy of its partner. The translation
must be nonzero and constant over the entire surface; rotating, shearing, and
nonmatching surface pairs are rejected at runtime.

@section p48_cap_geom_sec 7.1 Grid Generator Geometry Entries

@htmlinclude generated/capability_inventory_grid_generator_type.html

@subsection p48_cap_geom_box_sub box

@anchor p48_cap_geom_box

**Identity.** `grid.gen box`, or `grid.generator.grid_type: box`.

**What it does.** Generates a Cartesian block whose two bounding walls may each be shaped. The wall-normal coordinate blends between them, `y = y_lo(x, z) + eta * (y_hi(x, z) - y_lo(x, z))`, so wall clustering stays a fraction of the *local* height rather than of the declared bounds. Each wall is a piecewise height field built from segments laid end to end, each continuing from the height the previous reached.

**When to choose it.** Every box-like domain: plain channels and ducts, and anything whose walls are shaped. Because the wall is a list rather than a named shape, one geometry covers the backward- and forward-facing step, sudden expansion and contraction, open cavity, rib-roughened channel, wall-mounted obstacle, periodic hill, wavy wall, ramp and diffuser.

**Parameters it owns.** `--ncells-i`, `--ncells-j`, `--ncells-k`, `--bounds-x`, `--bounds-y`, `--bounds-z`, the wall fields `--wall-j-lo` and `--wall-j-hi`, their spanwise envelopes `--wall-j-lo-span` and `--wall-j-hi-span`, and the sinusoidal warp amplitudes `--amp-A`, `--amp-B`, `--amp-C`.

**Interactions.** Reads defaults from the `[box]` config section; `config/grids/` ships plane-channel, square-duct and isotropic-box profiles. Segment lengths must sum to the axis extent, because a wall that stops short leaves its remainder undefined. `origin` places the *realized* bounding box, so shaping a wall shrinks the block and moves it relative to what `bounds_y` alone suggests.

**Diagnostics.** A segment that turns the wall in less than one cell is refused: below that, a shorter segment produces the identical single-cell jump, so the refinement asked for would not exist. Segments under four cells are reported. `Max_Non_Orthogonality_deg` in the `.info` report carries the cost - measured on a 200-cell axis, a step spanning one cell reports 84.2 degrees, four cells 73.7, ten cells 55.7, twenty cells 36.6.

**Evidence.** Production exercised - the `examples/periodic_test` driven-channel and driven-duct cases generate their meshes with this type, and all of them reproduce byte-identically the meshes the retired `warp` type produced.

**Limitations.** The wall is a single-valued function of position by construction, so it cannot express an overhang, a detached body, or anything the flow passes underneath; those need an immersed boundary or multiple blocks, and neither is available in the runtime today. Walls may be shaped on one axis pair only: profiling two intersecting pairs at once needs three-dimensional transfinite interpolation, which is not implemented. No solve has yet been run on a stepped mesh, so metric quality at a resolved corner is reported but not validated. Experimental until a case exercises the shaped-wall path end to end.

@subsection p48_cap_geom_sweep_sub sweep

@anchor p48_cap_geom_sweep

**Identity.** `grid.gen sweep`, or `grid.generator.grid_type: sweep`.

**What it does.** Sweeps a cross-section along a piecewise centreline. Shape and path are independent, so any section travels any path. The frame is carried by parallel transport, which means each arc rotates the whole frame about its own turning axis and no twist accumulates about the tangent.

**When to choose it.** Ducts and pipes: straight or bent, single-plane or not. A U-bend is one arc of 180 degrees, an S-bend two arcs of opposite sign, and a nozzle or diffuser is any of these with a scale field along the path.

**Parameters it owns.** `--ncells-i`, `--ncells-j`, `--ncells-k`, `--cross-section`, `--side-lengths`, `--path`, `--cross-section-scale`.

**Interactions.** Reads defaults from the `[sweep]` config section; `config/grids/coarse_square_tube_curved.cfg` ships a bent square duct. Logical axes follow the two cross-section directions and then the path, which makes the result right-handed. An `arc` turns in the plane perpendicular to its axis and may not turn about an axis parallel to the current tangent, because an arc has no defined plane there; chain arcs about different axes to leave a plane.

**Diagnostics.** The `.info` report carries `Right_Handed` alongside `Jacobian_Sign_Consistent`, which a uniformly inverted mesh would satisfy. On an out-of-plane path the parallel-transported frame reports 0.6 degrees maximum non-orthogonality, where a frame built from a fixed up-vector would shear the section badly.

**Evidence.** Production exercised - `examples/bent_channel` generates its square duct and quarter turn with this type. It reproduces the 2.9 MB `.picgrid` that example used to ship to 5e-8, which is the precision that file was written at, and the mesh the retired `cpipe` type produced from the same parameters to 1.8e-15 over a domain of size 13 - arithmetic reordering rather than a difference in geometry.

**Limitations.** One block only. `circle` uses the square-to-disc map rather than an O-grid: `Metric.c` has a `cgrid` branch for a circumferential seam, but it is wired only to `programmatic_c`'s `cgrids` flag, and even there `src/grid.c` never builds anything but a Cartesian box - no path in the runtime today both produces real O-grid coordinates and can seam them. The disc covers the whole section with no hole and no axis singularity, but cell size varies with angle and is worst in the four diagonal regions - a 32x32 section reports 79.7 degrees maximum and 12.1 average non-orthogonality against 0.004 for a rectangle. That is adequate for a cylindrical domain and is not a wall-resolved pipe mesh, which needs a butterfly topology and therefore multiple blocks. Experimental until a case exercises it.


@section p48_cap_xsec_sec 7.2 Swept Cross-Section Entries

@htmlinclude generated/capability_inventory_grid_cross_section.html

@subsection p48_cap_xsec_rectangle_sub rectangle

@anchor p48_cap_xsec_rectangle

**Identity.** `grid.gen sweep --cross-section rectangle`, or `cross_section = rectangle`.

**What it does.** Sweeps the logical square itself, scaled by `--side-lengths`. Grid lines stay parallel to the section axes, so a straight rectangular sweep is an exactly Cartesian block.

**When to choose it.** Square and rectangular ducts, and anything where a circular section would change the mesh rather than the physics. It is the better-conditioned of the two: a 32x32 section reports 0.004 degrees maximum non-orthogonality.

**Parameters it owns.** `--side-lengths I J`, the two full widths.

**Interactions.** Sweeps along any `--path`, and scales along it with `--cross-section-scale`. The logical section axes run on [-0.5, 0.5], so `first_cell_i_*` and `first_cell_j_*` are fractions of the full side length.

**Diagnostics.** `Max_Non_Orthogonality_deg` stays near zero on a straight sweep; anything larger comes from the path, not the section.

**Evidence.** Production exercised - `examples/bent_channel` sweeps a `rectangle` section through its bend, and `config/grids/coarse_square_tube_curved.cfg` ships an unselected second one.

**Limitations.** Corners are square, so a rectangular duct has four concave corner lines where the boundary layers of two walls meet; that is the geometry, not a mesh defect, but it is where near-wall resolution has to be judged on both walls at once.

@subsection p48_cap_xsec_circle_sub circle

@anchor p48_cap_xsec_circle

**Identity.** `grid.gen sweep --cross-section circle`, or `cross_section = circle`.

**What it does.** Maps the logical square onto a disc by `x = u*sqrt(1 - v^2/2)`, `y = v*sqrt(1 - u^2/2)`. The square's edges land exactly on the circle and its corners on the 45-degree points, so one block covers the whole disc with no hole and no axis singularity. `--side-lengths` gives the two diameters.

**When to choose it.** Cylindrical and round-duct domains at laminar or moderate Reynolds number, where what matters is that the domain is round.

**Parameters it owns.** `--side-lengths I J`, read as the two diameters.

**Interactions.** The logical section axes run on [-1, 1] rather than [-0.5, 0.5], so `first_cell_i_*` and `first_cell_j_*` are fractions of the full diameter. Wall clustering thickens cells toward the boundary but not uniformly around it.

**Diagnostics.** `Max_Non_Orthogonality_deg` and `Max_Volume_Expansion_Ratio_Global` both report the diagonal-region cost; on a 32x32 section they read 79.7 degrees and 2.81 against 0.004 and 1.0 for a rectangle.

**Evidence.** Implemented only. Unit tests check that the map fills the disc, reaches its edge and leaves no hole; no shipped case selects it.

**Limitations.** Cell size varies with angle and is worst in the four diagonal regions: a 32x32 section reports 79.7 degrees maximum and 12.1 average non-orthogonality, against 0.004 for a rectangle. Composing a radial redistribution transform evens out wall-normal spacing but not azimuthal spacing. This is not a wall-resolved pipe mesh - that needs an O-grid annulus stitched to a core block, which is multi-block, and the runtime has no block-interface handler. The alternative single-block topology, an O-grid, is not offered: `Metric.c` has a `cgrid` branch for its circumferential seam, but nothing in the runtime today can both build real O-grid coordinates and reach that branch - it is wired only to `programmatic_c`, whose grid generator never builds anything but a Cartesian box.

@section p48_cap_wall_sec 7.3 Wall Height Field Segment Entries

@htmlinclude generated/capability_inventory_grid_wall_segment.html

A wall is these laid end to end along the streamwise axis. Every segment takes `len`, and every segment starts at the height the previous one reached, so a wall is continuous by construction. Lengths must sum to the axis extent. The first segment may carry `y` to set the datum.

@subsection p48_cap_wall_flat_sub flat

@anchor p48_cap_wall_flat

**Identity.** `flat:len=L[,y=Y]`.

**What it does.** Holds the current height across its length. The whole wall omitted is equivalent to one `flat` at the bounds datum.

**When to choose it.** The developed stretches between features: an inlet run before a step, the recovery after it, the gap between ribs.

**Parameters it owns.** `len`, and `y` on the first segment only.

**Interactions.** Omitting a wall entirely is equivalent to a single `flat` at the `--bounds-y` datum, which is how every shipped channel and duct configuration uses `box`.

**Diagnostics.** None of its own: a flat stretch contributes no non-orthogonality.

**Evidence.** Production exercised - the `examples/periodic_test` cases are flat-walled boxes.

**Limitations.** None beyond needing a length; a `flat` cannot itself set a height except as the first segment, where `y` establishes the datum.

@subsection p48_cap_wall_step_sub step

@anchor p48_cap_wall_step

**Identity.** `step:len=L,dy=D[,shape=smoothstep|cosine|linear]`.

**What it does.** Changes the wall height by `dy` across its own length, following a shape that is exactly flat at both ends. Sharpness *is* the length: a shorter segment is a sharper corner.

**When to choose it.** Backward- and forward-facing steps, sudden expansions and contractions, and the rising and falling faces of ribs and wall-mounted obstacles.

**Parameters it owns.** `len`, `dy`, and optionally `shape`.

**Interactions.** Length is the only sharpness control, so it interacts directly with `--ncells-i`: refining the axis is what lets a shorter step be represented. Two steps of opposite sign make a rib; a step and a `ramp` make an expansion into a diffuser.

**Diagnostics.** A sub-cell segment is refused by name. Under four cells it is reported to stderr. `Max_Non_Orthogonality_deg` carries the cost.

**Evidence.** Implemented only. Unit tests cover the refusal, the report, and that a shorter corner costs orthogonality; no shipped case selects a shaped wall.

**Limitations.** A segment spanning less than one cell is refused, because below that a shorter length produces the identical single-cell jump and the refinement asked for would not exist. Under four cells the corner is reported as coarse. `Max_Non_Orthogonality_deg` is the number to judge it by: on a 200-cell axis a one-cell corner reports 84.2 degrees, four cells 73.7, ten 55.7, twenty 36.6. A true zero-length corner is not expressible in one block at all.

@subsection p48_cap_wall_ramp_sub ramp

@anchor p48_cap_wall_ramp

**Identity.** `ramp:len=L,dy=D`.

**What it does.** Changes the height by `dy` linearly across its length.

**When to choose it.** Diffusers, contractions and inclined walls, where the wall angle is the quantity of interest and the kinks at each end are acceptable.

**Parameters it owns.** `len` and `dy`.

**Interactions.** Follows and precedes `flat` naturally; between two `arc` segments it forms a straight run with rounded ends.

**Diagnostics.** The end kinks appear as a local `Max_Non_Orthogonality_deg` spike independent of the ramp's own slope.

**Evidence.** Implemented only. Unit tests cover its continuity with neighbouring segments.

**Limitations.** Slope is discontinuous at both ends, unlike `step` and `arc`. That kink is a local non-orthogonality spike; use `arc` where the transition itself matters.

@subsection p48_cap_wall_arc_sub arc

@anchor p48_cap_wall_arc

**Identity.** `arc:len=L,dy=D[,shape=cosine|smoothstep|linear]`.

**What it does.** The same height change as `step`, defaulting to a cosine shape - a rounded transition rather than a squared-off one.

**When to choose it.** Where the curvature of the transition is part of the geometry rather than an artefact of meshing it: a filleted expansion, a rounded shoulder.

**Parameters it owns.** `len`, `dy`, and optionally `shape`.

**Interactions.** Shares the resolution floor with `step`, since both turn the wall across their own length. Choosing `shape=cosine` on a `step` and `shape=smoothstep` on an `arc` makes the two interchangeable.

**Diagnostics.** Same floor and reporting as `step`.

**Evidence.** Implemented only. Covered by the shared segment tests.

**Limitations.** Shares `step`'s resolution floor and reporting, and is symmetric about its midpoint: an asymmetric fillet needs two segments.

@subsection p48_cap_wall_sine_sub sine

@anchor p48_cap_wall_sine

**Identity.** `sine:len=L,amp=A[,cycles=N]`.

**What it does.** Oscillates about the entry height with amplitude `amp` over `cycles` full periods across the segment.

**When to choose it.** Wavy-wall channels, and any periodic corrugation where the wall is a single harmonic.

**Parameters it owns.** `len`, `amp`, and optionally `cycles`.

**Interactions.** For a streamwise-periodic wavy channel, `cycles` must be a whole number and the domain length must match the wavelength, or the seam will not close.

**Diagnostics.** `Max_Non_Orthogonality_deg` grows with amplitude over wavelength; there is no separate check on the corrugation slope.

**Evidence.** Implemented only. No shipped case selects it.

**Limitations.** Returns to the entry height only when `cycles` is a whole number; a fractional count leaves a step at the join with the next segment, which the continuity of the surrounding wall will then carry forward. Amplitude is not checked against the channel height, so a large one can close the domain.

@subsection p48_cap_wall_gaussian_sub gaussian

@anchor p48_cap_wall_gaussian

**Identity.** `gaussian:len=L,amp=A[,width=W,center=C]`.

**What it does.** Adds a Gaussian bump of height `amp` centred at `center` within the segment, with standard deviation `width`.

**When to choose it.** Isolated smooth obstacles and bumps, where the feature should decay rather than end.

**Parameters it owns.** `len`, `amp`, and optionally `width` and `center`.

**Interactions.** Give the segment several `width` of run either side of `center` so the bump has decayed before the join.

**Diagnostics.** The residual height at the segment ends is not reported; check the wall profile in the `--vts` output if the bump is narrow.

**Evidence.** Implemented only. No shipped case selects it.

**Limitations.** A Gaussian never reaches zero, so the wall does not return exactly to its entry height at the segment ends unless the segment is several `width` long; a narrow bump in a short segment leaves a small discontinuity at the join.

@subsection p48_cap_wall_hill_sub hill

@anchor p48_cap_wall_hill

**Identity.** `hill:len=L,height=H`.

**What it does.** A symmetric cos-squared rise and fall of height `height` across the segment. Value and slope both vanish at each end, so hills tile without a kink.

**When to choose it.** Periodic-hill configurations and any repeating smooth obstacle, where the wall must close on itself to be periodic.

**Parameters it owns.** `len` and `height`.

**Interactions.** Tiles without a kink, so a periodic-hill domain is one `hill` per wavelength and the streamwise axis can be declared periodic.

**Diagnostics.** Value and slope vanish at both ends by construction, so the joins contribute nothing to `Max_Non_Orthogonality_deg`.

**Evidence.** Implemented only. A unit test checks that the profile closes on itself; no shipped case selects it.

**Limitations.** This is a cos-squared profile, not the ERCOFTAC periodic-hill polynomial: comparisons against that benchmark's separation and reattachment data are not like-for-like. No shipped case exercises it.

@section p48_cap_path_sec 7.4 Centreline Path Segment Entries

@htmlinclude generated/capability_inventory_grid_path_segment.html

@subsection p48_cap_path_straight_sub straight

@anchor p48_cap_path_straight

**Identity.** `straight:len=L`.

**What it does.** Continues along the current tangent for `len`, carrying the frame unchanged.

**When to choose it.** Inlet and outlet runs, and the developed sections between bends. A path of one `straight` is a prismatic duct.

**Parameters it owns.** `len`.

**Interactions.** Carries the frame unchanged, so it neither creates nor corrects twist. `first_cell_k_*` and `stretch_k` distribute cells along it.

**Diagnostics.** Contributes no non-orthogonality of its own.

**Evidence.** Production exercised - both straight runs of `examples/bent_channel`'s path are this segment.

**Limitations.** None beyond needing a positive length.

@subsection p48_cap_path_arc_sub arc

@anchor p48_cap_path_arc

**Identity.** `arc:radius=R,deg=A[,axis=x|y|z]`.

**What it does.** Turns through `deg` on a constant radius, in the plane perpendicular to `axis` (default `z`). The whole frame rotates with the turn, which is parallel transport: no twist accumulates about the tangent.

**When to choose it.** Any bend. A U-bend is one arc of 180 degrees, an S-bend two arcs of opposite sign, and a path leaves a single plane by chaining arcs about different axes.

**Parameters it owns.** `radius`, `deg`, and optionally `axis`.

**Interactions.** Consecutive arcs about different axes take the path out of a plane, which is where parallel transport differs from a fixed up-vector. A negative `deg` reverses the turn, which is how an S-bend is written.

**Diagnostics.** `Max_Non_Orthogonality_deg` is the check on bend tightness; on an out-of-plane two-arc path it reads 0.6 degrees, so a large value means the radius is too small for the section.

**Evidence.** Production exercised for a single, default-axis arc - `examples/bent_channel`'s quarter turn. Chaining arcs about different axes to leave a plane, a negative `deg`, and the refusal of an arc about its own tangent are unit-tested only; no shipped case exercises any of them.

**Limitations.** Constant radius only, so a variable-radius or spiral bend needs several segments to approximate. An arc may not turn about an axis parallel to the current tangent, which has no defined plane and is refused. Curvature is not checked against the section size, so a bend radius comparable to the duct width will produce badly skewed cells on the inside of the turn; `Max_Non_Orthogonality_deg` is where that shows.

@section p48_cap_xform_sec 7.5 Placement and Similarity Transform Entries

@htmlinclude generated/capability_inventory_grid_transform.html

Transforms apply after the geometry map, in the order given, and compose on top of the `origin` placement rather than replacing it. Every kind takes `kind:key=value,key=value`.

@subsection p48_cap_xform_anchor_sub anchor

@anchor p48_cap_xform_anchor

**Identity.** `anchor:mode=bbox_min|bbox_max|bbox_center|centroid|origin,x=X,y=Y,z=Z`.

**What it does.** Moves the grid so the named reference point lands at the given position.

**When to choose it.** When the natural handle on the domain is not its minimum corner - centring a box on the origin, or putting a duct's centroid somewhere specific.

**Parameters it owns.** `mode`, and `x`, `y`, `z`.

**Interactions.** `origin` already applies `anchor:mode=bbox_min` before any transform list runs, so naming `anchor` again replaces that placement rather than adding to it.

**Diagnostics.** The realized bounding box is in the `.info` report as `X_Range`, `Y_Range` and `Z_Range`; compare against them to confirm the placement.

**Evidence.** Production exercised in its default form - every generated grid is placed by the implicit `bbox_min` anchor.

**Limitations.** The reference is computed from the realized geometry, so shaping a wall changes where `bbox_min` is; `origin` already anchors `bbox_min` by default, and a second `anchor` overrides that rather than composing with it.

@subsection p48_cap_xform_translate_sub translate

@anchor p48_cap_xform_translate

**Identity.** `translate:dx=X,dy=Y,dz=Z`.

**What it does.** Shifts every node by a fixed displacement.

**When to choose it.** Offsetting a grid from wherever placement put it, without reference to its extent.

**Parameters it owns.** `dx`, `dy`, `dz`.

**Interactions.** Composes with anything, and is the transform to reach for when the geometry's extent should not influence where it ends up.

**Diagnostics.** Visible directly in the reported ranges.

**Evidence.** Implemented only. Covered by unit tests for order-dependent composition.

**Limitations.** None; it is the one transform whose effect does not depend on the geometry.

@subsection p48_cap_xform_scale_sub scale

@anchor p48_cap_xform_scale

**Identity.** `scale:sx=X,sy=Y,sz=Z[,about=bbox_min|bbox_center|centroid|origin]`.

**What it does.** Scales each axis independently about the named reference point.

**When to choose it.** Adjusting a generated domain's proportions without editing its bounds, or converting units.

**Parameters it owns.** `sx`, `sy`, `sz`, and optionally `about`.

**Interactions.** Applied after wall shaping and stretching, so it rescales targeted first-cell sizes along with everything else. Supply `--length-ref` if the report should show what the solver receives.

**Diagnostics.** `Max_Aspect_Ratio` and the first-cell entries in the report both move under a non-uniform scale.

**Evidence.** Implemented only. Unit tests cover the reference-point behaviour and the zero-factor refusal.

**Limitations.** A zero factor is refused since it collapses the grid. Non-uniform factors change cell aspect ratios and can undo carefully targeted wall spacing; this is not the way to non-dimensionalize, which the launcher does once by `length_ref`.

@subsection p48_cap_xform_rotate_sub rotate

@anchor p48_cap_xform_rotate

**Identity.** `rotate:axis=x|y|z,deg=A[,about=bbox_min|bbox_center|centroid|origin]`.

**What it does.** Rotates the grid about an axis through the named reference point.

**When to choose it.** Orienting a generated domain to match an external geometry or an inflow direction.

**Parameters it owns.** `axis`, `deg`, and optionally `about`.

**Interactions.** Mutually exclusive with `--periodic`. Compose with `translate` to rotate about a point that is not one of the named references.

**Diagnostics.** `Jacobian_Sign_Consistent` and `Right_Handed` both stay true: a rotation is orientation-preserving.

**Evidence.** Implemented only. A unit test covers the periodic refusal.

**Limitations.** Refused on any grid with a declared periodic axis. The runtime's periodic reconstruction offsets only the matching Cartesian component across the seam and copies the other two, so it assumes an axis-aligned seam; whether an arbitrarily oriented periodic grid is handled correctly is unverified, and this refuses rather than producing one quietly.

@subsection p48_cap_xform_mirror_sub mirror

@anchor p48_cap_xform_mirror

**Identity.** `mirror:plane=xy|yz|zx[,about=...]`.

**What it does.** Reflects the grid across the named plane through the reference point.

**When to choose it.** Producing the opposite-handed version of an asymmetric geometry - a bend that turns the other way - without re-specifying it.

**Parameters it owns.** `plane`, and optionally `about`.

**Interactions.** Two mirrors restore the original handedness, as does composing one with an odd `permute`.

**Diagnostics.** `Right_Handed` flips to false. `Jacobian_Sign_Consistent` does not move, which is exactly why the two are reported separately.

**Evidence.** Implemented only. Covered by the handedness tests.

**Limitations.** Inverts handedness. `Jacobian_Sign_Consistent` cannot detect this, because a uniformly inverted grid agrees with itself; read `Right_Handed` instead. The solver repairs a uniformly left-handed grid by flipping its metric vectors, so this is a reporting concern rather than a failure.

@subsection p48_cap_xform_permute_sub permute

@anchor p48_cap_xform_permute

**Identity.** `permute:order=xyz|xzy|yxz|yzx|zxy|zyx`.

**What it does.** Reassigns which physical axis each coordinate becomes.

**When to choose it.** Putting the streamwise direction on the axis a case expects, without regenerating the geometry. It subsumes the orientation flag the retired bent-pipe geometries carried.

**Parameters it owns.** `order`.

**Interactions.** Replaces the orientation flag the retired bent-pipe geometries carried: reorienting a swept duct is a transform rather than a geometry parameter.

**Diagnostics.** An odd order flips `Right_Handed`; an even one leaves it true, and unit tests pin both.

**Evidence.** Implemented only. Covered by the handedness tests.

**Limitations.** An odd permutation inverts handedness, with the same reporting consequence as `mirror`; an even one does not. It permutes physical coordinates, not logical axes, so the mapping between logical index and wall direction is unchanged.


@section p48_related_sec 8. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 17_Workflow_Extensibility**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@ref p54_geometric_periodic "Periodic Boundaries and Driven Flows"**
