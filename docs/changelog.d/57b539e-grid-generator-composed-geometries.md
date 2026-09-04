- `generators/grid.gen` rebuilt around two composed geometries, `box` and `sweep`,
  replacing the three hardcoded ones (`cpipe`, `pipe`, `warp`).
  - `box` is a Cartesian block whose two bounding walls may each be a piecewise height
    field - `flat`, `step`, `ramp`, `arc`, `sine`, `gaussian`, `hill` segments laid end to
    end, each continuing from the height the previous reached. One representation reaches
    the backward- and forward-facing step, sudden expansion and contraction, open cavity,
    rib-roughened channel, wall-mounted obstacle, periodic hill, wavy wall, ramp and
    diffuser - no multi-block or immersed-boundary support was added; all of this is a
    single warped box. A spanwise envelope scales a wall's departure from its datum, which
    is what gives a wall-attached obstacle finite width. `box` subsumes `warp` exactly:
    every shipped `[warp]` config section regenerates byte-identical under `[box]` with
    nothing but the section renamed, which is the proof its migration was faithful.
  - `sweep` carries a cross-section (`rectangle` or the square-to-disc `circle` map,
    chosen over an O-grid because that needs a circumferential seam the runtime can only
    reach through an unselectable `cgrid` branch) along a piecewise centreline of
    `straight` and `arc` segments, using a parallel-transported frame so a path leaving
    one plane does not twist the section about its own tangent. It reproduces the retired
    `cpipe` geometry to 1.8e-15 over a domain of size 13 (arithmetic reordering, not a
    difference in geometry) and the 2.9 MB `.picgrid` `examples/bent_channel` used to ship
    to 5e-8, the precision that file was written at.
  - A geometry-independent placement and similarity transform list -
    `anchor`/`translate`/`scale`/`rotate`/`mirror`/`permute` - composes on top of the
    existing `origin` placement, in one `kind:key=value` grammar shared by wall segments,
    path segments and transforms alike, reading identically on the command line and in a
    `.cfg`.
  - Given the case's reference scales (or explicit `--length-ref`/`--nu`/`--velocity-ref`/
    `--re-tau`/`--u-tau`), the `.info` report now also states the grid in solver units and
    nominal wall units - reproducing a shipped channel config's hand-computed `y+ = 0.5`,
    `dx+ = 8.8`, `dz+ = 17.5` instead of leaving that arithmetic in an unchecked comment.
    Coordinates are never scaled by the generator; non-dimensionalization stays exactly
    where it already was.
  - New refusals rather than silent bad output: a wall segment turning the wall in less
    than one cell (below that a shorter segment produces the identical single-cell jump),
    node counts the multigrid hierarchy cannot coarsen, one-sided stretching on a periodic
    axis, and rotating a grid with a periodic axis (the runtime's periodic reconstruction
    assumes an axis-aligned seam). An unknown config key is now an error, not a warning -
    it previously left the built-in default silently in force.
  - Two pre-existing defects the rewrite surfaced: a single-node axis returned a bare
    array where every other path returns `(positions, factor)`, breaking the caller's
    unpack; and `Jacobian_Sign_Consistent` cannot detect a uniformly inverted grid, since
    such a grid agrees with itself. A separate `Right_Handed` report catches it, and
    showed that the retired `pipe` O-grid had always emitted a left-handed mesh - not
    fatal, since the solver already repairs a uniformly negative Jacobian at runtime, but
    previously invisible in the generator's own output.
  - `generators/` was invisible to `audit_family_census` and the review-packet's change
    router, so a change to any of the five bundled generators named no test target and no
    owning page; both now see it, and `grid.gen` has the freshness surface that lets
    `48_Grid_Generator_Guide` be flagged stale when the generator moves without it.
  - `examples/bent_channel` now generates its grid instead of shipping it, and
    `examples/decaying_isotropic_turbulence` gets its own periodic-cube config instead of
    overriding every field of a curved-tube template through `cli_args`. Together with
    deleting a duplicate and a derived `.vts` that should never have been tracked,
    `examples/` drops from 9.1 MB to 3.5 MB.
