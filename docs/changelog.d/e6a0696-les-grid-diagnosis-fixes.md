- Diagnosing the Humphrey bent-channel failure found and fixed defects in the LES filter
  widths, the Clark gradient term and the metric helpers, added the Vreman and WALE
  subgrid models and the Scotti filter width, and made the grid report state what
  predicts grid-scale oscillation.
  - The failure itself was numerical, not a solver or physics fault: central convection,
    which any LES model selects regardless of `central_diff`, leaves a two-cell
    streamwise oscillation undamped where the streamwise cell Reynolds number is large
    and the subgrid model supplies almost nothing, as in a laminar or transitional
    stretch. Solver reference page 08 now says so.
  - The Clark gradient term multiplied per-cell velocity differences, which already
    carry the cell spacing, by a squared cell extent as well. Its stress scaled as
    `Delta^4` and was smaller than intended by `Delta^2` - 1e-4 to 1e-6 on a production
    mesh - so any earlier run with `gradient_model.enabled` effectively ran without it.
    It now matches `(Delta^2/12) du/dx du/dx`, and such runs will behave differently.
  - The `geometric_mean` and `max_edge` filter widths used the Cartesian components of
    the cell diagonal, so an identical cell got up to nine times the eddy viscosity once
    a bend turned it. They now use the cell's extents along its own grid directions,
    and their widths change on curved grids.
  - A fixed 1e-12 floor on the determinant of the face-area vectors declared any
    wall-resolved cell singular, and the refusal was dropped, so selecting `max_edge` or
    `geometric_mean` on such a grid produced a NaN eddy viscosity. Degeneracy is now
    judged against the cell's own scale and the error propagates. The three helpers
    involved have no production caller left and are marked deprecated.
  - New: `filter_width: scotti` (Scotti, Meneveau and Lilly 1993), and the
    coefficient-free models `model: vreman` and `model: wale`, which vanish in pure shear
    and need no test filter, averaging or coefficient field. Vreman weights each grid
    direction by the cell's own edge. All are experimental and unvalidated against a
    reference flow.
  - The constant model accepts `filter_width` directly instead of through a PETSc
    passthrough. Each LES model now refuses the parameters of the others, including
    `dynamic_frequency` under the constant model, and `-les` is range-checked instead of
    being read into an uninitialized value when absent.
  - The grid generator report states the wall-unit length and what it was built from,
    since an `Re_tau` quoted on a half-width with a full-width `length_ref` halves every
    plus number, and adds the cell Reynolds number, the effective cell Reynolds number
    under an assumed `--nut-ratio`, and the time viscosity takes to clear a two-cell
    oscillation along each direction.
  - Every tracked C source and header now routes to a freshness surface or contract;
    41 of 75 routed nowhere. The test map gave `src/les.c` no `unit-les` target.
