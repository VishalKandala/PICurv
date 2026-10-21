- `grid.gen` no longer refuses `rotate` on a grid with `--periodic` axes: a rigidly rotated
  periodic grid has a seam translation in any direction, which the runtime validates and applies in
  full. Rotational periodicity (a periodic axis that curves) is still unsupported.
