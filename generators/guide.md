# Generators Guide

This directory contains standalone, extensible generators used by PICurv
workflows. Each generator can be invoked directly or selected through PICurv
configuration.

- `grid.gen`: structured PICGRID generation and legacy-grid conversion. Two composed
  geometries rather than a list of named shapes: `box` is a Cartesian block whose walls
  are piecewise height fields, and `sweep` carries a cross-section along a piecewise
  centreline. Wall fields, paths and the placement transforms share one
  `kind:key=value` list grammar that reads the same on the command line and in a `.cfg`.
  Given the case's reference scales it also reports the grid in solver and nominal wall
  units, and refuses meshes the runtime cannot use: multigrid-illegal node counts,
  a periodic seam joining unequal cells, and rotating a periodic grid.
- `profile.gen`: dimensional PICSLICE generation and field slicing.
- `ic.gen`: expression-driven and configurable spectral-random PETSc
  initial-condition vector generation. Select the latter with
  `generator: spectral_random_velocity`.
- `plot.gen`: normalized scalar-history, iterative-convergence, and spectrum
  report rendering.
- `convert_grid_from_legacy_to_picgrid.py`: legacy conversion command that
  delegates to `grid.gen legacy1d`.

Repository defaults resolve these canonical paths directly.
