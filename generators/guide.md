# Generators Guide

This directory contains standalone, extensible generators used by PICurv
workflows. Each generator can be invoked directly or selected through PICurv
configuration.

- `grid.gen`: structured PICGRID generation. Two composed geometries rather than a list
  of named shapes: `box` is a Cartesian block whose walls are piecewise height fields,
  and `sweep` carries a cross-section along a piecewise centreline. Wall fields, paths
  and the placement transforms share one
  `kind:key=value` list grammar that reads the same on the command line and in a `.cfg`.
  Given the case's reference scales it also reports the grid in solver and nominal wall
  units, and refuses meshes the runtime cannot use: multigrid-illegal node counts,
  a periodic seam joining unequal cells, and rotating a periodic grid.
- `profile.gen`: dimensional PICSLICE generation and field slicing.
- `ic.gen`: the initial-condition generators. `generator: spectral_velocity` writes a
  spectral random field in a periodic box, or with `wall_axes` a discrete-curl
  wall-compatible perturbation for a channel or duct; `generator: stability_mode_velocity`
  seeds a linear-stability mode. `generator: resampled_velocity`
  imports an external periodic field onto the grid.
- Every file here is a runnable `.gen` script. The code they share is in `picurv_cli/`:
  `periodic_spectral.py` (Fourier machinery for uniform periodic axes),
  `petsc_binary.py` (the PETSc binary Vec reader and writer, at the PetscInt width
  `--petsc-index-bytes` names), `picgrid.py` (the PICGRID grid format and reader) and
  `source_expressions.py` (the expression language);
  `periodic_spectral.py` also owns `require_numpy()`, the lazy NumPy import every script
  uses. Each script puts that directory on `sys.path` when it needs one.
- `plot.gen`: normalized scalar-history, iterative-convergence, and spectrum
  report rendering.

Repository defaults resolve these canonical paths directly.

- `spectra.gen`: shell spectra for periodic boxes, or actual selected plane/line
  FFTs on Cartesian grids; the conductor owns checkpoint discovery and output paths.
- `reduce.gen`: application reducers over field-statistics profiles - a channel's
  wall-unit profile, a duct's cross-section and bisectors, and per-step profile
  functionals (thicknesses, extremes). Reads only small CSVs and the PICGRID; the
  conductor's `reducers` post stage supplies viscosity, scales, window bounds and
  metrics paths from the run.
