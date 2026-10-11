- Analytical `TGV3D` can now run on a file-backed grid (`grid.mode: file`). The field is
  evaluated at the given grid's cell centres instead of the solver replacing the grid with
  its own `[0, 2π]² × [0, 0.2π]` box, which allows particle-interpolation verification on
  curvilinear (for example interior-warped) grids. `programmatic_c` behaves as before;
  `grid_gen` is still rejected.
