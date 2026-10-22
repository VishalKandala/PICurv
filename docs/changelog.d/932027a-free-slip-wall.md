- Experimental: boundary handler `free_slip` for type `WALL` (zero normal velocity, zero shear; no wall
  function or wall eddy-viscosity damping), supported by every momentum solver including Newton-Krylov. The
  `SYMMETRY` boundary type is retired: a symmetry plane is a free-slip wall, and a half channel with a
  free-slip centreline reproduces the full channel to round-off.
