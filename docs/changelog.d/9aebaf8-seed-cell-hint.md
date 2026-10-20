- Experimental: freshly seeded particles start their first location search in the cell they
  were placed in, rather than searching from scratch. This applies to `Volume`, to `Surface`
  and `SurfaceEdges` on ranks that hold the inlet, and to `PointSource`, whose particles go
  directly to the rank owning the point. Particle locations are unchanged. The initial settlement
  of 262,144 Volume particles on a 64^3 box dropped from 192 s to 18 s (debug build).
