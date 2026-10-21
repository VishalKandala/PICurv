- Experimental: `interpolation.method: FaceFlux` builds particle velocity from the host cell's face
  fluxes, so it is divergence-free where the solver's fluxes are and keeps tracer density uniform;
  it is first-order accurate. `Contra2Cart` now shares its cell solve with it, with unchanged results.
