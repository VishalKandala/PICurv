- A fresh particle run now scatters the particles' initial `Psi` to the Eulerian mean
  before the first step, as a restart already did; the first IEM update previously relaxed
  toward a mean that had never seen the initial state. `scalar_transport.iem_constant: 0`
  is accepted and switches micromixing off, carrying each particle's `Psi` bit-identically,
  and `search_metrics.csv` gains `lost_psi_sum`, the `Psi` removed with lost particles each
  step.
