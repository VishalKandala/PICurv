- Added experimental `channel_spectral_velocity` and `duct_spectral_velocity`
  initial-condition generators for Cartesian channels with two no-slip walls and
  rectangular ducts with four no-slip walls. Seeded spectral perturbations use a
  discrete curl, support stretched wall-normal grids, and preserve the requested
  bulk flux. They reuse the existing file-based velocity initialization and asset
  staging paths; no C solver changes or recurring perturbation injection are needed.
- Extended `spectra.gen` and the conductor with `plane_spectrum` and `line_spectrum`
  tasks for initial fields and saved checkpoints. Each task transforms one selected
  physical plane or line along uniform periodic axes, retaining signed wavenumbers,
  Cartesian component energies, sample coordinates, and a Parseval check. Parallel
  samples are not averaged; shell spectra remain available for periodic-box cases.
- Replaced the turbulent-channel example's coarse wall-modelled setup with a
  nominal Re_tau 180 candidate DNS startup using `grid.gen` and a stretched
  128 x 128 x 256-cell grid. Added a spectral duct startup configuration and updated
  the capability registry, generated inventories, and usage documentation. These
  configurations are experimental: successful initialization and short runtime
  checks do not establish transition, sustained turbulence, or grid convergence;
  full-scale cluster validation remains outstanding.
