- Experimental post-processing statistics ensembles and profiles: `post.yml -> field_statistics.ensembles`
  derives the moments of the particles within each cell (`source: particles`, fields `Psi` and
  `velocity`) or of a snapshot field (`source: eulerian`), and `field_statistics.homogeneous` writes
  per-step profile CSVs merged over homogeneous directions for windows and ensembles alike.
