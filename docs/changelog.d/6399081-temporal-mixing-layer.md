- Experimental: example `temporal_mixing_layer`, the 2D non-reacting mixing layer of Colucci et al. (1998), the
  particle FDF benchmark, with Picard and Newton Krylov profiles. Promoted to supported with it: Brownian increments
  keyed on seed, particle ID and step, verified on the cluster for diffusive restarts and across rank counts, and the
  restart staging check of carried statistics windows, with its rollback of a refused restart and the grid-input
  warning, verified on the cluster on the login node and through Slurm.
