- `P_nodal` and `Psi_nodal` are no longer diluted toward zero on the domain boundary:
  pressure's and the particle scalar's dummy cells now repeat the adjacent cell (zero
  normal gradient), filled by `UpdateDummyCells`, which now takes the field it fills.
  Analytical flows now synchronize periodic dummy cells as the solved flow does, which
  corrects periodic `UNIFORM_FLOW` and `TGV3D` velocity dummies. Page 57 records what a
  pressure boundary condition would involve.
