- Fixed: the dummy cell centres next to a periodic seam are now the cells across the seam
  shifted by the validated periodic translation, including along periodic edges and at corners,
  where the old spacing-based shift was zero and particle interpolation stopped with
  "Degenerate plane". Periodic cases without particles are unchanged.
