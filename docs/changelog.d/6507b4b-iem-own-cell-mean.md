- IEM micromixing now relaxes each particle's `Psi` toward its own cell's mean. It
  previously relaxed every particle toward zero: the update read the cell mean from a
  ghosted copy that nothing refreshed after the scatter, and indexed it one cell off from
  the scatter's storage convention. The effect was invisible while every particle started
  at `Psi = 0`, and appears with the first configured initial value.
