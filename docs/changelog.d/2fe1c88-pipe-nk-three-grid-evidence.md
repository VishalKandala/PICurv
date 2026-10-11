- Recorded `pipe-poiseuille-nk-three-grid-2026-10-10`: the Newton-Krylov solver reproduces
  Hagen-Poiseuille flow at second order on the swept circle with 9, 17 and 33 cells across.
  It is the solver's first analytical evidence and also backs the `sweep` generator and the
  `circle` cross-section. The dual-time solver did not converge on the 33-cell grid at the
  same time step, which page 24 now states.
