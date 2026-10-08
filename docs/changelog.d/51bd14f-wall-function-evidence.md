- The wall-function models (`log_law`, `werner`, `cabot`) now declare their unit-test
  coverage (`make unit-runtime`, `make unit-boundaries`), so review tooling routes a
  wall-model change to its tests. The Case Reference now documents the first-cell y+
  guard that was already in place: a run stops after ten consecutive steps with the
  wall-face mean y+ outside the selected law's range. The master template's
  wall-function comments now state correctly that `roughness_height` is rejected for
  `werner` and `cabot`, and that `cabot` uses the resolved pressure gradient.
