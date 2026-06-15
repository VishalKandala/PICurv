# Generators Guide

This directory contains standalone, extensible generators used by PICurv
workflows. Each generator can be invoked directly or selected through PICurv
configuration.

- `grid.gen`: structured PICGRID generation and legacy-grid conversion.
- `profile.gen`: dimensional PICSLICE generation and field slicing.
- `ic.gen`: expression-driven PETSc initial-condition vector generation.
- `plot.gen`: normalized time-history request rendering.
- `convert_grid_from_legacy_to_picgrid.py`: legacy conversion command that
  delegates to `grid.gen legacy1d`.

Repository defaults resolve these canonical paths directly.
