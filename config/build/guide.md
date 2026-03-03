# Build Config Guide

This directory stores Makefile configuration fragments for different build environments.

## Files

- `config.local.mk`: local workstation defaults.
- `config.cluster.mk`: cluster/HPC oriented defaults.
- `config.petsc.mk`: PETSc-related build settings.

## Usage

Select a profile via:

```bash
make SYSTEM=local
make SYSTEM=cluster
```

or equivalent through `picurv build` wrappers.

## Maintenance Notes

- keep compiler/MPI settings consistent with PETSc build,
- avoid hard-coding host-specific paths in shared configs,
- document any new `SYSTEM=` option in top-level build docs.
