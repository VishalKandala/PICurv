# Build Config Guide

This directory stores Makefile configuration fragments for different build environments. These fragments let PICurv keep one build workflow while still supporting local workstations and cluster/HPC toolchains with different compiler/MPI settings.

## Files and Intent

- `config.local.mk`: local workstation defaults for day-to-day development.
- `config.cluster.mk`: cluster/HPC-oriented defaults where modules/toolchains differ from local.
- `config.petsc.mk`: PETSc path/compiler integration details shared by both environments.

## Selection and Invocation

Select a profile with `SYSTEM=`:

```bash
make SYSTEM=local
make SYSTEM=cluster
```

The same profile selection is reachable through `picurv build` wrappers, so users do not need to invoke Make directly unless debugging build internals.

## Practical Guidance

- Keep compiler family compatible with your PETSc build (`gcc` PETSc with `gcc` project toolchain, etc.).
- Keep MPI wrappers consistent across PETSc and PICurv compilation.
- Avoid host-specific absolute paths in shared configs; use environment variables when possible.

## Maintenance Checklist

1. Document any new `SYSTEM=` target in top-level docs.
2. Verify `make doctor` on each supported profile.
3. When changing flags, confirm both `simulator` and `postprocessor` still build.
4. Update CI or cluster setup notes if profile assumptions changed.
