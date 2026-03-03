# PICurv Solver Profiles Library

## 1. Purpose

This directory stores reusable `solver.yml` profiles for numerical strategy.

Solver profiles control:
- operation mode (`solve`/`load`/`analytical`),
- momentum solver selection and tolerances,
- solver-specific blocks (for example `dual_time_picard_rk4`),
- pressure/multigrid configuration,
- PETSc passthrough options.

## 2. Included Profiles

- `Imp-MG-Standard.yml`: baseline implicit multigrid-oriented setup.

For full schema coverage, use:
- `examples/master_template/master_solver.yml`

## 3. How to Use

```bash
./bin/picurv run --solve -n 8 \
  --case my_study/case.yml \
  --solver config/solvers/Imp-MG-Standard.yml \
  --monitor config/monitors/Standard_Output.yml
```

## 4. Notes

Use `strategy.momentum_solver` with the exact values `Explicit RK4` or `Dual Time Picard RK4`.
