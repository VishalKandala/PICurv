# Configuration Guide

This directory contains reusable YAML profiles and build/runtime configuration assets. It is best treated as a configuration library: stable, versioned building blocks that users compose into case-specific workflows.

For CFD users, the key idea is separation of concerns. Instead of creating one monolithic YAML file, PICurv uses role-oriented contracts (`case`, `solver`, `monitor`, `post`, and optional `cluster`/`study`) so you can change numerical strategy without rewriting geometry definitions, or change post outputs without touching solver controls.

## Sub-guides

- `build/guide.md`
- `grids/guide.md`
- `monitors/guide.md`
- `postprocessors/guide.md`
- `runtime/guide.md`
- `schedulers/guide.md`
- `solvers/guide.md`
- `studies/guide.md`

## How To Use This Library Effectively

1. Start from an initialized example (`./bin/picurv init ...`) so contracts are already valid.
2. Replace only the role file you are actively experimenting with.
3. Validate after each change (`./bin/picurv validate ...`) before launching runs.
4. Promote stable reusable profiles back into `config/` so team workflows converge.

## Common CFD Iteration Patterns

- Change `solver.yml` while freezing `case.yml` to isolate numerics.
- Change `case.yml` grid/BC setup while keeping monitor/post constant for comparison.
- Reuse one monitor/post pair across many solver variants for fair diagnostics.

## Canonical Contract References

- https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
