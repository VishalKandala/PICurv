# Configuration Guide

This directory contains reusable YAML profiles and build/runtime configuration assets. It is best treated as a configuration library: stable, versioned building blocks that users compose into case-specific workflows.

For CFD users, the key idea is separation of concerns. Instead of creating one monolithic YAML file, PICurv uses role-oriented contracts (`case`, `solver`, `monitor`, `post`, and optional `cluster`/`study`) so you can change numerical strategy without rewriting geometry definitions, or change post outputs without touching solver controls.

Two repo-wide patterns are especially important in the current codebase:

- `solver.yml -> verification.sources.*` is reserved for verification-only injections/overrides when no cleaner end-to-end path exists. For example, `verification.sources.scalar` prescribes particle `Psi` from analytical truth and enables the runtime diagnostic `logs/scatter_metrics.csv` without changing ordinary production runs.
- `study.yml` supports either cross-product sweeps under `parameters:` or explicit coupled bundles under `parameter_sets:` when multiple overrides must move together.

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
4. For Slurm workflows, prefer `picurv run ... --cluster ... --no-submit` first, inspect `runs/<run_id>/scheduler/`, then use `picurv submit --run-dir ...`.
5. If a submitted Slurm run needs to be stopped, use `picurv cancel --run-dir ...` instead of relying on a separately tracked job ID. For solver jobs that should write a final off-cadence checkpoint first, use `picurv cancel --run-dir ... --stage solve --graceful`.
6. Generated Slurm solver jobs enable the runtime walltime guard by default; tune it in `cluster.yml -> execution.walltime_guard` only when the default warmup/headroom policy needs adjustment.
7. Keep `cluster.yml -> execution.extra_sbatch.signal` as the fallback path for preemption/termination signals or runs that may not reach the guard warmup window.
8. Promote stable reusable profiles back into `config/` so team workflows converge.
9. When designing verification studies, start from the verification pathway in `solver.yml` rather than embedding one-off study logic into production configs.

## Common CFD Iteration Patterns

- Change `solver.yml` while freezing `case.yml` to isolate numerics.
- Change `case.yml` grid/BC setup while keeping monitor/post constant for comparison.
- Reuse one monitor/post pair across many solver variants for fair diagnostics.
- Stage cluster jobs, inspect generated scheduler artifacts, and submit only after the rendered launcher/resources look correct.

## Canonical Contract References

- https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
