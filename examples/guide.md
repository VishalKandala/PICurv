# Examples Guide

This directory contains runnable case templates used by `./bin/picurv init` and reference configs for sweeps/clusters. The examples are designed as operational baselines, not just demos: each one exercises a different grid/physics pathway and can be used as a seed for production studies.

## Example Families

- `flat_channel/`: baseline first-run case using programmatic grid generation.
- `bent_channel/`: file-based curvilinear grid ingestion and curved-geometry BC behavior.
- `brownian_motion/`: analytical zero-flow particle diffusion verification workflow.
- `interpolation_test/`: TGV3D analytical flow interpolation accuracy test with 10k particles.
- `master_template/`: exhaustive reference templates for all config roles.

## How To Start A New Study

1. Initialize a starter case:
   - `./bin/picurv init flat_channel --dest my_case`
   - `./bin/picurv init bent_channel --dest my_case`
2. Validate copied configs:
   - `./bin/picurv validate --case ... --solver ... --monitor ... --post ...`
3. Run dry-run planning, then actual solve/post execution.

If project binaries are already built, `init` copies available executables into the case directory so it can run self-contained. The initializer also writes `.picurv-origin.json` and an inert `.picurv-execution.yml`, which enables case-local maintenance commands (`status-source`, `build`, `pull-source`, `sync-binaries`, `sync-config`) and gives each case a safe place for site-specific launcher overrides when needed.

## Composition Guidance

- Treat config roles as modular. You can often reuse a validated `solver.yml`, `monitor.yml`, or `post.yml` across multiple `case.yml` variants.
- Keep benchmark studies close to example defaults first, then perturb one dimension at a time.
- Promote stable custom profiles back into `config/` for team-level reuse.

## Related Docs

- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
