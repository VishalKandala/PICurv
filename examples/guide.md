# Examples Guide

This directory contains runnable case templates used by `./bin/picurv init` and reference configs for sweeps/clusters. The examples are designed as operational baselines, not just demos: each one exercises a different grid/physics pathway and can be used as a seed for production studies.

## Currency Guarantee

The documentation certification treats this directory as a starter-product surface, not incidental test data. The tracked starter-content contract inventories every top-level example directory and every example YAML. On each certification it verifies that every YAML belongs to a declared, valid composition, validates all declared case and study bundles, and runs `picurv init` for every example directory to confirm the copied files match their source template. `master_template/` and `periodic_test/` are reference collections: they are copied and their declared component bundles are validated, while runnable studies should start from a concrete case family.

### Adding or Changing an Example

Treat a new example as a supported user entry point. Add its directory to
`tests/tooling/starter_content_contract.json`; classify reference collections
there, list any non-role YAML as auxiliary, and declare every case/solver/
monitor/post or cluster/study composition it supplies. The audit rejects an
unlisted top-level directory or YAML file, so the documentation certification
cannot pass until the example, its usage guidance, and its validation contract
are updated together. Before committing, run:

```bash
python3 tests/tooling/audit_starter_content.py
```

## Example Families

- `flat_channel/`: baseline first-run case using programmatic grid generation.
- `bent_channel/`: file-based curvilinear grid ingestion and curved-geometry BC behavior.
- `brownian_motion/`: analytical zero-flow particle diffusion verification workflow.
- `search_robustness/`: runtime walking-search and migration robustness characterization, with Brownian baselines plus deterministic `UNIFORM_FLOW` migration stress cases.
- `drift_uniform_flow/`: analytical uniform-flow particle advection drift verification.
- `drift_diffusivity_gradient/`: analytical zero-flow plus verification-source diffusivity-gradient drift verification.
- `interpolation_test/`: TGV3D analytical flow interpolation accuracy test with 10k particles.
- `scatter_verification/`: verification-path prescribed scalar truth injection plus runtime particle-to-grid scatter metrics and grid-size sweeps.
- `decaying_isotropic_turbulence/`: triply periodic LES decay benchmark, and the worked example for online field statistics; isotropy makes the accumulated Reynolds stresses self-checking.
- `turbulent_channel/`: driven periodic plane channel at `Re_tau ~ 1000`, wall-modelled.
  The wall-model counterpart to `decaying_isotropic_turbulence`: that case exercises the
  subgrid closure where there is no wall, this one exercises what only exists because
  there is. Distinct from `periodic_test/driven_channel/les_retau180`, which is
  wall-resolved constant-Smagorinsky at `Re_tau = 180`; this is wall-modelled dynamic
  Smagorinsky with the coefficient averaged over the two homogeneous axes. Ships with a
  laminar initial condition and will not transition on its own - see its README.
- `periodic_test/`: reference collection of periodic-boundary configurations. Alongside the
  plain geometric smoke cases it holds the driven-periodic validation families
  `driven_channel/` (laminar exact verification, DNS at `Re_tau = 180` and `395`, and a
  constant-Smagorinsky LES repeat) and `driven_duct/` (square duct at `Re_b = 4410`,
  which sustains secondary flow of the second kind). Their grids are generated from
  checked-in `.cfg` files under `config/grids/` rather than shipped as `.picgrid`.
- `master_template/`: exhaustive reference templates for all config roles.

## How To Start A New Study

1. Initialize a starter case:
   - `./bin/picurv init flat_channel --dest my_case`
   - `./bin/picurv init bent_channel --dest my_case`
2. Validate copied configs:
   - `./bin/picurv validate --case ... --solver ... --monitor ... --post ...`
3. If the case uses `grid_gen` or generated inlet profiles, optionally inspect
   deterministic artifacts first:
   - `./bin/picurv precompute --case ... --output-dir precomputed/<name>`
4. Run dry-run planning, then actual solve/post execution.

For long or live runs, keep the full post-analysis window in `post.yml`, then use:

```bash
./bin/picurv run --post-process --continue --run-dir runs/<run_id> --post my_case/post.yml
```

PICurv will resume the same recipe from the first unfinished step, cap the launch to the current live solver frontier, and refuse a second concurrent post writer on the same run directory.

`init` creates the case directory with config files — no binaries are copied. Runtime executables (`simulator`, `postprocessor`) are resolved from the project `bin/` directory via PATH. The initializer also writes `.picurv-origin.json` and an inert `.picurv-execution.yml`, which enables maintenance commands (`status-source`, `build`, `pull-source`, `sync-binaries`, `sync-config`) and gives each case a safe place for site-specific launcher overrides when needed. To pin specific binary versions into a case, run `picurv sync-binaries`.

## Composition Guidance

- Treat config roles as modular. You can often reuse a validated `solver.yml`, `monitor.yml`, or `post.yml` across multiple `case.yml` variants.
- Keep benchmark studies close to example defaults first, then perturb one dimension at a time.
- Use `prescribed_flow.source.type: file` for inspected PICSLICE inlet profiles,
  `source.type: generated` for analytical profiles generated by `picurv`, and
  `source.type: field_slice` when reusing an old `ufield` slice as a new inlet.
- Promote stable custom profiles back into `<repo>/config/` for team-level reuse.
- The shipped implicit profiles enable residual-aware Picard-Jameson
  convergence and adaptive global pseudo-CFL rollback. Preserve those defaults
  for initial runs, then tune from the solver reference rather than changing
  several controller values at once.

## Related Docs

- https://vishalkandala.me/docs/picurv/02_Tutorial_Programmatic_Grid.html
- https://vishalkandala.me/docs/picurv/03_Tutorial_File-Based_Grid.html
- https://vishalkandala.me/docs/picurv/49_Workflow_Recipes_and_Config_Cookbook.html
- https://vishalkandala.me/docs/picurv/37_Sweep_Studies_Guide.html
- https://vishalkandala.me/docs/picurv/24_Dual_Time_Picard_Jameson_RK.html
