# Master Template Guide

This directory provides fully commented template files for all major YAML roles. These files are designed as contract references and onboarding aids, not direct production inputs.

A practical pattern is: initialize a runnable example first, then use these templates to discover advanced options and copy only the relevant blocks into your case-local configs.
For search and migration observability patterns, the dedicated reference bundle
is `examples/search_robustness/`.

## Template Files

- `master_case.yml`
- `master_solver.yml`
- `master_monitor.yml`
- `master_postprocessor.yml`
- `master_cluster.yml`
- `execution.example.yml`
- `master_study.yml`

## Suggested Workflow

1. Start from a template matching the role you want to extend.
2. Copy minimal required blocks into runnable case-local YAML.
3. Keep each role modular rather than merging all settings into one large file.
4. Validate with `picurv validate` after each change.
5. Execute with `picurv run` or `picurv sweep`.

Post-profile usage note:
- keep `master_postprocessor.yml -> run_control` as the full logical analysis window you want a recipe to represent. When you later run `picurv run --post-process --continue --run-dir ... --post ...`, PICurv resumes the same recipe from the first unfinished step instead of requiring manual `start_step` edits.
- on live solver runs, PICurv also caps each post launch to the highest fully available contiguous source prefix for that recipe.

Launcher-related note:

- `master_cluster.yml` is only for batch/scheduler policy and batch-specific overrides.
- generated Slurm solver jobs use the runtime walltime guard by default; override it in `master_cluster.yml -> execution.walltime_guard` only when the default policy is not appropriate.
- `picurv init` writes `.picurv-execution.yml` into new cases with inert defaults.
- for an existing case or repo-root site config, start from `execution.example.yml`.
- legacy `.picurv-local.yml` remains supported for local-only compatibility, but it does not feed generated batch jobs.
- run directory names are still generated automatically by `picurv` as `<case_basename>_<timestamp>`.

## Why This Matters For CFD Teams

- It reduces configuration drift between projects.
- It makes solver/monitor/post strategy comparisons explicit and reproducible.
- It lowers onboarding cost for users new to the PICurv contract model.

## Related Docs

- `master_template.md`
- `../search_robustness/search_robustness.md`
- https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
