# Master Template Guide

This directory provides fully commented template files for all major YAML roles. These files are designed as contract references and onboarding aids, not direct production inputs.

A practical pattern is: initialize a runnable example first, then use these templates to discover advanced options and copy only the relevant blocks into your case-local configs.

## Template Files

- `master_case.yml`
- `master_solver.yml`
- `master_monitor.yml`
- `master_postprocessor.yml`
- `master_cluster.yml`
- `master_study.yml`

## Suggested Workflow

1. Start from a template matching the role you want to extend.
2. Copy minimal required blocks into runnable case-local YAML.
3. Keep each role modular rather than merging all settings into one large file.
4. Validate with `picurv validate` after each change.
5. Execute with `picurv run` or `picurv sweep`.

## Why This Matters For CFD Teams

- It reduces configuration drift between projects.
- It makes solver/monitor/post strategy comparisons explicit and reproducible.
- It lowers onboarding cost for users new to the PICurv contract model.

## Related Docs

- `master_template.md`
- https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
