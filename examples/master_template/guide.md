# Master Template Guide

This directory provides fully commented template files for all major YAML roles.
Use these as schema references, then copy only needed sections into runnable case directories.

## Template Files

- `master_case.yml`
- `master_solver.yml`
- `master_monitor.yml`
- `master_postprocessor.yml`
- `master_cluster.yml`
- `master_study.yml`

## Suggested Workflow

1. start from template relevant to your role,
2. trim to minimum required fields,
3. validate with `pic.flow validate`,
4. run with `pic.flow run` or `pic.flow sweep`.

## Related Docs

- `master_template.md`
- https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
