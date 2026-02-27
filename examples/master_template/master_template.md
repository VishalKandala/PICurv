# PIC-Flow Master Configuration Templates

## 1. Overview

This directory contains fully commented reference templates for:

- `master_case.yml`
- `master_solver.yml`
- `master_monitor.yml`
- `master_postprocessor.yml`
- `master_cluster.yml`
- `master_study.yml`

These are reference contracts, not intended to be run directly.

## 2. Recommended Usage

1. Start with a runnable example:
   - `./scripts/pic.flow init flat_channel --dest my_study`
2. Modify study-local YAML files for your case.
3. Use master templates to discover advanced options and copy validated snippets.

## 3. Best Practices

- Keep reusable solver/monitor/post profiles in `config/solvers`, `config/monitors`, and `config/postprocessors`.
- Keep case-specific physics/geometries in study-local files.
- Prefer structured schema keys over passthrough flags when both exist.

## 4. Contract and Mapping Docs

- User contract: https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- Developer ingestion map: https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- Extension playbook: https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
- Workflow extensibility notes: https://vishalkandala.me/picurv-docs/17_Workflow_Extensibility.html
- Conductor CLI: https://vishalkandala.me/picurv-docs/05_The_Conductor_Script.html
