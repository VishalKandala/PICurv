# PICurv Master Configuration Templates

## 1. Overview

This directory contains fully commented reference templates for all major workflow roles:

- `master_case.yml`
- `master_solver.yml`
- `master_monitor.yml`
- `master_postprocessor.yml`
- `master_cluster.yml`
- `execution.example.yml`
- `master_study.yml`

These files are intended to document contract breadth and option interactions. They are not meant to be executed directly as-is.
For a concrete end-to-end observability example focused on particle settlement,
search effort, and migration behavior, use `examples/search_robustness/` as the
paired runtime reference.

## 2. Recommended Usage Pattern

1. Initialize a runnable example:
   - `./bin/picurv init flat_channel --dest my_study`
2. Keep edits in study-local YAML files.
3. Use master templates to discover advanced options and copy validated snippets.
4. Re-run `picurv validate` after each copied block.

For post-processing templates, treat `run_control` as the full logical window you want the recipe to represent. During `picurv run --post-process --continue --run-dir ... --post ...`, PICurv now computes the effective restart step internally, caps each launch to the currently available live-source frontier, and preserves single-writer safety on the run directory.

For shared site MPI quirks across login-node and batch runs, `picurv init` now creates `.picurv-execution.yml` in each new case with inert defaults. Existing cases or repo-root site configs can still start from `execution.example.yml`. Keep `master_cluster.yml` focused on batch scheduler policy and batch-only overrides. Legacy `.picurv-local.yml` still works for local-only compatibility.

## 3. How To Read A Master Template Efficiently

- Read the role purpose first (`case`, `solver`, `monitor`, `post`, `cluster`, `study`).
- Identify required fields and defaults.
- Mark optional blocks relevant to your scenario (for example restart, sweep, or cluster paths).
- Copy the smallest meaningful block and test before adding more options.

## 4. Best Practices

- Keep reusable solver/monitor/post profiles in `config/solvers`, `config/monitors`, and `config/postprocessors`.
- Keep geometry-specific and study-specific physics choices in case-local files.
- Prefer explicit schema keys over passthrough flags when both exist.
- Document non-default assumptions in the case directory for reproducibility.

## 5. Contract and Mapping Docs

- User contract: https://vishalkandala.me/picurv-docs/14_Config_Contract.html
- Developer ingestion map: https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- Extension playbook: https://vishalkandala.me/picurv-docs/16_Config_Extension_Playbook.html
- Workflow extensibility notes: https://vishalkandala.me/picurv-docs/17_Workflow_Extensibility.html
- Conductor CLI: https://vishalkandala.me/picurv-docs/05_The_Conductor_Script.html
- Search robustness example: ../search_robustness/search_robustness.md
