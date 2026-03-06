@page 41_Getting_Started_Index Getting Started

@anchor _Getting_Started_Index

This index is the recommended onboarding path for new users.
It is organized to minimize first-run failure risk:
build first, then run a known working case, then inspect outputs.
PICurv's YAML files are modular by design, so the goal is to learn how to reuse and mix profiles,
not to create one giant per-run configuration from scratch.

@tableofcontents

@section p41_when_sec 1. When to Use This Section

Use this section if you are in one of these states:

- first time cloning PICurv,
- returning after dependency/toolchain changes,
- onboarding someone to the repository and needing a deterministic starting path.

@section p41_distinction_sec 2. First Simulation Entry

The canonical first-run page is:

- **@subpage 02_Tutorial_Programmatic_Grid**
  complete walkthrough with artifact inspection and troubleshooting context.

@section p41_path_sec 3. Recommended Read Order

1. **@subpage 01_Installation**
2. **@subpage 02_Tutorial_Programmatic_Grid**
3. **@subpage 03_Tutorial_File-Based_Grid**
4. **@subpage 04_Visualization_Tutorial**
5. **@subpage 05_The_Conductor_Script**
6. **@subpage 49_Workflow_Recipes_and_Config_Cookbook**

@section p41_outputs_sec 4. Expected Outcomes After Completing This Path

You should be able to:

- build `bin/simulator` and `bin/postprocessor`,
- generate valid runtime control artifacts from YAML inputs,
- execute `run --solve --post-process` locally,
- inspect VTK outputs in ParaView,
- recombine `case.yml`, `solver.yml`, `monitor.yml`, and `post.yml` intentionally,
- map first-run failures to corrective actions.

@section p41_next_sec 5. Where to Go Next

- Operational run authoring: **@subpage 42_User_Guide_Index**
- Practical reusable config patterns: **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- Direct grid generation and wrapped `grid.mode: grid_gen`: **@subpage 48_Grid_Generator_Guide**
- Troubleshooting catalog: **@subpage 39_Common_Fatal_Errors**
- CI/smoke test contract: **@subpage 40_Testing_and_Quality_Guide**
- Cross-axis map view: **@subpage Documentation_Map**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Getting Started** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.

