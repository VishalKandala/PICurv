@page 41_Getting_Started_Index Getting Started

This index is the recommended onboarding path for new users.
It is organized to minimize first-run failure risk:
build first, then run a known working case, then inspect outputs.
PICurv's YAML files are modular by design, so the goal is to learn how to reuse and mix profiles,
not to create one giant per-run configuration from scratch.

@tableofcontents

@section when_sec 1. When to Use This Section

Use this section if you are in one of these states:

- first time cloning PICurv,
- returning after dependency/toolchain changes,
- onboarding someone to the repository and needing a deterministic starting path.

@section distinction_sec 2. First Simulation Entry

The canonical first-run page is:

- **@subpage 02_Tutorial_Programmatic_Grid**
  complete walkthrough with artifact inspection and troubleshooting context.

@section path_sec 3. Recommended Read Order

1. **@subpage 01_Installation**
2. **@subpage 02_Tutorial_Programmatic_Grid**
3. **@subpage 03_Tutorial_File-Based_Grid**
4. **@subpage 04_Visualization_Tutorial**
5. **@subpage 05_The_Conductor_Script**
6. **@subpage 49_Workflow_Recipes_and_Config_Cookbook**

@section outputs_sec 4. Expected Outcomes After Completing This Path

You should be able to:

- build `bin/simulator` and `bin/postprocessor`,
- generate valid runtime control artifacts from YAML inputs,
- execute `run --solve --post-process` locally,
- inspect VTK outputs in ParaView,
- recombine `case.yml`, `solver.yml`, `monitor.yml`, and `post.yml` intentionally,
- map first-run failures to corrective actions.

@section next_sec 5. Where to Go Next

- Operational run authoring: **@subpage 42_User_Guide_Index**
- Practical reusable config patterns: **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- Direct grid generation and wrapped `grid.mode: grid_gen`: **@subpage 48_Grid_Generator_Guide**
- Troubleshooting catalog: **@subpage 39_Common_Fatal_Errors**
- CI/smoke test contract: **@subpage 40_Testing_and_Quality_Guide**
- Cross-axis map view: **@subpage Documentation_Map**
