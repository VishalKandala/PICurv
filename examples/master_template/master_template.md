# PIC-Flow Master Configuration Templates

## 1. Overview

This directory contains master templates for each of the core PIC-Flow configuration files: `master_case.yml`, `master_solver.yml`,`master_postprocessor.yml` and `master_monitor.yml`.

These files are **not intended to be run directly**. They are heavily commented reference documents that showcase **every possible configuration option** that the `pic.flow` platform understands.

## 2. How to Use These Files

1.  **Reference Manual:** When you are unsure about a specific setting or want to see all available options for a section, open the relevant master template file. For example, to see all available multigrid settings, open `master_solver.yml` and look under the `pressure_solver.multigrid` section.

2.  **Snippet Library:** Use these files to copy-paste sections into your own configuration files. If you start with a simple `flat_channel` case but want to add turbulence modeling, you can copy the `models.physics.turbulence` section from `master_case.yml` into your study's `case.yml`.

3.  **Learning Tool:** Reading through the comments in these files is the best way to learn about the full capabilities of the PIC-Flow platform.

## 3. Best Practices & Workflow

The recommended workflow is to combine simple, working templates with the detailed master templates:

1.  **Start with `init`:** Always begin a new project by running `pic.flow init <template_name>` with a simple, working template like `flat_channel` or `bent_channel`. This gives you a guaranteed-to-run starting point.
2.  **Customize:** Edit the `case.yml` in your new study directory to match the specific physics you want to simulate.
3.  **Enhance with Master Templates:** When you need to add more advanced features (e.g., a custom PETSc solver setting, particle physics, or a different turbulence model), refer to the master templates, copy the relevant sections, and paste them into your study's configuration files.
4.  **Build Your Own Library:** As you create new, useful configurations, save them in the central `solver_profiles/` and `monitor_profiles/` directories so you can easily reuse them in future studies.