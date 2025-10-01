@page 05_The_Conductor_Script The Conductor Script: `pic-flow`

The `pic-flow` script is the primary user interface for the PICurv platform. It is a powerful Python-based "conductor" that automates the entire simulation workflow, from building the source code to launching simulations and post-processing results.

This guide serves as the main reference for all of its commands and command-line arguments. All commands are run from the root directory of your PICurv project.

@tableofcontents

@section usage_sec 1. General Usage

The script is invoked via `./bin/pic-flow` and uses a command-based structure, similar to `git`.

```bash
./bin/pic-flow [COMMAND] [ARGUMENTS...]
```

There are three main commands:
- **`init`**: Creates a new study directory from a template.
- **`build`**: Compiles the C solver and other tools.
- **`run`**: Executes a simulation or post-processing workflow.

You can get help at any time by running:
```bash
./bin/pic-flow --help
./bin/pic-flow run --help
```

@section init_command_sec 2. The `init` Command

The `init` command is used to bootstrap a new simulation study by copying a pre-configured template.

**Usage:**
```bash
./bin/pic-flow init <template_name> [OPTIONS]
```

**Arguments:**
- `<template_name>`: (Required) The name of a template directory located in `examples/`. Common choices are `flat_channel` or `bent_channel`.

**Options:**
- `--dest <new_directory_name>`: Specifies a name for the new study directory. If omitted, the directory will be named after the template.
- `--copy-binaries`: If this flag is present, the C executables (`picsolver`, `postprocessor`) will be physically copied into the new study directory. This creates a fully portable, self-contained study. By default, symbolic links are created instead.

**Example:**
```bash
# Create a new study named "my_channel_case" from the "flat_channel" template.
./bin/pic-flow init flat_channel --dest my_channel_case
```

@section build_command_sec 3. The `build` Command

The `build` command is a wrapper around the project's `Makefile`. It is used to compile all C source code.

**Usage:**
```bash
./bin/pic-flow build [make_arguments...]
```

Any arguments provided after `build` are passed directly to the `make` command.

**Examples:**
```bash
# Compile all executables (equivalent to 'make all')
./bin/pic-flow build

# Build only the postprocessor
./bin/pic-flow build postprocessor

# Clean the entire project of all build artifacts
./bin/pic-flow build clean-project

# Build using a specific system configuration for a cluster
./bin/pic-flow build SYSTEM=cluster
```

@section run_command_sec 4. The `run` Command

The `run` command is the most powerful command. It orchestrates the execution of the C-solver and post-processor based on the provided YAML configuration files.

**Usage:**
```bash
./bin/pic-flow run [STAGES] [INPUTS] [OPTIONS]
```

**Workflow Stages:**
You must specify at least one of the following stages:
- `--solve`: Executes the `picsolver` executable. This will create a new, timestamped output directory inside the `runs/` folder.
- `--post-process`: Executes the `postprocessor` executable. This can be run on the results of a `--solve` stage in the same command, or on a previously completed run directory using the `--run-dir` flag.

**Input Files:**
- `--case <path/to/case.yml>`: (Required for `--solve`) Path to the case definition file, which specifies the physics and geometry.
- `--solver <path/to/solver.yml>`: (Required for `--solve`) Path to the solver profile, which specifies the numerical strategy.
- `--monitor <path/to/monitor.yml>`: (Required for `--solve`) Path to the monitor profile, which controls I/O and logging.
- `--post <path/to/post.yml>`: (Required for `--post-process`) Path to the post-processing recipe file.
- `--run-dir <path/to/run_directory>`: (Required for standalone `--post-process`) Specifies an existing run directory (e.g., `runs/flat_channel_...`) to post-process. Not needed if running `--solve` and `--post-process` together.

**Execution Options:**
- `-n, --num-procs <integer>`: The number of MPI processes to use for the simulation. Defaults to `1` (serial execution).

**Examples:**

1.  **Run a solver simulation on 8 cores:**
    ```bash
    ./bin/pic-flow run --solve -n 8 \
        --case my_studies/case.yml \
        --solver my_studies/solver.yml \
        --monitor my_studies/monitor.yml
    ```

2.  **Run a solver simulation and immediately post-process the results:**
    ```bash
    ./bin/pic-flow run --solve --post-process -n 8 \
        --case my_studies/case.yml \
        --solver my_studies/solver.yml \
        --monitor my_studies/monitor.yml \
        --post my_studies/post-recipe.yml
    ```

3.  **Run post-processing on an existing simulation directory:**
    ```bash
    ./bin/pic-flow run --post-process \
        --run-dir runs/flat_channel_20240401-153000 \
        --post my_studies/another-post-recipe.yml
    ```

@section next_steps_sec 5. Next Steps

Now that you understand how to use the `pic-flow` conductor, the next step is to learn the details of the configuration files it consumes.

Proceed to the **@subpage 06_Case_Reference** to learn about all the parameters available in the main `case.yml` file.
