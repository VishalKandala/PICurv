# ==============================================================================
#
#                    Makefile for the PIC Solver Project
#
# ==============================================================================
#
#  Author:      Vishal Indivar Kandala
#  Version:     1.0 (Platform-Agnostic Release)
#  Date:        2025-08-27
#  License:     MIT License
#  Repository:  https://github.com/VishalKandala/PICurv 
#
# ------------------------------------------------------------------------------
#
# ## Description
#
#   This Makefile manages the compilation, linking, and execution of a C-based
#   parallel scientific computing project. It is engineered for portability,
#   relying on an external configuration system (`config.*.mk` files) to adapt
#   to different build environments (e.g., local development vs. HPC cluster).
#
#   Its core philosophy is to separate the build *logic* (in this file) from
#   the environment *configuration* (in the config files).
#
# ------------------------------------------------------------------------------
#
# ## Prerequisites
#
#   Before using this Makefile, ensure the following are installed and configured:
#
#   1. **Build Tools:** A standard C toolchain (e.g., GCC, Clang) and GNU Make.
#   2. **MPI Library:** An MPI implementation (e.g., OpenMPI, MPICH) that provides `mpicc`.
#   3. **PETSc Library:** A compiled PETSc distribution. The following environment
#      variables MUST be set correctly before invoking make:
#      - `PETSC_DIR`:  Absolute path to the PETSc root directory.
#      - `PETSC_ARCH`: PETSc build architecture name (e.g., `arch-linux-c-debug`).
#
# ------------------------------------------------------------------------------
#
# ## File Structure
#
#   This build system expects the following file structure:
#
#   .
#   ├── Makefile            (This file)
#   ├── config.local.mk     (Configuration for local development)
#   ├── config.cluster.mk   (Configuration for the HPC cluster)
#   ├── src/                (Directory for all .c source files)
#   ├── include/            (Directory for all .h header files)
#   └── ...
#
# ------------------------------------------------------------------------------
#
# ## Usage & Configuration
#
#   All commands are run from the project root directory.
#
#   ### Primary Targets
#
#   - **Build all executables:**
#     $ make
#     $ make SYSTEM=cluster
#
#   - **Run the main solver:**
#     $ make run
#
#   - **Clean the project:**
#     $ make clean-project
#
#   - **Build documentation:**
#     $ make build-docs
#
#   ### Build Customization (via command-line variables)
#
#   - **Specify the build system:** `SYSTEM=cluster` (default is `local`)
#   - **Change processor count for `make run`:** `NPROCS=8`
#   - **Pass arguments to the solver:** `ARGS="--input case.txt"`
#   - **Enable sanitizers (in debug mode):** `SAN=1`
#
#   **Example:** Run on the cluster with 64 processors and a specific input file:
#   $ make run SYSTEM=cluster NPROCS=64 ARGS="--input production.txt"
#
# ==============================================================================

# --- 1. Project Structure ---
# Define the directory layout for the project.
SRCDIR := src
INCDIR := include
OBJDIR := obj
BINDIR := bin
DOCSDIR := docs

# --- 2. System Configuration ---
# Select and include the appropriate configuration file based on the SYSTEM variable.
SYSTEM ?= local
CONFIG_FILE = config.$(SYSTEM).mk
ifneq ("$(wildcard $(CONFIG_FILE))","")
    include $(CONFIG_FILE)
else
    $(error Cannot find configuration file '$(CONFIG_FILE)'. Please create it or specify a valid SYSTEM.)
endif
$(info Building for system: $(SYSTEM_NAME))


# --- 3. Source & Object File Definitions ---
# Explicitly list the object files required for each final executable.
PICSOLVER_OBJS := $(addprefix $(OBJDIR)/, \
                 picsolver.o setup.o logging.o grid.o io.o Metric.o \
                 Boundaries.o wallfunction.o simulation.o walkingsearch.o \
                 ParticleSwarm.o ParticleMotion.o interpolation.o \
                 initialcondition.o rhs.o solvers.o implicitsolvers.o poisson.o)

POSTPROCESSOR_OBJS := $(addprefix $(OBJDIR)/, \
                     postprocessor.o setup.o logging.o grid.o io.o Metric.o \
                     Boundaries.o wallfunction.o postprocessing_kernels.o)

# --- 4. Executable Definitions ---
# Define the final paths for the compiled programs.
PICSOLVER_EXE     := $(BINDIR)/picsolver
POSTPROCESSOR_EXE := $(BINDIR)/postprocessor

# ==============================================================================
# --- 5. Build Targets & Rules ---
# ==============================================================================
.PHONY: all picsolver postprocessor dirs

## @target all
## @brief Default target. Builds all project executables.
all: picsolver postprocessor

## @target picsolver
## @brief Builds the `picsolver` executable.
picsolver: $(PICSOLVER_EXE)

## @target postprocessor
## @brief Builds the `postprocessor` executable.
postprocessor: $(POSTPROCESSOR_EXE)

# Explicit rule for linking the main solver.
$(PICSOLVER_EXE): $(PICSOLVER_OBJS) | dirs
	@echo "--- Linking Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)
	@echo "--- Build Complete: $(@) ---"

# Explicit rule for linking the post-processor.
$(POSTPROCESSOR_EXE): $(POSTPROCESSOR_OBJS) | dirs
	@echo "--- Linking Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)
	@echo "--- Build Complete: $(@) ---"


# Generic rule for compiling any .c file from SRCDIR into an object file in OBJDIR.
$(OBJDIR)/%.o: $(SRCDIR)/%.c | dirs
	@echo "--- Compiling: $< ---"
	$(CC_TO_USE) $(CFLAGS_TO_USE) -c $< -o $@

## @target dirs
## @brief (Internal) Ensures all necessary build directories exist.
dirs: 
	@mkdir -p $(OBJDIR) $(BINDIR) $(DOCSDIR)

# ==============================================================================
# --- 6. Execution, Auxiliary, & Cleanup Targets ---
# ==============================================================================
.PHONY: run build-docs open-docs tags clean-project cleanobj clean-project-docs clean-project-tags

## @target run
## @brief Runs the main solver using the system-specific MPI launcher.
run: $(PICSOLVER_EXE)
	$(MPI_LAUNCHER) -n $(NPROCS) $< $(ARGS)

## @target build-docs
## @brief Generates Doxygen documentation for the project.
build-docs: dirs
	@echo "==> Generating Doxygen documentation..."
	@doxygen docs/Doxyfile
	@echo "HTML documentation generated in docs/docs_build/html/index.html"

## @target open-docs
## @brief Opens the generated documentation in a web browser.
open-docs: build-docs
ifeq ($(shell uname),Darwin)
	@open docs/docs_build/html/index.html
else
	@xdg-open docs/docs_build/html/index.html || echo "Could not open browser automatically."
endif

## @target tags
## @brief Generates an Emacs TAGS file for code navigation.
tags:
	@echo "==> Generating TAGS file..."
	@find $(SRCDIR) $(INCDIR) -type f \( -name "*.c" -o -name "*.h" \) -print | etags -

## @target clean-project
## @brief Removes all build artifacts generated by this project.
clean-project: cleanobj clean-project-docs clean-project-tags
	@echo "Project cleaned."

## @target cleanobj
## @brief (Internal) Removes object files and executables.
cleanobj:
	@rm -rf $(OBJDIR) $(BINDIR)

## @target clean-project-docs
## @brief (Internal) Removes generated documentation.
clean-project-docs:
	@rm -rf docs/docs_build

## @target clean-project-tags
## @brief (Internal) Removes the TAGS file.
clean-project-tags:
	@rm -f TAGS