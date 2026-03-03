# ==============================================================================
#
#                    Makefile for the PICurv Project
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
#   relying on an external configuration system (`config/build/config.*.mk` files) to adapt
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
#      - `PETSC_ARCH`: Optional(required only for old style in-tree builds.) New style prefix installs don't use it.
#
# ------------------------------------------------------------------------------
#
# ## File Structure
#
#   This build system expects the following file structure:
#
#   .
#   ├── Makefile            (This file)
#   ├── config/build/config.local.mk     (Configuration for local development)
#   ├── config/build/config.cluster.mk   (Configuration for the HPC cluster)
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
SCRIPTDIR := scripts
TESTCDIR := tests/c
TESTOBJDIR := $(OBJDIR)/tests
TESTBINDIR := $(BINDIR)/tests
SMOKEDIR := tests/smoke
SMOKE_RUNNER := $(SMOKEDIR)/run_smoke.sh
TEST_NPROCS ?= 1

# --- 2. System Configuration ---
# Select and include the appropriate configuration file based on the SYSTEM variable.
SYSTEM ?= local
NO_CONFIG_GOALS := test test-python doctor install-check audit-ingress build-docs open-docs tags clean-project cleanobj clean-project-docs clean-project-tags clean-unit
NEEDS_BUILD_CONFIG := 1

ifneq ($(MAKECMDGOALS),)
ifeq ($(strip $(filter-out $(NO_CONFIG_GOALS),$(MAKECMDGOALS))),)
NEEDS_BUILD_CONFIG :=
endif
endif

ifdef NEEDS_BUILD_CONFIG
CONFIG_FILE := config/build/config.$(SYSTEM).mk
ifneq ("$(wildcard $(CONFIG_FILE))","")
    include $(CONFIG_FILE)
else
    $(error Cannot find configuration file '$(CONFIG_FILE)'. Please create it or specify a valid SYSTEM.)
endif
$(info Building for system: $(SYSTEM_NAME))
endif


# --- 3. Source & Object File Definitions ---
# Explicitly list the object files required for each final executable.
SIMULATOR_OBJS := $(addprefix $(OBJDIR)/, \
                 simulator.o setup.o logging.o grid.o io.o Metric.o AnalyticalSolutions.o\
                 Boundaries.o BC_Handlers.o wallfunction.o runloop.o walkingsearch.o BodyForces.o\
                 ParticleSwarm.o ParticleMotion.o ParticlePhysics.o interpolation.o \
                 initialcondition.o rhs.o solvers.o momentumsolvers.o poisson.o\
				 les.o  Filter.o)

POSTPROCESSOR_OBJS := $(addprefix $(OBJDIR)/, \
                     postprocessor.o setup.o logging.o grid.o io.o Metric.o AnalyticalSolutions.o\
                     Boundaries.o BC_Handlers.o wallfunction.o postprocessing_kernels.o vtk_io.o \
					 ParticleSwarm.o ParticleMotion.o interpolation.o walkingsearch.o \
					 particle_statistics.o)

# --- 4. Executable Definitions ---
# Define the final paths for the compiled programs.
SIMULATOR_EXE     := $(BINDIR)/simulator
POSTPROCESSOR_EXE := $(BINDIR)/postprocessor
CONDUCTOR_EXE     := $(BINDIR)/picurv
DOCTOR_EXE        := $(TESTBINDIR)/doctor
UNIT_GEOMETRY_EXE := $(TESTBINDIR)/unit_geometry
UNIT_SOLVER_EXE   := $(TESTBINDIR)/unit_solver
UNIT_PARTICLES_EXE := $(TESTBINDIR)/unit_particles
UNIT_IO_EXE       := $(TESTBINDIR)/unit_io
UNIT_POST_EXE     := $(TESTBINDIR)/unit_post
TEST_CFLAGS_TO_USE := $(CFLAGS_TO_USE) -I$(TESTCDIR)
TEST_SUPPORT_OBJ  := $(TESTOBJDIR)/test_support.o
DOCTOR_OBJ        := $(TESTOBJDIR)/test_install_check.o
UNIT_GEOMETRY_OBJ := $(TESTOBJDIR)/test_geometry.o
UNIT_SOLVER_OBJ   := $(TESTOBJDIR)/test_solver_kernels.o
UNIT_PARTICLES_OBJ := $(TESTOBJDIR)/test_particle_kernels.o
UNIT_IO_OBJ       := $(TESTOBJDIR)/test_io.o
UNIT_POST_OBJ     := $(TESTOBJDIR)/test_postprocessing.o
TEST_COMMON_OBJS  := $(sort $(filter-out $(OBJDIR)/simulator.o $(OBJDIR)/postprocessor.o,$(SIMULATOR_OBJS) $(POSTPROCESSOR_OBJS)))

# ==============================================================================
# --- 5. Build Targets & Rules ---
# ==============================================================================
.PHONY: all simulator postprocessor conductor dirs

## @target all
## @brief Default target. Builds all project executables.
all: simulator postprocessor conductor

## @target simulator
## @brief Builds the `simulator` executable.
simulator: $(SIMULATOR_EXE)

## @target postprocessor
## @brief Builds the `postprocessor` executable.
postprocessor: $(POSTPROCESSOR_EXE)

## @target conductor
## @brief Installs the `picurv` conductor script.
conductor: $(CONDUCTOR_EXE)

# Explicit rule for linking the main solver.
$(SIMULATOR_EXE): $(SIMULATOR_OBJS) | dirs
	@echo "--- Linking Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)
	@echo "--- Build Complete: $(@) ---"

# Explicit rule for linking the post-processor.
$(POSTPROCESSOR_EXE): $(POSTPROCESSOR_OBJS) | dirs
	@echo "--- Linking Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)
	@echo "--- Build Complete: $(@) ---"

# This rule copies the source script to the bin directory and makes it executable.
$(CONDUCTOR_EXE): $(SCRIPTDIR)/picurv | dirs
	@echo "--- Installing Conductor Script: $(@) ---"
	@cp $< $@
	@chmod +x $@

# Generic rule for compiling any .c file from SRCDIR into an object file in OBJDIR.
$(OBJDIR)/%.o: $(SRCDIR)/%.c | dirs
	@echo "--- Compiling: $< ---"
	$(CC_TO_USE) $(CFLAGS_TO_USE) -c $< -o $@

# Generic rule for compiling any test .c file from TESTCDIR into TESTOBJDIR.
$(TESTOBJDIR)/%.o: $(TESTCDIR)/%.c | dirs
	@echo "--- Compiling Test: $< ---"
	$(CC_TO_USE) $(TEST_CFLAGS_TO_USE) -c $< -o $@

$(DOCTOR_EXE): $(DOCTOR_OBJ) $(TEST_SUPPORT_OBJ) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_GEOMETRY_EXE): $(UNIT_GEOMETRY_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_SOLVER_EXE): $(UNIT_SOLVER_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_PARTICLES_EXE): $(UNIT_PARTICLES_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_IO_EXE): $(UNIT_IO_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_EXE): $(UNIT_POST_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

## @target dirs
## @brief (Internal) Ensures all necessary build directories exist.
dirs: 
	@mkdir -p $(OBJDIR) $(BINDIR) $(DOCSDIR) $(TESTOBJDIR) $(TESTBINDIR)

# ==============================================================================
# --- 6. Execution, Auxiliary, & Cleanup Targets ---
# ==============================================================================
.PHONY: run test test-python doctor doctor-runner install-check smoke unit unit-geometry unit-solver unit-particles unit-io unit-post ctest ctest-geometry ctest-solver ctest-particles ctest-io ctest-post check build-docs open-docs tags audit-ingress clean-project cleanobj clean-project-docs clean-project-tags clean-unit

## @target run
## @brief Runs the main solver using the system-specific MPI launcher.
run: $(SIMULATOR_EXE)
	$(MPI_LAUNCHER) -n $(NPROCS) $< $(ARGS)

## @target test-python
## @brief Runs the Python regression test suite.
test-python:
	@python3 -m pytest -q

## @target test
## @brief Backward-compatible alias for `test-python`.
test: test-python

## @target doctor
## @brief Validates PETSc provisioning by building and running a minimal PETSc smoke binary.
doctor:
	@if [ -z "$$PETSC_DIR" ]; then \
		echo "PICurv doctor: PETSC_DIR is not set."; \
		echo "Export PETSC_DIR (and PETSC_ARCH if your install requires it), then rerun 'make doctor'."; \
		exit 1; \
	fi
	@if [ ! -d "$$PETSC_DIR" ]; then \
		echo "PICurv doctor: PETSC_DIR points to a missing directory: $$PETSC_DIR"; \
		exit 1; \
	fi
	$(if $(findstring n,$(MAKEFLAGS)),@echo "$(MAKE) --no-print-directory doctor-runner SYSTEM=$(SYSTEM)",@$(MAKE) --no-print-directory doctor-runner SYSTEM=$(SYSTEM))

## @target doctor-runner
## @brief (Internal) Builds and runs the PETSc installation smoke binary.
doctor-runner: $(DOCTOR_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target install-check
## @brief Compatibility alias for `doctor`.
install-check: doctor

## @target unit-geometry
## @brief Runs the geometry and interpolation C unit tests.
unit-geometry: $(UNIT_GEOMETRY_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-solver
## @brief Runs the solver utility C unit tests.
unit-solver: $(UNIT_SOLVER_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-particles
## @brief Runs the particle kernel C unit tests.
unit-particles: $(UNIT_PARTICLES_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-io
## @brief Runs the I/O infrastructure C unit tests.
unit-io: $(UNIT_IO_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-post
## @brief Runs the post-processing/statistics C unit tests.
unit-post: $(UNIT_POST_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit
## @brief Runs the full isolated C unit/component suite.
unit: unit-geometry unit-solver unit-particles unit-io unit-post

## @target ctest
## @brief Compatibility alias for `unit`.
ctest: unit

## @target ctest-geometry
## @brief Compatibility alias for `unit-geometry`.
ctest-geometry: unit-geometry

## @target ctest-solver
## @brief Compatibility alias for `unit-solver`.
ctest-solver: unit-solver

## @target ctest-particles
## @brief Compatibility alias for `unit-particles`.
ctest-particles: unit-particles

## @target ctest-io
## @brief Compatibility alias for `unit-io`.
ctest-io: unit-io

## @target ctest-post
## @brief Compatibility alias for `unit-post`.
ctest-post: unit-post

## @target smoke
## @brief Runs lightweight executable-level simulator/postprocessor smoke checks.
smoke: simulator postprocessor
	@bash $(SMOKE_RUNNER) "$(SIMULATOR_EXE)" "$(POSTPROCESSOR_EXE)" "$(MPI_LAUNCHER)" "$(TEST_NPROCS)"

## @target check
## @brief Runs the full local validation sweep (Python, doctor, unit, smoke).
check:
	@$(MAKE) --no-print-directory test-python
	@$(MAKE) --no-print-directory doctor-runner SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory unit SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory smoke SYSTEM=$(SYSTEM)

## @target build-docs
## @brief Generates Doxygen documentation for the project.
build-docs: dirs
	@echo "==> Generating Doxygen documentation..."
	@mkdir -p logs
	@doxygen docs/Doxyfile
	@echo "HTML documentation generated in docs_build/html/index.html"
	@echo "Doxygen warnings log: logs/doxygen.warnings"

## @target open-docs
## @brief Opens the generated documentation in a web browser.
open-docs: build-docs
ifeq ($(shell uname),Darwin)
	@open docs_build/html/index.html
else
	@xdg-open docs_build/html/index.html || echo "Could not open browser automatically."
endif

## @target tags
## @brief Generates an Emacs TAGS file for code navigation.
tags:
	@echo "==> Generating TAGS file..."
	@find $(SRCDIR) $(INCDIR) -type f \( -name "*.c" -o -name "*.h" \) -print | etags -

## @target audit-ingress
## @brief Audits PETSc option ingress in setup/io against the maintained manifest.
audit-ingress:
	@python3 scripts/audit_ingress.py

## @target clean-project
## @brief Removes all build artifacts generated by this project.
clean-project: cleanobj clean-unit clean-project-docs clean-project-tags
	@echo "Project cleaned."

## @target cleanobj
## @brief (Internal) Removes object files and executables.
cleanobj:
	@echo "--- Removing object files and compiled executables"
	@rm -rf $(OBJDIR) 
	@rm -f $(SIMULATOR_EXE) $(POSTPROCESSOR_EXE) $(CONDUCTOR_EXE)

## @target clean-unit
## @brief Removes compiled C test objects and binaries.
clean-unit:
	@rm -rf $(TESTOBJDIR)
	@rm -rf $(TESTBINDIR)

## @target clean-project-docs
## @brief (Internal) Removes generated documentation.
clean-project-docs:
	@rm -rf docs_build

## @target clean-project-tags
## @brief (Internal) Removes the TAGS file.
clean-project-tags:
	@rm -f TAGS
.PHONY: show-config
show-config:
	@echo "SYSTEM=$(SYSTEM)"
	@echo "SYSTEM_NAME=$(SYSTEM_NAME)"
	@echo "PETSC_DIR=$(PETSC_DIR)"
	@echo "PETSC_ARCH=$(PETSC_ARCH)"
	@echo "CC_TO_USE=$(CC_TO_USE)"
	@echo "CFLAGS_TO_USE=$(CFLAGS_TO_USE)"
	@echo "LINKER_TO_USE=$(LINKER_TO_USE)"
	@echo "LIBS_TO_USE=$(LIBS_TO_USE)"
	@echo "MPI_LAUNCHER=$(MPI_LAUNCHER)"
