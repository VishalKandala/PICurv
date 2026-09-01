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
LOGDIR := logs
DOCSDIR := docs
PICURV_ENTRYPOINT := picurv_cli/picurv
TESTCDIR := tests/c
TESTOBJDIR := $(OBJDIR)/tests
TESTBINDIR := $(BINDIR)/tests
SMOKEDIR := tests/smoke
SMOKE_RUNNER := $(SMOKEDIR)/run_smoke.sh
DRIVEN_PERIODIC_RUNNER := $(SMOKEDIR)/run_driven_periodic_regression.sh
# Rank counts for the driven-periodic regression. 4 and 10 are the counts that
# previously disagreed on tracked-vs-true Poisson residuals; keep both.
DRIVEN_PERIODIC_NPROCS ?= 4 10
TEST_NPROCS ?= 1
TEST_MPI_NPROCS ?= 2
SMOKE_MPI_NPROCS ?= $(TEST_MPI_NPROCS)
SMOKE_MPI_MATRIX_NPROCS ?= 2 3
PY_COVERAGE_MIN ?= 40
C_COVERAGE_MIN ?= 70
BUILD_LOG := $(LOGDIR)/build.log
BUILD_WARNINGS_LOG := $(LOGDIR)/build.warnings.log
BUILD_AUDIT_GOALS ?= cleanobj clean-unit all unit

# --- 2. System Configuration ---
# Select and include the appropriate configuration file based on the SYSTEM variable.
SYSTEM ?= local
PORT ?= 8000
BIND ?= 127.0.0.1

NO_CONFIG_GOALS := test test-python coverage-python doctor install-check audit-ingress audit-docs-expansion audit-docs-site docs-inventory docs-topology docs-xref audit-capability audit-contracts audit-path-literals audit-family-census audit-page-types audit-inline-choices audit-subsystems audit-freshness audit-agent-setup sync-agent-skills docs-cli-reference audit-cli-reference scaffold review-packet build-docs preview-docs open-docs tags certify-docs certify-docs-fast install-git-hooks clean-project cleanobj clean-project-docs clean-project-tags clean-unit conductor
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

ifeq ($(COVERAGE),1)
CFLAGS_TO_USE += --coverage
LIBS_TO_USE += --coverage
endif

ifneq ($(wildcard $(OBJDIR)/*.gcno $(TESTOBJDIR)/*.gcno),)
LIBS_TO_USE += --coverage
endif

# --- 3. Source & Object File Definitions ---
# Explicitly list the object files required for each final executable.
SIMULATOR_OBJS := $(addprefix $(OBJDIR)/, \
                 simulator.o setup.o checksum.o statistics_moments.o statistics_target.o statistics_window.o statistics_accumulator.o statistics_config.o field_catalog.o particle_field_catalog.o logging.o grid.o io.o Metric.o AnalyticalSolutions.o\
                 Boundaries.o BC_Handlers.o wallfunction.o runloop.o walkingsearch.o BodyForces.o\
                 ParticleSwarm.o ParticleMotion.o ParticlePhysics.o interpolation.o \
                 initialcondition.o rhs.o solvers.o momentumsolvers.o momentum_newton_krylov.o poisson.o verification_sources.o\
				 les.o  Filter.o)

POSTPROCESSOR_OBJS := $(addprefix $(OBJDIR)/, \
                     postprocessor.o setup.o checksum.o statistics_moments.o statistics_target.o statistics_window.o statistics_accumulator.o statistics_config.o field_catalog.o particle_field_catalog.o logging.o grid.o io.o Metric.o AnalyticalSolutions.o\
                     Boundaries.o BC_Handlers.o wallfunction.o postprocessing_kernels.o vtk_io.o \
					 ParticleSwarm.o ParticleMotion.o interpolation.o walkingsearch.o \
					 particle_statistics.o verification_sources.o)

# --- 4. Executable Definitions ---
# Define the final paths for the compiled programs.
SIMULATOR_EXE     := $(BINDIR)/simulator
POSTPROCESSOR_EXE := $(BINDIR)/postprocessor
CONDUCTOR_EXE     := $(BINDIR)/picurv
DOCTOR_EXE        := $(TESTBINDIR)/doctor
UNIT_GEOMETRY_EXE := $(TESTBINDIR)/unit_geometry
UNIT_SETUP_EXE    := $(TESTBINDIR)/unit_setup
UNIT_SOLVER_EXE   := $(TESTBINDIR)/unit_solver
UNIT_NEWTON_KRYLOV_EXE := $(TESTBINDIR)/unit_momentum_newton_krylov
UNIT_NEWTON_FIXEDPOINT_EXE := $(TESTBINDIR)/unit_momentum_newton_boundary_fixedpoint
UNIT_PARTICLES_EXE := $(TESTBINDIR)/unit_particles
UNIT_STATISTICS_EXE := $(TESTBINDIR)/unit_statistics
UNIT_STATISTICS_TARGET_EXE := $(TESTBINDIR)/unit_statistics_target
UNIT_STATISTICS_WINDOW_EXE := $(TESTBINDIR)/unit_statistics_window
UNIT_STATISTICS_ACC_EXE := $(TESTBINDIR)/unit_statistics_accumulator
UNIT_STATISTICS_CFG_EXE := $(TESTBINDIR)/unit_statistics_config
UNIT_IO_EXE       := $(TESTBINDIR)/unit_io
UNIT_LOGGING_EXE  := $(TESTBINDIR)/unit_logging
UNIT_POST_EXE     := $(TESTBINDIR)/unit_post
UNIT_POST_VTK_EXE := $(TESTBINDIR)/unit_post_vtk
UNIT_POST_EULERIAN_VTK_MPI_EXE := $(TESTBINDIR)/unit_post_eulerian_vtk_mpi
UNIT_POST_PARTICLE_VTK_MPI_EXE := $(TESTBINDIR)/unit_post_particle_vtk_mpi
UNIT_POST_COMPUTE_MPI_EXE := $(TESTBINDIR)/unit_post_compute_mpi
UNIT_POSTPROCESSOR_EXE := $(TESTBINDIR)/unit_postprocessor
UNIT_POST_STATISTICS_EXE := $(TESTBINDIR)/unit_post_statistics
UNIT_GRID_EXE     := $(TESTBINDIR)/unit_grid
UNIT_METRIC_EXE   := $(TESTBINDIR)/unit_metric
UNIT_BOUNDARIES_EXE := $(TESTBINDIR)/unit_boundaries
UNIT_POISSON_RHS_EXE := $(TESTBINDIR)/unit_poisson_rhs
UNIT_LES_EXE      := $(TESTBINDIR)/unit_les
UNIT_RUNTIME_EXE := $(TESTBINDIR)/unit_runtime
UNIT_MPI_EXE := $(TESTBINDIR)/unit_mpi
UNIT_PERIODIC_DEV_EXE := $(TESTBINDIR)/unit_periodic_dev
TEST_CFLAGS_TO_USE := $(CFLAGS_TO_USE) -I$(TESTCDIR)
DEPFLAGS := -MMD -MP
TEST_SUPPORT_OBJ  := $(TESTOBJDIR)/test_support.o
DOCTOR_OBJ        := $(TESTOBJDIR)/test_install_check.o
UNIT_GEOMETRY_OBJ := $(TESTOBJDIR)/test_geometry.o
UNIT_SETUP_OBJ    := $(TESTOBJDIR)/test_setup_lifecycle.o
UNIT_SOLVER_OBJ   := $(TESTOBJDIR)/test_solver_kernels.o
UNIT_NEWTON_KRYLOV_OBJ := $(TESTOBJDIR)/test_momentum_newton_krylov.o
UNIT_NEWTON_FIXEDPOINT_OBJ := $(TESTOBJDIR)/test_momentum_newton_boundary_fixedpoint.o
UNIT_CANDIDATES_EXE := $(TESTBINDIR)/unit_momentum_candidates
UNIT_CANDIDATES_OBJ := $(TESTOBJDIR)/test_momentum_convective_candidates.o
UNIT_PARTICLES_OBJ := $(TESTOBJDIR)/test_particle_kernels.o
UNIT_STATISTICS_OBJ := $(TESTOBJDIR)/test_statistics_moments.o
UNIT_STATISTICS_TARGET_OBJ := $(TESTOBJDIR)/test_statistics_target.o
UNIT_STATISTICS_WINDOW_OBJ := $(TESTOBJDIR)/test_statistics_window.o
UNIT_STATISTICS_ACC_OBJ := $(TESTOBJDIR)/test_statistics_accumulator.o
UNIT_STATISTICS_CFG_OBJ := $(TESTOBJDIR)/test_statistics_config.o
UNIT_IO_OBJ       := $(TESTOBJDIR)/test_io.o
UNIT_LOGGING_OBJ  := $(TESTOBJDIR)/test_logging.o
UNIT_POST_OBJ     := $(TESTOBJDIR)/test_postprocessing.o
UNIT_POST_VTK_OBJ := $(TESTOBJDIR)/test_vtk_io.o
UNIT_POST_EULERIAN_VTK_MPI_OBJ := $(TESTOBJDIR)/test_post_eulerian_vtk_mpi.o
UNIT_POST_PARTICLE_VTK_MPI_OBJ := $(TESTOBJDIR)/test_post_particle_vtk_mpi.o
UNIT_POST_COMPUTE_MPI_OBJ := $(TESTOBJDIR)/test_post_compute_mpi.o
UNIT_POSTPROCESSOR_OBJ := $(TESTOBJDIR)/test_postprocessor.o
UNIT_POST_STATISTICS_OBJ := $(TESTOBJDIR)/test_statistics.o
TEST_POSTPROCESSOR_IMPL_OBJ := $(TESTOBJDIR)/postprocessor_no_main.o
UNIT_GRID_OBJ     := $(TESTOBJDIR)/test_grid.o
UNIT_METRIC_OBJ   := $(TESTOBJDIR)/test_metric.o
UNIT_BOUNDARIES_OBJ := $(TESTOBJDIR)/test_boundaries.o
UNIT_POISSON_RHS_OBJ := $(TESTOBJDIR)/test_poisson_rhs.o
UNIT_LES_OBJ      := $(TESTOBJDIR)/test_les.o
UNIT_RUNTIME_OBJ := $(TESTOBJDIR)/test_runtime_kernels.o
UNIT_MPI_OBJ := $(TESTOBJDIR)/test_mpi_kernels.o
UNIT_PERIODIC_DEV_OBJ := $(TESTOBJDIR)/test_periodic_dev.o
TEST_COMMON_OBJS  := $(sort $(filter-out $(OBJDIR)/simulator.o $(OBJDIR)/postprocessor.o,$(SIMULATOR_OBJS) $(POSTPROCESSOR_OBJS)))
DEPFILES := $(wildcard $(OBJDIR)/*.d) $(wildcard $(TESTOBJDIR)/*.d)

# ==============================================================================
# --- 5. Build Targets & Rules ---
# ==============================================================================
.PHONY: all simulator postprocessor conductor dirs FORCE

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

# This rule installs a small launcher so PATH-based access can use the managed
# Python environment created by bootstrap_install.sh while keeping
# picurv_cli/picurv as the stable source-tree entrypoint.
$(CONDUCTOR_EXE): FORCE $(PICURV_ENTRYPOINT) | dirs
	@echo "--- Installing Conductor Script: $(@) ---"
	@rm -f $@
	@{ \
	  echo '#!/usr/bin/env bash'; \
	  echo 'set -e'; \
	  echo 'PICURV_DIR="$$(cd "$$(dirname "$${BASH_SOURCE[0]}")/.." && pwd)"'; \
	  echo '_picurv_apply_managed_python_env() {'; \
	  echo '  if [ -f "$$PICURV_DIR/.picurv-python-env" ]; then'; \
	  echo '    . "$$PICURV_DIR/.picurv-python-env"'; \
	  echo '  fi'; \
	  echo '  export PYTHONNOUSERSITE=1'; \
	  echo '  unset PYTHONPATH'; \
	  echo '  if [ -n "$${PICURV_PYTHON_LD_LIBRARY_PATH:-}" ]; then'; \
	  echo '    export LD_LIBRARY_PATH="$$PICURV_PYTHON_LD_LIBRARY_PATH$${LD_LIBRARY_PATH:+:$$LD_LIBRARY_PATH}"'; \
	  echo '  fi'; \
	  echo '}'; \
	  echo 'if [ -x "$$PICURV_DIR/.picurv-venv/bin/python" ]; then'; \
	  echo '  _picurv_apply_managed_python_env'; \
	  echo '  exec "$$PICURV_DIR/.picurv-venv/bin/python" "$$PICURV_DIR/picurv_cli/picurv" "$$@"'; \
	  echo 'fi'; \
	  echo 'if [ -f "$$PICURV_DIR/.picurv-python" ]; then'; \
	  echo '  _picurv_python="$$(sed -n "1p" "$$PICURV_DIR/.picurv-python")"'; \
	  echo '  if [ -n "$$_picurv_python" ] && [ -x "$$_picurv_python" ]; then'; \
	  echo '    case "$$_picurv_python" in */.picurv-venv/bin/python) _picurv_apply_managed_python_env ;; esac'; \
	  echo '    exec "$$_picurv_python" "$$PICURV_DIR/picurv_cli/picurv" "$$@"'; \
	  echo '  fi'; \
	  echo 'fi'; \
	  echo 'if [ -n "$${PICURV_PYTHON:-}" ]; then'; \
	  echo '  exec "$$PICURV_PYTHON" "$$PICURV_DIR/picurv_cli/picurv" "$$@"'; \
	  echo 'fi'; \
	  echo 'exec python3 "$$PICURV_DIR/picurv_cli/picurv" "$$@"'; \
	} > $@
	@chmod +x $@

# Generic rule for compiling any .c file from SRCDIR into an object file in OBJDIR.
$(OBJDIR)/%.o: $(SRCDIR)/%.c | dirs
	@echo "--- Compiling: $< ---"
	$(CC_TO_USE) $(CFLAGS_TO_USE) $(DEPFLAGS) -c $< -o $@

# Generic rule for compiling any test .c file from TESTCDIR into TESTOBJDIR.
$(TESTOBJDIR)/%.o: $(TESTCDIR)/%.c | dirs
	@echo "--- Compiling Test: $< ---"
	$(CC_TO_USE) $(TEST_CFLAGS_TO_USE) $(DEPFLAGS) -c $< -o $@

# Test-only object build for postprocessor orchestration APIs without the executable main().
$(TEST_POSTPROCESSOR_IMPL_OBJ): $(SRCDIR)/postprocessor.c | dirs
	@echo "--- Compiling Test Support Source: $< (no main) ---"
	$(CC_TO_USE) $(TEST_CFLAGS_TO_USE) $(DEPFLAGS) -DPICURV_POSTPROCESSOR_NO_MAIN -c $< -o $@

$(DOCTOR_EXE): $(DOCTOR_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_GEOMETRY_EXE): $(UNIT_GEOMETRY_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_SETUP_EXE): $(UNIT_SETUP_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_SOLVER_EXE): $(UNIT_SOLVER_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

# The focused test includes a renamed copy for private-callback checks and also
# links the production object for public lifecycle/coverage checks.
$(UNIT_NEWTON_KRYLOV_EXE): $(UNIT_NEWTON_KRYLOV_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_NEWTON_FIXEDPOINT_EXE): $(UNIT_NEWTON_FIXEDPOINT_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Regression Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_CANDIDATES_EXE): $(UNIT_CANDIDATES_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_PARTICLES_EXE): $(UNIT_PARTICLES_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_STATISTICS_ACC_EXE): $(UNIT_STATISTICS_ACC_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_STATISTICS_CFG_EXE): $(UNIT_STATISTICS_CFG_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_STATISTICS_WINDOW_EXE): $(UNIT_STATISTICS_WINDOW_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_STATISTICS_TARGET_EXE): $(UNIT_STATISTICS_TARGET_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_STATISTICS_EXE): $(UNIT_STATISTICS_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_IO_EXE): $(UNIT_IO_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_LOGGING_EXE): $(UNIT_LOGGING_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_EXE): $(UNIT_POST_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_VTK_EXE): $(UNIT_POST_VTK_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_EULERIAN_VTK_MPI_EXE): $(UNIT_POST_EULERIAN_VTK_MPI_OBJ) $(TEST_POSTPROCESSOR_IMPL_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_PARTICLE_VTK_MPI_EXE): $(UNIT_POST_PARTICLE_VTK_MPI_OBJ) $(TEST_POSTPROCESSOR_IMPL_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_COMPUTE_MPI_EXE): $(UNIT_POST_COMPUTE_MPI_OBJ) $(TEST_POSTPROCESSOR_IMPL_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POSTPROCESSOR_EXE): $(UNIT_POSTPROCESSOR_OBJ) $(TEST_POSTPROCESSOR_IMPL_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POST_STATISTICS_EXE): $(UNIT_POST_STATISTICS_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_GRID_EXE): $(UNIT_GRID_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_METRIC_EXE): $(UNIT_METRIC_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_BOUNDARIES_EXE): $(UNIT_BOUNDARIES_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_POISSON_RHS_EXE): $(UNIT_POISSON_RHS_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_LES_EXE): $(UNIT_LES_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_RUNTIME_EXE): $(UNIT_RUNTIME_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_MPI_EXE): $(UNIT_MPI_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

$(UNIT_PERIODIC_DEV_EXE): $(UNIT_PERIODIC_DEV_OBJ) $(TEST_SUPPORT_OBJ) $(TEST_COMMON_OBJS) | dirs
	@echo "--- Linking Test Executable: $(@) ---"
	$(LINKER_TO_USE) -o $@ $^ $(LIBS_TO_USE)

## @target dirs
## @brief (Internal) Ensures all necessary build directories exist.
dirs: 
	@mkdir -p $(OBJDIR) $(BINDIR) $(DOCSDIR) $(TESTOBJDIR) $(TESTBINDIR)

# ==============================================================================
# --- 6. Execution, Auxiliary, & Cleanup Targets ---
# ==============================================================================
.PHONY: unit-momentum-candidates unit-newton-krylov unit-momentum-newton-boundary-fixedpoint run test test-python coverage coverage-python coverage-c doctor doctor-runner install-check smoke smoke-mpi smoke-mpi-matrix smoke-stress smoke-periodic smoke-periodic-dev smoke-driven-periodic unit unit-simulation unit-geometry unit-setup unit-solver unit-particles unit-statistics unit-statistics-target unit-statistics-window unit-statistics-accumulator unit-statistics-config unit-io unit-logging unit-post unit-post-eulerian-vtk-mpi unit-post-particle-vtk-mpi unit-post-compute-mpi unit-grid unit-metric unit-boundaries unit-poisson-rhs unit-les unit-runtime unit-mpi unit-periodic unit-periodic-dev ctest ctest-geometry ctest-setup ctest-solver ctest-particles ctest-io ctest-logging ctest-post ctest-grid ctest-metric ctest-boundaries ctest-poisson-rhs ctest-runtime ctest-mpi check check-mpi check-mpi-matrix check-full check-stress audit-build build-docs docs-xref open-docs tags audit-ingress audit-docs-expansion audit-docs-site preview-docs docs-inventory docs-topology audit-capability audit-contracts audit-path-literals audit-family-census audit-page-types audit-inline-choices audit-subsystems audit-freshness audit-agent-setup sync-agent-skills docs-cli-reference audit-cli-reference scaffold review-packet certify-docs certify-docs-fast install-git-hooks clean-project cleanobj clean-project-docs clean-project-tags clean-unit

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

## @target coverage-python
## @brief Runs dependency-free Python line-coverage gate for core runtime scripts.
coverage-python:
	@python3 tests/tooling/python_coverage_gate.py --min-line "$(PY_COVERAGE_MIN)" --output-dir "coverage/python"

## @target coverage-c
## @brief Rebuilds with gcov flags, runs unit+smoke plus MPI/stress validation layers, and enforces C line-coverage threshold.
coverage-c:
	@$(MAKE) --no-print-directory cleanobj clean-unit SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory COVERAGE=1 unit smoke unit-mpi smoke-mpi smoke-mpi-matrix smoke-stress SYSTEM=$(SYSTEM)
	@python3 tests/tooling/c_coverage_gate.py --src-dir "$(SRCDIR)" --obj-dir "$(OBJDIR)" --output-dir "coverage/c" --min-line "$(C_COVERAGE_MIN)"

## @target coverage
## @brief Runs Python and C coverage gates.
coverage: coverage-python coverage-c

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
	$(if $(or $(filter n,$(MAKEFLAGS)),$(findstring --just-print,$(MAKEFLAGS)),$(findstring --dry-run,$(MAKEFLAGS)),$(findstring --recon,$(MAKEFLAGS))),@echo "$(MAKE) --no-print-directory doctor-runner SYSTEM=$(SYSTEM)",@$(MAKE) --no-print-directory doctor-runner SYSTEM=$(SYSTEM))

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

## @target unit-setup
## @brief Runs setup, lifecycle, and cleanup C unit tests.
unit-setup: $(UNIT_SETUP_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-solver
## @brief Runs the solver utility C unit tests.
unit-solver: $(UNIT_SOLVER_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-newton-krylov
## @brief Runs focused matrix-free Newton--Krylov momentum tests on one and multiple ranks.
unit-newton-krylov: $(UNIT_NEWTON_KRYLOV_EXE)
	@$(MPI_LAUNCHER) -n 1 $< $(TEST_ARGS)
	@$(MPI_LAUNCHER) -n $(TEST_MPI_NPROCS) $< $(TEST_ARGS)

## @target unit-momentum-newton-boundary-fixedpoint
## @brief Runs the opt-in Newton--Krylov production integration regression on one and four ranks.
## @details Heavy fixture (full step-1 solve + step-2 solve + projection); kept out of
##          `unit`/`check`. The cheap seed-removal detector lives in unit-newton-krylov.
unit-momentum-newton-boundary-fixedpoint: $(UNIT_NEWTON_FIXEDPOINT_EXE)
	@status=0; \
	$(MPI_LAUNCHER) -n 1 $< $(TEST_ARGS) || status=$$?; \
	$(MPI_LAUNCHER) -n 4 $< $(TEST_ARGS) || status=$$?; \
	exit $$status

## @target unit-momentum-candidates
## @brief Runs the A4a convective-candidate study (FD Jacobian + RK; states A-C).
unit-momentum-candidates: $(UNIT_CANDIDATES_EXE)
	@ref=$$(mktemp /tmp/picurv_unit_momentum_candidates_ref.XXXXXX); \
	token=$$(date +%s).$$$$; \
	status=0; \
	$(MPI_LAUNCHER) -n 1 $< -candidate_ref_path $$ref -candidate_ref_token $$token || status=$$?; \
	if [ $$status -eq 0 ]; then \
	  $(MPI_LAUNCHER) -n $(TEST_MPI_NPROCS) $< -candidate_ref_path $$ref -candidate_ref_token $$token || status=$$?; \
	fi; \
	rm -f $$ref; \
	exit $$status

## @target unit-particles
## @brief Runs the particle kernel C unit tests.
unit-particles: $(UNIT_PARTICLES_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-statistics
## @brief Runs the centered-moment kernel C unit tests.
unit-statistics: $(UNIT_STATISTICS_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-statistics-target
## @brief Runs the spatial-target resolution C unit tests.
unit-statistics-target: $(UNIT_STATISTICS_TARGET_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-statistics-window
## @brief Runs the statistics window lifecycle and scheduling C unit tests.
unit-statistics-window: $(UNIT_STATISTICS_WINDOW_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-statistics-accumulator
## @brief Runs the per-window accumulator storage and application C unit tests.
unit-statistics-accumulator: $(UNIT_STATISTICS_ACC_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-statistics-config
## @brief Runs the field-statistics control ingress C unit tests.
unit-statistics-config: $(UNIT_STATISTICS_CFG_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-io
## @brief Runs the I/O infrastructure C unit tests.
unit-io: $(UNIT_IO_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-logging
## @brief Runs the logging infrastructure C unit tests.
unit-logging: $(UNIT_LOGGING_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-post
## @brief Runs post-processing kernels, VTK I/O, postprocessor orchestration, and statistics C unit tests.
unit-post: $(UNIT_POST_EXE) $(UNIT_POST_VTK_EXE) $(UNIT_POSTPROCESSOR_EXE) $(UNIT_POST_STATISTICS_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $(UNIT_POST_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $(UNIT_POST_VTK_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $(UNIT_POSTPROCESSOR_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $(UNIT_POST_STATISTICS_EXE)

## @target unit-grid
## @brief Runs the grid utility C unit tests.
unit-grid: $(UNIT_GRID_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-metric
## @brief Runs the metric kernel C unit tests.
unit-metric: $(UNIT_METRIC_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-boundaries
## @brief Runs boundary-system focused C unit tests.
unit-boundaries: $(UNIT_BOUNDARIES_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-poisson-rhs
## @brief Runs focused Poisson/RHS kernel C unit tests.
unit-poisson-rhs: $(UNIT_POISSON_RHS_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-les
## @brief Runs the LES subgrid-scale closure C unit tests.
unit-les: $(UNIT_LES_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-runtime
## @brief Runs focused runtime-kernel C unit tests.
unit-runtime: $(UNIT_RUNTIME_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-mpi
## @brief Runs dedicated multi-rank MPI-focused C unit tests.
unit-mpi: $(UNIT_MPI_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_MPI_NPROCS) $<
	@$(MAKE) --no-print-directory unit-post-compute-mpi SYSTEM=$(SYSTEM) TEST_MPI_NPROCS=$(TEST_MPI_NPROCS)
	@$(MAKE) --no-print-directory unit-post-eulerian-vtk-mpi SYSTEM=$(SYSTEM) TEST_MPI_NPROCS=$(TEST_MPI_NPROCS)
	@$(MAKE) --no-print-directory unit-post-particle-vtk-mpi SYSTEM=$(SYSTEM) TEST_MPI_NPROCS=$(TEST_MPI_NPROCS)

## @target unit-post-compute-mpi
## @brief Checks post-processing compute pipelines against analytic truth in serial and MPI.
unit-post-compute-mpi: $(UNIT_POST_COMPUTE_MPI_EXE)
	@$(MPI_LAUNCHER) -n 1 $<
	@$(MPI_LAUNCHER) -n $(TEST_MPI_NPROCS) $<

## @target unit-post-eulerian-vtk-mpi
## @brief Exhaustively compares serial and multirank instantaneous Eulerian VTK output.
unit-post-eulerian-vtk-mpi: $(UNIT_POST_EULERIAN_VTK_MPI_EXE)
	@tmpdir=$$(mktemp -d /tmp/picurv-eulerian-vtk-mpi.XXXXXX); \
	status=0; \
	mkdir -p $$tmpdir/serial $$tmpdir/parallel; \
	$(MPI_LAUNCHER) -n 1 $< -vtk_test_output_prefix $$tmpdir/serial/field || status=$$?; \
	if [ $$status -eq 0 ]; then \
	  $(MPI_LAUNCHER) -n $(TEST_MPI_NPROCS) $< -vtk_test_output_prefix $$tmpdir/parallel/field || status=$$?; \
	fi; \
	if [ $$status -eq 0 ] && ! cmp -s $$tmpdir/serial/field_00003.vts $$tmpdir/parallel/field_00003.vts; then \
	  echo "Serial and $(TEST_MPI_NPROCS)-rank Eulerian VTK files differ."; \
	  status=1; \
	fi; \
	if [ $$status -eq 0 ]; then \
	  echo "Serial and $(TEST_MPI_NPROCS)-rank Eulerian VTK files are byte-identical."; \
	fi; \
	rm -rf $$tmpdir; \
	exit $$status

## @target unit-post-particle-vtk-mpi
## @brief Exhaustively compares serial and multirank subsampled particle VTK output.
unit-post-particle-vtk-mpi: $(UNIT_POST_PARTICLE_VTK_MPI_EXE)
	@tmpdir=$$(mktemp -d /tmp/picurv-particle-vtk-mpi.XXXXXX); \
	status=0; \
	mkdir -p $$tmpdir/serial $$tmpdir/parallel; \
	$(MPI_LAUNCHER) -n 1 $< -vtk_test_output_prefix $$tmpdir/serial/particles || status=$$?; \
	if [ $$status -eq 0 ]; then \
	  $(MPI_LAUNCHER) -n $(TEST_MPI_NPROCS) $< -vtk_test_output_prefix $$tmpdir/parallel/particles || status=$$?; \
	fi; \
	if [ $$status -eq 0 ] && ! cmp -s $$tmpdir/serial/particles_00004.vtp $$tmpdir/parallel/particles_00004.vtp; then \
	  echo "Serial and $(TEST_MPI_NPROCS)-rank particle VTK files differ."; \
	  status=1; \
	fi; \
	if [ $$status -eq 0 ]; then \
	  echo "Serial and $(TEST_MPI_NPROCS)-rank particle VTK files are byte-identical."; \
	fi; \
	rm -rf $$tmpdir; \
	exit $$status

## @target unit-periodic
## @brief Runs focused geometric-periodic boundary tests.
unit-periodic: $(UNIT_PERIODIC_DEV_EXE)
	@$(MPI_LAUNCHER) -n $(TEST_NPROCS) $<

## @target unit-periodic-dev
## @brief Compatibility alias for `unit-periodic`.
unit-periodic-dev: unit-periodic

## @target unit-simulation
## @brief Runs the simulation-core debugging suites without the setup/post/IO layers.
unit-simulation: unit-boundaries unit-solver unit-newton-krylov unit-poisson-rhs unit-runtime unit-particles

## @target unit
## @brief Runs the full isolated C unit/component suite.
unit: unit-geometry unit-setup unit-solver unit-newton-krylov unit-particles unit-statistics unit-statistics-target unit-statistics-window unit-statistics-accumulator unit-statistics-config unit-io unit-logging unit-post unit-grid unit-metric unit-boundaries unit-poisson-rhs unit-les unit-runtime unit-periodic

## @target ctest
## @brief Compatibility alias for `unit`.
ctest: unit

## @target ctest-geometry
## @brief Compatibility alias for `unit-geometry`.
ctest-geometry: unit-geometry

## @target ctest-setup
## @brief Compatibility alias for `unit-setup`.
ctest-setup: unit-setup

## @target ctest-solver
## @brief Compatibility alias for `unit-solver`.
ctest-solver: unit-solver

## @target ctest-particles
## @brief Compatibility alias for `unit-particles`.
ctest-particles: unit-particles

## @target ctest-io
## @brief Compatibility alias for `unit-io`.
ctest-io: unit-io

## @target ctest-logging
## @brief Compatibility alias for `unit-logging`.
ctest-logging: unit-logging

## @target ctest-post
## @brief Compatibility alias for `unit-post`.
ctest-post: unit-post

## @target ctest-grid
## @brief Compatibility alias for `unit-grid`.
ctest-grid: unit-grid

## @target ctest-metric
## @brief Compatibility alias for `unit-metric`.
ctest-metric: unit-metric

## @target ctest-boundaries
## @brief Compatibility alias for `unit-boundaries`.
ctest-boundaries: unit-boundaries

## @target ctest-poisson-rhs
## @brief Compatibility alias for `unit-poisson-rhs`.
ctest-poisson-rhs: unit-poisson-rhs

## @target ctest-runtime
## @brief Compatibility alias for `unit-runtime`.
ctest-runtime: unit-runtime

## @target ctest-mpi
## @brief Compatibility alias for `unit-mpi`.
ctest-mpi: unit-mpi

## @target smoke
## @brief Runs executable-level smoke checks (template matrix + tiny flat/bent/brownian runtime flows).
smoke: simulator postprocessor conductor
	@bash $(SMOKE_RUNNER) "$(SIMULATOR_EXE)" "$(POSTPROCESSOR_EXE)" "$(MPI_LAUNCHER)" "$(TEST_NPROCS)"

## @target smoke-mpi
## @brief Runs executable-level multi-rank smoke checks (flat+bent tiny runtime flows).
smoke-mpi: simulator postprocessor conductor
	@bash $(SMOKE_RUNNER) "$(SIMULATOR_EXE)" "$(POSTPROCESSOR_EXE)" "$(MPI_LAUNCHER)" "$(SMOKE_MPI_NPROCS)"

## @target smoke-mpi-matrix
## @brief Runs multi-rank runtime smoke checks across a rank matrix (default 2 and 3).
smoke-mpi-matrix: simulator postprocessor conductor
	@set -e; \
	for n in $(SMOKE_MPI_MATRIX_NPROCS); do \
		echo "==> smoke-mpi-matrix rank=$$n"; \
		bash $(SMOKE_RUNNER) "$(SIMULATOR_EXE)" "$(POSTPROCESSOR_EXE)" "$(MPI_LAUNCHER)" "$$n"; \
	done

## @target smoke-stress
## @brief Runs opt-in medium-budget smoke stress checks that extend, but do not replace, the default smoke gate.
smoke-stress: simulator postprocessor conductor
	@bash $(SMOKE_RUNNER) "$(SIMULATOR_EXE)" "$(POSTPROCESSOR_EXE)" "$(MPI_LAUNCHER)" "$(TEST_NPROCS)" "stress"

## @target smoke-periodic
## @brief Runs the geometric-periodic runtime smoke harness.
smoke-periodic: simulator postprocessor conductor
	@bash $(SMOKE_RUNNER) "$(SIMULATOR_EXE)" "$(POSTPROCESSOR_EXE)" "$(MPI_LAUNCHER)" "$(TEST_NPROCS)" "periodic"

## @target smoke-periodic-dev
## @brief Compatibility alias for `smoke-periodic`.
smoke-periodic-dev: smoke-periodic

## @target smoke-driven-periodic
## @brief Asserts the multigrid coarse-solve residual invariant and the driven-flux contract.
## @details Runs at the rank counts in DRIVEN_PERIODIC_NPROCS (default "4 10"), two counts
##          that previously disagreed on tracked-vs-true Poisson residuals. Also checks that
##          the initial_flux target matches the analytic initial-condition flux and survives
##          a restart. Needs enough cores for the largest rank count.
smoke-driven-periodic: simulator conductor
	@bash $(DRIVEN_PERIODIC_RUNNER) "$(SIMULATOR_EXE)" "$(MPI_LAUNCHER)" $(DRIVEN_PERIODIC_NPROCS)

## @target check
## @brief Runs the full local validation sweep (Python, doctor, unit, smoke).
check:
	@$(MAKE) --no-print-directory test-python
	@$(MAKE) --no-print-directory doctor-runner SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory unit SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory smoke SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory smoke-periodic SYSTEM=$(SYSTEM)

## @target check-mpi
## @brief Runs `check` plus dedicated multi-rank MPI tests.
check-mpi: check
	@$(MAKE) --no-print-directory smoke-mpi SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory unit-mpi SYSTEM=$(SYSTEM)

## @target check-mpi-matrix
## @brief Runs `check` plus multi-rank matrix smoke and dedicated MPI unit tests.
check-mpi-matrix: check
	@$(MAKE) --no-print-directory smoke-mpi-matrix SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory unit-mpi SYSTEM=$(SYSTEM)

## @target check-full
## @brief Runs `check` plus all MPI validation layers (`unit-mpi`, `smoke-mpi`, `smoke-mpi-matrix`).
check-full: check
	@$(MAKE) --no-print-directory unit-mpi SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory smoke-mpi SYSTEM=$(SYSTEM)
	@$(MAKE) --no-print-directory smoke-mpi-matrix SYSTEM=$(SYSTEM)

## @target check-stress
## @brief Runs `check-full` plus the opt-in smoke stress tier.
check-stress: check-full
	@$(MAKE) --no-print-directory smoke-stress SYSTEM=$(SYSTEM)

## @target audit-build
## @brief Runs a clean compilation audit and writes full and warning-only logs under `logs/`.
audit-build:
	@mkdir -p $(LOGDIR)
	@rm -f $(BUILD_LOG) $(BUILD_WARNINGS_LOG)
	@echo "==> Running compilation audit ($(BUILD_AUDIT_GOALS))..."
	@/bin/bash -o pipefail -c '$(MAKE) --no-print-directory $(BUILD_AUDIT_GOALS) SYSTEM=$(SYSTEM) 2>&1 | tee "$(BUILD_LOG)"'
	@grep -E "warning:" "$(BUILD_LOG)" > "$(BUILD_WARNINGS_LOG)" || true
	@echo "Compilation log: $(BUILD_LOG)"
	@echo "Compilation warnings: $(BUILD_WARNINGS_LOG)"

## @target docs-xref
## @brief Generates Doxygen output and the optional dirty-byte-stamped source-reference index.
docs-xref: dirs
	@echo "==> Generating Doxygen documentation..."
	@mkdir -p $(LOGDIR)
	@command -v doxygen >/dev/null || { echo "Doxygen is required for docs-xref; install doxygen and rerun make docs-xref." >&2; exit 1; }
	@doxygen docs/Doxyfile
	@python3 tests/tooling/generate_xref_index.py --repo-root .

## @target build-docs
## @brief Generates Doxygen documentation and the optional review cross-reference index.
build-docs: docs-xref
	@python3 tests/tooling/stamp_docs_revision.py --html-dir docs_build/html
	@python3 tests/tooling/generate_doxygen_fallback_indexes.py --repo-root . --html-dir docs_build/html
	@python3 tests/tooling/inject_theme_sync.py --repo-root . --html-dir docs_build/html
	@echo "HTML documentation generated in docs_build/html/index.html"
	@echo "Doxygen warnings log: $(LOGDIR)/doxygen.warnings"

## @target certify-docs-fast
## @brief Validates documentation, public API comments, ingress, and all shipped configuration bundles for clean HEAD.
certify-docs-fast:
	@python3 tests/tooling/certify_documentation.py

## @target certify-docs
## @brief Runs the full commit-scoped documentation certification, including PETSc/MPI runtime validation.
certify-docs:
	@python3 tests/tooling/certify_documentation.py --runtime

## @target install-git-hooks
## @brief Configures this clone to use tracked Git hooks, including full certification before pushes to main.
install-git-hooks:
	@git config core.hooksPath .githooks
	@echo "Git hooks enabled from .githooks/ (main pushes run make certify-docs)."

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
	@python3 tests/tooling/audit_ingress.py

## @target audit-docs-expansion
## @brief Rejects generic documentation-expansion debris across the repository Markdown corpus.
audit-docs-expansion:
	@python3 tests/tooling/audit_generic_expansion.py

## @target docs-inventory
## @brief Regenerates the public capability inventory from the executable sources.
docs-inventory:
	@python3 tests/tooling/generate_capability_inventory.py

## @target audit-capability
## @brief Verifies capability parity across sources and Tier-2 documentation coverage.
audit-capability:
	@python3 tests/tooling/audit_capability_coverage.py

## @target docs-topology
## @brief Regenerates the run artifact topology snapshot from the CLI planner.
docs-topology:
	@python3 tests/tooling/extract_artifact_topology.py

## @target audit-contracts
## @brief Verifies enforced invariant contracts and reports tracked ones.
audit-contracts:
	@python3 tests/tooling/audit_contracts.py

## @target docs-cli-reference
## @brief Regenerates the CLI reference from the assembled argparse parser.
docs-cli-reference:
	@python3 tests/tooling/generate_cli_reference.py

## @target audit-cli-reference
## @brief Fails when the generated CLI reference no longer matches the live parser.
audit-cli-reference:
	@python3 tests/tooling/generate_cli_reference.py --check

## @target audit-freshness
## @brief Reports hard (blocking) and soft (advisory) documentation freshness suspicion.
audit-freshness:
	@python3 tests/tooling/audit_freshness.py

## @target attest-freshness
## @brief Records reviewed surfaces as current: `make attest-freshness ARGS="cli.parser"`.
attest-freshness:
	@python3 tests/tooling/audit_freshness.py --attest $(ARGS)

## @target audit-subsystems
## @brief Fails when a subsystem claims a status it has not documented.
audit-subsystems:
	@python3 tests/tooling/audit_subsystem_lifecycle.py

## @target audit-inline-choices
## @brief Rejects public closed choices written as unnamed inline literals.
audit-inline-choices:
	@python3 tests/tooling/audit_inline_choices.py

## @target audit-page-types
## @brief Fails when a published page has no declared document type.
audit-page-types:
	@python3 tests/tooling/audit_page_types.py

## @target audit-family-census
## @brief Reports public selector surfaces that no capability family covers.
audit-family-census:
	@python3 tests/tooling/audit_family_census.py

## @target audit-path-literals
## @brief Rejects unmanaged run-path literals in narrative prose.
audit-path-literals:
	@python3 tests/tooling/audit_path_literals.py

## @target audit-agent-setup
## @brief Verifies portable shared instructions, materialized skills, and local-settings hygiene.
audit-agent-setup:
	@python3 tests/tooling/audit_agent_setup.py

## @target sync-agent-skills
## @brief Materializes Claude skill copies from canonical `.agents/skills`, then audits parity.
sync-agent-skills:
	@python3 tests/tooling/audit_agent_setup.py --sync

## @target scaffold
## @brief Generates a documentation skeleton: `make scaffold ARGS="capability --family boundary.handler --value x"`.
scaffold:
	@python3 tests/tooling/scaffold_documentation.py $(ARGS)

## @target review-packet
## @brief Routes PAGE, CONTRACT, CAPABILITY, SUBSYSTEM, SURFACE, or CHANGED=working-tree.
review-packet:
	@python3 tests/tooling/review_packet.py \
		$(if $(PAGE),$(PAGE)) \
		$(if $(CONTRACT),--contract $(CONTRACT)) \
		$(if $(CAPABILITY),--capability $(CAPABILITY)) \
		$(if $(SUBSYSTEM),--subsystem $(SUBSYSTEM)) \
		$(if $(SURFACE),--surface $(SURFACE)) \
		$(if $(CHANGED),--changed $(CHANGED))

## @target preview-docs
## @brief Builds, gates, and serves the documentation site for local review.
## @details Serves the same artifact that gets published. Override the address with
##          `make preview-docs PORT=8080 BIND=0.0.0.0`.
preview-docs: build-docs
	@if [ -s $(LOGDIR)/doxygen.warnings ]; then \
		echo "Doxygen emitted warnings; refusing to preview:"; \
		cat $(LOGDIR)/doxygen.warnings; \
		exit 1; \
	fi
	@python3 tests/tooling/audit_capability_coverage.py
	@python3 tests/tooling/audit_docs_site.py
	@echo ""
	@echo "==> Serving documentation at http://localhost:$(PORT)/"
	@echo "    (bind $(BIND):$(PORT); press Ctrl-C to stop)"
	@echo ""
	@python3 -m http.server $(PORT) --bind $(BIND) --directory docs_build/html

## @target audit-docs-site
## @brief Verifies project-owned documentation URLs resolve and every navigation tab has a page.
audit-docs-site:
	@python3 tests/tooling/audit_docs_site.py

## @target clean-project
## @brief Removes all build artifacts generated by this project.
clean-project: cleanobj clean-unit clean-project-docs clean-project-tags
	@echo "Project cleaned."

## @target cleanobj
## @brief (Internal) Removes object files and compiled executables while preserving the conductor launcher.
cleanobj:
	@echo "--- Removing object files and compiled executables (preserving picurv launcher)"
	@rm -rf $(OBJDIR) 
	@rm -f $(SIMULATOR_EXE) $(POSTPROCESSOR_EXE)

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

-include $(DEPFILES)
