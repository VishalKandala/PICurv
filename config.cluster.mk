# ==============================================================================
#
#            Configuration for the HPC CLUSTER Environment
#
# ==============================================================================
#
# ## Description
#
#   This file defines all environment-specific variables for a production build
#   on a High-Performance Computing (HPC) cluster. It is optimized for maximum
#   performance and can be easily extended to link against cluster-specific
#   scientific libraries.
#
# ==============================================================================

# --- 1. System Identification ---
SYSTEM_NAME := HPC Cluster Environment

# --- 2. Prerequisite Checks & PETSc Integration ---
# Verify PETSC_DIR is set, then include a compatibility layer that supports
# both old-style (PETSC_ARCH) and new-style (installed prefix) PETSc installs.
ifndef PETSC_DIR
    $(error PETSC_DIR is not set. Did you forget to 'module load petsc' or export PETSC_DIR?)
endif
include config.petsc.mk

# --- 3. Automatic Build Type Detection ---
# Preserve original behavior when PETSC_ARCH exists; add minimal fallback when it doesn't.
ifdef PETSC_ARCH
ifeq (,$(findstring debug,$(PETSC_ARCH)))
    # If "debug" is NOT in the name, use aggressive, architecture-specific optimizations.
    $(info PETSc performance build detected. Using optimization flags: -O3 -march=native)
    OPT_FLAGS := -O3 -march=native
else
    # If "debug" is in the name, use flags for debugging on the cluster.
    $(info PETSc debug build detected. Using debug flags: -g -O0)
    OPT_FLAGS := -g -O0
endif
else
    # New-style install: no PETSC_ARCH exists. Choose the same cluster default you prefer.
    $(info PETSC_ARCH not set (new-style PETSc). Defaulting to: -O3 -march=native)
    OPT_FLAGS := -O3 -march=native
endif

# --- 4. Generic Variable Definitions for the Main Makefile ---
CC_TO_USE     := $(PCC)
CFLAGS_TO_USE := $(PCC_FLAGS) $(PETSC_CC_INCLUDES) -Wall $(OPT_FLAGS) -I$(INCDIR)
LINKER_TO_USE := $(CLINKER)
LIBS_TO_USE   := $(PETSC_LIB)

# --- 5. Execution Environment ---
# Use the cluster's job scheduler launcher (e.g., srun for Slurm).
MPI_LAUNCHER  := mpirun  -mca pml ucx   -mca btl '^uct,ofi'   -mca mtl '^ofi' # Change to 'srun' if using Slurm and preferred by the cluster.
NPROCS ?= 4

# --- 6. Platform-Specific Extensions ---
# This is where you can link against libraries provided by the cluster via modules.
#
# ## Example: Link against a cluster-provided high-performance HDF5 library.
# HDF5_PATH := /path/to/cluster/hdf5-module
# $(info Adding cluster-specific HDF5 library from $(HDF5_PATH))
# CFLAGS_TO_USE += -I$(HDF5_PATH)/include
# LIBS_TO_USE   += -L$(HDF5_PATH)/lib -lhdf5_parallel

