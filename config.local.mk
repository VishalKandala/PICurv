# ==============================================================================
#
#          Configuration for the LOCAL Development Environment
#
# ==============================================================================
#
# ## Description
#
#   This file defines all environment-specific variables for a local build.
#   It is optimized for debugging and rapid iteration. It automatically detects
#   whether the configured PETSc is a debug or optimized build and sets compiler
#   flags accordingly.
#
# ==============================================================================

# --- 1. System Identification ---
SYSTEM_NAME := Local Development Environment (e.g., WSL)

# --- 2. Prerequisite Checks & PETSc Integration ---
# Verify that PETSc environment variables are set, then include PETSc's makefiles.
ifndef PETSC_DIR
    $(error PETSC_DIR is not set. Please export it to your PETSc installation path.)
endif
ifndef PETSC_ARCH
    $(error PETSC_ARCH is not set. Please export your PETSc build architecture.)
endif
include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules
$(info Using PETSC_ARCH: $(PETSC_ARCH))

$(info ------------------------------------------------)
#$(info DEBUG: Value of PCC_FLAGS is [$(PCC_FLAGS)])
#$(info DEBUG: Value of CLINKER is [$(CLINKER)])
#$(info DEBUG: Value of PETSC_LIB is [$(PETSC_LIB)])
#$(info ------------------------------------------------)

# ==============================================================================

# --- 3. Automatic Build Type Detection ---
# Set optimization flags based on the name of the PETSc architecture.
ifeq (,$(findstring debug,$(PETSC_ARCH)))
    # If "debug" is NOT in the PETSC_ARCH name, assume an optimized build.
    $(info PETSc performance build detected. Using debug flags: -O0)
    OPT_FLAGS := -O0 # Change to -O3 for performance builds, for now -O0 for easier debugging.
else
    # If "debug" is in the name, use flags for debugging.
    $(info PETSc debug build detected. Using debug flags: -g -O0)
    OPT_FLAGS := -g -O0
endif

# --- 4. Generic Variable Definitions for the Main Makefile ---
# These variables are used by the generic rules in the main Makefile.
CC_TO_USE     := $(PCC)
# This is the corrected line. We manually add the missing PETSc include path.
# Use PETSc-provided include flags (covers $PETSC_DIR/include and arch/include)
CFLAGS_TO_USE := $(PCC_FLAGS) $(PETSC_CC_INCLUDES) -Wall $(OPT_FLAGS) -I$(INCDIR)
#CFLAGS_TO_USE := $(PCC_FLAGS) -Wall $(OPT_FLAGS) -I$(INCDIR)
LINKER_TO_USE := $(CLINKER)
LIBS_TO_USE   := $(PETSC_LIB)

# --- 5. Execution Environment ---
MPI_LAUNCHER  := mpiexec
NPROCS ?= 2

# --- 6. Optional Build Add-ons (controlled from command line) ---
# Enable with 'make SAN=1'
ifdef SAN
    $(info Enabling Address and Undefined Behavior Sanitizers)
    SANFLAGS      := -fsanitize=address,undefined -fno-omit-frame-pointer
    CFLAGS_TO_USE += $(SANFLAGS)
    LIBS_TO_USE   += $(SANFLAGS)
endif

# Enable with 'make TEC360HOME=/path/to/tecplot'
ifdef TEC360HOME
    $(info Enabling Tecplot I/O Support)
    CFLAGS_TO_USE += -I$(TEC360HOME)/include -DTECIO=1
    LIBS_TO_USE   += $(TEC360HOME)/lib/tecio64.a -lstdc++
endif