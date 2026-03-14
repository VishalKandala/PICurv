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

# --- 2. PETSc Integration ---
include config/build/config.petsc.mk
$(info ------------------------------------------------)
#$(info DEBUG: Value of PCC_FLAGS is [$(PCC_FLAGS)])
#$(info DEBUG: Value of CLINKER is [$(CLINKER)])
#$(info DEBUG: Value of PETSC_LIB is [$(PETSC_LIB)])
#$(info ------------------------------------------------)

# ==============================================================================
# --- 3. Automatic Build Type Detection ---
ifdef PETSC_ARCH
ifeq (,$(findstring debug,$(PETSC_ARCH)))
    $(info PETSc performance build detected. Using debug flags: -O0)
    OPT_FLAGS := -O0
else
    $(info PETSc debug build detected. Using debug flags: -g -O0)
    OPT_FLAGS := -g -O0
endif
else
    # New-style install: no PETSC_ARCH exists, so we choose a safe default.
    # This does NOT affect old-style behavior at all.
    $(info PETSC_ARCH not set (new-style PETSc). Defaulting to debug flags: -g -O0)
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
# Use PETSc's own mpiexec (matches the MPI that PETSc was built with).
# Falls back to system mpiexec when MPIEXEC is not set (e.g. pkg-config path).
MPI_LAUNCHER  := $(if $(MPIEXEC),$(MPIEXEC),mpiexec)
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
