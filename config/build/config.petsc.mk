# ==============================================================================
# config.petsc.mk — PETSc compatibility layer (old + new installs)
#
# Goal: make compilation/linking work on BOTH:
#   (A) Old-style in-tree PETSc builds that use PETSC_ARCH
#   (B) New-style installed-prefix PETSc (EasyBuild/Spack/system) with no PETSC_ARCH
#
# Design constraints:
#   - Do NOT set or modify optimization/debug flags (OPT_FLAGS etc.)
#   - Preserve old-style behavior exactly when PETSC_ARCH is set
#   - Prefer PETSc's own "variables" file when available (most faithful link line)
#   - Fall back to pkg-config only when needed
#
# Exports the same variable names your config files already use:
#   PCC, CLINKER, PCC_FLAGS, PETSC_CC_INCLUDES, PETSC_LIB
# ==============================================================================

ifndef PETSC_DIR
  $(error PETSC_DIR is not set. Please export it to your PETSc installation path.)
endif

# PETSC_ARCH is optional; required only for old-style in-tree builds
PETSC_ARCH ?=

# ------------------------------------------------------------------------------
# 1) Old-style: PETSC_ARCH set  -> keep original behavior
# ------------------------------------------------------------------------------
ifdef PETSC_ARCH
  include $(PETSC_DIR)/lib/petsc/conf/variables
  include $(PETSC_DIR)/lib/petsc/conf/rules
  $(info Using PETSc old-style (PETSC_ARCH): $(PETSC_ARCH))

else
  # ----------------------------------------------------------------------------
  # 2) New-style preferred: installed-prefix variables file (most faithful)
  #    Many packaged installs ship:
  #      $PETSC_DIR/lib/petsc/conf/variables   (or lib64)
  # ----------------------------------------------------------------------------
  PETSC_VARS_PREFIX := $(firstword \
    $(wildcard $(PETSC_DIR)/lib/petsc/conf/variables) \
    $(wildcard $(PETSC_DIR)/lib64/petsc/conf/variables))

  ifneq ($(strip $(PETSC_VARS_PREFIX)),)
    include $(PETSC_VARS_PREFIX)
    # rules is not required for your Makefile, but including it is harmless if present
    PETSC_RULES_PREFIX := $(firstword \
      $(wildcard $(PETSC_DIR)/lib/petsc/conf/rules) \
      $(wildcard $(PETSC_DIR)/lib64/petsc/conf/rules))
    ifneq ($(strip $(PETSC_RULES_PREFIX)),)
      include $(PETSC_RULES_PREFIX)
    endif

    $(info Using PETSc prefix variables: $(PETSC_VARS_PREFIX))

  else
    # ----------------------------------------------------------------------------
    # 3) Fallback: pkg-config (works when variables file is not shipped)
    # ----------------------------------------------------------------------------
    HAVE_PKG_CONFIG := $(shell command -v pkg-config >/dev/null 2>&1 && echo 1 || echo 0)
    ifneq ($(HAVE_PKG_CONFIG),1)
      $(error PETSC_ARCH is not set and no PETSc variables file found under $(PETSC_DIR)/lib*/petsc/conf/. Also pkg-config not found.)
    endif

    # Help pkg-config find PETSc inside PETSC_DIR
    export PKG_CONFIG_PATH := \
      $(PETSC_DIR)/lib/pkgconfig:$(PETSC_DIR)/lib64/pkgconfig:$(PKG_CONFIG_PATH)

    HAVE_PETSC_PC := $(shell pkg-config --exists petsc >/dev/null 2>&1 && echo 1 || echo 0)
    ifneq ($(HAVE_PETSC_PC),1)
      $(error PETSC_ARCH is not set, PETSc variables file not found, and pkg-config cannot find 'petsc'. Check PETSC_DIR or PKG_CONFIG_PATH.)
    endif

    # Mirror variable names your existing configs already use
    PCC               ?= mpicc
    CLINKER           ?= mpicc
    PCC_FLAGS         := $(shell pkg-config --cflags petsc)
    PETSC_CC_INCLUDES :=
    PETSC_LIB         := $(shell pkg-config --libs petsc)

    $(info Using PETSc via pkg-config (no PETSC_ARCH))
  endif
endif

# ------------------------------------------------------------------------------
# Final sanity check: if PETSC_LIB doesn't include libm on some pkg-config installs,
# your code may still need -lm. We do NOT force -lm here to avoid changing behavior.
# Prefer using PETSc variables file (above) for faithful dependency lists.
# If you hit 'pow@@GLIBC' unresolved, append -lm in the system config.
# ------------------------------------------------------------------------------

# Minimal diagnostics (safe)
# $(info PETSC_DIR=$(PETSC_DIR))
# $(info PETSC_ARCH=$(PETSC_ARCH))
# $(info PCC=$(PCC))
# $(info CLINKER=$(CLINKER))
# $(info PCC_FLAGS=$(PCC_FLAGS))
# $(info PETSC_CC_INCLUDES=$(PETSC_CC_INCLUDES))
# $(info PETSC_LIB=$(PETSC_LIB))

