# This Makefile builds the project executables, organizes the build process,
# and includes custom targets for cleaning and generating TAGS for code navigation.

# -----------------------------------------------------
# Directory Variables
# -----------------------------------------------------
SRCDIR    = src
INCDIR    = include
OBJDIR    = obj
BINDIR    = bin
SCRIPTDIR = scripts
DOCSDIR    = docs

# -----------------------------------------------------
# Compiler and Flags
# -----------------------------------------------------
CC        = mpicc
CFLAGS    = -Wall -g -O0 -I$(INCDIR)
SANFLAGS  = -O1 -fsanitize=address,undefined -fno-omit-frame-pointer

# Update LDFLAGS with correct rpath and rpath-link flags based on module paths
#LDFLAGS   = -Wl,-rpath,/sw/eb/sw/Hypre/2.28.0-foss-2022b \
#            -Wl,-rpath,/sw/eb/sw/SuperLU_DIST/8.1.2-foss-2022b \
#            -Wl,-rpath-link,/sw/eb/sw/Hypre/2.28.0-foss-2022b \
#            -Wl,-rpath-link,/sw/eb/sw/SuperLU_DIST/8.1.2-foss-2022b

#LDFLAGS   += -check-pointers=rw
LIBS      =
CLINKER   = $(CC)
PETSC_LIB =

ifdef TEC360HOME
CFLAGS   += -I${TEC360HOME}/include/ -DTECIO=1
LIBS     += ${TEC360HOME}/lib/tecio64.a -lstdc++
endif

ifdef SAN
CFLAGS += $(SANFLAGS)
LDFLAGS+= $(SANFLAGS)
endif

# -----------------------------------------------------
# PETSc Integration
# -----------------------------------------------------
ifdef PETSC_DIR
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
endif
# -----------------------------------------------------
# Source and Object Files
# -----------------------------------------------------
SOURCEC   = $(wildcard $(SRCDIR)/*.c)
OBJSC     = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCEC))
#------------------------------------------------------
# Specific object files for the 'testt' executable
#------------------------------------------------------
TESTT_OBJS = $(OBJDIR)/bcs.o $(OBJDIR)/bmv.o $(OBJDIR)/compgeom.o \
             $(OBJDIR)/ibm.o $(OBJDIR)/ibm_io.o $(OBJDIR)/init.o \
             $(OBJDIR)/main.o $(OBJDIR)/metrics.o $(OBJDIR)/poisson.o \
             $(OBJDIR)/rhs.o $(OBJDIR)/rheology.o $(OBJDIR)/variables.o \
             $(OBJDIR)/fsi.o $(OBJDIR)/implicitsolver.o $(OBJDIR)/fsi_move.o \
             $(OBJDIR)/solvers.o $(OBJDIR)/copepod.o $(OBJDIR)/fish.o \
             $(OBJDIR)/cstart.o $(OBJDIR)/spline.o $(OBJDIR)/les.o \
             $(OBJDIR)/k-omega.o $(OBJDIR)/wallfunction.o $(OBJDIR)/rhs2.o
#------------------------------------------------------
# Specific object files for the 'picsolver' executable (MINIMAL VERSION)
#------------------------------------------------------
PICSOLVER_OBJS = $(OBJDIR)/picsolver.o \
                 $(OBJDIR)/setup.o \
                 $(OBJDIR)/logging.o \
		 $(OBJDIR)/grid.o \
		 $(OBJDIR)/io.o \
		 $(OBJDIR)/Metric.o \
		 $(OBJDIR)/Boundaries.o \
		 $(OBJDIR)/wallfunction.o \
                 $(OBJDIR)/simulation.o \
                 $(OBJDIR)/walkingsearch.o \
                 $(OBJDIR)/ParticleSwarm.o \
                 $(OBJDIR)/ParticleMotion.o \
                 $(OBJDIR)/interpolation.o \
                 $(OBJDIR)/initialcondition.o $(OBJDIR)/rhs.o \
                 $(OBJDIR)/solvers.o  $(OBJDIR)/implicitsolvers.o \
                 $(OBJDIR)/poisson.o
# -----------------------------------------------------------
# Specific object files for the 'postprocessor' executable
# -----------------------------------------------------------
POSTPROCESSOR_OBJS = $(OBJDIR)/picsolver.o \
                     $(OBJDIR)/setup.o \
                     $(OBJDIR)/logging.o \
		             $(OBJDIR)/grid.o \
		             $(OBJDIR)/io.o \
		             $(OBJDIR)/Metric.o \
		             $(OBJDIR)/Boundaries.o \
		             $(OBJDIR)/wallfunction.o 
# Executable Names
# -----------------------------------------------------
TESTT_EXE         = $(BINDIR)/testt
DATA_EXE          = $(BINDIR)/data
INTTEST_EXE       = $(BINDIR)/inttest
ITFCSEARCH_EXE    = $(BINDIR)/itfcsearch
DATA_VTK_EXE      = $(BINDIR)/data_vtk
DATALIS_EXE       = $(BINDIR)/datalis
DATAFILE_EXE      = $(BINDIR)/datafile
SWARM_TEST_EXE    = $(BINDIR)/swarm_test
PICSOLVER_EXE     = $(BINDIR)/picsolver

# New Executable
POSTPROCESSOR_EXE   = $(BINDIR)/postprocessor

# -----------------------------------------------------
# Phony Targets
# -----------------------------------------------------
.PHONY: all cleanobj clean_all tags \
        testt data inttest \
        itfcsearch data_vtk datalis datafile swarm_test \
        postprocessor picsolver

# -----------------------------------------------------
# Default Target
# -----------------------------------------------------
all: testt data inttest \
     itfcsearch data_vtk datalis datafile swarm_test \
     postprocessor picsolver

# -----------------------------------------------------
# Create Necessary Directories
# -----------------------------------------------------
dirs: | $(OBJDIR) $(BINDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(BINDIR):
	@mkdir -p $(BINDIR)

# -----------------------------------------------------
# Compilation Rule
# -----------------------------------------------------
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@
#------------------------------------------------------
# -----------------------------------------------------
# 1. testt executable
# -----------------------------------------------------
#testt: dirs $(TESTT_EXE)
testt: dirs $(TESTT_EXE)

$(TESTT_EXE): $(TESTT_OBJS)
	@echo "================== DEBUGGING testt =================="
	@echo "Prerequisites for linking are ($^):"
	@echo $^
	@echo "====================================================="
	$(CLINKER) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

#$(TESTT_EXE): $(TESTT_OBJS)
#	$(CLINKER) $(CFLAGS) -o $@ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 2. data executable
# -----------------------------------------------------
data: dirs $(DATA_EXE)

$(DATA_EXE): $(OBJDIR)/variables.o $(OBJDIR)/compgeom.o $(OBJDIR)/data_ibm.o \
             $(OBJDIR)/ibm_io.o $(OBJDIR)/fsi.o $(OBJDIR)/fsi_move.o \
             $(OBJDIR)/fish.o $(OBJDIR)/data.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(PETSC_SNES_LIB) $(PETSC_TS_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 3. inttest (swarm_interp) executable
# -----------------------------------------------------
inttest: dirs $(INTTEST_EXE)

$(INTTEST_EXE): $(OBJDIR)/inttest.o $(OBJDIR)/interpolation.o $(OBJDIR)/walkingsearch.o \
                $(OBJDIR)/ParticleSwarm.o $(OBJDIR)/logging.o $(OBJDIR)/setup.o \
                $(OBJDIR)/AnalyticalSolution.o $(OBJDIR)/grid.o $(OBJDIR)/io.o $(OBJDIR)/ParticleMotion.o \
                $(OBJDIR)/Metric.o $(OBJDIR)/Boundaries.o $(OBJDIR)/BC_Handlers.o $(OBJDIR)/simulation.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 4. itfcsearch
# -----------------------------------------------------
itfcsearch: dirs $(ITFCSEARCH_EXE)

$(ITFCSEARCH_EXE): $(OBJDIR)/itfcsearch.o \
                   $(OBJDIR)/walkingsearch.o $(OBJDIR)/logging.o \
                   $(OBJDIR)/grid.o $(OBJDIR)/io.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 5. data_vtk
# -----------------------------------------------------
data_vtk: dirs $(DATA_VTK_EXE)

$(DATA_VTK_EXE): $(OBJDIR)/data_vtk.o $(OBJDIR)/logging.o \
                 $(OBJDIR)/io.o $(OBJDIR)/grid.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 6. datalis
# -----------------------------------------------------
datalis: dirs $(DATALIS_EXE)

$(DATALIS_EXE): $(OBJDIR)/datalis.o $(OBJDIR)/logging.o $(OBJDIR)/io.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 7. datafile
# -----------------------------------------------------
datafile: dirs $(DATAFILE_EXE)

$(DATAFILE_EXE): $(OBJDIR)/datafile.o $(OBJDIR)/logging.o $(OBJDIR)/io.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 8. swarm_test
# -----------------------------------------------------
swarm_test: dirs $(SWARM_TEST_EXE)

$(SWARM_TEST_EXE): $(OBJDIR)/swarm_test.o $(OBJDIR)/logging.o \
                   $(OBJDIR)/io.o $(OBJDIR)/grid.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# NEW: 9. postprocess executable
# -----------------------------------------------------
postprocessor: dirs $(POSTPROCESSOR_EXE)

$(POSTPROCESSOR_EXE): $(POSTPROCESSOR_OBJS)
    @echo "--- Builing Executable: $(@) ---"
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)
    @echo "--- Build Complete: $(@) ---"
# -----------------------------------------------------
# NEW: picsolver executable
# -----------------------------------------------------
picsolver: dirs $(PICSOLVER_EXE)

$(PICSOLVER_EXE): $(PICSOLVER_OBJS)
    @echo "--- Building Executable: $(@) ---"
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)
	@echo "--- Build Complete: $(@) ---"

# -----------------------------------------------------
# Custom Targets
# -----------------------------------------------------
.PHONY: tags
tags:
	find $(SRCDIR) $(INCDIR) $(SCRIPTDIR) -type f \( -name "*.c" -o -name "*.h" -o -name "*.py" \) -print | etags -f TAGS -

.PHONY: clean_tags
clean_tags:
	rm -f TAGS

cleanobj:
	rm -f $(OBJDIR)/*.o

clean_all: cleanobj
	rm -f $(BINDIR)/* TAGS

# -------------------------------
# Documentation (Doxygen) targets
# -------------------------------
DOXYGEN   ?= doxygen
DOCS_OUT  := $(DOCSDIR)/docs_build/html

.PHONY: docs open-docs clean-docs

# Ensure docs dir exists (in case it's missing)
$(DOCSDIR):
	@mkdir -p $(DOCSDIR)

docs: | $(DOCSDIR)
	@$(DOXYGEN) -v >/dev/null 2>&1 || { \
	  echo "Error: doxygen not found. Install it (and graphviz) then retry."; \
	  echo "  macOS:   brew install doxygen graphviz"; \
	  echo "  Ubuntu:  sudo apt-get install -y doxygen graphviz"; \
	  exit 1; }
	@echo "==> Generating Doxygen docs"
	@cd $(DOCSDIR) && $(DOXYGEN) Doxyfile
	@echo "HTML: $(DOCS_OUT)/index.html"

open-docs: docs
ifeq ($(shell uname),Darwin)
	@open $(DOCS_OUT)/index.html
else
	@xdg-open $(DOCS_OUT)/index.html || true
endif

clean-docs:
	@rm -rf $(DOCSDIR)/docs_build $(DOCSDIR)/doxygen.warnings
	@echo "Docs build artifacts removed."

# -------------------------------
# PETSc test rules (guarded)
# -------------------------------
# Include PETSc test rules only if PETSC_DIR is set and the file exists.
ifneq ($(strip $(PETSC_DIR)),)
  ifneq ("$(wildcard $(PETSC_DIR)/lib/petsc/conf/test)","")
    include $(PETSC_DIR)/lib/petsc/conf/test
  else
    $(warning PETSc test include not found at $(PETSC_DIR)/lib/petsc/conf/test; skipping)
  endif
else
  $(warning PETSC_DIR is not set; skipping PETSc test includes)
endif


# Include PETSc test rules
#include ${PETSC_DIR}/lib/petsc/conf/test
