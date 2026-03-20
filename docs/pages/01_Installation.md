@page 01_Installation Installation Guide

@anchor _Installation

This guide covers dependency setup, PETSc configuration, and binary build verification for PICurv.

@tableofcontents

@section p01_prereqs_sec 1. Prerequisites

Install these first:

- C compiler (`gcc` or `clang`)
- MPI implementation (`mpich` or `openmpi`)
- GNU Make
- Python 3.10+ + pip
- Git
- PETSc 3.20.3+ with `DMSwarm` support

@section p01_automated_sec 2. Automated Install (Recommended)

From the PICurv repo root:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
./scripts/bootstrap_install.sh
```

If PETSc is not installed yet, let the script build it:

```bash
./scripts/bootstrap_install.sh --install-petsc
```

The script installs system and Python dependencies, verifies PETSc/DMSwarm visibility, then builds:

- `bin/simulator` (compiled C solver)
- `bin/postprocessor` (compiled C post-processor)
- `bin/picurv` (symlink to `scripts/picurv`, the Python conductor)

@section p01_install_tools_sec 3. Install Base Toolchain (Manual Path)

Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran mpich git make pkg-config libx11-dev python3 python3-pip python3-venv
python3 -m pip install --upgrade pip
python3 -m pip install --user pyyaml numpy
```

RHEL/CentOS/Fedora:

```bash
sudo yum groupinstall -y "Development Tools"
sudo yum install -y mpich-devel python3 python3-pip git
python3 -m pip install --user pyyaml numpy
```

Optional (study plot generation):

```bash
python3 -m pip install --user matplotlib
```

@section p01_petsc_sec 4. Install PETSc

Recommended source install:

```bash
git clone -b v3.20.3 https://gitlab.com/petsc/petsc.git
cd petsc
```

Debug build example:

```bash
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            --download-fblaslapack --download-metis --download-parmetis \
            --with-dmswarm=1 --with-debugging=1
make all
make check
```

Optimized build example:

```bash
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            --download-fblaslapack --download-metis --download-parmetis \
            --with-dmswarm=1 --with-debugging=0 \
            --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3'
make all
make check
```

Official references:

- https://petsc.org/release/install/
- https://petsc.org/release/docs/manual/

@section p01_env_sec 5. Configure Environment Variables

Add to your shell profile (`~/.bashrc` or equivalent):

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
source /path/to/PICurv/etc/picurv.sh
```

The `etc/picurv.sh` script sets `PICURV_DIR`, adds `bin/` to your `PATH` for
compiled executables, and also exposes `scripts/` as a fallback so `picurv`
still resolves if `bin/picurv` is temporarily absent before a rebuild. It is
idempotent and safe to source multiple times.

Reload and verify:

```bash
source ~/.bashrc
echo "$PETSC_DIR"
echo "$PETSC_ARCH"
picurv --help
```

Verify PETSc has DMSwarm headers:

```bash
test -f "$PETSC_DIR/include/petscdmswarm.h" && echo "DMSwarm header found"
test -f "$PETSC_DIR/$PETSC_ARCH/include/petscconf.h" && echo "PETSc arch config found"
```

@section p01_clone_sec 6. Clone PICurv

```bash
git clone https://github.com/VishalKandala/PICurv.git
cd PICurv
```

@section p01_build_sec 7. Build with picurv

```bash
./scripts/picurv build
```

Expected binaries:

- `bin/simulator` (compiled C solver)
- `bin/postprocessor` (compiled C post-processor)
- `bin/picurv` (symlink → `scripts/picurv`)

Useful variants:

```bash
./scripts/picurv build clean-project
./scripts/picurv build SYSTEM=cluster
```

After `source etc/picurv.sh`, use `picurv` directly from any directory.

@section p01_verify_sec 8. Verify Installation

```bash
make doctor
```

Recommended sequence after a successful build:

1. `make doctor`
2. `make smoke` (stronger executable-level sanity check)
3. `make check` (full local validation sweep)
4. `make check-full` (comprehensive MPI-inclusive validation sweep)

What these prove:

- `make doctor` proves the local PETSc installation can build and run a minimal PETSc-backed binary.
- `make smoke` proves the compiled PICurv executables still launch.
- `make check` runs Python regressions plus PETSc-backed validation.
- `make check-full` additionally validates dedicated MPI unit coverage, fixed-size multi-rank smoke, and rank-matrix smoke in one pass.

What `make doctor` does not prove:

- it does not prove a full solver case is numerically correct
- it does not replace case-specific validation or convergence testing

@section p01_common_sec 9. Common Installation Failures

- `PETSC_DIR`/`PETSC_ARCH` unset or mismatched with built arch.
- MPI compiler wrappers unavailable in PATH.
- Python interpreter too old (PICurv expects Python 3.10+).
- PETSc configured without required downloaded dependencies.
- missing X11 dev library at link time (`cannot find -lX11` on Linux; install `libx11-dev`).
- stale object files after toolchain changes (use `clean-project`).

For runtime-level failures after successful build, see **@subpage 39_Common_Fatal_Errors**.

For the full testing model after installation, see **@subpage 40_Testing_and_Quality_Guide**.

@section p01_next_steps_sec 10. Next Steps

- First simulation walkthrough: **@subpage 02_Tutorial_Programmatic_Grid**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Installation Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
