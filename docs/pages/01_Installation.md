@page 01_Installation Installation Guide

This guide covers dependency setup, PETSc configuration, and binary build verification for PICurv.

@tableofcontents

@section prereqs_sec 1. Prerequisites

Install these first:

- C compiler (`gcc` or `clang`)
- MPI implementation (`mpich` or `openmpi`)
- GNU Make
- Python 3 + pip
- Git
- PETSc 3.20.3+ with `DMSwarm` support

@section install_tools_sec 2. Install Base Toolchain

Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran mpich python3 python3-pip git
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

@section petsc_sec 3. Install PETSc

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

@section env_sec 4. Configure Environment Variables

Add to your shell profile (`~/.bashrc` or equivalent):

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
```

Reload and verify:

```bash
source ~/.bashrc
echo "$PETSC_DIR"
echo "$PETSC_ARCH"
```

@section clone_sec 5. Clone PICurv

```bash
git clone https://github.com/VishalKandala/PICurv.git
cd PICurv
```

@section build_sec 6. Build with picurv

```bash
./scripts/picurv build
```

Expected binaries:

- `bin/simulator`
- `bin/postprocessor`

Useful variants:

```bash
./scripts/picurv build clean-project
./scripts/picurv build SYSTEM=cluster
```

After the build creates `bin/picurv`, prefer `./bin/picurv` for normal help, validation, and run commands.

@section verify_sec 7. Verify Installation

```bash
make doctor
```

Recommended sequence after a successful build:

1. `make doctor`
2. `make smoke` (stronger executable-level sanity check)
3. `make check` (full local validation sweep)

What these prove:

- `make doctor` proves the local PETSc installation can build and run a minimal PETSc-backed binary.
- `make smoke` proves the compiled PICurv executables still launch.
- `make check` runs Python regressions plus PETSc-backed validation.

What `make doctor` does not prove:

- it does not prove a full solver case is numerically correct
- it does not replace case-specific validation or convergence testing

@section common_sec 8. Common Installation Failures

- `PETSC_DIR`/`PETSC_ARCH` unset or mismatched with built arch.
- MPI compiler wrappers unavailable in PATH.
- PETSc configured without required downloaded dependencies.
- stale object files after toolchain changes (use `clean-project`).

For runtime-level failures after successful build, see **@subpage 39_Common_Fatal_Errors**.

For the full testing model after installation, see **@subpage 40_Testing_and_Quality_Guide**.

@section next_steps_sec 9. Next Steps

- First simulation walkthrough: **@subpage 02_Tutorial_Programmatic_Grid**
