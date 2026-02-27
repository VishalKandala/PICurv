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

@section build_sec 6. Build with pic.flow

```bash
./scripts/pic.flow build
```

Expected binaries:

- `bin/picsolver`
- `bin/postprocessor`

Useful variants:

```bash
./scripts/pic.flow build clean-project
./scripts/pic.flow build SYSTEM=cluster
```

@section verify_sec 7. Verify Installation

```bash
./scripts/pic.flow --help
./scripts/pic.flow run --help
./scripts/pic.flow sweep --help
./scripts/pic.flow validate --help
```

If these commands work and binaries exist, your installation is ready.

@section common_sec 8. Common Installation Failures

- `PETSC_DIR`/`PETSC_ARCH` unset or mismatched with built arch.
- MPI compiler wrappers unavailable in PATH.
- PETSc configured without required downloaded dependencies.
- stale object files after toolchain changes (use `clean-project`).

For runtime-level failures after successful build, see **@subpage 39_Common_Fatal_Errors**.

@section next_steps_sec 9. Next Steps

- Quick path: **@subpage 38_Start_Here_10_Minutes**
- Full first simulation: **@subpage 02_Tutorial_Programmatic_Grid**
