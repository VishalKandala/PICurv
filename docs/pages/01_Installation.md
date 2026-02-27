@page 01_Installation Installation Guide

This guide covers environment setup, dependency installation, and building PICurv with the `pic.flow` conductor.

@tableofcontents

@section prereqs_sec 1. Prerequisites

Install the following first:

- C compiler (`gcc` or `clang`)
- MPI implementation (`mpich` or `openmpi`)
- GNU Make
- Python 3
- Git
- PETSc (3.20.3 or newer, with `DMSwarm` support)

@section install_tools_sec 2. Install Base Toolchain

For Debian/Ubuntu:
```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran mpich python3 python3-pip git
python3 -m pip install --user pyyaml numpy
```

For RHEL/CentOS/Fedora:
```bash
sudo yum groupinstall -y "Development Tools"
sudo yum install -y mpich-devel python3 python3-pip git
python3 -m pip install --user pyyaml numpy

Optional (for sweep plot generation):
```bash
python3 -m pip install --user matplotlib
```
```

@section petsc_sec 3. Install PETSc

Recommended from source:
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

Official PETSc install documentation:
- https://petsc.org/release/install/
- https://petsc.org/release/docs/manual/

Optimized build example:
```bash
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            --download-fblaslapack --download-metis --download-parmetis \
            --with-dmswarm=1 --with-debugging=0 \
            --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3'
make all
make check
```

Set environment variables (`~/.bashrc` or equivalent):
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
```

Reload shell and verify:
```bash
source ~/.bashrc
echo "$PETSC_DIR"
echo "$PETSC_ARCH"
```

@section clone_sec 4. Clone PICurv

```bash
git clone https://github.com/VishalKandala/PICurv.git
cd PICurv
```

@section build_sec 5. Build with `pic.flow`

```bash
./scripts/pic.flow build
```

Expected artifacts:
- `bin/picsolver`
- `bin/postprocessor`

Useful variants:
```bash
./scripts/pic.flow build clean-project
./scripts/pic.flow build SYSTEM=cluster
```

@section verify_sec 6. Quick Sanity Checks

```bash
./scripts/pic.flow --help
./scripts/pic.flow run --help
./scripts/pic.flow sweep --help
```

@section next_steps_sec 7. Next Steps

Proceed to **@subpage 02_Tutorial_Programmatic_Grid** to run your first simulation.
