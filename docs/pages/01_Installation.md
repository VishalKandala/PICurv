@page 01_Installation Installation Guide

This guide provides step-by-step instructions for setting up your environment, installing dependencies, cloning the source code, and compiling the PICurv solver platform.

@tableofcontents

@section prereqs_sec 1. Prerequisites

Before you begin, ensure you have the following software installed on your system. These are essential for building and running the solver.

- **A C Compiler:** Such as GCC or Clang.
- **An MPI Implementation:** Required for parallel computing. Common choices are MPICH or OpenMPI.
- **Git:** For cloning the source code repository.
- **Python 3:** For running the `pic.flow` conductor script.
- **GNU Make:** The standard build automation tool.
- **PETSc (version 3.20.3 or newer):** This is the core dependency. PICurv relies on PETSc for its parallel data structures (DMDA, DMSwarm), linear solvers, and I/O.

@section install_sec 2. Installing Dependencies

This section describes how to install the required dependencies. Correctly installing PETSc is the most critical step.

@subsection install_tools_sec 2.1. Installing Build Tools, MPI, and Python

On most Linux distributions, you can install the essential tools using the system's package manager.

**For Debian/Ubuntu:**
```bash
sudo apt-get update
sudo apt-get install build-essential gfortran mpich python3 python3-pip git
pip3 install pyyaml numpy
```

**For CentOS/RHEL/Fedora:**
```bash
sudo yum groupinstall "Development Tools"
sudo yum install mpich-devel python3 python3-pip git
pip3 install pyyaml numpy
```

@subsection install_petsc_sec 2.2. Installing PETSc

We strongly recommend building PETSc from source to ensure all components required by PICurv are enabled.

1.  **Download PETSc:**
    Go to the [PETSc website](https://petsc.org/release/download/) and download a recent version (e.g., 3.20.3), or use git:
    ```bash
    git clone -b v3.20.3 https://gitlab.com/petsc/petsc.git
    cd petsc
    ```

2.  **Configure PETSc:**
    PICurv specifically requires `DMSwarm`. We recommend a configuration like the one below. Run this from the `petsc` source directory.

    *For a debugging build (recommended for development):*
    ```bash
    ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
                --download-fblaslapack --download-metis --download-parmetis \
                --with-dmswarm=1 \
                --with-debugging=1
    ```

    *For an optimized build (recommended for production runs):*
    ```bash
    ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
                --download-fblaslapack --download-metis --download-parmetis \
                --with-dmswarm=1 \
                --with-debugging=0 --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3'
    ```

3.  **Build and Test PETSc:**
    After configuration is complete, follow the instructions printed to the screen. This will look something like this:
    ```bash
    make all
    make check
    ```

@subsection verify_petsc_sec 2.3. Setting PETSc Environment Variables

Once PETSc is built, you **must** set two environment variables so that other applications (like PICurv) can find it.

1.  **Set Environment Variables:**
    Add the following lines to your shell startup file (e.g., `~/.bashrc`, `~/.zshrc`). **Replace `/path/to/your/petsc` with the actual path.** The `PETSC_ARCH` name will be printed on screen after you run `./configure`.

    ```bash
    # Example for a debug build on Linux
    export PETSC_DIR=/path/to/your/petsc
    export PETSC_ARCH=arch-linux-c-debug
    ```

2.  **Activate the Variables:**
    Source your startup file or open a new terminal:
    ```bash
    source ~/.bashrc
    ```

3.  **Confirm the Variables are Set:**
    Run `echo $PETSC_DIR` and `echo $PETSC_ARCH`. They should print the paths you just set. This is a critical step; PICurv's `Makefile` will fail if these are not set correctly.

@section get_code_sec 3. Getting the PICurv Code

With all prerequisites in place, you can now clone the PICurv repository.

```bash
git clone https://github.com/VishalKandala/PICurv.git
cd PICurv
```

@section build_sec 4. Compiling PICurv

The PICurv platform is managed by the `pic.flow` conductor script, which provides a simple interface to the underlying `Makefile`.

1.  **Navigate to the Root Directory:**
    Ensure you are in the top-level `PICurv/` directory.

2.  **Build the Solver and Postprocessor:**
    Use the `build` command of the `pic.flow` script.

    ```bash
    ./bin/pic-flow build
    ```
    This command will invoke the `Makefile` and compile all C source code. If the compilation is successful, you will find the executables `picsolver` and `postprocessor` in the `bin/` directory.

@subsection advanced_build_sec 4.1. Advanced Build Options

The `build` command can pass arguments directly to `make`. This is useful for cleaning the project or building for different systems.

- **Clean all build artifacts:**
  ```bash
  ./bin/pic-flow build clean-project
  ```

- **Build for a specific system (e.g., a cluster environment defined in `config.cluster.mk`):**
  ```bash
  ./bin/pic-flow build SYSTEM=cluster
  ```

@section next_steps_sec 5. Next Steps

Congratulations, you have successfully built the PICurv solver platform!

You are now ready to run your first simulation. Please proceed to the next tutorial: **@subpage 02_Tutorial_Programmatic_Grid**.
