@page getting_started Getting Started

This guide provides step-by-step instructions for setting up your environment, installing dependencies, cloning the source code, and compiling the PICurv solver.

@section prereqs_sec 1. Prerequisites

Before you begin, ensure you have the following software installed on your system. These are essential for building and running the solver.

- **A C Compiler:** Such as GCC or Clang.
- **An MPI Implementation:** Required for parallel computing. Common choices are MPICH or OpenMPI.
- **Git:** For cloning the source code repository.
- **Make:** The standard build automation tool.
- **PETSc (version 3.20.3 or newer):** This is the core dependency. PICurv relies on PETSc for its parallel data structures (DMDA, DMSwarm), linear and nonlinear solvers, and I/O.

@section install_sec 2. Installation and Verification

This section describes how to install the required dependencies.

@subsection install_tools_sec 2.1. Installing Build Tools and MPI

On most Linux distributions, you can install the compiler, Make, and an MPI library using the system's package manager.

**For Debian/Ubuntu:**
```bash
sudo apt-get update
sudo apt-get install build-essential gfortran mpich
```

**For CentOS/RHEL/Fedora:**
```bash
sudo yum groupinstall "Development Tools"
sudo yum install mpich-devel
```

**For macOS (using Homebrew):**
```bash
brew install gcc open-mpi
```

@subsection install_petsc_sec 2.2. Installing PETSc

Installing PETSc correctly is the most critical step. We recommend building it from source to ensure all required components are enabled.

1.  **Download PETSc:**
    Go to the [PETSc website](https://petsc.org/release/download/) and download a recent version (e.g., 3.20.3). Or, use git:
    ```bash
    git clone -b v3.20.3 https://gitlab.com/petsc/petsc.git
    cd petsc
    ```

2.  **Configure PETSc:**
    The PICurv solver specifically requires `DMSwarm`. We recommend a configuration like the one below. Run this from the `petsc` directory.

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
                --with-debugging=0 --COPTFLAGS='-O3 -march=native' --CXXOPTFLAGS='-O3 -march=native' --FOPTFLAGS='-O3 -march=native'
    ```

3.  **Build and Test PETSc:**
    After configuration is complete, follow the instructions printed to the screen.
    ```bash
    make all
    make check
    ```

@subsection verify_petsc_sec 2.3. Verifying the PETSc Installation

Once PETSc is built, you must set two environment variables so that other applications (like PICurv) can find it.

1.  **Set Environment Variables:**
    Add the following lines to your shell startup file (e.g., `~/.bashrc`, `~/.zshrc`). **Replace `/path/to/your/petsc` with the actual path.**
    ```bash
    export PETSC_DIR=/path/to/your/petsc
    export PETSC_ARCH=arch-linux-c-debug  # Or whatever arch was created during configure
    export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:${PATH}"
    ```
    Now, source your startup file or open a new terminal:
    ```bash
    source ~/.bashrc
    ```

2.  **Confirm the Variables are Set:**
    Run `echo $PETSC_DIR` and `echo $PETSC_ARCH`. They should print the paths you just set.

@section get_code_sec 3. Getting the PICurv Code

With all prerequisites in place, you can now clone the PICurv repository.

```bash
git clone https://github.com/your-username/picurv.git
cd picurv
```
*(Note: Replace the URL with your actual repository URL.)*

@section build_sec 4. Compiling PICurv

The provided `Makefile` is designed to automatically find your PETSc installation using the environment variables you set.

1.  **Navigate to the Root Directory:**
    Ensure you are in the top-level `picurv/` directory.

2.  **Clean Previous Builds (Optional but Recommended):**
    It's good practice to start with a clean slate.
    ```bash
    make clean_all
    ```

3.  **Build the Solver and Postprocessor:**
    The `Makefile` has targets for the main executables.
    ```bash
    make picsolver       # Builds the main CFD/particle solver
    make postprocessor   # Builds the data postprocessor utility
    ```
    If the compilation is successful, you will find the executables `picsolver` and `postprocessor` in the `bin/` directory.

@section next_steps_sec Next Steps

Congratulations, you have successfully built the PICurv solver!

You are now ready to run your first simulation. Please proceed to the **@ref user_guide** for instructions on setting up and running a test case.

```
