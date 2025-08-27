/** @mainpage 🌀 PICurv: A Parallel Particle-in-Cell Solver for Curvilinear LES
@brief PICurv is a parallel particle-in-cell (PIC) solver designed for low-Mach-number turbulent flows in complex geometries. Built using PETSc’s DMDA and DMSwarm infrastructure, it features:
- A fractional-step incompressible Navier–Stokes solver with dynamic Smagorinsky LES
- Support for body-fitted curvilinear grids and immersed boundaries (CURVIB)
- Fully MPI-parallel particle tracking with two-way Eulerian–Lagrangian coupling
- Trilinear interpolation (grid → particle) and conservative face-based projection (particle → grid)

The solver has been validated on canonical geometries such as curved pipes and demonstrates scalable performance for millions of particles on both structured Cartesian and curvilinear meshes.

@section keyfeatures 🔧 Key Features
Parallelized 3D particle-in-cell solver with PETSc
Low-Mach number incompressible Navier–Stokes with LES modeling
Two-way coupled Lagrangian particle transport
Trilinear interpolation and deposition kernels
Geometric flexibility via immersed boundary methods (CURVIB)
Modular architecture with plug-and-play routines for scalar transport, stochastic mixing, and particle diagnostics

@section gettingstarted 🚀 Getting Started 
- Dependencies: PETSc 3.20.3 or newer (Built with MPI and 'DMSWARM' support).
- Build:
To build the solver and postprocessor:

```bash
cd   $pic          # change to the pic directory.
make clean_all     # clean all object and executable files.
make picsolver       # Builds main solver executable `picsolver`
make postprocessor   # Builds postprocessor executable `postprocessor`
```

Ensure PETSc paths are configured correctly in your Makefile or environment.

@section directorystructure Directory Structure
```text
.
├── src/           # Source files (*.c)
├── include/       # Header files (*.h)
├── scripts/       # Python utilities
├── obj/           # Object files (*.o)
├── bin/           # Executables
├── docs/          # Docs files
├── test/          # Config directories for different test cases.
```

@section running ▶️ Running The Code

1.  **Create a new `run/` directory and `logs/`,`results/` sub-directories.**.
2.  **Place the `test/config` directory in `run/`
3.  **Create symbolic links** to the built executables:
    ```bash
    ln -s $pic/bin/picsolver 
    ln -s $pic/bin/postprocessor 
    ```
4.  Run the solver:
    ```bash
    mpirun -np 4 ./picsolver
    ```

*/
