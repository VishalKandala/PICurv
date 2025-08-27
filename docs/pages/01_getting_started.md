/**
@page getting_started Getting Started

This guide provides instructions for setting up your environment, building the solver, and running your first simulation.

@section dependencies Dependencies
- PETSc 3.20.3 or newer...

@section building Building the Code
To build the solver and postprocessor...
`make picsolver`

@section directory_structure Directory Structure
The project is organized as follows...
`
.
├── src/
...
`

@section first_run Running a Simulation
Create a new `run/` directory...
`mpirun -np 4 ./picsolver`

*/
