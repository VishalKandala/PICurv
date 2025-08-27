/**
@page user_guide User Guide

This section provides detailed information on how to configure and run simulations.

@section cli_options Command-Line Options
The solver is configured at runtime using command-line flags. The `CreateSimulationContext()` function parses these options. Key options include:
- `-np <count>`: Number of Lagrangian particles to initialize.
- `-steps_to_run <count>`: Number of time steps to execute.
- `-start_step <count>`: The step number to start or restart from.
- `-only_setup`: If present, the code will initialize everything and then exit without running the time loop. Useful for debugging setup issues.
- ... (Add other important flags here)

@section test_cases Example Test Cases
Detailed instructions for running cases in the `test/` directory.

@section output_files Output and Checkpointing
Description of the files generated in the `results/` directory.

*/
