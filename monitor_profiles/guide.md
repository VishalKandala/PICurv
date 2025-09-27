# PIC-Flow Monitor Profiles Library

## 1. Purpose of this Directory

This directory contains a library of reusable **Monitor Profiles** (`.yml` files). Each file defines a strategy for **observing and recording** a simulation run.

Changing a monitor profile **will never change the scientific result** of the simulation. It only changes what you see on the console, what performance data is collected, and how often results are saved to disk. This allows you to run the exact same case with different monitoring profiles for different purposes, such as production runs, debugging, or performance analysis.

A monitor profile controls:
- Console and log file verbosity (`logging`).
- Performance measurement and reporting (`profiling`).
- The frequency and location of data output (`io`).
- Detailed debugging of the numerical solvers (`solver_monitoring`).

## 2. Included Profiles

- **`Standard_Output.yml`**: The recommended profile for everyday runs. It provides high-level progress updates (`INFO` verbosity) and saves data at a reasonable interval.
- **`master_monitor_profile.yml`**: A heavily commented reference file showcasing **all possible** monitoring options. Use this file as a guide for creating your own custom profiles.

*(You should create other profiles as needed, such as `Debug.yml`, `Performance.yml`, `Frequent_Output.yml`.)*

## 3. How to Use

Specify which monitor profile to use with the `--monitor` flag:

```bash
# Example: Running a case with the standard monitor profile
./bin/pic.flow run \
    --case my_studies/my_case.yml \
    --solver ... \
    --monitor monitor_profiles/Standard_Output.yml
```

To debug the same case, you would simply create a `Debug.yml` profile and change the flag:
`--monitor monitor_profiles/Debug.yml`

## 4. Creating a New Monitor Profile

Create a new profile to suit a specific task.

**Example: Creating `Debug-Grid.yml` for targeted debugging:**
1.  **Copy `Standard_Output.yml`** to `Debug-Grid.yml`.
2.  **Edit `Debug-Grid.yml`:**
    -   Change `logging.verbosity` to `DEBUG`.
    -   Add the specific C function names you want to inspect to the `logging.enabled_functions` list (e.g., `SetupGridAndSolvers`).
    -   Decrease `io.data_output_frequency` to `10` to get more frequent snapshots.
3.  **Save and use:** Your new targeted debugging profile is ready.