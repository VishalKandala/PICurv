#!/bin/bash

# =============================================================================
#
#                    build.sh - Smart Build Wrapper
#
# =============================================================================
#
# ## Description
#
#   This script is a user-friendly front-end for the project's Makefile,
#   designed to be run from the `run/` directory.
#
#   It intelligently works around a known issue in the Makefile where the
#   default `make` command does not trigger a build. This script ensures that
#   running it without arguments correctly builds the project by explicitly
#   calling the `all` target.
#
# ## Behavior
#
#   - If run with NO arguments (`./build.sh`):
#     It executes `make all` in the code directory.
#
#   - If run WITH arguments (`./build.sh clean-project`):
#     It passes all arguments directly to make, executing `make clean-project`.
#
# =============================================================================

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
# Relative path to the directory containing the Makefile.
CODE_DIR="/root/PE/PICurv"

# --- Pre-flight Check ---
# Ensure the code directory actually exists.
if [ ! -d "$CODE_DIR" ]; then
    echo "Error: Code directory not found at '$CODE_DIR'."
    echo "Please ensure this script is in the 'run' directory and the 'code' directory is a sibling."
    exit 1
fi

# --- Main Logic ---
if [ $# -eq 0 ]; then
    # Case 1: No arguments were provided by the user.
    echo "--- No target specified. Defaulting to 'make all' for a full build. ---"
    make -C "$CODE_DIR" all
else
    # Case 2: The user provided arguments (e.g., 'clean-project').
    echo "--- Executing make in '$CODE_DIR' with specified arguments: '$*' ---"
    # "$@" expands to all the arguments the user provided.
    make -C "$CODE_DIR" "$@"
fi

echo "--- Make command finished successfully. ---"
