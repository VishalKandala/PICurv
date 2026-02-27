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

# --- Directory Discovery Helpers ---
is_project_root() {
    local dir="$1"
    [ -f "$dir/Makefile" ] && [ -d "$dir/src" ] && [ -d "$dir/include" ]
}

find_project_root_upwards() {
    local dir="$1"
    while [ "$dir" != "/" ]; do
        if is_project_root "$dir"; then
            echo "$dir"
            return 0
        fi
        dir="$(dirname "$dir")"
    done

    if is_project_root "/"; then
        echo "/"
        return 0
    fi

    return 1
}

# --- Configuration ---
# Resolve this script directory first (works even when called via symlink).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Allow explicit override for automation/use in external case directories.
CODE_DIR="${PICURV_CODE_DIR:-}"

# If not explicitly provided, try to discover project root automatically.
if [ -z "$CODE_DIR" ]; then
    CODE_DIR="$(find_project_root_upwards "$(pwd)" || true)"
fi
if [ -z "$CODE_DIR" ]; then
    CODE_DIR="$(find_project_root_upwards "$SCRIPT_DIR" || true)"
fi

# --- Pre-flight Check ---
if [ -z "$CODE_DIR" ] || ! is_project_root "$CODE_DIR"; then
    echo "Error: Could not locate project root automatically."
    echo "Searched upward from:"
    echo "  1) current directory: $(pwd)"
    echo "  2) script directory : $SCRIPT_DIR"
    echo "Expected to find: Makefile, src/, and include/."
    echo "If needed, set PICURV_CODE_DIR=/absolute/path/to/PICurv."
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
