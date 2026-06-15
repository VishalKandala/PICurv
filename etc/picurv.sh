#!/bin/bash
# PICurv environment setup.
# Source this file to add PICurv executables to your PATH.
#
#   source /path/to/PICurv/etc/picurv.sh
#
# Or add the line above to your ~/.bashrc for persistent access.

_picurv_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)"
if [ -n "$_picurv_dir" ]; then
    export PICURV_DIR="$_picurv_dir"
    if [ -x "$PICURV_DIR/.picurv-venv/bin/python" ]; then
        export PICURV_PYTHON="$PICURV_DIR/.picurv-venv/bin/python"
        export PICURV_PYTHON_IS_MANAGED=1
    elif [ -f "$PICURV_DIR/.picurv-python" ]; then
        _picurv_python="$(sed -n '1p' "$PICURV_DIR/.picurv-python" 2>/dev/null)"
        if [ -n "$_picurv_python" ] && [ -x "$_picurv_python" ]; then
            export PICURV_PYTHON="$_picurv_python"
            case "$_picurv_python" in
                */.picurv-venv/bin/python) export PICURV_PYTHON_IS_MANAGED=1 ;;
            esac
        fi
    fi
    # Keep bin/ ahead of picurv_cli/ for compiled executables, but also expose
    # picurv_cli/ so `picurv` remains available if bin/picurv is temporarily absent
    # after a pull/rebase and before `make conductor` recreates the launcher.
    for _picurv_path in "$PICURV_DIR/picurv_cli" "$PICURV_DIR/bin"; do
        case ":$PATH:" in
            *":$_picurv_path:"*) ;;  # already on PATH
            *) PATH="$_picurv_path:$PATH" ;;
        esac
    done
    export PATH
fi
unset _picurv_dir
unset _picurv_path
unset _picurv_python
