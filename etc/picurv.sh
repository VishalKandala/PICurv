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
    # Keep bin/ ahead of scripts/ for compiled executables, but also expose
    # scripts/ so `picurv` remains available if bin/picurv is temporarily absent
    # after a pull/rebase and before `make conductor` recreates the symlink.
    for _picurv_path in "$PICURV_DIR/scripts" "$PICURV_DIR/bin"; do
        case ":$PATH:" in
            *":$_picurv_path:"*) ;;  # already on PATH
            *) PATH="$_picurv_path:$PATH" ;;
        esac
    done
    export PATH
fi
unset _picurv_dir
unset _picurv_path
