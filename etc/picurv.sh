#!/bin/bash
# PICurv environment setup.
# Source this file to add picurv to your PATH.
#
#   source /path/to/PICurv/etc/picurv.sh
#
# Or add the line above to your ~/.bashrc for persistent access.

_picurv_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)"
if [ -n "$_picurv_dir" ]; then
    export PICURV_DIR="$_picurv_dir"
    case ":$PATH:" in
        *":$PICURV_DIR/bin:"*) ;;  # already on PATH
        *) export PATH="$PICURV_DIR/bin:$PATH" ;;
    esac
fi
unset _picurv_dir
