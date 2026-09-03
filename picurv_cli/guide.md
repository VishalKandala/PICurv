# PICurv CLI Package Guide

This directory contains the importable implementation of the PICurv conductor.

- `main.py`: command-line entrypoint function.
- `cli.py`: argument parser construction and command dispatch.
- `core.py`: current workflow implementation and compatibility surface.
- `storage/`: the run, study, and workspace storage lifecycle (see its own guide).

`picurv_cli/picurv` is the stable source-tree executable and `bin/picurv` selects
the configured Python runtime before invoking it. Further internal extraction
from `core.py` can proceed without changing those public entrypoints.
