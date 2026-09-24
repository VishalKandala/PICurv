# PICurv CLI Package Guide

This directory contains the importable implementation of the PICurv conductor.

- `main.py`: command-line entrypoint function.
- `cli.py`: argument parser construction and command dispatch.
- `core.py`: current workflow implementation and compatibility surface.
- `storage/`: the run, study, and workspace storage lifecycle (see its own guide).

`picurv_cli/picurv` is the stable source-tree executable and `bin/picurv` selects
the configured Python runtime before invoking it. Further internal extraction
from `core.py` can proceed without changing those public entrypoints.

Physical-time ParaView collections belong to the existing post workflow in
`core.py`: `normalize_post_paraview_series_config` owns policy,
`resolve_run_lineage` reads manifest ancestry, and `finalize_post_paraview_series`
indexes existing VTK output with checkpoint times. The post lock wrapper runs the
serial finalizer after successful field processing; caught-up invocations can
refresh the index without relaunching the field executable. No separate command
or lineage metadata file is required. See page 10 for the user contract and page 52
for artifact lifecycle behavior.
