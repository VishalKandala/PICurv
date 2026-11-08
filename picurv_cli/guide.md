# PICurv CLI Package Guide

This directory contains the importable implementation of the PICurv conductor.

- `main.py`: command-line entrypoint function.
- `cli.py`: argument parser construction and command dispatch.
- `core.py`: current workflow implementation and compatibility surface.
- `storage/`: the run, study, and workspace storage lifecycle (see its own guide).
- `source_expressions.py`: the source expression language's Python side - functions and
  variables, validation, constant substitution, region lowering, NumPy evaluation. `core.py`
  checks every source with it; `generators/ic.gen` evaluates one.
- `periodic_spectral.py`: Fourier machinery for uniform periodic axes, shared by `ic.gen`,
  `spectra.gen` and `reduce.gen`; also `require_numpy()`, the one lazy NumPy import (with
  its retry without other Pythons' site-packages) that `core.py`, `source_expressions.py`
  and every generator use.
- `petsc_binary.py`: the PETSc binary Vec envelope - the reader `profile.gen` and `spectra.gen`
  share, the writer `ic.gen` and `core.py` share, and the stricter validator `core.py` applies
  to every vector it stages for the solver.
- `picgrid.py`: the PICGRID text grid format - its coordinate precision, the header reader
  `core.py` uses while streaming a grid to the run, and the block reader `ic.gen`,
  `profile.gen`, `reduce.gen` and `spectra.gen` share.

The last four are plain modules, not package imports: `core.py` and the `generators/*.gen`
scripts put this directory on `sys.path` and import them by name, which also works when
`core.py` is loaded as a lone file, as the tests do.

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
