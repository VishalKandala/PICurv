
- New initial-condition generator `resampled_velocity` starts a periodic-box case from an
  existing velocity field: a published DNS, an experiment, another code's snapshot, or a field at
  another resolution. It reads raw binary, `.npy`, `.npz` or HDF5 (with `h5py` installed) through
  a declared layout, moves the field to the case's grid in Fourier space with an optional filter,
  and makes it divergence-free for PICurv. By default it stages face fluxes, which keep more of the
  field's energy than cell-centre velocities: 87% against 78% for the AGARD HOM02 DNS field on a
  64-cubed grid. Page 33 explains when to use it rather than `mode: file`.
- `mode: file` now refuses, before the solver starts, a field sized for a different grid, and
  names `resampled_velocity` as the way to bring such a field in.
- `picurv run` and `picurv sweep` report a failure to build a case's inputs as a `[FATAL]`
  message instead of a Python traceback.
