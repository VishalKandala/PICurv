- Added `io.paraview_series` to `post.pipeline`: PICurv writes a `.pvd` ParaView
  collection indexing existing VTK frames by checkpoint physical time after
  successful field postprocessing. `scope: lineage` follows the run manifest's
  restart ancestry, clips each parent at its child's fork step, and lets the
  child own a duplicate fork frame; `scope: run` indexes only the current run.
  Indexing is a presentation product over already-produced VTK output and does
  not change computational recipe identity. `picurv run --post-process
  --continue` can be repeated while a solver is active, appending newly
  available frames and atomically refreshing the collection each time.
