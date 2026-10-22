- Reducer friction velocities from `driven_flow` and `wall_model` average the metrics of the whole restart
  lineage, so a window carried into a branched run no longer fails. A post job whose postprocessor
  exits cleanly but leaves planned outputs unwritten, as a pinned older build does, now warns.
