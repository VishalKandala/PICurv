- The example directories and master templates were reviewed for stale
  configuration, and the channel and duct campaign cases were prepared for the
  cluster verification campaign. `plane_spectrum` and `line_spectrum` are now
  supported, on a measurement that matched a direct DFT of a known checkpoint to
  round-off. Several fixes change results or configuration.
  - Every shipped channel and duct grid, and the driven-periodic smoke fixture,
    used an odd cell count, so its multigrid coarsening was misaligned. They now use
    16, 32, 64, 128 or 256 cells and declare `mg_levels`, which makes `grid.gen`
    refuse a count that cannot coarsen.
  - The `msd_final` study metric read the last column of the MSD table
    (`frac_3sigma_pct`) instead of `MSD_total`. Explicit `log_regex` metrics now pass
    validation and must capture their value in a group.
  - The wall spectral seeds (`channel_spectral_velocity`, `duct_spectral_velocity`)
    now grow linearly from the wall instead of as the fourth power, which puts
    perturbation energy into the buffer layer. The same seed gives a different field.
  - Generated initial conditions and inlet profiles refuse `output_file`,
    `summary_json` and `spectrum_csv` instead of ignoring them. `validate` warns that
    `source_data.directory` is ignored instead of checking that it exists; omitting
    `source_data` keeps the recipe ID.
  - Spectra on a `programmatic_c` case are refused at validation instead of failing
    after the solve, and `picurv init` no longer drops a post recipe that has no
    pipeline.
  - The restart guides no longer tell users to set `eulerian_field_source: load` to
    continue a run; that mode replays stored fields and never solves.
  - New `les_wallmodel_retau1000` driven-channel case and `humphrey_bend_wm.cfg` bend
    grid; the driven channel and duct campaign cases seed with the wall spectral
    providers.
