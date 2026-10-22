- Experimental `reducers` post stage (`post.yml -> reducers`, `generators/reduce.gen`): `wall_profile`,
  `cross_section` and `profile_series` run after the field stage with their inputs resolved from the
  run, and `--only reducers` re-runs them alone. The driven channel and duct example tools are
  replaced by these tasks; `summarize --plot profile.*` now reads the `profile_series` histories.
  Profile CSVs print 17 significant digits.
