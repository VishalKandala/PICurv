- Driven periodic flows (`constant_flux` and `initial_flux`) now write
  `<run.analysis.metrics>/driven_flow.csv` each step: target and measured flux, bulk
  velocity, the controller correction, and the driving force per unit mass actually
  applied. Its time average gives the mean wall stress (u_tau^2 / h in a channel), so a
  driven run's friction velocity can be read back directly; `picurv summarize --plot
  driven_flow.driving_acceleration` draws it, and the driven-channel and driven-duct
  profile tools take it with `--driven-flow-csv`.
