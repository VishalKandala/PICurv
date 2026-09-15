- A queued job can no longer silently run a different build than it was staged with, and
  the commands that change the installation no longer report success when they changed
  nothing.
  - Every generated Slurm script - solve, post, sweep array, and sweep continuation - runs
    the executable's `--version` after module setup and before the launcher, and exits
    before any rank starts when the identity differs from the one read at staging.
    Previously a job resolved `bin/simulator` when it started, so a rebuild while it sat in
    the queue changed what ran.
  - `picurv run` and `picurv sweep` accept `--pin-executables` and `--no-pin-executables`,
    and `reproducibility.pin_executables: true` in `.picurv-workspace.yml` pins every run
    staged in that workspace. A pinned run carries its own `simulator` and `postprocessor`
    in `<run.config.bin>`, recorded in `active.json`; every later stage of the run launches
    those copies, and the manifest and software lock describe them. A continuation keeps
    the pin, `--pin-executables` re-pins it into the new configuration revision, and
    `--no-pin-executables` on a pinned run is refused. Sweep members each receive a copy
    launched through `$RUN_DIR`. Pinning is opt-in; without it the job-start check is the
    protection.
  - `picurv versions install` and `versions activate` pass make variables and options
    through the same path as `picurv build`, so
    `picurv versions install 0.1.0 SYSTEM=cluster` builds the cluster configuration
    instead of depending on `SYSTEM` being exported. A bare release resolves to its
    `v`-prefixed tag, which `versions activate` from a workspace pin needs. Success is
    reported only after both executables carry the identity of the commit just checked
    out.
  - `picurv pull-source` refuses a checkout detached by `versions install` or
    `versions activate`. Its multi-branch path used to pull the branches, restore the detached
    commit, and print success with the running code unchanged.
  - A binary whose `--version` fails, typically because its shared libraries are not
    loadable in the staging shell, is recorded with the exit status and first error line
    instead of "no build identity reported", and staging now warns when an identity cannot
    be read. A run manifest had recorded the latter for a binary that does report one.
  - `picurv init --pin-binaries` copies are used only when that case directory's own
    `picurv` is invoked, which `init` never sets up, so the `picurv` on `PATH` ignored
    them. The pages that recommended it for Slurm jobs now describe run-time pinning; the
    flag remains until the next release.
