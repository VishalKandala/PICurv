@page 67_Troubleshooting Troubleshooting by Symptom

@anchor _Troubleshooting

@pagemeta{How-to, Anyone with a failing run, Current behavior}

Start from what you observed. Each symptom below runs cheapest checks first, so you
spend a minute before you spend a cluster allocation.

This page complements **@subpage 39_Common_Fatal_Errors**, which is organized by
error message. Use that page when you have an exact string to search for; use this
one when you only have a behavior.

@tableofcontents

@section p67_validate_sec 1. Validation Fails Before Anything Runs

Good news: this is the cheapest failure mode, and it is the one the system is
designed to produce.

1. Read the reported key path and file — validation names both.
2. Check the value against the reference page for that role: **@subpage 07_Case_Reference**,
   **@subpage 08_Solver_Reference**, **@subpage 09_Monitor_Reference**,
   **@subpage 10_Post_Processing_Reference**.
3. If it is a selector value, check it against the generated inventory for that
   family — an unknown value is rejected rather than defaulted.
4. If it is a boundary condition, check the type/handler pairing at
   @ref p44_types_sec; not every handler is legal for every type.
5. Solver-specific blocks must match the selected solver. Setting a
   `dual_time_picard_jameson_rk` block while selecting `Newton Krylov` is an error,
   not an ignored stanza.

@section p67_noconverge_sec 2. The Solver Does Not Converge

Establish *where* it fails before changing anything.

1. **Is the Poisson side healthy?** Check the reported maximum divergence. If it is
   near machine zero, the pressure solve is fine and the problem is in momentum.
2. **Is it the momentum solve?** Under `Dual Time Picard Jameson RK`, look for
   "reached N total attempts without convergence". That counts total attempts,
   accepted plus rejected.
3. **Is `dtau` pinned at its floor?** A pseudo-CFL collapsed to `pseudo_cfl.minimum`
   with a residual ratio plateauing near 1 is the signature of a stalled controller
   rather than a hard case.
4. **Reduce the physical timestep** before touching solver internals. Many apparent
   solver failures are a timestep the physics will not support.
5. **Try a different solver.** @ref p08_cap_newton_krylov gives true Newton
   convergence where pseudo-time stalls.

@warning If this is a periodic wall-bounded case, read @ref p54_driven_limits_sub
first. A stall was recorded there on 2026-08-24 and **requires re-characterization at
current `HEAD`**; you may be reproducing a known open issue rather than a
configuration mistake.

@section p67_diverge_sec 3. The Solution Blows Up

1. Check the initial condition is physically sensible — visualize step 0 before
   blaming the solver.
2. Check inlet and outlet are consistent: an inlet with no outlet to conserve
   against will not behave.
3. For `Explicit RK4`, divergence within a few steps is the explicit stability limit;
   there is no recovery mechanism on that path. Reduce the timestep or change solver.
4. Check the boundary conditions actually applied — the startup banner lists the
   resolved handler per face. A face you thought was a wall may not be.

@section p67_particles_sec 4. Particles Are Lost or Behave Oddly

1. Check the reported lost-particle count per step. Steady loss usually means
   particles are leaving through a boundary you did not intend.
2. A freshly seeded particle has valid state: setup interpolates once at `t=0` before
   the loop. If early behavior looks wrong, suspect the seeding, not the ordering.
3. Isolate the mechanism using the verification family in
   **@subpage 65_Example_Catalog** — each isolates one term of the displacement
   equation, so a failure points at one thing.
4. For search or migration anomalies, **@subpage 53_Search_Robustness_Metrics_Reference**
   defines the metrics to look at.

@section p67_output_sec 5. Output Is Missing or Unreadable

1. Check the monitor configuration: output cadence controls whether anything is
   written at all.
2. Check the postprocessor's own stream log, written to **`<run.scheduler>/`** as
   `<run_id>_<prefix>.log` — not `<run.runtime_logs>/`, which carries the solver's runtime
   logging. Solver success does not imply postprocessor success.
3. Confirm the field you expect is one the run actually produced;
   **@subpage 56_Field_Identity_and_Layout_Catalog** is authoritative on what exists.
4. Remember `Ucat` is *derived* from `Ucont`. A visualization field is not the
   solver's evolved state — see **@subpage 20_Grid_Cell_Architecture_Guide** if a field looks inconsistent.

@section p67_restart_sec 6. Restart Is Rejected or Behaves Unexpectedly

1. Confirm you want the operation you asked for: `--continue` resumes the same run;
   `--restart-from` seeds a new one. See **@subpage 52_Run_Lifecycle_Guide**.
2. Check that what you changed is compatible with continuing — @ref p52_compat_sec.
3. Restart is **not** bit-exact. A small difference immediately after restart is
   expected and documented in **@subpage 29_Maintenance_Backlog**; it is a structural
   floor, not a tolerance you can tighten away.
4. For `initial_flux`, a restart restores the latched target from the checkpoint. A
   log line saying it re-measured instead means the checkpoint predates that metadata.

@section p67_cluster_sec 7. The Job Fails on a Cluster but Works Locally

1. Distinguish scheduler failure from solver failure: look at the scheduler `.err`
   file before the solver log.
2. Confirm binaries and paths resolve on the compute node, not just the login node.
3. Check the MPI rank count against the grid decomposition. A rank count that
   partitions the grid awkwardly can also limit how deeply the Poisson multigrid can
   coarsen, since each level must remain decomposable across the ranks it runs on.
   If the Poisson solve refuses your `mglevels` at one rank count but not another,
   this is usually why.
4. Stage without submitting and inspect the generated scripts before blaming the
   solver: **@subpage 36_Cluster_Run_Guide**.

@section p67_suspicious_sec 8. Results Look Numerically Suspicious

This is the hardest class, and the one where documentation helps least — but three
checks are worth doing before anything else:

1. **Check units.** Values are physical in YAML and non-dimensionalized before C.
   `target_flux` is a volumetric flux, not a velocity. A constant-factor error here
   produces plausible-looking wrong answers. See @ref p44_nondim_sec.
2. **Check the evidence for what you used.** **@subpage 66_Evidence_Matrix** shows
   what has actually been verified. An unverified capability is a candidate suspect.
3. **Check for known limitations** on the capability entry for each selector you set.

@section p67_les_sec 9. LES Results Look Wrong

The closure is quiet when it misbehaves: a wrong coefficient still produces a run that
completes and plots that look plausible. Turn the diagnostics on
(`les.diagnostics.enabled: true`) and read `<run.runtime_logs>/les_coefficient.csv`
before anything else.

**`cs_effective` is zero for the whole run.** Expected for the first two steps of a run
started from rest, because the procedure has no developed field to sample. Beyond that,
check that the model is `dynamic_smagorinsky`; the constant model reports the
coefficient you configured, not a measured one.

**`cs_effective` is far from 0.16-0.17 in isotropic turbulence.** Check the grid first,
not the model. That value assumes the grid cutoff sits in an inertial range; at low
`Re_lambda` a smaller coefficient is correct. Then check `limited_fraction` — if the
clip is binding, the ceiling is setting the answer.

**`cs_effective` is noisy or oscillating.** Almost always `averaging.mode: local`,
where the coefficient's denominator collapses wherever the resolved strain is briefly
small. If the flow has a homogeneous direction, use `averaging.mode: homogeneous`.

**`limited_fraction` is large.** The `max_cs` ceiling is doing the modelling. Either
the coefficient is genuinely diverging, which is a resolution or stability problem, or
the ceiling is set too low for the flow.

**`backscatter_fraction` is large and the run is over-dissipative.** The clipping modes
discard only the negative tail, which biases mean dissipation upward. Averaging over a
homogeneous set lets forward and reverse transfer cancel the way the physics intends;
`clipping.mode: none` keeps the sign and bounds the total viscosity instead.

**`nu_t_over_nu_mean` is far below one.** The model is contributing nothing. Check that
LES is actually enabled, and that the grid is coarse enough for a subgrid model to be
the right tool.

**Validation rejected `simpson_ik`.** That stencil averages over the central eta-plane
only and assumes xi and zeta are homogeneous. The check requires both declared
`PERIODIC`. Use `volume_weighted_box` on any other geometry.

More detail on every setting at **@subpage 07_Case_Reference**, and the reasoning at
**@subpage 72_LES_Turbulence_Closure**.

@section p67_related_sec 10. Related Documentation

- **@subpage 39_Common_Fatal_Errors** — organized by error message
- **@subpage 66_Evidence_Matrix** — what is verified and what is not
- **@subpage 40_Testing_and_Quality_Guide** — reproducing with the test suites
