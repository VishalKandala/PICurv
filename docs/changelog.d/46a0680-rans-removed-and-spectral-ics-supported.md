- The four generated initial-condition providers are now supported. `ic_gen` and
  `spectral_random_velocity` cite the measurement that already covered them, and the two
  wall-bounded providers were measured against their contract
  (`wall-spectral-ic-2026-09-22`): on stretched channel and duct grids the staged field
  carries the requested bulk velocity and perturbation RMS exactly, its mean equals the
  flux-normalized product parabola to 3e-15, the walls carry zero velocity, and the flux
  field the runtime reconstructs is divergence-free at round-off. They remain startup
  constructions: nothing here establishes sustained turbulence.
- **RANS is removed and returned to planned.** `models.physics.turbulence.rans` is
  refused at validation and `-rans` at setup; both name the status rather than ignoring
  the setting. Nothing behind the selector was ever implemented - the `K_Omega` fields
  were never allocated, the transport update in `FlowSolver` was commented out, and the
  function it called was defined nowhere - so a case that enabled it aborted on a null
  vector at the end of the first step. It was recorded known-defective on 2026-09-18;
  since no implementation had ever existed, removing the dead hooks returns it to
  `planned`, with the charter for a future closure on page 57 and the removed hooks
  listed in `src/guide.md`.
- Removed with it: the `K_Omega` fields and their catalog entries, the
  `FIELD_AVAILABILITY_RANS` flag, the two-component DM they alone used (never
  created), the `simCtx->rans` branches in the RHS, run loop, checkpoint I/O and
  Newton-Krylov scope check, and the `k_omega` capability value.
- **Checkpoint metadata:** `-checkpoint_rans` is no longer written or required. The
  format version is unchanged, so existing checkpoints still load; a checkpoint
  written now is not readable by a pre-removal binary, which requires that key.
- Wall functions are unaffected, being configured independently of any closure. Only
  the two pairing refusals that named RANS (`cabot` and `werner` under it) went with
  the selector.
- A subsystem may now move from `known-defective` to `planned`, with a recorded reason
  saying what was found and what was removed: a feature that never existed is not a
  feature that broke.
