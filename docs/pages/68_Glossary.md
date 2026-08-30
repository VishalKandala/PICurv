@page 68_Glossary Glossary

@anchor _Glossary

@pagemeta{Reference, All readers, Current behavior}

PICurv-specific terms, and the CFD terms this documentation uses in a particular way.
Deprecated spellings are listed with their canonical replacements so a search for the
old name lands here.

@tableofcontents

@section p68_fields_sec 1. Fields and Storage

**`Ucont`** — the contravariant volume fluxes stored on cell **faces**. This is the
quantity the solver evolves. Curvilinear codes evolve fluxes rather than Cartesian
velocities because the flux formulation keeps mass conservation exact on a distorted
grid.

**`Ucat`** — the Cartesian velocity at cell **centres**, reconstructed from `Ucont`.
Derived, not evolved. A visualization file showing `Ucat` is showing a reconstruction
of the solver's state, not the state itself.

**`Ucat_nodal`** — `Ucat` interpolated to grid nodes for visualization.

**`Aj`** — the inverse Jacobian of the curvilinear mapping, cell-centred.

**`Csi` / `Eta` / `Zet`** — face area vectors for the three computational directions
\f$\xi, \eta, \zeta\f$.

**`Nvert`** — per-cell classification marking blanked, immersed, or boundary-adjacent
cells.

**Shifted index architecture** — cell-centred data for geometric cell `i` is stored at
array index `i+1`. This is deliberate: it makes boundary handling symmetric. See
**@subpage 20_Grid_Cell_Architecture_Guide**.

@section p68_solver_sec 2. Solver Terms

**Pseudo-time / dual time** — an inner iteration driving each physical timestep to
convergence. `Dual Time Picard Jameson RK` iterates in pseudo-time; `Explicit RK4`
does not.

**Pseudo-CFL** — the dimensionless Courant number governing the pseudo-time step,
\f$\Delta\tau = \text{pseudo\_cfl} / \lambda_{max}\f$. Exclusive to the Picard-Jameson
solver. Not checkpointed: a restart begins from the configured initial value.

**Accepted versus attempted trial** — the pseudo-time controller may reject a trial
and roll back. `max_iterations` caps **accepted** iterations; a separate hard cap of
`3 x max_iterations` bounds total attempts.

**Fractional step** — the split that makes each timestep tractable: solve momentum
for an intermediate velocity, then project it onto a divergence-free field with a
pressure solve. See **@subpage 23_Fractional_Step_Method**.

**Projection** — the step that enforces incompressibility. The momentum solve alone
does not produce a divergence-free field.

@section p68_bc_sec 3. Boundary and Periodicity Terms

**Type versus handler** — `type` says what kind of boundary a face is (`wall`,
`inlet`, …); `handler` says how it is computed (`noslip`, `parabolic`, …). Not every
handler is legal for every type. See @ref p44_types_sec.

**Driven periodic** — a periodic pair plus a body force sustaining the flow against
wall friction. Without it, periodic wall-bounded flow decays.

**`target_flux`** — a **volumetric** flux (bulk velocity times cross-sectional area),
given in physical m³/s in YAML. Not a velocity. `picurv` divides it by
\f$U_{ref} L_{ref}^2\f$ before the runtime sees it.

**Seam plane** — the single plane where a periodic pair joins. `enforce_seam_flux`
adds a local correction there, on top of the global body force.

**Latent option** — present in the C enum but not exposed through validation. Not
selectable; listed so its absence is explicit rather than mysterious.

@section p68_particle_sec 4. Particle Terms

**Swarm** — the PETSc `DMSwarm` holding the Lagrangian particles.

**Settling** — locating each particle's host cell and migrating it to the owning MPI
rank after advection.

**Walking search** — the cell-location algorithm that walks from a particle's previous
host cell toward its new one, rather than searching globally.

**Scatter** — writing particle-carried quantities back onto the Eulerian grid. In the
current implementation the scatter path targets the particle-derived scalar `Psi`; it
does **not** feed particle momentum back into `Ucat` or `Ucont`. "Two-way coupling" is
therefore accurate for the scalar field and not established for momentum — prefer the
specific statement over the general phrase.

**Advect-then-interpolate** — the ordering of the particle step: a particle moves
using the velocity interpolated at the *previous* step, and interpolation closes the
step rather than opening it. See **@subpage 34_Particle_Model_Overview**.

@section p68_docs_sec 5. Documentation Terms

**Capability entry (Tier 2)** — the standard eight-part description of one selectable
value. Defined in **@subpage 64_Documentation_Extension_Framework**.

**Evidence facet** — an independent claim about what has been verified. Facets do not
form a ladder. See **@subpage 62_Capability_Status_Vocabulary**.

**Accepted spelling versus deprecated alias** — a spelling is a synonym with no
migration story (`smagorinsky` for `constant_smagorinsky`); a deprecated alias is a
retired name that should be migrated away from (`Dual Time Picard RK4`).

@section p68_deprecated_sec 6. Deprecated and Renamed Terms

| You may see | Current name | Notes |
|---|---|---|
| `Dual Time Picard RK4` | `Dual Time Picard Jameson RK` | Deprecated alias; still parses |
| `rk4_residual_noise_allowance_factor` | `jameson_residual_noise_allowance_factor` | Rename only |
| `apply_trim` | `enforce_seam_flux` | The new name says what it holds on target |
| `params.vector`, `params.velocity` | `vx`, `vy`, `vz` | Rejected, not silently accepted |
| `58_Turbulence_Statistics_Pipeline_Specification` | `58_Field_Statistics` | Page renamed |
| `/picurv-docs/<page>.html` | `/docs/picurv/<page>.html` | Old page URLs do not resolve |

@section p68_related_sec 7. Related Documentation

- **@subpage 56_Field_Identity_and_Layout_Catalog** — authoritative field definitions
- **@subpage 62_Capability_Status_Vocabulary** — status and evidence words
