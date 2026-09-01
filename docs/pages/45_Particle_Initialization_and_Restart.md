@page 45_Particle_Initialization_and_Restart Particle Initialization and Restart Guide

@anchor _Particle_Initialization_and_Restart

This page documents particle seeding, restart behavior, and early-step migration/settling logic.
It is written for both case authors and contributors working in `src/ParticleSwarm.c` and `src/ParticleMotion.c`.

@tableofcontents

@section p45_contract_sec 1. User Contract in `case.yml`

Particle controls live in:

```yaml
models:
  physics:
    particles:
      count: 50000
      init_mode: "Surface"            # Surface | Volume | PointSource | SurfaceEdges
      restart_mode: "init"            # init | load
      point_source:
        x: 0.5
        y: 0.5
        z: 0.5
```

Mapping to control flags:

- `count` -> `-numParticles`
- `init_mode` -> `-pinit`
- `restart_mode` -> `-particle_restart_mode`
- `point_source` -> `-psrc_x/-psrc_y/-psrc_z` (required when `init_mode` is `PointSource`)

Note: The interpolation method (`Trilinear` / `CornerAveraged`) is configured in `solver.yml`, not `case.yml`. See **@subpage 08_Solver_Reference** and **@subpage 27_Trilinear_Interpolation_and_Projection**.

@section p45_mode_values_sec 2. Accepted `init_mode` Values

`picurv` accepts these exact canonical strings:

- `Surface`
- `Volume`
- `PointSource`
- `SurfaceEdges`

Enum mapping in C (`ParticleInitializationType`):

- `0`: `PARTICLE_INIT_SURFACE_RANDOM`
- `1`: `PARTICLE_INIT_VOLUME`
- `2`: `PARTICLE_INIT_POINT_SOURCE`
- `3`: `PARTICLE_INIT_SURFACE_EDGES`

@section p45_entries_sec 3. Particle Initialization Mode Entries

@htmlinclude generated/capability_inventory_particle_init_mode.html

@subsection p45_cap_surface_sub Surface

@anchor p45_cap_surface

**Identity.** `particles.init_mode: Surface` -> `-pinit 0` ->
`ParticleInitializationType`.

**What it does.** Seeds particles across a bounding surface of the domain.

**When to choose it.** Studying what enters through a face - inlet seeding, or transport
from a boundary into the interior. Choose `Volume` when you want the interior populated
from the start.

**Parameters it owns.** The particle count; the surface is derived from the configured
inlet face.

**Interactions.** Needs an identifiable inlet face; the startup log reports which face
was chosen.

**Diagnostics.** "Inlet face for particle initialization identified as Face N" at
startup, and the per-step particle count.

**Evidence.** Implemented only.

**Limitations.** Concentrates particles at one boundary, so interior statistics take time
to become meaningful.

@subsection p45_cap_volume_sub Volume

@anchor p45_cap_volume

**Identity.** `particles.init_mode: Volume` -> `-pinit 1`.

**What it does.** Distributes particles through the domain volume.

**When to choose it.** Whenever you want interior statistics immediately - dispersion,
mixing, or any case where waiting for particles to arrive from a boundary wastes the run.

**Parameters it owns.** The particle count.

**Interactions.** Particle density follows the cell distribution, so a graded grid gives a
non-uniform physical density.

**Diagnostics.** Per-step particle count and the lost-particle counter.

**Evidence.** Production exercised - `examples/scatter_verification` seeds this way.

**Limitations.** No control over the spatial distribution beyond the grid itself.

@subsection p45_cap_pointsource_sub PointSource

@anchor p45_cap_pointsource

**Identity.** `particles.init_mode: PointSource` -> `-pinit 2`, with coordinates from
`-psrc_x`/`-psrc_y`/`-psrc_z`.

**What it does.** Releases all particles from a single specified point.

**When to choose it.** Plume and dispersion studies, and the verification cases where a
known release point makes an analytic comparison possible.

**Parameters it owns.** The point-source coordinates.

**Interactions.** The point must lie inside the domain; a point outside it produces
immediate particle loss.

**Diagnostics.** The lost-particle counter is the first signal of a misplaced source.

**Evidence.** Production exercised - `examples/drift_uniform_flow` and
`examples/brownian_motion` are built on this mode.

**Limitations.** All particles share one origin, so early statistics are highly
correlated.

@subsection p45_cap_surfaceedges_sub SurfaceEdges

@anchor p45_cap_surfaceedges

**Identity.** `particles.init_mode: SurfaceEdges` -> `-pinit 3`.

**What it does.** Seeds particles along the edges of a bounding surface rather than across
its face.

**When to choose it.** Exercising the locate-and-migrate machinery, where edges and
corners are the hard cases. It is a diagnostic seeding mode more than a physical one.

**Parameters it owns.** The particle count.

**Interactions.** Edge and corner cells are exactly where the walking search is most
likely to struggle, which is the point.

**Diagnostics.** Search metrics at **@subpage 53_Search_Robustness_Metrics_Reference**.

**Evidence.** Implemented only.

**Limitations.** Not a physically motivated distribution.

@section p45_mode_behavior_sec 4. Mode Behavior in C

Main setup flow:

1. @ref InitializeParticleSwarm creates the DMSwarm and fields.
2. `AssignInitialPropertiesToSwarm` seeds base particle state.
3. `PerformInitializedParticleSetup` settles/migrates and couples to Eulerian fields.

Mode details:

`Surface` (`0`):

- requires an identified INLET face from BC parsing,
- ranks not servicing inlet place particles at inlet-center fallback (`CMx_c/CMy_c/CMz_c`) before migration,
- `ReinitializeParticlesOnInletSurface` re-spreads particles on inlet partitions after first settlement.

`SurfaceEdges` (`3`):

- deterministic placement by particle ID on inlet-face lattice,
- if deterministic target is non-local after migration, code falls back to random inlet placement.

`Volume` (`1`):

- random logical coordinates inside locally owned cells,
- mapped to physical space via metric interpolation.

`PointSource` (`2`):

- all particles start at fixed coordinates `(psrc_x, psrc_y, psrc_z)`.

@section p45_restart_matrix_sec 5. Restart Behavior Matrix

`InitializeParticleSwarm` behavior depends on `StartStep` and restart mode:

| `StartStep` | `particle_restart_mode` | Behavior |
|---|---|---|
| `0` | any | initialize new population |
| `>0` | `init` | initialize new population in restarted flow |
| `>0` | `load` | load particle fields from restart files |

Operational note:

- For a run completed through step `N`, use `start_step: N`.
- Choose `restart_mode: load` to continue the existing particle swarm.
- Choose `restart_mode: init` to reseed a fresh particle population in the restarted flow field.
- Point `picurv` at the previous run with `--restart-from <previous_run_dir>` on the CLI (or `--continue` to resume the most recent run of the same case).

For loaded particles, fast migration path:

- @ref MigrateRestartParticlesUsingCellID uses stored cell IDs to migrate directly before full walking-search fallback.

@section p45_settle_sec 5. Early-Step Settlement and Coupling

For initialized particles (`StartStep == 0` path):

1. `LocateAllParticlesInGrid` performs location/migration.
2. surface modes call @ref ReinitializeParticlesOnInletSurface.
3. statuses are reset and location pass is repeated.
4. `InterpolateAllFieldsToSwarm` assigns flow fields at particle positions.
5. optional scatter updates Eulerian particle-derived fields.

For loaded particles:

1. `MigrateRestartParticlesUsingCellID` fast migration,
2. `LocateAllParticlesInGrid` resolves invalid/missing cases,
3. interpolation/scatter synchronize coupling state before stepping.

@section p45_cap_restart_sec 5.1 Particle Restart Mode Entries

@htmlinclude generated/capability_inventory_particle_restart_mode.html

@subsection p45_cap_restart_init_sub init

@anchor p45_cap_restart_init

**Identity.** `models.physics.particles.restart_mode: init` (the default).

**What it does.** Seeds the swarm from `init_mode` at the start of the run, regardless of whether Eulerian fields are being restarted. Particle state in a checkpoint is ignored.

**When to choose it.** Whenever the particles are the thing you are varying: re-seeding onto a restarted flow field lets you launch a fresh swarm into developed turbulence without re-running the spin-up. Also the correct choice for any first run.

**Parameters it owns.** None of its own. It selects which branch consumes `models.physics.particles.init_mode` and its parameters - see @ref p45_entries_sec.

**Interactions.** Independent of `eulerian_field_source`: a run may load fields and still init particles. Combined with a field restart it is the standard 'new swarm, developed flow' configuration.

**Diagnostics.** The startup banner reports the resolved particle restart mode and the seeded count. A swarm that seeds at the wrong size shows here, before the first step.

**Evidence.** Integration verified - `make unit-particles` covers the seeding branches.

**Limitations.** Any statistics accumulated by a previous run's particles are lost, because the particles they described no longer exist.

@subsection p45_cap_restart_load_sub load

@anchor p45_cap_restart_load

**Identity.** `models.physics.particles.restart_mode: load`.

**What it does.** Restores particle positions, velocities, and swarm fields from the checkpoint named by the restart, continuing the same particles rather than seeding new ones.

**When to choose it.** When particle history matters - dispersion, mean-squared displacement, residence time - and the trajectory must continue rather than restart. Required for any statistic accumulated across a restart boundary.

**Parameters it owns.** None of its own; the checkpoint supplies the state. The restart source is chosen by `--restart-from` or `--continue`, not by this key.

**Interactions.** Requires a checkpoint that actually contains particle state: a run whose particles were disabled writes none. The restart matrix in @ref p45_restart_matrix_sec gives the full combination table against `eulerian_field_source`.

**Diagnostics.** A missing or empty particle group in the checkpoint is a fatal startup error naming the step, not a silent reseed. The restored count is reported in the banner.

**Evidence.** Integration verified - `make unit-particles` covers the restore path.

**Limitations.** Restart is not bit-exact: boundary conditions are re-applied on load, which introduces a perturbation around 1e-7 regardless of solver tolerance. That is a floor on any claim of trajectory continuity across a restart.

@section p45_fields_sec 6. Swarm Fields Initialized at Startup

After position/PID/cell placeholders, initialization sets defaults for:

- `velocity` (vector),
- `weight` (vector),
- `Diffusivity` (scalar),
- `DiffusivityGradient` (vector),
- `Psi` (scalar).

Cell IDs start at `-1` until location confirms host cells.

@section p45_diagnostics_sec 7. Diagnostics and Sanity Checks

Check banner/log output for:

- selected particle initialization mode,
- identified inlet face for surface modes,
- migration/lost-particle counters after first steps.

Typical errors:

- no INLET face with surface modes,
- missing `point_source.{x,y,z}` for point source mode,
- restart mode not in `{init, load}`.

@section p45_extension_sec 8. Contributor Extension Points

If adding a new particle initialization mode:

1. add validation and mapping in `picurv_cli/core.py`,
2. extend `ParticleInitializationType` and parser wiring in C,
3. implement placement logic in `InitializeParticleBasicProperties` and any inlet reinit path,
4. update logging string mappings and tests,
5. update this page and related method references.

For the full selector extension checklist, see **@subpage 50_Modular_Selector_Extension_Guide**.

@section p45_refs_sec 9. Related Pages

- **@subpage 33_Initial_Conditions**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 52_Run_Lifecycle_Guide**
- **@subpage 26_Walking_Search_Method**
- **@subpage 27_Trilinear_Interpolation_and_Projection**
- **@subpage 39_Common_Fatal_Errors**
