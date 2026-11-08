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
      particles_per_cell: 10          # Surface | Volume; a PointSource takes count instead
      init_mode: "Surface"            # Surface | Volume | PointSource
      restart_mode: "init"            # init | load
      random_seed: 12345              # optional; seeds placement and Brownian draws
      point_source:
        x: 0.5
        y: 0.5
        z: 0.5
      fields:                         # optional; see 1.1
        Psi: {expression: 0.5}
```

Mapping to control flags:

- `particles_per_cell` -> `-numParticles`, resolved by the CLI from the grid: `Volume` seeds it
  times the cells of block 0, `Surface` it times the face cells of the first `INLET` face in
  face order (`-Xi +Xi -Eta +Eta -Zeta +Zeta`), the face the solver seeds. It is a mean: the
  particles are spread by volume or area, so on a graded grid small cells hold fewer, and the
  solver logs the per-cell minimum, mean and maximum. It may be fractional; the total is rounded.
- `count` -> `-numParticles`, for `PointSource` only, which has no cells; `Volume` and `Surface`
  refuse it
- `init_mode` -> `-pinit`
- `restart_mode` -> `-particle_restart_mode`
- `random_seed` -> `-particle_random_seed` (optional; integer `0`..`2147483647`, default `12345`)
- `point_source` -> `-psrc_x/-psrc_y/-psrc_z` (required when `init_mode` is `PointSource`); a
  physical position, divided by `length_ref` like the grid bounds it sits within
- `fields` -> `-particle_fields_*` (optional; see @ref p45_particle_values_sec)

Note: The interpolation method (`Trilinear` / `CornerAveraged` / `FaceFlux`) is configured in `solver.yml`, not `case.yml`. See **@subpage 08_Solver_Reference** and **@subpage 27_Trilinear_Interpolation_and_Projection**.

@section p45_particle_values_sec 1.1 Configured Particle Values (`fields`)

`fields` gives particle-carried fields a value when the population is created. Only fields
the particle carries from step to step can be set; every other particle field is
re-derived from the Eulerian fields or from particle location, so a value given to it
would not survive. Today that is `Psi`, the scalar that IEM micromixing relaxes (see
**@subpage 28_IEM_and_Statistical_Averaging**).

```yaml
sources:
  mixture: {expression: "0.5 + 0.5*tanh((y - y0)/delta)", constants: {y0: 1.0, delta: 0.1}}

models:
  physics:
    particles:
      particles_per_cell: 10
      init_mode: "Volume"
      fields:
        Psi: mixture                              # a name under sources:
        # Psi: {expression: 0.5}                  # or inline: a number
        # Psi: {expression: "0.5 + 0.5*sin(2*pi*xi)"}   # an expression
        # Psi:                                    # or regions over a background
        #   constants: {R: 0.1}
        #   expression:
        #     background: 0
        #     regions:
        #       - {shape: half_space, axis: x, at: 0.5, side: above, value: 1,
        #          edge: {profile: smoothstep, width: 0.02}}
        #       - {shape: ball, center: [0.5, 0.5, 0.5], radius: R, value: "uniform()"}
```

Each field's value is a source, as on an inlet or initial condition (@ref p33_source_ssec): the name of an entry
of the top-level `sources:`, or a mapping of `expression` and optional `constants`, or a `sum` of
sources. It is one value per particle, so it takes no `components`. An expression is a number, an
expression string, or `{background, regions}`, where each region in order paints its `value` over
what came before. `constants` may build on the case's top-level `constants:` and on `rho` and `nu`.

**Expressions** use the language of Eulerian sources
(**@subpage 33_Initial_Conditions**): `+ - * / % **`, comparisons (chains such as
`0.2 < x < 0.5` included), `and or not`, `abs sin cos tan tanh cosh sinh exp log erf sqrt floor
atan2 minimum maximum where(c, a, b)`, the series `sum(k, lo, hi, body)`, `besselj(n, x)
bessely(n, x) ber(x) bei(x)` (@ref p33_closed_form_ssec), and `pi`. A comparison or logical expression is 1 when true and 0 when
false. Numbers are decimal literals. A particle expression may use:

| Name | Value |
|---|---|
| `x`, `y`, `z` | particle position, physical |
| `t` | physical time of the event that creates the particle |
| `xi`, `eta`, `zeta` | arc-length fractions along the grid lines through the particle's cell, at the particle |
| `L_xi`, `L_eta`, `L_zeta` | those lines' lengths, physical |
| `pid` | particle ID |
| `uniform()`, `normal()` | a per-particle draw in [0, 1), or from N(0, 1) |

A draw is a pure function of `random_seed`, the particle ID, and the field; it does not
depend on the rank a particle lives on or the order particles are visited. Brownian
increments are drawn the same way, keyed on the step instead of the field, in a key space of
their own, so they never coincide with a field's draws; so are `Surface` and `Volume`
seeding points (section 4). `uniform(k)` and `normal(k)`
with a non-negative integer literal `k` name a stream shared across fields, so two fields
can draw the same number for one particle.

**Regions** are `half_space` (`axis`, `at`, `side: above | below`), `slab` (`axis`, `from`,
`to`), `box` (`min`, `max`), `ball` (`center`, `radius`), and `cylinder` (`axis`, `center`
in the other two coordinates, `radius`), all in physical coordinates. Without `edge`, a
region's boundary is sharp; `edge: {profile, width}` blends across a band of that width
centred on the boundary with `linear`, `smoothstep`, or `cosine`, reaching one half on the
boundary itself. Region numbers may name a constant.

**When it applies.** A value is applied once each particle's position is final: after
placement and location, and before interpolation and the first scatter, so the t=0 cell
mean already reflects it. A particle restored from a checkpoint (`restart_mode: load`)
keeps its saved value, and `fields` is ignored with a warning.

**Units.** Values are physical, in the field's units, and are divided by the field's
reference scale as they are stored (**@subpage 19_Nondimensionalization**). `Psi` is
dimensionless.

**Validation.** `picurv validate` refuses an unknown or re-derived field, an expression
outside the language, a region with a missing or unknown key, `fields` without particles or an inlet with particles,
and `fields` together with `solver.yml` `verification.sources.scalar`, which prescribes
`Psi` at every step. It warns when a `PointSource` value has no draw, because every
particle then starts at one point with one value.

**Diagnostics.** The launcher prints each lowered expression. After applying them, the
solver logs each field's particle count, mean, variance, minimum, and maximum (info
level), and `particle_initial_fields.csv` in the run's analysis directory records the same
per field and component, with a ten-bin histogram between the minimum and maximum.

**Evidence.** Analytically verified - `particle-field-sources-2026-11-06`: on a grid stretched with
ratio 2 on every axis, values set to the arc-length fractions `xi`, `eta`, `zeta` are exact to
2.2e-16 at 64 ranks.

Mapping to control flags: `-particle_fields_count`, then per field `-particle_fields_<i>_name`
and one `-particle_fields_<i>_expr_<c>` per component, holding the lowered expression.

@section p45_inject_sec 1.2 Injection at an Inlet

**Status: supported.** Injection itself is a small share of a step: 0.7% at 64 ranks and 1.5% at 256
with 1e8 particles. A step with injection still costs more, because particle location needs two
extra migration passes each step (GitHub #18).

An INLET face adds particles every step when its entry carries `particles`. Nothing needs
to be seeded, so a case without `particles_per_cell` starts flow-through from an empty domain, and two inlets
can carry different values:

```yaml
boundary_conditions:
- face: -Zeta
  type: INLET
  handler: prescribed_velocity
  source: {components: normal, expression: "6*y*(1-y)"}
  particles:
    particles_per_cell: 10               # mean per cell the entering fluid carries; a number or an expression in t
    fields: {Psi: {expression: "y"}}     # optional; the case's particles.fields otherwise
```

**Rate.** Each particle stands for an equal share of fluid, so particles enter in proportion
to the flow. `particles_per_cell` is a mean over the grid: with `N` cells of total volume `V`
it is the number density `n = particles_per_cell N / V`, the density a `Volume` seed of the
same value has. A face cell with volume flux `Q` receives on average `n Q dt` particles per
step and none where the flow leaves, which keeps that density under any velocity profile.
Expect a steady population of about `particles_per_cell N` after one flow-through time; on a
graded grid small cells hold fewer.
Nothing adds or removes particles to even this out, so recirculating or stagnant regions
fed slowly by the inlet can run short. `particles_per_cell_min` in `search_metrics.csv` shows it.

**Placement.** Within its face cell, a particle crosses at a point drawn in proportion to the
inward velocity the swarm's interpolation (`solver.yml` `interpolation.method`) gives there.
The particles already in the domain move by that velocity, so the layer they leave is thick
where the flow is fast and thin near a wall, and each face cell's particles refill it where
it opened. With `FaceFlux` the inward velocity is uniform across a face cell and the point
is uniform by area; `CornerAveraged` restages Ucat at the corners once per step for it. A
particle's value is evaluated on the face, where it crossed, at the step's time, as on a
seeded particle (@ref p45_particle_values_sec): `x y z t xi eta zeta L_xi L_eta L_zeta pid
uniform() normal()`, with `xi eta zeta` from the crossing cell. It then moves in by a uniform
share of its velocity times `dt`, so one step's particles do not travel as a sheet, and takes
the step's Brownian displacement with the diffusivity of its cell, which keeps the density
uniform up to the inlet. A face cell with flux but no inward velocity anywhere on it spreads
its particles by area and moves them at the cell's mean. The search that follows the particle
update locates it.

**Inlet policy.** An injecting inlet defaults to `particle_policy: reflect`. The inlet then
fixes the total flux entering, and a particle diffusing back across it stays in the domain.
With `remove` those particles are lost and the population settles below `n V`.

**Keyed draws and restart.** The number entering a face cell, each position, and each
identity follow from the seed, the face cell's position in the grid, and the step. Nothing is
stored, so a restart continues the same sequence, and the rank layout does not enter. The
flow still does: a face cell's mean count follows the flux, so round-off or a restart's flow
difference can change a count in the rare case the mean lies that close to the draw.

**Limits.** An injected ID is `2^40 + (step F + cell) 2^10 + k`, with `F` the grid's boundary
face cells and `k` numbering a face cell's particles in one step. A checkpoint stores IDs as
doubles, so a step is refused once an ID would pass `2^53`, and so is a step sending 1024 or
more particles through one face cell. Seeding at least `2^40` particles alongside an inlet is
refused. Inlets inject on block 0 only, and the statistics task `msd` is refused because it
measures displacement from one release. On a 32-bit-index build `logging.trace_particles` reads
IDs as 32-bit integers, so it cannot name an injected particle. An inlet's `fields` cannot be
combined with `verification.sources.scalar`, which prescribes `Psi` at every step.

**Diagnostics.** `particle_injection.csv` under `<run.analysis.metrics>` records each inlet's
flow rate, expected and injected particles and `Psi` added, every step. With `exited` and
`exited_psi_sum` in `search_metrics.csv` these close the particle and scalar budgets.

**Evidence.** Unit verified: `make unit-mpi` (`injection-follows-flux-and-is-keyed-on-the-grid`)
checks counts and identities against the formula at 2 and 4 ranks. Integration verified:
`make smoke` (`S4c`) fills an empty domain in an analytical uniform stream and closes the
particle and `Psi` budgets every step. Analytically verified - `particle-injection-2026-11-08`: in a
solved laminar channel with 1e7 particles at 64 ranks, the budget closes on every step of 1200, the
density is uniform from the inlet on within counting noise (near-inlet wall cells 0.999-1.014 against
noise 0.008), particle IDs are identical across two rank layouts and a restart, and a two-stream scalar
matches the erf reference (D_eff/D 1.0043 +- 0.0006 over two seeds against about 1.004). In an
analytical channel with 1e8 particles, IDs and positions are bit-identical across 4 x 4 x 16 and
2 x 4 x 32 ranks and across a restart at 2 x 4 x 32, with the budget exact.

Mapping to control flags: `-particle_inlet_count`, then per inlet `-particle_inlet_<n>_face`,
`-particle_inlet_<n>_particles_per_cell` and the field options under `-particle_inlet_<n>_fields_*`.

@section p45_mode_values_sec 2. Accepted `init_mode` Values

`picurv` accepts these exact canonical strings:

- `Surface`
- `Volume`
- `PointSource`

Enum mapping in C (`ParticleInitializationType`):

- `0`: `PARTICLE_INIT_SURFACE_RANDOM`
- `1`: `PARTICLE_INIT_VOLUME`
- `2`: `PARTICLE_INIT_POINT_SOURCE`

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

**Parameters it owns.** `particles_per_cell`, per face cell of the inlet face; the surface is
derived from the configured inlet face.

**Interactions.** Needs an identifiable inlet face; the startup log reports which face
was chosen. Particles are stratified over the inlet area (section 4): each inlet cell holds
its share of the count to within one particle, the point inside it is drawn uniform in area,
and every particle ID gets the same position at any rank layout. Particles start 1e-6 of a
cell inside the domain, in the cell they were placed in.

**Diagnostics.** "Inlet face for particle initialization identified as Face N" at
startup, and the per-step particle count.

**Evidence.** Analytically verified - `stratified-seeding-layout-2026-10-29`: positions, cells and IDs
bitwise equal by particle ID at 1, 1 x 1 x 2, 2 x 2 x 1, 4 x 1 x 1 and 2 x 2 x 2 ranks, on a
uniform and a graded inlet; every inlet cell held 25 or 26 of 16,000 particles against 25.6
expected, and the Kolmogorov-Smirnov statistic per axis was at most 0.35.
`surface-seeding-area-2026-10-29` established the uniform-by-area placement across rank
layouts and on a graded inlet. `stratified-seeding-layout-64rank-2026-10-29`: 1e6 Surface particles
bitwise equal by ID at 64 ranks in two layouts and on 1 rank.

**Limitations.** Concentrates particles at one boundary, so interior statistics take time
to become meaningful. The area of a curved face cell uses a 2x2 Gauss rule, exact only for
planar faces. Counts per cell are stratified, not Poisson-distributed.

@subsection p45_cap_volume_sub Volume

@anchor p45_cap_volume

**Identity.** `particles.init_mode: Volume` -> `-pinit 1`.

**What it does.** Distributes particles through the domain volume.

**When to choose it.** Whenever you want interior statistics immediately - dispersion,
mixing, or any case where waiting for particles to arrive from a boundary wastes the run.

**Parameters it owns.** `particles_per_cell`, the mean over the grid's cells.

**Interactions.** Particle density is uniform per unit physical volume on any grid.
Particles are stratified over the volume (section 4): each cell holds its share of the
count to within one particle, the point inside it is drawn uniform in volume by rejection
against the cell's largest Jacobian determinant, and every particle ID gets the same
position at any rank layout.

**Diagnostics.** Per-step particle count and the lost-particle counter.

**Evidence.** Production exercised - `examples/scatter_verification` seeds this way.
Analytically verified - `stratified-seeding-layout-2026-10-29`: 200,000 particles with positions, cells and IDs
bitwise equal by particle ID at five rank layouts on a uniform and a graded grid, and
345,600 on the periodic temporal mixing layer at 1, 2 and 4 ranks; on the uniform grid every
cell held 3 or 4 against 3.30 expected, and on the graded grid the density per tenth of x
stayed within 0.6% of uniform with Kolmogorov-Smirnov statistics at most 0.29.
`particle-seeding-restart-2026-09-21` first measured the uniformity and unique IDs.
`stratified-seeding-layout-64rank-2026-10-29`: 4e6 particles on a uniform and a graded grid bitwise
equal by ID at 64 ranks in two layouts and on 1 rank, no cell outside its share by more than one.
`particles-per-cell-2026-11-08`: `particles_per_cell: 10` on a graded grid of 524,288 cells at 64 ranks
seeded exactly 5,242,880 particles with unique IDs; counts per cell ran from 2 to 31 about a mean of
10, the spread volume weighting gives on a graded grid.

**Limitations.** No control over the spatial distribution beyond the grid itself. Counts
per cell are stratified, not Poisson-distributed.

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
`examples/brownian_motion` are built on this mode. Analytically verified - `particle-seeding-restart-2026-09-21`: every
particle starts exactly at the source.

**Limitations.** All particles share one origin, so early statistics are highly
correlated.

@section p45_mode_behavior_sec 4. Mode Behavior in C

Main setup flow:

1. @ref InitializeParticleSwarm creates the DMSwarm and fields.
2. `AssignInitialPropertiesToSwarm` seeds base particle state.
3. `PerformInitializedParticleSetup` settles/migrates and couples to Eulerian fields.

Mode details:

`Surface` (`0`) and `Volume` (`1`) stratify particles over a measure: the area of the inlet
face cells, or the volume of every cell.

- In global cell order (i fastest), with P(c) the measure of all cells before cell c, T the
  total and N the particle count, cell c seeds the particle IDs floor(N P(c) / T) to
  floor(N (P(c) + m_c) / T) - 1. Each cell holds its expected share to within one, and the
  IDs are 0 to N - 1.
- The measures are summed as integers, each cell's quantised against the largest cell, so
  P(c) and T, and with them each cell's particles, are identical at any rank layout. A rank
  seeds the particles of its own cells; for `Surface` a rank without part of the inlet seeds
  none, and `Surface` requires an identified INLET face from BC parsing.
- The point inside the cell is drawn uniform in its measure, by rejection, with every draw
  keyed on (`random_seed`, particle ID, attempt), and mapped to physical space via metric
  interpolation.

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

1. `LocateAllParticlesInGrid` performs location/migration, starting each walk in the
   cell the particle was seeded in.
2. @ref ParticleFieldPlanApply sets configured `fields` values, and
   @ref ParticleFieldPlanSummarize logs and records them.
3. `InterpolateAllFieldsToSwarm` assigns flow fields at particle positions.
4. the scatter updates Eulerian particle-derived fields.

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
Analytically verified - `particle-seeding-restart-2026-09-21`: a restart with `init` reseeded the run's own step-0
population and advected it by the carrier velocity to 5.5e-16.

**Limitations.** Any statistics accumulated by a previous run's particles are lost, because the particles they described no longer exist. Placement is keyed on `random_seed` and the particle ID, so with the seed unchanged a restart reseeds the same particles at the same positions the run started with, at any rank layout; change the seed for a statistically independent swarm.

@subsection p45_cap_restart_load_sub load

@anchor p45_cap_restart_load

**Identity.** `models.physics.particles.restart_mode: load`.

**What it does.** Restores particle positions, velocities, and swarm fields from the checkpoint named by the restart, continuing the same particles rather than seeding new ones.

**When to choose it.** When particle history matters - dispersion, mean-squared displacement, residence time - and the trajectory must continue rather than restart. Required for any statistic accumulated across a restart boundary.

**Parameters it owns.** None of its own; the checkpoint supplies the state. The restart source is chosen by `--restart-from` or `--continue`, not by this key.

**Interactions.** Requires a checkpoint that actually contains particle state: a run whose particles were disabled writes none. The restart matrix in @ref p45_restart_matrix_sec gives the full combination table against `eulerian_field_source`.

**Diagnostics.** A missing or empty particle group in the checkpoint is a fatal startup error naming the step, not a silent reseed. The restored count is reported in the banner.

**Evidence.** Integration verified - `make unit-particles` covers the restore path.
Analytically verified - `particle-seeding-restart-2026-09-21`: a 10-step run restarted for 10 more reproduced the
continuous 20-step swarm bitwise, lost particles included.

**Limitations.** The particle state itself restarts exactly, as the evidence shows on an
analytical carrier. Under a solved flow the Eulerian restart re-applies boundary
conditions, which perturbs the flow around 1e-7 regardless of solver tolerance, and the
particles inherit that; it is a floor on any claim of trajectory continuity across a
solved-flow restart. Brownian increments need no saved generator state: each is drawn from
`random_seed`, the particle ID and the step alone, the same keyed draw that `uniform()`
and `normal()` use. A restart with diffusion therefore continues the uninterrupted run's
increments, and a run on another rank count draws the same ones; `make smoke` checks both
bitwise on a zero-flow cloud. Production exercised - `keyed-brownian-continuation-2026-10-22`:
on a 200k-tracer channel at Sc 1 a restart is bitwise at the restart step and then differs
only by the fluid restart difference carried through the particle velocities (positions
5.3e-8 after 10 steps, against 0.037 with the earlier per-rank stream), and 64 and 144 ranks
agree to 9e-15. Runs made before this change drew from a per-rank stream and do not
reproduce.

@section p45_fields_sec 6. Swarm Fields Initialized at Startup

After position/PID/cell placeholders, initialization sets defaults for:

- `velocity` (vector),
- `weight` (vector),
- `Diffusivity` (scalar),
- `DiffusivityGradient` (vector),
- `Psi` (scalar).

These are the catalog defaults. Once positions are final, configured `fields` values
(@ref p45_particle_values_sec) replace the defaults of the fields they name.

Cell IDs start at `-1` until location confirms host cells.

@section p45_diagnostics_sec 7. Diagnostics and Sanity Checks

Check banner/log output for:

- selected particle initialization mode,
- identified inlet face for surface modes,
- migration/lost-particle counters after first steps.

Typical errors:

- no INLET face with surface modes,
- missing `point_source.{x,y,z}` for point source mode,
- restart mode not in `{init, load}`,
- a `fields` entry naming a field the runtime re-derives, or an expression outside the
  language (the error names the configuration path of the value).

The `[Particle IC]` log line and `particle_initial_fields.csv` show what the configured
values produced: a minimum and maximum that are equal where a spread was expected usually
means the expression does not depend on anything that varies across the particles.

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
- **@subpage 52_Run_Artifact_Lifecycle_Contract**
- **@subpage 26_Walking_Search_Method**
- **@subpage 27_Trilinear_Interpolation_and_Projection**
- **@subpage 39_Common_Fatal_Errors**
