@page 33_Initial_Conditions Initial Conditions

@anchor _Initial_Conditions

This page documents how PICurv initializes Eulerian velocity fields and particles at startup.
It covers both user-facing YAML inputs and the C implementation path that applies them.

@tableofcontents

@section p33_overview_sec 1. Where Initialization Happens

Startup sequence:

1. `picurv_cli/core.py` resolves the `source`.
2. A file or generator term is staged in the filename layout expected by @ref ReadFieldData.
3. @ref InitializeEulerianState chooses fresh solve, restart load, or analytical initialization.
4. On a fresh solve, @ref PopulateInitialUcont populates `Ucont`; existing finalization then applies
   boundary values and derives the remaining velocity state.
5. If particles are enabled, @ref InitializeParticleSwarm runs independently after Eulerian setup.

@section p33_euler_sec 2. Eulerian Field Initialization (`properties.initial_conditions`)

An initial condition is one `source` and optional `constraints`, the same definition a
@ref p44_cap_prescribed_velocity "prescribed_velocity" inlet takes. Its forms:

| `source` | Field |
|---|---|
| `{components, constants, expression}`, or a name under `sources:` | evaluated by the solver (@ref p33_source_ssec) |
| `{generator, params}` | written by a Python generator, staged (@ref p33_entries_sec) |
| `{type: file, field, path}` | read from a velocity file, staged |
| `{sum: [...]}` | expressions added to at most one staged field (@ref p33_source_ssec) |

The same language gives particle-carried fields their initial values: one expression per particle, which also reads `pid` and may draw `uniform()` and `normal()` (@ref p45_particle_values_sec).

A velocity file:

```yaml
properties:
  initial_conditions:
    source: {type: file, field: Ucat, path: initial_conditions/velocity.dat}
```

The input must be one PETSc binary `.dat` vector readable by @ref ReadFieldData. `Ucat` inputs
are converted by @ref Cart2Contra; `Ucont` inputs are used directly. The first implementation
supports one block only because it intentionally reuses the existing single-field
@ref ReadFieldData path.

The payload is physical: `Ucat` in velocity units, `Ucont` in volume flux per face. It is
divided by `velocity_ref` (or `velocity_ref length_ref^2`) as it is staged. A field saved by
an earlier PICurv run is in that run's solver units, so it declares the scales it was
written in, and staging rescales it to this case's:

```yaml
    source:
      type: file
      field: Ucat
      path: inputs/initial_conditions/velocity.dat
      source_case: runs/<earlier_run>/config/case.yml   # or: velocity_scale: 1.5
```

`velocity_scale` needs `length_scale` too for a `Ucont` payload, and `source_case` and the
explicit scales are mutually exclusive.

@subsection p33_source_ssec Source: a velocity written as expressions

@anchor p33_source

An expression source is given inline or by name from the case's top-level `sources:`:

```yaml
sources:
  duct_flow:
    components: curvilinear          # cartesian (default) | curvilinear
    constants: {U: 1.0}
    expression: {u_zeta: "36*U*xi*(1-xi)*eta*(1-eta)"}

properties:
  initial_conditions:
    source: duct_flow
    constraints: {scale: 1.0}        # optional
```

The runtime evaluates it at every cell centre, in physical units, over `x y z`, `t` (the start
time), the arc-length fractions `xi eta zeta` and the line lengths `L_xi L_eta L_zeta`, and
writes `Ucat`; @ref Cart2Contra then derives `Ucont`. A curvilinear source sets the fluxes
instead, as below. The fractions follow distance along each
grid line, so a profile written in them is the physical profile on a stretched grid, unlike the
index fractions a grid-line count would give.
The `normal` frame is a face's frame and is refused here.

`source` may instead be a `sum` of terms: sources (names or mappings) in one frame, which are
added, and at most one staged field: a `{generator, params}` term or a `{type: file, ...}` term.
The field is staged as it would be on its own, and the runtime adds the evaluated source to it;
under a staged `Ucont` the expressions must be curvilinear. A mean profile plus the turbulent
fluctuations of a wall seed:

```yaml
properties:
  initial_conditions:
    source:
      sum:
        - inflow                                  # a named source: the mean
        - generator: channel_spectral_velocity    # fluctuations only
          params: {seed: 12345, wall_axes: [j], streamwise_axis: k,
                   perturbation_rms: 0.1, initial_spectra: [...]}
    constraints: {bulk_velocity: 1.0, project: true}
```

`constraints` apply to the whole field, in order. `scale` multiplies it. `bulk_velocity` [m/s], a
number or an expression in `t`, means what it means on an inlet (@ref p44_cap_prescribed_velocity):
the source's expressions are rescaled, and anything staged is left as it is, so that the mean of the
whole equals it. Here the mean is over the cells with volume weights, along the source's own
direction at each cell; an inlet's is over its face with area weights, along the normal. It needs an
expression to rescale. `project`, which defaults to true for a sum and false otherwise, makes it divergence-free before the first step
with the time step's own pressure projection: one Poisson solve for the initial `Ucont`, whose
boundary fluxes are already set by the boundary conditions, and the correction of the interior
faces. The potential is discarded, so pressure keeps its initial value, and the solve's
convergence record is overwritten by the first step's. Without `project`, the first step's
projection makes the field divergence-free instead, after one momentum step on the field as
given. Projecting evens each plane's flux to the same value without changing `bulk_velocity`'s
volume mean on a periodic axis; between an inlet and an outlet, conservation sets each plane's
flux to the inlet's, so a `bulk_velocity` other than the inlet's does not survive it. A source's
`pressure` expression [Pa] sets the initial `P` at the cell centres, and is otherwise zero.

Every component reaches the velocity through the faces, evaluated at each face centre. A Cartesian
velocity is dotted with the face metric. A `curvilinear` component is the speed through the grid face
normal to its axis, on every target: the flux is that component times the face area. A `contravariant`
component is the volume flux through that face itself [m^3/s], so it writes `Ucont` directly:

```yaml
    source: {components: contravariant, expression: {u_zeta: "0.01*(1 + x*y)"}}
```

An initial condition takes the contravariant frame; an inlet and an analytical field do not, and
`bulk_velocity` reads a velocity, so it does not apply to one. A uniform
`u_zeta` is therefore discretely divergence-free on a curved grid whose zeta faces keep their area
along the line, and needs no projection.

@subsection p33_closed_form_ssec Closed-form profiles: series and special functions

A source's expressions, like particle values, may use a series and four special functions,
which the solver evaluates:

| Construct | Meaning |
|---|---|
| `sum(k, lo, hi, body)` | `body` summed for the integers `k = lo, ..., hi`, with `k` in scope only in `body`. The bounds are integer numbers (a constant inlines to one); at most 1000000 terms; series nest up to four deep, each with its own index |
| `besselj(n, x)`, `bessely(n, x)` | Bessel functions of the first and second kind, of integer order `n` (another order gives a value that is not finite, which stops the run) |
| `ber(x)`, `bei(x)` | Kelvin functions of order zero, `ber(x) + i bei(x) = J0(x e^{3 pi i/4})`, from their power series below `|x| = 30`, accurate to 1e-13 relative there, and their asymptotic expansion above it |

`picurv` validates these before launch and the solver evaluates them. Two exact laminar profiles,
in the form `tests/test_expression_language.py` and `make unit-particles` check:

```yaml
sources:
  duct_poiseuille:
    components: normal
    constants: {U: 1.0, beta: 1.0}
    expression: >-
      U*sum(k, 0, 99, (-1)**k/(2*k+1)**3
      *(1 - (exp(abs((2*k+1)*pi*beta*(eta - 0.5)) - (2*k+1)*pi*beta/2)
      + exp(-abs((2*k+1)*pi*beta*(eta - 0.5)) - (2*k+1)*pi*beta/2))/(1 + exp(-(2*k+1)*pi*beta)))
      *cos((2*k+1)*pi*(xi - 0.5)))
      / sum(k, 0, 99, 2/(pi*(2*k+1)**4)*(1 - 2*tanh((2*k+1)*pi*beta/2)/((2*k+1)*pi*beta)))
  womersley:
    components: normal
    constants: {U: 1.0, K: 1.0, omega: 6.283185307179586, alpha: 5.0, R: 0.5}
    expression: >-
      2*U*(1 - (x**2 + y**2)/R**2)
      + (K/omega)*((ber(alpha*sqrt(x**2 + y**2)/R)*bei(alpha) - bei(alpha*sqrt(x**2 + y**2)/R)*ber(alpha))
      /(ber(alpha)**2 + bei(alpha)**2)*cos(omega*t)
      + (1 - (ber(alpha*sqrt(x**2 + y**2)/R)*ber(alpha) + bei(alpha*sqrt(x**2 + y**2)/R)*bei(alpha))
      /(ber(alpha)**2 + bei(alpha)**2))*sin(omega*t))
```

`duct_poiseuille` is the fully developed flow in a rectangular duct whose section has width
along `xi` and height along `eta`, `beta` = height / width, at bulk velocity `U`: the
cosh-cosine series of the exact solution over its own area mean, both to 100 odd terms, with
`cosh(a)/cosh(b)` written so that it cannot overflow. On a square section it matches the
double sine series of the same solution to 1e-5 of `U`; written in the
arc-length fractions, it is the physical profile on a stretched grid. `womersley` is the flow in
a pipe of radius `R` along `z`, with its axis through `x = y = 0`: a Poiseuille mean of bulk `U` plus the
exact oscillatory part driven by `-(1/rho) dp/dz = K cos(omega t)` [m/s^2], with Womersley
number `alpha = R sqrt(omega/nu)`. It matches the complex solution
`Re{(K/(i omega))(1 - J0(i^{3/2} alpha s)/J0(i^{3/2} alpha)) e^{i omega t}}` to 1e-12, and at
small `alpha` it tends to the quasi-steady Poiseuille profile. A sweep does not place a pipe there by default; on an inlet, `r` replaces
`sqrt(x**2 + y**2)` and the profile holds wherever the pipe lies, with its centre on the face's. Either serves an inlet as
written and an initial condition with `components: curvilinear` and the expression under
`u_zeta`, on a grid whose `zeta` lines run along the flow.

@subsection p33_file_or_import_ssec File, Or Import?

Every Python provider ends in a `.dat` vector; what differs is who guarantees it fits the case.

| | `type: file` | `generator: resampled_velocity` |
|---|---|---|
| You supply | a vector already in PICurv's layout for this grid | a field from anywhere: another code, resolution, or layout |
| PICurv does | reads it verbatim | resamples and filters it onto this grid, projects it, writes the vector |
| Typical source | a field PICurv itself produced | a published DNS, an experiment, another code's snapshot |

A `type: file` vector sized for another grid is refused when it is staged, with a pointer to
@ref p33_cap_gen_resampled_velocity.

Your own generator:

```yaml
properties:
  initial_conditions:
    source:
      generator: ic_gen
      params:
        field: Ucat
        script: tools/custom_ic.py
        config_file: config/initial_conditions/custom.cfg
```

`picurv` chooses where the payload goes: it is always written to
`<run.inputs>/initial_condition/initial_condition.generated.dat`. The retired
destination keys `output_file`, `summary_json`, and `spectrum_csv` are refused
with an error rather than ignored.

The launcher runs the script, resolved against the case, as:

```text
python <script> -c <config_file> --field Ucat|Ucont --output <run.inputs>/initial_condition/initial_condition.generated.dat --grid <run.inputs>/grid/grid.run --length-ref <L> --velocity-ref <U> --petsc-index-bytes <4|8>
```

It receives the non-dimensional staged grid and the case's reference scales, and writes solver
units. `grid.mode: file` and `grid.mode: grid_gen` provide that grid directly; for single-block
`grid.mode: programmatic_c`, the launcher materializes a nondimensional
`<run.inputs>/grid/grid.run` bridge from scalar `programmatic_settings` first.
`picurv run --solve` materializes the result after grid preparation, and `picurv precompute --case ...`
materializes and stages the same artifact without running the solver.

Repository spectral provider, for a triply periodic box:

```yaml
properties:
  initial_conditions:
    source:
      generator: spectral_random_velocity
      params:
        field: Ucat
        seed: 12345
        spectrum: {type: k4_exponential, k0: 4.0, k_cut: 20.0}
        projection: {type: solenoidal, operator: picurv_discrete}
        normalization: {type: component_rms, target: 1.0}
        remove_mean: true
```

It draws a seeded random field, shapes it to the spectrum envelope
(`k4_exponential` is `E(k) ~ k^4 exp(-2 (k/k0)^2)`, nothing above `k_cut`), projects it
to be divergence-free under the chosen operator, and scales it to the requested
component RMS. `random.mean` and `normalization.target` are velocities and `k0`,
`k_cut` wavenumbers, all physical, and are converted before the provider runs on the
staged non-dimensional grid. It requires every face geometric-periodic and a fresh 3D run, and writes
`<run.analysis>/metrics/initial_condition_summary.json` and the staged spectrum beside
the field. `examples/decaying_isotropic_turbulence` configures it; its README covers the
keys. With `operator: picurv_discrete` the solver's own step-0 divergence is at round-off;
`initial-conditions-2026-09-21` measured 3.1e-14 against a flux scale of 4.3, component RMS 1.000000, the shell
spectrum within 4.7% (energy-weighted) of the envelope, and the same seed reproducing
the same field.

@section p33_entries_sec 3. Staged Field Entries

The `file` term and the expression config `ic_gen` once evaluated were checked against the field they claim to produce by
`initial-conditions-2026-09-21`: on a programmatic duct the staged interior face fluxes equal their
definitions to 3.5e-18 or exactly, and a `file` term reproduces the generated run it was taken from
bitwise. Step-0 `Ucat` next to walls and inlets reflects the boundary faces, because it is
reconstructed from the fluxes after boundary conditions apply.

@htmlinclude generated/capability_inventory_initial_target_field.html

@subsection p33_cap_field_ucat_sub Ucat

@anchor p33_cap_field_ucat

**Identity.** `field: Ucat` on a file or `ic_gen` term -> the Cartesian velocity field is what is
staged (internal code `0`).

**What it does.** Stages the Cartesian velocity, from which the contravariant fluxes are then derived.

**When to choose it.** The normal choice, and the one that matches how the field modes are
expressed - a uniform or parabolic profile is naturally stated in Cartesian components.

**Parameters it owns.** None; it selects which field the term holds.

**Interactions.** `Ucat` is the derived representation, so initializing it requires a
conversion to the evolved `Ucont` fluxes before the first step. See
**@subpage 20_Grid_Cell_Architecture_Guide**.

**Diagnostics.** Step-0 output shows the initialized field.

**Evidence.** Production exercised - `examples/flat_channel`. Analytically verified -
`initial-conditions-2026-09-21`: `Ucat` expressions staged the fluxes their Cartesian
values define.

**Limitations.** On a strongly curvilinear grid the conversion from a Cartesian profile to
fluxes is not exact in the sense of preserving the intended profile shape.

@subsection p33_cap_field_ucont_sub Ucont

@anchor p33_cap_field_ucont

**Identity.** `field: Ucont` on a file or `ic_gen` term -> the contravariant flux field is what is
staged (internal code `1`).

**What it does.** Stages the evolved flux variable directly,
skipping the Cartesian-to-flux conversion.

**When to choose it.** When you want exact control of the initial fluxes - for instance
setting a precise volumetric flux - rather than of the Cartesian velocity.

**Parameters it owns.** None.

**Interactions.** `Ucont` is what the solver actually evolves, so this path avoids one
conversion. The field modes are still expressed in the same components, which is what
makes this the less intuitive of the two.

**Diagnostics.** Step-0 output; check the derived `Ucat` looks as intended.

**Evidence.** Analytically verified - `initial-conditions-2026-09-21`: `ic_gen` `Ucont` expressions stage exactly,
and a `file` source of `Ucont` reproduces the generated run bitwise. No shipped example
selects it.

**Limitations.** Easier to get wrong than `Ucat`, because the values are fluxes rather
than velocities: a face flux carries the face area, so a velocity typed where a flux is
expected is off by that area. Only uniform Cartesian grids were checked.

@section p33_euler_modes_sec 4. C Runtime Entry Point

@ref PopulateInitialUcont is the fresh-solve entry point. It reads a staged field when `-ic_dir` names
one, through @ref ReadFieldData; when that field is `Ucat`, @ref Cart2Contra converts it to `Ucont`.
It then evaluates the source's expressions, through the same `PicurvVelocitySource` an inlet uses.
After that point, the existing finalization path treats every initial condition identically.

@section p33_restart_modes_sec 5. Authority and Restart Branches

In @ref InitializeEulerianState "InitializeEulerianState":

- `eulerian_field_source=solve` and `StartStep == 0` consumes `initial_conditions`.
- `eulerian_field_source=solve` and `StartStep > 0` consumes the restart path.
- `eulerian_field_source=load` consumes the existing load path.
- `eulerian_field_source=analytical` consumes the analytical initializer.

Thus `eulerian_field_source` supersedes `initial_conditions`; the launcher does not materialize
a configured file or `ic_gen` artifact when another source has authority.

Operational note:

- `StartStep` identifies the saved restart state being loaded, not the first new step to compute.
- If a run completed through step `N`, restart with `start_step: N`; the first newly advanced step will be `N+1`.

@section p33_particle_link_sec 6. Particle Initialization Relation

Particle initialization is configured in `case.yml -> models.physics.particles`, but executed by a separate subsystem.

For full particle mode and restart details, use:

- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**

@section p33_checks_sec 7. Practical Checks

After startup, confirm:

- the banner names the initial condition as the Eulerian source, and any staged file,
- no warnings about missing inlet face (for inlet-driven setups),
- first output step has non-empty `Ucat`/`Ucont` fields.

Common pitfalls:

- forgetting that initialization sets the interior only; boundary handlers then overwrite face values,
- providing a PETSc vector whose size does not match the target DM,
- comparing a velocity directly to `Ucont` without accounting for metric-face scaling.

@section p33_refs_sec 7. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 32_Analytical_Solutions**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 50_Modular_Selector_Extension_Guide**

@section p33_wall_spectral_sec Wall-bounded spectral seeds

@htmlinclude generated/capability_inventory_initial_file_generator.html

`channel_spectral_velocity` uses one wall axis; `duct_spectral_velocity` uses two.
The remaining axes are periodic, with driving allowed only on `streamwise_axis`, or the
streamwise axis runs from an inlet on its negative face to an outlet on its positive face. The
seed's modes are periodic along that axis either way; between an inlet and an outlet, project
the sum (@ref p33_source).
Both use `Ucat` on one Cartesian block. Uniform spacing is required along periodic
axes; wall-normal axes may stretch. `grid_gen` provides an inspectable staged grid.

```yaml
source:
  sum:
    - {components: cartesian, expression: {w: "6*eta*(1-eta)"}}
    - generator: channel_spectral_velocity
      params:
        field: Ucat
        seed: 12345
        wall_axes: [j]
        streamwise_axis: k
        perturbation_rms: 0.1
        wall_modes: 3
        spectrum: {type: k4_exponential, k0: 4.0, k_cut: 8.0}
        initial_spectra:
          - {task: plane_spectrum, axes: [i, k], fixed_indices: {j: 32}, subtract_mean: sample}
constraints: {bulk_velocity: 1.0, project: false}
```

These generators write the fluctuation alone; the mean is a source, normalized by `bulk_velocity`.
Parameters are physical: `perturbation_rms` is a velocity and the spectrum wavenumbers are per length, converted before the provider runs. `perturbation_rms` is
`sqrt(volume_mean(u'^2+v'^2+w'^2)/3)`, not the RMS of each component independently.
The spectral envelope shapes a seeded vector potential along periodic directions;
wall basis functions are `sin(n*pi*t)*sin(pi*t)`, evaluated at physical cell
centers: the least envelope that keeps the potential and its wall-normal derivative
zero at the wall, so the seed has the physical near-wall behaviour, wall-parallel
`u' ~ y` and wall-normal `v' ~ y^2`. A curl using the separable face-average divergence operators produces a
perturbation in the discrete divergence nullspace, including zero wall-face flux.
Periodic modes at or above one third of the cell count are removed. Wall mode count
plus one (the highest sine index added by the envelope) must fit below half the wall-axis cell count; stretched
meshes still need a physical resolution check. The envelope is not a prescribed final
isotropic velocity spectrum: the curl, wall functions, and stretching change it.

Perturbation net flux is zero, so the bulk velocity is the mean's. Match the driven boundary's
target flux to bulk velocity times cross-section. A curl of periodic modes plus a streamwise-uniform
mean is already divergence-free, so a periodic case need not project. Wall dummy values are odd reflections, periodic dummies wrap. Runtime retains
ownership of boundary enforcement and Cartesian/contravariant reconstruction.

`initial_spectra` is a required nonempty list using the plane/line task syntax from
@ref p10_spectra_sec. Initial measurements cannot use a statistics-window mean.
Staging measures the generated Ucat, before runtime reconstruction. Compare against
step-zero checkpoint spectra to distinguish initialization reconstruction from
subsequent physical evolution. The summary and spectra use the canonical metrics and
spectra directories, including through precomputed asset publication and reuse.
There is no per-timestep injection, and transition/sustained turbulence is not guaranteed.

@subsection p33_cap_gen_ic_gen_sub ic_gen

@anchor p33_cap_gen_ic_gen

**Identity.** `source: {generator: ic_gen}` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the staged-field runtime path.

**What it does.** Runs your own script, which writes a `Ucat` or `Ucont` vector for the staged grid.

**When to choose it.** A field that needs a program. A field with a formula is an expression source
(@ref p33_source_ssec), whose language this page and the solver share.

**Parameters it owns.** Requires script and config_file; field selects Ucat/Ucont; cli_args are appended.

**Interactions.** The script reads the staged PICGRID and writes the PETSc vector; a missing input or a
vector sized for another grid fails before the solver launches.

**Diagnostics.** The script's output is reported when it fails.

**Evidence.** Unit verified - `make test-python` covers the script contract and the CLI route.
`initial-conditions-2026-09-21` measured the expression config this generator once evaluated, whose
`Ucat` and `Ucont` fields the source frames reproduce: exactly on a flat channel for a test field with no variation along each face's normal, and to 3e-15 for
`Ucont` on a curved one. A `Ucat` expression differs there by 5e-7 in flux, because a source evaluates
at the face.

**Limitations.** Nothing is implied about a field a script writes.

@subsection p33_cap_gen_spectral_random_velocity_sub spectral_random_velocity

@anchor p33_cap_gen_spectral_random_velocity

**Identity.** `source: {generator: spectral_random_velocity}` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the staged-field runtime path.

**What it does.** The existing seeded triply periodic random spectral velocity.

**When to choose it.** Choose this for DIT; channel and duct providers supply wall-compatible seeds.

**Parameters it owns.** seed, random, spectrum, projection, normalization and remove_mean retain their existing contracts.

**Interactions.** Requires six geometric-periodic faces and a uniform Cartesian grid; restart remains authoritative.

**Diagnostics.** The IC summary reports selected-operator divergence and realized energy; shell spectra measure the generated field.

**Evidence.** Unit verified - `make test-python`. Production exercised in
`examples/decaying_isotropic_turbulence`. Analytically verified -
`initial-conditions-2026-09-21`: fluctuation RMS 1.000000 and kinetic energy 1.500000 as
requested, the solver's own step-0 divergence 3.1e-14 against a flux scale of 4.3 under
`operator: picurv_discrete`, the shell spectrum within 4.7% (energy-weighted) of the
`k4_exponential` envelope with nothing above `k_cut`, and identical seeds reproducing
identical fields.

**Limitations.** A continuum solenoidal field need not be discretely solenoidal - choose
`operator: picurv_discrete` when the solver's own divergence must vanish - and runtime
reconstruction changes the energy from the staged value, which the summary predicts
separately. A startup field, not developed turbulence.

@subsection p33_cap_gen_resampled_velocity_sub resampled_velocity

@anchor p33_cap_gen_resampled_velocity

**Identity.** `source: {generator: resampled_velocity}` dispatches through
`GENERATED_IC_PROVIDERS` into `generate_resampled_velocity()` in `generators/ic.gen`, then the
file IC runtime path; the Fourier resampling is `fourier_resample` in
`generators/periodic_spectral.py`.

**What it does.** Brings an external triply periodic velocity field onto this case's grid. It
reads the source, re-evaluates it at the target points from its Fourier series - truncating
when the grid is coarser, zero-padding when it is finer - optionally filters it, projects it to
be divergence-free under PICurv's discrete operator, and writes solver units. With `field: Ucont`
(the default) it samples face centres and stages volume fluxes, the representation PICurv evolves;
with `field: Ucat` it stages cell-centre velocities like the synthetic providers.

**When to choose it.** When the starting field exists already: a published DNS or experiment
field, a snapshot from another code, or a field at another resolution. `spectral_random_velocity`
synthesizes a field from a spectrum instead, and suits decay-law and parametric studies where the
initial spectrum is the variable. `mode: file` reads a vector already laid out for this grid,
verbatim; see @ref p33_file_or_import_ssec.

**Parameters it owns.** Required: `source_file` (a workspace input, normally a `reference-field`
import), `format` (`raw`, `npy`, `npz`, or `hdf5`), and `box_length`, the source's three periods
in physical length, which must equal the grid's. `cells` (samples along x, y, z) is required for
`raw` and checked otherwise. `layout`: `fastest_axis` (`x`, default, or `z`) and `components`
(`interleaved`, default, or `separate`); for `raw` also `dtype` (`float32` or `float64`),
`byte_order` (`little`, `big`, `native`) and `header_bytes`; for `npz` and `hdf5`, `datasets` (one
array name, or three for u, v, w). Optional: `sample_offset` (fraction of the source spacing at
which samples sit, 0 for nodes on the box origin and 0.5 for cell centres; default 0), `filter`
(`none`; `sharp` with `cutoff`, a wavenumber; `gaussian` or `box` with `width`, a length),
`projection` (`picurv_discrete`, default; `continuum`; `none`), `remove_mean` (default true), and
`field`. The source values are physical velocities, divided by `velocity_ref`.

**Interactions.** Requires six geometric-periodic faces, one block, and a uniform Cartesian grid;
wall-bounded sources are refused. Restart remains authoritative. The source file's checksum is part
of the asset identity, so changing its bytes builds a new object. `hdf5` needs `h5py`, which PICurv
does not install; without it the provider says how to add it.

**Diagnostics.** The IC summary reports the source, resampled, and projected kinetic energies, the
fraction the projection removed, the divergence before and after, and the energy, dissipation and
`re_lambda` of the field the runtime will hold; for `Ucont` that is the velocity reconstructed from
the staged fluxes (`measurement_state`). The shell spectrum is measured on that same field. A layout
that does not fit the file is refused with the byte or shape count it expected; a `box_length`
that differs from the grid is refused.

**Evidence.** Unit verified - `tests/test_resampled_velocity.py` (`make test-python`).
Analytically verified - `resampled-velocity-hom02-2026-10-03`: a band-limited field coarsened and
refined to round-off (1e-14) at the target cell and face centres, every container and layout
reading the same field, and the HOM02 DNS field (AGARD-AR-345) imported at 64-cubed with every shell
below the cutoff equal to the source's to 7e-16 before projection; the solver's step-1 energy
matched the provider's prediction within one step's decay, at divergence 3.5e-12.

**Limitations.** Periodic boxes only. The projection changes the field near the grid cutoff: with
`Ucat` it removed 2.6% of HOM02's energy at 64-cubed (21% of the k = 32 shell), with `Ucont` 0.075%.
The runtime then holds a second-order reconstruction of what was staged - 87% of the band's energy
for `Ucont` and 78% for `Ucat` on that field - so compare spectra only well below the cutoff.
`projection: continuum` keeps the source's modes but leaves the staged field discretely divergent.

@subsection p33_cap_gen_channel_spectral_velocity_sub channel_spectral_velocity

@anchor p33_cap_gen_channel_spectral_velocity

**Identity.** `source: {generator: channel_spectral_velocity}` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the staged-field runtime path.

**What it does.** A discrete-curl spectral perturbation vanishing at one wall pair, with no mean.

**When to choose it.** Choose for two no-slip walls; use duct_spectral_velocity for four walls.

**Parameters it owns.** The wall-bounded parameter block above applies with one wall_axes entry; initial_spectra is required.

**Interactions.** Single-block 3D Cartesian grid, one no-slip pair, Ucat file staging; the spanwise axis is periodic, and the streamwise axis is periodic or runs from an inlet to an outlet. As a `sum` term it supplies the fluctuations about a source's mean (@ref p33_source).

**Diagnostics.** Summary reports perturbation RMS, component means and discrete divergence; selected plane/line spectra measure the seed.

**Evidence.** Unit verified - `make test-python`. Production exercised in
`examples/turbulent_channel`. Analytically verified - `wall-spectral-ic-2026-09-24`, on a
stretched 32 x 32 x 64 channel when the generator owned a flux-normalized `4t(1-t)` mean, which
the source sum reproduces to 2e-13 at 32 cells: the perturbation RMS is exactly the requested 0.1 with zero component means, the
dummy layers are exact odd reflections so the wall faces carry zero velocity, the
runtime's step-0 flux field has a cell divergence of 1.4e-17 against a flux scale of
4.1e-2, and the wall-parallel perturbation grows linearly from the wall (first-to-second
cell RMS ratio 2.42, the discrete curl's value for a potential ~ y^2).

**Limitations.** A startup construction: amplitude and resolution do not establish
sustained turbulence, and the spectrum shapes the vector potential rather than the
velocity. Measured on a reduced grid on one rank, not at the shipped production size. The
measurement covers a periodic streamwise axis; the inlet-outlet axis is exercised only by
`tests/test_velocity_sources.py` and a one-period demonstration run.

@subsection p33_cap_gen_duct_spectral_velocity_sub duct_spectral_velocity

@anchor p33_cap_gen_duct_spectral_velocity

**Identity.** `source: {generator: duct_spectral_velocity}` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the staged-field runtime path.

**What it does.** A discrete-curl spectral perturbation vanishing at two wall pairs, with no mean.

**When to choose it.** Choose for rectangular/square ducts with four walls; channel_spectral_velocity has only one wall axis.

**Parameters it owns.** The wall-bounded parameter block above applies with two wall_axes entries; initial_spectra selects streamwise lines.

**Interactions.** Single-block 3D Cartesian grid, two no-slip pairs, Ucat file staging; the streamwise axis is periodic or runs from an inlet to an outlet. As a `sum` term it supplies the fluctuations about a source's mean (@ref p33_source).

**Diagnostics.** Summary reports perturbation RMS and discrete divergence; actual selected line spectra report sample energy and Parseval residual.

**Evidence.** Unit verified - `make test-python`. Production exercised in
`examples/periodic_test/driven_duct`. Analytically verified - `wall-spectral-ic-2026-09-24`,
on a doubly stretched 32 x 32 x 64 duct when the generator owned a flux-normalized product of
`4t(1-t)` parabolas, which the source sum reproduces to 4e-13 at 32 cells: the perturbation RMS is exactly the requested 0.1, all four walls carry zero
velocity through exact odd reflection, the runtime's step-0 flux field has a cell
divergence of 1.0e-17 against a flux scale of 2.7e-2, and near-wall growth is linear on
both wall axes (ratios 2.42 and 2.44).

**Limitations.** A startup construction establishes no sustained turbulence. Circular pipes and
immersed boundaries are outside this provider. Measured on a reduced grid on one rank. The
inlet-outlet streamwise axis is not covered by that measurement.

@section p33_stability_mode_sec Linear-stability mode seeds

`stability_mode_velocity` builds a parallel base flow `U(y)` along `streamwise_axis`, varying
across `cross_axis`, plus the most unstable inviscid (Rayleigh) eigenmodes at chosen streamwise
harmonics of the box, and optional broadband noise. It suits temporal free shear flows: a tanh
profile gives a mixing layer, `sech^2` a jet (sinuous and varicose modes), a velocity deficit a
wake. The third axis, `normal`, must be periodic like the streamwise one; the cross-stream faces
may be any non-periodic type.

```yaml
mode: generated
generator: stability_mode_velocity
params:
  field: Ucat
  streamwise_axis: k
  cross_axis: j
  base_profile: "0.5*tanh(2*y/0.05)"     # physical velocity of the physical cross coordinate
  modes:
    - {harmonic: 2, amplitude: 0.01}     # fundamental of a box two wavelengths long
    - {harmonic: 1, amplitude: 0.01, phase: 0.0}
  noise: {rms: 0.02, seed: 1, planar: true, spectrum: {type: white}, wall_modes: 8}
  # velocity_difference: 1.0             # default |U(y_max) - U(y_min)|
  # eigen_points: 401
```

For harmonic `n` the wavenumber is `alpha = 2 pi n / L_s`, `L_s` the streamwise length. The
eigenfunction solves `(U - c)(phi'' - alpha^2 phi) - U'' phi = 0` with `phi = 0` at both
cross-stream ends, by second-order differences on `eigen_points` uniform points, and the mode
with the largest growth rate `alpha c_i` is kept. Each mode is the discrete curl of
`phi e^{i(alpha s + phase)}` along the normal axis, so the field is in the discrete divergence
nullspace and has no normal velocity, and it is scaled so its streamwise fluctuation peaks at
`amplitude * velocity_difference`. `noise.rms` is the volume RMS of the noise as a fraction of the
velocity difference; it is the wall-bounded seed's divergence-free perturbation with the
cross-stream ends as its walls (zero velocity there). `noise.planar: true` keeps it in the
streamwise and cross-stream plane and uniform along the normal axis, which 2D mode requires.

The base profile may use only the cross-stream coordinate (`x`, `y` or `z` for `i`, `j`, `k`)
and is evaluated in physical units, then divided by the reference velocity. Streamwise and normal
axes must be uniform and axis-aligned with at least four cells; the cross-stream axis may stretch.

@subsection p33_cap_gen_stability_mode_velocity_sub stability_mode_velocity

@anchor p33_cap_gen_stability_mode_velocity

**Identity.** `source: {generator: stability_mode_velocity}` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the staged-field runtime path.

**What it does.** A parallel base flow plus Rayleigh eigenmodes at chosen streamwise harmonics and optional divergence-free noise (@ref p33_stability_mode_sec).

**When to choose it.** Temporal mixing layers, jets and wakes forced at their natural instability, and linear growth-rate checks of the solver.

**Parameters it owns.** `streamwise_axis`, `cross_axis`, `base_profile`, `modes` (`harmonic`, `amplitude`, `phase`), `velocity_difference`, `eigen_points`, `noise`, and the optional `initial_spectra` over the periodic axes.

**Interactions.** Single-block axis-aligned Cartesian grid; geometric periodic streamwise and normal axes; non-periodic cross-stream faces. In 2D mode the normal axis must be `i` and `noise.planar` true. Refuses a harmonic at or above a third of the streamwise cells and a profile with no unstable mode at a requested harmonic.

**Diagnostics.** The initial-condition summary reports, per mode, the wavenumber, complex phase speed, growth rate (solver and physical units) and the streamwise amplitude left at the cross-stream ends; a growth-rate scan over the resolvable harmonics and the fastest-growing one, for sizing the box; and the discrete divergence.

**Evidence.** Unit verified - `make test-python` checks Michalke's tanh result (most unstable at
`alpha = 0.4449`, growth rate 0.0949 for `U = (1 + tanh y)/2`), the neutral wavenumbers of the tanh
layer (1) and the Bickley jet's sinuous mode (2), a divergence below 1e-12 on uniform and stretched
cross-stream grids, zero normal velocity with planar noise, the amplitude normalization, the units,
the refusals and the conductor routing. Analytically verified - `stability-mode-growth-2026-10-21`:
a 2D-mode PICurv solve of the seeded tanh layer grows at 0.0882 and 0.0927 on 64 and 128 cross-stream
cells against the Rayleigh rate 0.0948, converging toward it (Richardson estimate 0.0948).

**Limitations.** Temporal and parallel only: no spatial (complex-wavenumber) modes, no inflow
forcing. Inviscid, so the modes do not satisfy no-slip at cross-stream walls; the summary reports
what is left there, negligible when the walls are far. Streamwise modes only: no spanwise or
oblique wavevectors. Axis-aligned Cartesian grids.
