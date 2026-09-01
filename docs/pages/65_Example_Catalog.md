@page 65_Example_Catalog Example Catalog

@anchor _Example_Catalog

@pagemeta{Reference, All readers, Cost column populated as runs are measured}

Every example directory shipped under `examples/`, what it demonstrates, and what
kind of evidence it provides. Start here rather than browsing the directory — a
listing tells you a case exists, not whether it is a tutorial, a verification test,
or a research campaign.

`examples/master_template` is intentionally excluded: it is a commented reference of
every available key, not a runnable case.

@tableofcontents

@section p65_reading_sec 1. How to Read This Catalog

**Kind** says what the example is *for*:

- **Tutorial** — meant to be followed end to end while learning. Cheap, robust, and
  documented step by step.
- **Verification** — isolates one mechanism and checks it against an exact or
  analytic answer. Small and fast; a failure means the code is wrong, not that the
  physics is hard.
- **Benchmark** — reproduces a published or in-tree reference result.
- **Characterization** — instruments behavior to produce measurements rather than to
  pass or fail.

**Evidence** uses the facet vocabulary from **@subpage 62_Capability_Status_Vocabulary**.

**Cost** is deliberately blank where it has not been measured on stated hardware. An
unmeasured cost is recorded as unknown rather than estimated, because a wrong runtime
estimate is worse than none when someone is sizing a cluster request.

@section p65_start_sec 2. Start Here

| Example | Kind | Demonstrates | Execution | Cost |
|---|---|---|---|---|
| `flat_channel` | Tutorial | Laminar channel on a programmatically generated grid; the canonical first run | serial / MPI | not yet measured |
| `bent_channel` | Tutorial | The same workflow driven by a **file-based** curvilinear grid | serial / MPI | not yet measured |

`flat_channel` is the case **@subpage 41_Getting_Started_Index** walks you through.
Run it before anything else on this page.

@section p65_verify_sec 3. Particle Verification Family

These isolate one term of the particle update each, which is why they are the right
place to look when particle behavior is suspect: each failure points at one mechanism.

| Example | Kind | Isolates | Evidence |
|---|---|---|---|
| `drift_uniform_flow` | Verification | Deterministic advection by interpolated velocity | Designed for analytical comparison; threshold not gated |
| `drift_diffusivity_gradient` | Verification | The `grad D * dt` drift term | Designed for analytical comparison; threshold not gated |
| `brownian_motion` | Verification | The stochastic displacement, against a diffusion solution | Designed for analytical comparison; threshold not gated |
| `interpolation_test` | Verification | Grid-to-particle interpolation accuracy | Designed for analytical comparison; threshold not gated |
| `scatter_verification` | Verification | The particle-to-grid scatter path | Designed for analytical comparison; threshold not gated |
| `search_robustness` | Characterization | Walking-search and migration behavior under stress | Benchmark characterized |

Together these address every term in the displacement equation documented at
**@subpage 34_Particle_Model_Overview**.

@note These cases are *designed* to compare against exact or analytic answers, but a
numerical acceptance threshold has not been run and gated as part of the current
documentation work. That is why **@subpage 66_Evidence_Matrix** records no
`analytical` facet for them: design intent is not evidence.

@section p65_flow_sec 4. Flow and Turbulence

| Example | Kind | Demonstrates | Status |
|---|---|---|---|
| `decaying_isotropic_turbulence` | Benchmark | LES decay in a triply periodic box, 64³ cells | Selects an experimental model — see below |
| `periodic_test/driven_channel` | Benchmark | Driven periodic channel; DNS and LES variants | See caveat below |
| `periodic_test/driven_duct` | Benchmark | Square duct with secondary flow of the second kind | See caveat below |

@note **The shipped LES configuration is experimental.**
`decaying_isotropic_turbulence` selects `model: dynamic_smagorinsky` with
`averaging.mode: homogeneous`, which on a triply periodic box gives one coefficient
for the whole domain and writes `Cs(t)` to `<run.runtime_logs>/les_coefficient.csv`. The
formulation is unit-tested but the coefficient magnitude has not been validated: this
case is the run that would settle it, and `Cs(t)` is expected to settle near
0.16-0.17. Until that is recorded, treat the magnitude as uncharacterized.

@warning The `periodic_test` campaigns carry an open status: pseudo-time momentum
convergence on periodic wall-bounded flow was observed to stall on 2026-08-24 and
**requires re-characterization at current `HEAD`** after the convergence-criterion
change. Read @ref p54_driven_limits_sub before planning a campaign around them. Use
the periodic boundary handlers documented at @ref p44_cap_geometric "geometric" and
the driven handlers; note the LES caveat above applies to their LES variants.

@section p65_authoring_sec 5. Using an Example as a Starting Point

The intended workflow is to copy and adapt rather than author from scratch:

```bash
./bin/picurv init flat_channel --dest my_case
```

Then change one thing at a time, re-validating between changes. **@subpage 49_Workflow_Recipes_and_Config_Cookbook** covers profile recombination, and
`examples/master_template` documents every available key with commentary — it is a
reference, not a runnable case.

@section p65_related_sec 6. Related Documentation

- **@subpage 41_Getting_Started_Index** — the guided first run
- **@subpage 40_Testing_and_Quality_Guide** — which examples CI exercises
- **@subpage 62_Capability_Status_Vocabulary** — what the evidence words mean
