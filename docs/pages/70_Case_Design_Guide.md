@page 70_Case_Design_Guide Case Design Guide

@anchor _Case_Design_Guide

@pagemeta{Explanation, Users building a new case, Decision guidance - not an executable recipe}

Building a case from scratch is a different journey from adapting an example, and
"read all four YAML references" is not an answer to it. This page walks the
**decisions** in the order they constrain each other, and links out for the exact
syntax at each step.

@note This is decision guidance, not an executable recipe. It does not end with a
working case, because the choices below are yours to make. For a path that does end
with a running simulation, use **@subpage 41_Getting_Started_Index** and then adapt
what it produced.

@warning Adapting a shipped example is almost always cheaper and safer. Start at
**@subpage 65_Example_Catalog** and only come back here when nothing shipped is close
enough to your problem.

@tableofcontents

@section p70_order_sec 1. The Order That Matters

Decisions constrain each other in one direction, so taking them out of order means
redoing work:

```text
physical problem
  -> geometry and grid
    -> boundary topology
      -> scaling and properties
        -> numerical method
          -> stability and convergence
            -> observability
              -> outputs and analysis
```

Geometry constrains which boundaries are possible; boundaries constrain which
scalings make sense; scaling constrains which timestep is stable. Choosing a solver
before you know your grid is choosing in the dark.

At each step: build the **smallest valid configuration first**, validate it, and only
then add complexity.

@section p70_geometry_sec 2. Geometry and Grid

**Minimum viable choice.** A programmatically generated box (`grid.mode:
programmatic_c`). It needs no external files and validates fastest.

**Alternatives.** A file-based curvilinear grid when the geometry is not a box, or
`grid_gen` when you want the shipped generator to build it. See
**@subpage 48_Grid_Generator_Guide**.

**Constraints to check now, not later.**

- Each axis should remain coarsenable for the multigrid depth you intend — see
  @ref p25_mg_sec if the Poisson solve later refuses your level count.
- Periodic directions impose a translation contract on the grid; a grid that is not
  truly periodic is rejected, not approximated. See @ref p54_grid_sec.

**Syntax.** @ref p07_grid_sec.

@section p70_boundaries_sec 3. Boundary Topology

Decide the *shape* of the problem before the numbers: which faces are walls, which
admit flow, which are homogeneous.

**Choose the type first**, then the handler — not every handler is legal for every
type. The compatibility matrix and the entries for both are at
@ref p44_types_sec and @ref p44_entries_sec.

**Common topologies.**

| Problem | Topology |
|---|---|
| Through-flow channel or duct | `inlet` + `outlet` + `wall` |
| Statistically homogeneous, decaying | `periodic` on the homogeneous axes with @ref p44_cap_geometric "geometric" |
| Sustained wall-bounded turbulence | `periodic` streamwise with a driven handler, `wall` on the others |

**The trap.** A periodic wall-bounded flow decays under wall friction unless it is
driven. If you want sustained turbulence, you need @ref p44_cap_constant_flux "constant_flux"
or @ref p44_cap_initial_flux "initial_flux", on **both** faces of the pair.

@section p70_scaling_sec 4. Scaling and Properties

Values in YAML are **physical**; `picurv` non-dimensionalizes them before the runtime
sees them. Getting this wrong produces plausible-looking wrong answers rather than
errors, which makes it the highest-risk step on this page.

- Set the reference scales deliberately; do not inherit them by accident.
- `target_flux` is a **volumetric** flux, not a velocity.
- The project-wide scaling conventions are in **@subpage 19_Nondimensionalization**;
  @ref p44_nondim_sec covers the boundary-condition specifics.

@section p70_method_sec 5. Numerical Method

**Minimum viable choice.** `Dual Time Picard Jameson RK`, the default production
solver. It tolerates far larger timesteps than the explicit path.

**Alternatives and why they exist.** Compare against siblings at @ref p08_entries_sec
before switching. `Newton Krylov` is the right escalation when pseudo-time stalls;
`Explicit RK4` currently has **no positive-path verification** (see
**@subpage 66_Evidence_Matrix**) and should not be a first choice.

**Turbulence model.** Use none.

@warning Both LES models are **known-defective** and unsafe for scientific use:
@ref p07_cap_les_constant_smagorinsky "constant_smagorinsky" is inert on fresh runs, and
@ref p07_cap_les_dynamic_smagorinsky "dynamic_smagorinsky" has no dynamic content. There is
currently no usable subgrid model in this tree. If your problem requires LES, that is a
blocker to resolve before designing the case, not a setting to choose.

@section p70_stability_sec 6. Stability and Convergence

1. Start with a conservative timestep and loosen it once the case runs.
2. Leave solver tolerances at their defaults for the first successful run. Tuning
   tolerances before you have a baseline means you cannot tell what your change did.
3. Run a handful of steps and confirm convergence behavior before committing to a
   long run. @ref p67_noconverge_sec covers what to look at.

@section p70_observability_sec 7. Observability

Configure this **before** the first expensive run, not after it fails:

- Logging cadence frequent enough to see divergence early, sparse enough not to
  dominate runtime.
- Checkpoint cadence that lets you restart without losing meaningful work.
- Field statistics if you need averages — starting them late costs a rerun. See
  **@subpage 58_Field_Statistics**.

**Syntax.** @ref p09_io_sec.

@section p70_outputs_sec 8. Outputs and Analysis

Decide what you will actually compute from the run, because it determines what must
be written. **@subpage 10_Post_Processing_Reference** covers the pipelines;
**@subpage 04_Visualization_Tutorial** covers getting pictures out.

@section p70_checklist_sec 9. Before The First Expensive Run

- [ ] `picurv validate` passes on all four roles.
- [ ] `--dry-run` output shows the artifacts and commands you expect.
- [ ] The case ran for a few steps locally and converged.
- [ ] Boundary handlers in the startup banner are the ones you intended.
- [ ] Units checked against @ref p44_nondim_sec.
- [ ] Output cadence will not exhaust your storage quota.
- [ ] Checkpoint cadence lets you resume without losing meaningful work.
- [ ] You know which capabilities you are relying on that are unverified
      (**@subpage 66_Evidence_Matrix**).

@section p70_endpoint_sec 10. Where This Ends

You finish this page with decisions, not with a case. To turn them into one:

1. `./bin/picurv init <closest example> --dest my_case` — start from the nearest
   shipped case rather than an empty file.
2. Apply your decisions one at a time, running `picurv validate` between each.
3. Work through the checklist in section 9.
4. Run a few steps locally before committing to anything expensive.

@section p70_related_sec 11. Related Documentation

- **@subpage 49_Workflow_Recipes_and_Config_Cookbook** — recombining profiles
- **@subpage 52_Run_Lifecycle_Guide** — running, restarting, archiving
- **@subpage 67_Troubleshooting** — when it does not work
