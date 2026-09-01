@page 72_LES_Turbulence_Closure LES Turbulence Closure

@anchor _LES_Turbulence_Closure

Large-eddy simulation resolves the motions the grid can carry and models the rest. The
modelling is concentrated in a single coefficient, and the whole of this page is about
where that coefficient comes from, what the code actually computes, and which choices
are yours.

The YAML keys are specified in @ref 07_Case_Reference; this page explains the physics
behind them. The runtime call order is in @ref 46_C_Runtime_Execution_Map.

@tableofcontents

@section p72_problem_sec 1. The Closure Problem

Filtering the incompressible Navier-Stokes equations at the grid scale leaves the
nonlinear term with a residue that cannot be evaluated from the filtered field alone:

    tau_ij = (u_i u_j)~ - u~_i u~_j

This is the **subgrid stress**: the momentum transfer done by eddies smaller than a
cell. It involves the unfiltered velocity, which the simulation never had. Closing the
equations means modelling it.

Smagorinsky's closure is the simplest one that works, and it is an eddy-viscosity
ansatz: pretend the unresolved motions act like an extra viscosity that drains energy
from the resolved field.

    tau_ij ~= -2 (Cs Delta)^2 |S~| S~_ij  =  -2 C Delta^2 |S~| S~_ij

with `Delta` the filter width, `S~_ij` the resolved strain rate, and
`|S~| = sqrt(2 S~_ij S~_ij)`. The eddy viscosity is `nu_t = C Delta^2 |S~|`, added to
the molecular viscosity in the momentum equations.

Everything then rests on the single number `C`. Note the notation carefully, because
the code follows it: **`C` is the coefficient that multiplies `Delta^2 |S~|`, which is
`Cs^2` in the classical form, not `Cs`**. The stored field holds `C`. Reported
diagnostics convert to `Cs`, which is the quantity with a value readers recognise.

@section p72_constant_sec 2. The Constant Model, and Why It Is Not Enough

@ref p07_cap_les_constant_smagorinsky "constant_smagorinsky" prescribes `Cs` and is
done. It is cheap, predictable, and correct wherever a single number is appropriate.

It is not appropriate in most places. The right value differs between isotropic
turbulence (Lilly's 0.16-0.17, valid only where the grid cutoff sits in an inertial
range), a wall-bounded flow (where the coefficient must vanish at the wall and no
constant does), a transitional region (where any nonzero value dissipates a flow that
is not yet turbulent), and a laminar pocket (where the model should switch itself off
and cannot).

The dynamic procedure exists to let the flow supply the number.

@section p72_germano_sec 3. The Germano Identity

The insight (Germano et al. 1991) is that although the simulation cannot see below its
grid, it can see between the grid scale and a coarser one. Take the already-filtered
field and filter it again with a wider **test filter** of width `Delta^ = r Delta`, and
three stresses come into view:

- `tau_ij` — the stress at the grid filter. Unknown; this is what is being modelled.
- `T_ij` — the stress at the combined grid-and-test filter. Also unknown.
- `L_ij = (u~_i u~_j)^ - u^_i u^_j` — the **Leonard stress**, computable entirely from
  the resolved field.

The identity connecting them is exact and purely algebraic, with no modelling in it:

    L_ij = T_ij - (tau_ij)^

`L_ij` is the stress carried by the eddies between `Delta` and `Delta^`: the band the
simulation can actually resolve and therefore measure. That is the whole trick.

Now apply the same Smagorinsky ansatz at both levels with the same coefficient:

    tau_ij ~= -2 C Delta^2  |S~| S~_ij
    T_ij   ~= -2 C Delta^^2 |S^| S^_ij

Substituting into the identity gives `L_ij ~= C M_ij` with

    M_ij = -2 Delta^2 ( alpha |S^| S^_ij  -  (|S~| S~_ij)^ ),    alpha = (Delta^/Delta)^2

@warning The two terms of `M_ij` are structurally different and must stay so. The first
filters `S~` and *then* forms the product, so it is a **product of filtered
quantities**. The second forms the product `|S~| S~_ij` pointwise on the fine field and
*then* filters the whole thing, so it is a **filter of a product**. Those are not equal,
and their difference is the entire information content of the dynamic procedure. Written
with both terms as the same expression, `M_ij` collapses to a fixed multiple of
`|S^| S^_ij`, and the result is a well-scaled, dimensionally correct, flow-responsive
coefficient that is **not** the Germano-Lilly coefficient. It cannot be cited as one.
This is exactly the shape the implementation once had, which is why
@ref GermanoModelTensor is one function with one call site and why
`tests/c/test_les.c` pins the two terms apart on a field where they provably differ.

@subsection p72_lilly_ssec 3.1 Lilly's Contraction

`L_ij = C M_ij` is six equations for one unknown. Lilly (1992) resolves the
overdetermination by least squares, minimising the residual over the tensor components:

    C = <L_ij M_ij> / <M_ij M_ij>

The angle brackets are not decoration. The derivation pulls `C` out of a filtering
operation, which is only legitimate if `C` is constant over whatever set the brackets
average. That assumption is the subject of the next section.

`M_ij` is made **deviatoric** before the contraction. The trace of `L_ij` is absorbed
by the pressure in an incompressible flow and carries no information about the
coefficient; projecting `M_ij` trace-free removes it from the contraction exactly,
rather than relying on the discrete divergence being small enough to ignore.

@section p72_averaging_sec 4. Averaging

@ref p07_les_avg_sec selects the set the two contractions are averaged over, and it
matters for three separate reasons.

**The derivation requires it.** Lilly's least squares assumes `C` is constant across
the averaging set. A pointwise coefficient was already inside the test filter used to
build `M_ij`, so the local form is internally inconsistent. Averaging over a region
where `C` genuinely is uniform restores the consistency.

**The denominator vanishes.** Wherever the resolved strain is briefly small, `M:M`
collapses and the pointwise quotient becomes a ratio of two small noisy numbers.
Averaging over a set large enough to contain structure fixes this at the source rather
than guarding it after the fact.

**Backscatter should cancel where it physically does.** Negative `L:M` is real. Averaged
over a homogeneous set, forward and reverse transfer combine the way the physics
intends. Clipped pointwise, only the negative tail is removed, and the mean dissipation
comes out biased upward.

Two implementation properties are worth knowing. The sums are **volume-weighted**, so a
refined region does not dominate an average merely by holding more cells. And the
numerator and denominator are averaged **separately** before dividing: the mean of the
quotients is a different number from the quotient of the means, and it is the second
that the closure defines. `tests/c/test_les.c` separates them explicitly.

@subsection p72_homogeneous_ssec 4.1 Deriving the Directions

`averaging.mode: homogeneous` takes its directions from the case's periodicity unless
they are named explicitly. This is not a shortcut; periodicity is already the
declaration that the flow repeats along an axis, which is the same claim homogeneity
makes.

The closure reads the resolved per-axis flags rather than re-inspecting the boundary
faces. The two cannot disagree: the loader derives the flags from the boundary files
and rejects any case whose opposite faces, or whose blocks, do not agree on
periodicity. Those flags are also what the DMDA itself was built from and what the
shared spatial target consults, so there is one answer to the question rather than
three places that each work it out. It means the two canonical cases need no extra configuration:

- A **triply periodic box** averages over all three directions and yields one
  coefficient per update. This is the arrangement for decaying isotropic turbulence,
  and `Cs(t)` is then a single readable curve rather than a noisy field.
- A **plane channel**, periodic in xi and zeta, averages over those two and retains
  eta, yielding a wall-normal `Cs(y)` profile. The coefficient falls toward the wall on
  its own, which is the behaviour a constant coefficient cannot reproduce.

Naming a direction the case does not declare periodic is accepted with a warning, since
homogeneity is a claim only the user can make.

@section p72_clipping_sec 5. Limiting and Backscatter

A negative `C` means `nu_t < 0`: energy flowing from the unresolved scales back to the
resolved ones. Backscatter is physical, and in a locally averaged dynamic model it
occurs at a large fraction of points in developed turbulence.

It is also destabilising if unbounded, so something has to limit it. The three modes at
@ref p07_les_clip_sec differ in what they limit and at what cost.

The important design point is **where** the bound is applied. The clipping modes bound
the coefficient's sign, which removes backscatter entirely and biases the mean
dissipation. `clipping.mode: none` instead keeps the signed coefficient and bounds the
**total viscosity**:

    nu + nu_t >= min_viscosity_ratio * nu

That constrains the quantity which actually has to stay positive for the momentum
operator to remain well posed, and it leaves the physics alone everywhere the bound is
not active. It is the better guard where backscatter matters, and it pairs naturally
with homogeneous averaging, where the averaged coefficient goes negative only when the
whole homogeneous set is backscattering on balance.

@section p72_widths_sec 6. Filter Widths

Two widths enter, and both are configurable because neither has a universally right
value on a curvilinear grid.

The **grid filter width** `Delta` is derived per cell by @ref p07_les_width_sec. The
default cube root of the cell volume is exact for a cube and progressively optimistic as
a cell is stretched, because a geometric mean is dominated by the short directions. On a
wall-normal channel mesh `max_edge` is the more defensible choice.

The **test-to-grid ratio** `Delta^/Delta` sets `alpha` and defaults to 2.0, giving
`alpha = 4`. The value is exposed because the effective width of the discrete
three-point stencil is not exactly twice the grid spacing: matching second moments for a
top-hat gives `sqrt(6) ~= 2.449`. Varying it is a standard sensitivity study, and the
ingress rejects any value at or below 1, where the two filters coincide and the
procedure has nothing to measure.

@section p72_gradient_sec 7. The Gradient Term and Mixed Models

`gradient_model.enabled` adds the Clark tensor-diffusivity term to the viscous flux. It
is the leading term of a Taylor expansion of the Leonard stress, and it is
structurally different from an eddy viscosity: it is not purely dissipative and it does
not require a coefficient to be measured.

Combined with a Smagorinsky closure this is the classical **mixed model** — the gradient
term supplies the correct stress structure while the eddy viscosity supplies the
dissipation the gradient term lacks. The two are configured independently because they
are independent terms.

@note The gradient term is available under both momentum solver families, reaching each
through the shared residual. Under Newton-Krylov the preconditioner does not model it,
so Krylov iteration counts can rise with it enabled; watch them when you turn it on.

@section p72_diagnostics_sec 8. Reading the Diagnostics

With `diagnostics.enabled`, the solver appends one row per step (or per `cadence` steps)
to `<run.runtime_logs>/les_coefficient.csv`. Three groups of columns answer different questions.

**Is the coefficient right?** `cs_effective` is the whole-domain
`sqrt(<L:M>/<M:M>)`, computed the same way regardless of averaging mode so it stays
comparable across modes. For decaying isotropic turbulence it should settle near
**0.16-0.17**. `cs_mean`, `coefficient_rms`, `coefficient_min`, and `coefficient_max`
describe the spread; under global averaging the spread should be zero by construction.

**How much is the model doing?** `nu_t_mean`, `nu_t_max`, and `nu_t_over_nu_mean` give
the eddy-viscosity level. `k_sgs_mean` is the modelled subgrid kinetic energy from the
Yoshizawa relation `k_sgs = 2 C_I Delta^2 |S~|^2`. Compared against the resolved
turbulent kinetic energy from @ref 58_Field_Statistics it measures how much of the
energy the grid resolves — Pope's criterion, which asks for the resolved share to exceed
roughly 80%. The two halves come from different pipelines on purpose: resolved turbulent
energy needs the mean subtracted correctly, which the statistics windows do and the LES
module has no business doing.

**What is the limiting costing?** `backscatter_fraction` and `limited_fraction` are the
volume fractions that were negative, and that the clip modified, **before** clipping.
No stored field preserves that state, which is why these columns exist. A `max_cs` that
never binds reads near zero.

@warning These are instantaneous volume statistics, not time averages. Window-averaged,
spatially resolved statistics of `Nu_t` and `CS` come from @ref 58_Field_Statistics
instead, which can target them like any other catalogued field.

@section p72_extending_sec 9. Extending the Closure

The module is built from kernels rather than from one procedure, so a new closure reuses
the parts it needs instead of restating them. The public ones are declared in
`include/les.h`:

- symmetric tensor algebra — @ref SymTensorDeviator, @ref SymTensorContract,
  @ref SymTensorSelfOuter, @ref SymTensorCombine — which knows nothing about turbulence;
- @ref StrainRateFromGradients and @ref ComputeCellFilterWidth, shared by every
  eddy-viscosity model;
- @ref LeonardStress and @ref GermanoModelTensor, the identity itself;
- @ref PicurvSpatialRatioAverage — declared in `include/statistics_target.h`, not here,
  because averaging fields over a spatial domain is not an LES concern; it turns a local
  coefficient into a plane- or volume-averaged one and needs no separate code path per
  mode;
- @ref ClipModelCoefficient and @ref EddyViscosityFromCoefficient.

**Another eddy-viscosity closure** (WALE, Vreman, the sigma model) is one function
consuming the strain rate and the filter width and producing `C`, plus the selector
wiring in @ref 50_Modular_Selector_Extension_Guide. Dispatch is an enum and a branch in
`FlowSolver`, the same pattern the momentum solvers use; there is no registry to
register with.

**A different averaging strategy**, such as the Lagrangian dynamic model of Meneveau,
Lund and Cabot (1996), replaces @ref PicurvSpatialRatioAverage with a relaxation
transported along pathlines. It is not implemented and no placeholder exists for it.

**A transported-PDF closure** such as a velocity filtered density function is a
different kind of model: it evolves particles rather than adding momentum diffusion, and
the Lagrangian machinery it needs is already in the tree. What it needs from here is
`nu_t`, `Delta`, and the subgrid kinetic energy — the first two are fields and the third
is @ref SubgridKineticEnergy. No accessor is provided in advance, because an interface
with no caller is a guess.

@section p72_status_sec 10. Status and What Would Change It

Both models are **experimental**. The formulation is implemented and unit-tested:
`tests/c/test_les.c` checks the model tensor against its closed form on constant strain,
pins the filtered product apart from the product of filtered factors, verifies the
procedure returns exactly zero on uniform flow, and checks that global averaging gives
one coefficient per block; `tests/c/test_mpi_kernels.c` checks that the averaging
reduction is independent of the decomposition.

None of that validates the coefficient's **magnitude**. The run that would is
`examples/decaying_isotropic_turbulence` with `averaging.mode: homogeneous`, where
`cs_effective` should settle near 0.16-0.17 and the energy spectrum should show no
pile-up at the grid cutoff. Until that is run and recorded, coefficient magnitudes from
this tree are uncharacterized.

@section p72_references_sec 11. References

- Smagorinsky, J. (1963). General circulation experiments with the primitive equations.
  *Monthly Weather Review* 91(3), 99-164.
- Clark, R. A., Ferziger, J. H., Reynolds, W. C. (1979). Evaluation of subgrid-scale
  models using an accurately simulated turbulent flow. *JFM* 91(1), 1-16.
- Yoshizawa, A. (1986). Statistical theory for compressible turbulent shear flows.
  *Physics of Fluids* 29(7), 2152-2164.
- Germano, M., Piomelli, U., Moin, P., Cabot, W. H. (1991). A dynamic subgrid-scale eddy
  viscosity model. *Physics of Fluids A* 3(7), 1760-1765.
- Lilly, D. K. (1992). A proposed modification of the Germano subgrid-scale closure
  method. *Physics of Fluids A* 4(3), 633-635.
- Meneveau, C., Lund, T. S., Cabot, W. H. (1996). A Lagrangian dynamic subgrid-scale
  model of turbulence. *JFM* 319, 353-385.

@section p72_related_sec 12. Related Pages

- @subpage 07_Case_Reference — the YAML keys and every selector value
- @subpage 21_Methods_Overview — where this closure sits among the governing methods
- @subpage 58_Field_Statistics — window-averaged statistics of the model's output fields
- @subpage 50_Modular_Selector_Extension_Guide — the checklist for adding a selector
- @subpage 46_C_Runtime_Execution_Map — when the closure runs within a timestep
