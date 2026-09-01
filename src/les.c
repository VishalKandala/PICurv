/**
 * @file les.c
 * @brief Implements the Large Eddy Simulation (LES) subgrid-scale closure.
 *
 * The closure is built from small kernels that each do one thing, so the physics is
 * readable at the call sites and testable without a running solve:
 *
 * - symmetric-tensor algebra (`SymTensor*`), which knows nothing about turbulence;
 * - strain rate and filter width, shared by every eddy-viscosity model;
 * - the Germano identity pieces, `LeonardStress()` and `GermanoModelTensor()`;
 * - the coefficient averaging, which is `PicurvSpatialRatioAverage()` from the shared
 *   spatial-target module: it turns the local dynamic model into the plane- or
 *   volume-averaged one without a separate code path per mode;
 * - limiting and eddy-viscosity assembly.
 *
 * Two routines drive them. `ComputeSmagorinskyConstant()` runs the Germano-Lilly
 * procedure for the dynamic model, and `ComputeEddyViscosityLES()` turns whichever
 * coefficient is in force into the eddy viscosity the momentum equations add to the
 * molecular one. `LogLESDiagnostics()` reports the coefficient's time history.
 *
 * The coefficient stored in `UserCtx::CS` is `C`, the factor multiplying
 * `Delta^2 |S|`, which is `Cs^2` in the classical notation and is signed when
 * backscatter is admitted. Conversion to the familiar `Cs` happens only where a
 * human reads the number.
 */

#include "les.h"
#include "statistics_target.h"

// A small constant to prevent division by zero in sensitive calculations.
static const double LES_EPSILON = 1.0e-12;

/*================================================================================*
 *                       INTERNAL (STATIC) FUNCTIONS                              *
 *================================================================================*/

/**
 * @brief Synchronizes the completed Smagorinsky coefficient field.
 */
static PetscErrorCode FinalizeSmagorinskyConstantField(UserCtx *user)
{
    const FieldId fields[] = {FIELD_ID_CS};

    PetscFunctionBeginUser;
    PetscCall(SynchronizePeriodicCellFields(user, 1, fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_CS));
    PetscFunctionReturn(0);
}

/**
 * @brief Resolves the iteration domain the closure computes over.
 *
 * Delegates to the shared spatial target so the closure, the statistics pipeline, and
 * anything else looping over cell-centred data agree on which indices are physical.
 * The domain excludes PETSc halo storage and the solver's boundary, dummy, and
 * duplicate-periodic indices; the last of those matters here, because counting a
 * duplicate plane in a spatial average would weight a real cell twice.
 *
 * The plan is resolved for `CS` as a representative cell-centred field. The bounds
 * depend on the layout and the block's periodicity, not on which field is being
 * iterated, so the scratch vectors this module allocates share them.
 */
static PetscErrorCode LESResolveDomain(UserCtx *user, SpatialTargetPlan *plan)
{
    PetscFunctionBeginUser;
    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS, PICURV_STATISTICS_MASK_FLUID, plan));
    PetscFunctionReturn(0);
}

/**
 * @brief Reports the block's per-axis periodicity in `(xi, eta, zeta)` order.
 *
 * Reads the resolved global flags rather than the boundary face configuration. The two
 * cannot disagree: the flags are derived from the boundary files, and the loader
 * rejects a case whose opposite faces or whose blocks do not agree on periodicity. The
 * flags are also what the DMDA itself was built from, and what the shared spatial
 * target consults, so taking them here keeps one source of truth.
 */
static void LESPeriodicAxes(const SimCtx *simCtx, PetscBool periodic[3])
{
    periodic[0] = (PetscBool)(simCtx->i_periodic != 0);
    periodic[1] = (PetscBool)(simCtx->j_periodic != 0);
    periodic[2] = (PetscBool)(simCtx->k_periodic != 0);
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorSelfOuter"
/**
 * @brief Implementation of \ref SymTensorSelfOuter().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorSelfOuter()
 */
SymTensor SymTensorSelfOuter(Cmpnts v)
{
    SymTensor t;

    t.xx = v.x * v.x;
    t.xy = v.x * v.y;
    t.xz = v.x * v.z;
    t.yy = v.y * v.y;
    t.yz = v.y * v.z;
    t.zz = v.z * v.z;
    return t;
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorCombine"
/**
 * @brief Implementation of \ref SymTensorCombine().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorCombine()
 */
SymTensor SymTensorCombine(PetscReal a, SymTensor x, PetscReal b, SymTensor y)
{
    SymTensor t;

    t.xx = a * x.xx + b * y.xx;
    t.xy = a * x.xy + b * y.xy;
    t.xz = a * x.xz + b * y.xz;
    t.yy = a * x.yy + b * y.yy;
    t.yz = a * x.yz + b * y.yz;
    t.zz = a * x.zz + b * y.zz;
    return t;
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorTrace"
/**
 * @brief Implementation of \ref SymTensorTrace().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorTrace()
 */
PetscReal SymTensorTrace(SymTensor t)
{
    return t.xx + t.yy + t.zz;
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorDeviator"
/**
 * @brief Implementation of \ref SymTensorDeviator().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorDeviator()
 */
SymTensor SymTensorDeviator(SymTensor t)
{
    const PetscReal isotropic = SymTensorTrace(t) / 3.0;

    t.xx -= isotropic;
    t.yy -= isotropic;
    t.zz -= isotropic;
    return t;
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorContract"
/**
 * @brief Implementation of \ref SymTensorContract().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorContract()
 */
PetscReal SymTensorContract(SymTensor a, SymTensor b)
{
    // The three off-diagonal components each stand for two entries of the full tensor.
    return a.xx * b.xx + a.yy * b.yy + a.zz * b.zz +
           2.0 * (a.xy * b.xy + a.xz * b.xz + a.yz * b.yz);
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorTimesVector"
/**
 * @brief Implementation of \ref SymTensorTimesVector().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorTimesVector()
 */
Cmpnts SymTensorTimesVector(SymTensor t, Cmpnts v)
{
    Cmpnts result;

    result.x = t.xx * v.x + t.xy * v.y + t.xz * v.z;
    result.y = t.xy * v.x + t.yy * v.y + t.yz * v.z;
    result.z = t.xz * v.x + t.yz * v.y + t.zz * v.z;
    return result;
}

#undef __FUNCT__
#define __FUNCT__ "SymTensorNormSq"
/**
 * @brief Implementation of \ref SymTensorNormSq().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SymTensorNormSq()
 */
PetscReal SymTensorNormSq(SymTensor t)
{
    return SymTensorContract(t, t);
}

/*================================================================================*
 *                        STRAIN RATE AND FILTER WIDTH                            *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "StrainRateFromGradients"
/**
 * @brief Implementation of \ref StrainRateFromGradients().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see StrainRateFromGradients()
 */
PetscErrorCode StrainRateFromGradients(Cmpnts dudx, Cmpnts dvdx, Cmpnts dwdx,
                                       SymTensor *strain, PetscReal *magnitude)
{
    SymTensor s;

    PetscFunctionBeginUser;

    s.xx = dudx.x;
    s.xy = 0.5 * (dudx.y + dvdx.x);
    s.xz = 0.5 * (dudx.z + dwdx.x);
    s.yy = dvdx.y;
    s.yz = 0.5 * (dvdx.z + dwdx.y);
    s.zz = dwdx.z;

    if (strain) *strain = s;
    // |S| = sqrt(2 S_ij S_ij), with the off-diagonal doubling handled by the contraction.
    if (magnitude) *magnitude = PetscSqrtReal(2.0 * SymTensorNormSq(s));

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeCellFilterWidth"
/**
 * @brief Implementation of \ref ComputeCellFilterWidth().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see ComputeCellFilterWidth()
 */
PetscErrorCode ComputeCellFilterWidth(LESFilterWidthModel model, PetscReal aj,
                                      Cmpnts csi, Cmpnts eta, Cmpnts zet, PetscReal *delta)
{
    double dx = 0.0, dy = 0.0, dz = 0.0;

    PetscFunctionBeginUser;
    PetscCheck(delta != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Filter width destination cannot be NULL.");
    PetscCheck(aj > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Cell Jacobian must be positive to define a filter width; received %g.", (double)aj);

    if (model == LES_FILTER_WIDTH_CUBE_ROOT_VOLUME) {
        // Cell volume alone: exact for a cube, and increasingly optimistic as the
        // cell is stretched, because it averages the three extents geometrically.
        *delta = PetscPowReal(1.0 / aj, 1.0 / 3.0);
        PetscFunctionReturn(0);
    }

    PetscCall(ComputeCellCharacteristicLengthScale(aj, csi, eta, zet, &dx, &dy, &dz));

    switch (model) {
    case LES_FILTER_WIDTH_MAX_EDGE:
        // The largest unresolved scale the cell can carry; conservative on stretched grids.
        *delta = PetscMax(dx, PetscMax(dy, dz));
        break;
    case LES_FILTER_WIDTH_GEOMETRIC_MEAN:
    default:
        *delta = PetscPowReal(PetscMax(dx * dy * dz, LES_EPSILON), 1.0 / 3.0);
        break;
    }

    PetscFunctionReturn(0);
}

/*================================================================================*
 *                        GERMANO IDENTITY KERNELS                                *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "LeonardStress"
/**
 * @brief Implementation of \ref LeonardStress().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see LeonardStress()
 */
SymTensor LeonardStress(Cmpnts velocity_filtered, SymTensor velocity_product_filtered)
{
    // L_ij = (u_i u_j)^ - u^_i u^_j : the stress carried between the two filter widths.
    return SymTensorCombine(1.0, velocity_product_filtered,
                            -1.0, SymTensorSelfOuter(velocity_filtered));
}

#undef __FUNCT__
#define __FUNCT__ "GermanoModelTensor"
/**
 * @brief Implementation of \ref GermanoModelTensor().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see GermanoModelTensor()
 */
SymTensor GermanoModelTensor(PetscReal delta, PetscReal alpha,
                             PetscReal strain_magnitude_filtered,
                             SymTensor strain_filtered,
                             SymTensor strain_product_filtered)
{
    // Term one: the model evaluated on the test-filtered strain. Both factors were
    // filtered separately, then multiplied.
    const SymTensor test_scale = SymTensorCombine(alpha * strain_magnitude_filtered,
                                                  strain_filtered, 0.0, strain_filtered);

    // Term two: the grid-scale model formed pointwise and then test-filtered. This is
    // the filter of a product, and it is not equal to the product of the filtered
    // factors above. That inequality carries the whole information content of the
    // dynamic procedure.
    const SymTensor grid_scale = strain_product_filtered;

    const SymTensor model = SymTensorCombine(-2.0 * delta * delta, test_scale,
                                             2.0 * delta * delta, grid_scale);

    // Trace-free, so the Leonard stress's trace cannot enter the contraction through
    // whatever discrete divergence the velocity field carries.
    return SymTensorDeviator(model);
}

/*================================================================================*
 *                     COEFFICIENT AVERAGING AND LIMITING                         *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "ResolveLESAveragingDirections"
/**
 * @brief Implementation of \ref ResolveLESAveragingDirections().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see ResolveLESAveragingDirections()
 */
PetscErrorCode ResolveLESAveragingDirections(UserCtx *user, PetscBool direction[3])
{
    const LESConfig *config = &user->simCtx->les_config;
    PetscInt         selected = 0;
    PetscBool        periodic[3];

    PetscFunctionBeginUser;
    PetscCheck(direction != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Averaging direction destination cannot be NULL.");
    LESPeriodicAxes(user->simCtx, periodic);

    for (PetscInt axis = 0; axis < 3; axis++) direction[axis] = PETSC_FALSE;

    switch (config->averaging_mode) {
    case LES_AVERAGING_LOCAL:
        // No averaging: the coefficient is whatever the single cell's contractions give.
        PetscFunctionReturn(0);

    case LES_AVERAGING_GLOBAL:
        for (PetscInt axis = 0; axis < 3; axis++) direction[axis] = PETSC_TRUE;
        PetscFunctionReturn(0);

    case LES_AVERAGING_HOMOGENEOUS:
    default:
        for (PetscInt axis = 0; axis < 3; axis++) {
            if (config->averaging_direction[axis]) {
                direction[axis] = PETSC_TRUE;
                selected++;
            }
        }
        // Naming no direction asks for the periodic axes, which is where the flow is
        // homogeneous in every case the solver can currently express.
        if (selected == 0) {
            for (PetscInt axis = 0; axis < 3; axis++) {
                if (periodic[axis]) {
                    direction[axis] = PETSC_TRUE;
                    selected++;
                }
            }
        }
        PetscCheck(selected > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                   "Homogeneous LES averaging was requested, but this block declares no "
                   "periodic axis and the configuration names no direction. Name the "
                   "homogeneous directions explicitly or select local averaging.");
        PetscFunctionReturn(0);
    }
}

#undef __FUNCT__
#define __FUNCT__ "ClipModelCoefficient"
/**
 * @brief Implementation of \ref ClipModelCoefficient().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see ClipModelCoefficient()
 */
PetscReal ClipModelCoefficient(PetscReal coefficient, const LESConfig *config, PetscBool *limited)
{
    PetscReal result = coefficient;

    switch (config->clip_mode) {
    case LES_CLIP_NONE:
        // Signed coefficient retained: backscatter survives, and stability is left to
        // the total-viscosity floor rather than to the coefficient's sign.
        break;
    case LES_CLIP_CLIP_NEGATIVE:
        if (result < 0.0) result = 0.0;
        break;
    case LES_CLIP_CLAMP:
    default:
        if (result < 0.0) result = 0.0;
        if (result > config->max_cs * config->max_cs) result = config->max_cs * config->max_cs;
        break;
    }

    if (limited) *limited = (PetscBool)(result != coefficient);
    return result;
}

#undef __FUNCT__
#define __FUNCT__ "EddyViscosityFromCoefficient"
/**
 * @brief Implementation of \ref EddyViscosityFromCoefficient().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see EddyViscosityFromCoefficient()
 */
PetscReal EddyViscosityFromCoefficient(PetscReal coefficient, PetscReal delta,
                                       PetscReal strain_magnitude,
                                       PetscReal molecular_viscosity,
                                       PetscReal min_viscosity_ratio)
{
    PetscReal nu_t = coefficient * delta * delta * strain_magnitude;

    // A negative eddy viscosity is physical up to the point where it would drive the
    // total viscosity non-positive and make the momentum operator ill posed.
    const PetscReal floor = (min_viscosity_ratio - 1.0) * molecular_viscosity;

    if (nu_t < floor) nu_t = floor;
    return nu_t;
}

#undef __FUNCT__
#define __FUNCT__ "SubgridKineticEnergy"
/**
 * @brief Implementation of \ref SubgridKineticEnergy().
 * @details Full API contract is documented with the header declaration in
 *          `include/les.h`.
 * @see SubgridKineticEnergy()
 */
PetscReal SubgridKineticEnergy(PetscReal yoshizawa_ci, PetscReal delta, PetscReal strain_magnitude)
{
    return 2.0 * yoshizawa_ci * delta * delta * strain_magnitude * strain_magnitude;
}

/*================================================================================*
 *                          DYNAMIC PROCEDURE INTERNALS                           *
 *================================================================================*/

/**
 * @brief Maps a cell index to the index whose strain rate represents it.
 *
 * The strain rate must be known one cell beyond the computed range, because the
 * test-filter stencil reaches there. At an interior index the answer is the index
 * itself. At the two layout-boundary planes it is not, and for opposite reasons.
 *
 * On a periodic axis, index 0 and index `m-1` are duplicates of the last and first
 * interior cells. A central difference evaluated at index 0 straddles the wrap and
 * returns nonsense, so the strain is evaluated instead at ghost index `-2`, which
 * addresses the same physical cell from the interior side. Index `m+1` plays the
 * same role at the other end. These are the ghost indices the periodic cell
 * synchronization itself copies from, and the DMDA carries the three ghost layers
 * they need whenever any axis is periodic.
 *
 * On a non-periodic axis the two planes are the solver's boundary layer, where a
 * central difference would reach outside the domain. The nearest interior cell's
 * strain is used instead, which is a zeroth-order extrapolation into a layer that
 * only ever contributes to a filter average.
 */
static PetscInt LESStrainSampleIndex(PetscInt index, PetscInt extent, PetscBool periodic)
{
    if (index == 0)          return periodic ? -2 : 1;
    if (index == extent - 1) return periodic ? extent + 1 : extent - 2;
    return index;
}

/**
 * @brief Fills the strain-rate scratch fields across the computed range and its halo.
 */
static PetscErrorCode ComputeStrainRateField(UserCtx *user, Cmpnts ***ucat,
                                             Vec strain_diagonal, Vec strain_offdiagonal,
                                             Vec strain_magnitude)
{
    DM                da = user->da, fda = user->fda;
    DMDALocalInfo     info;
    SpatialTargetPlan plan;
    PetscBool         periodic[3];
    Cmpnts         ***diagonal = NULL, ***offdiagonal = NULL;
    PetscReal      ***magnitude = NULL;

    PetscFunctionBeginUser;

    PetscCall(DMDAGetLocalInfo(da, &info));
    PetscCall(LESResolveDomain(user, &plan));
    LESPeriodicAxes(user->simCtx, periodic);

    PetscCall(DMDAVecGetArray(fda, strain_diagonal, &diagonal));
    PetscCall(DMDAVecGetArray(fda, strain_offdiagonal, &offdiagonal));
    PetscCall(DMDAVecGetArray(da, strain_magnitude, &magnitude));

    // One cell beyond the computed range in every direction, because that is how far
    // the 3x3x3 test-filter stencil reaches from the outermost computed cell.
    for (PetscInt k = plan.start[2] - 1; k <= plan.end[2]; k++)
    for (PetscInt j = plan.start[1] - 1; j <= plan.end[1]; j++)
    for (PetscInt i = plan.start[0] - 1; i <= plan.end[0]; i++) {
        const PetscInt si = LESStrainSampleIndex(i, info.mx, periodic[0]);
        const PetscInt sj = LESStrainSampleIndex(j, info.my, periodic[1]);
        const PetscInt sk = LESStrainSampleIndex(k, info.mz, periodic[2]);

        Cmpnts    dudx, dvdx, dwdx;
        SymTensor strain;
        PetscReal strain_abs;

        PetscCall(ComputeVectorFieldDerivatives(user, si, sj, sk, ucat, &dudx, &dvdx, &dwdx));
        PetscCall(StrainRateFromGradients(dudx, dvdx, dwdx, &strain, &strain_abs));

        diagonal[k][j][i].x    = strain.xx;
        diagonal[k][j][i].y    = strain.yy;
        diagonal[k][j][i].z    = strain.zz;
        offdiagonal[k][j][i].x = strain.xy;
        offdiagonal[k][j][i].y = strain.xz;
        offdiagonal[k][j][i].z = strain.yz;
        magnitude[k][j][i]     = strain_abs;
    }

    PetscCall(DMDAVecRestoreArray(fda, strain_diagonal, &diagonal));
    PetscCall(DMDAVecRestoreArray(fda, strain_offdiagonal, &offdiagonal));
    PetscCall(DMDAVecRestoreArray(da, strain_magnitude, &magnitude));

    PetscFunctionReturn(0);
}

/**
 * @brief Rebuilds a symmetric strain tensor from its two scratch storage vectors.
 */
static inline SymTensor LESStrainAt(Cmpnts diagonal, Cmpnts offdiagonal)
{
    SymTensor s;

    s.xx = diagonal.x;
    s.yy = diagonal.y;
    s.zz = diagonal.z;
    s.xy = offdiagonal.x;
    s.xz = offdiagonal.y;
    s.yz = offdiagonal.z;
    return s;
}

/*================================================================================*
 *                              RUNTIME ENTRY POINTS                              *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "ComputeSmagorinskyConstant"
/**
 * @brief Implementation of \ref ComputeSmagorinskyConstant().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/les.h`.
 * @see ComputeSmagorinskyConstant()
 */
PetscErrorCode ComputeSmagorinskyConstant(UserCtx *user)
{
    SimCtx          *simCtx = user->simCtx;
    const LESConfig *config = &simCtx->les_config;
    DM                da = user->da, fda = user->fda;
    SpatialTargetPlan plan;
    PetscBool         average_direction[3];
    PetscReal        alpha;

    Vec lStrainDiag, lStrainOff, lStrainAbs, lLM, lMM, lCsq;

    Cmpnts    ***ucat = NULL, ***strain_diag = NULL, ***strain_off = NULL;
    Cmpnts    ***csi = NULL, ***eta = NULL, ***zet = NULL;
    PetscReal ***strain_abs = NULL, ***nvert = NULL, ***aj = NULL;
    PetscReal ***lm = NULL, ***mm = NULL, ***csq = NULL, ***cs_out = NULL;

    LESDiagnosticsState diagnostics = {0.0, 0.0, 0.0, 0.0, 0.0, PETSC_TRUE};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    PetscCheck(simCtx->les == DYNAMIC_SMAGORINSKY, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "The dynamic procedure runs only for the dynamic Smagorinsky model. The "
               "constant model prescribes its coefficient from configuration and holds "
               "no coefficient field.");

    // The dynamic procedure needs a developed field to sample. Starting from rest it
    // has nothing to measure, so the first steps run with no subgrid dissipation.
    if (simCtx->step < 2 && simCtx->StartStep == 0) {
        PetscCall(VecSet(user->CS, 0.0));
        PetscCall(FinalizeSmagorinskyConstantField(user));
        user->les_diagnostics = (LESDiagnosticsState){0.0, 0.0, 0.0, 0.0, 0.0, PETSC_FALSE};
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "Holding the dynamic coefficient at zero for the opening steps (step=%d).\n",
                  simCtx->step);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    PetscCall(LESResolveDomain(user, &plan));
    PetscCall(ResolveLESAveragingDirections(user, average_direction));

    // alpha is the squared ratio of the two filter widths and sets the relative weight
    // of the test-scale and grid-scale terms in the model tensor.
    alpha = config->test_filter_width_ratio * config->test_filter_width_ratio;

    PetscCall(VecDuplicate(user->lUcat, &lStrainDiag));
    PetscCall(VecDuplicate(user->lUcat, &lStrainOff));
    PetscCall(VecDuplicate(user->lNvert, &lStrainAbs));
    PetscCall(VecDuplicate(user->lNvert, &lLM));
    PetscCall(VecDuplicate(user->lNvert, &lMM));
    PetscCall(VecDuplicate(user->lNvert, &lCsq));
    PetscCall(VecSet(lStrainDiag, 0.0));
    PetscCall(VecSet(lStrainOff, 0.0));
    PetscCall(VecSet(lStrainAbs, 0.0));
    PetscCall(VecSet(lLM, 0.0));
    PetscCall(VecSet(lMM, 0.0));
    PetscCall(VecSet(lCsq, 0.0));

    PetscCall(DMDAVecGetArray(fda, user->lUcat, &ucat));

    // --- 1. Strain rate everywhere the filter stencil will look -----------------
    PetscCall(ComputeStrainRateField(user, ucat, lStrainDiag, lStrainOff, lStrainAbs));

    PetscCall(DMDAVecGetArrayRead(fda, lStrainDiag, &strain_diag));
    PetscCall(DMDAVecGetArrayRead(fda, lStrainOff, &strain_off));
    PetscCall(DMDAVecGetArrayRead(da, lStrainAbs, &strain_abs));
    PetscCall(DMDAVecGetArrayRead(fda, user->lCsi, &csi));
    PetscCall(DMDAVecGetArrayRead(fda, user->lEta, &eta));
    PetscCall(DMDAVecGetArrayRead(fda, user->lZet, &zet));
    PetscCall(DMDAVecGetArrayRead(da, user->lNvert, &nvert));
    PetscCall(DMDAVecGetArrayRead(da, user->lAj, &aj));
    PetscCall(DMDAVecGetArray(da, lLM, &lm));
    PetscCall(DMDAVecGetArray(da, lMM, &mm));

    // --- 2. Germano contractions at every fluid cell ----------------------------
    for (PetscInt k = plan.start[2]; k < plan.end[2]; k++)
    for (PetscInt j = plan.start[1]; j < plan.end[1]; j++)
    for (PetscInt i = plan.start[0]; i < plan.end[0]; i++) {
        double    velocity_x[3][3][3], velocity_y[3][3][3], velocity_z[3][3][3];
        double    strain_abs_stencil[3][3][3], weights[3][3][3];
        SymTensor strain_stencil[3][3][3], product_stencil[3][3][3], velocity_product[3][3][3];
        Cmpnts    velocity_filtered;
        SymTensor strain_filtered, product_filtered, velocity_product_filtered;
        SymTensor leonard, model;
        PetscReal strain_abs_filtered, delta, volume;

        if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) {
            lm[k][j][i] = 0.0;
            mm[k][j][i] = 0.0;
            continue;
        }

        // --- 2a. Gather the 3x3x3 stencil -------------------------------------
        for (PetscInt r = -1; r <= 1; r++)
        for (PetscInt q = -1; q <= 1; q++)
        for (PetscInt p = -1; p <= 1; p++) {
            const PetscInt  R = r + 1, Q = q + 1, P = p + 1;
            const PetscInt  KK = k + r, JJ = j + q, II = i + p;
            const SymTensor strain = LESStrainAt(strain_diag[KK][JJ][II], strain_off[KK][JJ][II]);
            const PetscReal strain_abs_here = strain_abs[KK][JJ][II];

            velocity_x[R][Q][P] = ucat[KK][JJ][II].x;
            velocity_y[R][Q][P] = ucat[KK][JJ][II].y;
            velocity_z[R][Q][P] = ucat[KK][JJ][II].z;

            strain_stencil[R][Q][P]     = strain;
            strain_abs_stencil[R][Q][P] = strain_abs_here;

            // The grid-scale model formed pointwise, before any test filtering. Test
            // filtering this product is the term the model tensor cannot do without,
            // and it is not the same as multiplying the two filtered factors.
            product_stencil[R][Q][P] = SymTensorCombine(strain_abs_here, strain, 0.0, strain);

            velocity_product[R][Q][P] = SymTensorSelfOuter(ucat[KK][JJ][II]);
            weights[R][Q][P]          = aj[KK][JJ][II];
        }

        // --- 2b. Test-filter every quantity the identity needs ----------------
        velocity_filtered.x = ApplyLESTestFilter(config->test_filter_kernel, velocity_x, weights);
        velocity_filtered.y = ApplyLESTestFilter(config->test_filter_kernel, velocity_y, weights);
        velocity_filtered.z = ApplyLESTestFilter(config->test_filter_kernel, velocity_z, weights);
        strain_abs_filtered = ApplyLESTestFilter(config->test_filter_kernel, strain_abs_stencil, weights);
        PetscCall(ApplyLESTestFilterSymTensor(config->test_filter_kernel, strain_stencil, weights, &strain_filtered));
        PetscCall(ApplyLESTestFilterSymTensor(config->test_filter_kernel, product_stencil, weights, &product_filtered));
        PetscCall(ApplyLESTestFilterSymTensor(config->test_filter_kernel, velocity_product, weights, &velocity_product_filtered));

        // --- 2c. Leonard stress and model tensor -------------------------------
        leonard = LeonardStress(velocity_filtered, velocity_product_filtered);

        PetscCall(ComputeCellFilterWidth(config->filter_width_model, aj[k][j][i],
                                         csi[k][j][i], eta[k][j][i], zet[k][j][i], &delta));

        model = GermanoModelTensor(delta, alpha, strain_abs_filtered, strain_filtered, product_filtered);

        // --- 2d. Lilly's least-squares contraction -----------------------------
        // Stored pre-multiplied by cell volume. The average that follows is over
        // volume, not over cell count, so that a refined region does not dominate it
        // merely by holding more cells; carrying the weight here keeps the shared
        // averaging kernel free of any weighting policy of its own.
        volume = 1.0 / aj[k][j][i];
        lm[k][j][i] = volume * SymTensorContract(leonard, model);
        mm[k][j][i] = volume * SymTensorNormSq(model);
    }

    PetscCall(DMDAVecRestoreArrayRead(fda, lStrainDiag, &strain_diag));
    PetscCall(DMDAVecRestoreArrayRead(fda, lStrainOff, &strain_off));
    PetscCall(DMDAVecRestoreArray(da, lLM, &lm));
    PetscCall(DMDAVecRestoreArray(da, lMM, &mm));
    PetscCall(DMDAVecRestoreArray(fda, user->lUcat, &ucat));

    // --- 3. Average the two contractions, then divide ---------------------------
    // Both sums are formed before the division. Averaging the pointwise quotients
    // instead would not be the least-squares coefficient Lilly's closure defines.
    PetscCall(PicurvSpatialRatioAverage(user, &plan, lLM, lMM, NULL, average_direction,
                                        PetscObjectComm((PetscObject)da), lCsq, NULL));

    // --- 4. Limit the coefficient and record what limiting cost -----------------
    PetscCall(DMDAVecGetArrayRead(da, lLM, &lm));
    PetscCall(DMDAVecGetArrayRead(da, lMM, &mm));
    PetscCall(DMDAVecGetArrayRead(da, lCsq, &csq));
    PetscCall(DMDAVecGetArray(da, user->CS, &cs_out));

    for (PetscInt k = plan.start[2]; k < plan.end[2]; k++)
    for (PetscInt j = plan.start[1]; j < plan.end[1]; j++)
    for (PetscInt i = plan.start[0]; i < plan.end[0]; i++) {
        PetscBool limited = PETSC_FALSE;
        PetscReal weight, raw;

        if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) {
            cs_out[k][j][i] = 0.0;
            continue;
        }

        weight = 1.0 / aj[k][j][i];
        raw    = csq[k][j][i];

        cs_out[k][j][i] = ClipModelCoefficient(raw, config, &limited);

        // The two contraction fields already carry the cell volume from the loop
        // above, so they are summed as they stand rather than weighted twice.
        diagnostics.contraction_lm     += lm[k][j][i];
        diagnostics.contraction_mm     += mm[k][j][i];
        diagnostics.fluid_volume       += weight;
        if (raw < 0.0) diagnostics.backscatter_volume += weight;
        if (limited)   diagnostics.limited_volume     += weight;
    }

    PetscCall(DMDAVecRestoreArrayRead(da, lLM, &lm));
    PetscCall(DMDAVecRestoreArrayRead(da, lMM, &mm));
    PetscCall(DMDAVecRestoreArrayRead(da, lCsq, &csq));
    PetscCall(DMDAVecRestoreArray(da, user->CS, &cs_out));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lCsi, &csi));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lEta, &eta));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lZet, &zet));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lNvert, &nvert));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lAj, &aj));

    PetscCall(VecDestroy(&lStrainDiag));
    PetscCall(VecDestroy(&lStrainOff));
    PetscCall(VecDestroy(&lStrainAbs));
    PetscCall(VecDestroy(&lLM));
    PetscCall(VecDestroy(&lMM));
    PetscCall(VecDestroy(&lCsq));

    user->les_diagnostics = diagnostics;
    PetscCall(FinalizeSmagorinskyConstantField(user));

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeEddyViscosityLES"
/**
 * @brief Implementation of \ref ComputeEddyViscosityLES().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/les.h`.
 * @see ComputeEddyViscosityLES()
 */
PetscErrorCode ComputeEddyViscosityLES(UserCtx *user)
{
    SimCtx          *simCtx = user->simCtx;
    const LESConfig *config = &simCtx->les_config;
    DM                da = user->da, fda = user->fda;
    SpatialTargetPlan plan;
    const PetscReal  molecular_viscosity = 1.0 / simCtx->ren;

    // The constant model's coefficient is a number from the configuration, so it needs
    // no field, no ghost exchange, and no periodic synchronization.
    const PetscBool  coefficient_is_field = (PetscBool)(simCtx->les == DYNAMIC_SMAGORINSKY);
    const PetscReal  prescribed_coefficient = config->constant_cs * config->constant_cs;

    Cmpnts    ***ucat = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;
    PetscReal ***coefficient = NULL, ***nu_t = NULL, ***nvert = NULL, ***aj = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    PetscCall(LESResolveDomain(user, &plan));

    PetscCall(DMDAVecGetArrayRead(fda, user->lUcat, &ucat));
    PetscCall(DMDAVecGetArrayRead(fda, user->lCsi, &csi));
    PetscCall(DMDAVecGetArrayRead(fda, user->lEta, &eta));
    PetscCall(DMDAVecGetArrayRead(fda, user->lZet, &zet));
    PetscCall(DMDAVecGetArray(da, user->Nu_t, &nu_t));
    PetscCall(DMDAVecGetArrayRead(da, user->lNvert, &nvert));
    PetscCall(DMDAVecGetArrayRead(da, user->lAj, &aj));
    if (coefficient_is_field) PetscCall(DMDAVecGetArrayRead(da, user->lCs, &coefficient));

    for (PetscInt k = plan.start[2]; k < plan.end[2]; k++)
    for (PetscInt j = plan.start[1]; j < plan.end[1]; j++)
    for (PetscInt i = plan.start[0]; i < plan.end[0]; i++) {
        Cmpnts    dudx, dvdx, dwdx;
        PetscReal strain_magnitude, delta, model_coefficient;

        if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) {
            nu_t[k][j][i] = 0.0;
            continue;
        }

        PetscCall(ComputeVectorFieldDerivatives(user, i, j, k, (Cmpnts ***)ucat, &dudx, &dvdx, &dwdx));
        PetscCall(StrainRateFromGradients(dudx, dvdx, dwdx, NULL, &strain_magnitude));
        PetscCall(ComputeCellFilterWidth(config->filter_width_model, aj[k][j][i],
                                         csi[k][j][i], eta[k][j][i], zet[k][j][i], &delta));

        model_coefficient = coefficient_is_field ? coefficient[k][j][i] : prescribed_coefficient;

        nu_t[k][j][i] = EddyViscosityFromCoefficient(model_coefficient, delta, strain_magnitude,
                                                     molecular_viscosity, config->min_viscosity_ratio);
    }

    if (coefficient_is_field) PetscCall(DMDAVecRestoreArrayRead(da, user->lCs, &coefficient));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lUcat, &ucat));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lCsi, &csi));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lEta, &eta));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lZet, &zet));
    PetscCall(DMDAVecRestoreArray(da, user->Nu_t, &nu_t));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lNvert, &nvert));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lAj, &aj));

    {
        const FieldId periodic_fields[] = {FIELD_ID_NU_T};

        PetscCall(SynchronizePeriodicCellFields(user, 1, periodic_fields));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_NU_T));
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LogLESDiagnostics"
/**
 * @brief Implementation of \ref LogLESDiagnostics().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/les.h`.
 * @see LogLESDiagnostics()
 */
PetscErrorCode LogLESDiagnostics(UserCtx *user)
{
    SimCtx          *simCtx = user->simCtx;
    const LESConfig *config = &simCtx->les_config;
    DM                da = user->da, fda = user->fda;
    SpatialTargetPlan plan;
    MPI_Comm         comm;

    const PetscBool dynamic = (PetscBool)(simCtx->les == DYNAMIC_SMAGORINSKY);

    // Index 0..4: fluid volume, coefficient first and second moments, nu_t volume sum,
    // subgrid-energy volume sum. Index 5..6: the pre-clip fractions. Reduced together
    // so the diagnostic costs one collective rather than seven.
    PetscReal local_sum[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    PetscReal global_sum[7];
    PetscReal local_max[2] = {0.0, 0.0}, global_max[2];
    PetscReal local_min_coefficient = PETSC_MAX_REAL, global_min_coefficient;
    PetscReal contraction[2], global_contraction[2];

    Cmpnts    ***ucat = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;
    PetscReal ***coefficient = NULL, ***nu_t = NULL, ***nvert = NULL, ***aj = NULL;

    PetscFunctionBeginUser;

    if (!config->diagnostics_enabled) PetscFunctionReturn(0);
    if (config->diagnostics_cadence <= 0) PetscFunctionReturn(0);
    if (simCtx->step % config->diagnostics_cadence != 0) PetscFunctionReturn(0);

    PetscCall(LESResolveDomain(user, &plan));
    PetscCall(PetscObjectGetComm((PetscObject)da, &comm));

    PetscCall(DMDAVecGetArrayRead(fda, user->lUcat, &ucat));
    PetscCall(DMDAVecGetArrayRead(fda, user->lCsi, &csi));
    PetscCall(DMDAVecGetArrayRead(fda, user->lEta, &eta));
    PetscCall(DMDAVecGetArrayRead(fda, user->lZet, &zet));
    PetscCall(DMDAVecGetArrayRead(da, user->lNu_t, &nu_t));
    PetscCall(DMDAVecGetArrayRead(da, user->lNvert, &nvert));
    PetscCall(DMDAVecGetArrayRead(da, user->lAj, &aj));
    if (dynamic) PetscCall(DMDAVecGetArrayRead(da, user->lCs, &coefficient));

    for (PetscInt k = plan.start[2]; k < plan.end[2]; k++)
    for (PetscInt j = plan.start[1]; j < plan.end[1]; j++)
    for (PetscInt i = plan.start[0]; i < plan.end[0]; i++) {
        Cmpnts    dudx, dvdx, dwdx;
        PetscReal strain_magnitude, delta, weight, value;

        if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) continue;

        weight = 1.0 / aj[k][j][i];
        value  = dynamic ? coefficient[k][j][i]
                         : config->constant_cs * config->constant_cs;

        PetscCall(ComputeVectorFieldDerivatives(user, i, j, k, (Cmpnts ***)ucat, &dudx, &dvdx, &dwdx));
        PetscCall(StrainRateFromGradients(dudx, dvdx, dwdx, NULL, &strain_magnitude));
        PetscCall(ComputeCellFilterWidth(config->filter_width_model, aj[k][j][i],
                                         csi[k][j][i], eta[k][j][i], zet[k][j][i], &delta));

        local_sum[0] += weight;
        local_sum[1] += weight * value;
        local_sum[2] += weight * value * value;
        local_sum[3] += weight * nu_t[k][j][i];
        local_sum[4] += weight * SubgridKineticEnergy(config->yoshizawa_ci, delta, strain_magnitude);

        local_max[0] = PetscMax(local_max[0], value);
        local_max[1] = PetscMax(local_max[1], nu_t[k][j][i]);
        local_min_coefficient = PetscMin(local_min_coefficient, value);
    }

    if (dynamic) PetscCall(DMDAVecRestoreArrayRead(da, user->lCs, &coefficient));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lUcat, &ucat));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lCsi, &csi));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lEta, &eta));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lZet, &zet));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lNu_t, &nu_t));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lNvert, &nvert));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lAj, &aj));

    local_sum[5] = user->les_diagnostics.valid ? user->les_diagnostics.backscatter_volume : 0.0;
    local_sum[6] = user->les_diagnostics.valid ? user->les_diagnostics.limited_volume : 0.0;
    contraction[0] = user->les_diagnostics.valid ? user->les_diagnostics.contraction_lm : 0.0;
    contraction[1] = user->les_diagnostics.valid ? user->les_diagnostics.contraction_mm : 0.0;

    PetscCallMPI(MPI_Allreduce(local_sum, global_sum, 7, MPIU_REAL, MPI_SUM, comm));
    PetscCallMPI(MPI_Allreduce(local_max, global_max, 2, MPIU_REAL, MPI_MAX, comm));
    PetscCallMPI(MPI_Allreduce(&local_min_coefficient, &global_min_coefficient, 1, MPIU_REAL, MPI_MIN, comm));
    PetscCallMPI(MPI_Allreduce(contraction, global_contraction, 2, MPIU_REAL, MPI_SUM, comm));

    if (simCtx->rank == 0) {
        const PetscReal volume        = (global_sum[0] > 0.0) ? global_sum[0] : 1.0;
        const PetscReal mean          = global_sum[1] / volume;
        const PetscReal mean_square   = global_sum[2] / volume;
        const PetscReal variance      = PetscMax(mean_square - mean * mean, 0.0);
        const PetscReal molecular     = 1.0 / simCtx->ren;
        const PetscReal effective     = (global_contraction[1] > LES_EPSILON)
                                            ? (global_contraction[0] / global_contraction[1]) : 0.0;
        char  path[PETSC_MAX_PATH_LEN + 32];
        FILE *file = NULL;

        // Reported as Cs rather than as the stored coefficient, because that is the
        // number with a recognised value: near 0.16 for decaying isotropic turbulence.
        // A negative average means the box is backscattering on balance, and is written
        // as a negative Cs rather than hidden behind a square root.
        const PetscReal cs_effective = (effective >= 0.0) ? PetscSqrtReal(effective)
                                                          : -PetscSqrtReal(-effective);
        const PetscReal cs_mean      = (mean >= 0.0) ? PetscSqrtReal(mean) : -PetscSqrtReal(-mean);

        PetscCall(PetscSNPrintf(path, sizeof(path), "%s/les_coefficient.csv", simCtx->log_dir));
        file = fopen(path, "a");
        PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "Unable to open LES diagnostics file '%s'.", path);
        if (ftell(file) == 0) {
            fprintf(file,
                    "step,time,cs_effective,cs_mean,coefficient_mean,coefficient_rms,"
                    "coefficient_min,coefficient_max,nu_t_mean,nu_t_max,nu_t_over_nu_mean,"
                    "k_sgs_mean,backscatter_fraction,limited_fraction\n");
        }
        if (simCtx->continueMode && simCtx->step == simCtx->StartStep + 1) {
            fprintf(file, "# Continuation from step %" PetscInt_FMT "\n", simCtx->StartStep);
        }
        fprintf(file,
                "%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                (int)simCtx->step, (double)simCtx->ti,
                (double)cs_effective, (double)cs_mean,
                (double)mean, (double)PetscSqrtReal(variance),
                (double)global_min_coefficient, (double)global_max[0],
                (double)(global_sum[3] / volume), (double)global_max[1],
                (double)(global_sum[3] / volume / molecular),
                (double)(global_sum[4] / volume),
                (double)(global_sum[5] / volume), (double)(global_sum[6] / volume));
        PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close LES diagnostics file '%s'.", path);

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "  LES coefficient: Cs_effective=%.4f, nu_t/nu (mean)=%.3e, backscatter=%.1f%%\n",
                  (double)cs_effective, (double)(global_sum[3] / volume / molecular),
                  (double)(100.0 * global_sum[5] / volume));
    }

    PetscFunctionReturn(0);
}
