/**
 * @file test_les.c
 * @brief C unit tests for the LES subgrid-scale closure.
 *
 * The suite is layered the way the module is. Pure kernels are checked against hand
 * computed values with no PETSc objects at all; the Germano pieces are checked on
 * fields whose answers are known in closed form; and the driver routines are checked
 * on a small DMDA fixture.
 *
 * Two of the cases exist because of specific defects and are worth naming. The
 * "filtered-product" case pins the distinction between the filter of a product and
 * the product of filtered factors, which the model tensor once collapsed. The
 * "duplicate-plane" case pins the requirement that a spatial average count each
 * physical cell once on a periodic block, where the two layout boundary planes are
 * copies of interior cells.
 */

#include "test_support.h"

#include "les.h"
#include "setup.h"
#include "statistics_target.h"

/*================================================================================*
 *                              KERNEL-LEVEL TESTS                                *
 *================================================================================*/

/** @brief Tests the symmetric-tensor primitives against hand-computed values. */
static PetscErrorCode TestSymTensorAlgebra(void)
{
    const Cmpnts    v = {1.0, 2.0, 3.0};
    const SymTensor a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const SymTensor b = {2.0, 0.0, 0.0, 2.0, 0.0, 2.0};
    SymTensor       outer, deviator, combined;

    PetscFunctionBeginUser;

    outer = SymTensorSelfOuter(v);
    PetscCall(PicurvAssertRealNear(1.0, outer.xx, 1.0e-14, "self outer product should square the x component"));
    PetscCall(PicurvAssertRealNear(2.0, outer.xy, 1.0e-14, "self outer product should cross x and y"));
    PetscCall(PicurvAssertRealNear(6.0, outer.yz, 1.0e-14, "self outer product should cross y and z"));
    PetscCall(PicurvAssertRealNear(9.0, outer.zz, 1.0e-14, "self outer product should square the z component"));

    /* trace(a) = 1 + 4 + 6 = 11 */
    PetscCall(PicurvAssertRealNear(11.0, SymTensorTrace(a), 1.0e-14, "trace should sum the diagonal"));

    deviator = SymTensorDeviator(a);
    PetscCall(PicurvAssertRealNear(0.0, SymTensorTrace(deviator), 1.0e-13,
                                   "the deviator should be trace free"));
    PetscCall(PicurvAssertRealNear(a.xy, deviator.xy, 1.0e-14,
                                   "the deviator should leave off-diagonal components alone"));

    /* a:b doubles the off-diagonal terms; b is diagonal, so a:b = 2*(1 + 4 + 6) = 22. */
    PetscCall(PicurvAssertRealNear(22.0, SymTensorContract(a, b), 1.0e-13,
                                   "contraction against a diagonal tensor should scale its trace"));
    /* a:a = 1 + 16 + 36 + 2*(4 + 9 + 25) = 53 + 76 = 129 */
    PetscCall(PicurvAssertRealNear(129.0, SymTensorNormSq(a), 1.0e-13,
                                   "the squared norm must count off-diagonal components twice"));

    combined = SymTensorCombine(2.0, a, -1.0, b);
    PetscCall(PicurvAssertRealNear(0.0, combined.xx, 1.0e-14, "combine should scale and subtract"));
    PetscCall(PicurvAssertRealNear(4.0, combined.xy, 1.0e-14, "combine should scale off-diagonal terms"));

    PetscFunctionReturn(0);
}

/** @brief Tests strain-rate assembly and its magnitude for a known velocity gradient. */
static PetscErrorCode TestStrainRateFromGradients(void)
{
    /* du/dy = 2 and every other derivative zero: pure shear. */
    const Cmpnts dudx = {0.0, 2.0, 0.0};
    const Cmpnts dvdx = {0.0, 0.0, 0.0};
    const Cmpnts dwdx = {0.0, 0.0, 0.0};
    SymTensor    strain;
    PetscReal    magnitude;

    PetscFunctionBeginUser;
    PetscCall(StrainRateFromGradients(dudx, dvdx, dwdx, &strain, &magnitude));

    PetscCall(PicurvAssertRealNear(1.0, strain.xy, 1.0e-14,
                                   "pure shear should give S_xy = half the velocity gradient"));
    PetscCall(PicurvAssertRealNear(0.0, strain.xx, 1.0e-14, "pure shear should have no normal strain"));
    PetscCall(PicurvAssertRealNear(0.0, SymTensorTrace(strain), 1.0e-14,
                                   "an incompressible gradient should give a trace-free strain"));
    /* |S| = sqrt(2 S_ij S_ij) = sqrt(2 * 2 * 1) = 2 */
    PetscCall(PicurvAssertRealNear(2.0, magnitude, 1.0e-13,
                                   "the strain magnitude must count both off-diagonal entries"));

    PetscFunctionReturn(0);
}

/** @brief Tests that the three filter-width models separate on an anisotropic cell. */
static PetscErrorCode TestFilterWidthModelsSeparateOnStretchedCell(void)
{
    /* Cartesian-aligned metrics for a cell of size 4 x 1 x 1: the face area vectors
       are the products of the two spanning edges, and the Jacobian is 1/volume. */
    const Cmpnts    csi = {1.0, 0.0, 0.0};
    const Cmpnts    eta = {0.0, 4.0, 0.0};
    const Cmpnts    zet = {0.0, 0.0, 4.0};
    const PetscReal aj  = 0.25;
    PetscReal       cube_root = 0.0, geometric = 0.0, max_edge = 0.0;

    PetscFunctionBeginUser;

    PetscCall(ComputeCellFilterWidth(LES_FILTER_WIDTH_CUBE_ROOT_VOLUME, aj, csi, eta, zet, &cube_root));
    PetscCall(ComputeCellFilterWidth(LES_FILTER_WIDTH_GEOMETRIC_MEAN, aj, csi, eta, zet, &geometric));
    PetscCall(ComputeCellFilterWidth(LES_FILTER_WIDTH_MAX_EDGE, aj, csi, eta, zet, &max_edge));

    /* Volume is 4, so the cube-root model gives 4^(1/3). */
    PetscCall(PicurvAssertRealNear(PetscPowReal(4.0, 1.0 / 3.0), cube_root, 1.0e-12,
                                   "the cube-root model should depend only on cell volume"));
    PetscCall(PicurvAssertRealNear(PetscPowReal(4.0, 1.0 / 3.0), geometric, 1.0e-12,
                                   "the geometric mean of the extents should recover the same volume"));
    PetscCall(PicurvAssertRealNear(4.0, max_edge, 1.0e-12,
                                   "the max-edge model should report the long direction"));
    PetscCall(PicurvAssertBool((PetscBool)(max_edge > cube_root),
                               "a stretched cell should make the max-edge width the larger of the two"));

    PetscFunctionReturn(0);
}

/** @brief Tests that the Leonard stress vanishes on a uniform velocity field. */
static PetscErrorCode TestLeonardStressVanishesOnUniformFlow(void)
{
    const Cmpnts    velocity = {1.5, -0.5, 2.0};
    const SymTensor product  = SymTensorSelfOuter(velocity);
    SymTensor       leonard;

    PetscFunctionBeginUser;
    /* Filtering a constant returns the constant, so both filtered inputs are exact. */
    leonard = LeonardStress(velocity, product);

    PetscCall(PicurvAssertRealNear(0.0, SymTensorNormSq(leonard), 1.0e-24,
                                   "a uniform field carries no stress between the two filter widths"));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests the model tensor where the two filter terms provably coincide.
 *
 * On a constant strain field the test filter is the identity, so the filter of the
 * product equals the product of the filtered factors and the model tensor collapses
 * to a closed form, `-2 Delta^2 (alpha - 1) |S| S_ij`, made trace free. This is the
 * one configuration in which the corrected tensor and the collapsed expression that
 * preceded it agree, which is what makes it a clean analytic anchor.
 */
static PetscErrorCode TestGermanoModelTensorOnConstantStrain(void)
{
    const PetscReal delta = 0.5;
    const PetscReal alpha = 4.0;
    Cmpnts          dudx = {0.0, 2.0, 0.0}, dvdx = {0.0, 0.0, 0.0}, dwdx = {0.0, 0.0, 0.0};
    SymTensor       strain, product, model, expected;
    PetscReal       magnitude;

    PetscFunctionBeginUser;
    PetscCall(StrainRateFromGradients(dudx, dvdx, dwdx, &strain, &magnitude));

    product = SymTensorCombine(magnitude, strain, 0.0, strain);
    model   = GermanoModelTensor(delta, alpha, magnitude, strain, product);

    expected = SymTensorDeviator(SymTensorCombine(-2.0 * delta * delta * (alpha - 1.0) * magnitude,
                                                  strain, 0.0, strain));

    PetscCall(PicurvAssertRealNear(expected.xy, model.xy, 1.0e-13,
                                   "constant strain should give the closed-form model tensor"));
    PetscCall(PicurvAssertRealNear(0.0, SymTensorTrace(model), 1.0e-13,
                                   "the model tensor must be trace free"));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the model tensor uses the filtered product, not the filtered factors.
 *
 * Builds a stencil whose strain magnitude varies across it, so the test filter of
 * `|S| S_ij` genuinely differs from `|S|^ S^_ij`. The model tensor must be built from
 * the former. Collapsing it onto the latter yields `-2 Delta^2 (alpha - 1) |S|^ S^_ij`,
 * which is well scaled, responds to the flow, and is not the Germano-Lilly tensor;
 * this case fails if the two are ever conflated again.
 */
static PetscErrorCode TestGermanoModelTensorUsesFilteredProduct(void)
{
    const LESTestFilterKernel kernel = LES_TEST_FILTER_VOLUME_WEIGHTED_BOX;
    const PetscReal           delta = 1.0;
    const PetscReal           alpha = 4.0;
    double                    magnitude_stencil[3][3][3], weights[3][3][3];
    SymTensor                 strain_stencil[3][3][3], product_stencil[3][3][3];
    SymTensor                 strain_filtered, product_filtered, model, collapsed;
    PetscReal                 magnitude_filtered, separation;

    PetscFunctionBeginUser;

    for (PetscInt k = 0; k < 3; ++k)
    for (PetscInt j = 0; j < 3; ++j)
    for (PetscInt i = 0; i < 3; ++i) {
        /* A strain that grows across the stencil: filtering no longer commutes with
           forming the product, which is exactly the regime the identity relies on. */
        const PetscReal shear = 1.0 + 0.5 * (PetscReal)(i + j + k);
        Cmpnts          dudx = {0.0, 2.0 * shear, 0.0};
        Cmpnts          dvdx = {0.0, 0.0, 0.0};
        Cmpnts          dwdx = {0.0, 0.0, 0.0};
        SymTensor       strain;
        PetscReal       magnitude;

        PetscCall(StrainRateFromGradients(dudx, dvdx, dwdx, &strain, &magnitude));
        strain_stencil[k][j][i]    = strain;
        magnitude_stencil[k][j][i] = magnitude;
        product_stencil[k][j][i]   = SymTensorCombine(magnitude, strain, 0.0, strain);
        weights[k][j][i]           = 1.0;
    }

    magnitude_filtered = ApplyLESTestFilter(kernel, magnitude_stencil, weights);
    PetscCall(ApplyLESTestFilterSymTensor(kernel, strain_stencil, weights, &strain_filtered));
    PetscCall(ApplyLESTestFilterSymTensor(kernel, product_stencil, weights, &product_filtered));

    /* The filter of the product must not equal the product of the filtered factors. */
    separation = PetscAbsReal(product_filtered.xy - magnitude_filtered * strain_filtered.xy);
    PetscCall(PicurvAssertBool((PetscBool)(separation > 1.0e-6),
                               "a varying strain field must separate the filtered product from "
                               "the product of the filtered factors"));

    model     = GermanoModelTensor(delta, alpha, magnitude_filtered, strain_filtered, product_filtered);
    collapsed = SymTensorDeviator(SymTensorCombine(-2.0 * delta * delta * (alpha - 1.0) * magnitude_filtered,
                                                   strain_filtered, 0.0, strain_filtered));

    PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(model.xy - collapsed.xy) > 1.0e-6),
                               "the model tensor must be built from the filtered product, not from "
                               "the product of separately filtered quantities"));
    PetscFunctionReturn(0);
}

/** @brief Tests each limiting mode, including the sign that carries backscatter. */
static PetscErrorCode TestClipModelCoefficientModes(void)
{
    LESConfig config;
    PetscBool limited = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(LESConfigSetDefaults(&config));
    config.max_cs = 0.2; /* ceiling on Cs, so the coefficient ceiling is 0.04 */

    config.clip_mode = LES_CLIP_CLAMP;
    PetscCall(PicurvAssertRealNear(0.0, ClipModelCoefficient(-0.01, &config, &limited), 1.0e-14,
                                   "clamping should discard a negative coefficient"));
    PetscCall(PicurvAssertBool(limited, "clamping a negative coefficient counts as limiting"));
    PetscCall(PicurvAssertRealNear(0.04, ClipModelCoefficient(1.0, &config, &limited), 1.0e-14,
                                   "clamping should cap the coefficient at max_cs squared"));
    PetscCall(PicurvAssertRealNear(0.01, ClipModelCoefficient(0.01, &config, &limited), 1.0e-14,
                                   "clamping should pass an admissible coefficient through"));
    PetscCall(PicurvAssertBool((PetscBool)(!limited), "an admissible coefficient is not limited"));

    config.clip_mode = LES_CLIP_CLIP_NEGATIVE;
    PetscCall(PicurvAssertRealNear(0.0, ClipModelCoefficient(-0.01, &config, &limited), 1.0e-14,
                                   "clip_negative should discard a negative coefficient"));
    PetscCall(PicurvAssertRealNear(1.0, ClipModelCoefficient(1.0, &config, &limited), 1.0e-14,
                                   "clip_negative should impose no ceiling"));

    config.clip_mode = LES_CLIP_NONE;
    PetscCall(PicurvAssertRealNear(-0.01, ClipModelCoefficient(-0.01, &config, &limited), 1.0e-14,
                                   "no limiting should preserve the sign that carries backscatter"));
    PetscCall(PicurvAssertBool((PetscBool)(!limited), "an unlimited coefficient is not modified"));
    PetscFunctionReturn(0);
}

/** @brief Tests eddy-viscosity assembly and the total-viscosity floor. */
static PetscErrorCode TestEddyViscosityFloorBoundsTotalViscosity(void)
{
    const PetscReal molecular = 0.01;

    PetscFunctionBeginUser;

    /* C = 0.04, Delta = 2, |S| = 3  ->  nu_t = 0.04 * 4 * 3 = 0.48 */
    PetscCall(PicurvAssertRealNear(0.48, EddyViscosityFromCoefficient(0.04, 2.0, 3.0, molecular, 0.0),
                                   1.0e-13, "eddy viscosity should be the coefficient times Delta^2 |S|"));

    /* A backscattering coefficient may drive nu_t negative, but only down to the point
       where the total viscosity would vanish. */
    PetscCall(PicurvAssertRealNear(-molecular,
                                   EddyViscosityFromCoefficient(-100.0, 2.0, 3.0, molecular, 0.0),
                                   1.0e-14, "the floor should stop nu_t at minus the molecular value"));
    PetscCall(PicurvAssertRealNear(-0.5 * molecular,
                                   EddyViscosityFromCoefficient(-100.0, 2.0, 3.0, molecular, 0.5),
                                   1.0e-14, "a ratio of one half should keep half the molecular viscosity"));
    PetscCall(PicurvAssertRealNear(-0.006,
                                   EddyViscosityFromCoefficient(-0.0005, 2.0, 3.0, molecular, 0.0),
                                   1.0e-14, "a mild backscatter should survive the floor untouched"));
    PetscFunctionReturn(0);
}

/** @brief Tests the Yoshizawa subgrid kinetic energy against its closed form. */
static PetscErrorCode TestSubgridKineticEnergy(void)
{
    PetscFunctionBeginUser;
    /* 2 * 0.09 * 4 * 9 = 6.48 */
    PetscCall(PicurvAssertRealNear(6.48, SubgridKineticEnergy(0.09, 2.0, 3.0), 1.0e-13,
                                   "the Yoshizawa relation should scale with Delta^2 |S|^2"));
    PetscFunctionReturn(0);
}

/*================================================================================*
 *                            FIXTURE-LEVEL TESTS                                 *
 *================================================================================*/

/**
 * @brief Fills a ghosted local scalar field over the whole local array.
 *
 * Writes the halo as well as the owned cells, so a test can poison the layout
 * boundary planes and observe whether an averaging routine counts them.
 */
static PetscErrorCode FillLocalScalar(UserCtx *user, Vec local, PetscReal value)
{
    PetscFunctionBeginUser;
    PetscCall(VecSet(local, value));
    (void)user;
    PetscFunctionReturn(0);
}

/**
 * @brief Writes a scalar into every cell the closure treats as owned and interior.
 */
static PetscErrorCode SetInteriorScalar(UserCtx *user, Vec local, PetscReal value)
{
    DMDALocalInfo info;
    PetscReal  ***array = NULL;
    PetscInt      xs, xe, ys, ye, zs, ze;

    PetscFunctionBeginUser;
    PetscCall(DMDAGetLocalInfo(user->da, &info));
    xs = (info.xs == 0) ? 1 : info.xs;
    ys = (info.ys == 0) ? 1 : info.ys;
    zs = (info.zs == 0) ? 1 : info.zs;
    xe = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    ye = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    ze = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    PetscCall(DMDAVecGetArray(user->da, local, &array));
    for (PetscInt k = zs; k < ze; ++k)
    for (PetscInt j = ys; j < ye; ++j)
    for (PetscInt i = xs; i < xe; ++i) array[k][j][i] = value;
    PetscCall(DMDAVecRestoreArray(user->da, local, &array));
    PetscFunctionReturn(0);
}

/** @brief Reads one cell of a ghosted local scalar field. */
static PetscErrorCode ReadLocalScalar(UserCtx *user, Vec local, PetscInt i, PetscInt j, PetscInt k,
                                      PetscReal *value)
{
    PetscReal ***array = NULL;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArrayRead(user->da, local, &array));
    *value = array[k][j][i];
    PetscCall(DMDAVecRestoreArrayRead(user->da, local, &array));
    PetscFunctionReturn(0);
}

/**
 * @brief Declares periodic boundary pairs on the axes a test wants homogeneous.
 *
 * The shared fixture builds a periodic DMDA from the `SimCtx` flags but leaves the
 * boundary face configuration untouched. The periodic field synchronization reads the
 * faces rather than the flags, so a test that synchronizes has to declare them. The
 * closure's own periodicity questions are answered from the flags, so tests that only
 * ask those do not call this.
 */
static void DeclarePeriodicFaces(UserCtx *user, PetscBool xi, PetscBool eta, PetscBool zeta)
{
    const BCFace    negative[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace    positive[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};
    const PetscBool selected[3] = {xi, eta, zeta};

    for (PetscInt axis = 0; axis < 3; ++axis) {
        if (!selected[axis]) continue;
        user->boundary_faces[negative[axis]].mathematical_type = PERIODIC;
        user->boundary_faces[positive[axis]].mathematical_type = PERIODIC;
    }
}

/** @brief Tests that an empty direction set reproduces the pointwise local model. */
static PetscErrorCode TestAverageRatioLocalIsPointwise(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Vec          numerator, denominator, ratio;
    SpatialTargetPlan plan;
    PetscReal ***num = NULL;
    PetscReal    value = 0.0;
    const PetscBool none[3] = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->lNvert, 0.0));
    PetscCall(VecSet(user->lAj, 1.0));

    PetscCall(VecDuplicate(user->lNvert, &numerator));
    PetscCall(VecDuplicate(user->lNvert, &denominator));
    PetscCall(VecDuplicate(user->lNvert, &ratio));
    PetscCall(FillLocalScalar(user, denominator, 4.0));
    PetscCall(FillLocalScalar(user, ratio, 0.0));
    PetscCall(FillLocalScalar(user, numerator, 0.0));

    /* A numerator that varies cell by cell, so a pointwise result cannot be confused
       with an averaged one. */
    PetscCall(DMDAVecGetArray(user->da, numerator, &num));
    for (PetscInt k = 1; k < 6; ++k)
    for (PetscInt j = 1; j < 6; ++j)
    for (PetscInt i = 1; i < 6; ++i) num[k][j][i] = (PetscReal)(i + j + k);
    PetscCall(DMDAVecRestoreArray(user->da, numerator, &num));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvSpatialRatioAverage(user, &plan, numerator, denominator, NULL, none,
                                        PETSC_COMM_WORLD, ratio, NULL));

    PetscCall(ReadLocalScalar(user, ratio, 2, 3, 4, &value));
    PetscCall(PicurvAssertRealNear((2.0 + 3.0 + 4.0) / 4.0, value, 1.0e-13,
                                   "with no averaging direction the quotient should stay pointwise"));

    PetscCall(VecDestroy(&numerator));
    PetscCall(VecDestroy(&denominator));
    PetscCall(VecDestroy(&ratio));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that averaging divides summed numerators by summed denominators.
 *
 * The distinction matters: the mean of the quotients differs from the quotient of the
 * means, and Lilly's least-squares closure calls for the latter. Separating them needs
 * a denominator that varies too. Three cells carry 1/1 and two carry 6/2, so every
 * pointwise quotient is 1 or 3 and their mean is 9/5, while the ratio of the sums is
 * 15/7. Only the second answer can come from averaging the two fields first.
 */
static PetscErrorCode TestAverageRatioDividesSummedFields(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Vec          numerator, denominator, ratio;
    SpatialTargetPlan plan;
    PetscReal ***num = NULL, ***den = NULL;
    PetscReal    value = 0.0;
    const PetscBool all[3] = {PETSC_TRUE, PETSC_TRUE, PETSC_TRUE};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->lNvert, 0.0));
    PetscCall(VecSet(user->lAj, 1.0));

    PetscCall(VecDuplicate(user->lNvert, &numerator));
    PetscCall(VecDuplicate(user->lNvert, &denominator));
    PetscCall(VecDuplicate(user->lNvert, &ratio));
    PetscCall(FillLocalScalar(user, numerator, 0.0));
    PetscCall(FillLocalScalar(user, denominator, 0.0));
    PetscCall(FillLocalScalar(user, ratio, 0.0));

    PetscCall(DMDAVecGetArray(user->da, numerator, &num));
    PetscCall(DMDAVecGetArray(user->da, denominator, &den));
    for (PetscInt k = 1; k < 6; ++k)
    for (PetscInt j = 1; j < 6; ++j)
    for (PetscInt i = 1; i < 6; ++i) {
        const PetscBool low = (PetscBool)(i < 4);

        num[k][j][i] = low ? 1.0 : 6.0;
        den[k][j][i] = low ? 1.0 : 2.0;
    }
    PetscCall(DMDAVecRestoreArray(user->da, numerator, &num));
    PetscCall(DMDAVecRestoreArray(user->da, denominator, &den));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvSpatialRatioAverage(user, &plan, numerator, denominator, NULL, all,
                                        PETSC_COMM_WORLD, ratio, NULL));

    /* 15/7, the ratio of the sums. Averaging the pointwise quotients would give 9/5. */
    PetscCall(ReadLocalScalar(user, ratio, 2, 3, 4, &value));
    PetscCall(PicurvAssertRealNear(15.0 / 7.0, value, 1.0e-13,
                                   "the global average should divide summed numerators by summed "
                                   "denominators, not average the pointwise quotients"));
    PetscCall(ReadLocalScalar(user, ratio, 4, 2, 1, &value));
    PetscCall(PicurvAssertRealNear(15.0 / 7.0, value, 1.0e-13,
                                   "a global average should be the same number in every cell"));

    PetscCall(VecDestroy(&numerator));
    PetscCall(VecDestroy(&denominator));
    PetscCall(VecDestroy(&ratio));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that averaging over two directions leaves a profile along the third.
 */
static PetscErrorCode TestAverageRatioRetainsUnaveragedDirection(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Vec          numerator, denominator, ratio;
    SpatialTargetPlan plan;
    PetscReal ***num = NULL;
    PetscReal    at_low = 0.0, at_high = 0.0;
    /* Average over xi and zeta, retain eta: the channel-flow arrangement. */
    const PetscBool ik[3] = {PETSC_TRUE, PETSC_FALSE, PETSC_TRUE};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->lNvert, 0.0));
    PetscCall(VecSet(user->lAj, 1.0));

    PetscCall(VecDuplicate(user->lNvert, &numerator));
    PetscCall(VecDuplicate(user->lNvert, &denominator));
    PetscCall(VecDuplicate(user->lNvert, &ratio));
    PetscCall(FillLocalScalar(user, numerator, 0.0));
    PetscCall(FillLocalScalar(user, denominator, 0.0));
    PetscCall(FillLocalScalar(user, ratio, 0.0));
    PetscCall(SetInteriorScalar(user, denominator, 1.0));

    /* Depends on eta alone, so averaging over xi and zeta must leave it untouched. */
    PetscCall(DMDAVecGetArray(user->da, numerator, &num));
    for (PetscInt k = 1; k < 6; ++k)
    for (PetscInt j = 1; j < 6; ++j)
    for (PetscInt i = 1; i < 6; ++i) num[k][j][i] = (PetscReal)j;
    PetscCall(DMDAVecRestoreArray(user->da, numerator, &num));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvSpatialRatioAverage(user, &plan, numerator, denominator, NULL, ik,
                                        PETSC_COMM_WORLD, ratio, NULL));

    PetscCall(ReadLocalScalar(user, ratio, 2, 1, 3, &at_low));
    PetscCall(ReadLocalScalar(user, ratio, 4, 5, 2, &at_high));
    PetscCall(PicurvAssertRealNear(1.0, at_low, 1.0e-13,
                                   "averaging over xi and zeta should preserve the eta profile"));
    PetscCall(PicurvAssertRealNear(5.0, at_high, 1.0e-13,
                                   "each retained eta plane should keep its own average"));

    PetscCall(VecDestroy(&numerator));
    PetscCall(VecDestroy(&denominator));
    PetscCall(VecDestroy(&ratio));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that a spatial average ignores the periodic duplicate planes.
 *
 * On a periodic block, index 0 and index `m-1` hold copies of interior cells. Counting
 * them would weight those cells twice. The test poisons both planes and requires the
 * average not to move.
 */
static PetscErrorCode TestAverageRatioIgnoresPeriodicDuplicatePlanes(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Vec          numerator, denominator, ratio;
    SpatialTargetPlan plan;
    PetscReal ***num = NULL;
    DMDALocalInfo info;
    PetscReal    clean = 0.0, poisoned = 0.0;
    const PetscBool all[3] = {PETSC_TRUE, PETSC_TRUE, PETSC_TRUE};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 6, 6,
                                                         PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    DeclarePeriodicFaces(user, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE);
    PetscCall(DMDAGetLocalInfo(user->da, &info));
    PetscCall(VecSet(user->lNvert, 0.0));
    PetscCall(VecSet(user->lAj, 1.0));

    PetscCall(VecDuplicate(user->lNvert, &numerator));
    PetscCall(VecDuplicate(user->lNvert, &denominator));
    PetscCall(VecDuplicate(user->lNvert, &ratio));
    PetscCall(FillLocalScalar(user, numerator, 0.0));
    PetscCall(FillLocalScalar(user, denominator, 0.0));
    PetscCall(FillLocalScalar(user, ratio, 0.0));
    PetscCall(SetInteriorScalar(user, numerator, 2.0));
    PetscCall(SetInteriorScalar(user, denominator, 1.0));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvSpatialRatioAverage(user, &plan, numerator, denominator, NULL, all,
                                        PETSC_COMM_WORLD, ratio, NULL));
    PetscCall(ReadLocalScalar(user, ratio, 2, 2, 2, &clean));

    /* Poison both layout boundary planes in every direction. */
    PetscCall(DMDAVecGetArray(user->da, numerator, &num));
    for (PetscInt k = 0; k < info.mz; ++k)
    for (PetscInt j = 0; j < info.my; ++j)
    for (PetscInt i = 0; i < info.mx; ++i) {
        const PetscBool duplicate = (PetscBool)(i == 0 || j == 0 || k == 0 ||
                                                i == info.mx - 1 || j == info.my - 1 || k == info.mz - 1);

        if (duplicate) num[k][j][i] = 1.0e6;
    }
    PetscCall(DMDAVecRestoreArray(user->da, numerator, &num));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvSpatialRatioAverage(user, &plan, numerator, denominator, NULL, all,
                                        PETSC_COMM_WORLD, ratio, NULL));
    PetscCall(ReadLocalScalar(user, ratio, 2, 2, 2, &poisoned));

    PetscCall(PicurvAssertRealNear(2.0, clean, 1.0e-13,
                                   "the interior average should be the interior value"));
    PetscCall(PicurvAssertRealNear(clean, poisoned, 1.0e-13,
                                   "the duplicate planes must not contribute to a spatial average"));

    PetscCall(VecDestroy(&numerator));
    PetscCall(VecDestroy(&denominator));
    PetscCall(VecDestroy(&ratio));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Tests that homogeneous averaging falls back to the block's periodic axes. */
static PetscErrorCode TestHomogeneousAveragingDerivesPeriodicAxes(void)
{
    SimCtx   *simCtx = NULL;
    UserCtx  *user = NULL;
    PetscBool direction[3];

    PetscFunctionBeginUser;
    /* Periodic in xi and zeta only: the channel arrangement. The fixture sets the
       resolved periodicity flags, which are what the closure reads; the boundary faces
       are left alone here so the test cannot pass by consulting them instead. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 6, 6,
                                                         PETSC_TRUE, PETSC_FALSE, PETSC_TRUE));
    simCtx->les_config.averaging_mode = LES_AVERAGING_HOMOGENEOUS;

    PetscCall(ResolveLESAveragingDirections(user, direction));
    PetscCall(PicurvAssertBool(direction[0], "a periodic xi axis should be averaged over"));
    PetscCall(PicurvAssertBool((PetscBool)(!direction[1]),
                               "a wall-normal axis should keep its own coefficient profile"));
    PetscCall(PicurvAssertBool(direction[2], "a periodic zeta axis should be averaged over"));

    /* An explicit selection overrides the periodic default. */
    simCtx->les_config.averaging_direction[0] = PETSC_TRUE;
    simCtx->les_config.averaging_direction[1] = PETSC_FALSE;
    simCtx->les_config.averaging_direction[2] = PETSC_FALSE;
    PetscCall(ResolveLESAveragingDirections(user, direction));
    PetscCall(PicurvAssertBool(direction[0], "an explicit direction should be honoured"));
    PetscCall(PicurvAssertBool((PetscBool)(!direction[2]),
                               "an explicit selection should replace the periodic default"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Tests that local and global averaging ignore the configured direction list. */
static PetscErrorCode TestAveragingModesSelectTheirOwnDirections(void)
{
    SimCtx   *simCtx = NULL;
    UserCtx  *user = NULL;
    PetscBool direction[3];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    simCtx->les_config.averaging_direction[0] = PETSC_TRUE;

    simCtx->les_config.averaging_mode = LES_AVERAGING_LOCAL;
    PetscCall(ResolveLESAveragingDirections(user, direction));
    PetscCall(PicurvAssertBool((PetscBool)(!direction[0] && !direction[1] && !direction[2]),
                               "local averaging should select no direction at all"));

    simCtx->les_config.averaging_mode = LES_AVERAGING_GLOBAL;
    PetscCall(ResolveLESAveragingDirections(user, direction));
    PetscCall(PicurvAssertBool((PetscBool)(direction[0] && direction[1] && direction[2]),
                               "global averaging should select every direction"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the constant model builds its viscosity without a coefficient field.
 *
 * `UserCtx::CS` and `UserCtx::lCs` are deliberately left unallocated. If the constant
 * path ever reaches for them again this case crashes rather than silently reintroducing
 * a field of one repeated number.
 */
static PetscErrorCode TestConstantModelNeedsNoCoefficientField(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Cmpnts    ***ucat = NULL;
    PetscReal ***nu_t = NULL;
    const PetscReal constant_cs = 0.2;
    /* u = x gives |S| = sqrt(2); with Aj = 1 the filter width is 1. */
    const PetscReal expected = constant_cs * constant_cs * PetscSqrtReal(2.0);

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));
    simCtx->les = CONSTANT_SMAGORINSKY;
    simCtx->les_config.constant_cs = constant_cs;

    PetscCall(DMCreateGlobalVector(user->da, &user->Nu_t));
    PetscCall(DMCreateLocalVector(user->da, &user->lNu_t));
    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(VecSet(user->Nu_t, 0.0));

    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k)
    for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j)
    for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
        ucat[k][j][i].x = (PetscReal)i;
        ucat[k][j][i].y = 0.0;
        ucat[k][j][i].z = 0.0;
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));

    PetscCall(PicurvAssertBool((PetscBool)(user->lCs == NULL),
                               "the constant model should allocate no coefficient field"));
    PetscCall(ComputeEddyViscosityLES(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->Nu_t, &nu_t));
    PetscCall(PicurvAssertRealNear(expected, nu_t[2][2][2], 1.0e-9,
                                   "the constant model should read its coefficient from configuration"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Nu_t, &nu_t));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Tests that the dynamic procedure refuses to run for the constant model. */
static PetscErrorCode TestDynamicProcedureRejectsConstantModel(void)
{
    SimCtx        *simCtx = NULL;
    UserCtx       *user = NULL;
    PetscErrorCode call_ierr;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));
    simCtx->les = CONSTANT_SMAGORINSKY;
    simCtx->step = 5;
    simCtx->StartStep = 0;

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    call_ierr = ComputeSmagorinskyConstant(user);
    PetscCall(PetscPopErrorHandler());

    PetscCall(PicurvAssertBool((PetscBool)(call_ierr != 0),
                               "the dynamic procedure should reject a model that has no "
                               "coefficient field to fill"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the dynamic procedure returns a zero coefficient on uniform flow.
 *
 * Exercises the whole path: the strain precompute and its halo, the stencil gather,
 * every test filter, both contractions, the averaging, and the limiting. A uniform
 * field has no strain and no stress between the two filter widths, so the Germano
 * identity has nothing to fit and the coefficient must come out exactly zero.
 */
static PetscErrorCode TestDynamicProcedureVanishesOnUniformFlow(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Cmpnts    ***ucat = NULL;
    PetscReal ***coefficient = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    simCtx->les = DYNAMIC_SMAGORINSKY;
    simCtx->step = 5;
    simCtx->StartStep = 0;

    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    PetscCall(DMCreateLocalVector(user->da, &user->lCs));
    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(VecSet(user->lAj, 1.0));
    PetscCall(VecSet(user->lNvert, 0.0));

    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k)
    for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j)
    for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
        ucat[k][j][i].x = 1.25;
        ucat[k][j][i].y = -0.5;
        ucat[k][j][i].z = 3.0;
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));

    PetscCall(ComputeSmagorinskyConstant(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->CS, &coefficient));
    for (PetscInt k = 1; k < 6; ++k)
    for (PetscInt j = 1; j < 6; ++j)
    for (PetscInt i = 1; i < 6; ++i) {
        PetscCall(PicurvAssertRealNear(0.0, coefficient[k][j][i], 1.0e-20,
                                       "a uniform field should give the dynamic procedure "
                                       "nothing to fit"));
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->CS, &coefficient));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests the dynamic procedure end to end on a periodic sheared field.
 *
 * Runs the full procedure on a block that is periodic in all three directions, where
 * the strain precompute has to reach across the wrap correctly, and asks for global
 * averaging. Two things must hold: every cell receives the same finite coefficient,
 * because a global average is one number for the block; and the coefficient stays
 * inside the configured ceiling.
 */
static PetscErrorCode TestDynamicProcedureGlobalAverageIsUniform(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Cmpnts    ***ucat = NULL;
    PetscReal ***coefficient = NULL;
    PetscReal    reference = 0.0;
    const FieldId cell_fields[] = {FIELD_ID_UCAT};
    const PetscReal wavenumber = 2.0 * PETSC_PI / 5.0; /* five distinct cells per direction */

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 6, 6,
                                                         PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    DeclarePeriodicFaces(user, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE);
    simCtx->les = DYNAMIC_SMAGORINSKY;
    simCtx->step = 5;
    simCtx->StartStep = 0;
    simCtx->les_config.averaging_mode = LES_AVERAGING_GLOBAL;

    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    PetscCall(DMCreateLocalVector(user->da, &user->lCs));
    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(VecSet(user->lAj, 1.0));
    PetscCall(VecSet(user->lNvert, 0.0));

    /* A divergence-free periodic shear: enough structure that both contractions are
       nonzero, and smooth enough that the coefficient stays in range. */
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k)
    for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j)
    for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
        ucat[k][j][i].x =  PetscSinReal(wavenumber * (PetscReal)j);
        ucat[k][j][i].y =  PetscSinReal(wavenumber * (PetscReal)k);
        ucat[k][j][i].z =  PetscSinReal(wavenumber * (PetscReal)i);
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));

    PetscCall(ComputeSmagorinskyConstant(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->CS, &coefficient));
    reference = coefficient[1][1][1];
    PetscCall(PicurvAssertBool((PetscBool)(PetscIsNormalReal(reference) || reference == 0.0),
                               "the dynamic coefficient must be finite"));
    PetscCall(PicurvAssertBool((PetscBool)(reference >= 0.0 &&
                                           reference <= simCtx->les_config.max_cs *
                                                        simCtx->les_config.max_cs),
                               "clamping should keep the coefficient inside its ceiling"));
    for (PetscInt k = 1; k < 6; ++k)
    for (PetscInt j = 1; j < 6; ++j)
    for (PetscInt i = 1; i < 6; ++i) {
        PetscCall(PicurvAssertRealNear(reference, coefficient[k][j][i], 1.0e-13,
                                       "global averaging should give the block one coefficient"));
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->CS, &coefficient));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Entry point for the LES closure suite.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"sym-tensor-algebra", TestSymTensorAlgebra},
        {"strain-rate-from-gradients", TestStrainRateFromGradients},
        {"filter-width-models-separate-on-stretched-cell", TestFilterWidthModelsSeparateOnStretchedCell},
        {"leonard-stress-vanishes-on-uniform-flow", TestLeonardStressVanishesOnUniformFlow},
        {"germano-model-tensor-on-constant-strain", TestGermanoModelTensorOnConstantStrain},
        {"germano-model-tensor-uses-filtered-product", TestGermanoModelTensorUsesFilteredProduct},
        {"clip-model-coefficient-modes", TestClipModelCoefficientModes},
        {"eddy-viscosity-floor-bounds-total-viscosity", TestEddyViscosityFloorBoundsTotalViscosity},
        {"subgrid-kinetic-energy", TestSubgridKineticEnergy},
        {"average-ratio-local-is-pointwise", TestAverageRatioLocalIsPointwise},
        {"average-ratio-divides-summed-fields", TestAverageRatioDividesSummedFields},
        {"average-ratio-retains-unaveraged-direction", TestAverageRatioRetainsUnaveragedDirection},
        {"average-ratio-ignores-periodic-duplicate-planes", TestAverageRatioIgnoresPeriodicDuplicatePlanes},
        {"homogeneous-averaging-derives-periodic-axes", TestHomogeneousAveragingDerivesPeriodicAxes},
        {"averaging-modes-select-their-own-directions", TestAveragingModesSelectTheirOwnDirections},
        {"constant-model-needs-no-coefficient-field", TestConstantModelNeedsNoCoefficientField},
        {"dynamic-procedure-rejects-constant-model", TestDynamicProcedureRejectsConstantModel},
        {"dynamic-procedure-vanishes-on-uniform-flow", TestDynamicProcedureVanishesOnUniformFlow},
        {"dynamic-procedure-global-average-is-uniform", TestDynamicProcedureGlobalAverageIsUniform},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv LES closure tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-les", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
