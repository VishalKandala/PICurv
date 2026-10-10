/**
 * @file test_poisson_rhs.c
 * @brief C unit tests for Poisson, RHS, body-force, and diffusivity helpers.
 */

#include "test_support.h"

#include "poisson.h"
#include "rhs.h"
#include "setup.h"
#include "verification_sources.h"
/**
 * @brief Allocates Poisson/RHS support vectors required by the tests.
 */

static PetscErrorCode EnsurePoissonAndRhsVectors(UserCtx *user)
{
    PetscFunctionBeginUser;
    if (!user->Phi) PetscCall(DMCreateGlobalVector(user->da, &user->Phi));
    if (!user->lPhi) PetscCall(DMCreateLocalVector(user->da, &user->lPhi));
    if (!user->Aj) PetscCall(DMCreateGlobalVector(user->da, &user->Aj));
    if (!user->lAj) PetscCall(DMCreateLocalVector(user->da, &user->lAj));
    if (!user->Diffusivity) PetscCall(DMCreateGlobalVector(user->da, &user->Diffusivity));
    if (!user->lDiffusivity) PetscCall(DMCreateLocalVector(user->da, &user->lDiffusivity));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that pressure updates add the correction potential.
 */

static PetscErrorCode TestUpdatePressureAddsPhi(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(EnsurePoissonAndRhsVectors(user));
    PetscCall(VecSet(user->P, 1.0));
    PetscCall(VecSet(user->Phi, 0.25));

    PetscCall(UpdatePressure(user));
    PetscCall(PicurvAssertVecConstant(user->P, 1.25, 1.0e-12, "UpdatePressure should add Phi into P"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the Poisson RHS is zero for zero-divergence velocity.
 */

static PetscErrorCode TestPoissonRHSZeroDivergence(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec B = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(EnsurePoissonAndRhsVectors(user));

    simCtx->dt = 0.5;
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));

    PetscCall(VecDuplicate(user->P, &B));
    PetscCall(ComputePoissonRHS(user, B));
    PetscCall(PicurvAssertVecConstant(B, 0.0, 1.0e-12, "zero velocity divergence should produce zero Poisson RHS"));
    PetscCall(PicurvAssertRealNear(0.0, simCtx->poissonSourceImbalance, 1.0e-12, "global Poisson RHS sum should be zero"));

    PetscCall(VecDestroy(&B));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the body-force dispatcher across supported source terms.
 */

static PetscErrorCode TestComputeBodyForcesDispatcher(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec rct = NULL;
    Cmpnts ***rct_arr = NULL;
    Cmpnts ***l_csi = NULL;
    Cmpnts ***l_eta = NULL;
    Cmpnts ***l_zet = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecDuplicate(user->Ucont, &rct));
    PetscCall(VecZeroEntries(rct));

    user->boundary_faces[BC_FACE_NEG_X].face_id = BC_FACE_NEG_X;
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    user->boundary_faces[BC_FACE_POS_X].face_id = BC_FACE_POS_X;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    simCtx->bulkVelocityCorrection = 2.0;
    simCtx->dt = 1.0;
    simCtx->forceScalingFactor = 1.0;
    simCtx->drivingForceMagnitude = 0.0;
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMDAVecGetArray(user->fda, user->lCsi, &l_csi));
    PetscCall(DMDAVecGetArray(user->fda, user->lEta, &l_eta));
    PetscCall(DMDAVecGetArray(user->fda, user->lZet, &l_zet));
    l_csi[1][1][1].x = 1.0; l_csi[1][1][1].y = 0.0; l_csi[1][1][1].z = 0.0;
    l_eta[1][1][1].x = 0.0; l_eta[1][1][1].y = 1.0; l_eta[1][1][1].z = 0.0;
    l_zet[1][1][1].x = 0.0; l_zet[1][1][1].y = 0.0; l_zet[1][1][1].z = 1.0;
    PetscCall(DMDAVecRestoreArray(user->fda, user->lZet, &l_zet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lEta, &l_eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lCsi, &l_csi));

    PetscCall(ComputeBodyForces(user, rct));
    PetscCall(DMDAVecGetArrayRead(user->fda, rct, &rct_arr));
    PetscCall(PicurvAssertRealNear(1.5, simCtx->drivingForceMagnitude, 1.0e-12,
                                   "ComputeBodyForces should update the driven-flow controller magnitude"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].y, 1.0e-12, "ComputeBodyForces should keep y unchanged"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].z, 1.0e-12, "ComputeBodyForces should keep z unchanged"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, rct, &rct_arr));

    PetscCall(VecDestroy(&rct));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests Eulerian diffusivity for the molecular-only configuration.
 */

static PetscErrorCode TestComputeEulerianDiffusivityMolecularOnly(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***diff_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(EnsurePoissonAndRhsVectors(user));

    simCtx->ren = 2.0;
    simCtx->schmidt_number = 4.0;
    simCtx->les = PETSC_FALSE;
    PetscCall(VecZeroEntries(user->Diffusivity));

    PetscCall(ComputeEulerianDiffusivity(user));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diff_arr));
    PetscCall(PicurvAssertRealNear(0.125, diff_arr[1][1][1], 1.0e-12, "molecular diffusivity interior value"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diff_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that convection vanishes for a quiescent field.
 */

static PetscErrorCode TestConvectionZeroField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec conv = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecDuplicate(user->lUcont, &conv));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));

    PetscCall(Convection(user, user->lUcont, user->lUcat, conv));
    PetscCall(PicurvAssertVecConstant(conv, 0.0, 1.0e-12, "Convection should vanish for a quiescent field"));

    PetscCall(VecDestroy(&conv));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that viscous terms vanish for a uniform field.
 */

static PetscErrorCode TestViscousUniformField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec visc = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecDuplicate(user->lUcont, &visc));
    PetscCall(VecSet(user->Ucont, 1.0));
    PetscCall(VecSet(user->Ucat, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));

    PetscCall(Viscous(user, user->lUcont, user->lUcat, visc));
    PetscCall(PicurvAssertVecConstant(visc, 0.0, 1.0e-12, "Viscous should vanish for a uniform field"));

    PetscCall(VecDestroy(&visc));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Rescales the minimal fixture's unit metrics to cubic cells of edge h.
 */
static PetscErrorCode ScaleMinimalMetricsToSpacing(UserCtx *user, PetscReal h)
{
    Vec area_vectors[12] = {user->Csi, user->Eta, user->Zet, user->ICsi, user->IEta, user->IZet,
                            user->JCsi, user->JEta, user->JZet, user->KCsi, user->KEta, user->KZet};
    Vec local_area[12]   = {user->lCsi, user->lEta, user->lZet, user->lICsi, user->lIEta, user->lIZet,
                            user->lJCsi, user->lJEta, user->lJZet, user->lKCsi, user->lKEta, user->lKZet};
    Vec jacobians[4]     = {user->Aj, user->IAj, user->JAj, user->KAj};
    Vec local_jac[4]     = {user->lAj, user->lIAj, user->lJAj, user->lKAj};

    PetscFunctionBeginUser;
    for (PetscInt v = 0; v < 12; ++v) {
        PetscCall(VecScale(area_vectors[v], h * h));   /* a face of a cube of edge h */
        PetscCall(DMGlobalToLocalBegin(user->fda, area_vectors[v], INSERT_VALUES, local_area[v]));
        PetscCall(DMGlobalToLocalEnd(user->fda, area_vectors[v], INSERT_VALUES, local_area[v]));
    }
    for (PetscInt v = 0; v < 4; ++v) {
        PetscCall(VecSet(jacobians[v], 1.0 / (h * h * h)));
        PetscCall(DMGlobalToLocalBegin(user->da, jacobians[v], INSERT_VALUES, local_jac[v]));
        PetscCall(DMGlobalToLocalEnd(user->da, jacobians[v], INSERT_VALUES, local_jac[v]));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Returns the Clark contribution over the molecular viscous one, per unit x, at one cell.
 *
 * @details The field is u_x = x^2 on cubic cells of edge h, so du/dx = 2x and the exact
 *          gradient-model stress tau_11 = (h^2/12)(du/dx)^2 varies along x. Both the
 *          model's and viscosity's contributions leave Viscous() through the same flux
 *          differences, so their ratio does not depend on how Viscous() normalizes its
 *          output.
 */
static PetscErrorCode ClarkToMolecularRatioPerX(PetscReal h, PetscInt probe, PetscReal *ratio_per_x)
{
    SimCtx   *simCtx = NULL;
    UserCtx  *user = NULL;
    Vec       molecular = NULL, with_clark = NULL;
    Cmpnts ***ucat = NULL, ***visc0 = NULL, ***visc1 = NULL;
    const PetscInt n = 10, mid = n / 2;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, n, n, n));
    PetscCall(ScaleMinimalMetricsToSpacing(user, h));
    simCtx->ren = 1.0;
    simCtx->les = NO_LES_MODEL;

    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k)
    for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j)
    for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
        const PetscReal x = (i - 0.5) * h;   /* cell centre */
        ucat[k][j][i].x = x * x; ucat[k][j][i].y = 0.0; ucat[k][j][i].z = 0.0;
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(VecSet(user->lUcont, 0.0));

    PetscCall(VecDuplicate(user->lUcont, &molecular));
    PetscCall(VecDuplicate(user->lUcont, &with_clark));
    simCtx->les_gradient_model = 0;
    PetscCall(Viscous(user, user->lUcont, user->lUcat, molecular));
    simCtx->les_gradient_model = 1;
    PetscCall(Viscous(user, user->lUcont, user->lUcat, with_clark));

    PetscCall(DMDAVecGetArrayRead(user->fda, molecular, &visc0));
    PetscCall(DMDAVecGetArrayRead(user->fda, with_clark, &visc1));
    {
        const PetscReal v_molecular = visc0[mid][mid][probe].x;
        const PetscReal v_clark     = visc1[mid][mid][probe].x - v_molecular;
        PetscCheck(PetscAbsReal(v_molecular) > 0.0, PETSC_COMM_SELF, PETSC_ERR_PLIB,
                   "The molecular viscous term vanished at the probe; the fixture is wrong.");
        *ratio_per_x = (v_clark / v_molecular) / ((probe - 0.5) * h);
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, with_clark, &visc1));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, molecular, &visc0));
    PetscCall(VecDestroy(&with_clark));
    PetscCall(VecDestroy(&molecular));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the Clark gradient model scales as the square of the filter width.
 *
 * @details The gradient-model stress is tau_ij = (Delta_k^2/12) du_i/dx_k du_j/dx_k, so
 *          for a fixed physical field its force grows as Delta^2 relative to the
 *          molecular viscous force. Halving the spacing must therefore divide the
 *          Clark-to-molecular ratio at a given x by four. Fixtures with unit cells hide
 *          any error in the power of Delta, which is why the spacing is varied here.
 */
static PetscErrorCode TestClarkGradientModelScalesWithFilterWidthSquared(void)
{
    PetscReal coarse = 0.0, fine = 0.0;

    PetscFunctionBeginUser;
    /* Probe cells chosen so the two spacings give comparable x without being on a
       boundary: cell 5 at h = 0.2 and cell 5 at h = 0.1. The ratio is per unit x. */
    PetscCall(ClarkToMolecularRatioPerX(0.2, 5, &coarse));
    PetscCall(ClarkToMolecularRatioPerX(0.1, 5, &fine));
    PetscCall(PicurvAssertBool((PetscBool)(coarse < 0.0 && fine < 0.0),
                               "the gradient model removes momentum where du/dx grows along x"));
    PetscCall(PicurvAssertRealNear(4.0, coarse / fine, 1.0e-8,
                                   "halving the filter width divides the gradient-model force by four"));
    /* Exact value: the molecular flux carries both halves of the symmetric stress,
       2 nu du/dx per unit area, so its difference across a cell is 4 nu h^3 for
       u = x^2; the model's is -(2/3) h^5 x. Their ratio per unit x is -h^2/(6 nu),
       which pins the 1/12 coefficient as well as the power of Delta. */
    PetscCall(PicurvAssertRealNear(-0.2 * 0.2 / 6.0, coarse, 1.0e-10,
                                   "gradient-model force matches (Delta^2/12) du/dx du/dx at h = 0.2"));
    PetscCall(PicurvAssertRealNear(-0.1 * 0.1 / 6.0, fine, 1.0e-10,
                                   "gradient-model force matches (Delta^2/12) du/dx du/dx at h = 0.1"));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the full RHS remains zero without forcing on a quiescent field.
 */

static PetscErrorCode TestComputeRHSZeroFieldNoForcing(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(ComputeRHS(user, user->Rhs));
    PetscCall(PicurvAssertVecConstant(user->Rhs, 0.0, 1.0e-12, "ComputeRHS should remain zero without forcing on a quiescent field"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests diffusivity-gradient computation on a constant field.
 */

static PetscErrorCode TestComputeEulerianDiffusivityGradientConstantField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Diffusivity, 0.75));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Diffusivity, INSERT_VALUES, user->lDiffusivity));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Diffusivity, INSERT_VALUES, user->lDiffusivity));

    PetscCall(ComputeEulerianDiffusivityGradient(user));
    PetscCall(PicurvAssertVecConstant(user->DiffusivityGradient, 0.0, 1.0e-12, "constant diffusivity should yield zero diffusivity gradient"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests verification-driven linear diffusivity override and its gradient.
 */

static PetscErrorCode TestComputeEulerianDiffusivityVerificationLinearX(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***diff_arr = NULL;
    Cmpnts ***grad_arr = NULL;
    Cmpnts ***cent = NULL;
    const PetscReal gamma0 = 0.5;
    const PetscReal slope_x = 0.25;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "analytical", sizeof(simCtx->eulerianSource)));
    simCtx->verificationDiffusivity.enabled = PETSC_TRUE;
    PetscCall(PetscStrncpy(simCtx->verificationDiffusivity.mode, "analytical", sizeof(simCtx->verificationDiffusivity.mode)));
    PetscCall(PetscStrncpy(simCtx->verificationDiffusivity.profile, "LINEAR_X", sizeof(simCtx->verificationDiffusivity.profile)));
    simCtx->verificationDiffusivity.gamma0 = gamma0;
    simCtx->verificationDiffusivity.slope_x = slope_x;

    PetscCall(PicurvAssertBool(VerificationDiffusivityOverrideActive(simCtx),
                               "verification diffusivity override should report as active"));
    PetscCall(ComputeEulerianDiffusivity(user));
    PetscCall(ComputeEulerianDiffusivityGradient(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diff_arr));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Cent, &cent));
    PetscCall(PicurvAssertRealNear(gamma0 + slope_x * cent[2][2][2].x, diff_arr[2][2][2], 1.0e-12,
                                   "verification override should populate the linear diffusivity field"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diff_arr));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->DiffusivityGradient, &grad_arr));
    PetscCall(PicurvAssertRealNear(slope_x / user->IM, grad_arr[2][2][2].x, 1.0e-12,
                                   "linear verification diffusivity should yield the current finite-difference x gradient"));
    PetscCall(PicurvAssertRealNear(0.0, grad_arr[2][2][2].y, 1.0e-12,
                                   "linear verification diffusivity should keep y gradient zero"));
    PetscCall(PicurvAssertRealNear(0.0, grad_arr[2][2][2].z, 1.0e-12,
                                   "linear verification diffusivity should keep z gradient zero"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->DiffusivityGradient, &grad_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that Poisson matrix assembly produces a populated operator on a tiny Cartesian grid.
 */
static PetscErrorCode TestAssemblePoissonOperatorPopulatesRows(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt rows = 0, cols = 0;
    PetscInt ncols = 0;
    const PetscInt *col_idx = NULL;
    const PetscScalar *values = NULL;
    PetscInt interior_row = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(AssemblePoissonOperator(user));
    PetscCall(PicurvAssertBool((PetscBool)(user->A != NULL), "AssemblePoissonOperator should allocate the Poisson operator"));

    PetscCall(MatGetSize(user->A, &rows, &cols));
    PetscCall(PicurvAssertIntEqual(user->info.mx * user->info.my * user->info.mz, rows, "Poisson operator row count should match the DA node count"));
    PetscCall(PicurvAssertIntEqual(rows, cols, "Poisson operator should be square"));

    interior_row = (2 * user->info.my + 2) * user->info.mx + 2;
    PetscCall(MatGetRow(user->A, interior_row, &ncols, &col_idx, &values));
    PetscCall(PicurvAssertBool((PetscBool)(ncols > 1), "Interior Poisson row should couple to neighboring nodes"));
    PetscCall(PicurvAssertBool((PetscBool)(values != NULL), "Interior Poisson row should expose non-null coefficients"));
    PetscCall(MatRestoreRow(user->A, interior_row, &ncols, &col_idx, &values));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Fills every face metric and Jacobian of the minimal fixture with reproducible
 *        random values, so each face gradient carries all of its cross terms.
 */
static PetscErrorCode RandomizeMinimalFaceMetrics(UserCtx *user)
{
    Vec         face_metrics[9] = {user->ICsi, user->IEta, user->IZet, user->JCsi, user->JEta,
                                   user->JZet, user->KCsi, user->KEta, user->KZet};
    Vec         local_metrics[9] = {user->lICsi, user->lIEta, user->lIZet, user->lJCsi, user->lJEta,
                                    user->lJZet, user->lKCsi, user->lKEta, user->lKZet};
    Vec         jacobians[4] = {user->Aj, user->IAj, user->JAj, user->KAj};
    Vec         local_jac[4] = {user->lAj, user->lIAj, user->lJAj, user->lKAj};
    PetscRandom rnd;

    PetscFunctionBeginUser;
    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rnd));
    PetscCall(PetscRandomSetInterval(rnd, 0.5, 1.5));
    PetscCall(PetscRandomSetSeed(rnd, 2026));
    PetscCall(PetscRandomSeed(rnd));
    for (PetscInt v = 0; v < 9; ++v) {
        PetscCall(VecSetRandom(face_metrics[v], rnd));
        PetscCall(DMGlobalToLocalBegin(user->fda, face_metrics[v], INSERT_VALUES, local_metrics[v]));
        PetscCall(DMGlobalToLocalEnd(user->fda, face_metrics[v], INSERT_VALUES, local_metrics[v]));
    }
    for (PetscInt v = 0; v < 4; ++v) {
        PetscCall(VecSetRandom(jacobians[v], rnd));
        PetscCall(DMGlobalToLocalBegin(user->da, jacobians[v], INSERT_VALUES, local_jac[v]));
        PetscCall(DMGlobalToLocalEnd(user->da, jacobians[v], INSERT_VALUES, local_jac[v]));
    }
    PetscCall(PetscRandomDestroy(&rnd));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the operator is exactly the divergence of the projection's gradient.
 * @details With random non-orthogonal metrics and one solid cell, projecting a zero flux
 *          with a random Phi and forming the right-hand side of the result must reproduce
 *          -A Phi on every fluid row. The identity holds for any metric values only if the
 *          operator and the projection use the same face gradient, including its one-sided
 *          forms beside walls and solid cells.
 */
static PetscErrorCode TestOperatorAndProjectionShareOneFaceGradient(void)
{
    SimCtx            *simCtx = NULL;
    UserCtx           *user = NULL;
    Vec                B = NULL, APhi = NULL;
    PetscRandom        rnd;
    PetscReal       ***nvert;
    const PetscReal ***b, ***aphi;
    PetscReal          worst = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 7, 7, 7));
    PetscCall(EnsurePoissonAndRhsVectors(user));
    simCtx->dt = 0.3;
    PetscCall(RandomizeMinimalFaceMetrics(user));

    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    if (user->info.xs <= 3 && 3 < user->info.xs + user->info.xm &&
        user->info.ys <= 3 && 3 < user->info.ys + user->info.ym &&
        user->info.zs <= 3 && 3 < user->info.zs + user->info.zm) nvert[3][3][3] = 1.0;
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rnd));
    PetscCall(PetscRandomSetSeed(rnd, 7));
    PetscCall(PetscRandomSeed(rnd));
    PetscCall(VecSetRandom(user->Phi, rnd));
    PetscCall(PetscRandomDestroy(&rnd));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));

    PetscCall(AssemblePoissonOperator(user));
    PetscCall(ProjectVelocity(user));
    PetscCall(VecDuplicate(user->P, &B));
    PetscCall(ComputePoissonRHS(user, B));
    PetscCall(VecDuplicate(user->Phi, &APhi));
    PetscCall(MatMult(user->A, user->Phi, APhi));

    PetscCall(DMDAVecGetArrayRead(user->da, B, &b));
    PetscCall(DMDAVecGetArrayRead(user->da, APhi, &aphi));
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    for (PetscInt k = PetscMax(user->info.zs, 1); k < PetscMin(user->info.zs + user->info.zm, user->info.mz - 1); k++) {
        for (PetscInt j = PetscMax(user->info.ys, 1); j < PetscMin(user->info.ys + user->info.ym, user->info.my - 1); j++) {
            for (PetscInt i = PetscMax(user->info.xs, 1); i < PetscMin(user->info.xs + user->info.xm, user->info.mx - 1); i++) {
                if (nvert[k][j][i] > 0.1) continue;
                worst = PetscMax(worst, PetscAbsReal(b[k][j][i] + aphi[k][j][i]) / (1.0 + PetscAbsReal(aphi[k][j][i])));
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->da, APhi, &aphi));
    PetscCall(DMDAVecRestoreArrayRead(user->da, B, &b));
    PetscCall(PicurvAssertRealNear(0.0, worst, 1.0e-12,
                                   "the right-hand side of the projected flux should equal -A Phi on every fluid row"));

    PetscCall(VecDestroy(&APhi));
    PetscCall(VecDestroy(&B));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that projection leaves a zero pressure-correction field unchanged.
 */

static PetscErrorCode TestProjectionZeroPhiLeavesVelocityUnchanged(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Phi, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(ProjectVelocity(user));
    PetscCall(PicurvAssertVecConstant(user->Ucont, 0.0, 1.0e-12, "ProjectVelocity should leave a zero-velocity field unchanged when Phi is zero"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that projection applies the expected x-direction correction for a linear pressure field.
 */
static PetscErrorCode TestProjectionLinearPhiCorrectsVelocity(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***phi = NULL;
    Cmpnts ***ucont = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    simCtx->dt = COEF_TIME_ACCURACY;
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMDAVecGetArray(user->da, user->Phi, &phi));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                phi[k][j][i] = (PetscReal)i;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->Phi, &phi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(ProjectVelocity(user));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear(-1.0, ucont[2][2][2].x, 1.0e-10, "ProjectVelocity should subtract the x pressure gradient under identity metrics"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[2][2][2].y, 1.0e-10, "ProjectVelocity should leave the y component unchanged for an x-only gradient"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[2][2][2].z, 1.0e-10, "ProjectVelocity should leave the z component unchanged for an x-only gradient"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Adds a smooth divergent perturbation to every interior face flux.
 * @details Boundary faces are left alone, so the inflow still balances the outflow and
 *          the Neumann pressure problem stays solvable.
 */
static PetscErrorCode PerturbInteriorFaceFluxes(UserCtx *user, PetscReal amplitude)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***ucont;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = PetscMax(info.zs, 1); k < PetscMin(info.zs + info.zm, info.mz - 1); k++) {
        for (PetscInt j = PetscMax(info.ys, 1); j < PetscMin(info.ys + info.ym, info.my - 1); j++) {
            for (PetscInt i = PetscMax(info.xs, 1); i < PetscMin(info.xs + info.xm, info.mx - 1); i++) {
                if (i <= info.mx - 3) ucont[k][j][i].x += amplitude * PetscSinReal(1.3 * i + 0.7 * j + 0.4 * k);
                if (j <= info.my - 3) ucont[k][j][i].y += amplitude * PetscCosReal(0.9 * i + 1.1 * j + 0.3 * k);
                if (k <= info.mz - 3) ucont[k][j][i].z += amplitude * PetscSinReal(0.5 * i + 0.8 * j + 1.7 * k);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCONT));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the production multigrid Poisson solve and projection remove a divergence.
 * @details Runs PoissonSolver_Multigrid through the real setup path with a three-level hierarchy,
 *          then the same UpdatePressure/ProjectVelocity pair the time loop uses, and requires
 *          the projected field to be divergence-free to the solve tolerance.
 */
static PetscErrorCode TestPoissonSolverMultigridProjectsToDivergenceFree(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    PetscReal before;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(
        NULL, PETSC_FALSE,
        "-im 17\n-jm 17\n-km 17\n-mg_level 3\n"
        "-ps_ksp_rtol 1.0e-12\n-ps_ksp_atol 1.0e-14\n-ps_ksp_max_it 200\n",
        &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(PicurvAssertIntEqual(3, simCtx->usermg.mglevels, "fixture should build a three-level hierarchy"));

    PetscCall(PerturbInteriorFaceFluxes(user, 0.2));
    PetscCall(ComputeDivergence(user));
    before = simCtx->MaxDiv;
    PetscCall(PicurvAssertBool((PetscBool)(before > 1.0e-2), "perturbation should make the field divergent"));

    PetscCall(PoissonSolver_Multigrid(&simCtx->usermg));
    PetscCall(UpdatePressure(user));
    PetscCall(ProjectVelocity(user));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCONT));
    PetscCall(ComputeDivergence(user));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->MaxDiv < 1.0e-9 * before),
                               "multigrid solve plus projection should leave no divergence above the solve tolerance"));

    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Counts the occurrences of @p needle in a text file; a missing file counts zero.
 */
static PetscErrorCode CountFileOccurrences(const char *path, const char *needle, PetscInt *count)
{
    FILE *file;
    char  line[1024];

    PetscFunctionBeginUser;
    *count = 0;
    file = fopen(path, "r");
    if (!file) PetscFunctionReturn(0);
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, needle)) (*count)++;
    }
    fclose(file);
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the multigrid solver is built once and reused on later steps.
 * @details Two consecutive solves must use the same Krylov solver and operator, both must
 *          converge, and the convergence log must hold one header per step.
 */
static PetscErrorCode TestPoissonSolverMultigridReusesItsSolver(void)
{
    SimCtx            *simCtx = NULL;
    UserCtx           *user = NULL;
    char               tmpdir[PETSC_MAX_PATH_LEN], log_path[PETSC_MAX_PATH_LEN + 64];
    KSP                first_ksp;
    Mat                first_operator;
    KSPConvergedReason reason;
    PetscInt           headers = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(
        NULL, PETSC_FALSE,
        "-im 17\n-jm 17\n-km 17\n-mg_level 3\n"
        "-ps_ksp_rtol 1.0e-12\n-ps_ksp_atol 1.0e-14\n-ps_ksp_max_it 200\n",
        &simCtx, &user, tmpdir, sizeof(tmpdir)));

    simCtx->step = simCtx->StartStep + 1;
    PetscCall(PerturbInteriorFaceFluxes(user, 0.2));
    PetscCall(PoissonSolver_Multigrid(&simCtx->usermg));
    first_ksp = user->ksp;
    first_operator = user->A;
    PetscCall(PicurvAssertBool((PetscBool)(first_ksp != NULL), "the first solve should build and keep the solver"));
    PetscCall(KSPGetConvergedReason(user->ksp, &reason));
    PetscCall(PicurvAssertBool((PetscBool)(reason > 0), "the first solve should converge"));

    simCtx->step++;
    PetscCall(PerturbInteriorFaceFluxes(user, 0.1));
    PetscCall(PoissonSolver_Multigrid(&simCtx->usermg));
    PetscCall(PicurvAssertBool((PetscBool)(user->ksp == first_ksp), "the second solve should reuse the solver"));
    PetscCall(PicurvAssertBool((PetscBool)(user->A == first_operator), "the second solve should reuse the operator"));
    PetscCall(KSPGetConvergedReason(user->ksp, &reason));
    PetscCall(PicurvAssertBool((PetscBool)(reason > 0), "the reused solver should converge"));

    PetscCall(PetscSNPrintf(log_path, sizeof(log_path), "%s/Poisson_Solver_Convergence_History_Block_0.log", simCtx->log_dir));
    PetscCall(CountFileOccurrences(log_path, "--- Convergence for Timestep", &headers));
    PetscCall(PicurvAssertIntEqual(2, headers, "the convergence log should hold one header per solve"));

    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests the null space attached to the operator.
 * @details Removing it from a constant must leave zero everywhere, and from a linear field
 *          must leave zero mean over the interior with zero dummy layers.
 */
static PetscErrorCode TestPoissonNullSpaceRemovesTheInteriorMean(void)
{
    SimCtx       *simCtx = NULL;
    UserCtx      *user = NULL;
    char          tmpdir[PETSC_MAX_PATH_LEN];
    MatNullSpace  nullsp = NULL;
    Vec           x = NULL;
    PetscReal  ***xa, norm, sum;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(
        NULL, PETSC_FALSE, "-im 9\n-jm 9\n-km 9\n-mg_level 2\n", &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(PoissonSolver_Multigrid(&simCtx->usermg));
    PetscCall(MatGetNullSpace(user->A, &nullsp));
    PetscCall(PicurvAssertBool((PetscBool)(nullsp != NULL), "the operator should carry the Neumann null space"));

    PetscCall(VecDuplicate(user->Phi, &x));
    PetscCall(VecSet(x, 3.0));
    PetscCall(MatNullSpaceRemove(nullsp, x));
    PetscCall(VecNorm(x, NORM_INFINITY, &norm));
    PetscCall(PicurvAssertRealNear(0.0, norm, 1.0e-12, "removing the null space should annihilate a constant"));

    PetscCall(DMDAVecGetArray(user->da, x, &xa));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++)
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++)
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) xa[k][j][i] = (PetscReal)i;
    PetscCall(DMDAVecRestoreArray(user->da, x, &xa));
    PetscCall(MatNullSpaceRemove(nullsp, x));
    PetscCall(VecSum(x, &sum));
    PetscCall(PicurvAssertRealNear(0.0, sum, 1.0e-10, "the interior mean should vanish"));
    PetscCall(DMDAVecGetArray(user->da, x, &xa));
    if (user->info.xs == 0) PetscCall(PicurvAssertRealNear(0.0, xa[user->info.zs][user->info.ys][0], 1.0e-14, "dummy layers should be zero"));
    PetscCall(DMDAVecRestoreArray(user->da, x, &xa));

    PetscCall(VecDestroy(&x));
    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that an over-coarsened hierarchy stops the run at the Poisson solve.
 * @details Three levels on a 9-node grid leave a coarse operator the preconditioner cannot
 *          factor. The solve must return PETSC_ERR_NOT_CONVERGED rather than hand an
 *          unsolved Phi to the projection.
 */
static PetscErrorCode TestPoissonSolverMultigridRefusesAnOvercoarsenedHierarchy(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    PetscErrorCode solve_error;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(
        NULL, PETSC_FALSE,
        "-im 9\n-jm 9\n-km 9\n-mg_level 3\n"
        /* The shipped profiles solve the coarse level with a replicated direct factor. */
        "-ps_ksp_type fgmres\n-ps_pc_type mg\n"
        "-ps_mg_coarse_ksp_type preonly\n-ps_mg_coarse_pc_type redundant\n",
        &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(PerturbInteriorFaceFluxes(user, 0.2));

    PetscCall(PetscPushErrorHandler(PetscReturnErrorHandler, NULL));
    solve_error = PoissonSolver_Multigrid(&simCtx->usermg);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_NOT_CONVERGED, (PetscInt)solve_error,
                                   "an unfactorable coarse level should stop the solve, not reach the projection"));

    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-poisson-rhs PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"update-pressure-adds-phi", TestUpdatePressureAddsPhi},
        {"poisson-rhs-zero-divergence", TestPoissonRHSZeroDivergence},
        {"compute-body-forces-dispatcher", TestComputeBodyForcesDispatcher},
        {"compute-eulerian-diffusivity-molecular-only", TestComputeEulerianDiffusivityMolecularOnly},
        {"convection-zero-field", TestConvectionZeroField},
        {"viscous-uniform-field", TestViscousUniformField},
        {"clark-gradient-model-scales-with-filter-width-squared", TestClarkGradientModelScalesWithFilterWidthSquared},
        {"compute-rhs-zero-field-no-forcing", TestComputeRHSZeroFieldNoForcing},
        {"compute-eulerian-diffusivity-gradient-constant-field", TestComputeEulerianDiffusivityGradientConstantField},
        {"compute-eulerian-diffusivity-verification-linear-x", TestComputeEulerianDiffusivityVerificationLinearX},
        {"assemble-poisson-operator-populates-rows", TestAssemblePoissonOperatorPopulatesRows},
        {"operator-and-projection-share-one-face-gradient", TestOperatorAndProjectionShareOneFaceGradient},
        {"projection-zero-phi-leaves-velocity-unchanged", TestProjectionZeroPhiLeavesVelocityUnchanged},
        {"projection-linear-phi-corrects-velocity", TestProjectionLinearPhiCorrectsVelocity},
        {"poisson-solver-multigrid-projects-to-divergence-free", TestPoissonSolverMultigridProjectsToDivergenceFree},
        {"poisson-solver-multigrid-reuses-its-solver", TestPoissonSolverMultigridReusesItsSolver},
        {"poisson-null-space-removes-the-interior-mean", TestPoissonNullSpaceRemovesTheInteriorMean},
        {"poisson-solver-multigrid-refuses-an-overcoarsened-hierarchy", TestPoissonSolverMultigridRefusesAnOvercoarsenedHierarchy},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv Poisson/RHS tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-poisson-rhs", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
