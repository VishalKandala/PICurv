/**
 * @file test_postprocessing.c
 * @brief C unit tests for post-processing kernel functions.
 */

#include "test_support.h"

#include "postprocessing_kernels.h"
/**
 * @brief Tests specific turbulent kinetic energy computation on particle data.
 */

static PetscErrorCode TestComputeSpecificKE(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal (*vel_arr)[3] = NULL;
    PetscScalar *ske_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    vel_arr[0][0] = 1.0;
    vel_arr[0][1] = 2.0;
    vel_arr[0][2] = 2.0;
    vel_arr[1][0] = 0.0;
    vel_arr[1][1] = 3.0;
    vel_arr[1][2] = 4.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(ComputeSpecificKE(user, "velocity", "ske"));
    PetscCall(DMSwarmGetField(user->post_swarm, "ske", NULL, NULL, (void *)&ske_arr));
    PetscCall(PicurvAssertRealNear(4.5, ske_arr[0], 1.0e-12, "ComputeSpecificKE first particle"));
    PetscCall(PicurvAssertRealNear(12.5, ske_arr[1], 1.0e-12, "ComputeSpecificKE second particle"));
    PetscCall(DMSwarmRestoreField(user->post_swarm, "ske", NULL, NULL, (void *)&ske_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle displacement computation against reference positions.
 */

static PetscErrorCode TestComputeDisplacement(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal (*pos_arr)[3] = NULL;
    PetscScalar *disp_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "disp"));

    simCtx->psrc_x = 1.0;
    simCtx->psrc_y = 2.0;
    simCtx->psrc_z = 3.0;

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 1.0; pos_arr[0][1] = 2.0; pos_arr[0][2] = 3.0;
    pos_arr[1][0] = 4.0; pos_arr[1][1] = 6.0; pos_arr[1][2] = 3.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(ComputeDisplacement(user, "disp"));
    PetscCall(DMSwarmGetField(user->post_swarm, "disp", NULL, NULL, (void *)&disp_arr));
    PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(disp_arr[0]), 1.0e-12, "ComputeDisplacement first particle"));
    PetscCall(PicurvAssertRealNear(5.0, PetscRealPart(disp_arr[1]), 1.0e-12, "ComputeDisplacement second particle"));
    PetscCall(DMSwarmRestoreField(user->post_swarm, "disp", NULL, NULL, (void *)&disp_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests nodal averaging of scalar cell-centered data.
 */

static PetscErrorCode TestComputeNodalAverageScalar(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    const PetscScalar ***p_nodal_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->P, 7.0));
    PetscCall(VecSet(user->P_nodal, -1.0));

    PetscCall(ComputeNodalAverage(user, "P", "P_nodal"));

    PetscCall(DMDAVecGetArrayRead(user->da, user->P_nodal, (void *)&p_nodal_arr));
    PetscCall(PicurvAssertRealNear(7.0, PetscRealPart(p_nodal_arr[0][0][0]), 1.0e-12, "ComputeNodalAverage interior node"));
    PetscCall(PicurvAssertRealNear(-1.0, PetscRealPart(p_nodal_arr[user->KM][user->JM][user->IM]), 1.0e-12,
                                   "ComputeNodalAverage untouched non-physical boundary node"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P_nodal, (void *)&p_nodal_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the cell-centred Q-criterion can be averaged to nodes.
 *
 * Qcrit is computed at cell centres and a .vts carries point data only, so writing it
 * directly placed every value half a cell from its node. It now goes through the same
 * nodal average as P, which needs it catalogued for the ghost refresh.
 */
static PetscErrorCode TestComputeNodalAverageQcrit(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    const PetscScalar ***q_nodal_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->Qcrit, 3.0));
    PetscCall(VecSet(user->Qcrit_nodal, -1.0));

    PetscCall(ComputeNodalAverage(user, "Qcrit", "Qcrit_nodal"));

    PetscCall(DMDAVecGetArrayRead(user->da, user->Qcrit_nodal, (void *)&q_nodal_arr));
    PetscCall(PicurvAssertRealNear(3.0, PetscRealPart(q_nodal_arr[0][0][0]), 1.0e-12,
                                   "ComputeNodalAverage should average Qcrit onto nodes"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Qcrit_nodal, (void *)&q_nodal_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests normalization of one field relative to a reference field.
 */

static PetscErrorCode TestNormalizeRelativeField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt ref_idx = 0;
    PetscScalar ref_value = 0.0;
    PetscReal vmin = 0.0, vmax = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscCalloc1(1, &simCtx->pps));
    simCtx->pps->reference[0] = 1;
    simCtx->pps->reference[1] = 1;
    simCtx->pps->reference[2] = 1;

    PetscCall(VecSet(user->P, 10.0));
    ref_idx = simCtx->pps->reference[2] * (user->info.mx * user->info.my) +
              simCtx->pps->reference[1] * user->info.mx +
              simCtx->pps->reference[0];
    PetscCall(VecSetValue(user->P, ref_idx, 4.0, INSERT_VALUES));
    PetscCall(VecAssemblyBegin(user->P));
    PetscCall(VecAssemblyEnd(user->P));

    PetscCall(NormalizeRelativeField(user, "P"));

    PetscCall(VecGetValues(user->P, 1, &ref_idx, &ref_value));
    PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(ref_value), 1.0e-12, "NormalizeRelativeField reference value"));
    PetscCall(VecMin(user->P, NULL, &vmin));
    PetscCall(VecMax(user->P, NULL, &vmax));
    PetscCall(PicurvAssertRealNear(0.0, vmin, 1.0e-12, "NormalizeRelativeField minimum"));
    PetscCall(PicurvAssertRealNear(6.0, vmax, 1.0e-12, "NormalizeRelativeField maximum"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests dimensionalization of pressure data from nondimensional values.
 */

static PetscErrorCode TestDimensionalizePressureField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->scaling.P_ref = 3.0;

    PetscCall(VecSet(user->P, 2.0));
    PetscCall(DimensionalizeField(user, "P"));
    PetscCall(PicurvAssertVecConstant(user->P, 6.0, 1.0e-12, "DimensionalizeField should scale pressure by P_ref"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the per-step dimensionalization leaves the persistent grid alone.
 * @details The postprocessor runs DimensionalizeAllLoadedFields once per processed step.
 *          Scaling the coordinates there multiplied them by L_ref once per step; they are
 *          loaded once and must be scaled once, by the postprocessor loop.
 */
static PetscErrorCode TestDimensionalizeAllLoadedFieldsLeavesCoordinatesAlone(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec coordinates = NULL, before = NULL;
    PetscReal difference = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->scaling.L_ref = 2.0;
    simCtx->scaling.U_ref = 3.0;
    simCtx->scaling.P_ref = 5.0;
    PetscCall(DMGetCoordinates(user->da, &coordinates));
    PetscCall(VecDuplicate(coordinates, &before));
    PetscCall(VecCopy(coordinates, before));
    PetscCall(VecSet(user->Ucat, 1.0));

    PetscCall(DimensionalizeAllLoadedFields(user));
    PetscCall(DimensionalizeAllLoadedFields(user));

    PetscCall(VecAXPY(before, -1.0, coordinates));
    PetscCall(VecNorm(before, NORM_INFINITY, &difference));
    PetscCall(PicurvAssertRealNear(0.0, difference, 1.0e-14,
                                   "per-step dimensionalization must not rescale the grid coordinates"));
    PetscCall(PicurvAssertVecConstant(user->Ucat, 9.0, 1.0e-12,
                                      "each call scales the reloaded velocity by U_ref"));

    PetscCall(VecDestroy(&before));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that a dimensionalized Q-criterion carries (U_ref/L_ref)^2, not U_ref^2.
 * @details Ucat is already dimensional when Q is computed, but the metrics are not, so
 *          the result must be divided by L_ref^2.
 */
static PetscErrorCode TestComputeQCriterionDimensionalizedScalesByLengthSquared(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucat = NULL;
    Vec nondimensional = NULL;
    PostProcessParams *pps = NULL;
    PetscReal difference = 0.0, reference = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++)
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++)
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++)
                ucat[k][j][i] = (Cmpnts){-(PetscReal)j, (PetscReal)i, 0.0}; /* solid-body rotation */
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));

    PetscCall(ComputeQCriterion(user));
    PetscCall(VecDuplicate(user->Qcrit, &nondimensional));
    PetscCall(VecCopy(user->Qcrit, nondimensional));
    PetscCall(VecNorm(nondimensional, NORM_INFINITY, &reference));
    PetscCall(PicurvAssertBool((PetscBool)(reference > 0.1), "a rotating field should have a positive Q"));

    PetscCall(PetscNew(&pps));
    pps->dimensionalize = PETSC_TRUE;
    simCtx->pps = pps;
    simCtx->scaling.L_ref = 2.0;
    PetscCall(ComputeQCriterion(user));
    PetscCall(VecAXPY(nondimensional, -4.0, user->Qcrit));
    PetscCall(VecNorm(nondimensional, NORM_INFINITY, &difference));
    PetscCall(PicurvAssertRealNear(0.0, difference, 1.0e-12 * reference,
                                   "dimensionalized Q must be the nondimensional Q over L_ref^2"));

    simCtx->pps = NULL;
    PetscCall(PetscFree(pps));
    PetscCall(VecDestroy(&nondimensional));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests Q-criterion computation for a quiescent velocity field.
 */

static PetscErrorCode TestComputeQCriterionZeroFlow(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));

    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(ComputeQCriterion(user));
    PetscCall(PicurvAssertVecConstant(user->Qcrit, 0.0, 1.0e-12, "ComputeQCriterion should be zero for uniform flow"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-post PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"compute-specific-ke", TestComputeSpecificKE},
        {"compute-displacement", TestComputeDisplacement},
        {"compute-nodal-average-scalar", TestComputeNodalAverageScalar},
        {"compute-nodal-average-qcrit", TestComputeNodalAverageQcrit},
        {"normalize-relative-field", TestNormalizeRelativeField},
        {"dimensionalize-pressure-field", TestDimensionalizePressureField},
        {"compute-qcriterion-zero-flow", TestComputeQCriterionZeroFlow},
        {"dimensionalize-all-loaded-fields-leaves-coordinates-alone", TestDimensionalizeAllLoadedFieldsLeavesCoordinatesAlone},
        {"compute-qcriterion-dimensionalized-scales-by-length-squared", TestComputeQCriterionDimensionalizedScalesByLengthSquared},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv post-processing tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-post", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
