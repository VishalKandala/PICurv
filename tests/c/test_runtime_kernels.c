/**
 * @file test_runtime_kernels.c
 * @brief C unit tests for runtime, particle, wall, and walltime-guard helpers.
 */

#include "test_support.h"

#include "BC_Handlers.h"
#include "ParticleMotion.h"
#include "ParticlePhysics.h"
#include "ParticleSwarm.h"
#include "interpolation.h"
#include "initialcondition.h"
#include "les.h"
#include "runloop.h"
#include "setup.h"
#include "wallfunction.h"
#include "walkingsearch.h"

/**
 * @brief Synchronizes the minimal runtime fixture's global fields into their persistent local ghosts.
 */
static PetscErrorCode SyncRuntimeFieldGhosts(UserCtx *user)
{
    PetscFunctionBeginUser;
    PetscCall(DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Psi, INSERT_VALUES, user->lPsi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Psi, INSERT_VALUES, user->lPsi));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Diffusivity, INSERT_VALUES, user->lDiffusivity));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Diffusivity, INSERT_VALUES, user->lDiffusivity));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->DiffusivityGradient, INSERT_VALUES, user->lDiffusivityGradient));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->DiffusivityGradient, INSERT_VALUES, user->lDiffusivityGradient));
    PetscFunctionReturn(0);
}

/**
 * @brief Seeds one localized swarm particle with the cell, position, weight, and status data used by runtime tests.
 */
static PetscErrorCode SeedSingleParticle(UserCtx *user,
                                         PetscInt ci,
                                         PetscInt cj,
                                         PetscInt ck,
                                         PetscReal x,
                                         PetscReal y,
                                         PetscReal z,
                                         PetscReal wx,
                                         PetscReal wy,
                                         PetscReal wz,
                                         PetscInt status_value)
{
    PetscReal *positions = NULL;
    PetscReal *weights = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;

    PetscFunctionBeginUser;
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    positions[0] = x;
    positions[1] = y;
    positions[2] = z;
    weights[0] = wx;
    weights[1] = wy;
    weights[2] = wz;
    cell_ids[0] = ci;
    cell_ids[1] = cj;
    cell_ids[2] = ck;
    status[0] = status_value;

    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle distribution remainder handling across ranks.
 */

static PetscErrorCode TestDistributeParticlesRemainderHandling(void)
{
    PetscInt particles_per_rank = -1;
    PetscInt remainder = -1;

    PetscFunctionBeginUser;
    PetscCall(DistributeParticles(10, 0, 3, &particles_per_rank, &remainder));
    PetscCall(PicurvAssertIntEqual(4, particles_per_rank, "rank 0 should receive one remainder particle"));
    PetscCall(PicurvAssertIntEqual(1, remainder, "remainder should be reported correctly"));

    PetscCall(DistributeParticles(10, 2, 3, &particles_per_rank, &remainder));
    PetscCall(PicurvAssertIntEqual(3, particles_per_rank, "last rank should receive base particle count"));
    PetscCall(PicurvAssertIntEqual(1, remainder, "remainder should remain unchanged across ranks"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests basic particle-inside-bounding-box classification cases.
 */

static PetscErrorCode TestIsParticleInsideBoundingBoxBasicCases(void)
{
    BoundingBox bbox;
    Particle particle;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&bbox, sizeof(bbox)));
    PetscCall(PetscMemzero(&particle, sizeof(particle)));

    bbox.min_coords.x = 0.0;
    bbox.min_coords.y = 0.0;
    bbox.min_coords.z = 0.0;
    bbox.max_coords.x = 1.0;
    bbox.max_coords.y = 2.0;
    bbox.max_coords.z = 3.0;

    particle.loc.x = 0.25;
    particle.loc.y = 1.0;
    particle.loc.z = 2.5;
    PetscCall(PicurvAssertBool(IsParticleInsideBoundingBox(&bbox, &particle), "particle should be inside bounding box"));

    particle.loc.x = 1.5;
    PetscCall(PicurvAssertBool((PetscBool)!IsParticleInsideBoundingBox(&bbox, &particle), "particle should be outside bounding box"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle weight updates against expected ratios.
 */

static PetscErrorCode TestUpdateParticleWeightsComputesExpectedRatios(void)
{
    Particle particle;
    PetscReal distances[NUM_FACES] = {1.0, 3.0, 2.0, 2.0, 4.0, 1.0};
    PetscReal clamped[NUM_FACES] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&particle, sizeof(particle)));
    PetscCall(UpdateParticleWeights(distances, &particle));
    PetscCall(PicurvAssertRealNear(0.25, particle.weights.x, 1.0e-12, "x interpolation weight"));
    PetscCall(PicurvAssertRealNear(0.5, particle.weights.y, 1.0e-12, "y interpolation weight"));
    PetscCall(PicurvAssertRealNear(0.2, particle.weights.z, 1.0e-12, "z interpolation weight"));

    PetscCall(UpdateParticleWeights(clamped, &particle));
    PetscCall(PicurvAssertRealNear(0.5, particle.weights.x, 1.0e-12, "clamped x weight remains centered"));
    PetscCall(PicurvAssertRealNear(0.5, particle.weights.y, 1.0e-12, "clamped y weight remains centered"));
    PetscCall(PicurvAssertRealNear(0.5, particle.weights.z, 1.0e-12, "clamped z weight remains centered"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle position updates without Brownian forcing.
 */

static PetscErrorCode TestUpdateParticlePositionWithoutBrownianContribution(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Particle particle;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->dt = 0.25;

    PetscCall(PetscMemzero(&particle, sizeof(particle)));
    particle.loc.x = 1.0;
    particle.loc.y = -2.0;
    particle.loc.z = 3.0;
    particle.vel.x = 0.5;
    particle.vel.y = -1.0;
    particle.vel.z = 2.0;
    particle.diffusivitygradient.x = 0.1;
    particle.diffusivitygradient.y = 0.2;
    particle.diffusivitygradient.z = -0.3;
    particle.diffusivity = 0.0;

    PetscCall(UpdateParticlePosition(user, &particle));
    PetscCall(PicurvAssertRealNear(1.15, particle.loc.x, 1.0e-12, "x position update"));
    PetscCall(PicurvAssertRealNear(-2.2, particle.loc.y, 1.0e-12, "y position update"));
    PetscCall(PicurvAssertRealNear(3.425, particle.loc.z, 1.0e-12, "z position update"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle position updates driven only by diffusivity-gradient drift.
 */

static PetscErrorCode TestUpdateParticlePositionDiffusivityGradientOnly(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Particle particle;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->dt = 0.5;

    PetscCall(PetscMemzero(&particle, sizeof(particle)));
    particle.loc.x = 0.25;
    particle.loc.y = 0.5;
    particle.loc.z = 0.75;
    particle.vel.x = 0.0;
    particle.vel.y = 0.0;
    particle.vel.z = 0.0;
    particle.diffusivitygradient.x = 0.2;
    particle.diffusivitygradient.y = -0.1;
    particle.diffusivitygradient.z = 0.3;
    particle.diffusivity = 0.0;

    PetscCall(UpdateParticlePosition(user, &particle));
    PetscCall(PicurvAssertRealNear(0.35, particle.loc.x, 1.0e-12, "drift-only x position update"));
    PetscCall(PicurvAssertRealNear(0.45, particle.loc.y, 1.0e-12, "drift-only y position update"));
    PetscCall(PicurvAssertRealNear(0.90, particle.loc.z, 1.0e-12, "drift-only z position update"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests IEM relaxation updates for particle-carried fields.
 */

static PetscErrorCode TestUpdateParticleFieldIEMRelaxation(void)
{
    PetscReal psi = 1.0;
    const PetscReal dt = 0.5;
    const PetscReal diffusivity = 0.2;
    const PetscReal mean_val = 3.0;
    const PetscReal cell_vol = 8.0;
    const PetscReal c_model = 2.0;
    PetscReal unchanged = 7.0;

    PetscFunctionBeginUser;
    PetscCall(UpdateParticleField("Psi", dt, &psi, diffusivity, mean_val, cell_vol, c_model));
    PetscCall(PicurvAssertRealNear(
        mean_val + (1.0 - mean_val) * PetscExpReal(-(c_model * diffusivity / PetscPowReal(cell_vol, 0.6666667)) * dt),
        psi,
        1.0e-12,
        "IEM update should match analytical relaxation"));

    PetscCall(UpdateParticleField("UnrelatedField", dt, &unchanged, diffusivity, mean_val, cell_vol, c_model));
    PetscCall(PicurvAssertRealNear(7.0, unchanged, 1.0e-12, "unknown field should remain unchanged"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that non-Ucont requests do not modify interior field initialization.
 */

static PetscErrorCode TestSetInitialInteriorFieldIgnoresNonUcontRequest(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->Ucont, 7.0));

    PetscCall(SetInitialInteriorField(user, "P"));
    PetscCall(PicurvAssertVecConstant(user->Ucont, 7.0, 1.0e-12, "non-Ucont request should not modify Ucont"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests constant-profile interior initialization on a Z-direction inlet.
 */

static PetscErrorCode TestSetInitialInteriorFieldConstantProfileOnZInlet(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucont = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));
    user->GridOrientation = 1;
    user->identifiedInletBCFace = BC_FACE_NEG_Z;
    simCtx->FieldInitialization = 1;
    simCtx->InitialConstantContra.x = 0.0;
    simCtx->InitialConstantContra.y = 0.0;
    simCtx->InitialConstantContra.z = 2.0;
    PetscCall(VecSet(user->Ucont, 0.0));

    PetscCall(SetInitialInteriorField(user, "Ucont"));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear(0.0, ucont[1][1][1].x, 1.0e-12, "x contravariant component remains zero"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[1][1][1].y, 1.0e-12, "y contravariant component remains zero"));
    PetscCall(PicurvAssertRealNear(2.0, ucont[1][1][1].z, 1.0e-12, "z contravariant component should match constant inlet magnitude"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[0][1][1].z, 1.0e-12, "boundary/ghost cell should remain unchanged"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests direct interpolation from Eulerian fields to one localized swarm particle.
 */

static PetscErrorCode TestInterpolateAllFieldsToSwarmConstantFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***grad = NULL;
    PetscReal *velocity = NULL;
    PetscReal *diffusivity = NULL;
    PetscReal *diffusivity_gradient = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));
    PetscCall(VecSet(user->Ucat, 2.0));
    PetscCall(VecSet(user->Diffusivity, 0.25));

    PetscCall(DMDAVecGetArray(user->fda, user->DiffusivityGradient, &grad));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                grad[k][j][i].x = 0.1;
                grad[k][j][i].y = 0.2;
                grad[k][j][i].z = 0.3;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->DiffusivityGradient, &grad));
    PetscCall(SyncRuntimeFieldGhosts(user));
    PetscCall(SeedSingleParticle(user, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, ACTIVE_AND_LOCATED));

    PetscCall(InterpolateAllFieldsToSwarm(user));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));
    PetscCall(DMSwarmGetField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmGetField(user->swarm, "DiffusivityGradient", NULL, NULL, (void **)&diffusivity_gradient));
    PetscCall(PicurvAssertRealNear(2.0, velocity[0], 1.0e-12, "Interpolated velocity x should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(2.0, velocity[1], 1.0e-12, "Interpolated velocity y should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(2.0, velocity[2], 1.0e-12, "Interpolated velocity z should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(0.25, diffusivity[0], 1.0e-12, "Interpolated scalar diffusivity should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(0.1, diffusivity_gradient[0], 1.0e-12, "Interpolated diffusivity-gradient x component"));
    PetscCall(PicurvAssertRealNear(0.2, diffusivity_gradient[1], 1.0e-12, "Interpolated diffusivity-gradient y component"));
    PetscCall(PicurvAssertRealNear(0.3, diffusivity_gradient[2], 1.0e-12, "Interpolated diffusivity-gradient z component"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DiffusivityGradient", NULL, NULL, (void **)&diffusivity_gradient));
    PetscCall(DMSwarmRestoreField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the corner-averaged (legacy) interpolation path on constant fields.
 */

static PetscErrorCode TestInterpolateAllFieldsToSwarmCornerAveragedConstantFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***grad = NULL;
    PetscReal *velocity = NULL;
    PetscReal *diffusivity = NULL;
    PetscReal *diffusivity_gradient = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->interpolationMethod = INTERP_CORNER_AVERAGED;
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));
    PetscCall(VecSet(user->Ucat, 2.0));
    PetscCall(VecSet(user->Diffusivity, 0.25));

    PetscCall(DMDAVecGetArray(user->fda, user->DiffusivityGradient, &grad));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                grad[k][j][i].x = 0.1;
                grad[k][j][i].y = 0.2;
                grad[k][j][i].z = 0.3;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->DiffusivityGradient, &grad));
    PetscCall(SyncRuntimeFieldGhosts(user));
    PetscCall(SeedSingleParticle(user, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, ACTIVE_AND_LOCATED));

    PetscCall(InterpolateAllFieldsToSwarm(user));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));
    PetscCall(DMSwarmGetField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmGetField(user->swarm, "DiffusivityGradient", NULL, NULL, (void **)&diffusivity_gradient));
    PetscCall(PicurvAssertRealNear(2.0, velocity[0], 1.0e-12, "CornerAveraged: interpolated velocity x should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(2.0, velocity[1], 1.0e-12, "CornerAveraged: interpolated velocity y should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(2.0, velocity[2], 1.0e-12, "CornerAveraged: interpolated velocity z should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(0.25, diffusivity[0], 1.0e-12, "CornerAveraged: interpolated scalar diffusivity should match constant Eulerian field"));
    PetscCall(PicurvAssertRealNear(0.1, diffusivity_gradient[0], 1.0e-12, "CornerAveraged: interpolated diffusivity-gradient x component"));
    PetscCall(PicurvAssertRealNear(0.2, diffusivity_gradient[1], 1.0e-12, "CornerAveraged: interpolated diffusivity-gradient y component"));
    PetscCall(PicurvAssertRealNear(0.3, diffusivity_gradient[2], 1.0e-12, "CornerAveraged: interpolated diffusivity-gradient z component"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DiffusivityGradient", NULL, NULL, (void **)&diffusivity_gradient));
    PetscCall(DMSwarmRestoreField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle-to-grid scattering using known cell occupancy and scalar values.
 */

static PetscErrorCode TestScatterAllParticleFieldsToEulerFieldsAveragesPsi(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt *cell_ids = NULL;
    PetscReal *psi = NULL;
    PetscReal ***psi_grid = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    cell_ids[0] = 0; cell_ids[1] = 0; cell_ids[2] = 0;
    cell_ids[3] = 0; cell_ids[4] = 0; cell_ids[5] = 0;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));

    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    psi[0] = 1.0;
    psi[1] = 3.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));

    PetscCall(ScatterAllParticleFieldsToEulerFields(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->Psi, &psi_grid));
    PetscCall(PicurvAssertRealNear(2.0, psi_grid[1][1][1], 1.0e-12, "Scatter should average particle Psi values into the owning cell"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Psi, &psi_grid));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle counting by geometric cell IDs using the production +1 storage shift.
 */

static PetscErrorCode TestCalculateParticleCountPerCellCountsGlobalCellIDs(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt *cell_ids = NULL;
    PetscReal ***counts = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 3, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    cell_ids[0] = 0; cell_ids[1] = 0; cell_ids[2] = 0;
    cell_ids[3] = 0; cell_ids[4] = 0; cell_ids[5] = 0;
    cell_ids[6] = 1; cell_ids[7] = 0; cell_ids[8] = 0;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));

    PetscCall(CalculateParticleCountPerCell(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->ParticleCount, &counts));
    PetscCall(PicurvAssertRealNear(2.0, counts[1][1][1], 1.0e-12, "Two particles in cell (0,0,0) should accumulate at shifted index (1,1,1)"));
    PetscCall(PicurvAssertRealNear(1.0, counts[1][1][2], 1.0e-12, "One particle in cell (1,0,0) should accumulate at shifted index (2,1,1)"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->ParticleCount, &counts));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests localized particle-status reset behavior for restart of the location workflow.
 */

static PetscErrorCode TestResetAllParticleStatusesLeavesLostParticlesUntouched(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt *status = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 3, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    status[0] = ACTIVE_AND_LOCATED;
    status[1] = LOST;
    status[2] = NEEDS_LOCATION;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    PetscCall(ResetAllParticleStatuses(user));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(PicurvAssertIntEqual(NEEDS_LOCATION, status[0], "ACTIVE_AND_LOCATED particles should be reset to NEEDS_LOCATION"));
    PetscCall(PicurvAssertIntEqual(LOST, status[1], "LOST particles should remain LOST"));
    PetscCall(PicurvAssertIntEqual(NEEDS_LOCATION, status[2], "NEEDS_LOCATION particles should remain unchanged"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests direct removal of particles that leave every rank bounding box.
 */

static PetscErrorCode TestCheckAndRemoveOutOfBoundsParticlesRemovesEscapedParticle(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal *positions = NULL;
    PetscInt removed_local = 0;
    PetscInt removed_global = 0;
    PetscInt nlocal = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    positions[0] = 0.5; positions[1] = 0.5; positions[2] = 0.5;
    positions[3] = 9.0; positions[4] = 9.0; positions[5] = 9.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(CheckAndRemoveOutOfBoundsParticles(user, &removed_local, &removed_global, simCtx->bboxlist));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(PicurvAssertIntEqual(1, removed_local, "Exactly one particle should be removed as out-of-bounds on a single rank"));
    PetscCall(PicurvAssertIntEqual(1, removed_global, "Global out-of-bounds removal count should match the local single-rank result"));
    PetscCall(PicurvAssertIntEqual(1, nlocal, "One in-bounds particle should remain after out-of-bounds removal"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests direct removal of particles already marked LOST by the location workflow.
 */

static PetscErrorCode TestCheckAndRemoveLostParticlesRemovesLostEntries(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt *status = NULL;
    PetscInt removed_local = 0;
    PetscInt removed_global = 0;
    PetscInt nlocal = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 3, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    status[0] = ACTIVE_AND_LOCATED;
    status[1] = LOST;
    status[2] = LOST;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    PetscCall(CheckAndRemoveLostParticles(user, &removed_local, &removed_global));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(PicurvAssertIntEqual(2, removed_local, "Two LOST particles should be removed locally"));
    PetscCall(PicurvAssertIntEqual(2, removed_global, "Global LOST-particle removal count should match the local single-rank result"));
    PetscCall(PicurvAssertIntEqual(1, nlocal, "One non-LOST particle should remain"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests Brownian displacement generation against a duplicated seeded RNG stream.
 */

static PetscErrorCode TestCalculateBrownianDisplacementDeterministicSeed(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    Cmpnts first;
    Cmpnts second;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContext(NULL, PETSC_FALSE, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    simCtx->dt = 0.25;
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->BrownianMotionRNG != NULL),
                               "runtime setup path should initialize the Brownian RNG"));

    PetscCall(PetscRandomSetSeed(simCtx->BrownianMotionRNG, 12345));
    PetscCall(PetscRandomSeed(simCtx->BrownianMotionRNG));

    PetscCall(CalculateBrownianDisplacement(user, 0.5, &first));
    PetscCall(PetscRandomSetSeed(simCtx->BrownianMotionRNG, 12345));
    PetscCall(PetscRandomSeed(simCtx->BrownianMotionRNG));
    PetscCall(CalculateBrownianDisplacement(user, 0.5, &second));

    PetscCall(PicurvAssertRealNear(first.x, second.x, 1.0e-12, "Resetting the Brownian RNG seed should reproduce the x displacement"));
    PetscCall(PicurvAssertRealNear(first.y, second.y, 1.0e-12, "Resetting the Brownian RNG seed should reproduce the y displacement"));
    PetscCall(PicurvAssertRealNear(first.z, second.z, 1.0e-12, "Resetting the Brownian RNG seed should reproduce the z displacement"));

    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests swarm-wide particle position updates using the same transport path as the runtime loop.
 */
static PetscErrorCode TestUpdateAllParticlePositionsMovesSwarmEntries(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal *positions = NULL;
    PetscReal *velocities = NULL;
    PetscReal *diffusivity = NULL;
    Cmpnts *diffusivity_gradient = NULL;
    PetscReal *psi = NULL;
    PetscReal *weights = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;
    PetscInt64 *pid = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));
    simCtx->dt = 0.25;

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void **)&velocities));
    PetscCall(DMSwarmGetField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmGetField(user->swarm, "DiffusivityGradient", NULL, NULL, (void **)&diffusivity_gradient));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmGetField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));

    positions[0] = 0.20; positions[1] = 0.30; positions[2] = 0.40;
    velocities[0] = 0.40; velocities[1] = -0.20; velocities[2] = 0.10;
    diffusivity[0] = 0.0;
    diffusivity_gradient[0].x = 0.10;
    diffusivity_gradient[0].y = 0.20;
    diffusivity_gradient[0].z = -0.10;
    psi[0] = 0.5;
    weights[0] = 0.5; weights[1] = 0.5; weights[2] = 0.5;
    cell_ids[0] = 0; cell_ids[1] = 0; cell_ids[2] = 0;
    status[0] = ACTIVE_AND_LOCATED;
    pid[0] = 7;

    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmRestoreField(user->swarm, "DiffusivityGradient", NULL, NULL, (void **)&diffusivity_gradient));
    PetscCall(DMSwarmRestoreField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void **)&velocities));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(UpdateAllParticlePositions(user));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(PicurvAssertRealNear(0.325, positions[0], 1.0e-12, "UpdateAllParticlePositions should advect x"));
    PetscCall(PicurvAssertRealNear(0.300, positions[1], 1.0e-12, "UpdateAllParticlePositions should advect y"));
    PetscCall(PicurvAssertRealNear(0.400, positions[2], 1.0e-12, "UpdateAllParticlePositions should advect z"));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the location orchestrator fast path when a particle already carries a valid prior cell.
 */
static PetscErrorCode TestLocateAllParticlesInGridPriorCellFastPath(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal *positions = NULL;
    PetscReal *weights = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;
    PetscInt64 *pid = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMalloc1(simCtx->size, &user->RankCellInfoMap));
    PetscCall(GetOwnedCellRange(&user->info, 0, &user->RankCellInfoMap[0].xs_cell, &user->RankCellInfoMap[0].xm_cell));
    PetscCall(GetOwnedCellRange(&user->info, 1, &user->RankCellInfoMap[0].ys_cell, &user->RankCellInfoMap[0].ym_cell));
    PetscCall(GetOwnedCellRange(&user->info, 2, &user->RankCellInfoMap[0].zs_cell, &user->RankCellInfoMap[0].zm_cell));
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    positions[0] = 0.375; positions[1] = 0.375; positions[2] = 0.375;
    weights[0] = 0.5; weights[1] = 0.5; weights[2] = 0.5;
    cell_ids[0] = 1; cell_ids[1] = 1; cell_ids[2] = 1;
    status[0] = NEEDS_LOCATION;
    pid[0] = 11;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(LocateAllParticlesInGrid(user, simCtx->bboxlist));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(PicurvAssertIntEqual(1, cell_ids[0], "prior-cell fast path should preserve the i cell id"));
    PetscCall(PicurvAssertIntEqual(1, cell_ids[1], "prior-cell fast path should preserve the j cell id"));
    PetscCall(PicurvAssertIntEqual(1, cell_ids[2], "prior-cell fast path should preserve the k cell id"));
    PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status[0], "prior-cell fast path should mark the particle ACTIVE_AND_LOCATED"));
    PetscCall(PicurvAssertIntEqual(1, simCtx->searchMetrics.searchAttempts, "prior-cell fast path should record one search attempt"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.searchPopulation, "prior-cell fast path should record one input particle"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.searchLocatedCount, "prior-cell fast path should count one located particle"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.searchLostCount, "prior-cell fast path should not lose the particle"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.reSearchCount, "prior-cell fast path should not re-search on later passes"));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->searchMetrics.traversalStepsSum > 0), "prior-cell fast path should accumulate traversal steps"));
    PetscCall(PicurvAssertIntEqual(1, simCtx->searchMetrics.maxParticlePassDepth, "prior-cell fast path should report one settlement pass"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the guess-then-verify orchestrator path for a local particle with an unknown prior cell.
 */
static PetscErrorCode TestLocateAllParticlesInGridGuessPathResolvesLocalParticle(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal *positions = NULL;
    PetscReal *weights = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;
    PetscInt64 *pid = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMalloc1(simCtx->size, &user->RankCellInfoMap));
    PetscCall(GetOwnedCellRange(&user->info, 0, &user->RankCellInfoMap[0].xs_cell, &user->RankCellInfoMap[0].xm_cell));
    PetscCall(GetOwnedCellRange(&user->info, 1, &user->RankCellInfoMap[0].ys_cell, &user->RankCellInfoMap[0].ym_cell));
    PetscCall(GetOwnedCellRange(&user->info, 2, &user->RankCellInfoMap[0].zs_cell, &user->RankCellInfoMap[0].zm_cell));
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    positions[0] = 0.625; positions[1] = 0.625; positions[2] = 0.625;
    weights[0] = 0.5; weights[1] = 0.5; weights[2] = 0.5;
    cell_ids[0] = -1; cell_ids[1] = -1; cell_ids[2] = -1;
    status[0] = NEEDS_LOCATION;
    pid[0] = 22;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField(user->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(LocateAllParticlesInGrid(user, simCtx->bboxlist));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(PicurvAssertIntEqual(2, cell_ids[0], "guess-path location should resolve the i cell id"));
    PetscCall(PicurvAssertIntEqual(2, cell_ids[1], "guess-path location should resolve the j cell id"));
    PetscCall(PicurvAssertIntEqual(2, cell_ids[2], "guess-path location should resolve the k cell id"));
    PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status[0], "guess-path location should mark the particle ACTIVE_AND_LOCATED"));
    PetscCall(PicurvAssertIntEqual(1, simCtx->searchMetrics.searchAttempts, "guess-path location should perform one robust search"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.searchPopulation, "guess-path location should count one input particle"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.searchLocatedCount, "guess-path location should count one located particle"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.searchLostCount, "guess-path location should not lose the particle"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.reSearchCount, "guess-path location should not count later-pass re-searches"));
    PetscCall(PicurvAssertIntEqual(1, simCtx->searchMetrics.bboxGuessFallbackCount, "guess-path location should record one bbox fallback"));
    PetscCall(PicurvAssertIntEqual(0, simCtx->searchMetrics.bboxGuessSuccessCount, "guess-path local resolution should not count as remote bbox success"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Verifies that later settlement passes increment re-search metrics.
 */
static PetscErrorCode TestLocateParticleOrFindMigrationTargetCountsReSearch(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Particle particle;
    ParticleLocationStatus status = NEEDS_LOCATION;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&particle, sizeof(particle)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMalloc1(simCtx->size, &user->RankCellInfoMap));
    PetscCall(GetOwnedCellRange(&user->info, 0, &user->RankCellInfoMap[0].xs_cell, &user->RankCellInfoMap[0].xm_cell));
    PetscCall(GetOwnedCellRange(&user->info, 1, &user->RankCellInfoMap[0].ys_cell, &user->RankCellInfoMap[0].ym_cell));
    PetscCall(GetOwnedCellRange(&user->info, 2, &user->RankCellInfoMap[0].zs_cell, &user->RankCellInfoMap[0].zm_cell));

    simCtx->searchMetrics.currentSettlementPass = 2;
    particle.PID = 33;
    particle.cell[0] = 1;
    particle.cell[1] = 1;
    particle.cell[2] = 1;
    particle.loc.x = 0.375;
    particle.loc.y = 0.375;
    particle.loc.z = 0.375;
    particle.weights.x = 0.5;
    particle.weights.y = 0.5;
    particle.weights.z = 0.5;

    PetscCall(LocateParticleOrFindMigrationTarget(user, &particle, &status));

    PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status, "direct re-search test should locate the particle"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.searchAttempts, "direct re-search test should record one robust walk"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.reSearchCount, "direct re-search test should increment re_search_count on later passes"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.maxTraversalFailCount, "direct re-search test should not hit MAX_TRAVERSAL"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests no-slip and free-slip wall helper kernels.
 */

static PetscErrorCode TestWallNoSlipAndFreeSlipHelpers(void)
{
    Cmpnts wall_velocity = {0.0, 0.0, 0.0};
    Cmpnts reference_velocity = {2.0, 4.0, 6.0};
    Cmpnts boundary_velocity = {0.0, 0.0, 0.0};
    Cmpnts free_slip_reference = {2.0, 3.0, 4.0};

    PetscFunctionBeginUser;
    noslip(NULL, 2.0, 1.0, wall_velocity, reference_velocity, &boundary_velocity, 1.0, 0.0, 0.0);
    PetscCall(PicurvAssertRealNear(1.0, boundary_velocity.x, 1.0e-12, "no-slip interpolated x"));
    PetscCall(PicurvAssertRealNear(2.0, boundary_velocity.y, 1.0e-12, "no-slip interpolated y"));
    PetscCall(PicurvAssertRealNear(3.0, boundary_velocity.z, 1.0e-12, "no-slip interpolated z"));

    freeslip(NULL, 2.0, 1.0, wall_velocity, free_slip_reference, &boundary_velocity, 1.0, 0.0, 0.0);
    PetscCall(PicurvAssertRealNear(1.0, boundary_velocity.x, 1.0e-12, "free-slip interpolated normal component"));
    PetscCall(PicurvAssertRealNear(3.0, boundary_velocity.y, 1.0e-12, "free-slip tangential y preserved"));
    PetscCall(PicurvAssertRealNear(4.0, boundary_velocity.z, 1.0e-12, "free-slip tangential z preserved"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests wall-model scalar helper kernels.
 */

static PetscErrorCode TestWallModelScalarHelpers(void)
{
    const PetscReal expected_smooth_e = PetscExpReal(0.41 * 5.5);
    PetscReal e_coeff = 0.0;
    PetscReal utau = 0.0;
    PetscReal residual = 0.0;

    PetscFunctionBeginUser;
    e_coeff = E_coeff(0.1, 0.0, 1.0e-3);
    PetscCall(PicurvAssertRealNear(expected_smooth_e, e_coeff, 1.0e-10, "smooth-wall E coefficient"));

    utau = find_utau_hydset(1.0e-3, 1.0, 1.0e-2, 0.1, 0.0);
    PetscCall(PicurvAssertBool((PetscBool)(utau > 0.0), "friction velocity should remain positive"));
    residual = f_hydset(1.0e-3, 1.0, 1.0e-2, utau, 0.0);
    PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(residual) < 1.0e-5), "Newton solve residual should be small"));

    PetscCall(PicurvAssertRealNear(0.0, nu_t(0.0), 1.0e-12, "eddy viscosity ratio at wall"));
    PetscCall(PicurvAssertBool((PetscBool)(integrate_1(1.0e-3, 1.0e-2, 0.1, 0) > 0.0), "integral helper should be positive"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests closed-form and iterative wall-model velocity helpers against inverse reconstructions.
 */

static PetscErrorCode TestWallModelVelocityHelpers(void)
{
    const PetscReal kinematic_viscosity = 1.0e-3;
    const PetscReal wall_distance = 2.0e-2;
    const PetscReal target_velocity = 1.0;
    const PetscReal roughness_length = 1.0e-4;
    PetscReal utau_loglaw = 0.0;
    PetscReal utau_werner = 0.0;
    PetscReal utau_cabot = 0.0;
    PetscReal wall_shear_velocity = 0.0;
    PetscReal wall_shear_normal = 0.0;

    PetscFunctionBeginUser;
    utau_loglaw = find_utau_loglaw(target_velocity, wall_distance, roughness_length);
    PetscCall(PicurvAssertRealNear(target_velocity, u_loglaw(wall_distance, utau_loglaw, roughness_length), 1.0e-12,
                                   "simple log-law inversion should reconstruct the target velocity"));

    utau_werner = find_utau_Werner(kinematic_viscosity, target_velocity, wall_distance, 0.1);
    PetscCall(PicurvAssertBool((PetscBool)(utau_werner > 0.0), "Werner-Wengle friction velocity should remain positive"));
    PetscCall(PicurvAssertRealNear(target_velocity, u_Werner(kinematic_viscosity, wall_distance, utau_werner), 1.0e-6,
                                   "Werner-Wengle inversion should reconstruct the target velocity"));

    find_utau_Cabot(kinematic_viscosity, target_velocity, wall_distance, 0.1, 0.0, 0.0,
                    &utau_cabot, &wall_shear_velocity, &wall_shear_normal);
    PetscCall(PicurvAssertBool((PetscBool)(utau_cabot > 0.0), "Cabot friction velocity should remain positive"));
    PetscCall(PicurvAssertRealNear(target_velocity, u_Cabot(kinematic_viscosity, wall_distance, utau_cabot, 0.0, wall_shear_velocity), 1.0e-6,
                                   "Cabot inversion should reconstruct the target velocity when pressure gradient is zero"));
    PetscCall(PicurvAssertRealNear(0.0, wall_shear_normal, 1.0e-10,
                                   "zero normal pressure gradient should keep Cabot normal wall shear at zero"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the vector wall-function wrappers on a tangential reference flow.
 */
static PetscErrorCode TestWallFunctionVectorWrappers(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts wall_velocity = {0.0, 0.0, 0.0};
    Cmpnts reference_velocity = {0.0, 1.0, 0.0};
    Cmpnts boundary_velocity = {0.0, 0.0, 0.0};
    PetscReal friction_velocity = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->ren = 1000.0;

    wall_function(user, 2.0e-2, 1.0e-2, wall_velocity, reference_velocity, &boundary_velocity, &friction_velocity, 1.0, 0.0, 0.0);
    PetscCall(PicurvAssertRealNear(0.0, boundary_velocity.x, 1.0e-12, "Werner wall function should preserve zero normal velocity"));
    PetscCall(PicurvAssertBool((PetscBool)(boundary_velocity.y > 0.0 && boundary_velocity.y < 1.0), "Werner wall function should damp tangential velocity"));
    PetscCall(PicurvAssertBool((PetscBool)(friction_velocity > 0.0), "Werner wall function should compute positive friction velocity"));

    wall_function_loglaw(user, 1.0e-4, 2.0e-2, 1.0e-2, wall_velocity, reference_velocity, &boundary_velocity, &friction_velocity, 1.0, 0.0, 0.0);
    PetscCall(PicurvAssertRealNear(0.0, boundary_velocity.x, 1.0e-12, "log-law wall function should preserve zero normal velocity"));
    PetscCall(PicurvAssertBool((PetscBool)(boundary_velocity.y > 0.0 && boundary_velocity.y <= 1.0), "log-law wall function should keep tangential velocity bounded"));
    PetscCall(PicurvAssertBool((PetscBool)(friction_velocity > 0.0), "log-law wall function should compute positive friction velocity"));

    wall_function_Cabot(user, 1.0e-4, 2.0e-2, 1.0e-2, wall_velocity, reference_velocity, &boundary_velocity, &friction_velocity,
                        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10);
    PetscCall(PicurvAssertRealNear(0.0, boundary_velocity.x, 1.0e-12, "Cabot wall function should preserve zero normal velocity"));
    PetscCall(PicurvAssertBool((PetscBool)(boundary_velocity.y > 0.0 && boundary_velocity.y <= 1.0), "Cabot wall function should keep tangential velocity bounded"));
    PetscCall(PicurvAssertBool((PetscBool)(friction_velocity > 0.0), "Cabot wall function should compute positive friction velocity"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests driven-flow validation when no driven handlers are present.
 */

static PetscErrorCode TestValidateDrivenFlowConfigurationNoDrivenHandlers(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(Validate_DrivenFlowConfiguration(user));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the constant Smagorinsky model helper path.
 */

static PetscErrorCode TestComputeSmagorinskyConstantConstantModel(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    simCtx->step = 2;
    simCtx->StartStep = 0;
    simCtx->les = CONSTANT_SMAGORINSKY;
    simCtx->Const_CS = 0.17;

    PetscCall(ComputeSmagorinskyConstant(user));
    PetscCall(PicurvAssertVecConstant(user->CS, 0.17, 1.0e-12, "constant Smagorinsky branch should fill CS"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the shared minimal fixture mirrors the production DA contract.
 */

static PetscErrorCode TestMinimalFixtureMirrorsProductionDMLayout(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    DM coord_dm = NULL;
    PetscInt mx = 0, my = 0, mz = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));

    PetscCall(DMDAGetInfo(user->da, NULL, &mx, &my, &mz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
    PetscCall(PicurvAssertIntEqual(user->IM + 1, mx, "minimal fixture should size da with IM+1 nodes"));
    PetscCall(PicurvAssertIntEqual(user->JM + 1, my, "minimal fixture should size da with JM+1 nodes"));
    PetscCall(PicurvAssertIntEqual(user->KM + 1, mz, "minimal fixture should size da with KM+1 nodes"));
    PetscCall(PicurvAssertIntEqual(mx, user->info.mx, "user->info should be sourced from the production da"));
    PetscCall(PicurvAssertIntEqual(my, user->info.my, "user->info my should match the da dimensions"));
    PetscCall(PicurvAssertIntEqual(mz, user->info.mz, "user->info mz should match the da dimensions"));

    PetscCall(DMGetCoordinateDM(user->da, &coord_dm));
    PetscCall(PicurvAssertBool((PetscBool)(coord_dm == user->fda),
                               "minimal fixture should derive fda from the coordinate-DM path"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the shared swarm fixture registers the production field set.
 */

static PetscErrorCode TestMinimalFixtureRegistersProductionSwarmFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt bs = 0;
    void *field_ptr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "position", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(3, bs, "solver swarm should register particle position"));
    PetscCall(PicurvAssertBool((PetscBool)(field_ptr != NULL), "position field should be retrievable"));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(3, bs, "solver swarm should register particle velocity"));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(3, bs, "solver swarm should register DMSwarm_CellID"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "weight", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(3, bs, "solver swarm should register particle weight"));
    PetscCall(DMSwarmRestoreField(user->swarm, "weight", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "Diffusivity", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(1, bs, "solver swarm should register particle diffusivity"));
    PetscCall(DMSwarmRestoreField(user->swarm, "Diffusivity", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "DiffusivityGradient", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(3, bs, "solver swarm should register particle diffusivity gradients"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DiffusivityGradient", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "Psi", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(1, bs, "solver swarm should register particle scalar Psi"));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", &bs, NULL, &field_ptr));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(1, bs, "solver swarm should register particle location status"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", &bs, NULL, &field_ptr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests solver history-vector shifting between time levels.
 */

static PetscErrorCode TestUpdateSolverHistoryVectorsShiftsStates(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));

    PetscCall(VecSet(user->Ucont, 11.0));
    PetscCall(VecSet(user->Ucont_o, 7.0));
    PetscCall(VecSet(user->Ucont_rm1, 3.0));
    PetscCall(VecSet(user->Ucat, 5.0));
    PetscCall(VecSet(user->Ucat_o, -1.0));
    PetscCall(VecSet(user->P, 9.0));
    PetscCall(VecSet(user->P_o, -2.0));

    PetscCall(UpdateSolverHistoryVectors(user));

    PetscCall(PicurvAssertVecConstant(user->Ucont_o, 11.0, 1.0e-12, "Ucont_o should receive current Ucont"));
    PetscCall(PicurvAssertVecConstant(user->Ucont_rm1, 7.0, 1.0e-12, "Ucont_rm1 should receive prior Ucont_o"));
    PetscCall(PicurvAssertVecConstant(user->Ucat_o, 5.0, 1.0e-12, "Ucat_o should receive current Ucat"));
    PetscCall(PicurvAssertVecConstant(user->P_o, 9.0, 1.0e-12, "P_o should receive current P"));
    PetscCall(PicurvAssertVecConstant(user->lUcont_o, 11.0, 1.0e-12, "lUcont_o ghost sync should match Ucont_o"));
    PetscCall(PicurvAssertVecConstant(user->lUcont_rm1, 7.0, 1.0e-12, "lUcont_rm1 ghost sync should match Ucont_rm1"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests owned-cell range accounting on a single MPI rank.
 */

static PetscErrorCode TestGetOwnedCellRangeSingleRankAccounting(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt xs = -1, ys = -1, zs = -1;
    PetscInt xm = -1, ym = -1, zm = -1;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));

    PetscCall(GetOwnedCellRange(&user->info, 0, &xs, &xm));
    PetscCall(GetOwnedCellRange(&user->info, 1, &ys, &ym));
    PetscCall(GetOwnedCellRange(&user->info, 2, &zs, &zm));

    PetscCall(PicurvAssertIntEqual(0, xs, "single-rank x cell-start index"));
    PetscCall(PicurvAssertIntEqual(0, ys, "single-rank y cell-start index"));
    PetscCall(PicurvAssertIntEqual(0, zs, "single-rank z cell-start index"));
    PetscCall(PicurvAssertIntEqual(user->info.mx - 2, xm, "single-rank x owned cell count"));
    PetscCall(PicurvAssertIntEqual(user->info.my - 2, ym, "single-rank y owned cell count"));
    PetscCall(PicurvAssertIntEqual(user->info.mz - 2, zm, "single-rank z owned cell count"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests neighbor-rank discovery on a single MPI rank.
 */

static PetscErrorCode TestComputeAndStoreNeighborRanksSingleRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(ComputeAndStoreNeighborRanks(user));

    PetscCall(PicurvAssertIntEqual(MPI_PROC_NULL, user->neighbors.rank_xm, "single-rank xm neighbor should be null"));
    PetscCall(PicurvAssertIntEqual(MPI_PROC_NULL, user->neighbors.rank_xp, "single-rank xp neighbor should be null"));
    PetscCall(PicurvAssertIntEqual(MPI_PROC_NULL, user->neighbors.rank_ym, "single-rank ym neighbor should be null"));
    PetscCall(PicurvAssertIntEqual(MPI_PROC_NULL, user->neighbors.rank_yp, "single-rank yp neighbor should be null"));
    PetscCall(PicurvAssertIntEqual(MPI_PROC_NULL, user->neighbors.rank_zm, "single-rank zm neighbor should be null"));
    PetscCall(PicurvAssertIntEqual(MPI_PROC_NULL, user->neighbors.rank_zp, "single-rank zp neighbor should be null"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests parsing of positive runtime walltime metadata values.
 */

static PetscErrorCode TestRuntimeWalltimeGuardParsesPositiveSeconds(void)
{
    PetscReal seconds = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool(RuntimeWalltimeGuardParsePositiveSeconds("300", &seconds), "plain positive seconds should parse"));
    PetscCall(PicurvAssertRealNear(300.0, seconds, 1.0e-12, "parsed walltime seconds"));
    PetscCall(PicurvAssertBool(RuntimeWalltimeGuardParsePositiveSeconds(" 42.5 ", &seconds), "whitespace-wrapped decimal seconds should parse"));
    PetscCall(PicurvAssertRealNear(42.5, seconds, 1.0e-12, "parsed decimal walltime seconds"));
    PetscCall(PicurvAssertBool((PetscBool)!RuntimeWalltimeGuardParsePositiveSeconds("nope", &seconds), "non-numeric metadata should fail parsing"));
    PetscCall(PicurvAssertBool((PetscBool)!RuntimeWalltimeGuardParsePositiveSeconds("-10", &seconds), "negative metadata should fail parsing"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests walltime-guard estimator helper calculations.
 */

static PetscErrorCode TestRuntimeWalltimeGuardEstimatorHelpers(void)
{
    PetscReal ewma_fast = 0.0;
    PetscReal ewma_slow = 0.0;
    PetscReal conservative_fast = 0.0;
    PetscReal conservative_slow = 0.0;
    PetscReal required_headroom = 0.0;

    PetscFunctionBeginUser;
    ewma_fast = RuntimeWalltimeGuardUpdateEWMA(PETSC_TRUE, 4.0, 6.0, 0.5);
    ewma_slow = RuntimeWalltimeGuardUpdateEWMA(PETSC_TRUE, ewma_fast, 12.0, 0.5);
    conservative_fast = RuntimeWalltimeGuardConservativeEstimate(5.0, ewma_fast, 6.0);
    conservative_slow = RuntimeWalltimeGuardConservativeEstimate(5.0, ewma_slow, 12.0);
    required_headroom = RuntimeWalltimeGuardRequiredHeadroom(8.0, 2.0, conservative_slow);

    PetscCall(PicurvAssertRealNear(5.0, ewma_fast, 1.0e-12, "EWMA after moderate step"));
    PetscCall(PicurvAssertRealNear(8.5, ewma_slow, 1.0e-12, "EWMA after newer slow step"));
    PetscCall(PicurvAssertRealNear(6.0, conservative_fast, 1.0e-12, "conservative estimate tracks latest moderate step"));
    PetscCall(PicurvAssertRealNear(12.0, conservative_slow, 1.0e-12, "conservative estimate tracks newest slow step"));
    PetscCall(PicurvAssertRealNear(24.0, required_headroom, 1.0e-12, "required headroom scales with conservative estimate"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests runtime walltime-guard shutdown trigger decisions.
 */

static PetscErrorCode TestRuntimeWalltimeGuardTriggerDecision(void)
{
    PetscBool should_trigger = PETSC_FALSE;
    PetscReal required_headroom = 0.0;

    PetscFunctionBeginUser;
    should_trigger = RuntimeWalltimeGuardShouldTrigger(9, 10, 15.0, 5.0, 2.0, 6.0, 6.0, 6.0, &required_headroom);
    PetscCall(PicurvAssertBool((PetscBool)!should_trigger, "guard should not trigger before warmup completes"));

    should_trigger = RuntimeWalltimeGuardShouldTrigger(10, 10, 40.0, 5.0, 2.0, 10.0, 12.0, 14.0, &required_headroom);
    PetscCall(PicurvAssertBool((PetscBool)!should_trigger, "guard should not trigger when remaining walltime exceeds required headroom"));
    PetscCall(PicurvAssertRealNear(28.0, required_headroom, 1.0e-12, "required headroom after warmup"));

    should_trigger = RuntimeWalltimeGuardShouldTrigger(10, 10, 28.0, 5.0, 2.0, 10.0, 12.0, 14.0, &required_headroom);
    PetscCall(PicurvAssertBool(should_trigger, "guard should trigger when remaining walltime reaches required headroom"));
    PetscCall(PicurvAssertRealNear(28.0, required_headroom, 1.0e-12, "required headroom remains unchanged at trigger threshold"));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-runtime PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"distribute-particles-remainder-handling", TestDistributeParticlesRemainderHandling},
        {"is-particle-inside-bbox-basic-cases", TestIsParticleInsideBoundingBoxBasicCases},
        {"update-particle-weights-computes-expected-ratios", TestUpdateParticleWeightsComputesExpectedRatios},
        {"update-particle-position-without-brownian-contribution", TestUpdateParticlePositionWithoutBrownianContribution},
        {"update-particle-position-diffusivity-gradient-only", TestUpdateParticlePositionDiffusivityGradientOnly},
        {"update-particle-field-iem-relaxation", TestUpdateParticleFieldIEMRelaxation},
        {"set-initial-interior-field-ignores-non-ucont-request", TestSetInitialInteriorFieldIgnoresNonUcontRequest},
        {"set-initial-interior-field-constant-profile-on-z-inlet", TestSetInitialInteriorFieldConstantProfileOnZInlet},
        {"interpolate-all-fields-to-swarm-constant-fields", TestInterpolateAllFieldsToSwarmConstantFields},
        {"interpolate-all-fields-to-swarm-corner-averaged-constant-fields", TestInterpolateAllFieldsToSwarmCornerAveragedConstantFields},
        {"scatter-all-particle-fields-to-euler-fields-averages-psi", TestScatterAllParticleFieldsToEulerFieldsAveragesPsi},
        {"calculate-particle-count-per-cell-counts-global-cell-ids", TestCalculateParticleCountPerCellCountsGlobalCellIDs},
        {"reset-all-particle-statuses-leaves-lost-particles-untouched", TestResetAllParticleStatusesLeavesLostParticlesUntouched},
        {"check-and-remove-out-of-bounds-particles-removes-escaped-particle", TestCheckAndRemoveOutOfBoundsParticlesRemovesEscapedParticle},
        {"check-and-remove-lost-particles-removes-lost-entries", TestCheckAndRemoveLostParticlesRemovesLostEntries},
        {"calculate-brownian-displacement-deterministic-seed", TestCalculateBrownianDisplacementDeterministicSeed},
        {"update-all-particle-positions-moves-swarm-entries", TestUpdateAllParticlePositionsMovesSwarmEntries},
        {"locate-all-particles-in-grid-prior-cell-fast-path", TestLocateAllParticlesInGridPriorCellFastPath},
        {"locate-all-particles-in-grid-guess-path-resolves-local-particle", TestLocateAllParticlesInGridGuessPathResolvesLocalParticle},
        {"locate-particle-or-find-migration-target-counts-research", TestLocateParticleOrFindMigrationTargetCountsReSearch},
        {"wall-noslip-and-freeslip-helpers", TestWallNoSlipAndFreeSlipHelpers},
        {"wall-model-scalar-helpers", TestWallModelScalarHelpers},
        {"wall-model-velocity-helpers", TestWallModelVelocityHelpers},
        {"wall-function-vector-wrappers", TestWallFunctionVectorWrappers},
        {"validate-driven-flow-configuration-no-driven-handlers", TestValidateDrivenFlowConfigurationNoDrivenHandlers},
        {"compute-smagorinsky-constant-constant-model", TestComputeSmagorinskyConstantConstantModel},
        {"minimal-fixture-mirrors-production-dm-layout", TestMinimalFixtureMirrorsProductionDMLayout},
        {"minimal-fixture-registers-production-swarm-fields", TestMinimalFixtureRegistersProductionSwarmFields},
        {"update-solver-history-vectors-shifts-states", TestUpdateSolverHistoryVectorsShiftsStates},
        {"get-owned-cell-range-single-rank-accounting", TestGetOwnedCellRangeSingleRankAccounting},
        {"compute-and-store-neighbor-ranks-single-rank", TestComputeAndStoreNeighborRanksSingleRank},
        {"runtime-walltime-guard-parses-positive-seconds", TestRuntimeWalltimeGuardParsesPositiveSeconds},
        {"runtime-walltime-guard-estimator-helpers", TestRuntimeWalltimeGuardEstimatorHelpers},
        {"runtime-walltime-guard-trigger-decision", TestRuntimeWalltimeGuardTriggerDecision},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv runtime-kernel tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-runtime", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
