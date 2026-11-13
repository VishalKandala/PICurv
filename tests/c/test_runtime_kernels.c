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

static PetscErrorCode TestIsParticleInBoxBasicCases(void)
{
    BoundingBox bbox;
    Cmpnts      pos;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&bbox, sizeof(bbox)));

    bbox.min_coords.x = 0.0;
    bbox.min_coords.y = 0.0;
    bbox.min_coords.z = 0.0;
    bbox.max_coords.x = 1.0;
    bbox.max_coords.y = 2.0;
    bbox.max_coords.z = 3.0;

    pos.x = 0.25;
    pos.y = 1.0;
    pos.z = 2.5;
    PetscCall(PicurvAssertBool(IsParticleInBox(&bbox, &pos), "particle should be inside bounding box"));

    pos.x = 1.0;
    PetscCall(PicurvAssertBool(IsParticleInBox(&bbox, &pos), "particle on the box face should count as inside"));

    pos.x = 1.5;
    PetscCall(PicurvAssertBool((PetscBool)!IsParticleInBox(&bbox, &pos), "particle should be outside bounding box"));
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
    PetscCall(UpdateParticleField(PARTICLE_FIELD_ID_PSI, dt, &psi, diffusivity, mean_val, PetscPowReal(cell_vol, 0.6666667), c_model));
    PetscCall(PicurvAssertRealNear(
        mean_val + (1.0 - mean_val) * PetscExpReal(-(c_model * diffusivity / PetscPowReal(cell_vol, 0.6666667)) * dt),
        psi,
        1.0e-12,
        "IEM update should match analytical relaxation"));

    PetscCall(UpdateParticleField(PARTICLE_FIELD_ID_VELOCITY, dt, &unchanged, diffusivity, mean_val, PetscPowReal(cell_vol, 0.6666667), c_model));
    PetscCall(PicurvAssertRealNear(7.0, unchanged, 1.0e-12, "unknown field should remain unchanged"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests spatially varying Cartesian-field conversion to contravariant fluxes.
 */
static PetscErrorCode TestCart2ContraConvertsCartesianField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucat = NULL;
    Cmpnts ***ucont = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = 0; k < 6; k++)
        for (PetscInt j = 0; j < 6; j++)
            for (PetscInt i = 0; i < 6; i++)
                ucat[k][j][i] = (Cmpnts){(PetscReal)i, 2.0, 3.0};
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));
    PetscCall(Cart2Contra(user));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear(1.5, ucont[1][1][1].x, 1.0e-10,
                                   "Xi flux uses face-interpolated Ucat"));
    PetscCall(PicurvAssertRealNear(2.0, ucont[1][1][1].y, 1.0e-10,
                                   "Eta flux uses Cartesian y velocity"));
    PetscCall(PicurvAssertRealNear(3.0, ucont[1][1][1].z, 1.0e-10,
                                   "Zeta flux uses Cartesian z velocity"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Verifies periodic Ucat endpoints are repaired before Cart2Contra terminal faces. */
static PetscErrorCode TestCart2ContraUsesFinalizedPeriodicUcat(void)
{
    SimCtx *simCtx=NULL; UserCtx *user=NULL; Cmpnts ***ucat=NULL,***ucont=NULL; PetscInt mx;
    const FieldId cell_fields[] = {FIELD_ID_UCAT}, staggered_fields[] = {FIELD_ID_UCONT};
    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx,&user,4,4,4,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE));
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type=PERIODIC;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type=PERIODIC;
    mx=user->info.mx;
    PetscCall(VecSet(user->Ucat,0.0));
    PetscCall(DMDAVecGetArray(user->fda,user->Ucat,&ucat));
    for(PetscInt k=1;k<user->info.mz-1;k++) for(PetscInt j=1;j<user->info.my-1;j++)
        for(PetscInt i=1;i<mx-1;i++) ucat[k][j][i]=(Cmpnts){10.0+(PetscReal)i,2.0,3.0};
    ucat[2][2][0]=(Cmpnts){-1000.0,-1000.0,-1000.0};
    ucat[2][2][mx-1]=(Cmpnts){-2000.0,-2000.0,-2000.0};
    PetscCall(DMDAVecRestoreArray(user->fda,user->Ucat,&ucat));
    PetscCall(SynchronizePeriodicCellFields(user,1,cell_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));
    PetscCall(Cart2Contra(user));
    PetscCall(SynchronizePeriodicStaggeredFields(user,1,staggered_fields));
    PetscCall(DMDAVecGetArrayRead(user->fda,user->Ucont,&ucont));
    PetscCall(PicurvAssertRealNear(0.5*((10.0+mx-2)+(10.0+1)),ucont[2][2][mx-2].x,1e-12,"terminal Xi face uses corrected periodic Ucat endpoint"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda,user->Ucont,&ucont));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx,&user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests loading a staged Ucat file IC and converting it to Ucont.
 */
static PetscErrorCode TestPopulateInitialUcontLoadsStagedUcat(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucont = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(VecSet(user->Ucat, 2.5));
    PetscCall(PetscStrncpy(simCtx->_io_context_buffer, tmpdir, sizeof(simCtx->_io_context_buffer)));
    simCtx->current_io_directory = simCtx->_io_context_buffer;
    PetscCall(WriteFieldData(user, "ufield00000_0", user->Ucat, "dat"));
    simCtx->current_io_directory = NULL;

    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(VecZeroEntries(user->Ucont));
    simCtx->initialConditionStaged = PETSC_TRUE;
    simCtx->initialConditionField = IC_FIELD_UCAT;
    PetscCall(PetscStrncpy(simCtx->initialConditionDirectory, tmpdir, sizeof(simCtx->initialConditionDirectory)));
    PetscCall(PopulateInitialUcont(user));

    PetscCall(PicurvAssertVecConstant(user->Ucat, 2.5, 1.0e-12, "file IC should restore Ucat"));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear(2.5, ucont[1][1][1].x, 1.0e-12, "file Ucat IC should populate interior Xi flux"));
    PetscCall(PicurvAssertRealNear(2.5, ucont[1][1][1].y, 1.0e-12, "file Ucat IC should populate interior Eta flux"));
    PetscCall(PicurvAssertRealNear(2.5, ucont[1][1][1].z, 1.0e-12, "file Ucat IC should populate interior Zeta flux"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Tests LES output staging from local vectors on a multiply periodic DMDA. */
static PetscErrorCode TestWriteLESFieldsCopiesOwnedPeriodicValues(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char euler_dir[PETSC_MAX_PATH_LEN];
    char block_dir[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(
        &simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    PetscCall(DMCreateLocalVector(user->da, &user->lCs));
    PetscCall(DMCreateGlobalVector(user->da, &user->Nu_t));
    PetscCall(DMCreateLocalVector(user->da, &user->lNu_t));
    PetscCall(VecSet(user->CS, 0.0));
    PetscCall(VecSet(user->lCs, 2.5));
    PetscCall(VecSet(user->Nu_t, 0.0));
    PetscCall(VecSet(user->lNu_t, 7.5));
    simCtx->les = DYNAMIC_SMAGORINSKY;

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(euler_dir, sizeof(euler_dir), "%s/eulerian", tmpdir));
    PetscCall(PetscSNPrintf(block_dir, sizeof(block_dir), "%s/block_0000", euler_dir));
    PetscCall(PicurvEnsureDir(euler_dir));
    PetscCall(PicurvEnsureDir(block_dir));
    PetscCall(WriteSimulationFields(user, tmpdir));

    PetscCall(PicurvAssertVecConstant(user->CS, 2.5, 1.0e-12,
                                      "checkpoint output should copy owned periodic lCs values"));
    PetscCall(PicurvAssertVecConstant(user->Nu_t, 7.5, 1.0e-12,
                                      "checkpoint output should copy owned periodic lNu_t values"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Stages a real Cartesian Fourier mode through the production file-IC path.
 *
 * The returned maximum is the native face-interpolation/backward-difference
 * divergence.  On the uniform periodic fixture its symbol is
 * (sin(theta_x), sin(theta_y), sin(theta_z)).
 */
static PetscErrorCode StagedFourierModeMaxDivergence(UserCtx *user,
                                                     PetscInt mode_x,
                                                     PetscInt mode_y,
                                                     PetscInt mode_z,
                                                     Cmpnts amplitude,
                                                     const char *tmpdir,
                                                     PetscReal *max_divergence)
{
    SimCtx *simCtx = user->simCtx;
    Cmpnts ***ucat = NULL, ***ucont = NULL;
    PetscInt n_x = user->info.mx - 2, n_y = user->info.my - 2, n_z = user->info.mz - 2;
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};

    PetscFunctionBeginUser;
    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = 0; k < user->info.mz; ++k) {
        for (PetscInt j = 0; j < user->info.my; ++j) {
            for (PetscInt i = 0; i < user->info.mx; ++i) {
                PetscReal phase = 2.0 * PETSC_PI * (
                    (PetscReal)(mode_x * (i - 1)) / (PetscReal)n_x +
                    (PetscReal)(mode_y * (j - 1)) / (PetscReal)n_y +
                    (PetscReal)(mode_z * (k - 1)) / (PetscReal)n_z);
                PetscReal wave = PetscCosReal(phase);
                ucat[k][j][i] = (Cmpnts){amplitude.x * wave, amplitude.y * wave, amplitude.z * wave};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));

    PetscCall(PetscStrncpy(simCtx->_io_context_buffer, tmpdir, sizeof(simCtx->_io_context_buffer)));
    simCtx->current_io_directory = simCtx->_io_context_buffer;
    PetscCall(WriteFieldData(user, "ufield00000_0", user->Ucat, "dat"));
    simCtx->current_io_directory = NULL;
    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(VecZeroEntries(user->Ucont));
    PetscCall(PopulateInitialUcont(user));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCONT));

    *max_divergence = 0.0;
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont));
    for (PetscInt k = 1; k < user->info.mz - 1; ++k) {
        for (PetscInt j = 1; j < user->info.my - 1; ++j) {
            for (PetscInt i = 1; i < user->info.mx - 1; ++i) {
                PetscReal divergence =
                    ucont[k][j][i].x - ucont[k][j][i - 1].x +
                    ucont[k][j][i].y - ucont[k][j - 1][i].y +
                    ucont[k][j][i].z - ucont[k - 1][j][i].z;
                *max_divergence = PetscMax(*max_divergence, PetscAbsReal(divergence));
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont));
    PetscFunctionReturn(0);
}

/** @brief Cross-checks selected exact Fourier modes against PICurv's native discrete symbol. */
static PetscErrorCode TestStagedUcatFourierModesUsePicurvDiscreteSymbol(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    PetscReal max_divergence = 0.0, sx, sy, sz;
    PetscInt n;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(
        &simCtx, &user, 9, 9, 9, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    for (PetscInt face = 0; face < NUM_FACES; ++face)
        user->boundary_faces[face].mathematical_type = PERIODIC;
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    simCtx->initialConditionStaged = PETSC_TRUE;
    simCtx->initialConditionField = IC_FIELD_UCAT;
    PetscCall(PetscStrncpy(simCtx->initialConditionDirectory, tmpdir, sizeof(simCtx->initialConditionDirectory)));

    PetscCall(StagedFourierModeMaxDivergence(
        user, 1, 0, 0, (Cmpnts){0.0, 1.0, 0.0}, tmpdir, &max_divergence));
    PetscCall(PicurvAssertRealNear(0.0, max_divergence, 1.0e-12,
                                   "axis-aligned transverse staged mode"));

    n = user->info.mx - 2;
    sx = PetscSinReal(2.0 * PETSC_PI / (PetscReal)n);
    sy = PetscSinReal(4.0 * PETSC_PI / (PetscReal)n);
    sz = PetscSinReal(6.0 * PETSC_PI / (PetscReal)n);
    PetscCall(StagedFourierModeMaxDivergence(
        user, 1, 2, 0, (Cmpnts){sy, -sx, 0.0}, tmpdir, &max_divergence));
    PetscCall(PicurvAssertRealNear(0.0, max_divergence, 1.0e-12,
                                   "oblique first discrete transverse polarization"));
    PetscCall(StagedFourierModeMaxDivergence(
        user, 1, 2, 3,
        (Cmpnts){sz * sx, sz * sy, -(sx * sx + sy * sy)},
        tmpdir, &max_divergence));
    PetscCall(PicurvAssertRealNear(0.0, max_divergence, 1.0e-12,
                                   "oblique second discrete transverse polarization"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests loading a staged Ucont file IC without Cartesian conversion.
 */
static PetscErrorCode TestPopulateInitialUcontLoadsStagedUcont(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(VecSet(user->Ucont, 4.5));
    PetscCall(PetscStrncpy(simCtx->_io_context_buffer, tmpdir, sizeof(simCtx->_io_context_buffer)));
    simCtx->current_io_directory = simCtx->_io_context_buffer;
    PetscCall(WriteFieldData(user, "vfield00000_0", user->Ucont, "dat"));
    simCtx->current_io_directory = NULL;

    PetscCall(VecZeroEntries(user->Ucont));
    simCtx->initialConditionStaged = PETSC_TRUE;
    simCtx->initialConditionField = IC_FIELD_UCONT;
    PetscCall(PetscStrncpy(simCtx->initialConditionDirectory, tmpdir, sizeof(simCtx->initialConditionDirectory)));
    PetscCall(PopulateInitialUcont(user));

    PetscCall(PicurvAssertVecConstant(user->Ucont, 4.5, 1.0e-12, "file Ucont IC should restore Ucont directly"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
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
 * @brief FaceFlux: constant face fluxes give that velocity; other fields still interpolate trilinearly.
 */
static PetscErrorCode TestInterpolateAllFieldsToSwarmFaceFluxConstantFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal *velocity = NULL;
    PetscReal *diffusivity = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->interpolationMethod = INTERP_FACE_FLUX;
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));
    PetscCall(VecSet(user->Ucont, 2.0)); /* Unit face-area vectors: the fluxes are the velocity */
    PetscCall(VecSet(user->Ucat, -7.0));  /* Not read for velocity */
    PetscCall(VecSet(user->Diffusivity, 0.25));
    PetscCall(SyncRuntimeFieldGhosts(user));
    PetscCall(SeedSingleParticle(user, 1, 1, 1, 0.5, 0.5, 0.5, 0.2, 0.7, 0.4, ACTIVE_AND_LOCATED));

    PetscCall(InterpolateAllFieldsToSwarm(user));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));
    PetscCall(DMSwarmGetField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    for (PetscInt c = 0; c < 3; c++)
        PetscCall(PicurvAssertRealNear(2.0, velocity[c], 1.0e-12, "FaceFlux: velocity should come from the constant face fluxes"));
    PetscCall(PicurvAssertRealNear(0.25, diffusivity[0], 1.0e-12, "FaceFlux: diffusivity should still interpolate trilinearly"));
    PetscCall(DMSwarmRestoreField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief FaceFlux: divergence-free fluxes give a divergence-free velocity, zero at a wall face,
 *        and the Contra2Cart value at the cell centre.
 *
 * The fluxes come from a stream function psi(i, j) = i^2 j + 0.3 j^2 on cell edges, so every cell's
 * discrete flux divergence is zero and the j = 0 wall face carries no flux. The velocity is linear
 * in each weight, so one-sided differences of the weights give its divergence exactly.
 */
static PetscErrorCode TestFaceFluxVelocityIsDivergenceFree(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucont = NULL, ***ucat = NULL;
    PetscReal *pos = NULL, *wgt = NULL, *vel = NULL;
    PetscInt *cell = NULL, *status = NULL;
    const PetscReal d = 0.25, w0[3] = {0.3, 0.6, 0.45};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    simCtx->interpolationMethod = INTERP_FACE_FLUX;
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++)
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++)
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) {
                const PetscReal psi    = (PetscReal)(i * i * j) + 0.3 * j * j;
                const PetscReal psi_jm = (PetscReal)(i * i * (j - 1)) + 0.3 * (j - 1) * (j - 1);
                const PetscReal psi_im = (PetscReal)((i - 1) * (i - 1) * j) + 0.3 * j * j;
                ucont[k][j][i].x = psi - psi_jm;
                ucont[k][j][i].y = -(psi - psi_im);
                ucont[k][j][i].z = 0.0;
            }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(SyncRuntimeFieldGhosts(user));
    PetscCall(Contra2Cart(user));

    /* Particles in cell (2, 0, 2): at w0, w0 + d e_x, w0 + d e_y, w0 + d e_z, the centre, and on the j = 0 wall. */
    PetscCall(PicurvCreateSwarmPair(user, 6, "ske"));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmGetField(user->swarm, "weight", NULL, NULL, (void **)&wgt));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    for (PetscInt p = 0; p < 6; p++) {
        for (PetscInt c = 0; c < 3; c++) {
            wgt[3*p + c] = (p == 4) ? 0.5 : w0[c] + ((p == c + 1) ? d : 0.0);
            pos[3*p + c] = 0.5;
        }
        if (p == 5) wgt[3*p + 1] = 0.0;
        cell[3*p + 0] = 2; cell[3*p + 1] = 0; cell[3*p + 2] = 2;
        status[p] = ACTIVE_AND_LOCATED;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmRestoreField(user->swarm, "weight", NULL, NULL, (void **)&wgt));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&pos));

    PetscCall(InterpolateAllFieldsToSwarm(user));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void **)&vel));
    {
        const PetscReal div = (vel[3*1 + 0] - vel[0]) / d + (vel[3*2 + 1] - vel[1]) / d + (vel[3*3 + 2] - vel[2]) / d;
        const PetscReal scale = PetscAbsReal(vel[3*1 + 0] - vel[0]) / d;
        PetscCall(PicurvAssertBool((PetscBool)(scale > 1.0), "the test field should vary across the cell"));
        PetscCall(PicurvAssertRealNear(0.0, div, 1.0e-10 * scale, "FaceFlux velocity should be divergence-free"));
    }
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    for (PetscInt c = 0; c < 3; c++)
        PetscCall(PicurvAssertRealNear((&ucat[3][1][3].x)[c], vel[3*4 + c], 1.0e-12, "the cell-centre value should be Contra2Cart's"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(PicurvAssertRealNear(0.0, vel[3*5 + 1], 1.0e-12, "the wall-normal velocity should vanish on the wall face"));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void **)&vel));
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
 * @brief Tests removal of LOST and EXITED particles, their separate tallies, and that survivors keep the one-at-a-time removal order.
 */

static PetscErrorCode TestCheckAndRemoveLostParticlesRemovesLostEntries(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt *status = NULL;
    PetscInt removed_local = 0;
    PetscInt removed_global = 0;
    PetscInt nlocal = 0;
    PetscReal *psi = NULL;
    PetscReal removed_psi = 0.0;
    PetscInt exited_global = 0;
    PetscReal exited_psi = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 5, "ske"));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    psi[0] = 0.5;
    psi[1] = 2.0;
    psi[2] = 3.0;
    psi[3] = 7.0;
    psi[4] = 11.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    status[0] = ACTIVE_AND_LOCATED;
    status[1] = LOST;
    status[2] = EXITED;
    status[3] = ACTIVE_AND_LOCATED;
    status[4] = ACTIVE_AND_LOCATED;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    PetscCall(CheckAndRemoveLostParticles(user, &removed_local, &removed_global, &removed_psi, &exited_global, &exited_psi));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(PicurvAssertRealNear(2.0, removed_psi, 0.0, "the LOST particle carried Psi 2"));
    PetscCall(PicurvAssertRealNear(3.0, exited_psi, 0.0, "the EXITED particle carried Psi 3"));
    PetscCall(PicurvAssertIntEqual(1, removed_local, "One LOST particle should be removed locally"));
    PetscCall(PicurvAssertIntEqual(1, removed_global, "Global LOST-particle removal count should match the local single-rank result"));
    PetscCall(PicurvAssertIntEqual(1, exited_global, "One EXITED particle should be removed"));
    PetscCall(PicurvAssertIntEqual(3, nlocal, "Three located particles should remain"));

    /* Backward one-at-a-time removal moves the last particle into slot 2, then the new
     * last into slot 1: survivors end as Psi 0.5, 7, 11 (a forward scan would give 0.5, 11, 7). */
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(PicurvAssertRealNear(0.5, psi[0], 0.0, "survivor order must match one-at-a-time removal (slot 0)"));
    PetscCall(PicurvAssertRealNear(7.0, psi[1], 0.0, "survivor order must match one-at-a-time removal (slot 1)"));
    PetscCall(PicurvAssertRealNear(11.0, psi[2], 0.0, "survivor order must match one-at-a-time removal (slot 2)"));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that sorting by host cell orders points k-slowest, keeps ties and every field together, and puts unlocated points last.
 */
static PetscErrorCode TestSortParticlesByCellKeepsEachPointWhole(void)
{
    SimCtx    *simCtx = NULL;
    UserCtx   *user = NULL;
    PetscInt  *cell = NULL;
    PetscInt64 *pid = NULL;
    PetscReal *psi = NULL, *pos = NULL;
    /* Cells (i,j,k) per point; point 4 has no cell. Expected order by (k,j,i), stable:
     * p3 (0,0,0), p1 (1,0,0), p5 (1,0,0), p0 (2,1,0), p2 (0,0,1), p4 unlocated. */
    const PetscInt  cells_in[6][3] = {{2, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 0, 0}, {-1, -1, -1}, {1, 0, 0}};
    const PetscInt  expected[6] = {3, 1, 5, 0, 2, 4};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 6, "ske"));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&pos));
    for (PetscInt p = 0; p < 6; ++p) {
        for (int d = 0; d < 3; ++d) { cell[3 * p + d] = cells_in[p][d]; pos[3 * p + d] = 10.0 * p + d; }
        pid[p] = 100 + p;
        psi[p] = 0.5 + p;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));

    PetscCall(SortParticlesByCell(user->swarm));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&pos));
    for (PetscInt q = 0; q < 6; ++q) {
        const PetscInt p = expected[q];
        PetscCall(PicurvAssertIntEqual(100 + p, (PetscInt)pid[q], "points must come out in (k,j,i) cell order, ties and unlocated points stable"));
        PetscCall(PicurvAssertRealNear(0.5 + p, psi[q], 0.0, "Psi must move with its point"));
        PetscCall(PicurvAssertRealNear(10.0 * p + 2, pos[3 * q + 2], 0.0, "position must move with its point"));
        PetscCall(PicurvAssertIntEqual(cells_in[p][0], cell[3 * q], "CellID must move with its point"));
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that a Brownian increment is a function of (seed, PID, step) alone.
 *
 * The same key reproduces the increment whatever was drawn before it, which is what
 * makes a restart and a different rank layout draw what the uninterrupted run drew;
 * changing any one key changes it.
 */

static PetscErrorCode TestBrownianDisplacementIsKeyedOnSeedParticleAndStep(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts base, again, other_pid, other_step, other_seed, none;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->dt = 0.25;
    simCtx->particleRandomSeed = 777;
    simCtx->step = 12;
    PetscCall(CalculateBrownianDisplacement(user, 41, 0.5, &base));
    PetscCall(CalculateBrownianDisplacement(user, 42, 0.5, &other_pid));
    PetscCall(CalculateBrownianDisplacement(user, 41, 0.5, &again));
    simCtx->step = 13;
    PetscCall(CalculateBrownianDisplacement(user, 41, 0.5, &other_step));
    simCtx->step = 12;
    simCtx->particleRandomSeed = 778;
    PetscCall(CalculateBrownianDisplacement(user, 41, 0.5, &other_seed));
    PetscCall(CalculateBrownianDisplacement(user, 41, 0.0, &none));

    PetscCall(PicurvAssertBool((PetscBool)(base.x == again.x && base.y == again.y && base.z == again.z),
                               "the same (seed, PID, step) must reproduce the increment after other draws"));
    PetscCall(PicurvAssertBool((PetscBool)(base.x != other_pid.x && base.x != other_step.x && base.x != other_seed.x),
                               "a different PID, step or seed must change the increment"));
    PetscCall(PicurvAssertBool((PetscBool)(none.x == 0.0 && none.y == 0.0 && none.z == 0.0),
                               "zero diffusivity must give no increment"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that increments over many particles have zero mean, variance 2 D dt, and
 *        uncorrelated components.
 */

static PetscErrorCode TestBrownianDisplacementMoments(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    const PetscInt count = 40000;
    const PetscReal diffusivity = 0.5;
    PetscReal sum[3] = {0.0, 0.0, 0.0}, sum_sq[3] = {0.0, 0.0, 0.0}, sum_xy = 0.0, variance;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->dt = 0.25;
    simCtx->particleRandomSeed = 12345;
    simCtx->step = 3;
    variance = 2.0 * diffusivity * simCtx->dt;
    for (PetscInt pid = 0; pid < count; ++pid) {
        Cmpnts d;

        PetscCall(CalculateBrownianDisplacement(user, pid, diffusivity, &d));
        const PetscReal c[3] = {d.x, d.y, d.z};
        for (int a = 0; a < 3; ++a) { sum[a] += c[a]; sum_sq[a] += c[a] * c[a]; }
        sum_xy += d.x * d.y;
    }
    /* Five standard errors: the mean's is sqrt(var/N), the variance's sqrt(2/N) var. */
    for (int a = 0; a < 3; ++a) {
        PetscCall(PicurvAssertRealNear(0.0, sum[a] / count, 5.0 * PetscSqrtReal(variance / count), "Brownian increments must have zero mean"));
        PetscCall(PicurvAssertRealNear(variance, sum_sq[a] / count, 5.0 * PetscSqrtReal(2.0 / count) * variance, "Brownian increments must have variance 2 D dt"));
    }
    PetscCall(PicurvAssertRealNear(0.0, sum_xy / count, 5.0 * variance / PetscSqrtReal((PetscReal)count), "Brownian components must be uncorrelated"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
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
    PetscCall(SetupDomainCellDecompositionMap(user));
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
    PetscCall(SetupDomainCellDecompositionMap(user));
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
    PetscCall(SetupDomainCellDecompositionMap(user));

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

    PetscCall(LocateParticleOrFindMigrationTarget(user, NULL, &particle, &status));

    PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status, "direct re-search test should locate the particle"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.searchAttempts, "direct re-search test should record one robust walk"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.reSearchCount, "direct re-search test should increment re_search_count on later passes"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.maxTraversalFailCount, "direct re-search test should not hit MAX_TRAVERSAL"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that a particle beyond a reflecting -Xi face is mirrored back and located, including after a tangential walk.
 */
static PetscErrorCode TestLocateParticleOrFindMigrationTargetReflectsAtWall(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Particle particle;
    ParticleLocationStatus status = NEEDS_LOCATION;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(SetupDomainCellDecompositionMap(user));
    user->boundary_faces[BC_FACE_NEG_X].particle_policy = PARTICLE_BOUNDARY_REFLECT;
    simCtx->searchMetrics.currentSettlementPass = 1;

    /* Straight crossing of -Xi out of its previous cell (0,1,1). */
    PetscCall(PetscMemzero(&particle, sizeof(particle)));
    particle.PID = 41;
    particle.cell[0] = 0; particle.cell[1] = 1; particle.cell[2] = 1;
    particle.loc.x = -0.05; particle.loc.y = 0.375; particle.loc.z = 0.375;
    PetscCall(LocateParticleOrFindMigrationTarget(user, NULL, &particle, &status));
    PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status, "a reflected particle should be located"));
    PetscCall(PicurvAssertRealNear(0.05, particle.loc.x, 1.0e-14, "x should mirror about the -Xi plane"));
    PetscCall(PicurvAssertIntEqual(0, particle.cell[0], "the particle stays in the boundary cell column"));
    PetscCall(PicurvAssertIntEqual(1, particle.cell[1], "a straight crossing keeps its j cell"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)simCtx->searchMetrics.reflectionCount, "one reflection"));

    /* Crossing -Xi while also moving one cell in +j: the walk moves along the wall first. */
    PetscCall(PetscMemzero(&particle, sizeof(particle)));
    particle.PID = 42;
    particle.cell[0] = 0; particle.cell[1] = 1; particle.cell[2] = 1;
    particle.loc.x = -0.05; particle.loc.y = 0.625; particle.loc.z = 0.375;
    PetscCall(LocateParticleOrFindMigrationTarget(user, NULL, &particle, &status));
    PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status, "a tangential crossing should also be located"));
    PetscCall(PicurvAssertRealNear(0.05, particle.loc.x, 1.0e-14, "x should mirror about the -Xi plane"));
    PetscCall(PicurvAssertIntEqual(2, particle.cell[1], "the reflection happens in the boundary cell that holds the particle tangentially"));
    PetscCall(PicurvAssertIntEqual(2, (PetscInt)simCtx->searchMetrics.reflectionCount, "one reflection per particle"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)simCtx->searchMetrics.reflectionFallbackCount, "no particle reached the cap"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
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
    PetscReal utau_werner = 0.0;
    PetscReal utau_cabot = 0.0;
    PetscReal wall_shear_velocity = 0.0;
    PetscReal wall_shear_normal = 0.0;

    PetscFunctionBeginUser;
    utau_werner = find_utau_Werner(kinematic_viscosity, target_velocity, wall_distance, 0.1);
    PetscCall(PicurvAssertBool((PetscBool)(utau_werner > 0.0), "Werner-Wengle friction velocity should remain positive"));
    PetscCall(PicurvAssertRealNear(target_velocity, u_Werner(kinematic_viscosity, wall_distance, utau_werner), 1.0e-12,
                                   "Werner-Wengle inversion should reconstruct the target velocity"));

    /* The pointwise pair the solver uses must invert exactly on BOTH sides of
       y+ = 11.81. The previous inverse round-tripped only in the sublayer, which is the
       one range the sample above happens to sit in, and so passed over the mismatch. */
    {
        const PetscReal yplus[5] = {2.0, 11.0, 11.81, 30.0, 500.0};
        const PetscReal u_tau = 0.05;
        const PetscReal nu = kinematic_viscosity;

        for (PetscInt i = 0; i < 5; ++i) {
            const PetscReal y = yplus[i] * nu / u_tau;
            const PetscReal u = u_Werner(nu, y, u_tau);
            PetscCall(PicurvAssertRealNear(u_tau, find_utau_Werner(nu, u, y, 0.0), 1.0e-12 * u_tau,
                                           "find_utau_Werner must invert u_Werner exactly on both branches"));
        }
    }

    /* The cell-integral pair: correct for a cell-averaged input, not wired into the
       solver. Its round trip must still hold on both sides of its own switch. */
    {
        const PetscReal A = 8.3, Bexp = 1.0 / 7.0;
        const PetscReal u_m = 0.5 * kinematic_viscosity / wall_distance *
                              PetscPowReal(A, 2.0 / (1.0 - Bexp));
        const PetscReal samples[2] = {0.25 * u_m, 4.0 * u_m};

        for (PetscInt i = 0; i < 2; ++i) {
            const PetscReal tau = taw_Werner(kinematic_viscosity, samples[i], wall_distance);

            PetscCall(PicurvAssertBool((PetscBool)(tau > 0.0),
                                       "the explicit Werner-Wengle stress must be positive"));
            PetscCall(PicurvAssertRealNear(
                samples[i], u_Werner_explicit(kinematic_viscosity, wall_distance, tau),
                1.0e-9 + 1.0e-9 * samples[i],
                "the explicit Werner-Wengle pair must round-trip on both branches"));
        }

        /* The branches meet at the switch, so approaching it from either side must give
           the same stress. A threshold written through the friction velocity, as the
           iterative residual's is, cannot satisfy this. */
        PetscCall(PicurvAssertRealNear(
            taw_Werner(kinematic_viscosity, u_m * (1.0 - 1.0e-9), wall_distance),
            taw_Werner(kinematic_viscosity, u_m * (1.0 + 1.0e-9), wall_distance), 1.0e-9,
            "the Werner-Wengle branches must be continuous at the transition speed"));

        /* Nearer the wall is slower, whatever the branch. This is the property the
           previous mismatched pairing violated. */
        {
            const PetscReal tau = taw_Werner(kinematic_viscosity, 4.0 * u_m, wall_distance);

            PetscCall(PicurvAssertBool(
                (PetscBool)(u_Werner_explicit(kinematic_viscosity, 0.5 * wall_distance, tau) <
                            u_Werner_explicit(kinematic_viscosity, wall_distance, tau)),
                "a cell nearer the wall must be assigned the slower speed"));
        }
    }

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
 * @brief Tests that the wall laws' branches meet, so each reference speed has one friction velocity.
 *
 * The log law switched from u+ = y+ to its logarithm at a fixed y+ = 11.81, where the two
 * differ by about 2.4%, and Werner-Wengle at the same rounded value, where its branches
 * differ by about 2e-5. Either jump makes the speed-to-u_tau map non-monotone near the
 * switch and the residual non-smooth for a finite-difference Jacobian.
 */
static PetscErrorCode TestWallLawBranchesAreContinuous(void)
{
    const double nu = 1.0e-3, y = 1.0e-2;
    const double edge = 11.81;    /* only a probe location; the switches are computed */
    double previous_utau = 0.0;

    PetscFunctionBeginUser;
    /* Profiles: sample u+ across the old switch at a fine step; no step may jump. */
    for (int law = 0; law < 2; ++law) {
        double previous = -1.0;
        for (int n = 0; n <= 400; ++n) {
            const double yplus = edge - 1.0 + 2.0 * n / 400.0;
            const double u_tau = yplus * nu / y;
            const double u = (law == 0) ? u_hydset_roughness(nu, y, u_tau, 0.0) : u_Werner(nu, y, u_tau);
            const double uplus = u / u_tau;
            if (previous > 0.0) {
                PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(uplus - previous) < 1.0e-2),
                                           law == 0 ? "the log-law profile must not jump at its sublayer edge"
                                                    : "the Werner-Wengle profile must not jump at its sublayer edge"));
                PetscCall(PicurvAssertBool((PetscBool)(uplus >= previous - 1.0e-12),
                                           "u+ must not decrease with y+"));
            }
            previous = uplus;
        }
    }
    /* Inverse: across the band of U y / nu where the old log law had two roots, the
       friction velocity must increase with the reference speed and reproduce it. */
    for (int n = 0; n <= 200; ++n) {
        const double U = (130.0 + 15.0 * n / 200.0) * nu / y;
        const double u_tau = find_utau_hydset(nu, U, y, 1.0e-3, 0.0);
        PetscCall(PicurvAssertRealNear(U, u_hydset_roughness(nu, y, u_tau, 0.0), 1.0e-6 * U,
                                       "the log-law friction velocity must reproduce the speed it inverts"));
        PetscCall(PicurvAssertBool((PetscBool)(u_tau > previous_utau), "the log-law friction velocity must increase with speed"));
        previous_utau = u_tau;
    }
    /* Forward and inverse are one relation at every y+, smooth and rough, including above
       y+ = 300 where the profile once returned -1 while the inverse kept the log branch.
       y+ = 100 is the control inside the old range. */
    {
        const double yplus_samples[] = {100.0, 299.0, 301.0, 500.0, 2000.0};
        const double roughness[] = {0.0, 5.0e-3};
        for (int r = 0; r < 2; ++r) {
            for (int n = 0; n < 5; ++n) {
                const double u_tau = yplus_samples[n] * nu / y;
                const double U = u_hydset_roughness(nu, y, u_tau, roughness[r]);
                PetscCall(PicurvAssertBool((PetscBool)(U > 0.0),
                                           "the log-law profile must be a speed at every y+, not a sentinel"));
                PetscCall(PicurvAssertRealNear(u_tau, find_utau_hydset(nu, U, y, 1.0e-3, roughness[r]), 1.0e-6 * u_tau,
                                               "the log-law inverse must recover the friction velocity of its profile"));
            }
        }
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the log law's friction velocity stays positive on a rough wall.
 *
 * Newton from a small guess could step onto the linear branch's mirror root at -u_tau; a
 * rough wall made that common, and the sign then corrupted u_tau, its minimum and y+.
 */
static PetscErrorCode TestLogLawFrictionVelocityStaysPositiveOnRoughWalls(void)
{
    const double nu = 1.0 / 4.0e4;
    const double roughness[4] = {1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2};
    const double distance[3] = {1.0e-3, 1.0e-2, 3.0e-2};
    const double speed[6] = {1.0e-3, 1.0e-2, 0.1, 0.5, 1.0, 2.0};

    PetscFunctionBeginUser;
    for (int r = 0; r < 4; ++r)
        for (int d = 0; d < 3; ++d)
            for (int s = 0; s < 6; ++s) {
                const double u_tau = find_utau_hydset(nu, speed[s], distance[d], 1.0e-3, roughness[r]);
                PetscCall(PicurvAssertBool((PetscBool)(u_tau > 0.0 && isfinite(u_tau)),
                                           "the log-law friction velocity must be positive and finite"));
            }
    PetscFunctionReturn(0);
}

/**
 * @brief Tests Cabot's stress sign and its wall-layer integral at large y+.
 *
 * Cabot's momentum balance gives a signed stress; under an adverse pressure gradient strong
 * enough to exceed the reference speed's share it reverses, and the solver applies it in that
 * direction. Its wall-layer integrals were 30 uniform points, about a third
 * off at y+ = 1000; they must now agree with a resolved reference.
 */
static PetscErrorCode TestCabotSignedStressAndIntegralAccuracy(void)
{
    const double nu = 1.0e-3, y = 1.0e-2;
    double u_tau = 0.0, stress = 0.0, stress_normal = 0.0;

    PetscFunctionBeginUser;
    /* Integral against a resolved trapezoid in log-spaced y (no closed form). */
    for (int mode = 0; mode < 2; ++mode) {
        const double u_tau_probe = 1000.0 * nu / y;   /* y+ = 1000 at the end */
        const int    fine = 200000;
        double       reference_integral = 0.0, y_prev = 0.0, f_prev = 0.0;
        for (int n = 1; n <= fine; ++n) {
            const double yy = y * (double)n * (double)n / ((double)fine * (double)fine);
            const double f = ((mode == 0) ? 1.0 : yy) / ((1.0 + nu_t(yy * u_tau_probe / nu)) * nu);
            if (n == 1) f_prev = ((mode == 0) ? 1.0 : 0.0) / nu;
            reference_integral += 0.5 * (f + f_prev) * (yy - y_prev);
            y_prev = yy; f_prev = f;
        }
        PetscCall(PicurvAssertRealNear(reference_integral, integrate_1(nu, y, u_tau_probe, mode),
                                       2.0e-3 * PetscAbsReal(reference_integral),
                                       "Cabot's wall-layer integral must match a resolved reference at y+ = 1000"));
    }

    /* No pressure gradient: the momentum-balance stress opposes the flow. */
    find_utau_Cabot(nu, 1.0, 2.0 * y, 0.01, 0.0, 0.0, &u_tau, &stress, &stress_normal);
    PetscCall(PicurvAssertBool((PetscBool)(u_tau > 0.0 && stress > 0.0),
                               "with no pressure gradient Cabot's stress opposes the flow"));
    /* A strong adverse gradient along the flow reverses the wall-layer balance; the signed
       stress is what ApplyWallFunction orients the wall stress by. */
    find_utau_Cabot(nu, 1.0, 2.0 * y, 0.01, 1.0e3, 0.0, &u_tau, &stress, &stress_normal);
    PetscCall(PicurvAssertBool((PetscBool)(stress < 0.0),
                               "an adverse gradient that reverses the balance must reverse Cabot's stress"));
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

    PetscCall(UpdateSolverHistoryVectors(user, PETSC_FALSE));

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

    /* A non-periodic single-rank grid: no face, edge or corner neighbour, and itself at the centre. */
    for (PetscInt n = 0; n < 27; ++n) {
        PetscCall(PicurvAssertIntEqual(n == RANK_NEIGHBOR_INDEX(0, 0, 0) ? 0 : MPI_PROC_NULL, user->neighbors.rank[n],
                                       "single-rank neighbours should be null except the rank itself"));
    }
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));

    /* Periodic in x: the rank is its own neighbour across both x faces, and nothing else exists. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    PetscCall(ComputeAndStoreNeighborRanks(user));
    for (PetscInt n = 0; n < 27; ++n) {
        const PetscBool along_x = (PetscBool)(n == RANK_NEIGHBOR_INDEX(-1, 0, 0) || n == RANK_NEIGHBOR_INDEX(0, 0, 0) ||
                                              n == RANK_NEIGHBOR_INDEX(1, 0, 0));
        PetscCall(PicurvAssertIntEqual(along_x ? 0 : MPI_PROC_NULL, user->neighbors.rank[n],
                                       "x-periodic single-rank neighbours should wrap to the rank along x only"));
    }

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
        {"is-particle-in-box-basic-cases", TestIsParticleInBoxBasicCases},
        {"update-particle-weights-computes-expected-ratios", TestUpdateParticleWeightsComputesExpectedRatios},
        {"update-particle-position-without-brownian-contribution", TestUpdateParticlePositionWithoutBrownianContribution},
        {"update-particle-position-diffusivity-gradient-only", TestUpdateParticlePositionDiffusivityGradientOnly},
        {"update-particle-field-iem-relaxation", TestUpdateParticleFieldIEMRelaxation},
        {"cart2contra-converts-cartesian-field", TestCart2ContraConvertsCartesianField},
        {"cart2contra-uses-finalized-periodic-ucat", TestCart2ContraUsesFinalizedPeriodicUcat},
        {"populate-initial-ucont-loads-staged-ucat", TestPopulateInitialUcontLoadsStagedUcat},
        {"write-les-fields-copies-owned-periodic-values", TestWriteLESFieldsCopiesOwnedPeriodicValues},
        {"staged-ucat-fourier-modes-use-picurv-discrete-symbol", TestStagedUcatFourierModesUsePicurvDiscreteSymbol},
        {"populate-initial-ucont-loads-staged-ucont", TestPopulateInitialUcontLoadsStagedUcont},
        {"interpolate-all-fields-to-swarm-constant-fields", TestInterpolateAllFieldsToSwarmConstantFields},
        {"interpolate-all-fields-to-swarm-corner-averaged-constant-fields", TestInterpolateAllFieldsToSwarmCornerAveragedConstantFields},
        {"interpolate-all-fields-to-swarm-face-flux-constant-fields", TestInterpolateAllFieldsToSwarmFaceFluxConstantFields},
        {"face-flux-velocity-is-divergence-free", TestFaceFluxVelocityIsDivergenceFree},
        {"scatter-all-particle-fields-to-euler-fields-averages-psi", TestScatterAllParticleFieldsToEulerFieldsAveragesPsi},
        {"calculate-particle-count-per-cell-counts-global-cell-ids", TestCalculateParticleCountPerCellCountsGlobalCellIDs},
        {"reset-all-particle-statuses-leaves-lost-particles-untouched", TestResetAllParticleStatusesLeavesLostParticlesUntouched},
        {"check-and-remove-lost-particles-removes-lost-entries", TestCheckAndRemoveLostParticlesRemovesLostEntries},
        {"sort-particles-by-cell-keeps-each-point-whole", TestSortParticlesByCellKeepsEachPointWhole},
        {"brownian-displacement-keyed-on-seed-particle-step", TestBrownianDisplacementIsKeyedOnSeedParticleAndStep},
        {"brownian-displacement-moments", TestBrownianDisplacementMoments},
        {"update-all-particle-positions-moves-swarm-entries", TestUpdateAllParticlePositionsMovesSwarmEntries},
        {"locate-particle-or-find-migration-target-reflects-at-wall", TestLocateParticleOrFindMigrationTargetReflectsAtWall},
        {"locate-all-particles-in-grid-prior-cell-fast-path", TestLocateAllParticlesInGridPriorCellFastPath},
        {"locate-all-particles-in-grid-guess-path-resolves-local-particle", TestLocateAllParticlesInGridGuessPathResolvesLocalParticle},
        {"locate-particle-or-find-migration-target-counts-research", TestLocateParticleOrFindMigrationTargetCountsReSearch},
        {"wall-law-branches-are-continuous", TestWallLawBranchesAreContinuous},
        {"log-law-friction-velocity-stays-positive-on-rough-walls", TestLogLawFrictionVelocityStaysPositiveOnRoughWalls},
        {"cabot-signed-stress-and-integral-accuracy", TestCabotSignedStressAndIntegralAccuracy},
        {"wall-model-scalar-helpers", TestWallModelScalarHelpers},
        {"wall-model-velocity-helpers", TestWallModelVelocityHelpers},
        {"validate-driven-flow-configuration-no-driven-handlers", TestValidateDrivenFlowConfigurationNoDrivenHandlers},
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
