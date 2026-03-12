/**
 * @file test_runtime_kernels.c
 * @brief C test module for PICurv.
 */

#include "test_support.h"

#include "BC_Handlers.h"
#include "ParticleMotion.h"
#include "ParticlePhysics.h"
#include "ParticleSwarm.h"
#include "initialcondition.h"
#include "les.h"
#include "runloop.h"
#include "setup.h"
#include "wallfunction.h"
/**
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
 */

static PetscErrorCode TestUpdateSolverHistoryVectorsShiftsStates(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));

    PetscCall(DMCreateGlobalVector(user->fda, &user->Ucont_o));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Ucont_rm1));
    PetscCall(DMCreateLocalVector(user->fda, &user->lUcont_o));
    PetscCall(DMCreateLocalVector(user->fda, &user->lUcont_rm1));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Ucat_o));
    PetscCall(DMCreateGlobalVector(user->da, &user->P_o));

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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Entry point for this unit-test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"distribute-particles-remainder-handling", TestDistributeParticlesRemainderHandling},
        {"is-particle-inside-bbox-basic-cases", TestIsParticleInsideBoundingBoxBasicCases},
        {"update-particle-weights-computes-expected-ratios", TestUpdateParticleWeightsComputesExpectedRatios},
        {"update-particle-position-without-brownian-contribution", TestUpdateParticlePositionWithoutBrownianContribution},
        {"update-particle-field-iem-relaxation", TestUpdateParticleFieldIEMRelaxation},
        {"set-initial-interior-field-ignores-non-ucont-request", TestSetInitialInteriorFieldIgnoresNonUcontRequest},
        {"set-initial-interior-field-constant-profile-on-z-inlet", TestSetInitialInteriorFieldConstantProfileOnZInlet},
        {"wall-noslip-and-freeslip-helpers", TestWallNoSlipAndFreeSlipHelpers},
        {"wall-model-scalar-helpers", TestWallModelScalarHelpers},
        {"validate-driven-flow-configuration-no-driven-handlers", TestValidateDrivenFlowConfigurationNoDrivenHandlers},
        {"compute-smagorinsky-constant-constant-model", TestComputeSmagorinskyConstantConstantModel},
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
