#include "verification_sources.h"

#include "expression.h"
#include "logging.h"
#include "setup.h"

/** @brief Names a verification field's expression reads, in the order of its values. */
static const char *const kVerificationVariableNames[4] = {"x", "y", "z", "t"};

/**
 * @brief Compiles a verification field's expression and enables the field.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode VerificationFieldCompile(const char *text, VerificationFieldConfig *field)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PicurvExpressionDestroy(&field->expression));
    PetscCall(PicurvExpressionCompile(text, kVerificationVariableNames, 4, PETSC_FALSE, &field->expression));
    field->enabled = PETSC_TRUE;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Reads and compiles the configured verification fields.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode VerificationFieldsCreate(SimCtx *simCtx)
{
    char      text[4096];
    PetscBool found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-verification_diffusivity_expr", text, sizeof(text), &found));
    if (found) PetscCall(VerificationFieldCompile(text, &simCtx->verificationDiffusivity));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-verification_scalar_expr", text, sizeof(text), &found));
    if (found) PetscCall(VerificationFieldCompile(text, &simCtx->verificationScalar));
    PetscCheck(!simCtx->verificationDiffusivity.enabled || strcmp(simCtx->eulerianSource, "analytical") == 0,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               "verification diffusivity overrides require -euler_field_source \"analytical\".");
    PetscCheck(!simCtx->verificationScalar.enabled || strcmp(simCtx->eulerianSource, "analytical") == 0,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               "verification scalar overrides require -euler_field_source \"analytical\".");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Releases the compiled verification fields.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode VerificationFieldsDestroy(SimCtx *simCtx)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PicurvExpressionDestroy(&simCtx->verificationDiffusivity.expression));
    PetscCall(PicurvExpressionDestroy(&simCtx->verificationScalar.expression));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Evaluates a verification field at a solver-unit point and time.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode VerificationFieldEvaluate(const SimCtx *simCtx, const VerificationFieldConfig *field, Cmpnts point,
                                         PetscReal solver_time, PetscReal *value)
{
    const PetscReal L_ref = simCtx->scaling.L_ref;
    PetscReal values[4] = {point.x * L_ref, point.y * L_ref, point.z * L_ref, 0.0};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PicurvPhysicalTime(simCtx, solver_time, &values[3]));
    PetscCall(PicurvExpressionEvaluate(field->expression, values, NULL, value));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Reports whether the verification-only diffusivity override is enabled.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscBool VerificationDiffusivityOverrideActive(const SimCtx *simCtx)
{
    PROFILE_FUNCTION_BEGIN;
    if (!simCtx) { PROFILE_FUNCTION_END; return PETSC_FALSE; }
    PROFILE_FUNCTION_END;
    return simCtx->verificationDiffusivity.enabled;
}

/**
 * @brief Reports whether the verification-only scalar override is enabled.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscBool VerificationScalarOverrideActive(const SimCtx *simCtx)
{
    PROFILE_FUNCTION_BEGIN;
    if (!simCtx) { PROFILE_FUNCTION_END; return PETSC_FALSE; }
    PROFILE_FUNCTION_END;
    return simCtx->verificationScalar.enabled;
}

/**
 * @brief Fills the Eulerian diffusivity field from the configured verification profile.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode ApplyVerificationDiffusivityOverride(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx        *simCtx = user->simCtx;
    PetscReal    ***diff_arr = NULL;
    const Cmpnts ***cent = NULL;
    DMDALocalInfo  info = user->info;
    PetscReal      min_gamma = PETSC_MAX_REAL;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // Leave the production diffusivity untouched unless this verification mode is selected.
    if (!VerificationDiffusivityOverrideActive(simCtx)) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    ierr = DMDAVecGetArray(user->da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);

    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                PetscReal gamma;
                ierr = VerificationFieldEvaluate(simCtx, &simCtx->verificationDiffusivity, cent[k][j][i], simCtx->ti, &gamma); CHKERRQ(ierr);
                gamma /= simCtx->scaling.U_ref * simCtx->scaling.L_ref;
                diff_arr[k][j][i] = gamma;
                min_gamma = PetscMin(min_gamma, gamma);
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);

    {
        PetscReal global_min_gamma = 0.0;
        ierr = MPI_Allreduce(&min_gamma, &global_min_gamma, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD); CHKERRMPI(ierr);
        // Diffusivity must remain positive on every rank for the PDE operator to be valid.
        if (global_min_gamma <= 0.0) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "verification diffusivity override produced non-positive Gamma values (min=%g).",
                    (double)global_min_gamma);
        }
    }

    {
        const FieldId periodic_fields[] = {FIELD_ID_DIFFUSIVITY};
        // Make the analytic field coherent across periodic seams before ghost updates.
        ierr = SynchronizePeriodicCellFields(user, 1, periodic_fields); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY); CHKERRQ(ierr);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applied the verification diffusivity.\n");

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Fills the particle `Psi` field from the configured verification scalar profile.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode ApplyVerificationScalarOverrideToParticles(UserCtx *user)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    // The normal particle model owns Psi unless the verification override is enabled.
    if (!VerificationScalarOverrideActive(user->simCtx)) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    PetscCall(SetVerificationScalarOnParticles(user, PARTICLE_FIELD_ID_PSI));
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applied the verification scalar to particle field Psi.\n");
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SetVerificationScalarOnParticles().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/verification_sources.h`.
 * @see SetVerificationScalarOnParticles()
 */
PetscErrorCode SetVerificationScalarOnParticles(UserCtx *user, ParticleFieldId particle_field_id)
{
    PetscErrorCode ierr;
    PetscInt       nlocal = 0;
    PetscReal     *positions = NULL;
    PetscReal     *scalar_values = NULL;
    const ParticleFieldDescriptor *descriptor = NULL;
    const char     *swarm_field_name = NULL;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    if (!user->swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx->swarm is NULL.");

    ierr = ParticleFieldGetDescriptor(particle_field_id, &descriptor); CHKERRQ(ierr);
    PetscCheck(descriptor->components == 1 && descriptor->data_type == PETSC_REAL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP,
               "Verification scalar assignment requires a one-component PETSC_REAL particle field; '%s' has %" PetscInt_FMT " components and type %s.",
               descriptor->canonical_name, descriptor->components, PetscDataTypes[descriptor->data_type]);
    swarm_field_name = descriptor->canonical_name;

    ierr = DMSwarmGetLocalSize(user->swarm, &nlocal); CHKERRQ(ierr);
    if (nlocal == 0) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    ierr = DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(user->swarm, swarm_field_name, NULL, NULL, (void **)&scalar_values); CHKERRQ(ierr);

    for (PetscInt p = 0; p < nlocal; ++p) {
        PetscReal value = 0.0;
        const Cmpnts point = {positions[3 * p + 0], positions[3 * p + 1], positions[3 * p + 2]};
        ierr = VerificationFieldEvaluate(user->simCtx, &user->simCtx->verificationScalar, point, user->simCtx->ti, &value); CHKERRQ(ierr);
        scalar_values[p] = value;
    }

    ierr = DMSwarmRestoreField(user->swarm, swarm_field_name, NULL, NULL, (void **)&scalar_values); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&positions); CHKERRQ(ierr);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SetVerificationScalarAtCellCenters().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/verification_sources.h`.
 * @see SetVerificationScalarAtCellCenters()
 */
PetscErrorCode SetVerificationScalarAtCellCenters(UserCtx *user, Vec targetVec)
{
    PetscErrorCode  ierr;
    PetscReal     ***target = NULL;
    const Cmpnts  ***cent = NULL;
    DMDALocalInfo   info;
    PetscInt        xs, xe, ys, ye, zs, ze, mx, my, mz;
    PetscInt        lxs, lxe, lys, lye, lzs, lze;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    if (!targetVec) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "targetVec cannot be NULL.");

    info = user->info;
    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;
    lxs = (xs == 0) ? xs + 1 : xs; lxe = (xe == mx) ? xe - 1 : xe;
    lys = (ys == 0) ? ys + 1 : ys; lye = (ye == my) ? ye - 1 : ye;
    lzs = (zs == 0) ? zs + 1 : zs; lze = (ze == mz) ? ze - 1 : ze;

    ierr = VecSet(targetVec, 0.0); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, targetVec, &target); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);

    for (PetscInt k = lzs; k < lze; ++k) {
        for (PetscInt j = lys; j < lye; ++j) {
            for (PetscInt i = lxs; i < lxe; ++i) {
                PetscReal value = 0.0;
                ierr = VerificationFieldEvaluate(user->simCtx, &user->simCtx->verificationScalar, cent[k][j][i], user->simCtx->ti, &value); CHKERRQ(ierr);
                target[k][j][i] = value;
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, targetVec, &target); CHKERRQ(ierr);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
