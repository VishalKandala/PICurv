/**
 * @file AnalyticalSolutions.c
 * @brief Prescribes the Eulerian fields from the case's source instead of solving for them.
 *
 * @details With `eulerian_field_source: analytical`, the field is the initial condition's source,
 * evaluated at the step's time through the operations every target of a source shares: fluxes at
 * the faces (PicurvSourceAddFluxes()) and a bulk velocity (PicurvSourceCellBulkFactor()). A Cartesian
 * source also sets `Ucat` and an optional pressure at the cell centres, exactly; a curvilinear one
 * derives `Ucat` from its fluxes.
 *
 * The boundary values are this module's own `Ubcs`. After setting the state it runs the sequence the
 * solved flow is finalized with (see `FinalizePostProjectionCellFields`): `UpdateDummyCells` on
 * non-periodic faces, `SynchronizePeriodicCellFields` on periodic ones, `UpdateCornerNodes` for edges
 * and corners, and the periodic synchronization again because corner averaging can overwrite
 * periodic endpoints. It does not call that finalizer, which would refresh `Ubcs` from the boundary
 * handlers.
 */

#include "AnalyticalSolutions.h"     // The header that declares this file's functions
#include "expression.h"
#include "verification_sources.h"
#include "Metric.h"
#include <petscmath.h>               // Provides PETSc-compatible math functions like sin, cos, exp

/** @brief Cell-centred fields the analytical state finalizes, in the order the post-projection
 *         finalization uses for the solved flow. */
static const FieldId kAnalyticalCellFields[] = {FIELD_ID_UCAT, FIELD_ID_P};

/** @brief The prescribed field: the case's source, and what it reuses from step to step. */
struct AnalyticalField {
    PicurvVelocitySource source;   /* `-ic_source_*`. */
    PicurvExpression    *pressure; /* `-ic_source_pressure` [Pa], or NULL for zero pressure. */
    PetscBool            time_dependent;
    PetscBool            pointwise; /* Reads only x, y, z, t in the Cartesian frame, and varies in space. */
    PetscBool            evaluated;
    PetscReal            factor;   /* From the bulk velocity. */
    Vec                 *fraction; /* Per block. */
    Vec                 *length;
};

/** @brief Read the source once, and work out what it depends on. Collective. */
static PetscErrorCode AnalyticalFieldCreate(SimCtx *simCtx)
{
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    struct AnalyticalField *field;
    char      text[4096];
    PetscBool found = PETSC_FALSE, uses;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PetscNew(&field));
    PetscCall(PicurvVelocitySourceCreate("ic_source", PETSC_FALSE, PETSC_FALSE, &field->source));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-ic_source_pressure", text, sizeof(text), &found));
    if (found) {
        PetscCall(PicurvExpressionCompile(text, kPicurvSourceVariableNames, PICURV_SOURCE_VAR_COUNT, PETSC_FALSE, &field->pressure));
        PetscCall(PicurvExpressionUsesVariable(field->pressure, PICURV_SOURCE_VAR_T, &uses));
        field->time_dependent = uses;
    }
    field->time_dependent = (PetscBool)(field->time_dependent || field->source.time_dependent ||
                                        field->source.bulk_time_dependent);
    field->pointwise = (PetscBool)(field->source.frame == PICURV_SOURCE_FRAME_CARTESIAN && !field->source.bulk);
    {
        PetscBool varies = PETSC_FALSE;
        for (PetscInt c = 0; c < field->source.count; ++c) {
            for (PetscInt var = 0; var < PICURV_SOURCE_VAR_COUNT; ++var) {
                PetscCall(PicurvExpressionUsesVariable(field->source.component[c], var, &uses));
                if (!uses) continue;
                if (var >= PICURV_SOURCE_VAR_XI) field->pointwise = PETSC_FALSE;
                if (var <= PICURV_SOURCE_VAR_Z) varies = PETSC_TRUE;
            }
        }
        field->pointwise = (PetscBool)(field->pointwise && varies);
    }
    field->factor = 1.0;
    PetscCall(PetscCalloc2(simCtx->block_number, &field->fraction, simCtx->block_number, &field->length));
    for (PetscInt bi = 0; bi < simCtx->block_number; ++bi) {
        PetscCall(DMCreateGlobalVector(user_finest[bi].fda, &field->fraction[bi]));
        PetscCall(DMCreateGlobalVector(user_finest[bi].fda, &field->length[bi]));
        PetscCall(ComputeGridLineFractions(&user_finest[bi], field->fraction[bi], field->length[bi]));
    }
    simCtx->analyticalField = field;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref AnalyticalFieldDestroy().
 * @see AnalyticalFieldDestroy()
 */
PetscErrorCode AnalyticalFieldDestroy(SimCtx *simCtx)
{
    struct AnalyticalField *field = simCtx->analyticalField;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (field) {
        PetscCall(PicurvVelocitySourceDestroy(&field->source));
        PetscCall(PicurvExpressionDestroy(&field->pressure));
        for (PetscInt bi = 0; bi < simCtx->block_number; ++bi) {
            PetscCall(VecDestroy(&field->fraction[bi]));
            PetscCall(VecDestroy(&field->length[bi]));
        }
        PetscCall(PetscFree2(field->fraction, field->length));
        PetscCall(PetscFree(simCtx->analyticalField));
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref AnalyticalFieldSupportsInterpolationError().
 * @see AnalyticalFieldSupportsInterpolationError()
 */
PetscBool AnalyticalFieldSupportsInterpolationError(const SimCtx *simCtx)
{
    PROFILE_FUNCTION_BEGIN;
    PROFILE_FUNCTION_END;
    return (simCtx && simCtx->analyticalField) ? simCtx->analyticalField->pointwise : PETSC_FALSE;
}

/** @brief Set Ubcs on every physical boundary face to the source at the face centre. */
static PetscErrorCode SetBoundaryVelocityFromSource(UserCtx *user, const struct AnalyticalField *field, PetscInt block,
                                                    PetscReal physical_time)
{
    const DMDALocalInfo *info = &user->info;
    const PetscReal L_ref = user->simCtx->scaling.L_ref, U_ref = user->simCtx->scaling.U_ref;
    const PetscInt start[3] = {info->xs, info->ys, info->zs};
    const PetscInt end[3] = {info->xs + info->xm, info->ys + info->ym, info->zs + info->zm};
    const PetscInt nodes[3] = {info->mx, info->my, info->mz};
    const Cmpnts ***frac, ***len, ***csi, ***eta, ***zet, ***cent_face[3];
    const PetscReal ***aj;
    Cmpnts ***ubcs;
    PetscReal values[PICURV_SOURCE_VAR_COUNT];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(DMDAVecGetArrayRead(user->fda, field->fraction[block], &frac));
    PetscCall(DMDAVecGetArrayRead(user->fda, field->length[block], &len));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, &csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, &eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, &zet));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCentx, &cent_face[0]));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCenty, &cent_face[1]));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCentz, &cent_face[2]));
    PetscCall(DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs));
    values[PICURV_SOURCE_VAR_T] = physical_time;
    for (PetscInt normal = 0; normal < 3; ++normal) {
        const PetscInt a = (normal + 1) % 3, b = (normal + 2) % 3;
        for (PetscInt side = 0; side < 2; ++side) {
            const PetscInt layer = side ? end[normal] - 1 : start[normal];
            if (side ? end[normal] != nodes[normal] : start[normal] != 0) continue;
            for (PetscInt q = start[b]; q < end[b]; ++q)
            for (PetscInt p = start[a]; p < end[a]; ++p) {
                PetscInt f[3], c[3];
                Cmpnts edges[3], v;
                f[normal] = layer; f[a] = p; f[b] = q;
                /* Grid-line data lives on interior cells; rows of dummies read the nearest one. */
                c[normal] = side ? nodes[normal] - 2 : 1;
                c[a] = PetscMax(PetscMin(p, nodes[a] - 2), 1); c[b] = PetscMax(PetscMin(q, nodes[b] - 2), 1);
                PicurvSourceSetVariables(values, cent_face[normal][f[2]][f[1]][f[0]], frac[c[2]][c[1]][c[0]],
                                         len[c[2]][c[1]][c[0]], L_ref);
                values[PICURV_SOURCE_VAR_XI + normal] = side ? 1.0 : 0.0;
                PetscCall(ComputeCellEdgeVectors(aj[c[2]][c[1]][c[0]], csi[c[2]][c[1]][c[0]], eta[c[2]][c[1]][c[0]],
                                                 zet[c[2]][c[1]][c[0]], edges));
                PetscCall(PicurvVelocitySourceEvaluate(&field->source, values, edges, &v));
                ubcs[f[2]][f[1]][f[0]].x = field->factor * v.x / U_ref;
                ubcs[f[2]][f[1]][f[0]].y = field->factor * v.y / U_ref;
                ubcs[f[2]][f[1]][f[0]].z = field->factor * v.z / U_ref;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCentz, &cent_face[2]));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCenty, &cent_face[1]));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCentx, &cent_face[0]));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, field->length[block], &len));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, field->fraction[block], &frac));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/** @brief Set a block's Ucat and P at every interior cell centre from the source. */
static PetscErrorCode SetCellFieldsFromSource(UserCtx *user, const struct AnalyticalField *field, PetscInt block,
                                              PetscReal physical_time)
{
    const DMDALocalInfo *info = &user->info;
    const ScalingCtx *scaling = &user->simCtx->scaling;
    const PetscInt xs = info->xs, xe = info->xs + info->xm, mx = info->mx;
    const PetscInt ys = info->ys, ye = info->ys + info->ym, my = info->my;
    const PetscInt zs = info->zs, ze = info->zs + info->zm, mz = info->mz;
    const PetscBool set_velocity = (PetscBool)(field->source.frame != PICURV_SOURCE_FRAME_CURVILINEAR);
    const Cmpnts ***cent, ***frac, ***len;
    Cmpnts ***ucat;
    PetscReal ***p, values[PICURV_SOURCE_VAR_COUNT];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Cent, &cent));
    PetscCall(DMDAVecGetArrayRead(user->fda, field->fraction[block], &frac));
    PetscCall(DMDAVecGetArrayRead(user->fda, field->length[block], &len));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArray(user->da, user->P, &p));
    values[PICURV_SOURCE_VAR_T] = physical_time;
    for (PetscInt k = (zs == 0) ? 1 : zs; k < ((ze == mz) ? mz - 1 : ze); ++k)
    for (PetscInt j = (ys == 0) ? 1 : ys; j < ((ye == my) ? my - 1 : ye); ++j)
    for (PetscInt i = (xs == 0) ? 1 : xs; i < ((xe == mx) ? mx - 1 : xe); ++i) {
        PicurvSourceSetVariables(values, cent[k][j][i], frac[k][j][i], len[k][j][i], scaling->L_ref);
        if (set_velocity) {
            Cmpnts v;
            PetscCall(PicurvVelocitySourceEvaluate(&field->source, values, NULL, &v));
            ucat[k][j][i].x = field->factor * v.x / scaling->U_ref;
            ucat[k][j][i].y = field->factor * v.y / scaling->U_ref;
            ucat[k][j][i].z = field->factor * v.z / scaling->U_ref;
        }
        if (field->pressure) {
            PetscReal value;
            PetscCall(PicurvExpressionEvaluate(field->pressure, values, NULL, &value));
            p[k][j][i] = value / scaling->P_ref;
        } else {
            p[k][j][i] = 0.0;
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, field->length[block], &len));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, field->fraction[block], &frac));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref AnalyticalSolutionEngine().
 * @details The field is the case's source, evaluated at the step's time. As on an inlet, a source
 *          that does not read `t` is evaluated once.
 * @see AnalyticalSolutionEngine()
 */
PetscErrorCode AnalyticalSolutionEngine(SimCtx *simCtx)
{
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    struct AnalyticalField *field;
    PetscReal physical_time;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!simCtx->analyticalField) PetscCall(AnalyticalFieldCreate(simCtx));
    field = simCtx->analyticalField;
    if (field->evaluated && !field->time_dependent) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
    for (PetscInt bi = 0; bi < simCtx->block_number; ++bi) {
        UserCtx *user = &user_finest[bi];
        const FieldId staggered_fields[] = {FIELD_ID_UCONT};
        if (field->source.bulk) {
            PetscCall(VecSet(user->Ucat, 0.0));
            PetscCall(PicurvSourceCellBulkFactor(user, &field->source, field->fraction[bi], field->length[bi], physical_time, &field->factor));
        }
        PetscCall(VecSet(user->Ucont, 0.0));
        PetscCall(PicurvSourceAddFluxes(user, &field->source, field->fraction[bi], field->length[bi], physical_time, field->factor, PETSC_TRUE));
        PetscCall(SetBoundaryVelocityFromSource(user, field, bi, physical_time));
        PetscCall(SetCellFieldsFromSource(user, field, bi, physical_time));
        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered_fields));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCONT));
        if (field->source.frame == PICURV_SOURCE_FRAME_CURVILINEAR) PetscCall(Contra2Cart(user));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_P));
        PetscCall(UpdateDummyCells(user, FIELD_ID_UCAT));
        PetscCall(UpdateDummyCells(user, FIELD_ID_P));
        PetscCall(SynchronizePeriodicCellFields(user, 2, kAnalyticalCellFields));
        PetscCall(UpdateCornerNodes(user, FIELD_ID_UCAT));
        PetscCall(UpdateCornerNodes(user, FIELD_ID_P));
        PetscCall(SynchronizePeriodicCellFields(user, 2, kAnalyticalCellFields));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_P));
    }
    field->evaluated = PETSC_TRUE;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SetAnalyticalSolutionForParticles().
 * @details Only a pointwise source reaches here; see AnalyticalFieldSupportsInterpolationError().
 * @see SetAnalyticalSolutionForParticles()
 */
PetscErrorCode SetAnalyticalSolutionForParticles(Vec tempVec, SimCtx *simCtx)
{
    const struct AnalyticalField *field = simCtx->analyticalField;
    const ScalingCtx *scaling = &simCtx->scaling;
    PetscInt nLocal;
    PetscReal *data, values[PICURV_SOURCE_VAR_COUNT], physical_time;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(AnalyticalFieldSupportsInterpolationError(simCtx), PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "The analytical field cannot be evaluated at a particle: it needs a Cartesian source in x, y, z and t.");
    PetscCall(PetscArrayzero(values, PICURV_SOURCE_VAR_COUNT));
    PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
    values[PICURV_SOURCE_VAR_T] = physical_time;
    PetscCall(VecGetLocalSize(tempVec, &nLocal));
    PetscCall(VecGetArray(tempVec, &data));
    /* Interleaved [x0,y0,z0, x1,y1,z1, ...], replaced by the velocity there. */
    for (PetscInt i = 0; i < nLocal; i += 3) {
        Cmpnts v;
        values[PICURV_SOURCE_VAR_X] = data[i] * scaling->L_ref;
        values[PICURV_SOURCE_VAR_Y] = data[i + 1] * scaling->L_ref;
        values[PICURV_SOURCE_VAR_Z] = data[i + 2] * scaling->L_ref;
        PetscCall(PicurvVelocitySourceEvaluate(&field->source, values, NULL, &v));
        data[i] = v.x / scaling->U_ref;
        data[i + 1] = v.y / scaling->U_ref;
        data[i + 2] = v.z / scaling->U_ref;
    }
    PetscCall(VecRestoreArray(tempVec, &data));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SetAnalyticalScalarFieldOnParticles().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see SetAnalyticalScalarFieldOnParticles()
 */
PetscErrorCode SetAnalyticalScalarFieldOnParticles(UserCtx *user, ParticleFieldId particle_field_id)
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
               "Analytical scalar assignment requires a one-component PETSC_REAL particle field; '%s' has %" PetscInt_FMT " components and type %s.",
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
 * @brief Implementation of \ref SetAnalyticalScalarFieldAtCellCenters().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see SetAnalyticalScalarFieldAtCellCenters()
 */
PetscErrorCode SetAnalyticalScalarFieldAtCellCenters(UserCtx *user, Vec targetVec)
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
