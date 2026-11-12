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
 * The boundary values are this module's own `Ubcs`, so after setting the state it runs the part of
 * the solved flow's finalization that does not refresh them from the boundary handlers
 * (FinalizeCellFieldsFromBoundaryValues()).
 */

#include "AnalyticalSolutions.h"
#include "expression.h"
#include "Boundaries.h"
#include "Metric.h"
#include "logging.h"
#include "setup.h"

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
    PicurvSourceFaceSlots (*faces)[6]; /* Per block, what the source reads on each face. */
};

/** @brief Read the source once, and work out what it depends on. Collective. */
static PetscErrorCode AnalyticalFieldCreate(SimCtx *simCtx)
{
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    struct AnalyticalField *field;
    PetscBool uses;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PetscNew(&field));
    PetscCall(PicurvVelocitySourceCreate("ic_source", PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, &field->source));
    PetscCall(PicurvSourcePressureCreate("ic_source", PETSC_FALSE, &field->pressure));
    if (field->pressure) {
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
    PetscCall(PetscCalloc3(simCtx->block_number, &field->fraction, simCtx->block_number, &field->length,
                           simCtx->block_number, &field->faces));
    for (PetscInt bi = 0; bi < simCtx->block_number; ++bi) {
        PetscCall(DMCreateGlobalVector(user_finest[bi].fda, &field->fraction[bi]));
        PetscCall(DMCreateGlobalVector(user_finest[bi].fda, &field->length[bi]));
        PetscCall(ComputeGridLineFractions(&user_finest[bi], field->fraction[bi], field->length[bi]));
        for (PetscInt face = 0; face < 6; ++face) {
            PetscCall(PicurvSourceFaceSlotsCreate(&user_finest[bi], (BCFace)face, field->fraction[bi], field->length[bi],
                                                  &field->faces[bi][face]));
        }
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
            for (PetscInt face = 0; face < 6; ++face) PetscCall(PicurvSourceFaceSlotsDestroy(&field->faces[bi][face]));
        }
        PetscCall(PetscFree3(field->fraction, field->length, field->faces));
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
    PROFILE_HELPER_BEGIN;
    const PetscBool pointwise = (simCtx && simCtx->analyticalField) ? simCtx->analyticalField->pointwise : PETSC_FALSE;
    PROFILE_HELPER_END;
    return pointwise;
}

/** @brief Set Ubcs on every physical boundary face to the source at the face centre. */
static PetscErrorCode SetBoundaryVelocityFromSource(UserCtx *user, const struct AnalyticalField *field, PetscInt block,
                                                    PetscReal physical_time)
{
    const PetscReal U_ref = user->simCtx->scaling.U_ref;
    Cmpnts ***ubcs;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs));
    for (PetscInt face = 0; face < 6; ++face) {
        const PicurvSourceFaceSlots *slots = &field->faces[block][face];
        for (PetscInt q = 0; q < slots->elements.n[1]; ++q)
        for (PetscInt p = 0; p < slots->elements.n[0]; ++p) {
            const PetscInt s = p + slots->elements.n[0] * q;
            PetscInt f[3];
            Cmpnts v;
            /* A + face's boundary value sits one past its flux, in the dummy layer. */
            f[slots->elements.normal] = slots->elements.layer + face % 2; f[slots->elements.a] = slots->elements.lo[0] + p; f[slots->elements.b] = slots->elements.lo[1] + q;
            slots->values[s][PICURV_SOURCE_VAR_T] = physical_time;
            PetscCall(PicurvVelocitySourceEvaluate(&field->source, slots->values[s], NULL, slots->edges[s], &v));
            ubcs[f[2]][f[1]][f[0]].x = field->factor * v.x / U_ref;
            ubcs[f[2]][f[1]][f[0]].y = field->factor * v.y / U_ref;
            ubcs[f[2]][f[1]][f[0]].z = field->factor * v.z / U_ref;
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs));
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
        /* A curvilinear source's Ucat follows from its fluxes, below. */
        PetscCall(PicurvSourceSetCellFields(user, (field->source.frame != PICURV_SOURCE_FRAME_CURVILINEAR) ? &field->source : NULL,
                                            field->factor, field->pressure, field->fraction[bi], field->length[bi], physical_time));
        if (!field->pressure) PetscCall(VecSet(user->P, 0.0));
        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered_fields));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCONT));
        if (field->source.frame == PICURV_SOURCE_FRAME_CURVILINEAR) PetscCall(Contra2Cart(user));
        PetscCall(FinalizeCellFieldsFromBoundaryValues(user));
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
        PetscCall(PicurvVelocitySourceEvaluate(&field->source, values, NULL, NULL, &v));
        data[i] = v.x / scaling->U_ref;
        data[i + 1] = v.y / scaling->U_ref;
        data[i + 2] = v.z / scaling->U_ref;
    }
    PetscCall(VecRestoreArray(tempVec, &data));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
