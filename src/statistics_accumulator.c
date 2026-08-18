/**
 * @file statistics_accumulator.c
 * @brief Per-window PETSc accumulator storage and pointwise application.
 *
 * Full API contract is documented with the declarations in
 * `include/statistics_accumulator.h`.
 */

#include "statistics_accumulator.h"
#include "statistics_moments.h"
#include "statistics_target.h"
#include "field_catalog.h"

/** @brief Upper-triangular row-major component pairs for a three-vector self-product. */
static const PetscInt kProductFirst[6]  = {0, 0, 0, 1, 1, 2};
static const PetscInt kProductSecond[6] = {0, 1, 2, 1, 2, 2};

/**
 * @brief Implementation of \ref PicurvProductComponentCount().
 * @see PicurvProductComponentCount()
 */
PetscErrorCode PicurvProductComponentCount(PetscInt dof, PetscInt *count)
{
    PetscFunctionBeginUser;
    PetscCheck(count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Count output is required.");
    PetscCheck(dof == 1 || dof == 3, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Self-products are supported for scalar and three-vector fields, got dof %" PetscInt_FMT ".", dof);
    *count = (dof == 1) ? 1 : 6;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvCovarianceComponentCount().
 * @see PicurvCovarianceComponentCount()
 */
PetscErrorCode PicurvCovarianceComponentCount(PetscInt dof_a, PetscInt dof_b, PetscInt *count)
{
    PetscFunctionBeginUser;
    PetscCheck(count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Count output is required.");
    /* Phase 2 accepts scalar-scalar and vector-scalar pairs; vector-vector cross
     * products are an explicit non-goal. */
    PetscCheck((dof_a == 1 && dof_b == 1) || (dof_a == 3 && dof_b == 1) || (dof_a == 1 && dof_b == 3),
               PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Covariance is supported for scalar-scalar and vector-scalar pairs, got dof %" PetscInt_FMT
               " and %" PetscInt_FMT ".", dof_a, dof_b);
    *count = (dof_a == 3 || dof_b == 3) ? 3 : 1;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowStorageCreate().
 * @see PicurvWindowStorageCreate()
 */
PetscErrorCode PicurvWindowStorageCreate(UserCtx *user, const PicurvWindowDefinition *definition,
                                         PicurvWindowStorage *storage)
{
    PetscInt m2_total = 0;
    PetscInt cm_total = 0;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context, definition, and storage are required.");
    PetscCall(PetscMemzero(storage, sizeof(*storage)));
    storage->field_count = definition->field_count;
    storage->covariance_count = definition->covariance_count;

    /* Occupancy is per point and scalar regardless of what is accumulated. */
    PetscCall(VecDuplicate(user->P, &storage->count));
    PetscCall(VecSet(storage->count, 0.0));
    PetscCall(VecDuplicate(user->P, &storage->weight));
    PetscCall(VecSet(storage->weight, 0.0));
    PetscCall(VecDuplicate(user->P, &storage->weight_sq));
    PetscCall(VecSet(storage->weight_sq, 0.0));

    if (definition->field_count > 0) {
        PetscCall(PetscCalloc1((size_t)definition->field_count, &storage->mean));
        PetscCall(PetscCalloc1((size_t)definition->field_count + 1, &storage->m2_offset));
    }
    for (PetscInt f = 0; f < definition->field_count; ++f) {
        FieldView view;
        PetscInt  components = 0;

        PetscCall(FieldGetView(user, (FieldId)definition->fields[f].field_id, &view));
        /* The mean matches the source field's own layout, so it is duplicated
         * directly from it and inherits DM and decomposition. */
        PetscCall(VecDuplicate(view.global_vec, &storage->mean[f]));
        PetscCall(VecSet(storage->mean[f], 0.0));
        storage->m2_offset[f] = m2_total;
        if (definition->fields[f].want_second) {
            PetscCall(PicurvProductComponentCount(view.descriptor->dof, &components));
            m2_total += components;
        }
    }
    if (definition->field_count > 0) storage->m2_offset[definition->field_count] = m2_total;

    if (m2_total > 0) {
        PetscCall(PetscCalloc1((size_t)m2_total, &storage->m2));
        for (PetscInt c = 0; c < m2_total; ++c) {
            PetscCall(VecDuplicate(user->P, &storage->m2[c]));
            PetscCall(VecSet(storage->m2[c], 0.0));
        }
    }

    if (definition->covariance_count > 0) {
        PetscCall(PetscCalloc1((size_t)definition->covariance_count + 1, &storage->cm_offset));
    }
    for (PetscInt p = 0; p < definition->covariance_count; ++p) {
        FieldView view_a, view_b;
        PetscInt  components = 0;

        PetscCall(FieldGetView(user, (FieldId)definition->covariances[p].first, &view_a));
        PetscCall(FieldGetView(user, (FieldId)definition->covariances[p].second, &view_b));
        PetscCall(PicurvCovarianceComponentCount(view_a.descriptor->dof, view_b.descriptor->dof, &components));
        storage->cm_offset[p] = cm_total;
        cm_total += components;
    }
    if (definition->covariance_count > 0) storage->cm_offset[definition->covariance_count] = cm_total;

    if (cm_total > 0) {
        PetscCall(PetscCalloc1((size_t)cm_total, &storage->cm));
        for (PetscInt c = 0; c < cm_total; ++c) {
            PetscCall(VecDuplicate(user->P, &storage->cm[c]));
            PetscCall(VecSet(storage->cm[c], 0.0));
        }
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowStorageDestroy().
 * @see PicurvWindowStorageDestroy()
 */
PetscErrorCode PicurvWindowStorageDestroy(PicurvWindowStorage *storage)
{
    PetscFunctionBeginUser;
    if (storage == NULL) PetscFunctionReturn(0);
    if (storage->count) PetscCall(VecDestroy(&storage->count));
    if (storage->weight) PetscCall(VecDestroy(&storage->weight));
    if (storage->weight_sq) PetscCall(VecDestroy(&storage->weight_sq));
    if (storage->mean) {
        for (PetscInt f = 0; f < storage->field_count; ++f) {
            if (storage->mean[f]) PetscCall(VecDestroy(&storage->mean[f]));
        }
        PetscCall(PetscFree(storage->mean));
    }
    if (storage->m2) {
        const PetscInt total = storage->m2_offset ? storage->m2_offset[storage->field_count] : 0;
        for (PetscInt c = 0; c < total; ++c) {
            if (storage->m2[c]) PetscCall(VecDestroy(&storage->m2[c]));
        }
        PetscCall(PetscFree(storage->m2));
    }
    if (storage->cm) {
        const PetscInt total = storage->cm_offset ? storage->cm_offset[storage->covariance_count] : 0;
        for (PetscInt c = 0; c < total; ++c) {
            if (storage->cm[c]) PetscCall(VecDestroy(&storage->cm[c]));
        }
        PetscCall(PetscFree(storage->cm));
    }
    if (storage->m2_offset) PetscCall(PetscFree(storage->m2_offset));
    if (storage->cm_offset) PetscCall(PetscFree(storage->cm_offset));
    PetscCall(PetscMemzero(storage, sizeof(*storage)));
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowAccumulate().
 * @see PicurvWindowAccumulate()
 */
PetscErrorCode PicurvWindowAccumulate(UserCtx *user, const PicurvWindowDefinition *definition,
                                      PicurvWindowStorage *storage, PetscReal weight)
{
    SpatialTargetPlan plan;
    PetscReal ***nvert = NULL;
    PetscReal ***count_arr = NULL, ***weight_arr = NULL, ***weight_sq_arr = NULL;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context, definition, and storage are required.");
    PetscCheck(weight > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Accepted states carry a positive weight, got %g.", (double)weight);
    if (definition->field_count == 0) PetscFunctionReturn(0);

    /* Every requested field shares the pointwise cell-centered domain in Phase 2,
     * so the plan is resolved once rather than per field. */
    PetscCall(SpatialTargetPlanCreate(user, (FieldId)definition->fields[0].field_id,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));

    PetscCall(DMDAVecGetArrayRead(user->da, user->Nvert, &nvert));
    PetscCall(DMDAVecGetArray(user->da, storage->count, &count_arr));
    PetscCall(DMDAVecGetArray(user->da, storage->weight, &weight_arr));
    PetscCall(DMDAVecGetArray(user->da, storage->weight_sq, &weight_sq_arr));

    /* Pass one advances per-point occupancy, so every field pass afterwards can
     * recover the pre-state weight uniformly as (stored weight - this weight).
     * Occupancy is per point and per accepted state, never per field. */
    for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
        for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
            for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) continue;
                count_arr[k][j][i]     += 1.0;
                weight_arr[k][j][i]    += weight;
                weight_sq_arr[k][j][i] += weight * weight;
            }
        }
    }

    /* Pass two updates each requested field's mean and, when asked, its centered
     * self-product. A three-vector's self-product is the six symmetric co-moments
     * between component pairs, not three per-component variances, so the co-moment
     * kernel drives every product component including the diagonal. */
    for (PetscInt f = 0; f < definition->field_count; ++f) {
        FieldView view;
        PetscReal ***src_scalar = NULL, ***mean_scalar = NULL;
        Cmpnts    ***src_vector = NULL, ***mean_vector = NULL;
        PetscReal ***product[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
        PetscInt  dof = 0, components = 0;
        const PetscBool second = definition->fields[f].want_second;

        PetscCall(FieldGetView(user, (FieldId)definition->fields[f].field_id, &view));
        dof = view.descriptor->dof;
        if (dof == 1) {
            PetscCall(DMDAVecGetArrayRead(view.dm, view.global_vec, &src_scalar));
            PetscCall(DMDAVecGetArray(view.dm, storage->mean[f], &mean_scalar));
        } else {
            PetscCall(DMDAVecGetArrayRead(view.dm, view.global_vec, &src_vector));
            PetscCall(DMDAVecGetArray(view.dm, storage->mean[f], &mean_vector));
        }
        if (second) {
            PetscCall(PicurvProductComponentCount(dof, &components));
            for (PetscInt c = 0; c < components; ++c) {
                PetscCall(DMDAVecGetArray(user->da, storage->m2[storage->m2_offset[f] + c], &product[c]));
            }
        }

        for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
            for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
                for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                    PetscReal value[3] = {0.0, 0.0, 0.0};
                    PetscReal mean_old[3] = {0.0, 0.0, 0.0};
                    PetscReal prior_weight = 0.0;
                    PetscReal prior_weight_sq = 0.0;
                    PetscReal prior_count = 0.0;

                    if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) continue;

                    prior_weight    = weight_arr[k][j][i] - weight;
                    prior_weight_sq = weight_sq_arr[k][j][i] - weight * weight;
                    prior_count     = count_arr[k][j][i] - 1.0;

                    if (dof == 1) {
                        value[0] = src_scalar[k][j][i];
                        mean_old[0] = mean_scalar[k][j][i];
                    } else {
                        value[0] = src_vector[k][j][i].x;
                        value[1] = src_vector[k][j][i].y;
                        value[2] = src_vector[k][j][i].z;
                        mean_old[0] = mean_vector[k][j][i].x;
                        mean_old[1] = mean_vector[k][j][i].y;
                        mean_old[2] = mean_vector[k][j][i].z;
                    }

                    /* Products are updated before the means are written back,
                     * because the co-moment kernel needs the pre-update means. */
                    for (PetscInt c = 0; c < components; ++c) {
                        const PetscInt a = (dof == 1) ? 0 : kProductFirst[c];
                        const PetscInt b = (dof == 1) ? 0 : kProductSecond[c];
                        PicurvCoMomentState pair;

                        pair.count = prior_count;
                        pair.weight = prior_weight;
                        pair.weight_sq = prior_weight_sq;
                        pair.mean_x = mean_old[a];
                        pair.mean_y = mean_old[b];
                        pair.cm = product[c][k][j][i];
                        PetscCall(PicurvCoMomentStateUpdate(&pair, value[a], value[b], weight));
                        product[c][k][j][i] = pair.cm;
                    }

                    for (PetscInt c = 0; c < dof; ++c) {
                        PicurvMomentState moment;

                        moment.count = prior_count;
                        moment.weight = prior_weight;
                        moment.weight_sq = prior_weight_sq;
                        moment.mean = mean_old[c];
                        moment.m2 = 0.0;
                        PetscCall(PicurvMomentStateUpdate(&moment, value[c], weight));
                        if (dof == 1) mean_scalar[k][j][i] = moment.mean;
                        else if (c == 0) mean_vector[k][j][i].x = moment.mean;
                        else if (c == 1) mean_vector[k][j][i].y = moment.mean;
                        else mean_vector[k][j][i].z = moment.mean;
                    }
                }
            }
        }

        if (second) {
            for (PetscInt c = 0; c < components; ++c) {
                PetscCall(DMDAVecRestoreArray(user->da, storage->m2[storage->m2_offset[f] + c], &product[c]));
            }
        }
        if (dof == 1) {
            PetscCall(DMDAVecRestoreArrayRead(view.dm, view.global_vec, &src_scalar));
            PetscCall(DMDAVecRestoreArray(view.dm, storage->mean[f], &mean_scalar));
        } else {
            PetscCall(DMDAVecRestoreArrayRead(view.dm, view.global_vec, &src_vector));
            PetscCall(DMDAVecRestoreArray(view.dm, storage->mean[f], &mean_vector));
        }
    }

    PetscCall(DMDAVecRestoreArray(user->da, storage->weight_sq, &weight_sq_arr));
    PetscCall(DMDAVecRestoreArray(user->da, storage->weight, &weight_arr));
    PetscCall(DMDAVecRestoreArray(user->da, storage->count, &count_arr));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Nvert, &nvert));
    PetscFunctionReturn(0);
}
