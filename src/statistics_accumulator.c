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
#include "io.h"
#include "logging.h"

/** @brief Longest output list a recipe may request. */
#define STATISTICS_DERIVED_OUTPUT_LENGTH 256

/**
 * @brief Upper-triangular row-major component pairs for a three-vector self-product.
 *
 * This pair table is the single definition of the symmetric component order. The
 * diagonal positions and the two-axis component labels are both derived from it
 * rather than restated, because accumulation reads the table while derivation reads
 * what follows from it: a hand-written copy that drifted would mislabel every
 * Reynolds stress and take the root of the wrong component, with nothing failing.
 */
static const PetscInt kProductFirst[6]  = {0, 0, 0, 1, 1, 2};
static const PetscInt kProductSecond[6] = {0, 1, 2, 1, 2, 2};

/** @brief Axis labels indexing the pair table above. */
static const char *const kAxisName[3] = {"x", "y", "z"};

/**
 * @brief Internal helper: the product index carrying one component's own variance.
 * @details Local to this translation unit. Found by searching the pair table for the
 *          entry pairing a component with itself, so the diagonal cannot be stated
 *          separately from the order it belongs to.
 */
static PetscInt ProductDiagonalIndex(PetscInt component)
{
    for (PetscInt c = 0; c < 6; ++c) {
        if (kProductFirst[c] == component && kProductSecond[c] == component) return c;
    }
    return -1;
}

/**
 * @brief Internal helper: writes the two-axis label of one product component.
 * @details Local to this translation unit. Built from the pair table, so "xy" and the
 *          slot it names can never refer to different pairs.
 */
static PetscErrorCode ProductComponentLabel(PetscInt index, char *out, size_t size)
{
    PetscFunctionBeginUser;
    PetscCheck(index >= 0 && index < 6, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Product component %" PetscInt_FMT " is outside the symmetric set.", index);
    PetscCall(PetscSNPrintf(out, size, "%s%s", kAxisName[kProductFirst[index]],
                            kAxisName[kProductSecond[index]]));
    PetscFunctionReturn(0);
}

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
    /* Scalar-scalar and vector-scalar pairs only; a vector-vector cross product
     * would need a nine-component carrier that nothing allocates. */
    PetscCheck((dof_a == 1 && dof_b == 1) || (dof_a == 3 && dof_b == 1) || (dof_a == 1 && dof_b == 3),
               PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Covariance is supported for scalar-scalar and vector-scalar pairs, got dof %" PetscInt_FMT
               " and %" PetscInt_FMT ".", dof_a, dof_b);
    *count = (dof_a == 3 || dof_b == 3) ? 3 : 1;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvStatisticsComponentDM().
 * @see PicurvStatisticsComponentDM()
 */
PetscErrorCode PicurvStatisticsComponentDM(UserCtx *user, PetscInt components, DM *dm)
{
    PetscFunctionBeginUser;
    PetscCheck(user != NULL && dm != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context and DM output are required.");
    switch (components) {
    case 1: *dm = user->da;   break;
    case 3: *dm = user->fda;  break;
    case 6: *dm = user->fda6; break;
    default:
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "No block DM carries %" PetscInt_FMT " accumulator components.", components);
    }
    PetscCheck(*dm != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "The block DM for %" PetscInt_FMT " accumulator components was never created.",
               components);
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowStorageCreate().
 * @see PicurvWindowStorageCreate()
 */
PetscErrorCode PicurvWindowStorageCreate(UserCtx *user, const PicurvWindowDefinition *definition,
                                         PicurvWindowStorage *storage)
{
    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context, definition, and storage are required.");
    PetscCall(PetscMemzero(storage, sizeof(*storage)));
    storage->field_count = definition->field_count;
    storage->covariance_count = definition->covariance_count;

    /* Occupancy is per point and scalar regardless of what is accumulated. */
    PetscCall(DMCreateGlobalVector(user->da, &storage->count));
    PetscCall(VecSet(storage->count, 0.0));
    PetscCall(DMCreateGlobalVector(user->da, &storage->weight));
    PetscCall(VecSet(storage->weight, 0.0));
    PetscCall(DMCreateGlobalVector(user->da, &storage->weight_sq));
    PetscCall(VecSet(storage->weight_sq, 0.0));

    if (definition->field_count > 0) {
        PetscCall(PetscCalloc1((size_t)definition->field_count, &storage->mean));
        PetscCall(PetscCalloc1((size_t)definition->field_count, &storage->m2));
    }
    for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
        FieldView view;
        PetscInt  components = 0;
        DM        product_dm = NULL;

        PetscCall(FieldGetView(user, (FieldId)definition->fields[field_index].field_id, &view));
        /* The mean matches the source field's own layout, so it comes from that
         * field's DM and inherits its decomposition. */
        PetscCall(DMCreateGlobalVector(view.dm, &storage->mean[field_index]));
        PetscCall(VecSet(storage->mean[field_index], 0.0));
        if (!definition->fields[field_index].want_second) continue;
        PetscCall(PicurvProductComponentCount(view.descriptor->dof, &components));
        PetscCall(PicurvStatisticsComponentDM(user, components, &product_dm));
        PetscCall(DMCreateGlobalVector(product_dm, &storage->m2[field_index]));
        PetscCall(VecSet(storage->m2[field_index], 0.0));
    }

    if (definition->covariance_count > 0) {
        PetscCall(PetscCalloc1((size_t)definition->covariance_count, &storage->cm));
    }
    for (PetscInt pair_index = 0; pair_index < definition->covariance_count; ++pair_index) {
        FieldView view_a, view_b;
        PetscInt  components = 0;
        DM        pair_dm = NULL;

        PetscCall(FieldGetView(user, (FieldId)definition->covariances[pair_index].first, &view_a));
        PetscCall(FieldGetView(user, (FieldId)definition->covariances[pair_index].second, &view_b));
        PetscCall(PicurvCovarianceComponentCount(view_a.descriptor->dof, view_b.descriptor->dof, &components));
        PetscCall(PicurvStatisticsComponentDM(user, components, &pair_dm));
        PetscCall(DMCreateGlobalVector(pair_dm, &storage->cm[pair_index]));
        PetscCall(VecSet(storage->cm[pair_index], 0.0));
    }

    /* Report what this window costs and what it covers, once per block at setup.
     * The point count is a function of layout, dimensions, and periodicity alone, so
     * it tells an operator what "per point" will mean before any sample is taken. */
    {
        SpatialTargetPlan plan;
        PetscInt payloads = 0;
        PetscInt points = 0;

        PetscCall(PicurvWindowStoragePayloadCount(storage, &payloads));
        PetscCall(SpatialTargetPlanCreate(user, (FieldId)definition->fields[0].field_id,
                                          PICURV_STATISTICS_MASK_FLUID, &plan));
        PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &points));
        LOG_ALLOW(LOCAL, LOG_DEBUG,
                  "Statistics window '%s' block %d: %d accumulator vector(s) over %d point(s).\n",
                  definition->name, (int)user->_this, (int)payloads, (int)points);
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
    for (PetscInt field_index = 0; storage->mean && field_index < storage->field_count; ++field_index) {
        if (storage->mean[field_index]) PetscCall(VecDestroy(&storage->mean[field_index]));
    }
    for (PetscInt field_index = 0; storage->m2 && field_index < storage->field_count; ++field_index) {
        if (storage->m2[field_index]) PetscCall(VecDestroy(&storage->m2[field_index]));
    }
    for (PetscInt pair_index = 0; storage->cm && pair_index < storage->covariance_count; ++pair_index) {
        if (storage->cm[pair_index]) PetscCall(VecDestroy(&storage->cm[pair_index]));
    }
    if (storage->mean) PetscCall(PetscFree(storage->mean));
    if (storage->m2) PetscCall(PetscFree(storage->m2));
    if (storage->cm) PetscCall(PetscFree(storage->cm));
    PetscCall(PetscMemzero(storage, sizeof(*storage)));
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: locates the storage slot holding one field's running mean.
 * @details Local to this translation unit. A covariance member must also appear in
 *          the window's field list, because the co-moment update needs that field's
 *          running mean; this is where that requirement is enforced.
 */
static PetscErrorCode FindFieldSlot(const PicurvWindowDefinition *definition, PetscInt field_id,
                                    PetscInt *slot)
{
    PetscFunctionBeginUser;
    for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
        if (definition->fields[field_index].field_id == field_id) {
            *slot = field_index;
            PetscFunctionReturn(0);
        }
    }
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
            "Window '%s' requests a covariance over field '%s', which is not in its field list.",
            definition->name, FieldCanonicalName((FieldId)field_id));
}

/**
 * @brief Implementation of \ref PicurvWindowStoragePayloadCount().
 * @see PicurvWindowStoragePayloadCount()
 */
PetscErrorCode PicurvWindowStoragePayloadCount(const PicurvWindowStorage *storage, PetscInt *count)
{
    PetscFunctionBeginUser;
    PetscCheck(storage != NULL && count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Storage and count output are required.");
    *count = 3 + storage->field_count + storage->covariance_count;
    /* A field without a requested second moment holds no product vector, so the
     * enumeration skips it rather than emitting an empty payload. */
    for (PetscInt field_index = 0; storage->m2 && field_index < storage->field_count; ++field_index) {
        if (storage->m2[field_index]) *count += 1;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowStoragePayload().
 * @see PicurvWindowStoragePayload()
 */
PetscErrorCode PicurvWindowStoragePayload(UserCtx *user, const PicurvWindowDefinition *definition,
                                          const PicurvWindowStorage *storage, PetscInt index,
                                          PicurvStatisticsPayload *payload)
{
    PetscInt total = 0;
    PetscInt cursor = index;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL && payload != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context, definition, storage, and payload output are required.");
    PetscCall(PicurvWindowStoragePayloadCount(storage, &total));
    PetscCheck(index >= 0 && index < total, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Payload index %" PetscInt_FMT " is outside [0, %" PetscInt_FMT ").", index, total);
    PetscCall(PetscMemzero(payload, sizeof(*payload)));

    /* Occupancy first, then means, then products, then co-moments. The order is the
     * persistence contract: the manifest inventory, the writer, and the reader all
     * walk it identically. */
    if (cursor < 3) {
        static const char *const occupancy_name[3] = {"count", "weight", "weight_sq"};
        Vec occupancy[3];

        occupancy[0] = storage->count;
        occupancy[1] = storage->weight;
        occupancy[2] = storage->weight_sq;
        PetscCall(PetscStrncpy(payload->name, occupancy_name[cursor], sizeof(payload->name)));
        payload->vec = occupancy[cursor];
        payload->components = 1;
        payload->role = "occupancy";
        payload->layout = FieldLayoutName(FIELD_LAYOUT_CELL_CENTERED);
        PetscFunctionReturn(0);
    }
    cursor -= 3;

    if (cursor < storage->field_count) {
        const FieldDescriptor *descriptor = NULL;

        PetscCall(FieldGetDescriptor((FieldId)definition->fields[cursor].field_id, &descriptor));
        PetscCall(PetscSNPrintf(payload->name, sizeof(payload->name), "%s_mean",
                                descriptor->canonical_name));
        payload->vec = storage->mean[cursor];
        payload->components = descriptor->dof;
        payload->role = "mean";
        payload->layout = FieldLayoutName(descriptor->layout);
        PetscFunctionReturn(0);
    }
    cursor -= storage->field_count;

    {
        PetscInt product_count = 0;
        PetscInt seen = 0;

        for (PetscInt field_index = 0; storage->m2 && field_index < storage->field_count; ++field_index) {
            if (storage->m2[field_index]) ++product_count;
        }
        if (cursor < product_count) {
            for (PetscInt field_index = 0; field_index < storage->field_count; ++field_index) {
                const FieldDescriptor *descriptor = NULL;

                if (!storage->m2[field_index]) continue;
                if (seen++ != cursor) continue;
                PetscCall(FieldGetDescriptor((FieldId)definition->fields[field_index].field_id, &descriptor));
                PetscCall(PetscSNPrintf(payload->name, sizeof(payload->name), "%s_m2",
                                        descriptor->canonical_name));
                payload->vec = storage->m2[field_index];
                PetscCall(PicurvProductComponentCount(descriptor->dof, &payload->components));
                payload->role = "second_moment";
                payload->layout = FieldLayoutName(descriptor->layout);
                PetscFunctionReturn(0);
            }
        }
        cursor -= product_count;
    }

    {
        const FieldDescriptor *first = NULL;
        const FieldDescriptor *second = NULL;

        PetscCheck(cursor >= 0 && cursor < storage->covariance_count, PETSC_COMM_SELF, PETSC_ERR_PLIB,
                   "Payload index %" PetscInt_FMT " fell through the storage enumeration.", index);
        PetscCall(FieldGetDescriptor((FieldId)definition->covariances[cursor].first, &first));
        PetscCall(FieldGetDescriptor((FieldId)definition->covariances[cursor].second, &second));
        PetscCall(PetscSNPrintf(payload->name, sizeof(payload->name), "%s_%s_cm",
                                first->canonical_name, second->canonical_name));
        payload->vec = storage->cm[cursor];
        PetscCall(PicurvCovarianceComponentCount(first->dof, second->dof, &payload->components));
        payload->role = "co_moment";
        payload->layout = FieldLayoutName(first->layout);
    }
    PetscFunctionReturn(0);
}

/** @brief Output kinds a postprocessing recipe may request, in enumeration order. */
typedef enum {
    DERIVED_MEAN = 0,
    DERIVED_REYNOLDS_STRESS,
    DERIVED_RMS,
    DERIVED_TKE,
    DERIVED_FLUX,
    DERIVED_KIND_COUNT
} DerivedKind;

/** @brief Recipe spellings of the output kinds. */
static const char *const kDerivedKindName[DERIVED_KIND_COUNT] = {
    "mean", "reynolds_stress", "rms", "tke", "flux"
};

/**
 * @brief Internal helper: reports which output kinds a recipe requested.
 * @details Local to this translation unit.
 */
static PetscErrorCode ParseDerivedKinds(const char *outputs, PetscBool wanted[DERIVED_KIND_COUNT])
{
    char buffer[STATISTICS_DERIVED_OUTPUT_LENGTH];
    char *cursor = buffer;

    PetscFunctionBeginUser;
    for (PetscInt k = 0; k < DERIVED_KIND_COUNT; ++k) wanted[k] = PETSC_FALSE;
    PetscCheck(outputs != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output list is required.");
    PetscCall(PetscStrncpy(buffer, outputs, sizeof(buffer)));

    while (cursor && *cursor) {
        char *comma = strchr(cursor, ',');
        PetscBool matched = PETSC_FALSE;

        if (comma) *comma = '\0';
        TrimWhitespace(cursor);
        if (cursor[0] != '\0') {
            for (PetscInt k = 0; k < DERIVED_KIND_COUNT; ++k) {
                if (strcmp(cursor, kDerivedKindName[k])) continue;
                wanted[k] = PETSC_TRUE;
                matched = PETSC_TRUE;
                break;
            }
            PetscCheck(matched, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                       "Unknown field-statistics output '%s'. Available outputs are "
                       "mean, reynolds_stress, rms, tke, and flux.", cursor);
        }
        cursor = comma ? comma + 1 : NULL;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: counts the derived fields each kind contributes.
 * @details Local to this translation unit. A kind contributes nothing when the state
 *          it needs was never accumulated, so a recipe may name every output without
 *          having to know which window carries which moment.
 */
static PetscErrorCode DerivedKindExtent(const PicurvWindowDefinition *definition,
                                        const PicurvWindowStorage *storage,
                                        DerivedKind kind, PetscInt *count)
{
    PetscFunctionBeginUser;
    *count = 0;
    switch (kind) {
    case DERIVED_MEAN:
        *count = storage->field_count;
        break;
    case DERIVED_REYNOLDS_STRESS:
    case DERIVED_RMS:
        for (PetscInt field_index = 0; storage->m2 && field_index < storage->field_count; ++field_index) {
            const FieldDescriptor *descriptor = NULL;
            PetscInt components = 0;

            if (!storage->m2[field_index]) continue;
            PetscCall(FieldGetDescriptor((FieldId)definition->fields[field_index].field_id, &descriptor));
            PetscCall(PicurvProductComponentCount(descriptor->dof, &components));
            /* A stress tensor emits every component; an RMS emits only the
             * diagonal, because an off-diagonal has no square root to take. */
            *count += (kind == DERIVED_REYNOLDS_STRESS) ? components : descriptor->dof;
        }
        break;
    case DERIVED_TKE:
        /* Turbulent kinetic energy is the trace of a three-vector's stress tensor,
         * so it exists only for a vector field carrying a second moment. */
        for (PetscInt field_index = 0; storage->m2 && field_index < storage->field_count; ++field_index) {
            const FieldDescriptor *descriptor = NULL;

            if (!storage->m2[field_index]) continue;
            PetscCall(FieldGetDescriptor((FieldId)definition->fields[field_index].field_id, &descriptor));
            if (descriptor->dof == 3) *count += 1;
        }
        break;
    default:
        *count = storage->covariance_count;
        break;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowDerivedCount().
 * @see PicurvWindowDerivedCount()
 */
PetscErrorCode PicurvWindowDerivedCount(const PicurvWindowDefinition *definition,
                                        const PicurvWindowStorage *storage,
                                        const char *outputs, PetscInt *count)
{
    PetscBool wanted[DERIVED_KIND_COUNT];

    PetscFunctionBeginUser;
    PetscCheck(definition != NULL && storage != NULL && count != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Definition, storage, and count are required.");
    PetscCall(ParseDerivedKinds(outputs, wanted));
    *count = 0;
    for (PetscInt k = 0; k < DERIVED_KIND_COUNT; ++k) {
        PetscInt extent = 0;

        if (!wanted[k]) continue;
        PetscCall(DerivedKindExtent(definition, storage, (DerivedKind)k, &extent));
        *count += extent;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: takes a square root of a variance that may be barely negative.
 * @details Local to this translation unit. Centered accumulation can leave a variance
 *          a few ulps below zero when the signal is nearly constant. Clamping is
 *          confined to this one place, applies only under a root, and never touches
 *          stored state; a genuinely negative variance is a defect and is reported.
 */
static PetscErrorCode SafeStandardDeviation(PetscReal variance, const char *label, PetscReal *result)
{
    PetscFunctionBeginUser;
    if (variance >= 0.0) {
        *result = PetscSqrtReal(variance);
        PetscFunctionReturn(0);
    }
    PetscCheck(variance >= -PICURV_STATISTICS_VARIANCE_FLOOR, PETSC_COMM_SELF, PETSC_ERR_FP,
               "Derived variance for '%s' is %g, which is too negative to be floating-point "
               "cancellation; the accumulated state is inconsistent.", label, (double)variance);
    *result = 0.0;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: resolves which kind and member one derived index selects.
 * @details Local to this translation unit. Walks the same order `DerivedKindExtent`
 *          counts, so the enumeration and the count cannot disagree.
 */
static PetscErrorCode ResolveDerivedIndex(const PicurvWindowDefinition *definition,
                                          const PicurvWindowStorage *storage,
                                          const PetscBool wanted[DERIVED_KIND_COUNT],
                                          PetscInt index, DerivedKind *kind, PetscInt *offset)
{
    PetscFunctionBeginUser;
    for (PetscInt k = 0; k < DERIVED_KIND_COUNT; ++k) {
        PetscInt extent = 0;

        if (!wanted[k]) continue;
        PetscCall(DerivedKindExtent(definition, storage, (DerivedKind)k, &extent));
        if (index < extent) {
            *kind = (DerivedKind)k;
            *offset = index;
            PetscFunctionReturn(0);
        }
        index -= extent;
    }
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
            "Derived index is past the end of the requested output set.");
}

/**
 * @brief Implementation of \ref PicurvWindowDerive().
 * @see PicurvWindowDerive()
 */
PetscErrorCode PicurvWindowDerive(UserCtx *user, const PicurvWindowDefinition *definition,
                                  const PicurvWindowStorage *storage, const char *outputs,
                                  PetscInt index, Vec scalar_target, Vec vector_target,
                                  PicurvDerivedField *field)
{
    PetscBool wanted[DERIVED_KIND_COUNT];
    DerivedKind kind = DERIVED_MEAN;
    SpatialTargetPlan plan;
    PetscReal ***weight_arr = NULL, ***weight_sq_arr = NULL, ***count_arr = NULL;
    PetscScalar ****source = NULL, ****target = NULL;
    const FieldDescriptor *descriptor = NULL;
    DM source_dm = NULL;
    Vec source_vec = NULL;
    PetscInt offset = 0, slot = 0, member = 0;
    /* The source's component count is not the output's: a scalar RMS reads a
     * six-component tensor, so each needs its own DM. */
    PetscInt source_components = 0;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL && field != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context, definition, storage, and output are required.");
    PetscCheck(index >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Derived index %" PetscInt_FMT " is negative.", index);
    PetscCall(PetscMemzero(field, sizeof(*field)));
    /* The output list is parsed once and the extents walked once. Resolving the
     * index is what bounds-checks it, so counting first would repeat both. */
    PetscCall(ParseDerivedKinds(outputs, wanted));
    PetscCall(ResolveDerivedIndex(definition, storage, wanted, index, &kind, &offset));

    /* Locate the accumulator slot and component this index selects. */
    switch (kind) {
    case DERIVED_MEAN:
        slot = offset;
        PetscCall(FieldGetDescriptor((FieldId)definition->fields[slot].field_id, &descriptor));
        source_vec = storage->mean[slot];
        field->components = descriptor->dof;
        source_components = descriptor->dof;
        PetscCall(PetscSNPrintf(field->name, sizeof(field->name), "%s_%s_mean",
                                definition->name, descriptor->canonical_name));
        break;
    case DERIVED_TKE:
        for (slot = 0; slot < storage->field_count; ++slot) {
            if (!storage->m2[slot]) continue;
            PetscCall(FieldGetDescriptor((FieldId)definition->fields[slot].field_id, &descriptor));
            if (descriptor->dof != 3) continue;
            if (offset-- == 0) break;
        }
        source_vec = storage->m2[slot];
        field->components = 1;
        PetscCall(PicurvProductComponentCount(descriptor->dof, &source_components));
        PetscCall(PetscSNPrintf(field->name, sizeof(field->name), "%s_%s_tke",
                                definition->name, descriptor->canonical_name));
        break;
    case DERIVED_FLUX:
        slot = offset;
        {
            const FieldDescriptor *first = NULL;
            const FieldDescriptor *second = NULL;

            PetscCall(FieldGetDescriptor((FieldId)definition->covariances[slot].first, &first));
            PetscCall(FieldGetDescriptor((FieldId)definition->covariances[slot].second, &second));
            PetscCall(PicurvCovarianceComponentCount(first->dof, second->dof, &field->components));
            source_components = field->components;
            source_vec = storage->cm[slot];
            PetscCall(PetscSNPrintf(field->name, sizeof(field->name), "%s_%s_%s_flux",
                                    definition->name, first->canonical_name, second->canonical_name));
        }
        break;
    default:
        /* Stress and RMS both walk the fields carrying a product; stress emits every
         * symmetric component, RMS only the diagonal. */
        for (slot = 0; slot < storage->field_count; ++slot) {
            PetscInt extent = 0;

            if (!storage->m2[slot]) continue;
            PetscCall(FieldGetDescriptor((FieldId)definition->fields[slot].field_id, &descriptor));
            PetscCall(PicurvProductComponentCount(descriptor->dof, &extent));
            if (kind == DERIVED_RMS) extent = descriptor->dof;
            if (offset < extent) { member = offset; break; }
            offset -= extent;
        }
        source_vec = storage->m2[slot];
        field->components = 1;
        PetscCall(PicurvProductComponentCount(descriptor->dof, &source_components));
        if (kind == DERIVED_RMS) {
            PetscCall(PetscSNPrintf(field->name, sizeof(field->name), "%s_%s_rms%s",
                                    definition->name, descriptor->canonical_name,
                                    descriptor->dof == 1 ? "" : kAxisName[member]));
        } else if (descriptor->dof == 1) {
            PetscCall(PetscSNPrintf(field->name, sizeof(field->name), "%s_%s_variance",
                                    definition->name, descriptor->canonical_name));
        } else {
            char label[8];

            PetscCall(ProductComponentLabel(member, label, sizeof(label)));
            PetscCall(PetscSNPrintf(field->name, sizeof(field->name), "%s_%s_R_%s",
                                    definition->name, descriptor->canonical_name, label));
        }
        break;
    }

    PetscCall(SpatialTargetPlanCreate(user, (FieldId)definition->fields[0].field_id,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvStatisticsComponentDM(user, source_components, &source_dm));
    {
        Vec destination = (field->components == 1) ? scalar_target : vector_target;
        DM destination_dm = NULL;

        PetscCheck(destination != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                   "Derived field '%s' needs a %d-component destination.",
                   field->name, (int)field->components);
        PetscCall(PicurvStatisticsComponentDM(user, field->components, &destination_dm));
        PetscCall(VecZeroEntries(destination));

        PetscCall(DMDAVecGetArrayRead(user->da, storage->weight, &weight_arr));
        PetscCall(DMDAVecGetArrayRead(user->da, storage->weight_sq, &weight_sq_arr));
        PetscCall(DMDAVecGetArrayRead(user->da, storage->count, &count_arr));
        PetscCall(DMDAVecGetArrayDOFRead(source_dm, source_vec, &source));
        PetscCall(DMDAVecGetArrayDOF(destination_dm, destination, &target));

        for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
            for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
                for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                    PicurvCoMomentState pair;

                    /* A point the window never sampled has no average at all, so it
                     * is left at zero rather than divided by a zero weight. */
                    if (weight_arr[k][j][i] <= 0.0) continue;
                    pair.count = count_arr[k][j][i];
                    pair.weight = weight_arr[k][j][i];
                    pair.weight_sq = weight_sq_arr[k][j][i];
                    pair.mean_x = 0.0;
                    pair.mean_y = 0.0;

                    if (kind == DERIVED_MEAN) {
                        for (PetscInt c = 0; c < field->components; ++c) {
                            target[k][j][i][c] = source[k][j][i][c];
                        }
                    } else if (kind == DERIVED_TKE) {
                        /* The trace of the stress tensor, halved: components 0, 3
                         * and 5 are xx, yy and zz in the stored symmetric order. */
                        PetscReal trace = 0.0;

                        for (PetscInt c = 0; c < 3; ++c) {
                            pair.cm = source[k][j][i][ProductDiagonalIndex(c)];
                            trace += PicurvCoMomentStateCovariance(&pair);
                        }
                        target[k][j][i][0] = 0.5 * trace;
                    } else if (kind == DERIVED_RMS) {
                        PetscReal deviation = 0.0;

                        pair.cm = source[k][j][i][(descriptor->dof == 1)
                                                  ? 0 : ProductDiagonalIndex(member)];
                        PetscCall(SafeStandardDeviation(PicurvCoMomentStateCovariance(&pair),
                                                        field->name, &deviation));
                        target[k][j][i][0] = deviation;
                    } else if (kind == DERIVED_FLUX) {
                        for (PetscInt c = 0; c < field->components; ++c) {
                            pair.cm = source[k][j][i][c];
                            target[k][j][i][c] = PicurvCoMomentStateCovariance(&pair);
                        }
                    } else {
                        pair.cm = source[k][j][i][member];
                        target[k][j][i][0] = PicurvCoMomentStateCovariance(&pair);
                    }
                }
            }
        }

        PetscCall(DMDAVecRestoreArrayDOF(destination_dm, destination, &target));
        PetscCall(DMDAVecRestoreArrayDOFRead(source_dm, source_vec, &source));
        PetscCall(DMDAVecRestoreArrayRead(user->da, storage->count, &count_arr));
        PetscCall(DMDAVecRestoreArrayRead(user->da, storage->weight_sq, &weight_sq_arr));
        PetscCall(DMDAVecRestoreArrayRead(user->da, storage->weight, &weight_arr));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowSpatialMean().
 * @see PicurvWindowSpatialMean()
 */
PetscErrorCode PicurvWindowSpatialMean(UserCtx *user, const PicurvWindowDefinition *definition,
                                       const PicurvWindowStorage *storage, Vec field,
                                       PetscReal *mean)
{
    SpatialTargetPlan plan;
    const PetscReal ***values = NULL;
    PetscReal ***weight_arr = NULL;
    PetscReal local_sum = 0.0;
    PetscReal local_count = 0.0;
    PetscReal totals[2] = {0.0, 0.0};
    PetscReal reduced[2] = {0.0, 0.0};

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL && field != NULL && mean != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context, definition, storage, field, and output are required.");
    *mean = 0.0;
    if (definition->field_count == 0) PetscFunctionReturn(0);

    PetscCall(SpatialTargetPlanCreate(user, (FieldId)definition->fields[0].field_id,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(DMDAVecGetArrayRead(user->da, field, &values));
    PetscCall(DMDAVecGetArray(user->da, storage->weight, &weight_arr));
    for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
        for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
            for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                /* A point with no accumulated weight holds a zero that means "never
                 * measured", so it is excluded from both the sum and the count. */
                if (weight_arr[k][j][i] <= 0.0) continue;
                local_sum += values[k][j][i];
                local_count += 1.0;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, storage->weight, &weight_arr));
    PetscCall(DMDAVecRestoreArrayRead(user->da, field, &values));

    totals[0] = local_sum;
    totals[1] = local_count;
    PetscCallMPI(MPI_Allreduce(totals, reduced, 2, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD));
    if (reduced[1] > 0.0) *mean = reduced[0] / reduced[1];
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowValidFractionRange().
 * @see PicurvWindowValidFractionRange()
 */
PetscErrorCode PicurvWindowValidFractionRange(UserCtx *user, const PicurvWindowDefinition *definition,
                                              const PicurvWindowStorage *storage, PetscInt sample_count,
                                              PetscReal *minimum, PetscReal *maximum)
{
    SpatialTargetPlan plan;
    PetscReal ***nvert = NULL;
    PetscReal ***count_arr = NULL;
    PetscReal local_min = PETSC_MAX_REAL;
    PetscReal local_max = 0.0;
    PetscReal reduced_min = 0.0, reduced_max = 0.0;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && definition != NULL && storage != NULL &&
               minimum != NULL && maximum != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context, definition, storage, and outputs are required.");
    *minimum = 1.0;
    *maximum = 0.0;
    if (definition->field_count == 0 || sample_count <= 0) PetscFunctionReturn(0);

    PetscCall(SpatialTargetPlanCreate(user, (FieldId)definition->fields[0].field_id,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Nvert, &nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, storage->count, &count_arr));
    for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
        for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
            for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                /* The mask is evaluated at the current state, but a point excluded
                 * now may have contributed earlier, so its stored count is what
                 * decides its fraction rather than its present eligibility. */
                const PetscReal fraction = count_arr[k][j][i] / (PetscReal)sample_count;

                local_min = PetscMin(local_min, fraction);
                local_max = PetscMax(local_max, fraction);
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, storage->count, &count_arr));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Nvert, &nvert));

    /* A rank owning no targeted point must not drag the minimum down, so its
     * sentinel is neutral under the reduction rather than zero. */
    PetscCallMPI(MPI_Allreduce(&local_min, &reduced_min, 1, MPIU_REAL, MPIU_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&local_max, &reduced_max, 1, MPIU_REAL, MPIU_MAX, PETSC_COMM_WORLD));
    *minimum = (reduced_min == PETSC_MAX_REAL) ? 1.0 : reduced_min;
    *maximum = reduced_max;
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

    /* Every requested field shares the pointwise cell-centered domain, so the plan
     * is resolved once rather than per field. */
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

    /* Pass two updates every cross-field co-moment. It runs before any mean is
     * written back, because a co-moment needs the pre-update mean of *both*
     * members, and pass three overwrites them field by field. */
    for (PetscInt pair = 0; pair < definition->covariance_count; ++pair) {
        FieldView view_a, view_b;
        PetscScalar ****src_a = NULL, ****mean_a = NULL;
        PetscScalar ****src_b = NULL, ****mean_b = NULL;
        PetscScalar ****co_moment = NULL;
        DM pair_dm = NULL;
        PetscInt slot_a = 0, slot_b = 0, dof_a = 0, dof_b = 0, components = 0;

        PetscCall(FindFieldSlot(definition, definition->covariances[pair].first, &slot_a));
        PetscCall(FindFieldSlot(definition, definition->covariances[pair].second, &slot_b));
        PetscCall(FieldGetView(user, (FieldId)definition->covariances[pair].first, &view_a));
        PetscCall(FieldGetView(user, (FieldId)definition->covariances[pair].second, &view_b));
        dof_a = view_a.descriptor->dof;
        dof_b = view_b.descriptor->dof;
        PetscCall(PicurvCovarianceComponentCount(dof_a, dof_b, &components));
        PetscCall(PicurvStatisticsComponentDM(user, components, &pair_dm));

        /* A component-indexed view serves every degree of freedom, so the loop
         * below needs no scalar-versus-vector branching. */
        PetscCall(DMDAVecGetArrayDOFRead(view_a.dm, view_a.global_vec, &src_a));
        PetscCall(DMDAVecGetArrayDOFRead(view_a.dm, storage->mean[slot_a], &mean_a));
        PetscCall(DMDAVecGetArrayDOFRead(view_b.dm, view_b.global_vec, &src_b));
        PetscCall(DMDAVecGetArrayDOFRead(view_b.dm, storage->mean[slot_b], &mean_b));
        PetscCall(DMDAVecGetArrayDOF(pair_dm, storage->cm[pair], &co_moment));

        for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
            for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
                for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                    if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) continue;
                    for (PetscInt c = 0; c < components; ++c) {
                        /* A scalar member contributes its single value against every
                         * component of a vector member, which is what makes a
                         * vector-scalar covariance a three-component object. */
                        const PetscInt component_a = (dof_a == 3) ? c : 0;
                        const PetscInt component_b = (dof_b == 3) ? c : 0;
                        PicurvCoMomentState state;

                        state.count = count_arr[k][j][i] - 1.0;
                        state.weight = weight_arr[k][j][i] - weight;
                        state.weight_sq = weight_sq_arr[k][j][i] - weight * weight;
                        state.mean_x = mean_a[k][j][i][component_a];
                        state.mean_y = mean_b[k][j][i][component_b];
                        state.cm = co_moment[k][j][i][c];
                        PetscCall(PicurvCoMomentStateUpdate(&state,
                                                            src_a[k][j][i][component_a],
                                                            src_b[k][j][i][component_b], weight));
                        co_moment[k][j][i][c] = state.cm;
                    }
                }
            }
        }

        PetscCall(DMDAVecRestoreArrayDOF(pair_dm, storage->cm[pair], &co_moment));
        PetscCall(DMDAVecRestoreArrayDOFRead(view_b.dm, storage->mean[slot_b], &mean_b));
        PetscCall(DMDAVecRestoreArrayDOFRead(view_b.dm, view_b.global_vec, &src_b));
        PetscCall(DMDAVecRestoreArrayDOFRead(view_a.dm, storage->mean[slot_a], &mean_a));
        PetscCall(DMDAVecRestoreArrayDOFRead(view_a.dm, view_a.global_vec, &src_a));
    }

    /* Pass three updates each requested field's mean and, when asked, its centered
     * self-product. A three-vector's self-product is the six symmetric co-moments
     * between component pairs, not three per-component variances, so the co-moment
     * kernel drives every product component including the diagonal. */
    for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
        FieldView view;
        PetscScalar ****src = NULL, ****mean = NULL, ****product = NULL;
        DM product_dm = NULL;
        PetscInt dof = 0, components = 0;
        const PetscBool second = definition->fields[field_index].want_second;

        PetscCall(FieldGetView(user, (FieldId)definition->fields[field_index].field_id, &view));
        dof = view.descriptor->dof;
        PetscCall(DMDAVecGetArrayDOFRead(view.dm, view.global_vec, &src));
        PetscCall(DMDAVecGetArrayDOF(view.dm, storage->mean[field_index], &mean));
        if (second) {
            PetscCall(PicurvProductComponentCount(dof, &components));
            PetscCall(PicurvStatisticsComponentDM(user, components, &product_dm));
            PetscCall(DMDAVecGetArrayDOF(product_dm, storage->m2[field_index], &product));
        }

        for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k) {
            for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j) {
                for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) {
                    PetscReal prior_weight = 0.0;
                    PetscReal prior_weight_sq = 0.0;
                    PetscReal prior_count = 0.0;

                    if (!SpatialTargetPlanMaskAllows(&plan, nvert[k][j][i])) continue;

                    prior_weight    = weight_arr[k][j][i] - weight;
                    prior_weight_sq = weight_sq_arr[k][j][i] - weight * weight;
                    prior_count     = count_arr[k][j][i] - 1.0;

                    /* Products are updated before the means are written back,
                     * because the co-moment kernel needs the pre-update means. */
                    for (PetscInt c = 0; c < components; ++c) {
                        const PetscInt a = (dof == 1) ? 0 : kProductFirst[c];
                        const PetscInt b = (dof == 1) ? 0 : kProductSecond[c];
                        PicurvCoMomentState state;

                        state.count = prior_count;
                        state.weight = prior_weight;
                        state.weight_sq = prior_weight_sq;
                        state.mean_x = mean[k][j][i][a];
                        state.mean_y = mean[k][j][i][b];
                        state.cm = product[k][j][i][c];
                        PetscCall(PicurvCoMomentStateUpdate(&state, src[k][j][i][a],
                                                            src[k][j][i][b], weight));
                        product[k][j][i][c] = state.cm;
                    }

                    for (PetscInt c = 0; c < dof; ++c) {
                        PicurvMomentState moment;

                        moment.count = prior_count;
                        moment.weight = prior_weight;
                        moment.weight_sq = prior_weight_sq;
                        moment.mean = mean[k][j][i][c];
                        moment.m2 = 0.0;
                        PetscCall(PicurvMomentStateUpdate(&moment, src[k][j][i][c], weight));
                        mean[k][j][i][c] = moment.mean;
                    }
                }
            }
        }

        if (second) PetscCall(DMDAVecRestoreArrayDOF(product_dm, storage->m2[field_index], &product));
        PetscCall(DMDAVecRestoreArrayDOF(view.dm, storage->mean[field_index], &mean));
        PetscCall(DMDAVecRestoreArrayDOFRead(view.dm, view.global_vec, &src));
    }

    PetscCall(DMDAVecRestoreArray(user->da, storage->weight_sq, &weight_sq_arr));
    PetscCall(DMDAVecRestoreArray(user->da, storage->weight, &weight_arr));
    PetscCall(DMDAVecRestoreArray(user->da, storage->count, &count_arr));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Nvert, &nvert));
    PetscFunctionReturn(0);
}
