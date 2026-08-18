/**
 * @file statistics_window.c
 * @brief Window lifecycle, scheduling, and weighting for the field-statistics pipeline.
 *
 * Full API contract is documented with the declarations in
 * `include/statistics_window.h`.
 */

#include "statistics_window.h"
#include "variables.h"
#include "logging.h"
#include "statistics_accumulator.h"
#include "field_catalog.h"
#include "checksum.h"

/** @brief Tolerance for treating two physical times as the same instant. */
#define WINDOW_TIME_EPSILON 1.0e-12

/**
 * @brief Implementation of \ref PicurvWindowStateName().
 * @see PicurvWindowStateName()
 */
const char *PicurvWindowStateName(PicurvWindowState state)
{
    switch (state) {
    case PICURV_WINDOW_PENDING:  return "pending";
    case PICURV_WINDOW_ACTIVE:   return "active";
    case PICURV_WINDOW_COMPLETE: return "complete";
    default:                     return "unknown";
    }
}

/**
 * @brief Implementation of \ref PicurvWindowInit().
 * @see PicurvWindowInit()
 */
PetscErrorCode PicurvWindowInit(PicurvWindow *window, const PicurvWindowDefinition *definition)
{
    PetscFunctionBeginUser;
    PetscCheck(window != NULL && definition != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Window and definition are required.");
    PetscCheck(definition->name[0] != '\0', PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Window name must not be empty.");
    if (definition->cadence_kind == PICURV_CADENCE_STEP) {
        PetscCheck(definition->step_cadence > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                   "Window '%s' needs a positive step cadence, got %" PetscInt_FMT ".",
                   definition->name, definition->step_cadence);
    } else {
        PetscCheck(definition->time_cadence > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                   "Window '%s' needs a positive time cadence, got %g.",
                   definition->name, (double)definition->time_cadence);
    }
    PetscCheck(!definition->bounded || definition->end_time > definition->start_time,
               PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Window '%s' must end after it starts.", definition->name);

    PetscCall(PetscMemzero(window, sizeof(*window)));
    window->definition = *definition;
    window->state = PICURV_WINDOW_PENDING;
    /* Anchoring the interval origin at the requested start, rather than at the
     * first accepted state, is what lets the first sample represent the interval
     * back to the requested bound instead of silently dropping it. */
    window->effective_start = definition->start_time;
    window->effective_end = definition->start_time;
    window->last_accepted_time = definition->start_time;
    window->activation_step = -1;
    window->last_event_step = PETSC_MIN_INT;
    window->next_time_target = 0;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: reports whether a state is due under the window's schedule.
 * @details Local to this translation unit. Advances no state.
 */
static PetscBool WindowStateIsDue(const PicurvWindow *window, PetscInt step, PetscReal time)
{
    if (window->definition.cadence_kind == PICURV_CADENCE_STEP) {
        /* The activation step is itself due; subsequent ones follow the stride. */
        if (window->activation_step < 0) return PETSC_TRUE;
        return (PetscBool)(((step - window->activation_step) % window->definition.step_cadence) == 0);
    }
    /* Time cadence: targets sit on the absolute grid effective_start + k*cadence, so
     * the schedule cannot drift as dt varies. A state is due once it reaches the
     * next outstanding target. */
    {
        const PetscReal target = window->effective_start +
                                 (PetscReal)window->next_time_target * window->definition.time_cadence;
        return (PetscBool)(time >= target - WINDOW_TIME_EPSILON);
    }
}

/**
 * @brief Internal helper: advances the time-cadence target past the accepted time.
 * @details Local to this translation unit. A single large step may overshoot several
 *          targets; they are consumed at once, because the accepted state's weight is
 *          the actual elapsed interval and so already accounts for the whole span.
 */
static void WindowAdvanceTimeTarget(PicurvWindow *window, PetscReal time)
{
    if (window->definition.cadence_kind != PICURV_CADENCE_TIME) return;
    while (window->effective_start +
           (PetscReal)window->next_time_target * window->definition.time_cadence
           <= time + WINDOW_TIME_EPSILON) {
        window->next_time_target += 1;
    }
}

/**
 * @brief Implementation of \ref PicurvWindowOfferState().
 * @see PicurvWindowOfferState()
 */
PetscErrorCode PicurvWindowOfferState(PicurvWindow *window, PetscInt step, PetscReal time,
                                      PetscBool *accepted, PetscReal *weight)
{
    PetscReal interval = 0.0;
    PetscReal right_edge = 0.0;
    PetscBool closes = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCheck(window != NULL && accepted != NULL && weight != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Window, acceptance flag, and weight output are required.");
    *accepted = PETSC_FALSE;
    *weight = 0.0;

    if (window->state == PICURV_WINDOW_COMPLETE) PetscFunctionReturn(0);
    /* A completed state is accepted at most once, so a step already seen is
     * rejected outright rather than double counted. */
    if (step == window->last_event_step) PetscFunctionReturn(0);
    if (time < window->definition.start_time - WINDOW_TIME_EPSILON) PetscFunctionReturn(0);

    /* First observation: a window resumed later than its requested start cannot
     * claim represented time it never saw, so the origin moves forward. */
    if (window->state == PICURV_WINDOW_PENDING) {
        if (time > window->definition.start_time + WINDOW_TIME_EPSILON) {
            window->effective_start = time;
            window->last_accepted_time = time;
        }
        window->state = PICURV_WINDOW_ACTIVE;
        window->activation_step = step;
        window->effective_end = window->effective_start;
        window->next_time_target = 0;
        WindowAdvanceTimeTarget(window, window->effective_start);
    }

    if (!WindowStateIsDue(window, step, time)) {
        /* An active off-schedule state changes no scientific state. */
        PetscFunctionReturn(0);
    }
    window->last_event_step = step;

    /* Clip the final represented interval to the requested bound. */
    right_edge = time;
    if (window->definition.bounded && time >= window->definition.end_time - WINDOW_TIME_EPSILON) {
        right_edge = window->definition.end_time;
        closes = PETSC_TRUE;
    }

    interval = right_edge - window->last_accepted_time;
    WindowAdvanceTimeTarget(window, time);

    if (interval <= WINDOW_TIME_EPSILON) {
        /* A zero-length interval is not a sample. It still anchors the origin, which
         * is what makes initial-state handling identical under both weightings. */
        window->last_accepted_time = right_edge;
        window->effective_end = right_edge;
        if (closes) window->state = PICURV_WINDOW_COMPLETE;
        PetscFunctionReturn(0);
    }

    *weight = (window->definition.weighting == PICURV_WEIGHTING_SAMPLE) ? 1.0 : interval;
    *accepted = PETSC_TRUE;

    window->last_accepted_time = right_edge;
    window->effective_end = right_edge;
    window->sample_count += 1;
    window->total_weight += *weight;
    window->represented_time += interval;
    if (closes) window->state = PICURV_WINDOW_COMPLETE;
    PetscFunctionReturn(0);
}

/** @brief Stable names of the hashed property groups, in serialization order. */
static const char *const kHashGroupNames[PICURV_WINDOW_HASH_GROUP_COUNT] = {
    "name", "start_time", "weighting", "cadence", "fields", "covariances", "mask", "target"
};

/**
 * @brief Implementation of \ref PicurvWindowHashGroupName().
 * @see PicurvWindowHashGroupName()
 */
const char *PicurvWindowHashGroupName(PetscInt group)
{
    if (group < 0 || group >= PICURV_WINDOW_HASH_GROUP_COUNT) return "unknown";
    return kHashGroupNames[group];
}

/**
 * @brief Internal helper: orders request indices so listing order cannot change the hash.
 * @details Local to this translation unit. Insertion sort over at most
 *          `PICURV_WINDOW_MAX_REQUESTS` entries, keyed on a primary and secondary
 *          integer so the same routine serves both fields and covariance pairs.
 */
static void SortRequestOrder(PetscInt count, const PetscInt *primary, const PetscInt *secondary,
                             PetscInt *order)
{
    for (PetscInt n = 0; n < count; ++n) order[n] = n;
    for (PetscInt n = 1; n < count; ++n) {
        const PetscInt candidate = order[n];
        PetscInt m = n - 1;

        while (m >= 0 &&
               (primary[order[m]] > primary[candidate] ||
                (primary[order[m]] == primary[candidate] &&
                 secondary[order[m]] > secondary[candidate]))) {
            order[m + 1] = order[m];
            --m;
        }
        order[m + 1] = candidate;
    }
}

/**
 * @brief Implementation of \ref PicurvWindowComputeHash().
 * @see PicurvWindowComputeHash()
 */
PetscErrorCode PicurvWindowComputeHash(const PicurvWindowDefinition *definition,
                                       char digest_hex[65],
                                       char group_digest_hex[][PICURV_WINDOW_HASH_GROUP_LENGTH])
{
    PicurvSHA256Context whole;
    PetscInt field_primary[PICURV_WINDOW_MAX_REQUESTS];
    PetscInt field_secondary[PICURV_WINDOW_MAX_REQUESTS];
    PetscInt field_order[PICURV_WINDOW_MAX_REQUESTS];
    PetscInt pair_primary[PICURV_WINDOW_MAX_REQUESTS];
    PetscInt pair_secondary[PICURV_WINDOW_MAX_REQUESTS];
    PetscInt pair_order[PICURV_WINDOW_MAX_REQUESTS];

    PetscFunctionBeginUser;
    PetscCheck(definition != NULL && digest_hex != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Definition and digest output are required.");
    PetscCheck(definition->field_count >= 0 && definition->field_count <= PICURV_WINDOW_MAX_REQUESTS,
               PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Window '%s' requests %" PetscInt_FMT " fields; at most %d are supported.",
               definition->name, definition->field_count, PICURV_WINDOW_MAX_REQUESTS);
    PetscCheck(definition->covariance_count >= 0 &&
               definition->covariance_count <= PICURV_WINDOW_MAX_REQUESTS,
               PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Window '%s' requests %" PetscInt_FMT " covariances; at most %d are supported.",
               definition->name, definition->covariance_count, PICURV_WINDOW_MAX_REQUESTS);

    for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
        field_primary[field_index] = definition->fields[field_index].field_id;
        field_secondary[field_index] = 0;
    }
    SortRequestOrder(definition->field_count, field_primary, field_secondary, field_order);
    for (PetscInt pair_index = 0; pair_index < definition->covariance_count; ++pair_index) {
        const PetscInt a = definition->covariances[pair_index].first;
        const PetscInt b = definition->covariances[pair_index].second;

        /* Covariance is symmetric, so the pair is canonicalized here rather than
         * requiring the user to list its members in catalog order. */
        pair_primary[pair_index] = PetscMin(a, b);
        pair_secondary[pair_index] = PetscMax(a, b);
    }
    SortRequestOrder(definition->covariance_count, pair_primary, pair_secondary, pair_order);

    PicurvSHA256Init(&whole);
    for (PetscInt group = 0; group < PICURV_WINDOW_HASH_GROUP_COUNT; ++group) {
        char text[1024];
        char group_digest[65];
        PicurvSHA256Context part;
        size_t used = 0;
        size_t length = 0;

        text[0] = '\0';
        switch (group) {
        case 0:
            PetscCall(PetscSNPrintf(text, sizeof(text), "name=%s\n", definition->name));
            break;
        case 1:
            PetscCall(PetscSNPrintf(text, sizeof(text), "start_time=%.17g\n",
                                    (double)definition->start_time));
            break;
        case 2:
            PetscCall(PetscSNPrintf(text, sizeof(text), "weighting=%s\n",
                                    definition->weighting == PICURV_WEIGHTING_SAMPLE
                                        ? "sample" : "physical_time"));
            break;
        case 3:
            if (definition->cadence_kind == PICURV_CADENCE_STEP) {
                PetscCall(PetscSNPrintf(text, sizeof(text), "cadence=step:%" PetscInt_FMT "\n",
                                        definition->step_cadence));
            } else {
                PetscCall(PetscSNPrintf(text, sizeof(text), "cadence=time:%.17g\n",
                                        (double)definition->time_cadence));
            }
            break;
        case 4:
            for (PetscInt n = 0; n < definition->field_count; ++n) {
                const PicurvWindowFieldRequest *request = &definition->fields[field_order[n]];
                const char *name = FieldCanonicalName((FieldId)request->field_id);

                PetscCheck(name != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                           "Window '%s' requests unknown field id %" PetscInt_FMT ".",
                           definition->name, request->field_id);
                PetscCall(PetscStrlen(text, &used));
                PetscCall(PetscSNPrintf(text + used, sizeof(text) - used, "field=%s:%s\n",
                                        name, request->want_second ? "first,second" : "first"));
            }
            break;
        case 5:
            for (PetscInt n = 0; n < definition->covariance_count; ++n) {
                const PetscInt index = pair_order[n];
                const char *first = FieldCanonicalName((FieldId)pair_primary[index]);
                const char *second = FieldCanonicalName((FieldId)pair_secondary[index]);

                PetscCheck(first != NULL && second != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                           "Window '%s' requests a covariance over an unknown field id.",
                           definition->name);
                PetscCall(PetscStrlen(text, &used));
                PetscCall(PetscSNPrintf(text + used, sizeof(text) - used, "covariance=%s,%s\n",
                                        first, second));
            }
            break;
        case 6:
            /* Exactly one mask is resolved, so this group hashes a constant. It
             * exists so that a future mask key extends this group instead of
             * renumbering every group after it and invalidating saved state. */
            PetscCall(PetscSNPrintf(text, sizeof(text), "mask=fluid\n"));
            break;
        default:
            PetscCall(PetscSNPrintf(text, sizeof(text), "target=pointwise\n"));
            for (PetscInt n = 0; n < definition->field_count; ++n) {
                const PicurvWindowFieldRequest *request = &definition->fields[field_order[n]];
                const FieldDescriptor *descriptor = NULL;

                PetscCall(FieldGetDescriptor((FieldId)request->field_id, &descriptor));
                PetscCall(PetscStrlen(text, &used));
                PetscCall(PetscSNPrintf(text + used, sizeof(text) - used, "layout=%s:%s\n",
                                        descriptor->canonical_name,
                                        FieldLayoutName(descriptor->layout)));
            }
            break;
        }

        PetscCall(PetscStrlen(text, &length));
        PicurvSHA256Update(&whole, text, length);
        if (group_digest_hex) {
            PicurvSHA256Init(&part);
            PicurvSHA256Update(&part, text, length);
            PicurvSHA256FinalHex(&part, group_digest);
            PetscCall(PetscStrncpy(group_digest_hex[group], group_digest,
                                   PICURV_WINDOW_HASH_GROUP_LENGTH));
        }
    }
    PicurvSHA256FinalHex(&whole, digest_hex);
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowFirstHashDifference().
 * @see PicurvWindowFirstHashDifference()
 */
PetscErrorCode PicurvWindowFirstHashDifference(const PicurvWindowDefinition *definition,
                                               const char *saved_group_digests,
                                               PetscInt *group)
{
    char digest[65];
    char current[PICURV_WINDOW_HASH_GROUP_COUNT][PICURV_WINDOW_HASH_GROUP_LENGTH];
    const char *cursor = saved_group_digests;

    PetscFunctionBeginUser;
    PetscCheck(definition != NULL && group != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Definition and group output are required.");
    *group = -1;
    if (!saved_group_digests) PetscFunctionReturn(0);
    PetscCall(PicurvWindowComputeHash(definition, digest, current));

    for (PetscInt index = 0; index < PICURV_WINDOW_HASH_GROUP_COUNT; ++index) {
        const char *comma = strchr(cursor, ',');
        const size_t length = comma ? (size_t)(comma - cursor) : strlen(cursor);
        size_t expected = 0;

        PetscCall(PetscStrlen(current[index], &expected));
        /* Every well-formed digest is fixed width, so a short segment means the
         * saved list is malformed rather than that this property changed. A
         * malformed list localizes to nothing: reporting the wrong property would
         * be worse than reporting none. */
        if (length != expected) PetscFunctionReturn(0);
        if (strncmp(cursor, current[index], length)) {
            *group = index;
            PetscFunctionReturn(0);
        }
        if (!comma) PetscFunctionReturn(0);
        cursor = comma + 1;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvWindowProgress().
 * @see PicurvWindowProgress()
 */
PetscReal PicurvWindowProgress(const PicurvWindow *window)
{
    PetscReal span = 0.0;

    if (window == NULL || !window->definition.bounded) return 0.0;
    span = window->definition.end_time - window->effective_start;
    if (span <= 0.0) return 0.0;
    return PetscMin(1.0, window->represented_time / span);
}

/**
 * @brief Implementation of \ref FieldStatisticsIsActive().
 * @see FieldStatisticsIsActive()
 */
PetscBool FieldStatisticsIsActive(const SimCtx *simCtx)
{
    if (!simCtx || !simCtx->fieldStatisticsEnabled) return PETSC_FALSE;
    return (PetscBool)(simCtx->fieldStatisticsWindowCount > 0 &&
                       simCtx->fieldStatisticsWindows != NULL);
}

/**
 * @brief Implementation of \ref FieldStatisticsUpdateWindows().
 * @see FieldStatisticsUpdateWindows()
 */
PetscErrorCode FieldStatisticsUpdateWindows(SimCtx *simCtx, PetscInt step, PetscReal time)
{
    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");
    if (!FieldStatisticsIsActive(simCtx)) PetscFunctionReturn(0);

    for (PetscInt window_index = 0; window_index < simCtx->fieldStatisticsWindowCount; ++window_index) {
        PicurvWindow *window = &simCtx->fieldStatisticsWindows[window_index];
        const PicurvWindowState before = window->state;
        PetscBool accepted = PETSC_FALSE;
        PetscReal weight = 0.0;

        PetscCall(PicurvWindowOfferState(window, step, time, &accepted, &weight));

        if (before == PICURV_WINDOW_PENDING && window->state == PICURV_WINDOW_ACTIVE) {
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "Statistics window '%s' active from t=%.6g.\n",
                      window->definition.name, (double)window->effective_start);
        }
        if (accepted) {
            /* Apply the accepted weight to every block's accumulators. Windows share
             * the source state but own independent accumulator state. */
            for (PetscInt bi = 0; bi < simCtx->block_number; ++bi) {
                UserCtx *user = &simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user[bi];
                if (!user->fieldStatisticsStorage) continue;
                PetscCall(PicurvWindowAccumulate(user, &window->definition,
                                                 &user->fieldStatisticsStorage[window_index], weight));
            }
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                      "Statistics window '%s' accepted step %d at t=%.6g with weight %.6g "
                      "(samples=%d, represented=%.6g).\n",
                      window->definition.name, step, (double)time, (double)weight,
                      window->sample_count, (double)window->represented_time);
        } else if (window->state == PICURV_WINDOW_ACTIVE) {
            /* An active window that took nothing was either off schedule or offered a
             * zero-length interval. Both are expected, so this is the level that
             * explains a sample count without implying anything is wrong. */
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                      "Statistics window '%s' did not sample step %d at t=%.6g "
                      "(last accepted t=%.6g, samples=%d).\n",
                      window->definition.name, step, (double)time,
                      (double)window->last_accepted_time, window->sample_count);
        }
        if (before != PICURV_WINDOW_COMPLETE && window->state == PICURV_WINDOW_COMPLETE) {
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "Statistics window '%s' complete: %d sample(s), total weight %.6g, "
                      "represented time %.6g.\n",
                      window->definition.name, window->sample_count,
                      (double)window->total_weight, (double)window->represented_time);
            if (window->sample_count == 0) {
                LOG_ALLOW(GLOBAL, LOG_WARNING,
                          "Statistics window '%s' completed without accepting any sample.\n",
                          window->definition.name);
            } else {
                UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

                /* A point that was never valid carries no average at all, which a
                 * reader of the mean field cannot see. Report it once, at the only
                 * moment the final coverage is known. The reduction is collective and
                 * runs once per window, not per step. */
                if (user && user->fieldStatisticsStorage) {
                    PetscReal lowest = 1.0, highest = 0.0;

                    PetscCall(PicurvWindowValidFractionRange(user, &window->definition,
                                                             &user->fieldStatisticsStorage[window_index],
                                                             window->sample_count, &lowest, &highest));
                    if (lowest <= 0.0) {
                        LOG_ALLOW(GLOBAL, LOG_WARNING,
                                  "Statistics window '%s' completed with points that were never "
                                  "valid; their mean and moments are undefined.\n",
                                  window->definition.name);
                    } else {
                        LOG_ALLOW(GLOBAL, LOG_INFO,
                                  "Statistics window '%s' per-point valid fraction spans [%.3f, %.3f].\n",
                                  window->definition.name, (double)lowest, (double)highest);
                    }
                }
            }
        }
    }
    PetscFunctionReturn(0);
}
