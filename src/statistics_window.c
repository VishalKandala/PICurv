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
 * @brief Implementation of \ref FieldStatisticsUpdateWindows().
 * @see FieldStatisticsUpdateWindows()
 */
PetscErrorCode FieldStatisticsUpdateWindows(SimCtx *simCtx, PetscInt step, PetscReal time)
{
    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");
    if (!simCtx->fieldStatisticsEnabled || simCtx->fieldStatisticsWindowCount <= 0) {
        PetscFunctionReturn(0);
    }
    PetscCheck(simCtx->fieldStatisticsWindows != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Field statistics are enabled but no window storage exists.");

    for (PetscInt w = 0; w < simCtx->fieldStatisticsWindowCount; ++w) {
        PicurvWindow *window = &simCtx->fieldStatisticsWindows[w];
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
                                                 &user->fieldStatisticsStorage[w], weight));
            }
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                      "Statistics window '%s' accepted step %d at t=%.6g with weight %.6g "
                      "(samples=%d, represented=%.6g).\n",
                      window->definition.name, step, (double)time, (double)weight,
                      window->sample_count, (double)window->represented_time);
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
            }
        }
    }
    PetscFunctionReturn(0);
}
