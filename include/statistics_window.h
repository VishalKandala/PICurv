/**
 * @file statistics_window.h
 * @brief Window lifecycle, scheduling, and weighting for the field-statistics pipeline.
 *
 * Implements the window semantics fixed in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification sections 2 and 4:
 * right-rectangle weighting, final-interval clipping, and the rule that a state
 * representing a zero-length interval is not a sample.
 *
 * This module decides **whether** a completed state is accepted and **what weight**
 * it carries. It holds no PETSc objects and performs no field accumulation; the
 * caller applies the returned weight through the moment kernels.
 */

#ifndef PICURV_STATISTICS_WINDOW_H
#define PICURV_STATISTICS_WINDOW_H

#include <petscsys.h>

struct SimCtx;

/** @brief Maximum stored length of a window name, including the terminator. */
#define PICURV_WINDOW_NAME_LENGTH 64

/** @brief Lifecycle state of one window. */
typedef enum {
    PICURV_WINDOW_PENDING = 0,  /**< Requested start not yet reached. */
    PICURV_WINDOW_ACTIVE,       /**< Accepting due states. */
    PICURV_WINDOW_COMPLETE      /**< Bounded end reached; accepts nothing further. */
} PicurvWindowState;

/** @brief How an accepted state's weight is determined. */
typedef enum {
    PICURV_WEIGHTING_SAMPLE = 0,    /**< Equal weight per accepted state. */
    PICURV_WEIGHTING_PHYSICAL_TIME  /**< Weight is the represented interval. */
} PicurvWeighting;

/** @brief Which schedule selects due states. Exactly one is used. */
typedef enum {
    PICURV_CADENCE_STEP = 0,  /**< Every n completed steps from activation. */
    PICURV_CADENCE_TIME       /**< First state at or past each nominal time target. */
} PicurvCadenceKind;

/** @brief The scientifically immutable definition of one window. */
typedef struct {
    char              name[PICURV_WINDOW_NAME_LENGTH];
    PetscReal         start_time;    /**< Requested start. */
    PetscReal         end_time;      /**< Requested end; ignored when @c bounded is false. */
    PetscBool         bounded;       /**< False for an open-ended window. */
    PicurvWeighting   weighting;
    PicurvCadenceKind cadence_kind;
    PetscInt          step_cadence;  /**< Used when cadence_kind is step; must be positive. */
    PetscReal         time_cadence;  /**< Used when cadence_kind is time; must be positive. */
} PicurvWindowDefinition;

/** @brief Runtime state of one window. */
typedef struct PicurvWindow {
    PicurvWindowDefinition definition;
    PicurvWindowState      state;
    PetscReal              effective_start;      /**< Origin of the first represented interval. */
    PetscReal              effective_end;        /**< End of the last represented interval. */
    PetscReal              last_accepted_time;   /**< Right edge of the last represented interval. */
    PetscInt               sample_count;
    PetscReal              total_weight;
    PetscReal              represented_time;     /**< Physical time the window covers. */
    PetscInt               activation_step;      /**< Step at which the window became active. */
    PetscInt               last_event_step;      /**< Guards against a step being offered twice. */
    PetscInt               next_time_target;     /**< k in effective_start + k*time_cadence. */
} PicurvWindow;

/**
 * @brief Validates a definition and initializes a window to the pending state.
 * @param[out] window     Window to initialize.
 * @param[in]  definition Requested definition; copied into the window.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for a non-positive cadence,
 *         an empty name, or a bounded window whose end does not exceed its start.
 */
PetscErrorCode PicurvWindowInit(PicurvWindow *window, const PicurvWindowDefinition *definition);

/**
 * @brief Offers one completed state to a window and reports the decision.
 *
 * Applies the interval convention in full: the state carries the interval ending
 * at it, measured from the previous accepted state or from the effective start;
 * a zero-length interval is not a sample; and a bounded window clips its final
 * interval to the requested end and then completes.
 *
 * When @p accepted is returned true the window's bookkeeping has already been
 * advanced, and the caller applies @p weight through the moment kernels. When it
 * is false the window is scientifically unchanged.
 *
 * Offering the same step twice is rejected, so a completed state cannot be
 * counted more than once.
 *
 * @param[in,out] window   Window to offer the state to.
 * @param[in]     step     Completed step number.
 * @param[in]     time     Physical time of the completed state.
 * @param[out]    accepted Whether the state became a sample.
 * @param[out]    weight   Weight to apply; zero when not accepted.
 * @return Zero on success, or a PETSc error for a null argument.
 */
PetscErrorCode PicurvWindowOfferState(PicurvWindow *window, PetscInt step, PetscReal time,
                                      PetscBool *accepted, PetscReal *weight);

/**
 * @brief Reports the fraction of a bounded window's span that has been represented.
 * @param[in] window Window to query.
 * @return Value in [0,1] for a bounded window, or zero for an open one.
 */
PetscReal PicurvWindowProgress(const PicurvWindow *window);

/**
 * @brief Returns a stable human-readable name for a window state.
 * @param[in] state Window lifecycle state.
 * @return Static string; never NULL.
 */
const char *PicurvWindowStateName(PicurvWindowState state);

/**
 * @brief Offers one completed state to every configured window.
 *
 * Called once per completed step from the runloop. Each due window advances its
 * own bookkeeping independently; windows share the source state but never share
 * accumulator state. Does nothing when field statistics are disabled or no
 * window is configured, which is the case until configuration ingress exists.
 *
 * @param[in,out] simCtx Simulation context carrying the window array.
 * @param[in]     step   Completed step number.
 * @param[in]     time   Physical time of the completed state.
 * @return Zero on success, or a PETSc error propagated from a window update.
 */
PetscErrorCode FieldStatisticsUpdateWindows(struct SimCtx *simCtx, PetscInt step, PetscReal time);

#endif /* PICURV_STATISTICS_WINDOW_H */
