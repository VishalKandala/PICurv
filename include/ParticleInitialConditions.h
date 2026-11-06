/**
 * @file ParticleInitialConditions.h
 * @brief Configured initial values of particle-carried fields.
 *
 * A plan binds each configured particle field to one compiled expression per component.
 * Expressions are written in physical units: `x`, `y`, `z` are physical coordinates,
 * `t` is physical time, and a value is divided by its field's catalog reference scale
 * when written. The plan is applied to a contiguous local range of swarm entries, so the
 * t=0 population and, later, particles injected at a boundary share one path.
 *
 * The expressions are compiled and evaluated by the shared engine in expression.h.
 */

#ifndef PARTICLE_INITIAL_CONDITIONS_H
#define PARTICLE_INITIAL_CONDITIONS_H

#include "variables.h"
#include "particle_field_catalog.h"
#include "expression.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Compiled initial values for every configured particle field. */
typedef struct ParticleFieldPlan ParticleFieldPlan;

/** @brief When and why a plan is applied. */
typedef struct {
    PetscReal  physical_time;  /**< Time, in seconds, the evaluated particles appear at. */
    PetscInt64 salt;           /**< Event identity keying random draws: 0 for the t=0 population. */
} ParticleFieldEvent;

/**
 * @brief Read a plan from the options database.
 * @details Reads `-<prefix>_count`, and for each field `i`, `-<prefix>_<i>_name` and one
 *          expression per component, `-<prefix>_<i>_expr_<c>`. Every named field must be a
 *          particle-carried field of the catalog (capability `USER_INITIALIZE`).
 * @param[in]  prefix Option prefix without its dash: `particle_fields` for the t=0 population,
 *                    `particle_inlet_<n>_fields` for an inlet.
 * @param[out] plan Compiled plan, or NULL when no field is configured.
 * @return Zero on success; `PETSC_ERR_ARG_WRONG` for an unknown or non-settable field or a
 *         malformed expression.
 */
PetscErrorCode ParticleFieldPlanCreate(const char *prefix, ParticleFieldPlan **plan);

/**
 * @brief Write the plan's values into swarm entries `[first, end)` of this rank.
 * @details Each rank evaluates its own entries. Every rank calls it even with an empty range: the
 *          first call that reads a grid line's fraction or length is collective, computing them once
 *          for the plan. Those expressions need each particle's cell and weights. A value that is not finite is an error naming the
 *          particle and its position.
 * @param[in,out] user  Block owning the swarm.
 * @param[in]     plan  Plan to apply; NULL does nothing.
 * @param[in]     event Time and identity of the event producing these particles.
 * @param[in]     first First local entry.
 * @param[in]     end   One past the last local entry.
 * @return Zero on success.
 */
PetscErrorCode ParticleFieldPlanApply(UserCtx *user, ParticleFieldPlan *plan,
                                      const ParticleFieldEvent *event, PetscInt first, PetscInt end);

/**
 * @brief Log and record the realized initial values of every planned field.
 * @details Collective. Reports count, mean, variance, minimum, maximum, and a ten-bin
 *          histogram per field component, in physical units, to the log and to
 *          `particle_initial_fields.csv` in the run's metrics directory.
 * @param[in] user Block owning the swarm.
 * @param[in] plan Plan whose fields are summarized; NULL does nothing.
 * @return Zero on success.
 */
PetscErrorCode ParticleFieldPlanSummarize(UserCtx *user, const ParticleFieldPlan *plan);

/**
 * @brief Free a plan and its compiled expressions.
 * @param[in,out] plan Plan to free; set to NULL.
 * @return Zero on success.
 */
PetscErrorCode ParticleFieldPlanDestroy(ParticleFieldPlan **plan);

#ifdef __cplusplus
}
#endif

#endif /* PARTICLE_INITIAL_CONDITIONS_H */
