/**
 * @file ParticleInitialConditions.h
 * @brief Configured initial values of particle-carried fields, and the expression
 *        language that defines them.
 *
 * A plan binds each configured particle field to one compiled expression per component.
 * Expressions are written in physical units: `x`, `y`, `z` are physical coordinates,
 * `t` is physical time, and a value is divided by its field's catalog reference scale
 * when written. The plan is applied to a contiguous local range of swarm entries, so the
 * t=0 population and, later, particles injected at a boundary share one path.
 *
 * The language is the one `generators/ic.gen` validates and evaluates for Eulerian
 * expressions; `tests/fixtures/expression_conformance.txt` holds cases both
 * implementations must agree on.
 */

#ifndef PARTICLE_INITIAL_CONDITIONS_H
#define PARTICLE_INITIAL_CONDITIONS_H

#include "variables.h"
#include "particle_field_catalog.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief A compiled expression: postfix code for a small stack evaluator. */
typedef struct PicurvExpression PicurvExpression;

/** @brief What a random draw is keyed on, besides the particle and its stream. */
typedef struct {
    PetscInt64 seed;   /**< Run seed (`-particle_random_seed`). */
    PetscInt64 pid;    /**< Particle identity; draws are pure functions of it. */
    PetscInt64 salt;   /**< Event that produced the particle: 0 for the t=0 population. */
    PetscInt64 field;  /**< Field whose expression draws, the default stream's namespace. */
} PicurvExpressionDrawKey;

/**
 * @brief Compile one expression against a list of variable names.
 * @details The grammar is Python's expression syntax restricted to numbers, the given
 *          names, `pi`, arithmetic (`+ - * / % **`), comparisons (chained ones included),
 *          `and`/`or`/`not`, and the functions `abs cos exp maximum minimum sin sqrt tan
 *          where`, plus `uniform` and `normal` when @p allow_random is true. Comparisons
 *          and logical operators yield 1 or 0; `%` takes the sign of the divisor; `**`
 *          binds tighter than a unary minus on its left and groups to the right.
 * @param[in]  text         Expression text.
 * @param[in]  names        Variable names, in the order evaluation supplies their values.
 * @param[in]  name_count   Number of names.
 * @param[in]  allow_random Whether `uniform` and `normal` may be called.
 * @param[out] expression   Compiled expression, owned by the caller.
 * @return Zero on success; `PETSC_ERR_ARG_WRONG` for a syntax error or an unknown name
 *         or function, with the offending position in the message.
 */
PetscErrorCode PicurvExpressionCompile(const char *text, const char *const *names, PetscInt name_count,
                                       PetscBool allow_random, PicurvExpression **expression);

/**
 * @brief Evaluate a compiled expression.
 * @param[in]  expression Compiled expression.
 * @param[in]  values     One value per compiled name, in compile order.
 * @param[in]  key        Draw key, required when the expression draws; may be NULL otherwise.
 * @param[out] result     Value of the expression; may be non-finite, which the caller judges.
 * @return Zero on success.
 */
PetscErrorCode PicurvExpressionEvaluate(const PicurvExpression *expression, const PetscReal *values,
                                        const PicurvExpressionDrawKey *key, PetscReal *result);

/**
 * @brief Free a compiled expression.
 * @param[in,out] expression Expression to free; set to NULL.
 * @return Zero on success.
 */
PetscErrorCode PicurvExpressionDestroy(PicurvExpression **expression);

/** @brief Compiled initial values for every configured particle field. */
typedef struct ParticleFieldPlan ParticleFieldPlan;

/** @brief When and why a plan is applied. */
typedef struct {
    PetscReal  physical_time;  /**< Time, in seconds, the evaluated particles appear at. */
    PetscInt64 salt;           /**< Event identity keying random draws: 0 for the t=0 population. */
} ParticleFieldEvent;

/**
 * @brief Read a plan from the options database.
 * @details Reads `-particle_fields_count`, and for each field `i`,
 *          `-particle_fields_<i>_name` and one expression per component,
 *          `-particle_fields_<i>_expr_<c>`. Every named field must be a particle-carried
 *          field of the catalog (capability `USER_INITIALIZE`).
 * @param[out] plan Compiled plan, or NULL when no field is configured.
 * @return Zero on success; `PETSC_ERR_ARG_WRONG` for an unknown or non-settable field or a
 *         malformed expression.
 */
PetscErrorCode ParticleFieldPlanCreate(ParticleFieldPlan **plan);

/**
 * @brief Write the plan's values into swarm entries `[first, end)` of this rank.
 * @details Not collective: each rank evaluates its own entries, and the domain bounds for
 *          `xn`, `yn`, `zn` come from the replicated rank bounding boxes. A value that is not
 *          finite is an error naming the particle and its position.
 * @param[in,out] user  Block owning the swarm.
 * @param[in]     plan  Plan to apply; NULL does nothing.
 * @param[in]     event Time and identity of the event producing these particles.
 * @param[in]     first First local entry.
 * @param[in]     end   One past the last local entry.
 * @return Zero on success.
 */
PetscErrorCode ParticleFieldPlanApply(UserCtx *user, const ParticleFieldPlan *plan,
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
