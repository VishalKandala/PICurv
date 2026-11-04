/**
 * @file expression.h
 * @brief The expression language shared by every configured field, and the keyed hashes
 *        its random draws and other reproducible draws are built from.
 *
 * Callers compile an expression against their own list of variable names and supply the
 * values in that order at evaluation. `generators/ic.gen` validates the same language;
 * `tests/fixtures/expression_conformance.txt` holds cases both implementations must agree on.
 */

#ifndef PICURV_EXPRESSION_H
#define PICURV_EXPRESSION_H

#include <stdint.h>
#include <petscsys.h>

#include "variables.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Disjoint key spaces of the keyed particle draws, so draws made for different
 *        purposes never coincide for one particle.
 */
typedef enum {
    PICURV_DRAW_SPACE_FIELD = 0,        /**< An expression's default, per-field stream. */
    PICURV_DRAW_SPACE_NAMED_STREAM = 1, /**< `uniform(k)`/`normal(k)`: a stream shared across fields. */
    PICURV_DRAW_SPACE_BROWNIAN = 2,     /**< Brownian increments; the salt is the step. */
    PICURV_DRAW_SPACE_PLACEMENT = 3     /**< Seeding points; the salt is the rejection attempt. */
} PicurvDrawSpace;

/**
 * @brief splitmix64 finalizer: a bijective mix of 64 bits.
 * @param[in] z Bits to mix.
 * @return The mixed bits.
 */
uint64_t PicurvMix64(uint64_t z);

/**
 * @brief Uniform on [0, 1) from the top 53 bits of a hash.
 * @param[in] hash Mixed bits.
 * @return A value in [0, 1).
 */
PetscReal PicurvUnitFromHash(uint64_t hash);

/**
 * @brief Hash that keys a particle draw: a pure function of its arguments, so a draw does
 *        not depend on the rank a particle lives on, the order particles are visited, or
 *        what was drawn before it.
 * @param[in] seed   Run seed (`-particle_random_seed`).
 * @param[in] space  Purpose of the draw.
 * @param[in] stream Stream within the space.
 * @param[in] pid    Particle identity.
 * @param[in] salt   Event that distinguishes repeated draws for one particle.
 * @return 64 mixed bits; derive further independent values with PicurvMix64().
 */
uint64_t PicurvParticleDrawHash(PetscInt64 seed, PicurvDrawSpace space, uint64_t stream,
                                PetscInt64 pid, PetscInt64 salt);

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
 *          `and`/`or`/`not`, and the functions `abs atan2 cos cosh erf exp floor log
 *          maximum minimum sin sinh sqrt tan tanh where`, plus `uniform` and `normal` when
 *          @p allow_random is true. Comparisons and logical operators yield 1 or 0; `%`
 *          takes the sign of the divisor; `**` binds tighter than a unary minus on its left
 *          and groups to the right.
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

/**
 * @brief Whether a compiled expression reads one of its variables.
 * @param[in]  expression Compiled expression.
 * @param[in]  index      Position of the variable in the compile-time name list.
 * @param[out] uses       True when some instruction loads that variable.
 * @return Zero on success.
 */
PetscErrorCode PicurvExpressionUsesVariable(const PicurvExpression *expression, PetscInt index, PetscBool *uses);

/* ------------------------------------------------------------------------- */
/*                         Prescribed velocity sources                        */
/* ------------------------------------------------------------------------- */

/** @brief Variables a prescribed velocity source may use, in the order of `values[]`. */
typedef enum {
    PICURV_SOURCE_VAR_X = 0,   /**< Physical x of the point [m]. */
    PICURV_SOURCE_VAR_Y,       /**< Physical y [m]. */
    PICURV_SOURCE_VAR_Z,       /**< Physical z [m]. */
    PICURV_SOURCE_VAR_T,       /**< Physical time [s]. */
    PICURV_SOURCE_VAR_XI,      /**< Arc-length fraction along the xi grid line, 0 at the - face. */
    PICURV_SOURCE_VAR_ETA,     /**< Arc-length fraction along the eta grid line. */
    PICURV_SOURCE_VAR_ZETA,    /**< Arc-length fraction along the zeta grid line. */
    PICURV_SOURCE_VAR_L_XI,    /**< Physical length of that xi grid line [m]. */
    PICURV_SOURCE_VAR_L_ETA,   /**< Physical length of that eta grid line [m]. */
    PICURV_SOURCE_VAR_L_ZETA,  /**< Physical length of that zeta grid line [m]. */
    PICURV_SOURCE_VAR_COUNT
} PicurvSourceVariable;

/** @brief How a source's values are read. */
typedef enum {
    PICURV_SOURCE_FRAME_CARTESIAN = 0,   /**< Three values: x, y, z velocity components. */
    PICURV_SOURCE_FRAME_CURVILINEAR = 1, /**< Three values along the unit xi, eta, zeta grid-line directions. */
    PICURV_SOURCE_FRAME_NORMAL = 2       /**< One value: the speed into the domain through a face. */
} PicurvSourceFrame;

/** @brief Where a velocity source is applied, which names its options. */
typedef enum {
    PICURV_SOURCE_TARGET_INITIAL = 0, /**< The Eulerian initial condition: `-ic_source_*`. */
    PICURV_SOURCE_TARGET_INLET = 1    /**< An inlet face: `-inlet_source_<id>_*`; the normal frame is allowed. */
} PicurvSourceTarget;

/** @brief A prescribed velocity: one compiled expression per value, and the frame they are in. */
typedef struct {
    PicurvSourceFrame frame;
    PetscInt          count;          /**< 1 for the normal frame, otherwise 3. */
    PicurvExpression *component[3];
    PetscBool         time_dependent; /**< Some component reads `t`. */
} PicurvVelocitySource;

/**
 * @brief Read and compile a velocity source from the options database.
 * @details Reads the target's `_frame` option (`cartesian`, `curvilinear`, or, on an
 *          inlet, `normal`) and one expression per value, `_expr_<c>`, written in physical
 *          units over the names of #PicurvSourceVariable. Random draws are not available.
 * @param[in]  target Target whose options are read.
 * @param[in]  id     Inlet source index (`-inlet_source_<id>_*`); ignored for the initial condition.
 * @param[out] source Compiled source; release it with PicurvVelocitySourceDestroy().
 * @return Zero on success; `PETSC_ERR_USER_INPUT` for a missing option or unknown frame,
 *         `PETSC_ERR_ARG_WRONG` for a malformed expression.
 */
PetscErrorCode PicurvVelocitySourceCreate(PicurvSourceTarget target, PetscInt id, PicurvVelocitySource *source);

/**
 * @brief Evaluate a source at one point.
 * @param[in]  source Compiled source.
 * @param[in]  values One value per #PicurvSourceVariable.
 * @param[in]  edges  Cell edge vectors along xi, eta, zeta (any length); read only in the
 *                    curvilinear frame.
 * @param[out] value  Physical Cartesian velocity, or for the normal frame the speed in `x`.
 * @return Zero on success; `PETSC_ERR_FP` for a non-finite value.
 */
PetscErrorCode PicurvVelocitySourceEvaluate(const PicurvVelocitySource *source, const PetscReal *values,
                                            const Cmpnts edges[3], Cmpnts *value);

/**
 * @brief Free a source's compiled expressions.
 * @param[in,out] source Source to release; left empty.
 * @return Zero on success.
 */
PetscErrorCode PicurvVelocitySourceDestroy(PicurvVelocitySource *source);

/** @brief Names of #PicurvSourceVariable, in order. */
extern const char *const kPicurvSourceVariableNames[PICURV_SOURCE_VAR_COUNT];


#ifdef __cplusplus
}
#endif

#endif /* PICURV_EXPRESSION_H */
