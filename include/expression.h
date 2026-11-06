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
    PICURV_SOURCE_FRAME_NORMAL = 2,      /**< One value: the speed into the domain through a face. */
    PICURV_SOURCE_FRAME_CONTRAVARIANT = 3 /**< Three values: the volume flux through the xi, eta, zeta faces; cells only. */
} PicurvSourceFrame;

/**
 * @brief A prescribed velocity, for any target: one compiled expression per value, the frame
 *        they are in, and the bulk velocity it is scaled to, if any. A target differs only in
 *        the points it evaluates the source at and the weight it takes a mean with.
 */
typedef struct {
    PicurvSourceFrame frame;
    PetscInt          count;          /**< 1 for the normal frame, otherwise 3; 0 for staged data alone. */
    PicurvExpression *component[3];
    PetscBool         time_dependent; /**< Some component reads `t`. */
    PicurvExpression *bulk;           /**< Target mean speed [m/s] in `t`, or NULL. */
    PetscBool         bulk_time_dependent;
} PicurvVelocitySource;

/**
 * @brief Read and compile a velocity source from the options database.
 * @details Reads `-<prefix>_frame` (`cartesian`, `curvilinear`, on a face `normal`, otherwise `contravariant`), one
 *          expression per value, `-<prefix>_expr_<c>`, written in physical units over the names
 *          of #PicurvSourceVariable, and the optional `-<prefix>_bulk_velocity`, an expression
 *          in `t`. Random draws are not available.
 * @param[in]  prefix   Option prefix without its dash, e.g. `ic_source` or `inlet_source_0`.
 * @param[in]  on_face  Whether the target is a face, where the normal frame means something.
 * @param[in]  has_data Whether the target also has staged data, which makes the expressions
 *                      optional; without a frame, data on a face is a normal speed.
 * @param[out] source   Compiled source; release it with PicurvVelocitySourceDestroy().
 * @return Zero on success; `PETSC_ERR_USER_INPUT` for a missing option or unknown frame,
 *         `PETSC_ERR_ARG_WRONG` for a malformed expression.
 */
PetscErrorCode PicurvVelocitySourceCreate(const char *prefix, PetscBool on_face, PetscBool has_data,
                                          PicurvVelocitySource *source);

/**
 * @brief Factor on a source's expressions that makes a target's whole mean equal its bulk velocity.
 * @details Solves `factor * expressions + data = bulk`, over the target's points and weight.
 *          Collective.
 * @param[in]  source        Source with a bulk velocity.
 * @param[in]  physical_time Time the bulk velocity is read at [s].
 * @param[in]  local         This rank's weighted sums, physical: the expressions' speed, the
 *                           staged data's speed along it, and the weight.
 * @param[in]  where         The target, for the error messages.
 * @param[out] factor        Factor on the expressions.
 * @return Zero on success; `PETSC_ERR_ARG_WRONGSTATE` for no weight or a zero mean.
 */
PetscErrorCode PicurvVelocitySourceBulkFactor(const PicurvVelocitySource *source, PetscReal physical_time,
                                              const PetscReal local[3], const char *where, PetscReal *factor);

/**
 * @brief Fill a source's position and grid-line variables at one point.
 * @param[out] values   Array of #PICURV_SOURCE_VAR_COUNT; `t` is left to the caller.
 * @param[in]  point    Solver-unit position.
 * @param[in]  fraction Grid-line fractions there, from ComputeGridLineFractions().
 * @param[in]  length   Solver-unit grid-line lengths there.
 * @param[in]  L_ref    Length reference [m].
 */
static inline void PicurvSourceSetVariables(PetscReal *values, Cmpnts point, Cmpnts fraction, Cmpnts length,
                                            PetscReal L_ref)
{
    values[PICURV_SOURCE_VAR_X] = point.x * L_ref;
    values[PICURV_SOURCE_VAR_Y] = point.y * L_ref;
    values[PICURV_SOURCE_VAR_Z] = point.z * L_ref;
    values[PICURV_SOURCE_VAR_XI] = fraction.x;
    values[PICURV_SOURCE_VAR_ETA] = fraction.y;
    values[PICURV_SOURCE_VAR_ZETA] = fraction.z;
    values[PICURV_SOURCE_VAR_L_XI] = length.x * L_ref;
    values[PICURV_SOURCE_VAR_L_ETA] = length.y * L_ref;
    values[PICURV_SOURCE_VAR_L_ZETA] = length.z * L_ref;
}

/** @brief What a source reads at each boundary-face slot this rank owns, gathered once. */
typedef struct {
    PetscInt   normal;  /**< Logical axis of the face normal. */
    PetscInt   a, b;    /**< The two tangential axes. */
    PetscInt   layer;   /**< Index along the normal of the face's metric and flux. */
    PetscInt   lo[2];   /**< First tangential indices of the owned slots. */
    PetscInt   n[2];    /**< Owned slots along the two tangential axes; 0 off the face. */
    PetscReal (*values)[PICURV_SOURCE_VAR_COUNT]; /**< Variables per slot; `t` is the caller's. */
    Cmpnts    (*edges)[3];                        /**< Directions per slot, for the curvilinear frame. */
} PicurvSourceFaceSlots;

/**
 * @brief Gather each owned slot of a boundary face: face-centre position, the grid-line fractions
 *        and lengths of the interior cell behind it, and its directions. A curvilinear component
 *        is the speed through the grid face, as on every target, so along the face's own axis
 *        the direction follows the face metric rather than the cell's edge.
 * @param[in]  user     Block with current metric and face-centre ghosts.
 * @param[in]  face_id  Face.
 * @param[in]  fraction Grid-line fractions, from ComputeGridLineFractions().
 * @param[in]  length   Grid-line lengths, from ComputeGridLineFractions().
 * @param[out] slots    The slots; release with PicurvSourceFaceSlotsDestroy().
 * @return PETSc error code.
 */
PetscErrorCode PicurvSourceFaceSlotsCreate(UserCtx *user, BCFace face_id, Vec fraction, Vec length,
                                           PicurvSourceFaceSlots *slots);

/**
 * @brief Release what PicurvSourceFaceSlotsCreate() allocated.
 * @param[in,out] slots The slots.
 * @return PETSc error code.
 */
PetscErrorCode PicurvSourceFaceSlotsDestroy(PicurvSourceFaceSlots *slots);

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
 * @brief Factor on a source's expressions that makes the volume mean of the whole field, along
 *        the source's own direction at each cell, equal its bulk velocity. Collective.
 * @param[in]  user          Block; its Ucat holds whatever is staged.
 * @param[in]  source        Source with a bulk velocity.
 * @param[in]  fraction      Grid-line fractions, from ComputeGridLineFractions().
 * @param[in]  length        Grid-line lengths, from ComputeGridLineFractions().
 * @param[in]  physical_time Time [s].
 * @param[out] factor        Factor on the expressions.
 * @return PETSc error code.
 */
PetscErrorCode PicurvSourceCellBulkFactor(UserCtx *user, const PicurvVelocitySource *source, Vec fraction, Vec length,
                                      PetscReal physical_time, PetscReal *factor);

/**
 * @brief Add a source's flux to every owned face of Ucont, evaluating it at each face centre.
 * @details A Cartesian velocity is dotted with the face metric; a curvilinear component is the
 *          speed through the face, times its area; a contravariant component is the flux itself.
 * @param[in,out] user           Block whose Ucont receives the fluxes.
 * @param[in]     source         Source with expressions.
 * @param[in]     fraction       Grid-line fractions, from ComputeGridLineFractions().
 * @param[in]     length         Grid-line lengths, from ComputeGridLineFractions().
 * @param[in]     physical_time  Time [s].
 * @param[in]     factor         Factor on the expressions.
 * @param[in]     boundary_faces Whether the boundary - faces are set too; each cell's + face always is.
 * @return PETSc error code.
 */
PetscErrorCode PicurvSourceAddFluxes(UserCtx *user, const PicurvVelocitySource *source, Vec fraction, Vec length,
                                     PetscReal physical_time, PetscReal factor, PetscBool boundary_faces);

/**
 * @brief Free a source's compiled expressions and bulk velocity.
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
