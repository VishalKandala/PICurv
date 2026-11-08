/**
 * @file expression.h
 * @brief The expression language shared by every configured field, the keyed hashes its
 *        random draws and other reproducible draws are built from, and the prescribed velocity
 *        sources written in it: their compilation, the face slots a boundary evaluates them at,
 *        the face fluxes and cell values they set, and their bulk-velocity scaling.
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
 * @brief Disjoint key spaces of the keyed draws, so draws made for different purposes never
 *        coincide for one particle, grid cell, or face cell.
 */
typedef enum {
    PICURV_DRAW_SPACE_FIELD = 0,        /**< An expression's default, per-field stream. */
    PICURV_DRAW_SPACE_NAMED_STREAM = 1, /**< `uniform(k)`/`normal(k)`: a stream shared across fields. */
    PICURV_DRAW_SPACE_BROWNIAN = 2,     /**< Brownian increments; the salt is the step. */
    PICURV_DRAW_SPACE_PLACEMENT = 3,    /**< Seeding points; the salt is the rejection attempt. */
    PICURV_DRAW_SPACE_INLET = 4,        /**< Injection at an inlet: stream 0 counts per face cell, stream 1 distance in. */
    PICURV_DRAW_SPACE_GRID_FIELD = 5,   /**< A velocity or pressure source's default stream, per component. */
    PICURV_DRAW_SPACE_GRID_NAMED_STREAM = 6, /**< `uniform(k)`/`normal(k)` in such a source: shared across components. */
    PICURV_DRAW_SPACE_INFLOW = 7        /**< Synthetic inflow noise: stream per block and component, salt the step. */
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
 * @brief Two independent standard normal deviates from a hash, by Box-Muller from two uniforms it yields.
 * @param[in]  hash   Hash bits; distinct bits give independent pairs.
 * @param[out] normal The pair.
 */
void PicurvNormalPairFromHash(uint64_t hash, PetscReal normal[2]);

/**
 * @brief A standard normal deviate from a hash: the first of PicurvNormalPairFromHash()'s pair.
 * @param[in] hash Hash bits; distinct bits give independent deviates.
 * @return The deviate.
 */
PetscReal PicurvNormalFromHash(uint64_t hash);

/**
 * @brief Hash that keys a draw: a pure function of its arguments, so a draw does not depend
 *        on the rank that makes it, the order points are visited, or what was drawn before it.
 * @param[in] seed   Run seed.
 * @param[in] space  Purpose of the draw.
 * @param[in] stream Stream within the space.
 * @param[in] id     Global identity of the particle, grid cell, or face cell.
 * @param[in] salt   Event that distinguishes repeated draws for one identity.
 * @return 64 mixed bits; derive further independent values with PicurvMix64().
 */
uint64_t PicurvDrawHash(PetscInt64 seed, PicurvDrawSpace space, uint64_t stream,
                                PetscInt64 id, PetscInt64 salt);

/** @brief A compiled expression: postfix code for a small stack evaluator. */
typedef struct PicurvExpression PicurvExpression;

/** @brief What a random draw is keyed on, besides its stream. */
typedef struct {
    PetscInt64 seed;   /**< Run seed: `-particle_random_seed`, or `-source_random_seed` on the grid. */
    PetscInt64 pid;    /**< Global identity of the particle, or on the grid of the cell or face cell. */
    PetscInt64 salt;   /**< Event: for a particle what produced it, 0 for the t=0 population; on a
                            face the step; in a cell #PICURV_DRAW_SALT_INITIAL_FIELD. */
    PetscInt64 field;  /**< Field or component whose expression draws, the default stream's namespace. */
    PetscBool  grid;   /**< The draw is for a grid point, not a particle. */
} PicurvExpressionDrawKey;

/** @brief Salt of a cell's draws: the field is set once, and no step has this number. */
#define PICURV_DRAW_SALT_INITIAL_FIELD ((PetscInt64)-1)

/**
 * @brief Compile one expression against a list of variable names.
 * @details The grammar is Python's expression syntax restricted to numbers, the given
 *          names, `pi`, arithmetic (`+ - * / % **`), comparisons (chained ones included),
 *          `and`/`or`/`not`, the functions `abs atan2 bei ber besselj bessely cos cosh erf
 *          exp floor log maximum minimum sin sinh sqrt tan tanh where`, the series
 *          `sum(k, lo, hi, body)` over integer bounds, plus `uniform` and `normal` when
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
 * @brief Whether an expression draws random numbers.
 * @param[in]  expression Compiled expression.
 * @param[out] draws      True when it calls `uniform` or `normal`.
 * @return Zero on success.
 */
PetscErrorCode PicurvExpressionDraws(const PicurvExpression *expression, PetscBool *draws);

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
    PICURV_SOURCE_VAR_CELL_COUNT, /**< Variables a source reads on every target. */
    PICURV_SOURCE_VAR_R = PICURV_SOURCE_VAR_CELL_COUNT, /**< On a face, distance from the face's centre [m]. */
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
    PetscBool         time_dependent; /**< Some component reads `t` or draws, so a face re-evaluates it each step. */
    PicurvExpression *bulk;           /**< Target mean speed [m/s] in `t`, or NULL. */
    PetscBool         bulk_time_dependent;
    PetscBool         bulk_is_flow_rate; /**< `bulk` is a volumetric flow rate [m^3/s] through a face instead. */
} PicurvVelocitySource;

/**
 * @brief Read and compile a velocity source from the options database.
 * @details Reads `-<prefix>_frame` (`cartesian`, `curvilinear`, on a face `normal`, otherwise `contravariant`), one
 *          expression per value, `-<prefix>_expr_<c>`, written in physical units over the names
 *          of #PicurvSourceVariable, and the optional `-<prefix>_bulk_velocity`, an expression
 *          in `t`; on a face, without a bulk velocity, the optional `-<prefix>_flow_rate`, a
 *          volumetric flow rate [m^3/s] in `t`, instead.
 * @param[in]  prefix   Option prefix without its dash, e.g. `ic_source` or `face_source_0`.
 * @param[in]  on_face  Whether the target is a face, where the normal frame means something.
 * @param[in]  has_data Whether the target also has staged data, which makes the expressions
 *                      optional; without a frame, data on a face is a normal speed.
 * @param[in]  allow_draws Whether the expressions may call `uniform` and `normal`; a target
 *                      evaluated at arbitrary points has nothing to key them on.
 * @param[out] source   Compiled source; release it with PicurvVelocitySourceDestroy().
 * @return Zero on success; `PETSC_ERR_USER_INPUT` for a missing option or unknown frame,
 *         `PETSC_ERR_ARG_WRONG` for a malformed expression.
 */
PetscErrorCode PicurvVelocitySourceCreate(const char *prefix, PetscBool on_face, PetscBool has_data,
                                          PetscBool allow_draws, PicurvVelocitySource *source);

/**
 * @brief Read and compile a source's optional pressure, `-<prefix>_pressure`, an expression in Pa
 *        over the cell variables of #PicurvSourceVariable.
 * @param[in]  prefix      Option prefix without its dash, e.g. `ic_source`.
 * @param[in]  allow_draws Whether the expression may call `uniform` and `normal`.
 * @param[out] pressure    Compiled expression, or NULL when the option is absent; release it with
 *                         PicurvExpressionDestroy().
 * @return Zero on success; `PETSC_ERR_ARG_WRONG` for a malformed expression.
 */
PetscErrorCode PicurvSourcePressureCreate(const char *prefix, PetscBool allow_draws, PicurvExpression **pressure);

/**
 * @brief Factor on a source's expressions that makes a target's whole mean equal its bulk velocity.
 * @details Solves `factor * expressions + data = bulk`, over the target's points and weight.
 *          Collective.
 * @param[in]  source        Source with a bulk velocity.
 * @param[in]  physical_time Time the bulk velocity is read at [s].
 * @param[in]  local         This rank's weighted sums, physical: the expressions' speed, the
 *                           staged data's speed along it, and the weight, an area [m^2] for a
 *                           flow rate.
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

/**
 * @brief Global identity of a boundary face cell, unique over the six faces of a block and
 *        independent of the decomposition.
 * @param[in] nodes Global node counts along i, j, k.
 * @param[in] face  Face.
 * @param[in] p     Zero-based cell index along the first tangential axis.
 * @param[in] q     Zero-based cell index along the second.
 * @return The identity, from 0 below `2 (nx ny + ny nz + nz nx)` for `nx ny nz` cells.
 */
static inline PetscInt64 PicurvFaceCellId(const PetscInt nodes[3], BCFace face, PetscInt p, PetscInt q)
{
    const PetscInt64 cells[3] = {nodes[0] - 2, nodes[1] - 2, nodes[2] - 2};
    PetscInt64 id = p + cells[((PetscInt)face / 2 == 0) ? 1 : 0] * q;

    for (PetscInt f = 0; f < (PetscInt)face; ++f) id += cells[(f / 2 + 1) % 3] * cells[(f / 2 + 2) % 3];
    return id;
}

/**
 * @brief Draw identity of a grid point: its block in the top 16 bits, above an identity unique
 *        within the block (a cell's natural index, or PicurvFaceCellId()).
 * @param[in] block Block index.
 * @param[in] id    Identity within the block, below 2^48.
 * @return The identity.
 */
static inline PetscInt64 PicurvGridDrawId(PetscInt block, PetscInt64 id)
{
    return ((PetscInt64)block << 48) + id;
}

/**
 * @brief Draw key of a boundary face cell at one step, the key of a face source's draws.
 * @param[in] seed  Run seed, `-source_random_seed`.
 * @param[in] block Block index.
 * @param[in] nodes Global node counts along i, j, k.
 * @param[in] face  Face.
 * @param[in] p     Zero-based cell index along the face's first tangential axis.
 * @param[in] q     Zero-based cell index along the second.
 * @param[in] step  Step the draws are made for, their salt.
 * @return The key; its `field` is 0, for PicurvVelocitySourceEvaluate() to set per component.
 */
static inline PicurvExpressionDrawKey PicurvFaceDrawKey(PetscInt64 seed, PetscInt block, const PetscInt nodes[3],
                                                        BCFace face, PetscInt p, PetscInt q, PetscInt64 step)
{
    const PicurvExpressionDrawKey key = {seed, PicurvGridDrawId(block, PicurvFaceCellId(nodes, face, p, q)), step, 0,
                                         PETSC_TRUE};
    return key;
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
 * @param[in]  key    The point's draw key, whose `field` this sets to each component; may be
 *                    NULL for a source that does not draw.
 * @param[in]  edges  Cell edge vectors along xi, eta, zeta (any length); read only in the
 *                    curvilinear frame.
 * @param[out] value  Physical Cartesian velocity, or for the normal frame the speed in `x`.
 * @return Zero on success; `PETSC_ERR_FP` for a non-finite value.
 */
PetscErrorCode PicurvVelocitySourceEvaluate(const PicurvVelocitySource *source, const PetscReal *values,
                                            const PicurvExpressionDrawKey *key, const Cmpnts edges[3], Cmpnts *value);

/**
 * @brief Set a block's Ucat, its P, or both at every interior cell centre this rank owns, from a
 *        source evaluated there. Dummy cells are left to the caller's boundary pass.
 * @details A cell's draws are keyed as the face fluxes it sets are, the pressure on the stream
 *          after the velocity's three components.
 * @param[in,out] user          Block with current cell centres; its Ucat and P receive the values.
 * @param[in]     velocity      Source of Cartesian components, physical, or NULL to leave Ucat.
 * @param[in]     factor        Factor on the velocity, from a bulk velocity, or 1.
 * @param[in]     pressure      Expression over the cell variables in Pa, or NULL to leave P.
 * @param[in]     fraction      Grid-line fractions, from ComputeGridLineFractions().
 * @param[in]     length        Grid-line lengths, from ComputeGridLineFractions().
 * @param[in]     physical_time Time [s].
 * @return PETSc error code.
 */
PetscErrorCode PicurvSourceSetCellFields(UserCtx *user, const PicurvVelocitySource *velocity, PetscReal factor,
                                         const PicurvExpression *pressure, Vec fraction, Vec length,
                                         PetscReal physical_time);

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
