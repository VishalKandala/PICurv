/**
 * @file field_catalog.h
 * @brief Authoritative identities and storage metadata for persistent Eulerian fields.
 *
 * Field identifiers describe conceptual fields at compile time. They do not own,
 * allocate, or modify PETSc objects. A FieldView resolves one identifier against
 * an existing UserCtx after the normal setup lifecycle has created its DMs and
 * vectors.
 */

#ifndef FIELD_CATALOG_H
#define FIELD_CATALOG_H

#include <stddef.h>

#include "variables.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Compile-time identity for a catalogued Eulerian field. */
typedef enum {
    FIELD_ID_INVALID = -1,
    FIELD_ID_COORDINATES = 0,
    FIELD_ID_UCAT,
    FIELD_ID_UCONT,
    FIELD_ID_UCONT_O,
    FIELD_ID_UCONT_RM1,
    FIELD_ID_P,
    FIELD_ID_NU_T,
    FIELD_ID_CS,
    FIELD_ID_U_TAU,
    FIELD_ID_NU_WALL,
    FIELD_ID_DIFFUSIVITY,
    FIELD_ID_DIFFUSIVITY_GRADIENT,
    FIELD_ID_CSI,
    FIELD_ID_ETA,
    FIELD_ID_ZET,
    FIELD_ID_NVERT,
    FIELD_ID_AJ,
    FIELD_ID_CENT,
    FIELD_ID_GRID_SPACE,
    FIELD_ID_CENTX,
    FIELD_ID_CENTY,
    FIELD_ID_CENTZ,
    FIELD_ID_ICSI,
    FIELD_ID_IETA,
    FIELD_ID_IZET,
    FIELD_ID_JCSI,
    FIELD_ID_JETA,
    FIELD_ID_JZET,
    FIELD_ID_KCSI,
    FIELD_ID_KETA,
    FIELD_ID_KZET,
    FIELD_ID_IAJ,
    FIELD_ID_JAJ,
    FIELD_ID_KAJ,
    FIELD_ID_PHI,
    FIELD_ID_PSI,
    FIELD_ID_NVERT_O,
    FIELD_ID_PARTICLE_COUNT,
    FIELD_ID_CELL_SCALAR_AT_CORNER,
    FIELD_ID_CELL_VECTOR_AT_CORNER,
    FIELD_ID_POST_SCALAR,
    FIELD_ID_POST_VECTOR,
    FIELD_ID_QCRIT,
    FIELD_ID_COUNT
} FieldId;

/** @brief Logical storage topology of a field. */
typedef enum {
    FIELD_LAYOUT_NODE_CENTERED = 0,
    FIELD_LAYOUT_CELL_CENTERED,
    FIELD_LAYOUT_I_FACE,
    FIELD_LAYOUT_J_FACE,
    FIELD_LAYOUT_K_FACE,
    FIELD_LAYOUT_COMPONENT_STAGGERED
} FieldLayout;

/** @brief UserCtx DM family used to store a field. */
typedef enum {
    FIELD_DM_DA = 0,
    FIELD_DM_FDA,
    FIELD_DM_COORDINATES
} FieldDMKind;

/** @brief Extra repair required after the normal PETSc global-to-local scatter. */
typedef enum {
    FIELD_SYNC_STANDARD = 0,
    FIELD_SYNC_I_FACE,
    FIELD_SYNC_J_FACE,
    FIELD_SYNC_K_FACE,
    FIELD_SYNC_COMPONENT_STAGGERED
} FieldSyncClass;

/** @brief Conditions controlling when a field can have runtime storage. */
typedef enum {
    FIELD_AVAILABILITY_ALWAYS       = 0u,
    FIELD_AVAILABILITY_FINEST_LEVEL = 1u << 0,
    FIELD_AVAILABILITY_TURBULENCE   = 1u << 1,
    FIELD_AVAILABILITY_LES_DYNAMIC  = 1u << 2,
    FIELD_AVAILABILITY_PARTICLES    = 1u << 4,
    FIELD_AVAILABILITY_WALL_MODEL   = 1u << 5
} FieldAvailability;

/** @brief Operations supported by a catalog entry. */
typedef enum {
    FIELD_CAPABILITY_NONE                      = 0u,
    FIELD_CAPABILITY_GHOST_UPDATE              = 1u << 0,
    FIELD_CAPABILITY_PERIODIC_CELL_SYNC        = 1u << 1,
    FIELD_CAPABILITY_PERIODIC_FACE_SYNC        = 1u << 2,
    FIELD_CAPABILITY_PERIODIC_STAGGERED_SYNC   = 1u << 3,
    FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT   = 1u << 4,
    FIELD_CAPABILITY_CHECKPOINT                = 1u << 5
} FieldCapabilities;

/** @brief How a field's physical dimension is known. */
typedef enum {
    FIELD_DIMENSION_FIXED = 0,       /**< The exponents are the field's dimension. */
    FIELD_DIMENSION_NOT_A_QUANTITY,  /**< Integer bookkeeping: identities, flags, ranks. */
    FIELD_DIMENSION_FROM_SOURCE      /**< Staging storage: carries whatever was written into it. */
} FieldDimensionKind;

/**
 * @brief Physical dimension as exponents of the reference length, velocity, and density.
 *
 * A solver value of dimension (a, b, c) becomes physical when multiplied by
 * `L_ref^a U_ref^b rho_ref^c`. The launcher records input dimensions with the same
 * triple (`picurv_cli/core.py`, `INPUT_QUANTITIES`), so input and output conversion use
 * one vocabulary. See docs/pages/19_Nondimensionalization.md.
 */
typedef struct {
    FieldDimensionKind kind;
    signed char        length;
    signed char        velocity;
    signed char        density;
} FieldDimension;

/* Named dimensions for catalog entries. Each expands to a brace initializer, so it can
 * be passed through an entry macro as one argument. */
#define FIELD_DIM_DIMENSIONLESS        {FIELD_DIMENSION_FIXED, 0, 0, 0}
#define FIELD_DIM_LENGTH               {FIELD_DIMENSION_FIXED, 1, 0, 0}
#define FIELD_DIM_AREA                 {FIELD_DIMENSION_FIXED, 2, 0, 0}
#define FIELD_DIM_INVERSE_VOLUME       {FIELD_DIMENSION_FIXED, -3, 0, 0}
#define FIELD_DIM_TIME                 {FIELD_DIMENSION_FIXED, 1, -1, 0}
#define FIELD_DIM_VELOCITY             {FIELD_DIMENSION_FIXED, 0, 1, 0}
#define FIELD_DIM_VOLUME_FLUX          {FIELD_DIMENSION_FIXED, 2, 1, 0}
#define FIELD_DIM_DIFFUSIVITY          {FIELD_DIMENSION_FIXED, 1, 1, 0}
#define FIELD_DIM_PRESSURE             {FIELD_DIMENSION_FIXED, 0, 2, 1}
#define FIELD_DIM_INVERSE_TIME_SQUARED {FIELD_DIMENSION_FIXED, -2, 2, 0}
#define FIELD_DIM_NOT_A_QUANTITY       {FIELD_DIMENSION_NOT_A_QUANTITY, 0, 0, 0}
#define FIELD_DIM_FROM_SOURCE          {FIELD_DIMENSION_FROM_SOURCE, 0, 0, 0}

/** @brief Immutable metadata for one field identity. */
typedef struct {
    FieldId           id;
    const char       *canonical_name;
    const char       *alias_1;
    const char       *alias_2;
    PetscInt          dof;
    FieldDMKind       dm_kind;
    FieldLayout       layout;
    FieldSyncClass    sync_class;
    unsigned int      availability;
    unsigned int      capabilities;
    size_t            global_vec_offset;
    size_t            local_vec_offset;
    FieldDimension    dimension;
} FieldDescriptor;

/** @brief Non-owning runtime objects resolved for one field and UserCtx. */
typedef struct {
    const FieldDescriptor *descriptor;
    DM                     dm;
    Vec                    global_vec;
    Vec                    local_vec;
} FieldView;

/**
 * @brief Return immutable metadata for a valid field identifier.
 * @param[in]  field_id    Compile-time field identity.
 * @param[out] descriptor  Catalog entry owned by the field-catalog module.
 * @return Zero on success; PETSc error for an invalid ID or null output.
 */
PetscErrorCode FieldGetDescriptor(FieldId field_id, const FieldDescriptor **descriptor);

/**
 * @brief Resolve a user-facing field name once into its typed identity.
 * @param[in]  field_name Canonical name or registered alias.
 * @param[out] field_id   Resolved field identity.
 * @return Zero on success; PETSc unknown-type error for an unregistered name.
 */
PetscErrorCode FieldIdFromName(const char *field_name, FieldId *field_id);

/**
 * @brief Look a name up without treating an unknown name as an error.
 * @details For callers that consult more than one catalog, such as a name that may be an
 *          Eulerian or a particle field.
 * @param[in]  field_name Canonical name or registered alias.
 * @param[out] field_id   Resolved identity, or `FIELD_ID_INVALID` when not found.
 * @param[out] found      Whether the name is registered.
 * @return Zero on success; PETSc error only for null arguments.
 */
PetscErrorCode FieldTryIdFromName(const char *field_name, FieldId *field_id, PetscBool *found);

/**
 * @brief Return the factor that turns a solver value of one dimension into physical units.
 * @param[in]  scaling   Reference scales of the run.
 * @param[in]  dimension Dimension to scale.
 * @param[out] scale     `L_ref^a U_ref^b rho_ref^c` for a fixed dimension `(a, b, c)`.
 * @return Zero on success; `PETSC_ERR_ARG_WRONGSTATE` for a dimension that is not fixed,
 *         since bookkeeping and staging storage have no scale of their own.
 */
PetscErrorCode FieldDimensionReferenceScale(const ScalingCtx *scaling, FieldDimension dimension,
                                            PetscReal *scale);

/**
 * @brief Write a printable form of a dimension, such as `L^2 U` or `dimensionless`.
 * @param[in]  dimension Dimension to describe.
 * @param[out] label     Destination buffer.
 * @param[in]  length    Buffer length.
 * @return Zero on success.
 */
PetscErrorCode FieldDimensionLabel(FieldDimension dimension, char *label, size_t length);

/**
 * @brief Return the canonical printable name for an ID.
 * @param[in] field_id Field identity.
 * @return Canonical name, or "InvalidField" for an invalid ID.
 */
const char *FieldCanonicalName(FieldId field_id);

/**
 * @brief Return a stable printable label for a field layout.
 * @param[in] layout Layout enum to describe.
 * @return Catalog-owned label, or `Invalid-Layout` for an invalid enum.
 */
const char *FieldLayoutName(FieldLayout layout);

/**
 * @brief Resolve the existing DM and global/local vectors for one field.
 *
 * This function never creates storage. Optional fields whose setup conditions
 * were not enabled retain valid descriptors but return PETSC_ERR_ARG_WRONGSTATE
 * because their runtime vectors are absent.
 *
 * @param[in]  user      Existing per-grid context.
 * @param[in]  field_id  Field identity to resolve.
 * @param[out] view      Non-owning runtime view.
 * @return Zero on success or a PETSc argument/state error.
 */
PetscErrorCode FieldGetView(UserCtx *user, FieldId field_id, FieldView *view);

#ifdef __cplusplus
}
#endif

#endif /* FIELD_CATALOG_H */
