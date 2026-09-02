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
    FIELD_ID_K_OMEGA,
    FIELD_ID_K_OMEGA_O,
    FIELD_ID_CELL_SCALAR_AT_CORNER,
    FIELD_ID_CELL_VECTOR_AT_CORNER,
    FIELD_ID_POST_SCALAR,
    FIELD_ID_POST_VECTOR,
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
    FIELD_DM_FDA2,
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
    FIELD_AVAILABILITY_RANS         = 1u << 3,
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
