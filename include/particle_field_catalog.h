/**
 * @file particle_field_catalog.h
 * @brief Typed identities and metadata for persistent solver-particle fields.
 *
 * Particle fields intentionally use a catalog separate from Eulerian fields.
 * They are registered on a DMSwarm, may use integer storage, and do not own a
 * global/local Vec pair in UserCtx.  Text is translated only where PETSc or a
 * user/file interface requires a registered DMSwarm name.
 */

#ifndef PARTICLE_FIELD_CATALOG_H
#define PARTICLE_FIELD_CATALOG_H

#include <petscdmswarm.h>

#include "field_catalog.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Compile-time identity for a persistent solver-particle field. */
typedef enum {
    PARTICLE_FIELD_ID_INVALID = -1,
    PARTICLE_FIELD_ID_POSITION = 0,
    PARTICLE_FIELD_ID_VELOCITY,
    PARTICLE_FIELD_ID_CELL_ID,
    PARTICLE_FIELD_ID_WEIGHT,
    PARTICLE_FIELD_ID_DIFFUSIVITY,
    PARTICLE_FIELD_ID_DIFFUSIVITY_GRADIENT,
    PARTICLE_FIELD_ID_PSI,
    PARTICLE_FIELD_ID_LOCATION_STATUS,
    PARTICLE_FIELD_ID_PID,
    PARTICLE_FIELD_ID_RANK,
    PARTICLE_FIELD_ID_COUNT
} ParticleFieldId;

/** @brief Owner of a field's DMSwarm registration. */
typedef enum {
    PARTICLE_FIELD_REGISTRATION_PICURV = 0,
    PARTICLE_FIELD_REGISTRATION_PETSC
} ParticleFieldRegistration;

/** @brief Operations supported by a particle-field catalog entry. */
typedef enum {
    PARTICLE_FIELD_CAPABILITY_NONE               = 0u,
    PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE = 1u << 0,
    PARTICLE_FIELD_CAPABILITY_MODEL_UPDATE       = 1u << 1,
    PARTICLE_FIELD_CAPABILITY_EULERIAN_SCATTER   = 1u << 2,
    PARTICLE_FIELD_CAPABILITY_CHECKPOINT          = 1u << 3
} ParticleFieldCapabilities;

/** @brief Immutable metadata for one persistent particle field. */
typedef struct {
    ParticleFieldId           id;
    const char               *canonical_name;
    const char               *alias_1;
    const char               *alias_2;
    PetscInt                  components;
    PetscDataType             data_type;
    ParticleFieldRegistration registration;
    unsigned int              capabilities;
    PetscReal                 default_real_value;
    FieldId                   eulerian_scatter_target;
} ParticleFieldDescriptor;

/**
 * @brief Return immutable metadata for a valid particle field ID.
 * @param[in]  field_id Typed particle-field identity.
 * @param[out] descriptor Catalog-owned immutable descriptor.
 * @return Zero on success; a PETSc argument error for invalid input.
 */
PetscErrorCode ParticleFieldGetDescriptor(ParticleFieldId field_id,
                                          const ParticleFieldDescriptor **descriptor);

/**
 * @brief Resolve a canonical name or registered alias at a text-ingress boundary.
 * @param[in]  field_name Canonical particle name or registered alias.
 * @param[out] field_id Resolved typed identity.
 * @return Zero on success; `PETSC_ERR_ARG_UNKNOWN_TYPE` for unknown text.
 */
PetscErrorCode ParticleFieldIdFromName(const char *field_name, ParticleFieldId *field_id);

/**
 * @brief Return the canonical PETSc DMSwarm name for an ID.
 * @param[in] field_id Typed particle-field identity.
 * @return Catalog-owned canonical name or `InvalidParticleField`.
 */
const char *ParticleFieldName(ParticleFieldId field_id);

#ifdef __cplusplus
}
#endif

#endif /* PARTICLE_FIELD_CATALOG_H */
