/**
 * @file particle_field_catalog.c
 * @brief Persistent solver-particle-field catalog implementation.
 */

#include "particle_field_catalog.h"

#define PARTICLE_FIELD_ENTRY(field_id_, name_, alias1_, alias2_, components_, type_, registration_, capabilities_, default_, scatter_target_) \
    [field_id_] = {                                                                                                                         \
        field_id_, name_, alias1_, alias2_, components_, type_, registration_, capabilities_, default_, scatter_target_                    \
    }

static const ParticleFieldDescriptor gParticleFieldCatalog[PARTICLE_FIELD_ID_COUNT] = {
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_POSITION, "position", "ParticlePosition", NULL,
                         3, PETSC_REAL, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT, 0.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_VELOCITY, "velocity", "ParticleVelocity", NULL,
                         3, PETSC_REAL, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE |
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT, 0.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_CELL_ID, "DMSwarm_CellID", "CellID", NULL,
                         3, PETSC_INT, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT, 0.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_WEIGHT, "weight", "ParticleWeight", NULL,
                         3, PETSC_REAL, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE |
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT, 0.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_DIFFUSIVITY, "Diffusivity", NULL, NULL,
                         1, PETSC_REAL, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE, 1.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_DIFFUSIVITY_GRADIENT, "DiffusivityGradient", NULL, NULL,
                         3, PETSC_REAL, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE, 1.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_PSI, "Psi", NULL, NULL,
                         1, PETSC_REAL, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE |
                         PARTICLE_FIELD_CAPABILITY_MODEL_UPDATE |
                         PARTICLE_FIELD_CAPABILITY_EULERIAN_SCATTER |
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT,
                         0.0, FIELD_ID_PSI),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_LOCATION_STATUS, "DMSwarm_location_status", "Migration Status", NULL,
                         1, PETSC_INT, PARTICLE_FIELD_REGISTRATION_PICURV,
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT, 0.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_PID, "DMSwarm_pid", "pid", "Particle ID",
                         1, PETSC_INT64, PARTICLE_FIELD_REGISTRATION_PETSC,
                         PARTICLE_FIELD_CAPABILITY_CHECKPOINT, 0.0, FIELD_ID_INVALID),
    PARTICLE_FIELD_ENTRY(PARTICLE_FIELD_ID_RANK, "DMSwarm_rank", "rank", NULL,
                         1, PETSC_INT, PARTICLE_FIELD_REGISTRATION_PETSC,
                         PARTICLE_FIELD_CAPABILITY_NONE, 0.0, FIELD_ID_INVALID)
};

_Static_assert(sizeof(gParticleFieldCatalog) / sizeof(gParticleFieldCatalog[0]) == PARTICLE_FIELD_ID_COUNT,
               "Particle field catalog must contain one entry per ParticleFieldId.");

/**
 * @brief Validates and returns one immutable particle-field descriptor.
 * @see ParticleFieldGetDescriptor()
 */
PetscErrorCode ParticleFieldGetDescriptor(ParticleFieldId field_id,
                                          const ParticleFieldDescriptor **descriptor)
{
    PetscFunctionBeginUser;
    PetscCheck(descriptor != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Particle field descriptor output cannot be NULL.");
    PetscCheck(field_id >= 0 && field_id < PARTICLE_FIELD_ID_COUNT,
               PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Invalid ParticleFieldId value %d.", (int)field_id);
    PetscCheck(gParticleFieldCatalog[field_id].id == field_id,
               PETSC_COMM_SELF, PETSC_ERR_PLIB,
               "Particle field catalog entry %d is not initialized consistently.", (int)field_id);

    *descriptor = &gParticleFieldCatalog[field_id];
    PetscFunctionReturn(0);
}

/**
 * @brief Resolves particle field text once by searching canonical names and aliases.
 * @see ParticleFieldIdFromName()
 */
PetscErrorCode ParticleFieldIdFromName(const char *field_name, ParticleFieldId *field_id)
{
    PetscFunctionBeginUser;
    PetscCheck(field_name != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Particle field name cannot be NULL.");
    PetscCheck(field_id != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "ParticleFieldId output cannot be NULL.");

    for (PetscInt index = 0; index < PARTICLE_FIELD_ID_COUNT; ++index) {
        const ParticleFieldDescriptor *descriptor = &gParticleFieldCatalog[index];
        PetscBool match = PETSC_FALSE;

        PetscCall(PetscStrcasecmp(field_name, descriptor->canonical_name, &match));
        if (!match && descriptor->alias_1) PetscCall(PetscStrcasecmp(field_name, descriptor->alias_1, &match));
        if (!match && descriptor->alias_2) PetscCall(PetscStrcasecmp(field_name, descriptor->alias_2, &match));
        if (match) {
            *field_id = descriptor->id;
            PetscFunctionReturn(0);
        }
    }

    *field_id = PARTICLE_FIELD_ID_INVALID;
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Particle field name '%s' is not registered in the persistent particle field catalog.",
            field_name);
}

/**
 * @brief Returns the canonical PETSc name without performing name lookup.
 * @see ParticleFieldName()
 */
const char *ParticleFieldName(ParticleFieldId field_id)
{
    if (field_id < 0 || field_id >= PARTICLE_FIELD_ID_COUNT) return "InvalidParticleField";
    if (gParticleFieldCatalog[field_id].id != field_id) return "InvalidParticleField";
    return gParticleFieldCatalog[field_id].canonical_name;
}

#undef PARTICLE_FIELD_ENTRY
