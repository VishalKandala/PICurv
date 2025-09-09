#ifndef POSTPROCESSING_KERNELS_H
#define POSTPROCESSING_KERNELS_H

#include "variables.h"
#include "logging.h"
#include "io.h" // For the UpdateLocalGhosts function prototype

// Function prototypes for post-processing kernels
PetscErrorCode ComputeNodalAverage(UserCtx* user, const char* in_field_name, const char* out_field_name);
PetscErrorCode ComputeQCriterion(UserCtx* user);
PetscErrorCode NormalizeRelativeField(UserCtx* user, const char* relative_field_name);

// Add more post-processing kernel prototypes as needed

// ===========================================================================
// Particle Post-Processing Kernels
// ===========================================================================

/**
 * @brief Computes the specific kinetic energy (KE per unit mass) for each particle.
 *
 * This kernel calculates SKE = 0.5 * |velocity|^2. It requires that the
 * velocity field exists and will populate the specific kinetic energy field.
 * The output field must be registered before this kernel is called.
 *
 * @param user           The UserCtx containing the DMSwarm.
 * @param velocity_field The name of the input vector field for particle velocity.
 * @param ske_field      The name of the output scalar field to store specific KE.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeSpecificKE(UserCtx* user, const char* velocity_field, const char* ske_field);

#endif // POSTPROCESSING_KERNELS_H