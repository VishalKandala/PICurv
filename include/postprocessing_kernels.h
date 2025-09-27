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
//  =========================================================================
//  Dimensionalization Kernels
//  =========================================================================
/**
 * @brief Scales a specified field from non-dimensional to dimensional units in-place.
 *
 * This function acts as a dispatcher. It takes the string name of a field,
 * identifies the corresponding PETSc Vec object and the correct physical
 * scaling factor (e.g., U_ref for velocity, P_ref for pressure), and then
 * performs an in-place VecScale operation. It correctly handles the different
 * physical dimensions of Cartesian velocity vs. contravariant volume flux.
 *
 * @param[in,out] user        The UserCtx containing the PETSc Vecs to be modified.
 * @param[in]     field_name  The case-insensitive string name of the field to dimensionalize
 *                            (e.g., "Ucat", "P", "Ucont", "Coordinates", "ParticlePosition", "ParticleVelocity").
 * @return PetscErrorCode
 */
PetscErrorCode DimensionalizeField(UserCtx *user, const char *field_name);

/**
 * @brief Orchestrates the dimensionalization of all relevant fields loaded from a file.
 *
 * This function is intended to be called in the post-processor immediately after
 * all solver output has been read into memory. It calls DimensionalizeField() for each of the core
 * physical quantities to convert the entire loaded state from non-dimensional to
 * dimensional units, preparing it for analysis and visualization.
 *
 * @param[in,out] user The UserCtx containing all the fields to be dimensionalized.
 * @return PetscErrorCode
 */
PetscErrorCode DimensionalizeAllLoadedFields(UserCtx *user);

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