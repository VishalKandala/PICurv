/**
 * @file ParticleSwarm.h
 * @brief Header file for Particle Swarm management functions.
 *
 * This file contains declarations of functions responsible for creating, managing,
 * initializing, migrating, and printing particle swarms within a simulation using PETSc's DMSwarm.
 */

#ifndef PARTICLE_SWARM_H
#define PARTICLE_SWARM_H

// Include necessary headers
#include <petsc.h>        // PETSc library header
#include <petscdmswarm.h> // PETSc DMSwarm header
#include <stdbool.h>
#include <math.h>
#include "variables.h"       // Common type definitions
#include "particle_field_catalog.h"
#include "logging.h"      // Logging macros and definitions
#include "walkingsearch.h"
#include "Metric.h"
#include "io.h"
// --------------------- Function Declarations ---------------------

/**
 * @brief Creates the particle swarm, registers its fields and sizes it for this rank's share.
 *
 * Volume and Surface seeding give each rank the particles its owned cells seed (by volume, or
 * by area on the inlet face, the stratified split of the seeding table); PointSource divides
 * the count evenly, the first ranks taking one more. Collective.
 *
 * @param[in,out] user          Finest-level context; receives the swarm.
 * @param[in]     numParticles  Total number of particles to create across all MPI processes.
 * @param[out]    particlesPerProcess Number of particles assigned to the local rank.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles, PetscInt *particlesPerProcess);

/**
 * @brief Initializes the DMSwarm object within the UserCtx structure.
 *
 * This function creates the DMSwarm, sets its type and dimension, and configures basic swarm properties.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSwarm(UserCtx* user);

/**
 * @brief Registers a swarm field without finalizing registration.
 *
 * This function calls DMSwarmRegisterPetscDatatypeField for the given field,
 * but does not finalize the registration. The finalization is deferred until
 * all fields have been registered.
 *
 * @param swarm      [in]  The DMSwarm object.
 * @param fieldName  [in]  Name of the field to register.
 * @param fieldDim   [in]  Dimension of the field (1 for scalar, 3 for vector, etc.).
 * @param dtype      [in]  The datatype of the swarm field being registered.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RegisterSwarmField(DM swarm, const char *fieldName, PetscInt fieldDim, PetscDataType dtype);

/**
 * @brief Registers necessary particle fields within the DMSwarm.
 *
 * This function registers fields such as position, velocity, CellID, and weight for each particle.
 *
 * @param[in,out] swarm   The DMSwarm object managing the particle swarm.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RegisterParticleFields(DM swarm);

/**
 * @brief Initializes all particle properties in the swarm.
 *
 * Surface seeding requires an INLET face. The basic properties (physical position, particle
 * ID and the seeded cell) are set by `InitializeParticleBasicProperties`: Volume and Surface
 * seeding stratify the particles over the cells' volume or inlet area and key every placement
 * draw on (seed, particle ID), so the population does not depend on the rank layout. Every
 * field the particle catalog marks for default initialization (velocity, weight, Diffusivity,
 * DiffusivityGradient and Psi) then takes its catalog default. Collective.
 *
 * @param[in,out] user               Pointer to the `UserCtx` structure.
 * @param[in]     particlesPerProcess Number of particles assigned to this MPI process.
 *
 * @return PetscErrorCode            Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialPropertiesToSwarm(UserCtx* user,
                                              PetscInt particlesPerProcess);

/**
 * @brief Distributes particles evenly across MPI processes, handling any remainders.
 *
 * This function calculates the number of particles each MPI process should handle,
 * distributing the remainder particles to the first few ranks if necessary.
 *
 * @param[in]     numParticles        Total number of particles to create across all MPI processes.
 * @param[in]     rank                MPI rank of the current process.
 * @param[in]     size                Total number of MPI processes.
 * @param[out]    particlesPerProcess Number of particles assigned to the current MPI process.
 * @param[out]    remainder           Remainder particles when dividing numParticles by size.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DistributeParticles(PetscInt numParticles, PetscMPIInt rank, PetscMPIInt size, PetscInt* particlesPerProcess, PetscInt* remainder);

/**
 * @brief Initializes a Particle struct with data from DMSwarm fields.
 *
 * This helper function populates a Particle structure using data retrieved from DMSwarm fields.
 * `PIDs`, `positions`, `cellIndices` and `LocStatus` are required. Pass NULL for an optional
 * field the caller does not need; it then takes a default: weights (1, 1, 1), velocity and
 * diffusivity gradient zero, diffusivity 1, psi 0. It logs nothing: the stage that uses a
 * particle reports its values.
 *
 * @param[in]     i            Index of the particle in the DMSwarm.
 * @param[in]     PIDs         Pointer to the array of particle IDs.
 * @param[in]     weights      Pointer to the array of particle weights.
 * @param[in]     positions    Pointer to the array of particle positions.
 * @param[in]     cellIndices  Pointer to the array of particle cell indices.
 * @param[in]     velocities   Pointer to the array of particle velocities.
 * @param[in]     LocStatus    Pointer to the array of cell location status indicators.
 * @param[in]     diffusivity  Pointer to the array of particle diffusivities.
 * @param[in]     diffusivitygradient Pointer to the array of particle diffusivity gradients.
 * @param[in]     psi          Pointer to the array of particle psi values.
 * @param[out]    particle     Pointer to the Particle struct to initialize.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode UnpackSwarmFields(PetscInt i, const PetscInt64 *PIDs, const PetscReal *weights,
                const PetscReal *positions, const PetscInt *cellIndices,
                PetscReal *velocities,PetscInt *LocStatus,PetscReal *diffusivity, Cmpnts *diffusivitygradient, PetscReal *psi, Particle *particle);

/**
 * @brief Updates DMSwarm data arrays from a Particle struct.
 *
 * This function writes data from the `Particle` struct back into the raw DMSwarm arrays.
 * It is robust: if any array pointer is NULL, that specific field is skipped.
 * This allows selective updating (e.g., update position but not velocity).
 *
 * @param[in]     i            Index of the particle in the local swarm arrays.
 * @param[in]     particle     Pointer to the Particle struct containing updated data.
 * @param[in,out] positions    (Optional) Array of particle positions (size 3*n).
 * @param[in,out] velocities   (Optional) Array of particle velocities (size 3*n).
 * @param[in,out] weights      (Optional) Array of particle weights (size 3*n).
 * @param[in,out] cellIndices  (Optional) Array of particle cell indices (size 3*n).
 * @param[in,out] status       (Optional) Array of location status (size 1*n).
 * @param[in,out] diffusivity  (Optional) Array of diffusivity values (size 1*n).
 * @param[in,out] diffusivitygradient (Optional) Array of diffusivity gradient values (size 3*n).
 * @param[in,out] psi          (Optional) Array of scalar Psi values (size 1*n).
 *
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode UpdateSwarmFields(PetscInt i, const Particle *particle,
                                 PetscReal *positions, 
                                 PetscReal *velocities,
                                 PetscReal *weights, 
                                 PetscInt  *cellIndices, 
                                 PetscInt  *status,
                                 PetscReal *diffusivity,
                                 Cmpnts *diffusivitygradient,
                                 PetscReal *psi);                
              
/**
 * @brief Updates a particle's interpolation weights based on distances to cell faces.
 *
 * This function computes interpolation weights using distances to the six
 * cell faces (`d`) and updates the `weight` field of the provided particle.
 *
 * @param[in]  d        Pointer to an array of distances to the six cell faces.
 * @param[out] particle Pointer to the Particle structure whose weights are to be updated.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode UpdateParticleWeights(PetscReal *d, Particle *particle);

/**
 * @brief High-level particle initialization orchestrator for a simulation run.
 *
 * This routine drives end-to-end swarm setup from the top-level simulation context:
 * creation/registration of particle fields, initial placement according to configured
 * mode, initial localization on the Eulerian grid, and startup interpolation needed
 * before entering the main run loop.
 *
 * @param[in,out] simCtx Master simulation context containing all blocks and run settings.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeParticleSwarm(SimCtx *simCtx);

/**
 * @brief Read the inlets that inject particles from the options database.
 * @details Reads `-particle_inlet_count`, and for each inlet `n`, `-particle_inlet_<n>_face`,
 *          `-particle_inlet_<n>_particles_per_cell` (the mean per cell of the grid, an expression in `t`) and the
 *          field plan under `particle_inlet_<n>_fields`.
 * @param[in,out] simCtx Context receiving the inlets; none when the count is absent or zero.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ParticleInletsCreate(SimCtx *simCtx);

/**
 * @brief Release what ParticleInletsCreate() and InjectParticlesAtInlets() allocated.
 * @param[in,out] simCtx Context holding the inlets.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ParticleInletsDestroy(SimCtx *simCtx);

/**
 * @brief Add the particles that entered through each injecting inlet during this step.
 * @details Collective. A face cell receives on average `n (u.A) dt` particles, none where the flow
 *          leaves, so the entering fluid has number density `n` whatever its velocity profile.
 *          Counts, positions and identities come from keyed draws on the grid position, the step
 *          and the particle, so they carry no state and do not depend on the rank layout. A
 *          particle is placed just inside the face, by area and inward velocity, and its field
 *          values are set there, where it crossed. It then moves in by a uniform share of `dt`
 *          times the velocity the configured interpolation method gives at that point (the face
 *          cell's mean inward velocity where the cell has none) plus this step's Brownian kick
 *          from the face cell's diffusivity. The new particles are left `NEEDS_LOCATION` for the
 *          settlement that follows. Writes `particle_injection.csv`.
 * @param[in,out] user Finest-level context of the block holding the swarm; reads `lUcont`, `lCoor`,
 *                     the face-normal metric, `lDiffusivity`, and the interpolated velocity's
 *                     sources (`lUcat`, or `lUcont` and the metrics for FaceFlux).
 * @return PetscErrorCode 0 on success; an error on a face that is not an INLET, an invalid number
 *         density, or a particle-ID or particle-count limit being reached.
 */
PetscErrorCode InjectParticlesAtInlets(UserCtx *user);

#endif // PARTICLE_SWARM_H
