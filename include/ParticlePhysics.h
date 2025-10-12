/**
 * @file ParticlePhysics.h
 * @brief Header file for Particle related physics modules.
 *
 * This file contains declarations of functions responsible for solving for fields in the particle swarms within a simulation using PETSc's DMSwarm.
 */

 #ifndef PARTICLE_PHYSICS_H
 #define PARTICLE_PHYSICS_H

// Include necessary headers
#include <petsc.h>        // PETSc library header
#include <petscdmswarm.h> // PETSc DMSwarm header
#include <stdbool.h>
#include <petscsys.h>     // For PetscRealloc
#include <math.h>
#include "variables.h"       // Common type definitions
#include "logging.h"      // Logging macros and definitions
#include "walkingsearch.h"  // Walking search function for particle migration 

/**
 * @brief Updates a single particle's field based on its state and the specified field name.
 *
 * This function serves as a switchboard for various particle property calculations.
 * Given a particle's data, it computes a new value for the specified field.
 *
 * Currently supported fields:
 *  - "Psi": Calculates the particle's kinetic energy (0.5 * |v|^2).
 *
 * @param[in]     fieldName   The name of the field to update (e.g., "Psi").
 * @param[in]     t           The current simulation time.
 * @param[in]     pos         The particle's physical position.
 * @param[in]     vel         The particle's velocity.
 * @param[in,out] psi_io      A pointer to the particle's current Psi value (used for accumulation if needed).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateParticleField(const char *fieldName,
                                                 PetscReal t,
                                                 Cmpnts pos,
                                                 Cmpnts vel,
                                                 PetscReal *psi_io);


/**
 * @brief Loops over all local particles and updates a specified field.
 *
 * This function orchestrates the update of a single particle field across the entire
 * local swarm. It gets access to the necessary particle data arrays and calls the
 * `UpdateParticleField` kernel for each particle.
 *
 * @param[in,out] user      Pointer to the UserCtx containing the swarm and simulation context.
 * @param[in]     fieldName The name of the field to update (e.g., "Psi").
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateFieldForAllParticles(UserCtx *user, const char *fieldName);


/**
 * @brief Orchestrates the update of all physical properties for particles.
 *
 * This function serves as the top-level entry point for updating particle-specific
 * physical quantities after their position and the surrounding fluid velocity are known.
 * It calls a sequence of more specific update routines for each property.
 *
 * For example, it can be configured to update:
 *  - Particle kinetic energy (stored in "Psi")
 *  - Particle accumulated shear stress
 *  - Particle temperature
 *  - etc.
 *
 * @param[in,out] user  Pointer to the UserCtx containing the swarm and simulation context.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateAllParticleFields(UserCtx *user);

#endif // PARTICLE_PHYSICS_H