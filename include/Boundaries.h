#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <petscsystypes.h>

// Include additional headers
#include "variables.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "interpolation.h"  // Interpolation routines
#include "AnalyticalSolution.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "BC_Handlers.h"    // Boundary Handlers 
#include "wallfunction.h"   //  wall functions for LES
//================================================================================
//
//                        PUBLIC SYSTEM-LEVEL FUNCTIONS
//
// These are the main entry points for interacting with the boundary system.
//
//================================================================================

/**
 * @brief Initializes the entire boundary system.
 * @param user The main UserCtx struct.
 * @param bcs_filename The path to the boundary conditions configuration file.
 */
PetscErrorCode BoundarySystem_Create(UserCtx *user, const char *bcs_filename);

/**
 * @brief Executes one full boundary condition update cycle for a time step.
 * @param user The main UserCtx struct.
 */
PetscErrorCode BoundarySystem_ExecuteStep(UserCtx *user);

/**
 * @brief Cleans up and destroys all boundary system resources.
 * @param user The main UserCtx struct.
 */
PetscErrorCode BoundarySystem_Destroy(UserCtx *user);

/**
 * @brief Determines if the current MPI rank owns any part of the globally defined inlet face,
 *        making it responsible for placing particles on that portion of the surface.
 *
 * The determination is based on the rank's owned nodes (from `DMDALocalInfo`) and
 * the global node counts, in conjunction with the `user->identifiedInletBCFace`.
 * A rank can service an inlet face if it owns the cells adjacent to that global boundary
 * and has a non-zero extent (owns cells) in the tangential dimensions of that face.
 *
 * @param user Pointer to the UserCtx structure, containing `identifiedInletBCFace`.
 * @param info Pointer to the DMDALocalInfo for the current rank's DA (node-based).
 * @param IM_nodes_global Global number of nodes in the I-direction (e.g., user->IM + 1 if user->IM is cell count).
 * @param JM_nodes_global Global number of nodes in the J-direction.
 * @param KM_nodes_global Global number of nodes in the K-direction.
 * @param[out] can_service_inlet_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                                   services (part of) the inlet, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CanRankServiceInletFace(UserCtx *user, const DMDALocalInfo *info,
                                              PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
                                              PetscBool *can_service_inlet_out);

/**
 * @brief Determines if the current MPI rank owns any part of a specified global face.
 *
 * This function is a general utility for parallel boundary operations. It checks if the
 * local domain of the current MPI rank is adjacent to a specified global boundary face.
 * A rank "services" a face if it owns the cells adjacent to that global boundary and has
 * a non-zero extent (i.e., owns at least one cell) in the tangential dimensions of that face.
 *
 * @param info              Pointer to the DMDALocalInfo for the current rank's DA.
 * @param face_id           The specific global face (e.g., BC_FACE_NEG_Z) to check.
 * @param[out] can_service_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                           services the face, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CanRankServiceFace(const DMDALocalInfo *info, BCFace face_id, PetscBool *can_service_out);


/**
 * @brief Assuming the current rank services the inlet face, this function selects a random
 *        cell (owned by this rank on that face) and random logical coordinates within that cell,
 *        suitable for placing a particle on the inlet surface.
 *
 * It is the caller's responsibility to ensure CanRankServiceInletFace returned true.
 *
 * @param user Pointer to UserCtx.
 * @param info Pointer to DMDALocalInfo for the current rank (node-based).
 * @param xs_gnode, ys_gnode, zs_gnode Local starting node indices (incl. ghosts) for the rank's DA.
 * @param IM_nodes_global, JM_nodes_global, KM_nodes_global Global node counts.
 * @param rand_logic_i_ptr, rand_logic_j_ptr, rand_logic_k_ptr Pointers to RNGs for logical coords.
 * @param[out] ci_metric_lnode_out, cj_metric_lnode_out, ck_metric_lnode_out Local node indices of the selected cell's origin (these are local to the rank's DA including ghosts).
 * @param[out] xi_metric_logic_out, eta_metric_logic_out, zta_metric_logic_out Logical coords [0,1] within the cell.
 * @return PetscErrorCode
 */
PetscErrorCode GetRandomCellAndLogicOnInletFace(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank, // Local starting node index (with ghosts) of the rank's DA patch
    PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out);

PetscErrorCode TranslateModernBCsToLegacy(UserCtx *user);

/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation/interface
   -1:  wallfunction
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   3:   slip wall/symmetry
   5:	Inlet
   4:	Outlet
   6:   farfield
   7:   periodic
   8:   Characteristic BC
   9:   Analytical Vortex
   10:  Oulet Junction
   11:  Annulus
   12:  Ogrid
   13:  Rheology
   14:  Outlet with Interface
  
*/

PetscErrorCode InflowFlux(UserCtx *user);

PetscErrorCode OutflowFlux(UserCtx *user);

PetscErrorCode FormBCS(UserCtx *user);

/**
 * @brief Acts as a temporary bridge to the legacy boundary condition implementation.
 *
 * This function is the key to our integration strategy. It matches the signature
 * of the modern `BoundarySystem_ExecuteStep` function that SetEulerianFields
 * expects to call.
 *
 * However, instead of containing new handler-based logic, it simply calls the
 * monolithic legacy `FormBCS` function. This allows the modern orchestrator to
 * drive the old solver logic without modification.
 *
 * @param user The UserCtx for a single block.
 * @return PetscErrorCode
 */
PetscErrorCode BoundarySystem_ExecuteStep_Legacy(UserCtx *user);

#endif // BOUNDARIES_H
