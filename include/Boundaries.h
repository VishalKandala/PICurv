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
 * @param IM_nodes_global Global number of nodes in the I-direction (e.g., user->IM + 1 if user->IM is cell count).
 * @param JM_nodes_global Global number of nodes in the J-direction.
 * @param KM_nodes_global Global number of nodes in the K-direction.
 * @param face_id           The specific global face (e.g., BC_FACE_NEG_Z) to check.
 * @param[out] can_service_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                           services the face, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CanRankServiceFace(const DMDALocalInfo *info, PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
                                  BCFace face_id, PetscBool *can_service_out);

/**
 * @brief Places particles in a deterministic grid/raster pattern on a specified domain face.
 *
 * This function creates a set of equidistant, parallel lines of particles near the four
 * edges of the face specified by user->identifiedInletBCFace. The number of lines drawn
 * from each edge is hardcoded within this function (default is 2).
 *
 * For example, if grid_layers=2 on face BC_FACE_NEG_X, the function will create particle lines at:
 * - y ~ 0*dy, y ~ 1*dy (parallel to the Z-axis, starting from the J=0 edge)
 * - y ~ y_max, y ~ y_max-dy (parallel to the Z-axis, starting from the J=max edge)
 * - z ~ 0*dz, z ~ 1*dz (parallel to the Y-axis, starting from the K=0 edge)
 * - z ~ z_max, z ~ z_max-dz (parallel to the Y-axis, starting from the K=max edge)
 *
 * The particle's final position is set just inside the target cell face to ensure it is
 * correctly located. The total number of particles (simCtx->np) is distributed as evenly
 * as possible among all generated lines.
 *
 * The function includes extensive validation to stop with an error if the requested grid
 * placement is geometrically impossible (e.g., in a 2D domain or if layers would overlap).
 * It also issues warnings for non-fatal but potentially unintended configurations.
 *
 * @param user Pointer to UserCtx, which must contain a valid identifiedInletBCFace.
 * @param info Pointer to DMDALocalInfo for the current rank's grid layout.
 * @param xs_gnode_rank, ys_gnode_rank, zs_gnode_rank Local starting node indices (incl. ghosts) for the rank's DA.
 * @param IM_cells_global, JM_cells_global, KM_cells_global Global cell counts.
 * @param particle_global_id The unique global ID of the particle being placed (from 0 to np-1).
 * @param[out] ci_metric_lnode_out Local I-node index of the selected cell's origin.
 * @param[out] cj_metric_lnode_out Local J-node index of the selected cell's origin.
 * @param[out] ck_metric_lnode_out Local K-node index of the selected cell's origin.
 * @param[out] xi_metric_logic_out Logical xi-coordinate [0,1] within the cell.
 * @param[out] eta_metric_logic_out Logical eta-coordinate [0,1] within the cell.
 * @param[out] zta_metric_logic_out Logical zta-coordinate [0,1] within the cell.
 * @param[out] placement_successful_out PETSC_TRUE if the point belongs to this rank, PETSC_FALSE otherwise.
 * @return PetscErrorCode
 */
PetscErrorCode GetDeterministicFaceGridLocation(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank,
    PetscInt IM_cells_global, PetscInt JM_cells_global, PetscInt KM_cells_global,
    PetscInt64 particle_global_id,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out,
    PetscBool *placement_successful_out);


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
PetscErrorCode GetRandomCellAndLogicalCoordsOnInletFace(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank, // Local starting node index (with ghosts) of the rank's DA patch
    PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out);

PetscErrorCode TranslateModernBCsToLegacy(UserCtx *user);

/**
 * @brief Applies inlet boundary conditions based on the modern BC handling system.
 *
 * This function iterates through all 6 domain faces. For each face identified as an
 * INLET, it applies the velocity profile specified by its assigned handler and
 * parameters (e.g., 'constant_velocity' with vx,vy,vz or 'parabolic' with u_max).
 *
 * It calculates the contravariant flux (Ucont), Cartesian velocity on the face (Ubcs),
 * and the staggered Cartesian velocity (Ucat). It also computes the total incoming
 * flux and area across all MPI ranks.
 *
 * @param user The main UserCtx struct containing the BC configuration and PETSc objects.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InflowFlux(UserCtx *user); 

/**
 * @brief Calculates the total outgoing flux through all OUTLET faces for reporting.
 *
 * NOTE: In a mixed modern/legacy environment, this function is for DIAGNOSTICS ONLY.
 * It reads the contravariant velocities and calculates the total flux passing through
 * faces marked as OUTLET. It does NOT apply any boundary conditions itself, as that
 * is still the responsibility of the legacy FormBCS function.
 *
 * @param user The main UserCtx struct containing BC config and PETSc vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode OutflowFlux(UserCtx *user);

// TO BE FIXED
PetscErrorCode FormBCS(UserCtx *user);

#endif // BOUNDARIES_H
