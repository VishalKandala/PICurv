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
#include "field_catalog.h"     // Typed Eulerian field identities
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "interpolation.h"  // Interpolation routines
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
 * @brief (Public) Validates the consistency and compatibility of the parsed boundary condition system.
 *
 * This function is the main entry point for all boundary condition validation. It should be
 * called from the main setup sequence AFTER the configuration file has been parsed by
 * `ParseAllBoundaryConditions` but BEFORE any `BoundaryCondition` handler objects are created.
 *
 * It acts as a dispatcher, calling specialized private sub-validators for different complex
 * BC setups (like driven flow) to ensure the combination of `mathematical_type` and `handler_type`
 * across all six faces is physically and numerically valid. This provides a "fail-fast"
 * mechanism to prevent users from running improperly configured simulations.
 *
 * @param user The UserCtx for a single block, containing the populated `boundary_faces` configuration.
 * @return PetscErrorCode 0 on success, non-zero PETSc error code on failure.
 */
PetscErrorCode BoundarySystem_Validate(UserCtx *user);

/**
 * @brief (Private) Creates and configures a specific BoundaryCondition handler object.
 *
 * This function acts as a factory. Based on the requested handler_type, it allocates
 * a BoundaryCondition object and populates it with the correct set of function
 * pointers corresponding to that specific behavior.
 *
 * @param handler_type The specific handler to create (e.g., BC_HANDLER_WALL_NOSLIP).
 * @param[out] new_bc_ptr  A pointer to where the newly created BoundaryCondition
 *                         object's address will be stored.
 * @return PetscErrorCode 0 on success.
 */

PetscErrorCode BoundaryCondition_Create(BCHandlerType handler_type, BoundaryCondition **new_bc_ptr);

/**
 * @brief Initializes the entire boundary system.
 *
 * @param[in,out] user Finest-level block context receiving parsed face configuration.
 * @param bcs_filename Path to the generated boundary-condition definition file.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode BoundarySystem_Initialize(UserCtx *user, const char *bcs_filename);

/**
 * @brief Propagates boundary condition configuration from finest to all coarser multigrid levels.
 *
 * Coarser levels need BC type information for geometric operations (e.g., periodic corrections)
 * but do NOT need full handler objects since timestepping only occurs at the finest level.
 * This function copies the boundary_faces configuration down the hierarchy.
 *
 * @param simCtx The master SimCtx containing the multigrid hierarchy
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode PropagateBoundaryConfigToCoarserLevels(SimCtx *simCtx);

/**
 * @brief Executes one full boundary condition update cycle for a time step.
 *
 * @param[in,out] user Block context whose boundary handlers update target values and fluxes.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode BoundarySystem_ExecuteStep(UserCtx *user);

/**
 * @brief (Private) A lightweight execution engine that calls the UpdateUbcs() method on all relevant handlers.
 *
 * This function's sole purpose is to re-evaluate the target boundary values (`ubcs`) for
 * flow-dependent boundary conditions (e.g., Symmetry, Outlets) after the interior
 * velocity field has changed, such as after the projection step.
 *
 * It operates based on a "pull" model: it iterates through all boundary handlers and
 * executes their `UpdateUbcs` method only if the handler has provided one. This makes the
 * system extensible, as new flow-dependent handlers can be added without changing this
 * engine. Handlers for fixed boundary conditions (e.g., a wall with a constant velocity)
 * will have their `UpdateUbcs` pointer set to `NULL` and will be skipped automatically.
 *
 * @note This function is a critical part of the post-projection refresh. It intentionally
 *       does NOT modify `ucont` and does NOT perform flux balancing.
 *
 * @param user The main UserCtx struct.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode BoundarySystem_RefreshUbcs(UserCtx *user);

/**
 * @brief Cleans up and destroys all boundary system resources.
 *
 * @param[in,out] user Block context whose boundary handlers and temporary state are released.
 * @return PetscErrorCode 0 on success.
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
 * For example, if grid_layers=2 on face BC_FACE_NEG_X, the function will create particle lines at:
 * - y ~ 0*dy, y ~ 1*dy (parallel to the Z-axis, starting from the J=0 edge)
 * - y ~ y_max, y ~ y_max-dy (parallel to the Z-axis, starting from the J=max edge)
 * - z ~ 0*dz, z ~ 1*dz (parallel to the Y-axis, starting from the K=0 edge)
 * - z ~ z_max, z ~ z_max-dz (parallel to the Y-axis, starting from the K=max edge)
 * The particle's final position is set just inside the target cell face to ensure it is
 * correctly located. The total number of particles (simCtx->np) is distributed as evenly
 * as possible among all generated lines.
 * The function includes extensive validation to stop with an error if the requested grid
 * placement is geometrically impossible (e.g., in a 2D domain or if layers would overlap).
 * It also issues warnings for non-fatal but potentially unintended configurations.
 *
 * @param user Inlet-boundary context that defines the target face and grid layers.
 * @param info Local ownership and ghost-range information.
 * @param xs_gnode_rank Global xi node index at this rank's owned lower corner.
 * @param ys_gnode_rank Global eta node index at this rank's owned lower corner.
 * @param zs_gnode_rank Global zeta node index at this rank's owned lower corner.
 * @param IM_cells_global Global number of xi cells.
 * @param JM_cells_global Global number of eta cells.
 * @param KM_cells_global Global number of zeta cells.
 * @param particle_global_id Global particle ordinal used for deterministic placement.
 * @param[out] ci_metric_lnode_out Local xi metric-node index of the chosen cell.
 * @param[out] cj_metric_lnode_out Local eta metric-node index of the chosen cell.
 * @param[out] ck_metric_lnode_out Local zeta metric-node index of the chosen cell.
 * @param[out] xi_metric_logic_out Logical xi coordinate within the chosen cell.
 * @param[out] eta_metric_logic_out Logical eta coordinate within the chosen cell.
 * @param[out] zta_metric_logic_out Logical zeta coordinate within the chosen cell.
 * @param[out] placement_successful_out PETSC_TRUE when this rank owns a valid placement.
 * @return PetscErrorCode 0 on success.
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
 * @param xs_gnode_rank Local i-start node index (including ghosts) for this rank.
 * @param ys_gnode_rank Local j-start node index (including ghosts) for this rank.
 * @param zs_gnode_rank Local k-start node index (including ghosts) for this rank.
 * @param IM_nodes_global Global node count in i.
 * @param JM_nodes_global Global node count in j.
 * @param KM_nodes_global Global node count in k.
 * @param rand_logic_i_ptr RNG handle for sampling local logical xi.
 * @param rand_logic_j_ptr RNG handle for sampling local logical eta.
 * @param rand_logic_k_ptr RNG handle for sampling local logical zta.
 * @param[out] ci_metric_lnode_out Local i node index of selected cell origin.
 * @param[out] cj_metric_lnode_out Local j node index of selected cell origin.
 * @param[out] ck_metric_lnode_out Local k node index of selected cell origin.
 * @param[out] xi_metric_logic_out Logical xi coordinate in [0,1].
 * @param[out] eta_metric_logic_out Logical eta coordinate in [0,1].
 * @param[out] zta_metric_logic_out Logical zta coordinate in [0,1].
 * @return PetscErrorCode
 */
PetscErrorCode GetRandomCellAndLogicalCoordsOnInletFace(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank, // Local starting node index (with ghosts) of the rank's DA patch
    PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out);

/**
 * @brief Classification of one staggered momentum row (location + component).
 *
 * @see ClassifyMomentumRow() for the meaning of each member and for the
 *      single-source-of-truth contract these values participate in.
 */
typedef enum {
    MOM_ROW_PHYSICAL = 0,       /**< Independent unknown governed by the momentum equation. */
    MOM_ROW_FIXED_CONDITIONED,  /**< Strong Dirichlet row; the value comes from ApplyBoundaryConditions(). */
    MOM_ROW_FIXED_HOMOGENEOUS,  /**< Dummy/tangential row carrying no unknown at all. */
    MOM_ROW_PERIODIC_DUPLICATE  /**< Duplicate of a wrapped representative row (see @p ri, @p rj, @p rk). */
} MomentumRowType;

/**
 * @brief Single source of truth for "which staggered momentum rows are unknowns".
 *
 * Every consumer of the momentum system must agree on which rows the solver is
 * responsible for, and every consumer must derive that answer from this function
 * rather than restating the index arithmetic locally. Three independent
 * restatements previously disagreed, and the disagreement was silent: the
 * residual assembly skipped the periodic duplicate column at index 0 while
 * nothing zeroed it, so `ComputeTotalResidual()`'s BDF term accumulated there
 * without bound and the reported residual norm stopped describing the state.
 *
 * The classification depends only on @p user->info, the configured boundary
 * types, and the queried index; it reads no field data and performs no
 * communication, so it is safe to call inside assembly loops.
 *
 * Periodicity of an axis is taken from that axis's NEGATIVE face, matching
 * `ComputeRHS()` and `TransferPeriodicStaggeredFieldByDirection()`. A periodic
 * axis is expected to carry PERIODIC on both of its faces.
 *
 * Callers act on the classification differently, and both actions are correct:
 *   - residual/pseudo-time consumers (`EnforceRHSBoundaryConditions()`) zero
 *     every non-physical row, because the value there is imposed immediately
 *     afterwards by the boundary sweep or the periodic synchronisation;
 *   - the matrix-free Newton path substitutes an explicit equation instead
 *     (`F = X - U_conditioned`, `F = X`, `F = X_dup - X_rep`), because a zeroed
 *     row would leave a zero Jacobian row.
 *
 * @param[in]  user      Block context supplying `info` and `boundary_faces`.
 * @param[in]  i         Location index along xi.
 * @param[in]  j         Location index along eta.
 * @param[in]  k         Location index along zeta.
 * @param[in]  component Staggered component of the row (0 = xi, 1 = eta, 2 = zeta).
 * @param[out] ri        Representative xi index; equals @p i unless the row wraps.
 * @param[out] rj        Representative eta index; equals @p j unless the row wraps.
 * @param[out] rk        Representative zeta index; equals @p k unless the row wraps.
 * @return The row classification. Only #MOM_ROW_PHYSICAL denotes an unknown.
 */
MomentumRowType ClassifyMomentumRow(UserCtx *user, PetscInt i, PetscInt j, PetscInt k,
                                    PetscInt component, PetscInt *ri, PetscInt *rj, PetscInt *rk);

/**
 * @brief Reports whether a momentum row is masked out by the solid-cell field.
 *
 * `ComputeRHS()` zeroes the residual at a solid cell and at the staggered rows whose
 * downstream neighbour is solid, so those rows carry no equation. A solver that only
 * marches on the residual needs nothing more: a zero residual means no update. A solver
 * that assembles a matrix and solves `F(X) = 0` does, because such a row has a zero
 * Jacobian row and a zero column, leaving the unknown undetermined. Those rows must be
 * constrained like any other row that carries no unknown.
 *
 * The condition mirrors the residual's own masking exactly. If one changes, so must the
 * other, or the assembled operator stops matching the residual it preconditions.
 *
 * @param[in] nvert     Ghosted solid-cell field, or NULL to apply no masking.
 * @param[in] i         Location index along xi.
 * @param[in] j         Location index along eta.
 * @param[in] k         Location index along zeta.
 * @param[in] component Staggered component of the row (0 = xi, 1 = eta, 2 = zeta).
 * @return `PETSC_TRUE` when the row is masked and therefore carries no unknown.
 */
PetscBool MomentumRowIsSolidMasked(const PetscReal ***nvert, PetscInt i, PetscInt j,
                                   PetscInt k, PetscInt component);

/**
 * @brief Zeroes every momentum RHS row that does not carry an independent unknown.
 *
 * The set of such rows is not restated here: each owned location and component is
 * asked of ClassifyMomentumRow(), and anything other than #MOM_ROW_PHYSICAL is
 * zeroed. That covers, without enumerating them,
 *
 *   - strong Dirichlet rows on non-periodic faces, so the time-stepping scheme
 *     cannot alter the values `ApplyBoundaryConditions()` has just set;
 *   - dummy layers at the far index of every axis, which hold no unknown; and
 *   - periodic duplicate columns, whose value the next
 *     `SynchronizePeriodicStaggeredFields()` copies from the wrapped master.
 *
 * The last case is the one that must not be skipped. `ComputeRHS()` leaves the
 * transverse components of a periodic duplicate column untouched, so a row left
 * unzeroed here retains its previous contents while `ComputeTotalResidual()`
 * adds the BDF term on top of them on every call. The residual norm then grows
 * by |dU|/dt per evaluation regardless of the state, and no pseudo-time
 * iteration can reduce it.
 *
 * Call immediately after the RHS vector is fully assembled (spatial + temporal
 * terms) and before it is used in a time-stepping update.
 *
 * @param user The UserCtx for the specific block being computed.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode EnforceRHSBoundaryConditions(UserCtx *user); 

/**
 * @brief Synchronizes periodic endpoint cells for a list of cell-centered fields.
 *
 * The fields are first communicated from global to local storage. Each periodic
 * direction is then transferred in i-j-k order, with an intermediate ghost
 * refresh after every active direction so periodic edges and corners inherit the
 * values established by earlier directions. Only global duplicate planes in active
 * periodic directions are repaired; non-periodic directions are untouched. The
 * routine is a no-op, including no local refresh, when every direction is
 * nonperiodic. During active periodic synchronization it internally refreshes the
 * local vectors, but it is not a general replacement for `UpdateLocalGhosts()`.
 *
 * Supported fields are selected by
 * `FIELD_CAPABILITY_PERIODIC_CELL_SYNC` in the field catalog.
 *
 * @param user The main UserCtx struct.
 * @param num_fields The number of entries in `field_ids`.
 * @param field_ids The cell-centered fields to synchronize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SynchronizePeriodicCellFields(UserCtx *user, PetscInt num_fields, const FieldId field_ids[]);

/**
 * @brief Synchronizes persistent fields belonging to one face family.
 *
 * The function performs deterministic I/J/K directional passes with an
 * intermediate ghost refresh after each active periodic direction. It updates
 * persistent global seam/dummy values only; face-specific local stencil repair
 * remains a separate operation.
 *
 * @param user The main UserCtx struct.
 * @param face_direction Face family shared by every field (`'i'`, `'j'`, or `'k'`).
 * @param[in] num_fields Count of registered face fields.
 * @param field_ids Registered persistent face-field identities.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SynchronizePeriodicFaceFields(UserCtx *user, char face_direction, PetscInt num_fields, const FieldId field_ids[]);

/**
 * @brief Synchronizes persistent component-staggered vector fields.
 *
 * The function performs deterministic I/J/K endpoint transfers with an
 * intermediate ghost refresh after every active periodic direction. Currently
 * `Ucont` is the only registered component-staggered field.
 *
 * @param user The main UserCtx struct.
 * @param num_fields Number of entries in `field_ids`.
 * @param field_ids Registered component-staggered field identities.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SynchronizePeriodicStaggeredFields(UserCtx *user, PetscInt num_fields,
                                                  const FieldId field_ids[]);

/**
 * @brief Repairs the outer adjacent periodic ghosts used by QUICK cell stencils.
 *
 * The supplied local vectors must already contain a current PETSc periodic
 * ghost exchange. The vector and scalar fields are repaired two logical cells
 * across each active periodic seam so QUICK's `i-1/i+2` equivalents are valid.
 *
 * @param user Main block context containing periodic boundary metadata.
 * @param local_vector_field Ghosted three-component cell-centered field.
 * @param local_scalar_field Ghosted scalar cell-centered field.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode PreparePeriodicQuickStencilFields(UserCtx *user, Vec local_vector_field,
                                                 Vec local_scalar_field);

/**
 * @brief Synchronizes one local-only component-staggered periodic work field.
 *
 * This helper communicates locally computed owned entries, establishes the
 * normal-component periodic endpoint values, and communicates once more.
 *
 * @param user Main block context containing periodic boundary metadata.
 * @param local_field Ghosted local component-staggered vector.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SynchronizePeriodicLocalStaggeredField(UserCtx *user, Vec local_field);

/**
 * @brief (Orchestrator) Updates all metric-related fields in the local ghost cell regions for periodic boundaries.
 *
 * This function synchronizes cell-centered `Aj` and the persistent I/J/K metric
 * face families through the canonical MPI-safe synchronizers.
 *
 * @param user The main UserCtx struct.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyMetricsPeriodicBCs(UserCtx *user);

/**
 * @brief Applies periodic boundary conditions by copying data across domain boundaries for all relevant fields.
 *
 * This is the canonical periodic orchestrator for geometric consistency. It updates
 * `Ucat`, `P`, and `Nvert` through the generic cell synchronizer and updates
 * staggered `Ucont` through the component-staggered synchronizer.
 *
 * Future extension rule: add new periodic variables by extending the existing field
 * string dispatchers and invoking them from this orchestrator.
 *
 * @param user The main UserCtx struct.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyPeriodicBCs(UserCtx *user);

/**
 * @brief Updates the dummy cells (ghost nodes) on the faces of the local domain for NON-PERIODIC boundaries.
 *
 * This function's role is to apply a second-order extrapolation to set the ghost
 * cell values based on the boundary condition value (stored in `ubcs`) and the
 * first interior cell.
 *
 * NOTE: This function deliberately IGNORES periodic boundaries. It is part of a
 * larger workflow where `ApplyPeriodicBCs` handles periodic faces first.
 *
 * CRITICAL DETAIL: This function uses shrunken loop ranges (lxs, lxe, etc.) to
 * intentionally update only the flat part of the faces, avoiding the edges and

 * corners. The edges and corners are then handled separately by `UpdateCornerNodes`.
 * This precisely replicates the logic of the original FormBCS function.
 *
 * @param user The main UserCtx struct containing all necessary data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateDummyCells(UserCtx *user);

/**
 * @brief Updates the corner and edge ghost nodes of the local domain by averaging.
 *
 * This function should be called AFTER the face ghost nodes are finalized by both
 * `ApplyPeriodicBCs` and `UpdateDummyCells`. It resolves the values at shared
 * edges and corners by averaging the values of adjacent, previously-computed
 * ghost nodes.
 *
 * The logic is generic and works correctly regardless of the boundary types on
 * the adjacent faces (e.g., it will correctly average a periodic face neighbor
 * with a wall face neighbor).
 *
 * @param user The main UserCtx struct containing all necessary data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateCornerNodes(UserCtx *user);

/**
 * @brief Applies wall function modeling to near-wall velocities for all wall-type boundaries.
 *
 * This function implements log-law wall functions to model the near-wall velocity profile
 * without fully resolving the viscous sublayer. It is applicable to ALL wall-type boundaries
 * regardless of their specific boundary condition (no-slip, moving wall, slip, etc.), as
 * determined by the mathematical_type being WALL.
 *
 * MATHEMATICAL BACKGROUND:
 * Wall functions bridge the gap between the wall (y=0) and the first computational cell
 * center by using empirical log-law relationships:
 *   - Viscous sublayer (y+ < 11.81): u+ = y+
 *   - Log-law region (y+ > 11.81): u+ = (1/κ) * ln(E * y+)
 * where u+ = u/u_τ, y+ = y*u_τ/ν, κ = 0.41 (von Karman constant), E = exp(κB)
 *
 * IMPLEMENTATION DETAILS:
 * Unlike standard boundary conditions that set ghost cell values, wall functions:
 *   1. Read velocity from the SECOND interior cell (i±2, j±2, k±2)
 *   2. Compute wall shear stress using log-law
 *   3. Modify velocity at the FIRST interior cell (i±1, j±1, k±1) 
 *   4. Keep ghost cell boundary values (ubcs, ucont) at zero
 *
 * WORKFLOW:
 *   - Called from ApplyBoundaryConditions after standard BC application
 *   - Operates on ucat (Cartesian velocity)
 *   - Updates ustar (friction velocity field) for diagnostics/turbulence models
 *   - Ghost cells remain zero; UpdateDummyCells handles extrapolation afterward
 *
 * GEOMETRIC QUANTITIES:
 *   sb = wall-normal distance from wall to first interior cell center
 *   sc = wall-normal distance from wall to second interior cell center  
 *   These are computed from cell Jacobians (aj) and face area vectors
 *
 * APPLICABILITY:
 *   - Requires simCtx->wallfunction = true
 *   - Only processes faces where mathematical_type == WALL
 *   - Skips solid-embedded cells (nvert >= 0.1)
 *
 * @param user The UserCtx containing all simulation state and geometry
 * @return PetscErrorCode 0 on success
 *
 * @note This function modifies interior cell velocities, NOT ghost cells
 * @note Wall roughness (ks) is currently set to 1e-16 (smooth wall)
 * @see wall_function_loglaw() in wallfunction.c for the actual log-law implementation
 * @see noslip() in wallfunction.c for the initial linear interpolation
 */
PetscErrorCode ApplyWallFunction(UserCtx *user);

/**
 * @brief Finalizes cell-centered fields after the projection step.
 *
 * This function completes the cell-centered state derived from the final,
 * divergence-free `Ucont` produced by `Projection`. It fills non-periodic
 * `Ucat` dummy faces, synchronizes periodic `Ucat` and `P` endpoints, resolves
 * edges and corners, and refreshes the corresponding local vectors.
 *
 * This function is fundamentally different from `ApplyBoundaryConditions`: it
 * does NOT modify `Ucont`, reapply wall functions, or rerun the full physical
 * boundary-condition workflow.
 *
 * WORKFLOW:
 * 1. Refreshes local `Ucat` and any flow-dependent `Ubcs` targets.
 * 2. Fills non-periodic dummy faces and establishes periodic cell endpoints.
 * 3. Resolves edges/corners, restores exact periodic relationships, and refreshes
 *    local `Ucat` and `P`.
 *
 * @param user The main UserCtx struct, containing all simulation state.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode FinalizePostProjectionCellFields(UserCtx *user);

/**
 * @brief Main boundary-condition orchestrator executed during solver timestepping.
 *
 * This routine performs the full BC workflow for the current block, including
 * dynamic boundary refresh, periodic transfer, dummy/corner updates, and optional
 * wall-function corrections in the same order expected by the runtime solver.
 * It may iterate boundary updates to enforce coupled boundary dependencies.
 *
 * @param user The main UserCtx struct containing field vectors and boundary system state.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyBoundaryConditions(UserCtx *user);

#endif // BOUNDARIES_H
