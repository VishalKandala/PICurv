/**
 * @file walkingsearch.h
 * @brief Header file for particle location functions using the walking search algorithm.
 *
 * This file contains declarations of functions that implement the walking search algorithm
 * to locate particles within a computational grid. It includes functions for computing
 * distances to cell faces, determining particle positions relative to cells, and updating
 * traversal parameters.
 */

#ifndef WALKINGSEARCH_H
#define WALKINGSEARCH_H

// Include necessary headers
#include <petsc.h>
#include <stdbool.h>
#include <math.h>
#include "variables.h"   // Common type definitions
#include "logging.h"  // Logging macros and definitions
#include "setup.h"    //  utility functions 

// --------------------- Function Declarations ---------------------

/**
 * @brief Estimates a characteristic length of the cell for threshold scaling.
 *
 * For a hexahedral cell with vertices cell->vertices[0..7], we approximate
 * the cell size by some measure—e.g. average edge length or diagonal.
 * @param[in]  cell A pointer to the Cell structure
 * @param[out] cellSize A pointer to a PetscReal where the characteristic size is stored.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
 PetscErrorCode GetCellCharacteristicSize(const Cell *cell, PetscReal *cellSize);

/**
 * @brief Computes the signed distances from a point to each face of a cubic cell.
 *
 * This function calculates the signed distances from a specified point `p` to each of the six
 * faces of a cubic cell. The cell is defined by its eight vertices, and the distances are
 * stored in the memory location pointed to by `d`.
 *
 * @param[in]  p         The target point in 3D space.
 * @param[in]  cell      Pointer to a `Cell` structure that defines the cubic cell via its vertices.
 * @param[out] d         A pointer to an array of six `PetscReal` values where the computed signed distances will be stored.
 * @param[in]  threshold A `PetscReal` value that specifies the minimum distance below which a computed distance is considered to be zero.
 * @param[out] normals   Optional, may be NULL: an array of six `Cmpnts` receiving each face plane's unit outward normal, in the same order as `d`.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CalculateDistancesToCellFaces(const Cmpnts p, const Cell *cell, PetscReal *d, const PetscReal threshold, Cmpnts *normals);

/**
 * @brief Signed distances from a point to precomputed cell face planes.
 *
 * @param[in]  p         The target point.
 * @param[in]  planes    Face planes from MetricBuildCellFacePlanes().
 * @param[out] d         Six signed distances, positive inside, indexed by `Face`; below `threshold` set to 0.
 * @param[in]  threshold Absolute distance below which a distance is zero.
 * @param[out] normals   Optional, may be NULL: each face's unit outward normal, in the order of `d`.
 *
 * @return PetscErrorCode 0 on success, or PETSC_ERR_USER for a degenerate face; faces are checked
 *         in BACK, FRONT, BOTTOM, TOP, LEFT, RIGHT order and the first degenerate one is reported.
 */
PetscErrorCode CalculateDistancesToCellFacePlanes(const Cmpnts p, const CellFacePlanes *planes, PetscReal *d,
                                                  const PetscReal threshold, Cmpnts *normals);

/**
 * @brief Classifies a point based on precomputed face distances.
 *
 * Given an array of six distances d[NUM_FACES] from the point to each face
 * of a hexahedral cell, this function determines:
 *    0 => inside
 *    1 => on a single face
 *    2 => on an edge (2 faces = 0)
 *    3 => on a corner (3+ faces = 0)
 *   -1 => outside
 *
 * @param[in]  d        Array of six face distances.
 * @param[out] result   Pointer to an integer classification: {0,1,2,3,-1}.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
 PetscErrorCode DeterminePointPosition(PetscReal *d, PetscInt *result);

/**
 * @brief Initializes traversal parameters for locating a particle.
 *
 * This function sets the initial cell indices based on the local grid boundaries
 * and initializes the traversal step counter.
 *
 * @param[in]  user            Pointer to the user-defined context containing grid information.
 * @param[in]  particle        Pointer to the Particle structure containing its location.
 * @param[out] idx             Pointer to store the initial i-index of the cell.
 * @param[out] idy             Pointer to store the initial j-index of the cell.
 * @param[out] idz             Pointer to store the initial k-index of the cell.
 * @param[out] traversal_steps Pointer to store the initial traversal step count.
 *
 * @return PetscErrorCode     Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeTraversalParameters(UserCtx *user, Particle *particle, PetscInt *idx, PetscInt *idy, PetscInt *idz, PetscInt *traversal_steps);

/**
 * @brief Checks if the current cell indices are within the local grid boundaries.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 * @param[in]  idx        The i-index of the current cell.
 * @param[in]  idy        The j-index of the current cell.
 * @param[in]  idz        The k-index of the current cell.
 * @param[out] is_within  Pointer to a PetscBool that will be set to PETSC_TRUE if within bounds, else PETSC_FALSE.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within);

/**
 * @brief Retrieves the coordinates of the eight vertices of the current cell.
 *
 * Vertices 0 to 7 are the nodes (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k), (i,j+1,k+1),
 * (i+1,j+1,k+1), (i+1,j,k+1) and (i,j,k+1) of cell (i,j,k) = (`idx`,`idy`,`idz`), the order
 * MetricBuildCellFacePlanes() reads.
 *
 * @param[in]  user  Pointer to the user-defined context containing grid information.
 * @param[in]  idx   The i-index of the current cell.
 * @param[in]  idy   The j-index of the current cell.
 * @param[in]  idz   The k-index of the current cell.
 * @param[out] cell  Pointer to a Cell structure to store the cell's vertices.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RetrieveCurrentCell(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, Cell *cell);

/**
 * @brief Classifies a point against a cell's precomputed face planes.
 *
 * Computes the signed distances to the six face planes, zeroing those below
 * `threshold * planes->size`, then classifies them with DeterminePointPosition().
 *
 * @param[in]  planes    Face planes with `size` set, e.g. from MetricGatherCellFacePlanes().
 * @param[out] d         Six signed face distances, thresholded at `threshold * planes->size`.
 * @param[in]  p         The particle position.
 * @param[out] position  0 inside, 1 to 3 on a face, edge or corner, -1 outside.
 * @param[in]  threshold Base threshold, scaled by the cell size.
 *
 * @return PetscErrorCode 0 on success, or PETSC_ERR_USER for a degenerate face.
 */
PetscErrorCode EvaluateParticlePositionInPlanes(const CellFacePlanes *planes, PetscReal *d, const Cmpnts p,
                                                PetscInt *position, const PetscReal threshold);


/**
 * @brief Finds the MPI rank that owns a given global cell index.
 * @ingroup DomainInfo
 *
 * The DMDA decomposition is a tensor product of per-axis ranges, with rank
 * i + m (j + n k) on process grid m x n x p, so the owner is found by a binary search of
 * the range starts along each axis in the pre-computed decomposition map
 * (`user->RankCellInfoMap`): O(log m + log n + log p), not O(ranks). It is the definitive
 * method for resolving cell ownership in the "Walk and Handoff" migration algorithm.
 *
 * If the provided indices are outside the range of any rank (e.g., negative or
 * beyond the global domain), the function will not find an owner and `owner_rank`
 * will be set to -1.
 *
 * @param[in]  user       Pointer to the UserCtx structure, which must contain the
 *                        initialized `RankCellInfoMap` and `num_ranks`.
 * @param[in]  i          Global i-index of the cell to find.
 * @param[in]  j          Global j-index of the cell to find.
 * @param[in]  k          Global k-index of the cell to find.
 * @param[out] owner_rank Pointer to a `PetscMPIInt` where the resulting owner rank will
 *                        be stored. It is set to -1 if no owner is found.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode FindOwnerOfCell(UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscMPIInt *owner_rank);

/**
 * @brief Mirrors a point about the plane of one face of a cell.
 *
 * The plane is the one the walking search tested the point against, so a point outside the
 * face by |d| ends inside it by |d|: p' = p + 2 d n, with d < 0 outside. Only the position changes.
 *
 * @param[in]     planes The cell's face planes (MetricGatherCellFacePlanes() or MetricBuildCellFacePlanes()).
 * @param[in]     face The face index (LEFT, RIGHT, BOTTOM, TOP, FRONT or BACK).
 * @param[in,out] p    The point to mirror.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code (including a degenerate face).
 */
PetscErrorCode ReflectParticleAcrossCellFace(const CellFacePlanes *planes, PetscInt face, Cmpnts *p);

/**
 * @brief Moves a point to just inside a cell along the segment toward the cell centroid.
 *
 * The terminal action when reflections reach their cap. Plane distances are linear along the
 * segment, so the entry parameter is the largest d_p / (d_p - d_c) over the faces the point is
 * outside of; the point is placed a relative 1e-6 of the remaining segment beyond it.
 *
 * @param[in]     cell   The cell to place the point in.
 * @param[in,out] p      The point; unchanged when it cannot be placed.
 * @param[out]    placed PETSC_FALSE when the centroid is not inside every face plane the point
 *                       is outside of (a cell warped that far), in which case `p` is unchanged.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode PlaceParticleInsideCell(const Cell *cell, Cmpnts *p, PetscBool *placed);

/**
 * @brief Locates a particle's host cell or identifies its migration target using a robust walk search.
 * @ingroup ParticleLocation
 *
 * This is the core search engine. It starts from a guess cell and walks through
 * the grid. It returns a definitive, actionable status indicating the outcome:
 * - `ACTIVE_AND_LOCATED`: The particle was found in a cell on the current rank.
 * - `MIGRATING_OUT`: The particle was found to belong to another rank. `particle->destination_rank` is set.
 * - `EXITED`: The particle lies beyond a domain-boundary face of the boundary cell the walk reached,
 *   with every face it is outside of on the domain boundary and at least one of them using the
 *   `remove` particle policy. When all of them use `reflect`, the particle is instead mirrored
 *   about the most-crossed face (`ReflectParticleAcrossCellFace`) and the search continues; after
 *   `REFLECTION_CAP` reflections it is placed just inside the cell (`PlaceParticleInsideCell`).
 * - `NEEDS_LOCATION` (with `cell` reset to -1): a pass-1 walk from the previous-step cell ended more
 *   than `BLOCKED_WALK_GUARD_CELLS` cell sizes beyond the boundary, so the walk was most likely blocked
 *   by a non-convex boundary; the caller searches the particle again from scratch.
 * - `LOST`: The search failed for a particle believed to be inside the domain.
 *
 * This function is globally aware and can walk across MPI rank boundaries. It contains
 * robust checks for global domain boundaries and a tie-breaker for numerically "stuck"
 * particles on cell faces. A face lies on the domain boundary when the cell index is at the
 * extreme of that axis, the same test the boundary-condition handlers use.
 *
 * @param[in]     user         Pointer to the UserCtx containing all grid and domain info.
 * @param[in]     planes       The block's primal face-plane arrays, already open with
 *                             `MetricGetFacePlaneArrays()`, or NULL to open and close them here.
 *                             A caller searching many particles opens them once and passes them.
 * @param[in,out] particle     Pointer to the Particle struct. Its fields are updated based
 *                             on the search outcome.
 * @param[out]    status_out   The final, actionable status of the particle after the search.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *
 * @note Testing status:
 *       This function is exercised indirectly by the orchestrator tests, but
 *       direct branch pinning for boundary clamp, ghost-region handoff,
 *       tie-breaker, `LOST`, and `MIGRATING_OUT` outcomes is still part of the
 *       next-gap test backlog.
 */
struct FacePlaneArrays; /* Metric.h; not included here, since Metric.h includes setup.h, which includes this header */
PetscErrorCode LocateParticleOrFindMigrationTarget(UserCtx *user,
                                                        const struct FacePlaneArrays *planes,
                                                        Particle *particle,
                                                        ParticleLocationStatus *status_out);

/**
 * @brief Reports the outcome of one particle's location search: its consolidated `locate` line.
 *
 * Emitted through LOG_PARTICLE(), so it prints for a traced particle or at VERBOSE with this
 * function allow-listed. Lost and exited particles are counted in `search_metrics.csv`.
 *
 * @param simCtx          Simulation context (traced IDs, step).
 * @param particle        Particle whose location attempt has completed.
 * @param status          Final locate/migration status assigned to the particle.
 * @param start_cell      Cell the walk started from.
 * @param traversal_steps Number of cell-to-cell traversal steps attempted.
 * @param reflections     Number of wall reflections applied during the walk.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ReportSearchOutcome(const SimCtx *simCtx, const Particle *particle,
                                   ParticleLocationStatus status, const PetscInt start_cell[3],
                                   PetscInt traversal_steps, PetscInt reflections);

/**
 * @brief Updates the cell indices based on the signed distances to each face.
 *
 * This function modifies the cell indices (`idx`, `idy`, `idz`) to move towards the direction
 * where the particle is likely to be located, based on positive distances indicating
 * that the particle is outside in that particular direction.
 *
 * @param[in]  d    An array of six `PetscReal` values representing the signed distances to each face:
 *                  - d[LEFT]: Left Face
 *                  - d[RIGHT]: Right Face
 *                  - d[BOTTOM]: Bottom Face
 *                  - d[TOP]: Top Face
 *                  - d[FRONT]: Front Face
 *                  - d[BACK]: Back Face
 * @param[out] idx  Pointer to the i-index of the cell to be updated.
 * @param[out] idy  Pointer to the j-index of the cell to be updated.
 * @param[out] idz  Pointer to the k-index of the cell to be updated.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode UpdateCellIndicesBasedOnDistances( PetscReal d[NUM_FACES], PetscInt *idx, PetscInt *idy, PetscInt *idz);

#endif // WALKINGSEARCH_H
