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
 * @brief Computes the signed distance from a point to the plane approximating a quadrilateral face.
 *
 * This function calculates the signed distance from a given point `p_target` to the plane
 * approximating a quadrilateral face defined by points `v1`, `v2`, `v3`, and `v4`.
 * The plane's initial normal is determined using two edge vectors emanating from `v1`
 * (i.e., `v2-v1` and `v4-v1`).
 *
 * **Normal Orientation:**
 * The initial normal's orientation is checked against the cell's interior.
 * A vector from the face centroid to the cell centroid (`vec_face_to_cell_centroid`) is computed.
 * If the dot product of the initial normal and `vec_face_to_cell_centroid` is positive,
 * it means the initial normal is pointing towards the cell's interior (it's an inward normal
 * relative to the cell). In this case, the normal is flipped to ensure it points outward.
 *
 * **Signed Distance Convention:**
 * The signed distance is the projection of the vector from `p_target` to the face's centroid
 * onto the (now guaranteed) outward-pointing plane's normal vector.
 *   - `d_signed < 0`: `p_target` is "outside" and "beyond" the face, in the direction of the outward normal.
 *                     This is the primary case indicating the particle has crossed this face.
 *   - `d_signed > 0`: `p_target` is "outside" but "behind" the face (on the same side as the cell interior
 *                     relative to this face plane).
 *   - `d_signed = 0`: `p_target` is on the face plane (within threshold).
 *
 * @param[in]  v1            First vertex defining the face (used as origin for initial normal calculation).
 * @param[in]  v2            Second vertex defining the face (used for first edge vector: v2-v1).
 * @param[in]  v3            Third vertex defining the face (used for face centroid calculation).
 * @param[in]  v4            Fourth vertex defining the face (used for second edge vector: v4-v1).
 * @param[in]  cell_centroid The centroid of the 8-vertex cell to which this face belongs.
 * @param[in]  p_target      The point from which the distance to the plane is calculated.
 * @param[out] d_signed      Pointer to store the computed signed distance.
 * @param[in]  threshold     The threshold below which the absolute distance is considered zero.
 *
 * @note
 * - The order of vertices `v1, v2, v3, v4` is used for calculating the face centroid.
 *   The order of `v1, v2, v4` is used for the initial normal calculation.
 * - The `cell_centroid` is crucial for robustly determining the outward direction of the normal.
 *
 * @return PetscErrorCode Returns 0 on success, PETSC_ERR_ARG_NULL if `d_signed` is NULL,
 *                        or PETSC_ERR_USER if a degenerate plane (near-zero normal) is detected.
 */
PetscErrorCode ComputeSignedDistanceToPlane(const Cmpnts v1, const Cmpnts v2, const Cmpnts v3, const Cmpnts v4,
                                            const Cmpnts cell_centroid,const Cmpnts p_target,
					    PetscReal *d_signed, const PetscReal threshold);

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
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CalculateDistancesToCellFaces(const Cmpnts p, const Cell *cell, PetscReal *d, const PetscReal threshold);

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
 * @brief Retrieves the coordinates of the eight vertices of a cell based on grid indices.
 *
 * This function populates the `cell` structure with the coordinates of the eight vertices
 * of a hexahedral cell given its grid indices `(idx, idy, idz)`.
 *
 * @param[in]  coor   Pointer to a 3D array of `Cmpnts` structures representing the grid coordinates.
 * @param[in]  idx    The i-index of the cell.
 * @param[in]  idy    The j-index of the cell.
 * @param[in]  idz    The k-index of the cell.
 * @param[out] cell   Pointer to a `Cell` structure to store the cell's vertices.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode GetCellVerticesFromGrid(Cmpnts ***coor, PetscInt idx, PetscInt idy, PetscInt idz, Cell *cell);

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
 * @brief Determines the spatial relationship of a particle relative to a cubic cell.
 *
 * This function evaluates whether a particle located at a specific point `p` in 3D space
 * is positioned inside the cell, on the boundary of the cell, or outside the cell. The
 * determination is based on the signed distances from the particle to each of the six
 * faces of the cell. The function utilizes the computed distances to ascertain the particle's
 * position with respect to the cell boundaries, considering a threshold to account for
 * floating-point precision.
 *
 * @param[in]  cell      A pointer to a `Cell` structure that defines the cubic cell via its
 *                       vertices. The cell's geometry is essential for accurately computing
 *                       the distances to each face.
 * @param[in]  d         A pointer to an array of six `PetscReal` values that store the
 *                       signed distances from the particle to each face of the cell. These
 *                       distances are typically computed using the `CalculateDistancesToCellFaces`
 *                       function.
 * @param[in]  p         The location of the particle in 3D space, represented by the `Cmpnts`
 *                       structure. This point is the reference for distance calculations to the
 *                       cell's faces.
 * @param[out] position  A pointer to an integer that will be set based on the particle's position
 *                       relative to the cell:
 *                       - `0`: The particle is inside the cell.
 *                       - `1`: The particle is on the boundary of the cell.
 *                       - `-1`: The particle is outside the cell.
 * @param[in]  threshold A `PetscReal` value that defines the minimum distance below which a
 *                       computed distance is considered to be zero. This threshold helps in
 *                       mitigating inaccuracies due to floating-point arithmetic, especially
 *                       when determining if the particle lies exactly on the boundary.
 *
 * @return PetscErrorCode Returns `0` if the function executes successfully. If an error occurs,
 *                        a non-zero error code is returned, indicating the type of failure.
 *
 * @note
 * - It is assumed that the `d` array has been properly allocated and contains valid distance
 *   measurements before calling this function.
 * - The function relies on `CalculateDistancesToCellFaces` to accurately compute the signed
 *   distances to each face. Any inaccuracies in distance calculations can affect the
 *   determination of the particle's position.
 * - The `threshold` parameter should be chosen based on the specific precision requirements
 *   of the application to balance between sensitivity and robustness against floating-point
 *   errors.
 * - The function includes a debug statement that prints the face distances, which can be useful
 *   for verifying the correctness of distance computations during development or troubleshooting.
 */
 PetscErrorCode EvaluateParticlePosition(const Cell *cell, PetscReal *d, const Cmpnts p, PetscInt *position, const PetscReal threshold);

/**
 * @brief Locates the grid cell containing a particle using a walking search.
 * @ingroup ParticleLocation
 *
 * This function implements a walking search algorithm to find the specific cell
 * (identified by global indices i, j, k) that encloses the particle's physical
 * location (`particle->loc`).
 *
 * The search starts from an initial guess cell (either the particle's previously known
 * cell or the corner of the local process domain) and iteratively steps to adjacent
 * cells based on the signed distances from the particle to the faces of the current
 * cell.
 *
 * Upon successful location (`position == 0`), the particle's `cell` field is updated
 * with the found indices (i,j,k), and the corresponding interpolation `weights`
 * are calculated and stored using the distances (`d`) relative to the final cell.
 *
 * Handles particles exactly on cell boundaries by attempting a tie-breaker if the
 * search gets stuck oscillating between adjacent cells due to the boundary condition.
 *
 * @param[in]  user     Pointer to the UserCtx structure containing grid information (DMDA, coordinates)
 *                      and domain boundaries.
 * @param[in,out] particle Pointer to the Particle structure. Its `loc` field provides the
 *                      target position. On successful return, its `cell` and `weights`
 *                      fields are updated. On failure, `cell` is set to `{-1, -1, -1}`
 *                      and `weights` to `{0.0, 0.0, 0.0}`.
 *
 * @return PetscErrorCode 0 on success. Non-zero error codes may indicate issues during
 *                        coordinate access, distance calculation, or other internal errors.
 *                        A return code of 0 does not guarantee the particle was found;
 *                        check `particle->cell[0] >= 0` afterward.
 *
 * @note Relies on helper functions like `InitializeTraversalParameters`, `CheckCellWithinLocalGrid`,
 *       `RetrieveCurrentCell`, `EvaluateParticlePosition`, `UpdateCellIndicesBasedOnDistances`,
 *       `UpdateParticleWeights`, and `FinalizeTraversal`.
 * @warning Ensure `particle->loc` is set correctly before calling.
 * @warning The function may fail to find the particle if it lies outside the domain accessible
 *          by the current process (including ghost cells) or if `MAX_TRAVERSAL` steps are exceeded.
 */
PetscErrorCode LocateParticleInGrid(UserCtx *user, Particle *particle);

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
 * @param[in]  info DMDALocalInfo structure that holds local & global domain bounds.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode UpdateCellIndicesBasedOnDistances( PetscReal d[NUM_FACES], PetscInt *idx, PetscInt *idy, PetscInt *idz, DMDALocalInfo *info);


/**
 * @brief Finalizes the traversal by reporting the results.
 *
 * This function prints the outcome of the traversal, indicating whether the particle
 * was found within a cell or not, and updates the particle's cell indices accordingly.
 *
 * @param[in]  user           Pointer to the user-defined context containing grid information.
 * @param[out] particle       Pointer to the Particle structure to update with cell indices.
 * @param[in]  traversal_steps The number of traversal steps taken.
 * @param[in]  cell_found      Flag indicating whether the particle was found within a cell.
 * @param[in]  idx             The i-index of the found cell.
 * @param[in]  idy             The j-index of the found cell.
 * @param[in]  idz             The k-index of the found cell.
 *
 * @return PetscErrorCode     Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeTraversal(UserCtx *user, Particle *particle, PetscInt traversal_steps, PetscBool cell_found, PetscInt idx, PetscInt idy, PetscInt idz);


/**
 * @brief Finds the MPI rank that owns a given global cell index.
 * @ingroup DomainInfo
 *
 * This function performs a linear search through the pre-computed decomposition map
 * (`user->RankCellInfoMap`) to determine which process is responsible for the cell
 * with global indices (i, j, k). It is the definitive method for resolving cell
 * ownership in the "Walk and Handoff" migration algorithm.
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
 * @brief Locates a particle's host cell or identifies its migration target using a robust walk search.
 * @ingroup ParticleLocation
 *
 * This is the core search engine. It starts from a guess cell and walks through
 * the grid. It returns a definitive, actionable status indicating the outcome:
 * - `ACTIVE_AND_LOCATED`: The particle was found in a cell on the current rank.
 * - `MIGRATING_OUT`: The particle was found to belong to another rank. `particle->destination_rank` is set.
 * - `LOST`: The search failed to find the particle within the global domain.
 *
 * This function is globally aware and can walk across MPI rank boundaries. It contains
 * robust checks for global domain boundaries and a tie-breaker for numerically "stuck"
 * particles on cell faces.
 *
 * @param[in]     user         Pointer to the UserCtx containing all grid and domain info.
 * @param[in,out] particle     Pointer to the Particle struct. Its fields are updated based
 *                             on the search outcome.
 * @param[out]    status_out   The final, actionable status of the particle after the search.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode LocateParticleOrFindMigrationTarget_TEST(UserCtx *user,
                                                        Particle *particle,
                                                        ParticleLocationStatus *status_out);

/**
 * @brief Logs the final outcome of the particle location search.
 */
PetscErrorCode ReportSearchOutcome(const Particle *particle,
                                          ParticleLocationStatus status,
				   PetscInt traversal_steps);

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
 * @param[in]  info DMDALocalInfo structure that holds local & global domain bounds.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode UpdateCellIndicesBasedOnDistancesTEST( PetscReal d[NUM_FACES], PetscInt *idx, PetscInt *idy, PetscInt *idz);

#endif // WALKINGSEARCH_H
