#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>
#include <petscerror.h>
#include <petscsys.h>

// Include additional headers
#include "variables.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "setup.h"          // Setup functions that are used across the codebase

// Macros and constants
#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

/**
 * @brief Generic macro to call the appropriate interpolation function based on the field type.
 *
 * This macro will select either the scalar or vector interpolation function based on
 * the type of the 'field' pointer. It also performs a compile-time check that 'centfield'
 * is of the same type as 'field'. If the types are not the same, a compile-time error is produced.
 *
 * Usage:
 *   InterpolateFieldFromCornerToCenter(field, centfield, user);
 */
#define InterpolateFieldFromCornerToCenter(field, centfield, user)                         \
  ( (void)sizeof(char[1 - 2*!!(!__builtin_types_compatible_p(typeof(field), typeof(centfield)))]),\
    _Generic((field),                                                                       \
      PetscReal ***: InterpolateFieldFromCornerToCenter_Scalar,                             \
      Cmpnts ***:    InterpolateFieldFromCornerToCenter_Vector                                \
    )(field, centfield, user) )

/*
#define InterpolateFieldFromCenterToCorner(centfield, field, info)         \
  _Generic((centfield),                                                    \
    PetscReal ***: InterpolateFieldFromCenterToCorner_Scalar,                \
    Cmpnts ***:    InterpolateFieldFromCenterToCorner_Vector                  \
  )(centfield, field, info)
*/
  
/**
 * @brief Macro to dispatch to the correct scalar or vector center-to-corner function
 *        based on a runtime block size variable.
 *
 * This macro uses a ternary operator to inspect the runtime value of 'blockSize' and
 * select the appropriate implementation. It expects the input and output pointers to be void*
 * and handles casting them to the correct, strongly-typed pointers.
 *
 * @param blockSize     The runtime block size (1 for scalar, 3 for vector).
 * @param centfield_ptr The input void* pointer to the 3D cell-centered data array.
 * @param corner_ptr    The output void* pointer to the 3D corner data array.
 * @param user_ctx      The UserCtx structure.
 */
#define InterpolateFieldFromCenterToCorner(blockSize, centfield_ptr, corner_ptr, user_ctx) \
    ( (blockSize) == 1 ? \
        InterpolateFieldFromCenterToCorner_Scalar_Petsc((PetscReal***)(centfield_ptr), (PetscReal***)(corner_ptr), (user_ctx)) : \
        InterpolateFieldFromCenterToCorner_Vector_Petsc((Cmpnts***)(centfield_ptr), (Cmpnts***)(corner_ptr), (user_ctx)) \
    )


/**
 * @brief A type-generic macro that interpolates a field from corner nodes to all face centers.
 *
 * This macro uses C11's _Generic feature to dispatch to the appropriate underlying
 * function based on the type of the input `corner_arr`.
 *
 *   - If `corner_arr` is of type `PetscReal***`, it calls `InterpolateCornerToFaceCenter_Scalar`.
 *   - If `corner_arr` is of type `Cmpnts***`, it calls `InterpolateCornerToFaceCenter_Vector`.
 *
 * @param corner_arr  Ghosted node-centered array (global indexing). The type of this
 *                    argument determines which function is called.
 * @param faceX_arr   Local array for X-faces.
 * @param faceY_arr   Local array for Y-faces.
 * @param faceZ_arr   Local array for Z-faces.
 * @param user_ctx    User context containing DMDA 'fda'.
 */
#define InterpolateCornerToFaceCenter(corner_arr, faceX_arr, faceY_arr, faceZ_arr, user_ctx) \
    _Generic((corner_arr),                                                                  \
        PetscReal***: InterpolateCornerToFaceCenter_Scalar,                                 \
        Cmpnts***:    InterpolateCornerToFaceCenter_Vector                                  \
    )(corner_arr, faceX_arr, faceY_arr, faceZ_arr, user_ctx)

/**
 * @brief Macro that calls either the scalar or vector piecewise interpolation function
 *        based on the type of the `fieldPtr` parameter (3D array).
 *
 * Usage example:
 *
 *   // For scalar:
 *   PetscReal ***fieldScal;
 *   PetscReal outVal;
 *   PieceWiseLinearInterpolation(fieldName, fieldScal, i, j, k, &outVal);
 *
 *   // For vector:
 *   Cmpnts ***fieldVec;
 *   Cmpnts vec;
 *   PieceWiseLinearInterpolation(fieldName, fieldVec, i, j, k, &vec);
 */
#define PieceWiseLinearInterpolation(fieldName, fieldPtr, i, j, k, outPtr) \
  _Generic((fieldPtr),                                                     \
    PetscReal ***: PieceWiseLinearInterpolation_Scalar,                    \
    Cmpnts      ***: PieceWiseLinearInterpolation_Vector                   \
  )(fieldName, fieldPtr, i, j, k, outPtr)

/**
 * @brief Macro that calls either the scalar or vector trilinear interpolation function
 *        based on the type of the `fieldPtr` parameter (3D array).
 *
 * Usage example:
 * 
 *   PetscReal result;
 *   Cmpnts vec;
 *   
 *   // For scalars:
 *   TrilinearInterpolation(fieldName, fieldScal, i, j, k, a1, a2, a3, &result);
 *
 *   // For vectors:
 *   TrilinearInterpolation(fieldName, fieldVec, i, j, k, a1, a2, a3, &vec);
 */
#define TrilinearInterpolation(fieldName, fieldPtr, i, j, k, a1, a2, a3, outPtr)           \
  _Generic((fieldPtr),                                                                     \
    PetscReal ***: TrilinearInterpolation_Scalar,                                          \
    Cmpnts      ***: TrilinearInterpolation_Vector                                         \
  )(fieldName, fieldPtr, i, j, k, a1, a2, a3, outPtr)

// Function declarations

/**
 * @brief Computes the trilinear interpolated scalar at a given point.
 *
 * @param[in]  fieldName A string representing the name of the scalar field (e.g., "temperature").
 * @param[in]  fieldScal 3D array of the field from a DMDA (indexed as [k][j][i]),
 *                      each cell a PetscReal.
 * @param[in]  i, j, k   Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1, a2, a3 Normalized coordinates within the cell ([0,1] range).
 * @param[out] val       Pointer to a PetscReal that will store the interpolated scalar.
 *
 * This function uses the standard 8-corner trilinear formula via `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same interface.
 */
PetscErrorCode TrilinearInterpolation_Scalar(
    const char   *fieldName,
    PetscReal  ***fieldScal,
    PetscInt      i,
    PetscInt      j,
    PetscInt      k,
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    PetscReal    *val);

/**
 * @brief Computes the trilinear interpolated vector (e.g., velocity) at a given point.
 *
 * @param[in]  fieldName  A string representing the name of the vector field (e.g., "velocity").
 * @param[in]  fieldVec   3D array of the field from a DMDA (indexed as [k][j][i]),
 *                        each cell of type Cmpnts.
 * @param[in]  i, j, k    Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1, a2, a3 Normalized coordinates within the cell ([0,1] range).
 * @param[out] vec        Pointer to a Cmpnts struct that will store the interpolated vector (x, y, z).
 *
 * This function uses the standard 8-corner trilinear formula via `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same interface.
 */
PetscErrorCode TrilinearInterpolation_Vector(
    const char   *fieldName,
    Cmpnts     ***fieldVec,
    PetscInt      i,
    PetscInt      j,
    PetscInt      k,
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    Cmpnts       *vec);

/**
 * @brief Interpolates a cell-centered field (scalar or vector) onto DMSwarm particles,
 *        using a robust, PETSc-idiomatic two-stage process.
 *
 * This function first converts the cell-centered input data to corner-node data,
 * storing this intermediate result in a PETSc Vec to correctly handle the communication
 * of ghost-point information across parallel ranks. It then performs a final trilinear
 * interpolation from the ghosted corner data to each particle's location.
 *
 * Workflow:
 *   1.  Create temporary PETSc Vecs (`cornerGlobal`, `cornerLocal`) to manage the
 *       intermediate corner-node data, using the existing nodal DMDA (`user->fda`).
 *   2.  Call a dispatch macro that uses the runtime block size (`bs`) to select the
 *       correct underlying center-to-corner function, writing results into `cornerGlobal`.
 *   3.  Perform a ghost-point exchange (`DMGlobalToLocal`) to transfer the boundary
 *       data from `cornerGlobal` into the ghost regions of `cornerLocal`.
 *   4.  Loop over all local particles. For each particle:
 *       a. Convert its global cell index to a local index relative to the ghosted array.
 *       b. Check if the particle's interpolation stencil is fully contained within
 *          the owned+ghost region. If not, log a warning and set the result to zero.
 *       c. Perform the final trilinear interpolation using the ghosted `cornerLocal` data.
 *   5.  Restore all PETSc objects to prevent memory leaks.
 *
 * @param[in]  user                     User context with DMDA, DMSwarm, etc.
 * @param[in]  fieldGlobal_cellCentered Vec (from fda) with cell-centered data (e.g., Ucat).
 * @param[in]  fieldName                Human-readable field name for logging (e.g., "Ucat").
 * @param[in]  swarmOutFieldName        Name of the DMSwarm field where interpolation results go.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateEulerFieldToSwarm(
    UserCtx    *user,
    Vec         fieldLocal_cellCentered,
    const char *fieldName,
    const char *swarmOutFieldName);


/**
 * @brief Interpolates all relevant fields from the DMDA to the DMSwarm.
 *
 * Currently, it interpolates:
 *   - user->Ucat (vector field) into the DMSwarm field "swarmVelocity".
 *
 * To add more fields, duplicate the call to InterpolateOneFieldOverSwarm and provide:
 *   - The global Vec for that field (e.g. user->Tcat for temperature),
 *   - A human-readable field name (for logging),
 *   - A DMSwarm output field name (e.g. "swarmTemperature").
 *
 * @param[in,out] user Pointer to a UserCtx containing:
 *                     - user->da (DM for the grid),
 *                     - user->swarm (DMSwarm for particles),
 *                     - user->Ucat (Vec for the vector field),
 *                     - possibly more fields like user->Tcat, user->Pcat, etc.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateAllFieldsToSwarm(UserCtx *user);

/**
 * @brief Interpolates particle velocities using trilinear interpolation.
 *
 * @param[in] user Pointer to the user-defined context containing grid and swarm information.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */

PetscErrorCode InterpolateParticleVelocities(UserCtx *user);

/**
 * @brief Safely interpolate a scalar field from corner nodes (from the coordinate DM)
 *        to cell centers (from the cell-centered DM) using the provided UserCtx.
 *
 * For each cell center in the physical region of the cell-centered DM (fda), this function
 * averages the up to 8 surrounding scalar values from the coordinate DM (da). On boundaries,
 * where fewer corners are available, a partial average is computed.
 *
 * The coordinate DM (da) is built on corners (IM+1 x JM+1 x KM+1) while the cell-centered DM (fda)
 * covers the physical cells (IM x JM x KM). Index offsets are adjusted via DMDAGetLocalInfo.
 *
 * @param[in]  field     3D array of corner-based scalar data (from user->da).
 * @param[out] centfield 3D array for the interpolated cell-center scalar data (for user->fda).
 * @param[in]  user      User context containing:
 *                       - da  : DM for the coordinate (corner) data.
 *                       - fda : DM for the cell-centered data.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Scalar(
    PetscReal ***field,
    PetscReal ***centfield,
    UserCtx *user);

/**
 * @brief Safely interpolate a vector field from corner nodes (from the coordinate DM)
 *        to cell centers (from the cell-centered DM) using the provided UserCtx.
 *
 * For each cell center in the physical region of the cell-centered DM (fda), this function
 * averages the 8 surrounding corner values from the coordinate DM (da). The coordinate DM
 * (da) is built on corners (IM+1 x JM+1 x KM+1) while the cell-centered DM (fda) covers
 * the physical cells (IM x JM x KM). Index offsets are adjusted using DMDAGetLocalInfo.
 *
 * @param[in]  field     3D array of corner-based vector data (from user->da).
 * @param[out] centfield 3D array for the interpolated cell-center vector data (for user->fda).
 * @param[in]  user      User context containing:
 *                       - da  : DM for the coordinate (corner) data.
 *                       - fda : DM for the cell-centered data.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Vector(
    Cmpnts ***field,
    Cmpnts ***centfield,
    UserCtx *user);

/**
 * @brief Interpolates a scalar field from cell centers to corner nodes.
 *
 * This function estimates the value of a scalar field at each grid node by averaging
 * the values from the cell centers of the cells surrounding that node (up to 8).
 * It handles physical boundaries by averaging only the available adjacent cells.
 *
 * Assumes input `centfield_arr` is from a ghosted local vector associated with `user->da` (DOF=1, s=2)
 * and output `field_arr` is from a ghosted local vector also associated with `user->da` (DOF=1, s=2).
 * Input array uses GLOBAL cell indices, output array uses GLOBAL node indices.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of scalar data at cell centers,
 *                            accessed via GLOBAL cell indices (k=0..KM-1, j=0..JM-1, i=0..IM-1).
 * @param[out] field_arr      Output: 3D array (ghosted) where interpolated node values are stored,
 *                            accessed via GLOBAL node indices (k=0..KM, j=0..JM, i=0..IM).
 * @param[in]  user           User context containing DMDA information (da).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar(PetscReal ***field_arr,
                                                  PetscReal ***centfield_arr,
                                                  UserCtx *user);

/**
 * @brief Interpolates a vector field from cell centers to corner nodes.
 *
 * This function estimates the value of a vector field at each grid node by averaging
 * the vector values from the cell centers of the cells surrounding that node (up to 8).
 * It handles physical boundaries by averaging only the available adjacent cells.
 *
 * Assumes input `centfield_arr` is from a ghosted local vector (e.g., representing ucat,
 * stored using node-indexing convention) and output `field_arr` is a ghosted local
 * vector associated with `user->fda` (DOF=3, s=2), accessed using global node indices.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of vector data conceptually at cell centers,
 *                            accessed via GLOBAL indices respecting the storage convention
 *                            (e.g., `ucat[k][j][i]` uses node index `i` but represents cell `C(i,j,k)` for interior).
 * @param[out] field_arr      Output: 3D array (ghosted) where interpolated node values are stored,
 *                            accessed via GLOBAL node indices (k=0..KM, j=0..JM, i=0..IM).
 * @param[in]  user           User context containing DMDA information (da and fda).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector(Cmpnts ***field_arr,
                                                  Cmpnts ***centfield_arr,
                                                  UserCtx *user);

/**
 * @brief Interpolates a vector field from cell centers to corner nodes.
 *
 * This version is adapted to write directly into a ghosted local array obtained from DMDAVecGetArray(),
 * which allows using GLOBAL indices for writing to the OWNED portion of the array.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of cell-centered data, accessed via GLOBAL indices.
 * @param[out] corner_arr     Output: 3D array (ghosted) where interpolated node values are stored,
 *                            also accessed via GLOBAL indices for the owned part.
 * @param[in]  user           User context containing DMDA information.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector_Petsc(
    Cmpnts ***centfield_arr, /* Input: Ghosted local array from Vec (read) */
    Cmpnts ***corner_arr,    /* Output: Ghosted local array from Vec (write) */
    UserCtx *user);

/**
 * @brief Interpolates a scalar field from cell centers to corner nodes.
 *
 * This version is adapted to write directly into a ghosted local array obtained from DMDAVecGetArray(),
 * which allows using GLOBAL indices for writing to the OWNED portion of the array.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of scalar data at cell centers, accessed via GLOBAL indices.
 * @param[out] corner_arr     Output: 3D array (ghosted) where interpolated node values are stored,
 *                            also accessed via GLOBAL indices for the owned part.
 * @param[in]  user           User context containing DMDA information.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar_Petsc(
    PetscReal ***centfield_arr, /* Input: Ghosted local array from Vec (read) */
    PetscReal ***corner_arr,    /* Output: Ghosted local array from Vec (write) */
    UserCtx *user);

/**
 * @brief Determines the target Eulerian DM and expected DOF for scattering a given particle field.
 * @ingroup scatter_module
 *
 * Based on hardcoded rules mapping particle field names to user context DMs (da/fda).
 * This function encapsulates the policy of where different fields should be scattered.
 *
 * @param[in] user             Pointer to the UserCtx containing da and fda.
 * @param[in] particleFieldName Name of the particle field (e.g., "P", "Ucat").
 * @param[out] targetDM        Pointer to store the determined target DM (da or fda).
 * @param[out] expected_dof    Pointer to store the expected DOF (1 or 3) for this field.
 *
 * @return PetscErrorCode Returns 0 on success, PETSC_ERR_ARG_UNKNOWN if field name is not recognized,
 *         or PETSC_ERR_ARG_NULL for NULL inputs.
 */
PetscErrorCode GetScatterTargetInfo(UserCtx *user, const char *particleFieldName,
                                    DM *targetDM, PetscInt *expected_dof);


/**
 * @brief Accumulates a particle field (scalar or vector) into a target grid sum vector.
 * @ingroup scatter_module_internal
 *
 * This function iterates through local particles, identifies their cell using the
 * "DMSwarm_CellID" field, and adds the particle's field value (`particleFieldName`)
 * to the corresponding cell location in the `gridSumVec`. It handles both scalar
 * (DOF=1) and vector (DOF=3) fields automatically based on the DOF of `gridSumDM`.
 *
 * IMPORTANT: The caller must ensure `gridSumVec` is zeroed before calling this
 * function if a fresh sum calculation is desired.
 *
 * @param[in] swarm           The DMSwarm containing particles.
 * @param[in] particleFieldName Name of the field on the particles (must match DOF).
 * @param[in] gridSumDM       The DMDA associated with `gridSumVec`. Its DOF determines
 *                            how many components are accumulated.
 * @param[in,out] gridSumVec  The Vec (associated with `gridSumDM`) to accumulate sums into.
 *
 * @return PetscErrorCode 0 on success. Errors if fields don't exist or DMs are incompatible.
 */
PetscErrorCode AccumulateParticleField(DM swarm, const char *particleFieldName,
                                       DM gridSumDM, Vec gridSumVec);

/**
 * @brief Normalizes a grid vector of sums by a grid vector of counts to produce an average.
 * @ingroup scatter_module_internal
 *
 * Calculates avgVec[i] = sumVec[i] / countVec[i] for each component of each
 * OWNED cell where countVec[i] > 0. Sets avgVec[i] = 0 otherwise.
 * Handles both scalar (DOF=1) and vector (DOF=3) data fields based on `dataDM`.
 * Uses basic `VecGetArray`/`VecGetArrayRead` and manual index calculation.
 *
 * @param[in] countDM    The DMDA associated with `countVec` (must have DOF=1).
 * @param[in] countVec   The Vec containing particle counts per cell (read-only).
 * @param[in] dataDM     The DMDA associated with `sumVec` and `avgVec` (must have DOF=1 or DOF=3).
 * @param[in] sumVec     The Vec containing the accumulated sums per cell (read-only).
 * @param[in,out] avgVec The Vec where the calculated averages will be stored (overwritten).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode NormalizeGridVectorByCount(DM countDM, Vec countVec,
                                          DM dataDM, Vec sumVec, Vec avgVec);

/**
 * @brief Scatters a particle field (scalar or vector) to the corresponding Eulerian field average.
 * @ingroup scatter_module
 *
 * This is the main user-facing function. It determines the target Eulerian DM
 * based on the `particleFieldName`, validates the provided `eulerFieldAverageVec`
 * against the target DM, and then orchestrates the scatter operation by calling
 * the internal helper function `ScatterParticleFieldToEulerField_Internal`.
 * The final averaged result is stored IN-PLACE in `eulerFieldAverageVec`.
 *
 * @param[in] user                 Pointer to UserCtx containing da, fda, swarm, ParticleCount.
 * @param[in] particleFieldName    Name of the field in the DMSwarm (e.g., "P", "Ucat").
 * @param[in,out] eulerFieldAverageVec Pre-created Vec associated with the correct target DM
 *                                 (implicitly da or fda). Result stored here.
 *
 * @return PetscErrorCode 0 on success. Errors on NULL input, unrecognized field name,
 *         or incompatible target vector.
 */
PetscErrorCode ScatterParticleFieldToEulerField(UserCtx *user,
                                                const char *particleFieldName,
                                                Vec eulerFieldAverageVec);

/**
 * @brief Scatters a predefined set of particle fields to their corresponding Eulerian fields.
 * @ingroup scatter_module
 *
 * This convenience function calls the unified `ScatterParticleFieldToEulerField`
 * for a standard set of fields ("P", potentially others). It assumes the target
 * Eulerian Vec objects (e.g., `user->P`, `user->Ucat`) exist in the UserCtx structure
 * and are correctly associated with their respective DMs (`user->da` or `user->fda`).
 * It zeros the target Vecs before scattering.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing all required DMs,
 *                     Vecs (`ParticleCount`, target Eulerian fields like `P`, `Ucat`), and `swarm`.
 *
 * @return PetscErrorCode 0 on success. Errors if prerequisites (like ParticleCount)
 *         are missing or if underlying scatter calls fail.
 */
PetscErrorCode ScatterAllParticleFieldsToEulerFields(UserCtx *user);

/**
 * @brief Interpolates a scalar field from corner nodes to all face centers.
 *
 * This routine computes the average of the four corner-node values
 * defining each face of a hexahedral cell:
 *   - X-faces (perpendicular to X): face between (i-1,i) in X-dir
 *   - Y-faces (perpendicular to Y): face between (j-1,j) in Y-dir
 *   - Z-faces (perpendicular to Z): face between (k-1,k) in Z-dir
 *
 * @param[in]  corner_arr  Ghosted node-centered array (global indexing) from user->fda.
 * @param[out] faceX_arr   Local array for X-faces sized [zm][ym][xm+1].
 * @param[out] faceY_arr   Local array for Y-faces sized [zm][ym+1][xm].
 * @param[out] faceZ_arr   Local array for Z-faces sized [zm+1][ym][xm].
 * @param[in]  user        User context containing DMDA 'fda' and GetOwnedCellRange.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateCornerToFaceCenter_Scalar(
    PetscReal ***corner_arr,
    PetscReal ***faceX_arr,
    PetscReal ***faceY_arr,
    PetscReal ***faceZ_arr,
    UserCtx *user);

/**
 * @brief Interpolates a vector field from corner nodes to all face centers.
 *
 * Identical to the scalar version, except it averages each component of the
 * Cmpnts struct at the four corner-nodes per face.
 *
 * @param[in]  corner_arr  Ghosted 3-component array (global node indices).
 * @param[out] faceX_arr   Local array of Cmpnts for X-faces sized [zm][ym][xm+1].
 * @param[out] faceY_arr   Local array of Cmpnts for Y-faces sized [zm][ym+1][xm].
 * @param[out] faceZ_arr   Local array of Cmpnts for Z-faces sized [zm+1][ym][xm].
 * @param[in]  user        User context containing DMDA 'fda'.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateCornerToFaceCenter_Vector(
    Cmpnts ***corner_arr,
    Cmpnts ***faceX_arr,
    Cmpnts ***faceY_arr,
    Cmpnts ***faceZ_arr,
    UserCtx *user);

#endif // INTERPOLATION_H

