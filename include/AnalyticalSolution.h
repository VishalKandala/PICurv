#ifndef ANALYTICALSOLUTION_H
#define ANALYTICALSOLUTION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include additional headers
#include "variables.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "setup.h"          // Setup Module required for array allocation and deallocation
#include "interpolation.h"  // All the different interpolation routines required 

/**
 * @brief Generic macro to select the appropriate local Cartesian field setter based on field type.
 *
 * This macro dispatches to:
 *   - SetLocalCartesianField_Vector if the field value pointer is of type Cmpnts*
 *   - SetLocalCartesianField_Scalar if the field value pointer is of type PetscReal*
 *
 * @param fieldName  Field name for logging.
 * @param fieldValue Pointer to the cell's field value (either Cmpnts* or PetscReal*).
 * @param coor       Pointer to the coordinate (Cmpnts*) to be used in the computation.
 */
#define SetLocalCartesianField(fieldName, fieldValue, coor, FieldInitialization, Constant) _Generic((fieldValue), \
    Cmpnts*:    SetLocalCartesianField_Vector, \
    PetscReal*: SetLocalCartesianField_Scalar \
    )(fieldName, fieldValue, coor,FieldInitialization,Constant)



// function declarations //


/**
 * @brief Sets the local Cartesian scalar field based on input coordinates.
 *
 * This function computes the scalar field value by combining the sine of the input
 * coordinate values. In this example, the scalar field is computed as the sum of the
 * sine functions of the x, y, and z coordinates.
 *
 * @param[in]     fieldName   Pointer to a string representing the field name (for logging purposes).
 * @param[in,out] scalarField Pointer to the PetscReal where the computed scalar field value will be stored.
 * @param[in]     coor        Pointer to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
PetscErrorCode SetLocalCartesianField_Scalar(const char *fieldName, PetscReal *scalarField, Cmpnts *coor,PetscInt FieldInitialization, PetscReal Constant);

/**
 * @brief Sets the local Cartesian vector field based on input coordinates.
 *
 * This function computes the vector components by applying the sine function to the
 * input coordinate values along the x, y, and z directions.
 *
 * @param[in]     fieldName Pointer to a string representing the field name (for logging purposes).
 * @param[in,out] vecField  Pointer to the `Cmpnts` structure where the computed vector field will be stored.
 * @param[in]     coor      Pointer to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
PetscErrorCode SetLocalCartesianField_Vector(const char *fieldName, Cmpnts *vecField, Cmpnts *coor, PetscInt FieldInitialization, PetscReal Constant);

/**
 * @brief Sets an analytical Cartesian field (scalar or vector) for cell centers based on a field name.
 *
 * This function looks up the field within the user context by comparing the provided field name
 * against the supported names:
 *   - Vector fields: "Ucat" and "Ucont"
 *   - Scalar fields: "P" and "nvert"
 *
 * If the field is found, the function verifies that the DM block size matches the expected value
 * (3 for vector fields and 1 for scalar fields), retrieves local DM information for both the coordinate DM (da)
 * and the cell-centered DM (fda), interpolates the corner-based coordinates (from da) to cell centers,
 * and then updates the field using the generic helper macro SetLocalCartesianField. Interior cells are
 * updated with the interpolated coordinates, and boundary cells are updated using the original coordinate data.
 *
 * If the field name is not found in the user context, the function throws an error.
 *
 * @param[in]  user      Pointer to the UserCtx structure containing:
 *                         - da: DM for coordinate (corner) data.
 *                         - fda: DM for cell-centered data.
 *                         - Ucat: vector field (Cmpnts ***)
 *                         - Ucont: vector field (Cmpnts ***)
 *                         - P: scalar field (PetscReal ***)
 *                         - nvert: scalar field (PetscReal ***)
 * @param[in]  fieldName Name of the field to update.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetAnalyticalCartesianField(UserCtx *user, const char *fieldName);

/**
 * @brief Sets an analytical contravariant field (Ucont) at nodes, where each component
 *        is evaluated using the physical coordinates of its respective face center.
 *
 * This function:
 * 1. Interpolates nodal physical coordinates to the centers of i-faces, j-faces, and k-faces
 *    for the cells relevant to the current rank, storing them in temporary local arrays.
 * 2. For each node (i,j,k) owned by the current rank (where Ucont components are stored):
 *    a. Retrieves the pre-calculated physical coordinate of the i-face center.
 *    b. Evaluates the analytical function for U^1 using this i-face center coordinate.
 *    c. Retrieves the pre-calculated physical coordinate of the j-face center.
 *    d. Evaluates the analytical function for U^2 using this j-face center coordinate.
 *    e. Retrieves the pre-calculated physical coordinate of the k-face center.
 *    f. Evaluates the analytical function for U^3 using this k-face center coordinate.
 *    g. Stores these (U^1, U^2, U^3) into user->Ucont[k][j][i].
 *
 * @param[in]  user      Pointer to the UserCtx structure.
 * @param[in]  fieldName Name of the field to update (currently only "Ucont").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetAnalyticalContravariantField(UserCtx *user, const char *fieldName);


/**
 * @brief Applies the analytical solution to the position vector.
 *
 * This function updates each entry in the provided PETSc vector by computing its sine,
 * thereby replacing each position with sin(position).
 *
 * @param tempVec The PETSc Vec containing particle positions which will be used to store the velocities.
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode SetAnalyticalSolution(Vec tempVec, PetscInt FieldInitialization);

/**
 * @brief Applies analytical boundary conditions to a specified global vector or scalar field.
 *
 * This function acts as a dispatcher based on the fieldName. It determines the
 * appropriate DM and Vec from the UserCtx and calls field-specific helper routines
 * (e.g., ApplyAnalyticalBC_Vector, ApplyAnalyticalBC_Scalar) to set boundary values
 * according to a predefined analytical function (e.g., sin(coord)).
 *
 * This should be called AFTER setting interior values and BEFORE updating local ghosts.
 *
 * @param[in] user      Pointer to the UserCtx containing DMs, coordinates, and field Vecs.
 * @param[in] fieldName Name of the field to apply BCs to ("Ucat", "Ucont", "P", "Nvert").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ApplyAnalyticalBC(UserCtx *user, const char *fieldName);

/**
 * @brief Sets the initial values for the INTERIOR of a specified Eulerian field.
 *
 * This function initializes the interior nodes of `Ucont` based on a profile selected
 * by `user->FieldInitialization`. It explicitly skips any node that lies on a global
 * boundary, as those values are set by the Boundary System's `Initialize` methods.
 *
 * The initialization is directional, aligned with the primary INLET face that was
 * identified by the parser. This ensures the initial flow is physically meaningful.
 *
 * Supported `user->FieldInitialization` profiles for "Ucont":
 *  - 0: Zero Velocity. All interior components of Ucont are set to 0.
 *  - 1: Constant Normal Velocity. The contravariant velocity component normal to the
 *       inlet direction is set such that the physical velocity normal to those grid
 *       planes is a constant `uin`. Other contravariant components are zero.
 *  - 2: Poiseuille Normal Velocity. The contravariant component normal to the
 *       inlet direction is set with a parabolic profile.
 *
 * @param user      The main UserCtx struct, containing all simulation data and configuration.
 * @param fieldName A string ("Ucont" or "P") identifying which field to initialize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetInitialInteriorField(UserCtx *user, const char *fieldName);

#endif // ANALYTICALSOLUTION_H
