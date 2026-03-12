#ifndef BC_HANDLERS_H
#define BC_HANDLERS_H

#include "variables.h"
#include "logging.h"

//================================================================================
//              VALIDATORS
//================================================================================

/**
 * @brief (Private) Validates all consistency rules for a driven flow (channel/pipe) setup.
 *
 * This function enforces a strict set of rules to ensure a driven flow simulation is
 * configured correctly. It is called by the main `BoundarySystem_Validate` dispatcher.
 *
 * The validation rules are checked in a specific order:
 * 1.  Detect if any `DRIVEN_` handler is active. If not, the function returns immediately.
 * 2.  Ensure that no `INLET`, `OUTLET`, or `FARFIELD` boundary conditions exist anywhere in the
 *     domain, as they are physically incompatible with a pressure-driven flow model.
 * 3.  Verify that both faces in the driven direction are of `mathematical_type PERIODIC`.
 * 4.  Verify that both faces in the driven direction use the exact same `DRIVEN_` handler type.
 *
 * @param user The UserCtx for a single block.
 * @return PetscErrorCode 0 on success, non-zero PETSc error code on failure.
 */
PetscErrorCode Validate_DrivenFlowConfiguration(UserCtx *user);


//================================================================================
//
//              HANDLER "CONSTRUCTOR" FUNCTION DECLARATIONS
//
// Each function is responsible for populating a BoundaryCondition struct
// with the correct function pointers for its specific behavior. They are
// implemented in BC_Handlers.c and called by the factory in Boundaries.c.
//
//================================================================================

/**
 * @brief Configures a BoundaryCondition object to behave as a no-slip, stationary wall.
 *
 * @param bc A
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_WallNoSlip(BoundaryCondition *bc);

/**
 * @brief Configures a BoundaryCondition object to behave as a constant velocity inlet.
 *
 * @param bc Parameter `bc` passed to `Create_InletConstantVelocity()`.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_InletConstantVelocity(BoundaryCondition *bc);

/**
 * @brief Configures a BoundaryCondition object for a parabolic inlet profile.
 *
 * The constructed handler computes inlet velocity as a profile function of
 * transverse coordinates, typically used for laminar channel/pipe initialization.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_InletParabolicProfile(BoundaryCondition *bc);

/**
 * @brief Configures a BoundaryCondition object for conservative outlet treatment.
 *
 * The constructed handler applies outlet updates that preserve the solver's
 * global mass/flux consistency assumptions.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_OutletConservation(BoundaryCondition *bc);

/**
 * @brief Configures a BoundaryCondition object for geometric periodic coupling.
 *
 * This constructor wires periodic boundary callbacks that exchange values across
 * opposite faces according to the configured periodic directions.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_PeriodicGeometric(BoundaryCondition *bc);

/**
 * @brief Configures a BoundaryCondition object for periodic driven-flow forcing.
 *
 * This constructor wires the periodic callbacks that enforce a prescribed driving
 * strategy (for example constant target flux) on a periodic direction pair.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_PeriodicDrivenConstant(BoundaryCondition *bc);

#endif // BC_HANDLERS_H
