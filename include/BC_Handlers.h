#ifndef BC_HANDLERS_H
#define BC_HANDLERS_H

#include "variables.h"
#include "logging.h"


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
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 */
PetscErrorCode Create_WallNoSlip(BoundaryCondition *bc);

/**
 * @brief Configures a BoundaryCondition object to behave as a constant velocity inlet.
 */
PetscErrorCode Create_InletConstantVelocity(BoundaryCondition *bc);

PetscErrorCode Create_InletParabolicProfile(BoundaryCondition *bc);

PetscErrorCode Create_OutletConservation(BoundaryCondition *bc);



#endif // BC_HANDLERS_H
