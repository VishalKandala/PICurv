#ifndef BODYFORCES_H
#define BODYFORCES_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h"
#include "Metric.h"

/**
 * @brief Applies a momentum source term to drive flow in a periodic channel or pipe.
 *
 * This function is the "engine" of the driven flow control system. It operates by:
 * 1.  Introspecting the boundary condition handlers to see if a `DRIVEN_` flow
 *     handler is active on any face. This determines if a driven flow is enabled
 *     and in which direction ('X', 'Y', or 'Z').
 * 2.  If a driven flow is active, it reads the `bulkVelocityCorrection` value that
 *     was computed by the handler's `PreStep` method and stored in the `SimCtx`.
 * 3.  It translates this velocity correction into a momentum source term.
 * 4.  It adds this source term to the appropriate component of the contravariant
 *     RHS vector (`Rct`) for all fluid cells in the domain.
 *
 * If no driven flow handler is found, this function does nothing.
 *
 * @param user The UserCtx containing the simulation state for a single block.
 * @param Rct  The PETSc Vec for the contravariant RHS, which will be modified in-place.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeDrivenChannelFlowSource(UserCtx *user, Vec Rct);

#endif // BODYFORCES_H