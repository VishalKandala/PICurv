#ifndef BODYFORCES_H
#define BODYFORCES_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h"
#include "Metric.h"

/**
 * @file BodyForces.h
 * @brief Momentum source terms added to the contravariant RHS.
 *
 * @section bf_contract Contract for a body force
 *
 * Every body force here is invoked from `ComputeBodyForces()` (`src/rhs.c`),
 * which is the single extension point. A new force should:
 *
 *  1. take `(UserCtx *user, Vec Rct)`,
 *  2. detect for itself whether it is active and return early if not,
 *  3. **accumulate** into `Rct` with `+=`, never assign.
 *
 * @section bf_state RULE: DO NOT ADVANCE PER-TIMESTEP STATE HERE
 *
 * `ComputeBodyForces()` is called from `ComputeRHS()`, which runs **once per
 * residual evaluation** - that is, once per Jameson RK stage under the Picard
 * solver and once per Newton residual evaluation (including every finite
 * difference probe) under Newton-Krylov. It is emphatically NOT called once per
 * physical timestep.
 *
 * Any force that carries state across calls - a filter, a ramp, a moving
 * average, an integral controller term - must therefore gate its update on
 * `simCtx->step` and reuse the resolved value for the rest of that step:
 *
 * @code
 *     if (simCtx->myForceStep != simCtx->step) {
 *         ... advance the state ...
 *         simCtx->myForceStep = simCtx->step;
 *     }
 * @endcode
 *
 * Advancing it unconditionally makes the applied force depend on how many
 * residual evaluations preceded it. That is history dependence, and it breaks
 * two things at once: `MomentumNewtonKrylov_FormResidual()` requires `F(X)` to
 * be a deterministic function of the trial vector alone, and the Picard
 * shadow-Jacobian estimate assumes body forces are a constant forcing with zero
 * velocity Jacobian.
 *
 * This is not hypothetical. The driven-flow smoothing EMA in
 * `ComputeDrivenChannelFlowSource()` had exactly this defect: the applied force
 * walked 0.5, 0.75, 0.875 ... of the way toward its target across evaluations
 * within a single timestep. `tests/smoke/run_driven_periodic_regression.sh`
 * asserts the force is piecewise constant per step; extend it when adding a
 * stateful force.
 */

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