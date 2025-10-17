#ifndef SOLVERS_H
#define SOLVERS_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "implicitsolvers.h"
#include "rhs.h"
#include "logging.h"
#include "poisson.h"
#include "setup.h"
#include "les.h"
/*================================================================================*
 *                      HIGH-LEVEL SOLVER ORCHESTRATOR                            *
 *================================================================================*/

/**
 * @brief Orchestrates a single time step of the Eulerian fluid solver.
 *
 * This is the refactored, high-level entry point for advancing the fluid state
 * from time t_n to t_{n+1}. It takes the master SimCtx as its primary argument.
 *
 * @param simCtx The master simulation context, containing all solver settings,
 *               multigrid structures, and data vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode FlowSolver(SimCtx *simCtx);

#endif // SOLVERS_H
