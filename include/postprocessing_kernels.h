#ifndef POSTPROCESSING_KERNELS_H
#define POSTPROCESSING_KERNELS_H

#include "variables.h"
#include "logging.h"
#include "io.h" // For the UpdateLocalGhosts function prototype

// Function prototypes for post-processing kernels

/**
 * @brief Interpolates a cell-centered field to nodal locations using local stencil averaging.
 *
 * The kernel reads the input field by name, computes nodal values, and stores the
 * output in the named destination field. Both fields must already exist in the
 * current `UserCtx`.
 *
 * @param[in,out] user           Block-level context that owns the source and destination vectors.
 * @param[in]     in_field_name  Name of the input field to sample.
 * @param[in]     out_field_name Name of the output field to populate.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeNodalAverage(UserCtx* user, const char* in_field_name, const char* out_field_name);

/**
 * @brief Populates the layout boundary of a field that was written on the interior only.
 *
 * @details This is **not** a boundary condition. It enforces no physics and satisfies no
 *          equation. It replicates the layout convention that boundary-condition
 *          application already gives solver state fields, so that stencil kernels —
 *          `ComputeNodalAverage()` above all — read defined values instead of the
 *          structural zeros an interior-only producer leaves behind.
 *
 *          Fields written by `PicurvWindowDerive()` and `ComputeQCriterion()` cover only
 *          the physical interior, so their layout boundary holds zeros that mean "never
 *          written" rather than "zero". A node on the domain edge averages four such
 *          entries with four real cells and comes out at roughly half its true value.
 *          Solver state fields never need this: `Ucat`, `P`, `Nu_t`, and `CS` all carry
 *          boundary values written when their boundary conditions were applied.
 *
 *          **Only periodic directions are handled.** There the correct value exists and
 *          is exact — the low dummy plane repeats the last physical plane and the high
 *          dummy plane repeats the first. On a non-periodic face the correct value
 *          depends on both the quantity and the boundary type, and no single convention
 *          serves a staging buffer that carries stresses, pressure, and eddy viscosity in
 *          turn; nothing is written there rather than something invented. The design for
 *          the non-periodic case is recorded in @ref 60_Field_Statistics_Planned_Extensions.
 *
 *          Operates on the global vector, so it must run **before** the caller's
 *          `UpdateLocalGhosts()`, which then carries the written values into every local
 *          ghost region. Rank interfaces are that scatter's responsibility and are never
 *          touched here.
 *
 *          Multi-block interface boundaries are out of scope.
 *
 * @param[in]     user       Block context supplying the DMDA layout and periodicity.
 * @param[in,out] global     Global vector whose layout boundary is populated.
 * @param[in]     components Degrees of freedom carried: 1 or 3.
 * @return Zero on success, or `PETSC_ERR_ARG_NULL` for a null argument, or
 *         `PETSC_ERR_ARG_OUTOFRANGE` for an unsupported component count.
 */
PetscErrorCode ExtendToLayoutBoundary(UserCtx *user, Vec global, PetscInt components);

/**
 * @brief Computes the Q-criterion diagnostic from the local velocity-gradient tensor.
 *
 * This kernel evaluates rotational versus strain-rate dominance and writes the
 * result into the configured Q-criterion output vector for visualization and flow
 * feature identification.
 *
 * @param[in,out] user Block-level context containing velocity fields and target output storage.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeQCriterion(UserCtx* user);

/**
 * @brief Normalizes pressure using the value at the configured logical grid point.
 *
 * The owning rank reads the reference value, shares it collectively, and every
 * rank subtracts it in-place from its distributed portion of the field.
 *
 * @param[in,out] user                Block-level context containing pressure and reference configuration.
 * @param[in]     relative_field_name Name of the field to normalize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode NormalizeRelativeField(UserCtx* user, const char* relative_field_name);

// Add more post-processing kernel prototypes as needed
//  =========================================================================
//  Dimensionalization Kernels
//  =========================================================================
/**
 * @brief Scales a specified field from non-dimensional to dimensional units in-place.
 *
 * This function acts as a dispatcher. It takes the string name of a field,
 * identifies the corresponding PETSc Vec object and the correct physical
 * scaling factor (e.g., U_ref for velocity, P_ref for pressure), and then
 * performs an in-place VecScale operation. It correctly handles the different
 * physical dimensions of Cartesian velocity vs. contravariant volume flux.
 *
 * @param[in,out] user        The UserCtx containing the PETSc Vecs to be modified.
 * @param[in]     field_name  The case-insensitive string name of the field to dimensionalize
 *                            (e.g., "Ucat", "P", "Ucont", "Coordinates", "ParticlePosition", "ParticleVelocity").
 * @return PetscErrorCode
 */
PetscErrorCode DimensionalizeField(UserCtx *user, const char *field_name);

/**
 * @brief Orchestrates the dimensionalization of all relevant fields loaded from a file.
 *
 * This function is intended to be called in the post-processor immediately after
 * all solver output has been read into memory. It calls DimensionalizeField() for each of the core
 * physical quantities to convert the entire loaded state from non-dimensional to
 * dimensional units, preparing it for analysis and visualization.
 *
 * @param[in,out] user The UserCtx containing all the fields to be dimensionalized.
 * @return PetscErrorCode
 */
PetscErrorCode DimensionalizeAllLoadedFields(UserCtx *user);

// ===========================================================================
// Particle Post-Processing Kernels
// ===========================================================================


/**
 * @brief Computes the specific kinetic energy (KE per unit mass) for each particle.
 *
 * This kernel calculates SKE = 0.5 * |velocity|^2. It requires that the
 * velocity field exists and will populate the specific kinetic energy field.
 * The output field must be registered before this kernel is called.
 *
 * @param user           The UserCtx containing the DMSwarm.
 * @param velocity_field The name of the input vector field for particle velocity.
 * @param ske_field      The name of the output scalar field to store specific KE.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeSpecificKE(UserCtx* user, const char* velocity_field, const char* ske_field);

/**
 * @brief Computes the displacement magnitude |r_i - r_0| for each particle (per-particle VTK kernel).
 *
 * Reference point r_0 = (simCtx->psrc_x, psrc_y, psrc_z).  Writes the scalar displacement to
 * post_swarm[disp_field].  This is a visualisation kernel only — use ComputeParticleMSD from
 * particle_statistics.h for quantitative global statistics.
 *
 * @param user       The UserCtx containing the DMSwarms.
 * @param disp_field Name of the output scalar field in post_swarm.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeDisplacement(UserCtx *user, const char *disp_field);

/**
 * @brief Derives one accumulated statistic and converts it to nodal values.
 *
 * The counterpart of ComputeNodalAverage() for accumulated window state: it takes an
 * index into a window's requested output set, normalizes the centered state behind
 * it, and leaves the result in the shared post-processing staging field ready for
 * output. The staging vectors are reused between calls, so a caller must consume the
 * result before requesting the next one.
 *
 * With `global_operations.dimensionalize` set, the result carries physical units: the
 * source field's reference scale raised to the power the derived kind carries: a mean
 * and an RMS are linear in it, a Reynolds stress and a turbulent kinetic energy
 * quadratic. A co-moment flux relates two possibly different fields, so its factor is
 * the product of their scales rather than one squared. That is what the per-field
 * scaling table alone cannot express, and why one blanket velocity factor would have
 * been wrong for three of the five kinds.
 *
 * @param[in]  user         Block context holding the accumulators and staging fields.
 * @param[in]  window_index Window whose state is derived.
 * @param[in]  outputs      Comma-separated output kinds the recipe requested.
 * @param[in]  output_index Index into that output set.
 * @param[out] out_name     Name of the derived field, window qualified.
 * @param[in]  name_size    Capacity of @p out_name.
 * @param[out] out_vec      Nodal vector holding the result; borrowed, not owned.
 * @param[out] out_components Component count of the result.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode ComputeWindowStatisticNodal(UserCtx *user, PetscInt window_index,
                                           const char *outputs, PetscInt output_index,
                                           char *out_name, size_t name_size,
                                           Vec *out_vec, PetscInt *out_components);

/**
 * @brief Appends one convergence row for an accumulated window to its CSV history.
 *
 * The counterpart of ComputeParticleMSD() for Eulerian window state: it reduces the
 * window to the few numbers that answer whether it has run long enough — sample
 * count, total weight, represented time, the per-point valid-fraction range, and the
 * mean turbulent kinetic energy — and appends them as one row per processed step.
 * No single field snapshot can answer that question, which is why the history exists
 * beside the field output rather than instead of it.
 *
 * The mean is taken over the fluid cells the window actually sampled. Cells outside
 * the target domain, and cells a moving mask excluded, hold zeros that mean "never
 * measured"; averaging over them would scale the result down by the fraction of the
 * grid the window never covered.
 *
 * @param[in] user          Block context holding the accumulators and staging fields.
 * @param[in] window_index  Window to summarize.
 * @param[in] output_prefix Output path prefix; the window name and `.csv` are appended.
 * @param[in] ti            Step being processed, recorded as the row's key.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode ComputeWindowStatisticsSummary(UserCtx *user, PetscInt window_index,
                                              const char *output_prefix, PetscInt ti);

#endif // POSTPROCESSING_KERNELS_H
