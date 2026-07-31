#ifndef METRIC_H
#define METRIC_H

// Include necessary headers
#include <petsc.h>
#include "variables.h"   // Common type definitions
#include "logging.h"  // Logging macros and definitions
#include <stdlib.h>
#include "io.h"
#include "setup.h"    // For SetDMDAProcLayout

/**
 * @brief Maps a logical point inside one hexahedral cell to physical space.
 *
 * Uses trilinear interpolation of the eight cell vertices. Logical coordinates
 * are measured from the cell's lower logical corner and normally lie in [0, 1].
 *
 * @param user Grid context used for cell indexing.
 * @param X Ghosted node-coordinate array for the block.
 * @param i Cell index in the xi direction.
 * @param j Cell index in the eta direction.
 * @param k Cell index in the zeta direction.
 * @param xi Local xi coordinate within the cell.
 * @param eta Local eta coordinate within the cell.
 * @param zta Local zeta coordinate within the cell.
 * @param[out] Xp Physical Cartesian position evaluated at the logical point.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode MetricLogicalToPhysical(UserCtx *user, const Cmpnts ***X,
                                       PetscInt i,PetscInt j,PetscInt k,
                                       PetscReal xi,PetscReal eta,PetscReal zta,
                                       Cmpnts *Xp);

/**
 * @brief Collects the eight node coordinates of one logical hexahedral cell.
 *
 * Vertex ordering matches the trilinear metric kernels in this module.
 *
 * @param user Grid context used for cell indexing.
 * @param X Ghosted node-coordinate array for the block.
 * @param i Cell index in the xi direction.
 * @param j Cell index in the eta direction.
 * @param k Cell index in the zeta direction.
 * @param[out] V Eight physical vertex coordinates in metric-kernel order.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode MetricGetCellVertices(UserCtx *user,const Cmpnts ***X,
                                     PetscInt i,PetscInt j,PetscInt k,
                                     Cmpnts V[8]);

/**
 * @brief Evaluates the trilinear mapping Jacobian and its determinant in one cell.
 *
 * The returned matrix is \f$\partial(x,y,z)/\partial(\xi,\eta,\zeta)\f$;
 * callers use its determinant to reject degenerate or inverted cells.
 *
 * @param user Grid context used for cell indexing.
 * @param X Ghosted node-coordinate array for the block.
 * @param i Cell index in the xi direction.
 * @param j Cell index in the eta direction.
 * @param k Cell index in the zeta direction.
 * @param xi Local xi coordinate within the cell.
 * @param eta Local eta coordinate within the cell.
 * @param zta Local zeta coordinate within the cell.
 * @param[out] J Mapping Jacobian in row-major Cartesian/logical form.
 * @param[out] detJ Determinant of J.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode MetricJacobian(UserCtx *user,const Cmpnts ***X,
                              PetscInt i,PetscInt j,PetscInt k,
                              PetscReal xi,PetscReal eta,PetscReal zta,
                              PetscReal J[3][3],PetscReal *detJ);

/**
 * @brief Converts a Cartesian velocity vector to contravariant logical components.
 *
 * Uses the inverse of the physical-to-logical Jacobian; `detJ` must describe
 * the same mapping as `J` and must be nonzero.
 *
 * @param J Mapping Jacobian at the evaluation point.
 * @param detJ Determinant of J.
 * @param u Cartesian velocity components.
 * @param[out] uc Contravariant xi/eta/zeta velocity components.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode MetricVelocityContravariant(const PetscReal J[3][3],
                                           PetscReal detJ,
                                           const PetscReal u[3],PetscReal uc[3]);

/**
 * @brief Computes the unit normal vectors and areas of the three faces of a computational cell.
 *
 * Given the metric vectors (csi, eta, zet), this function calculates the geometric
 * properties of the cell faces aligned with the i, j, and k directions.
 *
 * @param csi Covariant xi metric vector.
 * @param eta Covariant eta metric vector.
 * @param zet Covariant zeta metric vector.
 * @param[out] ni Unit normal of the xi-normal face.
 * @param[out] nj Unit normal of the eta-normal face.
 * @param[out] nk Unit normal of the zeta-normal face.
 * @param[out] Ai Area of the xi-normal face.
 * @param[out] Aj Area of the eta-normal face.
 * @param[out] Ak Area of the zeta-normal face.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CalculateFaceNormalAndArea(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3], double *Ai, double *Aj, double *Ak);

/**
 * @brief Inverts the 3x3 covariant metric tensor to obtain the contravariant metric tensor.
 *
 * In curvilinear coordinates, the input matrix `g` contains the dot products of the
 * covariant basis vectors (e.g., g_ij = e_i . e_j). Its inverse, `G`, is the
 * contravariant metric tensor, which is essential for transforming vectors and tensors
 * between coordinate systems.
 *
 * @param covariantTensor   Input: A 3x3 matrix representing the covariant metric tensor.
 * @param[out] contravariantTensor Inverse metric tensor written in place.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InvertCovariantMetricTensor(double covariantTensor[3][3], double contravariantTensor[3][3]);

/**
 * @brief Computes characteristic length scales (dx, dy, dz) for a curvilinear cell.
 *
 * For a non-uniform, non-orthogonal cell, there is no single "dx". This function
 * computes an effective length scale in each Cartesian direction based on the cell
 * volume and the areas of its faces.
 *
 * @param ajc Cell Jacobian/volume metric at the cell center.
 * @param csi Covariant xi metric vector.
 * @param eta Covariant eta metric vector.
 * @param zet Covariant zeta metric vector.
 * @param[out] dx Effective physical length associated with xi variation.
 * @param[out] dy Effective physical length associated with eta variation.
 * @param[out] dz Effective physical length associated with zeta variation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeCellCharacteristicLengthScale(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz);

/**
 * @brief Builds translated periodic images for cell centers and grid spacing.
 *
 * PETSc wraps field indices but does not translate coordinates. This routine
 * applies the validated per-axis geometric translation to wrapped center
 * coordinates and refreshes the associated local grid spacing.
 *
 * @param user The UserCtx containing grid and field data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyPeriodicCorrectionsToCellCentersAndSpacing(UserCtx *user);

/**
 * @brief Builds translated periodic images for i-face centers (Centx).
 *
 * Every active periodic axis can affect the ghost images of Centx. This must
 * be called after Centx is computed and before it is used for metric calculations.
 *
 * @param user The UserCtx containing grid and field data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyPeriodicCorrectionsToIFaceCenter(UserCtx *user);

/**
 * @brief Builds translated periodic images for j-face centers (Centy).
 *
 * Every active periodic axis can affect the ghost images of Centy. This must
 * be called after Centy is computed and before it is used for metric calculations.
 *
 * @param user The UserCtx containing grid and field data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyPeriodicCorrectionsToJFaceCenter(UserCtx *user);

/**
 * @brief Builds translated periodic images for k-face centers (Centz).
 *
 * Every active periodic axis can affect the ghost images of Centz. This must
 * be called after Centz is computed and before it is used for metric calculations.
 *
 * @param user The UserCtx containing grid and field data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyPeriodicCorrectionsToKFaceCenter(UserCtx *user);

/**
 * @brief Computes the primary face metric components (Csi, Eta, Zet), including
 *        boundary extrapolation, and stores them in the corresponding global Vec
 *        members of the UserCtx structure (user->Csi, user->Eta, user->Zet).
 *
 *        This is a self-contained routine that performs the following steps:
 *        1. Obtains local ghosted nodal coordinates using DMGetCoordinatesLocal.
 *        2. Calculates metrics for INTERIOR faces where finite difference stencils are valid.
 *        3. EXTRAPOLATES metrics for faces on the physical domain boundaries by copying
 *           from the nearest computed interior face.
 *        4. Assembles the global `user->Csi`, `user->Eta`, `user->Zet` Vecs.
 *        5. Updates the local ghosted `user->lCsi`, `user->lEta`, `user->lZet` Vecs.
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 *
 * @note
 *  - This function is a complete "compute and make ready" unit for Csi, Eta, and Zet.
 *  - It's recommended to call `VecZeroEntries` on user->Csi, Eta, Zet before this
 *    if they might contain old data.
 */
PetscErrorCode ComputeFaceMetrics(UserCtx *user);

/**
 * @brief Calculates the cell-centered inverse Jacobian determinant (1/J) for INTERIOR cells
 *        and stores it in `user->Aj`. This version includes boundary extrapolation.
 *
 *        Nodal coordinates are obtained internally.
 *        Refer to previous Doxygen comments for details on physical locations and
 *        storage convention (`aj_arr[k_n][j_n][i_n]` for cell `C(i_n-1,j_n-1,k_n-1)`).
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeCellCenteredJacobianInverse(UserCtx *user);

/**
 * @brief Verify and, when consistently inverted, repair the **right-handed** metric basis (`Csi`, `Eta`, `Zet`) and a
 *        **positive Jacobian** (`Aj`) over the whole domain.
 *
 * The metric-generation kernels are completely algebraic, so they will happily
 * deliver a *left-handed* basis if the mesh file enumerates nodes in the
 * opposite ζ-direction.  
 * This routine makes the orientation explicit and—if needed—repairs it
 * **once per run**:
 *
 * | Step | Action |
 * |------|--------|
 * | 1 | Compute global `Aj_min`, `Aj_max`.                          |
 * | 2 | **Mixed signs** (`Aj_min < 0 && Aj_max > 0`) &rarr; abort: the mesh is topologically inconsistent. |
 * | 3 | **All negative** (`Aj_max < 0`) &rarr; flip <br>`Csi`, `Eta`, `Zet`, `Aj` & update local ghosts. |
 * | 4 | Store `user->orientation = ±1` so BC / IC routines can apply sign-aware logic if they care about inlet direction. |
 *
 * @param[in,out] user  Fully initialised #UserCtx that already contains  
 *                      `Csi`, `Eta`, `Zet`, `Aj`, their **local** ghosts, and
 *                      valid distributed DMs.
 *
 * @return `0` on success or a PETSc error code on failure.
 *
 * @note  Call **immediately after** `ComputeCellCenteredJacobianInverse()` and
 *        before any routine that differentiates or applies BCs.
 *
 * @note Author metadata intentionally omitted in API docs.
 */
PetscErrorCode CheckAndFixGridOrientation(UserCtx *user);

/**
 * @brief Computes the physical location of cell centers and the spacing between them.
 *
 * This function calculates two key geometric properties from the nodal coordinates:
 * 1.  `Cent`: A vector field storing the (x,y,z) coordinates of the center of each grid cell.
 * 2.  `GridSpace`: A vector field storing the physical distance between adjacent
 *     cell centers in the i, j, and k computational directions.
 *
 * It is a direct adaptation of the corresponding logic from the legacy `FormMetrics`.
 *
 * @param user The UserCtx for a specific grid level. The function populates `user->Cent` and `user->GridSpace`.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ComputeCellCentersAndSpacing(UserCtx *user);

/**
 * @brief Computes metrics centered on constant-i faces (i-faces).
 *
 * This function calculates the metric terms (ICsi, IEta, IZet) and the inverse
 * Jacobian (IAj) located at the center of each i-face. The stencils use
 * i-face-centered coordinates (`Centx`) which must be computed first.
 * The logic is a direct adaptation of the legacy FormMetrics function.
 *
 * @param user The UserCtx for a specific grid level. Populates user->ICsi, etc.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ComputeIFaceMetrics(UserCtx *user);

/**
 * @brief Computes metrics centered on constant-j faces (j-faces).
 *
 * This function calculates the metric terms (`JCsi`, `JEta`, `JZet`) and the
 * inverse Jacobian (`JAj`) located at the geometric center of each constant-j
 * face. This is a critical step for staggered-grid finite difference schemes.
 *
 * The process is a direct and faithful refactoring of the corresponding logic
 * from the legacy `FormMetrics` function:
 * 1.  It first calculates the physical (x,y,z) coordinates of the center of
 *     each i-face and stores them in the `user->Centy` vector.
 * 2.  It then uses a boundary-aware, second-order finite difference stencil on
 *     the `Centy` field to compute the derivatives (e.g., d(x)/d(csi)).
 *     - Central differences are used in the grid interior.
 *     - One-sided differences are used at the physical domain boundaries.
 * 3.  Finally, these derivatives are used to compute the final metric terms and
 *     the inverse Jacobian, which are stored in their respective `Vec` objects.
 *
 * @param user The UserCtx for a specific grid level. This function populates
 *             the `user->JCsi`, `user->JEta`, `user->JZet`, and `user->JAj` vectors.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ComputeJFaceMetrics(UserCtx *user);

/**
 * @brief Computes metrics centered on constant-k faces (k-faces).
 *
 * This function calculates the metric terms (`KCsi`, `KEta`, `KZet`) and the
 * inverse Jacobian (`KAj`) located at the geometric center of each constant-j
 * face. This is a critical step for staggered-grid finite difference schemes.
 *
 * The process is a direct and faithful refactoring of the corresponding logic
 * from the legacy `FormMetrics` function:
 * 1.  It first calculates the physical (x,y,z) coordinates of the center of
 *     each i-face and stores them in the `user->Centz` vector.
 * 2.  It then uses a boundary-aware, second-order finite difference stencil on
 *     the `Centz` field to compute the derivatives (e.g., d(x)/d(csi)).
 *     - Central differences are used in the grid interior.
 *     - One-sided differences are used at the physical domain boundaries.
 * 3.  Finally, these derivatives are used to compute the final metric terms and
 *     the inverse Jacobian, which are stored in their respective `Vec` objects.
 *
 * @param user The UserCtx for a specific grid level. This function populates
 *             the `user->KCsi`, `user->KEta`, `user->KZet`, and `user->KAj` vectors.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ComputeKFaceMetrics(UserCtx *user);

/**
 * @brief Performs a diagnostic check on the divergence of the face area metric vectors.
 *
 * For a closed cell, the sum of the face area vectors should be zero (Gauss's
 * divergence theorem). This function computes a measure of this divergence and
 * reports the maximum value over the domain. A small value indicates a
 * well-formed grid. This is a direct adaptation of the legacy function.
 *
 * @param user The UserCtx for a specific grid level (typically the finest).
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ComputeMetricsDivergence(UserCtx *user);

/**
 * @brief Computes the max-min values of the grid metrics.
 *
 * This function serves as a diagnostic tool to assess the quality of the grid
 * metrics. It calculates the bounds of the face metrics (Csi, Eta, Zet).
 *
 * @param user The UserCtx, containing all necessary grid data.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeMetricNorms(UserCtx *user);

/**
 * @brief Orchestrates the calculation of all grid metrics.
 *
 * This function iterates through every UserCtx in the multigrid and multi-block
 * hierarchy. For each context, it calls a series of modern, modular helper
 * functions to compute the face metrics (Csi, Eta, Zet), the cell-centered
 * inverse Jacobian (Aj), and to validate the grid's orientation.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode
 */
PetscErrorCode CalculateAllGridMetrics(SimCtx *simCtx);
#endif /* METRIC_H */
