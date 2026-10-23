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
 * @deprecated Unreachable from any executable: its only caller is
 *             @ref ComputeCellCharacteristicLengthScale, which has none. Kept for now with
 *             its unit tests; do not add new callers - use @ref ComputeCellDirectionalExtents or
 *             @ref ComputeCellEdgeVectors, which do not depend on the cell's orientation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CalculateFaceNormalAndArea(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3], double *Ai, double *Aj, double *Ak);

/**
 * @brief Inverts the 3x3 matrix of face-area vectors to obtain the covariant directions.
 *
 * The caller fills the rows with `csi`, `eta` and `zet` - the face-area vectors, which
 * are the contravariant basis scaled by the cell Jacobian. The inverse's columns are the
 * covariant (tangent) directions `dx/dxi`, `dx/deta`, `dx/dzeta` divided by that
 * Jacobian, which is what `CalculateFaceNormalAndArea()` normalizes.
 *
 * Degeneracy is judged relative to the rows' own magnitude, because the determinant of
 * area vectors scales as the square of the cell volume: an absolute floor would reject a
 * well-formed cell purely for being small.
 *
 * @param covariantTensor   Input: rows are the three face-area vectors.
 * @param[out] contravariantTensor Inverse written in place.
 * @deprecated Unreachable from any executable: its only caller is
 *             @ref CalculateFaceNormalAndArea, which is reached only from the deprecated
 *             @ref ComputeCellCharacteristicLengthScale. Kept for now with its unit tests;
 *             do not add new callers - use @ref ComputeCellDirectionalExtents or
 *             @ref ComputeCellEdgeVectors, which do not depend on the cell's orientation.
 * @return PetscErrorCode 0 on success; `PETSC_ERR_MAT_LU_ZRPVT` when the rows are
 *         coplanar or one of them vanishes.
 */
PetscErrorCode InvertCovariantMetricTensor(double covariantTensor[3][3], double contravariantTensor[3][3]);

/**
 * @brief Computes characteristic length scales (dx, dy, dz) for a curvilinear cell.
 *
 * For a non-uniform, non-orthogonal cell, there is no single "dx". This function
 * computes an effective length scale in each Cartesian direction based on the cell
 * volume and the areas of its faces: the three outputs are the Cartesian components of
 * the cell diagonal, so they change when an otherwise identical cell is rotated. Use
 * @ref ComputeCellDirectionalExtents for a measure of the cell itself.
 *
 * @param ajc Cell Jacobian/volume metric at the cell center.
 * @param csi Covariant xi metric vector.
 * @param eta Covariant eta metric vector.
 * @param zet Covariant zeta metric vector.
 * @param[out] dx Effective physical length associated with xi variation.
 * @param[out] dy Effective physical length associated with eta variation.
 * @param[out] dz Effective physical length associated with zeta variation.
 * @deprecated No production caller. Its outputs are Cartesian projections of the cell
 *             diagonal, so they change when an identical cell is rotated; the LES filter
 *             widths and the Clark term stopped using it for that reason. Kept for now
 *             with its unit tests; do not add new callers - use @ref ComputeCellDirectionalExtents or
 *             @ref ComputeCellEdgeVectors, which do not depend on the cell's orientation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeCellCharacteristicLengthScale(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz);

/**
 * @brief Computes a cell's extent along each of its own grid directions.
 *
 * Each extent is the cell volume divided by the area of the face pair it crosses, which
 * is the distance between those two faces. Unlike
 * @ref ComputeCellCharacteristicLengthScale, whose outputs are the Cartesian components
 * of the cell diagonal and so change when the same cell is rotated, these depend on the
 * cell alone: a filter width built from them is the same for a cell in a straight
 * section and for an identical cell turned by a bend.
 *
 * @param ajc Inverse cell volume at the cell center.
 * @param csi Face-area vector of the xi faces.
 * @param eta Face-area vector of the eta faces.
 * @param zet Face-area vector of the zeta faces.
 * @param[out] l_xi   Extent across the xi faces.
 * @param[out] l_eta  Extent across the eta faces.
 * @param[out] l_zeta Extent across the zeta faces.
 * @return PetscErrorCode 0 on success; `PETSC_ERR_ARG_OUTOFRANGE` for a non-positive
 *         Jacobian or a face of zero area.
 */
PetscErrorCode ComputeCellDirectionalExtents(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet,
                                             double *l_xi, double *l_eta, double *l_zeta);

/**
 * @brief Computes a cell's edge vectors along its three grid directions.
 *
 * Each is the covariant basis vector `dx/dxi_k` over one index step: the physical
 * displacement from one cell to the next along that grid line, with both its direction
 * and its length. Projecting a velocity gradient onto it gives the change in velocity
 * across the cell in that direction, which is what a direction-resolved subgrid model
 * needs and what does not depend on how the grid is oriented.
 *
 * @param ajc Inverse cell volume at the cell center.
 * @param csi Face-area vector of the xi faces.
 * @param eta Face-area vector of the eta faces.
 * @param zet Face-area vector of the zeta faces.
 * @param[out] edges Edge vectors along xi, eta and zeta, in that order.
 * @return PetscErrorCode 0 on success; `PETSC_ERR_ARG_OUTOFRANGE` for a non-positive Jacobian.
 */
PetscErrorCode ComputeCellEdgeVectors(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet,
                                      Cmpnts edges[3]);

/**
 * @brief Builds translated periodic images for the dummy cell centers of periodic axes.
 *
 * PETSc wraps field indices but does not translate coordinates. This routine sets the dummy
 * slot before the first cell of each periodic axis to the last cell's center minus the
 * validated seam translation, and the slot after the last cell to the first cell's center
 * plus it, axis by axis with a ghost refresh in between, so edge and corner dummies carry the
 * translation of every periodic axis they cross. Grid spacing is not modified.
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
 * | 3 | **All negative** (`Aj_max < 0`) &rarr; abort: the grid is left-handed. Negating the metric vectors is not a valid repair (a mirrored duct solved that way was wrong by a factor of forty), so the grid must be renumbered along one logical axis instead. |
 * | 4 | Store `user->GridOrientation = +1`, the only orientation the solver accepts. |
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

/**
 * @brief Builds one face plane: the unit outward normal of (v2 - v1) x (v4 - v1) and the vertex mean.
 *
 * The single definition of a cell face plane: ComputeSignedDistanceToPlane() and the
 * precomputed planes both use it, so their distances are bit-identical. The normal points
 * away from `cell_centroid`. A degenerate face (normal magnitude below 1e-12) is reported
 * through `degenerate` and left unnormalised.
 *
 * @param[in]  v1            First face vertex; the normal's origin.
 * @param[in]  v2            Second face vertex; the normal's first edge is v2 - v1.
 * @param[in]  v3            Third face vertex; it enters only the centroid.
 * @param[in]  v4            Fourth face vertex; the normal's second edge is v4 - v1.
 * @param[in]  cell_centroid Centroid of the cell the face belongs to, which fixes the outward side.
 * @param[out] plane         The face plane.
 * @param[out] degenerate    Whether the normal magnitude is below 1e-12.
 * @return PETSC_TRUE when the orientation is ambiguous (the face plane contains the centroid offset).
 */
PetscBool MetricFacePlane(const Cmpnts v1, const Cmpnts v2, const Cmpnts v3, const Cmpnts v4,
                          const Cmpnts cell_centroid, FacePlane *plane, PetscBool *degenerate);

/**
 * @brief Builds the six face planes of a hexahedron, as the particle search tests them.
 *
 * Face vertex orders and orientation are those CalculateDistancesToCellFaces() has always
 * used; degenerate faces set their bit in `degenerate`. `size` is set to 0. This per-cell form
 * serves callers holding only a Cell; MetricComputeFacePlanes() precomputes the shared per-face form.
 *
 * @param[in]  cell   The cell, vertices in GetCellVerticesFromGrid() order.
 * @param[out] planes The cell's face planes.
 */
void MetricBuildCellFacePlanes(const Cell *cell, CellFacePlanes *planes);

/**
 * @brief Builds the block's particle face planes once, as catalogued per-face fields.
 *
 * Each face carries one plane, shared by the two cells on either side: a unit normal pointing
 * toward +i, +j or +k and the mean of its four vertices, built by MetricFacePlane() in the vertex
 * order of the minus face (LEFT, BOTTOM, BACK) of the cell above it and oriented by that cell
 * (by the cell below where only it is local). Grid-cell faces come from the node coordinates and
 * are stored in the I/J/K-face layouts (`PlaneNormal*`, `PlaneCentroid*`: the face whose lower
 * corner is node (i,j,k) sits at [k+1][j+1][i] for I, [k+1][j][i+1] for J, [k][j+1][i+1] for K);
 * dual-cell faces come from lCent and sit at their lower-corner cell centre (`DualPlaneNormal*`,
 * `DualPlaneCentroid*`). `SearchCellSize` holds each cell's GetCellCharacteristicSize() at its
 * shifted cell-centred index. A face with no plane (degenerate, or beyond the ghosted region) has a
 * zero normal. Owned values are the global fields; the locals are their ordinary ghost update
 * (UpdateLocalGhosts(), the plain scatter: no periodic face repair). Collective. Built by
 * CalculateAllGridMetrics() on the finest level of a particle run, after the coordinates and lCent
 * are final; a no-op once built. Coordinates are fixed thereafter (there are no moving grids); a
 * grid that changes must MetricDestroyFacePlanes() and rebuild.
 *
 * @param[in,out] user Block whose planes are built; MetricDestroyFacePlanes() frees them.
 * @return PetscErrorCode 0 on success, or a PETSc error.
 */
PetscErrorCode MetricComputeFacePlanes(UserCtx *user);

/**
 * @brief Frees the face-plane fields MetricComputeFacePlanes() built, if any.
 * @param[in,out] user Block whose planes are freed; their Vecs are reset to NULL.
 * @return PetscErrorCode 0 on success, or a PETSc error.
 */
PetscErrorCode MetricDestroyFacePlanes(UserCtx *user);

/** @brief Read-only flat arrays of one family of face planes (grid cells or dual cells). */
typedef struct FacePlaneArrays {
    PetscBool        dual;        ///< Dual-cell planes (from lCent) rather than grid-cell planes.
    const Cmpnts    *normal[3];   ///< Unit normals toward +i, +j, +k; zero where there is no plane.
    const Cmpnts    *centroid[3]; ///< Face vertex means.
    const PetscReal *size;        ///< Grid-cell sizes; NULL for dual cells.
} FacePlaneArrays;

/**
 * @brief Opens a family's face-plane arrays for reading; MetricComputeFacePlanes() must have run.
 * @param[in]  user   Block whose planes are read.
 * @param[in]  dual   PETSC_TRUE for the dual-cell planes, PETSC_FALSE for the grid-cell planes.
 * @param[out] arrays The open arrays; close with MetricRestoreFacePlaneArrays().
 * @return PetscErrorCode 0 on success, or PETSC_ERR_ARG_WRONGSTATE if the planes are not built.
 */
PetscErrorCode MetricGetFacePlaneArrays(UserCtx *user, PetscBool dual, FacePlaneArrays *arrays);

/**
 * @brief Closes arrays opened by MetricGetFacePlaneArrays().
 * @param[in]     user   Block whose planes were read.
 * @param[in,out] arrays The arrays to close.
 * @return PetscErrorCode 0 on success, or a PETSc error.
 */
PetscErrorCode MetricRestoreFacePlaneArrays(UserCtx *user, FacePlaneArrays *arrays);

/**
 * @brief Assembles one cell's six outward face planes from the shared per-face planes.
 *
 * A cell's plus faces use the stored normals and its minus faces their exact negation; a face
 * whose stored normal is zero sets its `degenerate` bit. For a grid cell `(i, j, k)` is its origin
 * node; for a dual cell, the lower-corner cell centre.
 *
 * @param[in]  arrays Open plane arrays.
 * @param[in]  info   The block's DMDA local info, for its ghosted extents.
 * @param[in]  i      Cell origin i index.
 * @param[in]  j      Cell origin j index.
 * @param[in]  k      Cell origin k index.
 * @param[out] planes The cell's planes, with `size` from `SearchCellSize` (0 for dual cells).
 */
void MetricGatherCellFacePlanes(const FacePlaneArrays *arrays, const DMDALocalInfo *info,
                                PetscInt i, PetscInt j, PetscInt k, CellFacePlanes *planes);

/**
 * @brief Index of point (i, j, k) in a block's ghosted local arrays (coordinates, lCent, planes).
 * @param[in] info The block's DMDA local info, for its ghosted extents.
 * @param[in] i    Global i index, within the ghosted range.
 * @param[in] j    Global j index, within the ghosted range.
 * @param[in] k    Global k index, within the ghosted range.
 * @return The flat index, x fastest.
 */
static inline PetscInt MetricCellFacePlanesIndex(const DMDALocalInfo *info, PetscInt i, PetscInt j, PetscInt k)
{
    const PetscInt di = i - info->gxs, dj = j - info->gys, dk = k - info->gzs;
    return (dk * info->gym + dj) * info->gxm + di;
}
#endif /* METRIC_H */
