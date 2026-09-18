#ifndef LES_H
#define LES_H

#include "variables.h"
#include "logging.h"
#include "Metric.h"
#include "setup.h"
#include "Filter.h"
#include "Boundaries.h"

/*================================================================================*
 *              SYMMETRIC TENSOR KERNELS (no LES-specific knowledge)              *
 *================================================================================*/

/**
 * @brief Forms a symmetric tensor from a vector's outer product with itself.
 *
 * Used to build the resolved velocity products the Leonard stress subtracts.
 *
 * @param[in] v Vector to square.
 * @return The tensor `v_i v_j`.
 */
SymTensor SymTensorSelfOuter(Cmpnts v);

/**
 * @brief Forms the linear combination `a*x + b*y`.
 *
 * One primitive covers addition, subtraction, and scaling, so tensor algebra at the
 * call sites reads as algebra instead of six repeated component expressions.
 *
 * @param[in] a Scalar multiplying @p x.
 * @param[in] x First tensor.
 * @param[in] b Scalar multiplying @p y.
 * @param[in] y Second tensor.
 * @return The combined tensor.
 */
SymTensor SymTensorCombine(PetscReal a, SymTensor x, PetscReal b, SymTensor y);

/**
 * @brief Returns the trace `t_kk`.
 * @param[in] t Tensor to trace.
 * @return Sum of the three diagonal components.
 */
PetscReal SymTensorTrace(SymTensor t);

/**
 * @brief Removes the isotropic part, returning `t_ij - (1/3) delta_ij t_kk`.
 *
 * The Germano model tensor must be deviatoric so that the trace of the Leonard
 * stress, which the incompressible pressure absorbs, cannot enter the contraction
 * through discrete divergence error.
 *
 * @param[in] t Tensor to project.
 * @return The trace-free part of @p t.
 */
SymTensor SymTensorDeviator(SymTensor t);

/**
 * @brief Contracts two symmetric tensors as `a_ij b_ij`.
 *
 * Only six components are stored, so the three off-diagonal products are counted
 * twice. Centralising that factor is the point: written out at each call site it is
 * a silent error whenever one term is missed.
 *
 * @param[in] a First tensor.
 * @param[in] b Second tensor.
 * @return The full nine-term contraction.
 */
PetscReal SymTensorContract(SymTensor a, SymTensor b);

/**
 * @brief Applies a symmetric tensor to a vector, returning `t_ij v_j`.
 *
 * The three stored off-diagonal components each stand for two entries of the full
 * tensor, so a caller writing this out by hand has six chances to reference the wrong
 * one. Centralising it is the same reason `SymTensorContract()` exists.
 *
 * @param[in] t Tensor to apply.
 * @param[in] v Vector to apply it to.
 * @return The resulting vector.
 */
Cmpnts SymTensorTimesVector(SymTensor t, Cmpnts v);

/**
 * @brief Returns the squared Frobenius norm `t_ij t_ij`.
 * @param[in] t Tensor to measure.
 * @return The self-contraction of @p t.
 */
PetscReal SymTensorNormSq(SymTensor t);

/*================================================================================*
 *                        STRAIN RATE AND FILTER WIDTH                            *
 *================================================================================*/

/**
 * @brief Builds the strain-rate tensor and its magnitude from a velocity gradient.
 *
 * `S_ij = 0.5 (du_i/dx_j + du_j/dx_i)` and `|S| = sqrt(2 S_ij S_ij)`. This is the
 * single definition of both quantities; the dynamic procedure and the eddy-viscosity
 * assembly consume it rather than each expanding the algebra themselves.
 *
 * @param[in]  dudx  Gradient of the x velocity component.
 * @param[in]  dvdx  Gradient of the y velocity component.
 * @param[in]  dwdx  Gradient of the z velocity component.
 * @param[out] strain Resolved strain-rate tensor; may be NULL if only the magnitude is wanted.
 * @param[out] magnitude Strain-rate magnitude `|S|`; may be NULL.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode StrainRateFromGradients(Cmpnts dudx, Cmpnts dvdx, Cmpnts dwdx,
                                       SymTensor *strain, PetscReal *magnitude);

/**
 * @brief Computes one cell's grid filter width under the selected width model.
 *
 * ::LES_FILTER_WIDTH_CUBE_ROOT_VOLUME uses the cell volume alone and so ignores
 * anisotropy; the other models take the cell's extents along its own grid directions
 * from @ref ComputeCellDirectionalExtents, so they follow the cell's shape and are
 * unchanged by its orientation.
 *
 * @param[in]  model Width model to apply.
 * @param[in]  aj    Inverse cell volume (the Jacobian determinant) at the cell.
 * @param[in]  csi   Xi-direction metric vector at the cell.
 * @param[in]  eta   Eta-direction metric vector at the cell.
 * @param[in]  zet   Zeta-direction metric vector at the cell.
 * @param[out] delta Resolved filter width.
 * @return PetscErrorCode 0 on success, or `PETSC_ERR_ARG_OUTOFRANGE` for a
 *         non-positive Jacobian.
 */
PetscErrorCode ComputeCellFilterWidth(LESFilterWidthModel model, PetscReal aj,
                                      Cmpnts csi, Cmpnts eta, Cmpnts zet, PetscReal *delta);

/*================================================================================*
 *                        GERMANO IDENTITY KERNELS                                *
 *================================================================================*/

/**
 * @brief Forms the Leonard stress `L_ij = (u_i u_j)^ - u^_i u^_j`.
 *
 * `L_ij` is the stress carried by the scales between the grid filter and the test
 * filter. It is computable entirely from the resolved field, which is what makes the
 * dynamic procedure possible.
 *
 * @param[in] velocity_filtered        Test-filtered resolved velocity `u^_i`.
 * @param[in] velocity_product_filtered Test filter of the products `(u_i u_j)^`.
 * @return The Leonard stress tensor.
 */
SymTensor LeonardStress(Cmpnts velocity_filtered, SymTensor velocity_product_filtered);

/**
 * @brief Forms the deviatoric Germano model tensor `M_ij`.
 *
 * `M_ij = -2 Delta^2 ( alpha |S^| S^_ij - (|S| S_ij)^ )`, made trace-free.
 *
 * The two terms are structurally different and must stay so. The first is the model
 * evaluated on the test-filtered strain, the product of separately filtered
 * quantities. The second is the grid-level model formed pointwise and then
 * test-filtered, that is, the filter of a product. Their difference is the entire
 * information content of the dynamic procedure: collapsing them into the same
 * expression yields a well-scaled, flow-responsive, and wrong coefficient that is
 * not the Germano-Lilly coefficient and cannot be cited as one.
 *
 * @param[in] delta                Grid filter width at the cell.
 * @param[in] alpha                Squared test-to-grid filter width ratio.
 * @param[in] strain_magnitude_filtered Test-filtered strain magnitude `|S|^`.
 * @param[in] strain_filtered      Test-filtered strain tensor `S^_ij`.
 * @param[in] strain_product_filtered Test filter of the product `(|S| S_ij)^`.
 * @return The deviatoric model tensor.
 */
SymTensor GermanoModelTensor(PetscReal delta, PetscReal alpha,
                             PetscReal strain_magnitude_filtered,
                             SymTensor strain_filtered,
                             SymTensor strain_product_filtered);

/*================================================================================*
 *                     COEFFICIENT AVERAGING AND LIMITING                         *
 *================================================================================*/

/**
 * @brief Resolves which logical directions the dynamic coefficient is averaged over.
 *
 * The resulting selector is passed to `PicurvSpatialRatioAverage()`, which performs
 * the reduction; this routine only decides which directions it spans.
 *
 * ::LES_AVERAGING_LOCAL selects none, ::LES_AVERAGING_GLOBAL selects all three, and
 * ::LES_AVERAGING_HOMOGENEOUS uses the directions the configuration names, defaulting
 * to the block's periodic axes when it names none. Periodicity is read from the
 * boundary face configuration, which is the same source the periodic field
 * synchronization uses.
 *
 * @param[in]  user      Block context supplying configuration and boundary faces.
 * @param[out] direction Per-direction averaging flags in `(xi, eta, zeta)` order.
 * @return PetscErrorCode 0 on success, or `PETSC_ERR_ARG_WRONGSTATE` when
 *         homogeneous averaging is requested but no direction can be resolved.
 */
PetscErrorCode ResolveLESAveragingDirections(UserCtx *user, PetscBool direction[3]);

/**
 * @brief Applies the configured admissible range to one model coefficient.
 *
 * The coefficient carries the sign of `<L:M>`; a negative value is backscatter,
 * energy returning from the unresolved scales, not an error. ::LES_CLIP_NONE keeps
 * it, leaving stability to the total-viscosity floor applied downstream.
 *
 * @param[in]  coefficient Raw coefficient from the Germano contraction.
 * @param[in]  config      LES configuration supplying the mode and ceiling.
 * @param[out] limited     Set to `PETSC_TRUE` when the value was modified; may be NULL.
 * @return The coefficient after limiting.
 */
PetscReal ClipModelCoefficient(PetscReal coefficient, const LESConfig *config, PetscBool *limited);

/**
 * @brief Builds the eddy viscosity from a model coefficient.
 *
 * `nu_t = C Delta^2 |S|`, floored so that `nu + nu_t >= min_viscosity_ratio * nu`.
 * The floor constrains the quantity that must stay positive for the momentum
 * operator to remain well posed, which is why backscatter is limited here rather
 * than by forcing the coefficient positive.
 *
 * @param[in] coefficient         Model coefficient, signed.
 * @param[in] delta               Grid filter width at the cell.
 * @param[in] strain_magnitude    Resolved strain magnitude `|S|`.
 * @param[in] molecular_viscosity Molecular kinematic viscosity.
 * @param[in] min_viscosity_ratio Lower bound on total viscosity as a multiple of molecular.
 * @return The admissible eddy viscosity.
 */
PetscReal EddyViscosityFromCoefficient(PetscReal coefficient, PetscReal delta,
                                       PetscReal strain_magnitude,
                                       PetscReal molecular_viscosity,
                                       PetscReal min_viscosity_ratio);

/**
 * @brief WALE eddy viscosity (Nicoud & Ducros 1999) from the resolved velocity gradient.
 *
 * `nu_t = (C_w Delta)^2 (S^d_ij S^d_ij)^(3/2) / [(S_ij S_ij)^(5/2) + (S^d_ij S^d_ij)^(5/4)]`,
 * where `S^d` is the traceless symmetric part of the squared gradient. It vanishes in pure
 * shear, so a laminar shear layer gets none, and it falls as the cube of wall distance
 * approaching a wall without a damping function. It needs no coefficient field, no test
 * filter, and no averaging, so it is valid on geometry with no homogeneous direction.
 *
 * @param[in] dudx Gradient of the x velocity component.
 * @param[in] dvdx Gradient of the y velocity component.
 * @param[in] dwdx Gradient of the z velocity component.
 * @param[in] delta Grid filter width at the cell.
 * @param[in] coefficient Model constant `C_w`.
 * @return The eddy viscosity, zero where the gradient carries no model invariant.
 */
PetscReal WALEEddyViscosity(Cmpnts dudx, Cmpnts dvdx, Cmpnts dwdx,
                            PetscReal delta, PetscReal coefficient);

/**
 * @brief Vreman eddy viscosity (Vreman 2004), weighting each grid direction by its own spacing.
 *
 * `nu_t = c sqrt(B_beta / (alpha_ij alpha_ij))`, with `alpha_ij = du_j/dx_i` and
 * `beta_ij = sum_m d_i^(m) d_j^(m)`, where `d^(m)` is the change in velocity across the
 * cell along grid direction `m` - `edges[m] . grad(u)`. On an aligned Cartesian grid that
 * is Vreman's `Delta_m^2 alpha_mi alpha_mj`; on a curvilinear grid it is the same quantity
 * measured along the cell's own edges, so the model is resolved per direction rather than
 * through one scalar width, and unchanged by the grid's orientation. It vanishes in pure
 * shear and needs no coefficient field, test filter, or averaging.
 *
 * @param[in] dudx Gradient of the x velocity component.
 * @param[in] dvdx Gradient of the y velocity component.
 * @param[in] dwdx Gradient of the z velocity component.
 * @param[in] edges Cell edge vectors along xi, eta, zeta, from @ref ComputeCellEdgeVectors.
 * @param[in] coefficient Model constant `c`.
 * @return The eddy viscosity, zero where `B_beta` vanishes.
 */
PetscReal VremanEddyViscosity(Cmpnts dudx, Cmpnts dvdx, Cmpnts dwdx,
                              const Cmpnts edges[3], PetscReal coefficient);

/**
 * @brief Returns the modelled subgrid kinetic energy at a cell.
 *
 * Uses the Yoshizawa relation `k_sgs = 2 C_I Delta^2 |S|^2`. Reported as a
 * diagnostic: compared against the resolved turbulent kinetic energy from the
 * field-statistics pipeline it measures how much of the energy the grid resolves.
 *
 * @param[in] yoshizawa_ci     Yoshizawa constant `C_I`.
 * @param[in] delta            Grid filter width at the cell.
 * @param[in] strain_magnitude Resolved strain magnitude `|S|`.
 * @return The modelled subgrid kinetic energy.
 */
PetscReal SubgridKineticEnergy(PetscReal yoshizawa_ci, PetscReal delta, PetscReal strain_magnitude);

/*================================================================================*
 *                              RUNTIME ENTRY POINTS                              *
 *================================================================================*/

/**
 * @brief Computes the dynamic Smagorinsky coefficient field for one block.
 *
 * Runs the Germano-Lilly procedure: test-filters the resolved field, forms the
 * Leonard stress and the model tensor, contracts them, averages the two
 * contractions over the configured directions, limits the quotient, and stores it in
 * `UserCtx::CS`.
 *
 * The stored quantity is the coefficient `C` that multiplies `Delta^2 |S|` in the
 * eddy viscosity, that is, `Cs^2` rather than `Cs`. It is signed under
 * ::LES_CLIP_NONE. Reported values are converted to `Cs` at the diagnostic boundary,
 * where the familiar target near 0.16 applies.
 *
 * Only ::DYNAMIC_SMAGORINSKY reaches this routine; the constant model prescribes its
 * coefficient from configuration and allocates no field. Callers honour
 * `LESConfig::dynamic_frequency`, and the caller must have refreshed `UserCtx::lUcat`
 * and its periodic images beforehand.
 *
 * @param[in,out] user The user context for a single computational block.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeSmagorinskyConstant(UserCtx *user);

/**
 * @brief Computes the turbulent eddy viscosity for one block.
 *
 * Evaluates `nu_t = C Delta^2 |S|` at every fluid cell and stores it in
 * `UserCtx::Nu_t`. The coefficient comes from `UserCtx::lCs` for the dynamic model
 * and directly from `LESConfig::constant_cs` for the constant model, which therefore
 * needs no coefficient field at all.
 *
 * @param[in,out] user The user context for a single computational block.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeEddyViscosityLES(UserCtx *user);

/**
 * @brief Appends one row of LES coefficient statistics to the run's log directory.
 *
 * Writes `les_coefficient.csv`: the effective coefficient reported as `Cs`, its
 * spatial spread, eddy-viscosity levels, the modelled subgrid energy, and the
 * fractions of the domain that were backscattering or limited before this step's
 * clipping. Those last two describe the pre-clipping state, which no stored field
 * preserves.
 *
 * Values are instantaneous volume statistics, not time averages; window averages of
 * the same model's fields come from the field-statistics pipeline instead.
 *
 * @param[in] user The user context for a single computational block.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode LogLESDiagnostics(UserCtx *user);

#endif // LES_H
