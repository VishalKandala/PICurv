/**
 * @file field_catalog.c
 * @brief Authoritative persistent-Eulerian-field catalog and runtime view resolver.
 */

#include "field_catalog.h"

#define FIELD_NO_VEC_OFFSET ((size_t)-1)

#define FIELD_ENTRY(field_id_, name_, alias1_, alias2_, dof_, dm_, layout_, sync_, availability_, capabilities_, global_member_, local_member_) \
    [field_id_] = {                                                                                                                     \
        field_id_, name_, alias1_, alias2_, dof_, dm_, layout_, sync_, availability_, capabilities_,                                  \
        offsetof(UserCtx, global_member_), offsetof(UserCtx, local_member_)                                                             \
    }

#define FIELD_COORDINATE_ENTRY(field_id_, name_, dof_, layout_, capabilities_) \
    [field_id_] = {                                                               \
        field_id_, name_, NULL, NULL, dof_, FIELD_DM_COORDINATES, layout_,         \
        FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS, capabilities_,             \
        FIELD_NO_VEC_OFFSET, FIELD_NO_VEC_OFFSET                                   \
    }

static const FieldDescriptor gFieldCatalog[FIELD_ID_COUNT] = {
    FIELD_COORDINATE_ENTRY(FIELD_ID_COORDINATES, "Coordinates", 3, FIELD_LAYOUT_NODE_CENTERED,
                           FIELD_CAPABILITY_GHOST_UPDATE),
    FIELD_ENTRY(FIELD_ID_UCAT, "Ucat", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, Ucat, lUcat),
    FIELD_ENTRY(FIELD_ID_UCONT, "Ucont", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_COMPONENT_STAGGERED,
                FIELD_SYNC_COMPONENT_STAGGERED, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_STAGGERED_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, Ucont, lUcont),
    FIELD_ENTRY(FIELD_ID_UCONT_O, "Ucont_o", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_COMPONENT_STAGGERED,
                FIELD_SYNC_COMPONENT_STAGGERED, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE, Ucont_o, lUcont_o),
    FIELD_ENTRY(FIELD_ID_UCONT_RM1, "Ucont_rm1", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_COMPONENT_STAGGERED,
                FIELD_SYNC_COMPONENT_STAGGERED, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_CHECKPOINT, Ucont_rm1, lUcont_rm1),
    FIELD_ENTRY(FIELD_ID_P, "P", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, P, lP),
    FIELD_ENTRY(FIELD_ID_NU_T, "Nu_t", "Eddy Viscosity", NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_TURBULENCE,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, Nu_t, lNu_t),
    /* Holds the model coefficient C that multiplies Delta^2 |S|, which is Cs^2 in the
       classical notation and is signed once backscatter is admitted. Exists only for
       the dynamic model; the constant model prescribes its coefficient from config. */
    FIELD_ENTRY(FIELD_ID_CS, "CS", "Cs", NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_LES_DYNAMIC,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, CS, lCs),
    /* Wall-model friction velocity, nonzero only in the first interior cell of a WALL
       face. It is the quantity a wall model is scored against, so it is checkpointed and
       exposed to postprocessing rather than being recomputed from the corrected velocity,
       which no longer carries the law that produced it. */
    FIELD_ENTRY(FIELD_ID_U_TAU, "Utau", "Friction Velocity", NULL, 1, FIELD_DM_DA,
                FIELD_LAYOUT_CELL_CENTERED, FIELD_SYNC_STANDARD,
                FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_WALL_MODEL,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, Friction_Velocity, lFriction_Velocity),
    /* The wall model's effective eddy viscosity at its own wall face. A wall-resolved
       run has no such face and carries zero here; the field exists so the viscous flux
       can deliver the modelled stress without re-deriving the wall distance and
       tangential speed that produced it. Derived state, so not checkpointed. */
    FIELD_ENTRY(FIELD_ID_NU_WALL, "NuWall", "Wall Eddy Viscosity", NULL, 1, FIELD_DM_DA,
                FIELD_LAYOUT_CELL_CENTERED, FIELD_SYNC_STANDARD,
                FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_WALL_MODEL,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC,
                Nu_Wall, lNu_Wall),
    FIELD_ENTRY(FIELD_ID_DIFFUSIVITY, "Diffusivity", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC, Diffusivity, lDiffusivity),
    FIELD_ENTRY(FIELD_ID_DIFFUSIVITY_GRADIENT, "DiffusivityGradient", NULL, NULL, 3, FIELD_DM_FDA,
                FIELD_LAYOUT_CELL_CENTERED, FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE, DiffusivityGradient, lDiffusivityGradient),
    FIELD_ENTRY(FIELD_ID_CSI, "Csi", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_I_FACE,
                FIELD_SYNC_I_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, Csi, lCsi),
    FIELD_ENTRY(FIELD_ID_ETA, "Eta", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_J_FACE,
                FIELD_SYNC_J_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, Eta, lEta),
    FIELD_ENTRY(FIELD_ID_ZET, "Zet", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_K_FACE,
                FIELD_SYNC_K_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, Zet, lZet),
    FIELD_ENTRY(FIELD_ID_NVERT, "Nvert", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC |
                FIELD_CAPABILITY_CHECKPOINT, Nvert, lNvert),
    FIELD_ENTRY(FIELD_ID_AJ, "Aj", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC, Aj, lAj),
    FIELD_ENTRY(FIELD_ID_CENT, "Cent", "Center-Coordinates", NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS, FIELD_CAPABILITY_GHOST_UPDATE, Cent, lCent),
    FIELD_ENTRY(FIELD_ID_GRID_SPACE, "GridSpace", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS, FIELD_CAPABILITY_GHOST_UPDATE, GridSpace, lGridSpace),
    FIELD_ENTRY(FIELD_ID_CENTX, "Centx", "X-Face-Centers", NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_I_FACE,
                FIELD_SYNC_I_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC |
                FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT, Centx, lCentx),
    FIELD_ENTRY(FIELD_ID_CENTY, "Centy", "Y-Face-Centers", NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_J_FACE,
                FIELD_SYNC_J_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC |
                FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT, Centy, lCenty),
    FIELD_ENTRY(FIELD_ID_CENTZ, "Centz", "Z-Face-Centers", NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_K_FACE,
                FIELD_SYNC_K_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC |
                FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT, Centz, lCentz),
    FIELD_ENTRY(FIELD_ID_ICSI, "ICsi", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_I_FACE,
                FIELD_SYNC_I_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, ICsi, lICsi),
    FIELD_ENTRY(FIELD_ID_IETA, "IEta", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_I_FACE,
                FIELD_SYNC_I_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, IEta, lIEta),
    FIELD_ENTRY(FIELD_ID_IZET, "IZet", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_I_FACE,
                FIELD_SYNC_I_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, IZet, lIZet),
    FIELD_ENTRY(FIELD_ID_JCSI, "JCsi", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_J_FACE,
                FIELD_SYNC_J_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, JCsi, lJCsi),
    FIELD_ENTRY(FIELD_ID_JETA, "JEta", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_J_FACE,
                FIELD_SYNC_J_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, JEta, lJEta),
    FIELD_ENTRY(FIELD_ID_JZET, "JZet", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_J_FACE,
                FIELD_SYNC_J_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, JZet, lJZet),
    FIELD_ENTRY(FIELD_ID_KCSI, "KCsi", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_K_FACE,
                FIELD_SYNC_K_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, KCsi, lKCsi),
    FIELD_ENTRY(FIELD_ID_KETA, "KEta", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_K_FACE,
                FIELD_SYNC_K_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, KEta, lKEta),
    FIELD_ENTRY(FIELD_ID_KZET, "KZet", NULL, NULL, 3, FIELD_DM_FDA, FIELD_LAYOUT_K_FACE,
                FIELD_SYNC_K_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, KZet, lKZet),
    FIELD_ENTRY(FIELD_ID_IAJ, "IAj", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_I_FACE,
                FIELD_SYNC_I_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, IAj, lIAj),
    FIELD_ENTRY(FIELD_ID_JAJ, "JAj", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_J_FACE,
                FIELD_SYNC_J_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, JAj, lJAj),
    FIELD_ENTRY(FIELD_ID_KAJ, "KAj", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_K_FACE,
                FIELD_SYNC_K_FACE, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_FACE_SYNC, KAj, lKAj),
    FIELD_ENTRY(FIELD_ID_PHI, "Phi", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_ALWAYS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_PERIODIC_CELL_SYNC, Phi, lPhi),
    FIELD_ENTRY(FIELD_ID_PSI, "Psi", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_PARTICLES,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_CHECKPOINT, Psi, lPsi),
    FIELD_ENTRY(FIELD_ID_NVERT_O, "Nvert_o", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE, Nvert_o, lNvert_o),
    FIELD_ENTRY(FIELD_ID_PARTICLE_COUNT, "ParticleCount", NULL, NULL, 1, FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_PARTICLES,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_CHECKPOINT, ParticleCount, lParticleCount),
    FIELD_ENTRY(FIELD_ID_K_OMEGA, "K_Omega", NULL, NULL, 2, FIELD_DM_FDA2, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_RANS,
                FIELD_CAPABILITY_GHOST_UPDATE | FIELD_CAPABILITY_CHECKPOINT, K_Omega, lK_Omega),
    FIELD_ENTRY(FIELD_ID_K_OMEGA_O, "K_Omega_o", NULL, NULL, 2, FIELD_DM_FDA2, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL | FIELD_AVAILABILITY_RANS,
                FIELD_CAPABILITY_GHOST_UPDATE, K_Omega_o, lK_Omega_o),
    /* Corner-staging workspace. Node-centered by construction: the interpolation
     * writes cell-centered data onto grid corners. Not checkpointed, since it is
     * transient scratch rebuilt on every conversion. */
    FIELD_ENTRY(FIELD_ID_CELL_SCALAR_AT_CORNER, "CellScalarAtCorner", NULL, NULL, 1,
                FIELD_DM_DA, FIELD_LAYOUT_NODE_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE, CellScalarAtCorner, lCellScalarAtCorner),
    FIELD_ENTRY(FIELD_ID_CELL_VECTOR_AT_CORNER, "CellVectorAtCorner", NULL, NULL, 3,
                FIELD_DM_FDA, FIELD_LAYOUT_NODE_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE, CellVectorAtCorner, lCellVectorAtCorner),
    /* Post-processing staging fields. A derived statistic is config-counted and so
     * has no compile-time offset of its own; staging the result here lets the
     * existing ghost, nodal-average, and logging paths address it by name instead
     * of each growing a second, view-based entry point. */
    FIELD_ENTRY(FIELD_ID_POST_SCALAR, "PostScalar", NULL, NULL, 1,
                FIELD_DM_DA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE, PostScalar, lPostScalar),
    FIELD_ENTRY(FIELD_ID_POST_VECTOR, "PostVector", NULL, NULL, 3,
                FIELD_DM_FDA, FIELD_LAYOUT_CELL_CENTERED,
                FIELD_SYNC_STANDARD, FIELD_AVAILABILITY_FINEST_LEVEL,
                FIELD_CAPABILITY_GHOST_UPDATE, PostVector, lPostVector)
};

_Static_assert(sizeof(gFieldCatalog) / sizeof(gFieldCatalog[0]) == FIELD_ID_COUNT,
               "Field catalog must contain one entry per FieldId.");

/**
 * @brief Returns the immutable catalog descriptor for one typed field identity.
 * @details Validates the enum range and the table's ID/index invariant before
 *          exposing the catalog-owned descriptor.
 * @see FieldGetDescriptor()
 */
PetscErrorCode FieldGetDescriptor(FieldId field_id, const FieldDescriptor **descriptor)
{
    PetscFunctionBeginUser;
    PetscCheck(descriptor != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Field descriptor output cannot be NULL.");
    PetscCheck(field_id >= 0 && field_id < FIELD_ID_COUNT, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Invalid FieldId value %d.", (int)field_id);
    PetscCheck(gFieldCatalog[field_id].id == field_id, PETSC_COMM_SELF, PETSC_ERR_PLIB,
               "Field catalog entry %d is not initialized consistently.", (int)field_id);

    *descriptor = &gFieldCatalog[field_id];
    PetscFunctionReturn(0);
}

/**
 * @brief Resolves a canonical field name or registered alias to a typed identity.
 * @details Name comparison is intentionally confined to ingress; numerical
 *          consumers receive the resolved @ref FieldId.
 * @see FieldIdFromName()
 */
PetscErrorCode FieldIdFromName(const char *field_name, FieldId *field_id)
{
    PetscFunctionBeginUser;
    PetscCheck(field_name != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Field name cannot be NULL.");
    PetscCheck(field_id != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "FieldId output cannot be NULL.");

    for (PetscInt index = 0; index < FIELD_ID_COUNT; ++index) {
        const FieldDescriptor *descriptor = &gFieldCatalog[index];
        PetscBool match = PETSC_FALSE;

        PetscCall(PetscStrcasecmp(field_name, descriptor->canonical_name, &match));
        if (!match && descriptor->alias_1) PetscCall(PetscStrcasecmp(field_name, descriptor->alias_1, &match));
        if (!match && descriptor->alias_2) PetscCall(PetscStrcasecmp(field_name, descriptor->alias_2, &match));
        if (match) {
            *field_id = descriptor->id;
            PetscFunctionReturn(0);
        }
    }

    *field_id = FIELD_ID_INVALID;
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Field name '%s' is not registered in the Eulerian field catalog.", field_name);
}

/**
 * @brief Returns catalog-owned printable text for a field identity.
 * @details Invalid enum values return a non-null diagnostic sentinel string.
 * @see FieldCanonicalName()
 */
const char *FieldCanonicalName(FieldId field_id)
{
    if (field_id < 0 || field_id >= FIELD_ID_COUNT) return "InvalidField";
    if (gFieldCatalog[field_id].id != field_id) return "InvalidField";
    return gFieldCatalog[field_id].canonical_name;
}

/**
 * @brief Converts layout metadata to the stable diagnostic label used by field loggers.
 * @see FieldLayoutName()
 */
const char *FieldLayoutName(FieldLayout layout)
{
    switch (layout) {
        case FIELD_LAYOUT_NODE_CENTERED:       return "Node-Centered";
        case FIELD_LAYOUT_CELL_CENTERED:       return "Cell-Centered";
        case FIELD_LAYOUT_I_FACE:              return "I-Face";
        case FIELD_LAYOUT_J_FACE:              return "J-Face";
        case FIELD_LAYOUT_K_FACE:              return "K-Face";
        case FIELD_LAYOUT_COMPONENT_STAGGERED: return "Component-Staggered";
        default:                               return "Invalid-Layout";
    }
}

/**
 * @brief Binds one descriptor to the PETSc objects already owned by a UserCtx.
 * @details Resolves DM and Vec handles without allocating, referencing, or
 *          destroying PETSc storage.
 * @see FieldGetView()
 */
PetscErrorCode FieldGetView(UserCtx *user, FieldId field_id, FieldView *view)
{
    const FieldDescriptor *descriptor = NULL;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "UserCtx cannot be NULL when resolving a field view.");
    PetscCheck(view != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Field view output cannot be NULL.");
    PetscCall(FieldGetDescriptor(field_id, &descriptor));

    view->descriptor = descriptor;
    view->dm = NULL;
    view->global_vec = NULL;
    view->local_vec = NULL;

    switch (descriptor->dm_kind) {
        case FIELD_DM_DA:
            view->dm = user->da;
            break;
        case FIELD_DM_FDA:
            view->dm = user->fda;
            break;
        case FIELD_DM_FDA2:
            view->dm = user->fda2;
            break;
        case FIELD_DM_COORDINATES:
            view->dm = user->fda;
            PetscCheck(user->da != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                       "Coordinate DM is unavailable for field '%s'.", descriptor->canonical_name);
            PetscCall(DMGetCoordinates(user->da, &view->global_vec));
            PetscCall(DMGetCoordinatesLocal(user->da, &view->local_vec));
            break;
        default:
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
                    "Field '%s' has an invalid DM selector.", descriptor->canonical_name);
    }

    if (descriptor->dm_kind != FIELD_DM_COORDINATES) {
        view->global_vec = *(Vec *)((char *)user + descriptor->global_vec_offset);
        view->local_vec = *(Vec *)((char *)user + descriptor->local_vec_offset);
    }

    PetscCheck(view->dm != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "DM for field '%s' is unavailable in this UserCtx.", descriptor->canonical_name);
    PetscCheck(view->global_vec != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Global vector for field '%s' is unavailable in this UserCtx.", descriptor->canonical_name);
    PetscCheck(view->local_vec != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Local vector for field '%s' is unavailable in this UserCtx.", descriptor->canonical_name);

    PetscFunctionReturn(0);
}

#undef FIELD_COORDINATE_ENTRY
#undef FIELD_ENTRY
#undef FIELD_NO_VEC_OFFSET
