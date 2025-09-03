/**
 * @file variables.h
 * @author Seokkoo Kang (original author) and contributors
 * @brief Main header file for a complex fluid dynamics solver.
 *
 * This file defines the core data structures, global variables, and function
 * prototypes for a CFD application. It combines features for a curvilinear,
 * multi-block, immersed boundary (IBM) solver with fluid-structure interaction (FSI),
 * a modern particle tracking system, a modular boundary condition framework,
 * and post-processing utilities. It heavily utilizes the PETSc library.
 */

#ifndef VARIABLES_H
#define VARIABLES_H

/*================================================================================*
 *                        PETSC & SYSTEM LIBRARIES                                *
 *================================================================================*/

// --- PETSc Includes ---
#include "petscvec.h"
#include "petscdmda.h"
#include "petscksp.h"
#include "petscsnes.h"

// --- Standard C/C++ Includes ---
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// --- C++ Specific Includes ---
#if defined(__cplusplus)
#include <vector>
#include <algorithm>
#include <cassert>
#endif

// Define a C-compatible interface for C++ compilers
#ifdef __cplusplus
extern "C" {
#endif

/*================================================================================*
 *                       MACROS & GLOBAL VARIABLES                                *
 *================================================================================*/

/// Coefficient controlling the temporal accuracy scheme (e.g., 1.5 for BDF2).
#define COEF_TIME_ACCURACY 1.5
/// Defines the maximum number of data fields for VTK point data.
#define MAX_POINT_DATA_FIELDS 10

// --- Global Configuration Flags (defined in a .c file) ---
extern int i_periodic, j_periodic, k_periodic;  ///< Periodicity flags for the main grid.
extern int ii_periodic, jj_periodic, kk_periodic; ///< Periodicity flags for a secondary grid/context.
extern int les;                                 ///< Flag to enable Large Eddy Simulation.
extern int wallfunction;                        ///< Flag to enable wall function models.
extern int rans;                                ///< Flag to enable Reynolds-Averaged Navier-Stokes models.
extern int i_homo_filter, j_homo_filter, k_homo_filter; ///< Flags for applying a homogeneous filter in each direction.
extern int testfilter_ik, testfilter_1d;        ///< Flags for specific test filters.
extern int clark;                               ///< Flag for Clark model formulation.


/*================================================================================*
 *                     CORE DATA STRUCTURES & TYPEDEFS                            *
 *================================================================================*/

//--------------------------------------------------------------------------------
//                 1. FORWARD DECLARATIONS & BASIC TYPES
//--------------------------------------------------------------------------------

// --- Forward Declarations ---
// These declarations inform the compiler that these struct types exist, allowing
// pointers to them to be used before their full content is defined. This is
// essential for breaking circular dependencies between complex structs.
typedef struct BC_Param_s BC_Param;
typedef struct BoundaryCondition_s BoundaryCondition;
typedef struct BoundaryFaceConfig_s BoundaryFaceConfig;
typedef struct UserCtx UserCtx;

// --- Foundational Geometric and Data Types ---

/** @brief A 3D point or vector with PetscScalar components. */
typedef struct {
    PetscScalar x, y, z;
} Cmpnts;

/** @brief A 2D point or vector with PetscScalar components. */
typedef struct {
    PetscScalar x, y;
} Cmpnts2;

/** @brief A 2D vector of PETSc real numbers (for geometry/coordinates). */
typedef struct {
  PetscReal	x, y;
} Cpt2D;

/** @brief Represents a single point in a time-varying flow waveform. */
typedef struct {
    PetscReal t, f;
} FlowWave;

/** @brief Legacy boundary condition data. Likely superseded by the newer BC system. */
typedef struct {
  Vec Ubcs;   ///< Boundary condition velocity values. (Comment: "An ugly hack, waste of memory")
  Vec Uch;    ///< Characteristic velocity for boundary conditions.
} BCS;


//--------------------------------------------------------------------------------
//                 2. PARTICLE LOCATION SYSTEM ENUMS
//--------------------------------------------------------------------------------

/**
 * @enum  ParticleLocationStatus
 * @brief Defines the state of a particle with respect to its location and migration
 *        status during the iterative SettleParticles() process.
 *
 * This enum is used to control the logic within the main particle relocation
 * orchestrator, ensuring that each particle is processed correctly and efficiently.
 * It is registered as a PETSC_INT field in the DMSwarm.
 */
typedef enum {
    /**
     * @brief The particle's host cell is unknown or invalid.
     * This is the initial state for newly arrived particles on a rank after migration,
     * or for particles whose location search failed in a previous step. Any particle
     * in this state is a candidate for processing in the current relocation pass.
     */
    NEEDS_LOCATION,

    /**
     * @brief The particle has been successfully located in a cell owned by the current rank.
     * This state indicates that the particle is settled for the current timestep and
     * does not need further location checks or migration. The orchestrator can skip
     * this particle in subsequent passes of the do-while loop.
     */
    ACTIVE_AND_LOCATED,

    /**
     * @brief The particle has been identified as needing to move to another rank.
     * This status is set by the sending rank during the identification phase. The particle's
     * `destination_rank` field will be valid. This state is temporary and primarily used
     * to signal that the particle should be included in the MPI migration step.
     */
    MIGRATING_OUT,

    /**
     * @brief The particle could not be located and is considered lost.
     * This status is set when the walk search exceeds its maximum traversal limit or
     * walks outside the global domain boundaries. Particles in this state are candidates
     * for removal by the `CheckAndRemoveOutOfBoundsParticles` function.
     */
    LOST,

    /**
     * @brief The default state for a newly created particle before any processing.
     * While `NEEDS_LOCATION` is functionally similar for the first pass, having a distinct
     * initial state can be useful for debugging and initialization logic.
     */
    UNINITIALIZED

} ParticleLocationStatus;

/** @brief Enumerates the six faces of a cubic cell for distance calculations. */
typedef enum {
    LEFT = 0,    /**< Left face (x-) */
    RIGHT,       /**< Right face (x+) */
    BOTTOM,      /**< Bottom face (y-) */
    TOP,         /**< Top face (y+) */
    FRONT,       /**< Front face (z-) */
    BACK,        /**< Back face (z+) */
    NUM_FACES    /**< Total number of faces */
} Face;


//--------------------------------------------------------------------------------
//                 3. PARTICLE LOCATION SYSTEM STRUCTS
//--------------------------------------------------------------------------------

/** @brief Defines a 3D axis-aligned bounding box. */
typedef struct {
    Cmpnts min_coords; ///< Minimum x, y, z coordinates of the bounding box.
    Cmpnts max_coords; ///< Maximum x, y, z coordinates of the bounding box.
} BoundingBox;

/** @brief Defines the vertices of a single hexahedral grid cell. */
typedef struct {
    Cmpnts vertices[8]; ///< Coordinates of the eight vertices of the cell.
} Cell;

/** @brief Defines a particle's core properties for Lagrangian tracking. */
typedef struct {
    PetscInt64 PID;     ///< Unique Particle ID.
    PetscInt cell[3];   ///< Computational indices (i, j, k) of the cell containing the particle.
    Cmpnts loc;         ///< Physical location (x,y,z) of the particle.
    Cmpnts vel;         ///< Physical velocity (vx,vy,vz) of the particle.
    Cmpnts weights;     ///< Interpolation weights within its host cell.

    ParticleLocationStatus location_status;  ///< Current state in the location/migration process.
    PetscMPIInt destination_rank;            ///< Target rank for migration (only valid if status MIGRATING_OUT)
} Particle;

/** @brief Stores the MPI ranks of neighboring subdomains. */
typedef struct {
    PetscMPIInt rank_xm, rank_xp; // Neighbors at -x, +x
    PetscMPIInt rank_ym, rank_yp; // Neighbors at -y, +y
    PetscMPIInt rank_zm, rank_zp; // Neighbors at -z, +z
} RankNeighbors;

/**
 * @struct RankCellInfo
 * @brief A lean struct to hold the global cell ownership range for a single MPI rank.
 * This is used to build a map of the entire domain's decomposition, allowing any
 * rank to quickly look up which rank owns a specific global cell index.
 */
typedef struct {
    PetscInt xs_cell, ys_cell, zs_cell; // Global starting cell indices (inclusive)
    PetscInt xm_cell, ym_cell, zm_cell; // Number of owned cells in each direction
} RankCellInfo;


/** @brief Information needed to migrate a single particle between MPI ranks. */
typedef struct {
    PetscInt local_index;
    PetscInt target_rank;
} MigrationInfo;


//--------------------------------------------------------------------------------
//                 4. BOUNDARY CONDITION SYSTEM ENUMS
//--------------------------------------------------------------------------------

/** @brief Identifies the six logical faces of a structured computational block. */
typedef enum {
    BC_FACE_NEG_X = 0, BC_FACE_POS_X = 1,
    BC_FACE_NEG_Y = 2, BC_FACE_POS_Y = 3,
    BC_FACE_NEG_Z = 4, BC_FACE_POS_Z = 5
} BCFace;

/** @brief Defines the general mathematical/physical category of a boundary. */
typedef enum {
    UNDEFINED = 0, NOGRAD, WALL, SYMMETRY, INLET,
    OUTLET, FARFIELD, PERIODIC, INTERFACE
} BCType;

/** @brief Defines the specific computational "strategy" for a boundary handler. */
typedef enum {
    BC_HANDLER_UNDEFINED = 0,BC_HANDLER_NOGRAD_COPY_GHOST,
    BC_HANDLER_WALL_NOSLIP, BC_HANDLER_WALL_MOVING,
    BC_HANDLER_SYMMETRY_PLANE,BC_HANDLER_INLET_PARABOLIC,
    BC_HANDLER_INLET_CONSTANT_VELOCITY, BC_HANDLER_INLET_PULSANTILE_FLUX, BC_HANDLER_INLET_DEVELOPED_PROFILE,
    BC_HANDLER_OUTLET_CONSERVATION, BC_HANDLER_OUTLET_PRESSURE,
    BC_HANDLER_FARFIELD_NONREFLECTING,
    BC_HANDLER_PERIODIC, BC_HANDLER_INTERFACE_OVERSET
} BCHandlerType;


//--------------------------------------------------------------------------------
//               5. BOUNDARY CONDITION SYSTEM STRUCTS
//--------------------------------------------------------------------------------

/** @brief A node in a linked list for storing key-value parameters from the bcs.dat file. */
struct BC_Param_s {
    char *key;
    char *value;
    struct BC_Param_s *next;
};

/** @brief Provides execution context for a boundary condition handler. */
typedef struct {
    UserCtx  *user;             ///< Access to all global simulation data.
    BCFace    face_id;          ///< The geometric face (0-5) being processed.
    const PetscReal *global_inflow_sum;  ///< Pointer to total domain inflow for the current step.
    const PetscReal *global_outflow_sum; ///< Pointer to total measured domain outflow for the current step.
} BCContext;

/** @brief The "virtual table" struct for a boundary condition handler object. */
struct BoundaryCondition_s {
    BCHandlerType type;
    void         *data;
    PetscErrorCode (*Initialize)(BoundaryCondition *self, BCContext *ctx);
    PetscErrorCode (*PreStep)(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow, PetscReal *local_outflow);
    PetscErrorCode (*Apply)(BoundaryCondition *self, BCContext *ctx);
    PetscErrorCode (*PlaceSource)(BoundaryCondition *self, BCContext *ctx, ...);
    PetscErrorCode (*Destroy)(BoundaryCondition *self);
};

/** @brief Holds the complete configuration for one of the six boundary faces. */
struct BoundaryFaceConfig_s {
    BCFace             face_id;
    BCType             mathematical_type;
    BCHandlerType      handler_type;
    BC_Param           *params;
    BoundaryCondition *handler;
};

//--------------------------------------------------------------------------------
//               6. IMBERSED BOUNDARY METHOD (IBM) & FSI STRUCTS
//--------------------------------------------------------------------------------

/**
 * @brief Holds interpolation and distance information for a single point
 *        in the Immersed Boundary Method (IBM).
 */
typedef struct {
  // --- Primary Interpolation Data (Grid-to-Surface) ---
  PetscInt	i1, j1, k1;   ///< Grid indices of the 1st interpolation stencil point.
  PetscInt	i2, j2, k2;   ///< Grid indices of the 2nd interpolation stencil point.
  PetscInt	i3, j3, k3;   ///< Grid indices of the 3rd interpolation stencil point.
  PetscReal	cr1, cr2, cr3;///< Interpolation coefficients for the stencil points.
  PetscReal	d_i;          ///< Distance to the interception point on grid cells.
  PetscInt	imode;        ///< Interception mode/scheme identifier.

  // --- Shortest Distance Data (Point-to-Surface) ---
  PetscInt	ni, nj, nk;   ///< Grid indices of the closest point on the solid surface.
  PetscReal	d_s;          ///< Shortest distance to a solid surface.
  Cmpnts	pmin;         ///< Coordinates of the closest point.
  PetscInt	cell;         ///< Index of the cell containing the closest point.
  PetscReal	cs1, cs2, cs3;///< Coefficients related to the closest surface point.

  // --- Secondary Interpolation Data Set 1 ---
  PetscInt	i11, j11, k11, i22, j22, k22, i33, j33, k33;
  PetscReal	cr11, cr22, cr33;
  PetscReal	d_ii;
  PetscInt	iimode;
  PetscReal	cs11, cs22, cs33;

  // --- Secondary Interpolation Data Set 2 ---
  PetscInt	ii1, jj1, kk1, ii2, jj2, kk2, ii3, jj3, kk3;
  PetscReal	ct1, ct2, ct3;
  PetscInt	smode;

  // --- Secondary Interpolation Data Set 3 ---
  PetscInt	ii11, jj11, kk11, ii22, jj22, kk22, ii33, jj33, kk33;
  PetscReal	ct11, ct22, ct33;
  PetscReal	d_ss;
  PetscInt	ssmode;
} IBMInfo;

/** @brief Represents a collection of nodes forming a surface for the IBM. */
typedef struct {
  PetscInt	nbnumber;     ///< Body number identifier.
  PetscInt	n_v;          ///< Number of vertices in the surface mesh.
  PetscInt	n_elmt;       ///< Number of elements (triangles) in the surface mesh.
  PetscInt	*nv1, *nv2, *nv3; ///< Node indices defining each triangle.
  PetscReal	*nf_x, *nf_y, *nf_z; ///< Normal vector for each element.
  PetscReal	*x_bp, *y_bp, *z_bp;   ///< Current coordinates of surface nodes.
  PetscReal	*x_bp0, *y_bp0, *z_bp0; ///< Initial coordinates of surface nodes.
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o; ///< Old (previous step) coordinates.
  PetscReal *cf,*cfsum,*CCp;
  PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270];
  Cmpnts	*u, *uold, *urm1;
  PetscReal     *dA ;         ///< Area of an element.
  PetscReal     *nt_x, *nt_y, *nt_z; ///< Tangent direction vector.
  PetscReal     *ns_x, *ns_y, *ns_z; ///< Azimuthal direction vector.
  PetscReal     *cent_x,*cent_y,*cent_z;
  PetscReal      *pres, *tau0, *tauN;
  PetscReal      *Bvel_u, *Bvel_v, *Bvel_w;
  PetscReal       x_min,x_max,y_min,y_max,z_min,z_max;
  Cmpnts *qvec;
  PetscReal *radvec;
} IBMNodes;


/** @brief Represents a collection of nodes forming a volume for the IBM. */
typedef struct {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;
  PetscInt	*nv1, *nv2, *nv3, *nv4;
  PetscReal	*x_bp, *y_bp, *z_bp;
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
  Cmpnts	*u, *uold, *urm1;
  PetscReal      V;
  PetscReal     *dV0 ;
  PetscReal     *cent_x,*cent_y,*cent_z;
  PetscReal     x_c,y_c,z_c;
  PetscReal     J[3][3];
  PetscReal     I_inv[3][3];
} IBMVNodes;

/** @brief A generic C-style linked list node for integers. */
typedef struct node{
  PetscInt Node;
  struct node *next;
} node;

/** @brief Head of a generic C-style linked list. */
typedef struct list{
  node *head;
} List;

/** @brief A generic C-style linked list node for indices. */
typedef struct list_node {
  PetscInt	index;
  struct list_node *next;
} Node_List;

/** @brief A linked list node containing IBM interpolation info. */
typedef struct IBMListNode {
  IBMInfo ibm_intp;
  struct IBMListNode* next;
} IBMListNode;

/** @brief Head of a linked list for IBM data. */
typedef struct IBMList {
  IBMListNode *head;
} IBMList;

/** @brief Information about a surface element for FSI calculations. */
typedef struct {
  PetscReal     P;    //Press on the surface elmt
  PetscInt      n_P; //number of Press Pts on the elmt
  PetscReal     Tow_ws, Tow_wt; //wall shear stress of the elmt
  PetscReal     Tow_wn; // normal stress
  PetscInt      Clsnbpt_i,Clsnbpt_j,Clsnbpt_k; //Closest Near Bndry Pt to each surf elmt
  PetscInt      icell,jcell,kcell;
  PetscInt      FoundAroundcell;
  PetscInt      Need3rdPoint;
} SurfElmtInfo;

/** @brief Holds all data related to the state and motion of a body in FSI. */
typedef struct {
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
  PetscReal    F_x_old,F_y_old,F_z_old; //Forces & Area
  PetscReal    F_x_real,F_y_real,F_z_real; //Forces & Area
  PetscReal    M_x,M_y,M_z; // Moments
  PetscReal    M_x_old,M_y_old,M_z_old; //Forces & Area
  PetscReal    M_x_real,M_y_real,M_z_real; //Forces & Area
  PetscReal    M_x_rm2,M_y_rm2,M_z_rm2; //Forces & Area
  PetscReal    M_x_rm3,M_y_rm3,M_z_rm3; //Forces & Area
  PetscReal    x_c,y_c,z_c; // center of rotation(mass)
  PetscReal    a_c[3];  //initial center of rotataion
  PetscReal    Mdpdn_x, Mdpdn_y,Mdpdn_z;
  PetscReal    Mdpdn_x_old, Mdpdn_y_old,Mdpdn_z_old;
  PetscReal    Power; //power
  PetscReal    clone;
  PetscInt     pbc[3];
  PetscReal    I_inv[3][3];
  PetscReal    L_n[3],L_o[3],L_r[3]; //Angular momentum
  PetscReal    alpha[3];
  PetscReal    acc[3];
  PetscReal    R[3][3],q[4],q_r[4];
  PetscReal    dS[6],dS_o[6],atk,atk_o; // Aitkin's iteration
  SurfElmtInfo  *elmtinfo;
  IBMInfo       *fsi_intp;
  PetscReal    Max_xbc,Min_xbc;
  PetscReal    Max_ybc,Min_ybc;
  PetscReal    Max_zbc,Min_zbc;
  PetscInt     CV_ys,CV_ye,CV_zs,CV_ze;
} FSInfo;

/** @brief Defines prescribed body motion using splines. */
typedef struct {
  PetscInt   n_time, n_midp, n_subit;
  PetscReal  *x_midp, *y_midp, *z_midp;
  PetscReal  *x_com, *y_com, *head_ang;
  PetscReal  *s1,*s2,*s3;
  PetscReal  *st1,*st2,*st3;
  Mat        Mphi;
  PetscReal  xmin,xmax,ymin,ymax,zmin,zmax;
} Cstart;

//--------------------------------------------------------------------------------
//               7. MULTIGRID AND POST-PROCESSING STRUCTS
//--------------------------------------------------------------------------------

/** @brief Context for Multigrid operations. */
typedef struct {
  UserCtx  *user;      ///< Pointer to the main application context for this level.
  PetscInt thislevel;  ///< The current grid level (0 is finest).
  DM       packer;     ///< DMDA for packing/unpacking fields.
} MGCtx;

/** @brief User-level context for managing the entire multigrid hierarchy. */
typedef struct {
  PetscInt  mglevels;   ///< Total number of multigrid levels.
  PetscInt  thislevel;  ///< Current level being operated on.
  PetscBool isc, jsc, ksc; ///< Local starting indices (legacy?).
  MGCtx     *mgctx;     ///< Array of multigrid contexts, one for each level.
  DM        packer;     ///< DM for field packing.
  SNES      snespacker; ///< SNES context for the packed system.
} UserMG;

/** @brief Parameters controlling post-processing output actions. */
typedef struct {
    PetscInt  startTime;       /* start of time loop */
    PetscInt  endTime;         /* end of time loop */
    PetscInt  timeStep;        /* increment of time loop */
    char      eulerianExt[8];  /* "vts", "vtk", etc. for grid data */
    char      particleExt[8];  /* "vtp", "vtk", etc. for particle data */
    char      eulerianPrefix[20]; /* directory to save euler fields in */
    char      particlePrefix[20]; /* directory to save particle fields in */
    PetscBool outputParticles; /* whether to write particle data or not */
} PostProcessParams;

/** @brief Enumerates the type of VTK file to be written. */
typedef enum {
    VTK_STRUCTURED, ///< For .vts (StructuredGrid) files from Eulerian data.
    VTK_POLYDATA    ///< For .vtp (PolyData) files from Lagrangian/particle data.
} VTKFileType;

/**
 * @brief Metadata structure for VTK output (.vts or .vtp files).
 * This struct is used for both Eulerian grid data and Lagrangian (particle) data.
 */
typedef struct _n_VTKMetaData {
    VTKFileType  fileType;      ///< Type of VTK file to generate.
    PetscInt          mx, my, mz;   ///< Grid dimensions (for structured).
    PetscInt          nnodes;       ///< Total number of structured grid nodes.
    PetscInt          npoints;      ///< Number of particle points (for polydata).
    PetscScalar      *coords;       ///< Coordinates array (3 components per point/node).
    PetscScalar      *scalarField;      ///< Pointer to scalar data array.
    const char       *scalarFieldName;  ///< Name of the scalar field.
    PetscInt          numScalarFields;  ///< Number of scalar fields.
    PetscScalar      *vectorField;      ///< Interleaved x,y,z vector field data.
    const char       *vectorFieldName;  ///< Name of the vector field (e.g., "Velocity").
    PetscInt          numVectorFields;  ///< Number of vector fields.
    PetscInt         *connectivity; ///< Connectivity array for cells/lines.
    PetscInt         *offsets;      ///< Offsets array for variable-size cells.
} VTKMetaData;


/*================================================================================*
 *                        MAIN APPLICATION CONTEXT                                *
 *================================================================================*/

/**
 * @brief User-defined context containing all simulation data, parameters, and
 *        configurations. This is the primary data structure passed throughout
 *        the application.
 */
typedef struct UserCtx {
    // --- Grid, Geometry & Parallelization ---
    DM da;                      ///< DMDA for scalar fields (e.g., Pressure, Nvert).
    DM fda;                     ///< DMDA for primary 3-component vector fields (e.g., Velocity).
    DM fda2;                    ///< DMDA for 2-component vector fields (e.g., k-omega).
    DMDALocalInfo info;         ///< Cached local grid info for the current rank.
    AO ao;                      ///< Application Ordering context for renumbering.
    PetscInt IM, JM, KM;        ///< Global grid dimensions (number of cells in i,j,k).
    PetscReal Max_X, Max_Y, Max_Z; ///< Physical maximum bounds of the grid.
    PetscReal Min_X, Min_Y, Min_Z; ///< Physical minimum bounds of the grid.
    PetscReal xMin,yMin,zMin;   ///< Physical minimum bounds of the grid.
    PetscReal xMax,yMax,zMax;   ///< Physical maximum bounds of the grid.
    PetscReal rx, ry, rz;       ///< Grid stretching ratios.
    PetscInt nblk;              ///< Number of grid blocks in the simulation.
    BoundingBox bbox;           ///< Bounding box for the local processor's grid domain.
    BoundingBox global_domain_bbox; ///< Bounding box for the entire global domain.
    RankNeighbors neighbors;    ///< MPI ranks of neighboring subdomains.
    PetscInt GridOrientation; ///< Integerer that determines Whether grid is right-handed (+1) or left-handed (-1) i.e normals are facing outward or inward.
    PetscInt isc, jsc, ksc;     ///< Starting corner indices for the local subdomain.
    PetscInt _this;             ///< Legacy block index, for multi-block contexts.
    PetscInt cgrid, df;         ///< Grid flags (e.g., coarse grid, degrees of freedom).
    PetscErrorCode aotopetsc;   ///< Legacy error code for AO to PETSc mapping.
    PetscBool assignedA;        ///< Flag indicating if matrix A has been assigned.

    // --- Simulation Parameters & State ---
    PetscReal ren;              ///< Reynolds number.
    PetscReal dt;               ///< Time step size.
    PetscReal st;               ///< Strouhal number.
    PetscReal cfl;              ///< Courant-Friedrichs-Lewy number.
    PetscReal vnn;              ///< Von Neumann number.
    PetscReal ratio;            ///< A generic ratio parameter.
    PetscInt step;              ///< Current timestep index.
    PetscReal ti;               ///< Current simulation time.
    PetscReal cdisx, cdisy, cdisz; ///< Central differencing dissipation parameters.
    PetscInt FieldInitialization; ///< Flag controlling how fields are initially set.
    PetscInt LoggingFrequency;  ///< Frequency for detailed logging.

    // --- Boundary Conditions & Interfacing ---
    BCS Bcs;                    ///< Legacy boundary condition structure.
    PetscInt bctype[6];         ///< Legacy integer array for BC types per face.
    PetscInt inttype[7];        ///< Legacy integer array for interface types.
    PetscInt itfcptsnumber;     ///< Number of interface points.
    PetscInt *itfcI, *itfcJ, *itfcK; ///< Indices of interface points.
    PetscInt *itfchostI, *itfchostJ, *itfchostK, *itfchostB; ///< Host cell data for interpolation.
    PetscReal *itfchostx, *itfchosty, *itfchostz; ///< Host cell coordinates.
    PetscInt ip[10],jp[10],kp[10],dispx[10],dispy[10],dispz[10],dispnn[10]; ///< Legacy probing/interface data.
    PetscInt *idx_from;         ///< Legacy index mapping data.
    BoundaryFaceConfig boundary_faces[6]; ///< The new, primary BC configuration array.
    BCType    face_bc_types[6];
    PetscBool inletFaceDefined;           ///< Legacy flag for particle system compatibility.
    BCFace    identifiedInletBCFace;      ///< Legacy field for particle system compatibility. 

    // --- Primary Flow Fields (Global & Local Views) ---
    Vec Ucont, lUcont; Vec Vcont, lVcont; Vec Ucat, lUcat; Vec Wcat, lWcat;
    Vec P, lP; Vec Phi, lPhi; Vec Nvert, lNvert; Vec GridSpace, lGridSpace;

    // --- Time-Stepping & Solver Workspace Fields ---
    Vec Ucont_o, lUcont_o; Vec Ucat_o; Vec P_o; Vec Nvert_o, lNvert_o;
    Vec Ucont_rm1, lUcont_rm1; Vec Rhs; Vec dUcont; Vec pUcont; Vec DUold;
    Vec Forcing; Vec Dt; Vec psuedot; Vec lUstar;

    // --- Pressure-Poisson System ---
    Mat A, C; Mat MR, MP; KSP ksp; MatNullSpace nullsp;
    Vec B; Vec Rhsp; Vec X, R;
    PetscBool multinullspace; PetscInt *KSKE;
    int	rhs_count; Vec Gid, Gidm; Vec Phi2, B2, Ucont2;
    int	local_Phi2_size, p_global_begin, reduced_p_size;

    // --- Grid Metrics (Global and Local) ---
    Vec Cent; Vec Centx, Centy, Centz; Vec lCent;
    Vec Csi, Eta, Zet, Aj; Vec ICsi, IEta, IZet, IAj;
    Vec JCsi, JEta, JZet, JAj; Vec KCsi, KEta, KZet, KAj;
    Vec lCsi, lEta, lZet, lAj; Vec lICsi, lIEta, lIZet, lIAj;
    Vec lJCsi, lJEta, lJZet, lJAj; Vec lKCsi, lKEta, lKZet, lKAj;
    Vec lNFace; Vec lMAreaCsi, lMAreaEta, lMAreaZet;

    // --- Turbulence Modeling (LES/RANS) ---
    Vec Nu_t, lNu_t; Vec CS, lCs; Vec K_Omega, lK_Omega;
    Vec K_Omega_o, lK_Omega_o; Vec Distance; Vec lF1; Cmpnts2 **komega_plane;

    // --- Statistical Averaging ---
    Vec Ucat_sum; Vec Ucat_cross_sum; Vec Ucat_square_sum; Vec P_sum; Vec USQ;

    // --- Immersed Boundary Method (IBM) ---
    PetscInt ibmnumber; IBMInfo  *ibm_intp; IBMNodes *ibm; IBMList *ibmlist;

    // --- Multigrid Hierarchy ---
    PetscInt thislevel; PetscInt mglevels; DM *da_f, *da_c;
    UserCtx *user_f, *user_c; Vec Ucont_MG; Vec *lNvert_c;

    // --- Nonlinear Solver (SNES) ---
    SNES snes; Mat J, PJ, J2; Vec RFC;
    ISColoring iscoloring; MatFDColoring matfdcoloring;

    // --- Transport of a Scalar 'Q' ---
    Vec Qnew; Vec Ql; Vec shrr; Vec ItfcQ, lItfcQ; Vec nhostQ; PetscReal LVOUT;

    // --- Inflow & General Interface Data ---
    FlowWave *inflow; PetscInt number_flowwave; Vec inletU;
    PetscReal r[101], tin[101], uinr[101][1001];
    Vec Itfc, lItfc; Vec lItfcKO; Vec nhostU;

    // --- Diagnostics & Post-Processing ---
    PetscReal utaw, utaw1; PetscReal utawsum, utaw1sum;
    PetscReal FluxInSum, FluxOutSum, FluxIntfcSum, FluxIntpSum;
    PetscReal AreaOutSum, AreaIntfcSum;
    Vec lSx, lSy, lSz, lS; Vec lLM, lMM, lNM;
    PostProcessParams pp;       ///< Parameters controlling post-processing.

    // --- Particle System ---
    DM swarm;                   ///< DMSwarm object for particle data.
    PetscMPIInt *miglist;      ///< List of ranks for particle migration.
    RankCellInfo *RankCellInfoMap; ///< List of RankCellInfo objects which contain the cells owned by each rank.
    PetscInt NumberofParticles; ///< Total number of particles in the simulation.
    Vec ParticleCount;          ///< Eulerian field to count particles per cell.
    PetscInt ParticleInitialization; ///< Flag controlling how particles are initially placed.
  
    // --- Modern System Additions ---
    PetscInt n_ranks;           ///< Total number of MPI ranks.

    // --- Initial Condition Parameters ---.
    Cmpnts    InitialConstantContra;   ///< A constant contravariant velocity vector.
    PetscReal InitialConstantPressure; ///< A constant pressure value.
    PetscReal InitialConstantNvert;    ///< A constant nvert value.
    
    // --- Simulation Mode Flags ---
    PetscBool averaging;        ///< Flag to enable statistical averaging.
    PetscBool les;              ///< Flag to enable Large Eddy Simulation.
    PetscBool rans;             ///< Flag to enable Reynolds-Averaged Navier-Stokes.
  
} UserCtx;

/*================================================================================*
 *                    FUNCTION PROTOTYPES                                         *
 *================================================================================*/

// --- Linked List Helpers ---
void initlist(List *ilist);
void insertnode(List *ilist, PetscInt Node);
void destroy(List *ilist);
void InitIBMList(IBMList *ilist);
void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp);
void DestroyIBMList(IBMList *ilist);

// --- IBM: Geometry, Search, and Interpolation ---
PetscErrorCode ibm_read_tecplot(IBMNodes *ibm, PetscInt nt, PetscBool flg);
PetscErrorCode ibm_read_ucd(IBMNodes *ibm, PetscInt ibi);
PetscErrorCode ibm_search_advanced(UserCtx *user, IBMNodes *ibm, PetscInt ibi);
PetscErrorCode ibm_interpolation_advanced(UserCtx *user, IBMNodes *ibm, PetscInt ibi, PetscInt Add_dUndt);

// --- FSI, Body Motion, and Collision ---
PetscErrorCode FsiInitialize(PetscInt n_elmt, FSInfo *fsi, PetscInt ibi);
PetscErrorCode Struc_Solver(UserMG *usermg, IBMNodes *ibm, FSInfo *fsi, Cstart *cstart, PetscInt itr_sc, PetscInt tistart, PetscBool *DoSCLoop);
PetscErrorCode Calc_forces_SI(FSInfo *FSinfo, UserCtx *user, IBMNodes *ibm, PetscInt ti, PetscInt ibi, PetscInt bi);
PetscErrorCode CollisionDetectionOfParticles(UserCtx *user, IBMNodes *ibm,FSInfo *fsi);

// --- Solvers and Time Integration ---
PetscErrorCode Flow_Solver(UserMG *usermg, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode ImplicitMomentumSolver(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Implicit_SNES(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);
PetscErrorCode Projection(UserCtx *user);
void PoissonSolver_Hypre(UserCtx *user, IBMNodes *ibm, IBMInfo *ibminfo);
PetscErrorCode PoissonSolver_MG(UserMG *usermg);

// --- Turbulence Modeling ---
void K_Omega_IC(UserCtx *user);
void Solve_K_Omega(UserCtx *user);
void Compute_Smagorinsky_Constant_1(UserCtx *user, Vec Ucont, Vec Ucat);
void wall_function(UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, double nx, double ny, double nz);

// --- Numerical Kernels & Discretization ---
void Compute_du_dxyz (double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc, double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz, double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz );
PetscErrorCode FormFunction1(UserCtx *user, Vec Rhs);
PetscErrorCode Divergence(UserCtx *user);
//PetscErrorCode Contra2Cart(UserCtx *user);

// --- Grid, Metrics, and Boundary Conditions ---
PetscErrorCode FormBCS(UserCtx *user);
PetscErrorCode Block_Interface_U(UserCtx *user);
PetscErrorCode InletRead(UserCtx *user);

// --- Multigrid ---
PetscErrorCode MG_Initial(UserMG *usermg, IBMNodes *ibm);
PetscErrorCode MG_Finalize(UserMG *usermg);
PetscErrorCode Implicit_DMMG(UserMG *usermg);

// --- I/O and Diagnostics ---
PetscErrorCode ibm_surface_VTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti);
PetscErrorCode Ucont_P_Binary_Output(UserCtx *user);
PetscErrorCode KE_Output(UserCtx *user);

// --- Add prototypes for new systems here (Particle location, BC handlers, etc.) ---


#ifdef __cplusplus
}
#endif

#endif // VARIABLES_H
