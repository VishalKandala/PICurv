/**
 * @file variables.h
 * @author Vishal Kandala
 * @brief Main header file for a complex fluid dynamics solver.
 *
 * This file defines the core data structures, global variables, and function
 * prototypes for a CFD application. It combines features for a curvilinear,
 * multi-block, immersed boundary (IBM) solver with fluid-structure interaction (FSI),
 * a modern particle tracking system, a modular boundary condition framework,
 * and post-processing utilities. It heavily utilizes the PETSc library.
 *
 * *** REVISION NOTE ***
 * This version introduces a central `SimulationContext` struct to encapsulate all
 * simulation-wide configuration, replacing the legacy system of global variables.
 * The `UserCtx` struct has been refactored to hold only data specific to a
 * single grid level and now contains a back-pointer to the `SimulationContext`.
 * All struct definitions now follow the `typedef struct Name { ... } Name;` convention.
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
#include "petscdmswarm.h"

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

/* --- NOTE: All legacy global 'extern' variables have been removed. --- */
/* They are now members of the SimulationContext struct. */


/*================================================================================*
 *                     CORE DATA STRUCTURES & TYPEDEFS                            *
 *================================================================================*/

//--------------------------------------------------------------------------------
//                 1. FORWARD DECLARATIONS & BASIC TYPES
//--------------------------------------------------------------------------------

// --- Forward Declarations ---
// These declarations allow pointers to these types before their full definition.
typedef struct SimCtx SimCtx;
typedef struct UserCtx UserCtx;
typedef struct BC_Param_s BC_Param; /* Retains _s for linked list safety */
typedef struct BoundaryCondition BoundaryCondition;
typedef struct BoundaryFaceConfig BoundaryFaceConfig;
typedef struct IBMNodes IBMNodes;
typedef struct IBMVNodes IBMVNodes;
typedef struct FSInfo FSInfo;
typedef struct UserMG UserMG;
typedef struct BoundingBox BoundingBox;
typedef struct Cell Cell;
typedef struct Particle Particle;
typedef struct RankNeighbors RankNeighbors;
typedef struct RankCellInfo RankCellInfo;
typedef struct MigrationInfo MigrationInfo;
typedef struct Cstart Cstart;
typedef struct MGCtx MGCtx;
typedef struct PostProcessParams PostProcessParams;
typedef struct VTKFieldInfo VTKFieldInfo;
typedef struct VTKMetaData VTKMetaData;
typedef struct IBMInfo IBMInfo;
typedef struct SurfElmtInfo SurfElmtInfo;
typedef struct ScalingCtx; 

// --- Foundational Geometric and Data Types ---

/** @brief A 3D point or vector with PetscScalar components. */
typedef struct Cmpnts {
    PetscScalar x, y, z;
} Cmpnts;

/** @brief A 2D point or vector with PetscScalar components. */
typedef struct Cmpnts2 {
    PetscScalar x, y;
} Cmpnts2;

/** @brief A 2D vector of PETSc real numbers (for geometry/coordinates). */
typedef struct Cpt2D {
  PetscReal	x, y;
} Cpt2D;

/** @brief Represents a single point in a time-varying flow waveform. */
typedef struct FlowWave {
    PetscReal t, f;
} FlowWave;

/** @brief Legacy boundary condition data. Likely superseded by the newer BC system. */
typedef struct BCS {
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
 */
typedef enum {
    NEEDS_LOCATION,
    ACTIVE_AND_LOCATED,
    MIGRATING_OUT,
    LOST,
    UNINITIALIZED
} ParticleLocationStatus;

/** @brief Enumerates the six faces of a cubic cell for distance calculations. */
typedef enum {
    LEFT = 0, RIGHT, BOTTOM, TOP, FRONT, BACK, NUM_FACES
} Face;


//--------------------------------------------------------------------------------
//                 3. PARTICLE LOCATION SYSTEM STRUCTS
//--------------------------------------------------------------------------------

/** @brief Defines a 3D axis-aligned bounding box. */
typedef struct BoundingBox {
    Cmpnts min_coords; ///< Minimum x, y, z coordinates of the bounding box.
    Cmpnts max_coords; ///< Maximum x, y, z coordinates of the bounding box.
} BoundingBox;

/** @brief Defines the vertices of a single hexahedral grid cell. */
typedef struct Cell {
    Cmpnts vertices[8]; ///< Coordinates of the eight vertices of the cell.
} Cell;

/** @brief Defines a particle's core properties for Lagrangian tracking. */
typedef struct Particle {
    PetscInt64 PID;
    PetscInt cell[3];
    Cmpnts loc;
    Cmpnts vel;
    Cmpnts weights;
    ParticleLocationStatus location_status;
    PetscMPIInt destination_rank;
} Particle;

/** @brief Stores the MPI ranks of neighboring subdomains. */
typedef struct RankNeighbors {
    PetscMPIInt rank_xm, rank_xp;
    PetscMPIInt rank_ym, rank_yp;
    PetscMPIInt rank_zm, rank_zp;
} RankNeighbors;

/** @brief A lean struct to hold the global cell ownership range for a single MPI rank. */
typedef struct RankCellInfo {
    PetscInt xs_cell, ys_cell, zs_cell;
    PetscInt xm_cell, ym_cell, zm_cell;
} RankCellInfo;

/** @brief Information needed to migrate a single particle between MPI ranks. */
typedef struct MigrationInfo {
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

/** @brief Defines the general mathematical/physical Category of a boundary. */
typedef enum {
    WALLFUNCTION       = -1,
    INTERFACE          = 0,
    WALL               = 1 , 
    MOVING_WALL        = 2,
    SYMMETRY           = 3,
    OUTLET             = 4,
    INLET              = 5, 
    FARFIELD           = 6,
    PERIODIC           = 7,
    CHARACTERISTIC_BC  = 8,
    ANALYTICAL_VORTEX  = 9,
    JUNCTION           = 10,
    ANNULUS            = 11,
    OGRID              = 12,
    RHEOLOGY           = 13,
    // Note:  Legacy 14 can be a JUNCTION with a specific handler.
} BCType;

/** @brief Defines the specific computational "strategy" for a boundary handler. */
typedef enum {
    BC_HANDLER_UNDEFINED = 0,
    BC_HANDLER_WALL_NOSLIP = 1, 
    BC_HANDLER_WALL_MOVING = 2,
    BC_HANDLER_SYMMETRY_PLANE = 3,
    BC_HANDLER_INLET_CONSTANT_VELOCITY = 4,
    BC_HANDLER_INLET_PARABOLIC = 5,
    BC_HANDLER_INLET_PROFILE_FROM_FILE = 6,
    BC_HANDLER_INLET_INTERP_FROM_FILE = 7,
    BC_HANDLER_INLET_PULSATILE_FLUX = 8,
    BC_HANDLER_FARFIELD_NONREFLECTING = 9,
    BC_HANDLER_OUTLET_CONSERVATION = 10, 
    BC_HANDLER_OUTLET_PRESSURE = 11,
    BC_HANDLER_PERIODIC = 12, 
    BC_HANDLER_INTERFACE_OVERSET = 13
  } BCHandlerType;

  
typedef enum{
    BC_PRIORITY_UNDEFINED = -1,
    BC_PRIORITY_INLET = 0,
    BC_PRIORITY_FARFIELD = 1,
    BC_PRIORITY_WALL = 2,
    BC_PRIORITY_OUTLET =3
} BCPriorityType;

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
typedef struct BCContext {
    UserCtx  *user;
    BCFace    face_id;
    const PetscReal *global_inflow_sum;
    const PetscReal *global_farfield_inflow_sum;
    const PetscReal *global_farfield_outflow_sum;
    const PetscReal *global_outflow_sum;
} BCContext;

/** @brief The "virtual table" struct for a boundary condition handler object. */
typedef struct BoundaryCondition {
    BCHandlerType type;
    BCPriorityType priority;
    void         *data;
    PetscErrorCode (*Initialize)(BoundaryCondition *self, BCContext *ctx);
    PetscErrorCode (*PreStep)(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow, PetscReal *local_outflow);
    PetscErrorCode (*Apply)(BoundaryCondition *self, BCContext *ctx);
    PetscErrorCode (*PostStep)(BoundaryCondition *self, BCContext *ctx, ...);
    PetscErrorCode (*Destroy)(BoundaryCondition *self);
} BoundaryCondition;

/** @brief Holds the complete configuration for one of the six boundary faces. */
typedef struct BoundaryFaceConfig {
    BCFace             face_id;
    BCType             mathematical_type;
    BCHandlerType      handler_type;
    BC_Param           *params;
    BoundaryCondition *handler;
} BoundaryFaceConfig;

//--------------------------------------------------------------------------------
//               6. IMBERSED BOUNDARY METHOD (IBM) & FSI STRUCTS
//--------------------------------------------------------------------------------

/** @brief Holds interpolation and distance information for a single IBM point. */
typedef struct IBMInfo {
  PetscInt	i1, j1, k1, i2, j2, k2, i3, j3, k3;
  PetscReal	cr1, cr2, cr3;
  PetscReal	d_i;
  PetscInt	imode;
  PetscInt	ni, nj, nk;
  PetscReal	d_s;
  Cmpnts	pmin;
  PetscInt	cell;
  PetscReal	cs1, cs2, cs3;
  PetscInt	i11, j11, k11, i22, j22, k22, i33, j33, k33;
  PetscReal	cr11, cr22, cr33;
  PetscReal	d_ii;
  PetscInt	iimode;
  PetscReal	cs11, cs22, cs33;
  PetscInt	ii1, jj1, kk1, ii2, jj2, kk2, ii3, jj3, kk3;
  PetscReal	ct1, ct2, ct3;
  PetscInt	smode;
  PetscInt	ii11, jj11, kk11, ii22, jj22, kk22, ii33, jj33, kk33;
  PetscReal	ct11, ct22, ct33;
  PetscReal	d_ss;
  PetscInt	ssmode;
} IBMInfo;

/** @brief Represents a collection of nodes forming a surface for the IBM. */
typedef struct IBMNodes {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;
  PetscInt	*nv1, *nv2, *nv3;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscReal	*x_bp, *y_bp, *z_bp;
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
  PetscReal *cf, *cfsum, *CCp;
  PetscReal	x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270];
  Cmpnts	*u, *uold, *urm1;
  PetscReal *dA, *nt_x, *nt_y, *nt_z, *ns_x, *ns_y, *ns_z;
  PetscReal *cent_x, *cent_y, *cent_z;
  PetscReal *pres, *tau0, *tauN, *Bvel_u, *Bvel_v, *Bvel_w;
  PetscReal x_min, x_max, y_min, y_max, z_min, z_max;
  Cmpnts    *qvec;
  PetscReal *radvec;
} IBMNodes;

/** @brief Represents a collection of nodes forming a volume for the IBM. */
typedef struct IBMVNodes {
  PetscInt	nbnumber;
  PetscInt	n_v, n_elmt;
  PetscInt	*nv1, *nv2, *nv3, *nv4;
  PetscReal	*x_bp, *y_bp, *z_bp;
  PetscReal	*x_bp0, *y_bp0, *z_bp0;
  PetscReal	*x_bp_o, *y_bp_o, *z_bp_o;
  Cmpnts	*u, *uold, *urm1;
  PetscReal V, *dV0;
  PetscReal *cent_x, *cent_y, *cent_z;
  PetscReal x_c, y_c, z_c;
  PetscReal J[3][3], I_inv[3][3];
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
typedef struct SurfElmtInfo {
  PetscReal     P;
  PetscInt      n_P;
  PetscReal     Tow_ws, Tow_wt, Tow_wn;
  PetscInt      Clsnbpt_i, Clsnbpt_j, Clsnbpt_k;
  PetscInt      icell, jcell, kcell;
  PetscInt      FoundAroundcell, Need3rdPoint;
} SurfElmtInfo;

/** @brief Holds all data related to the state and motion of a body in FSI. */
typedef struct FSInfo {
  PetscReal    S_new[6], S_old[6], S_real[6], S_realm1[6];
  PetscReal    S_ang_n[6], S_ang_o[6], S_ang_r[6], S_ang_rm1[6];
  PetscReal    red_vel, damp, mu_s;
  PetscReal    F_x, F_y, F_z, A_tot;
  PetscReal    F_x_old, F_y_old, F_z_old;
  PetscReal    F_x_real, F_y_real, F_z_real;
  PetscReal    M_x, M_y, M_z;
  PetscReal    M_x_old, M_y_old, M_z_old;
  PetscReal    M_x_real, M_y_real, M_z_real;
  PetscReal    M_x_rm2, M_y_rm2, M_z_rm2;
  PetscReal    M_x_rm3, M_y_rm3, M_z_rm3;
  PetscReal    x_c, y_c, z_c;
  PetscReal    a_c[3];
  PetscReal    Mdpdn_x, Mdpdn_y, Mdpdn_z;
  PetscReal    Mdpdn_x_old, Mdpdn_y_old, Mdpdn_z_old;
  PetscReal    Power, clone;
  PetscInt     pbc[3];
  PetscReal    I_inv[3][3];
  PetscReal    L_n[3], L_o[3], L_r[3];
  PetscReal    alpha[3], acc[3];
  PetscReal    R[3][3], q[4], q_r[4];
  PetscReal    dS[6], dS_o[6], atk, atk_o;
  SurfElmtInfo *elmtinfo;
  IBMInfo      *fsi_intp;
  PetscReal    Max_xbc, Min_xbc, Max_ybc, Min_ybc, Max_zbc, Min_zbc;
  PetscInt     CV_ys, CV_ye, CV_zs, CV_ze;
} FSInfo;

/** @brief Defines prescribed body motion using splines. */
typedef struct Cstart {
  PetscInt   n_time, n_midp, n_subit;
  PetscReal  *x_midp, *y_midp, *z_midp;
  PetscReal  *x_com, *y_com, *head_ang;
  PetscReal  *s1, *s2, *s3, *st1, *st2, *st3;
  Mat        Mphi;
  PetscReal  xmin, xmax, ymin, ymax, zmin, zmax;
} Cstart;

//--------------------------------------------------------------------------------
//               7. MULTIGRID AND POST-PROCESSING STRUCTS
//--------------------------------------------------------------------------------

/** @brief Context for Multigrid operations. */
typedef struct MGCtx {
  UserCtx  *user;
  PetscInt thislevel;
  DM       packer;
} MGCtx;

/** @brief User-level context for managing the entire multigrid hierarchy. */
typedef struct UserMG {
  PetscInt  mglevels;
  PetscInt  thislevel;
  PetscBool isc, jsc, ksc;
  MGCtx     *mgctx;
  DM        packer;
  SNES      snespacker;
} UserMG;

#define MAX_PIPELINE_LENGTH 1024
#define MAX_FIELD_LIST_LENGTH 1024
#define MAX_FILENAME_LENGTH 256
/// Defines the maximum number of data fields for VTK point data.
#define MAX_POINT_DATA_FIELDS 20 
#define MAX_VTK_FIELD_NAME_LENGTH 64 ///< Maximum length for VTK field names.
/**
 * @brief Holds all configuration parameters for a post-processing run.
 *        This is an enhanced version combining command-line and file-based settings.
 */
typedef struct PostProcessParams {
    // -- Source Directory --- (For Data)
    char source_dir[PETSC_MAX_PATH_LEN];
  
    // --- Time Controls (can be set by command line or file) ---
    PetscInt startTime;
    PetscInt endTime;
    PetscInt timeStep;
    PetscBool outputParticles;

    // --- Configuration primarily from the .cfg file ---
    char process_pipeline[MAX_PIPELINE_LENGTH];
    char output_fields_instantaneous[MAX_FIELD_LIST_LENGTH];
    char output_fields_averaged[MAX_FIELD_LIST_LENGTH];
    char output_prefix[MAX_FILENAME_LENGTH];
    char particle_pipeline[MAX_PIPELINE_LENGTH];
    char particle_fields[MAX_FIELD_LIST_LENGTH];
    char particle_output_prefix[MAX_FILENAME_LENGTH];
    PetscInt particle_output_freq;

    // --- Legacy settings ---
    char eulerianExt[8]; // from original PostProcessParams (repurposed for PreCheckAndResize() as the input file extension.)
    char particleExt[8]; // from original PostProcessParams

    // --- Processing Parameters ---
    PetscInt reference[3]; // Reference point for normalizing any field against.
    
} PostProcessParams;

/**
 * @brief Stores all necessary information for a single data array in a VTK file.
 */
typedef struct VTKFieldInfo {
    char        name[MAX_VTK_FIELD_NAME_LENGTH]; // Name of the data field
    PetscInt    num_components; // 1 for scalar, 3 for vector
    PetscScalar* data;          // Pointer to the raw data array (must be freed by caller)
} VTKFieldInfo;

/** @brief Enumerates the type of VTK file to be written. */
typedef enum {
    VTK_STRUCTURED,
    VTK_POLYDATA
} VTKFileType;

typedef struct VTKMetaData {
    VTKFileType     fileType;
    PetscInt        mx, my, mz;
    PetscInt        npoints;
    PetscScalar*    coords;
    VTKFieldInfo    point_data_fields[MAX_POINT_DATA_FIELDS];
    PetscInt        num_point_data_fields;
    PetscInt*       connectivity;
    PetscInt*       offsets;
} VTKMetaData;

/**
 * @brief Defines the execution mode of the application.
 */
typedef enum {
    EXEC_MODE_SOLVER,          // The application is running as the main fluid solver.
    EXEC_MODE_POSTPROCESSOR,   // The application is running as the post-processor.
    EXEC_MODE_UNKNOWN          // Default/error state.
} ExecutionMode;

//-------------------------------------------------------------------------------
//             8. SCALING AND DIMENSIONAL ANALYSIS STRUCTS
//-------------------------------------------------------------------------------
typedef struct ScalingCtx{
  PetscReal L_ref;
  PetscReal U_ref;
  PetscReal rho_ref;
  PetscReal P_ref;
}ScalingCtx;
//-------------------------------------------------------------------------------
/*================================================================================*
 *                        MAIN APPLICATION CONTEXT                                *
 *================================================================================*/

/**
 * @brief The master context for the entire simulation.
 *
 * This struct encapsulates all global configuration flags, physical parameters,
 * simulation control settings, and top-level data objects. It replaces the
 * numerous global variables from the legacy codebase, providing a single,
 * explicit source of truth for the simulation's setup.
 */
typedef struct SimCtx {

    //================ Group 1: Parallelism & MPI Information ================
    PetscMPIInt rank;
    PetscMPIInt size;

    //================ Group 2: Simulation Control, Time, and I/O ================

    PetscInt  step;
    PetscReal ti;
    PetscInt  StartStep;
    PetscInt  StepsToRun;
    PetscInt  tiout;
    PetscReal StartTime;
    PetscReal dt;
    PetscBool OnlySetup;
    PetscViewer logviewer;
    PetscInt    OutputFreq;
    ExecutionMode exec_mode;
    char eulerianSource[64];
    char restart_dir[PETSC_MAX_PATH_LEN];
    char output_dir[PETSC_MAX_PATH_LEN];
    char euler_subdir[PETSC_MAX_PATH_LEN];
    char particle_subdir[PETSC_MAX_PATH_LEN];
    char log_dir[PETSC_MAX_PATH_LEN];
    char _io_context_buffer[PETSC_MAX_PATH_LEN]; // Persistent store for I/O context strings.
    char *current_io_directory; // Pointer into the above buffer.

    //================ Group 3: High-Level Physics & Model Selection Flags ================
    PetscInt  immersed, movefsi, rotatefsi, sediment, rheology;
    PetscInt  invicid, TwoD, thin, moveframe, rotateframe, blank;
    PetscInt  dgf_x, dgf_y, dgf_z, dgf_ax, dgf_ay, dgf_az;
  
    //================ Group 4: Specific Simulation Case Flags ================
    PetscInt  cop, fish, fish_c, fishcyl, eel, pizza, turbine, Pipe, wing, hydro, MHV, LV;

    //================ Group 5: Solver & Numerics Parameters ================
    PetscInt  implicit, implicit_type, imp_MAX_IT;
    PetscReal imp_atol, imp_rtol, imp_stol;
    PetscInt  mglevels,mg_MAX_IT, mg_idx, mg_preItr, mg_poItr;
    PetscInt  poisson;
    PetscReal poisson_tol;
    PetscInt  STRONG_COUPLING,central;
    PetscReal ren, st,cfl, vnn, cdisx, cdisy, cdisz;
    PetscInt  FieldInitialization; 
    Cmpnts    InitialConstantContra;
    

    //================ Group 6: Physical & Geometric Parameters ================
    PetscInt  NumberOfBodies;
    PetscReal Flux_in, angle,max_angle;
    PetscReal CMx_c, CMy_c, CMz_c;
    ScalingCtx scaling;
    PetscReal  wall_roughness_height;

    //================ Group 7: Grid, Domain, and Boundary Condition Settings ================
    PetscInt  block_number, inletprofile, grid1d, Ogrid, channelz;
    PetscInt  i_periodic, j_periodic, k_periodic, blkpbc, pseudo_periodic;
    PetscBool generate_grid;        
    PetscReal grid_rotation_angle; 
    PetscReal Croty, Crotz; 
    char grid_file[PETSC_MAX_PATH_LEN];
    PetscInt da_procs_x, da_procs_y, da_procs_z;
    PetscInt num_bcs_files;
    char **bcs_files;
    PetscReal  FluxInSum, FluxOutSum,Fluxsum,FarFluxInSum,FarFluxOutSum;
    PetscReal  AreaInSum, AreaOutSum;
    PetscReal  U_bc;
    PetscInt   ccc;
    PetscReal  ratio;
  
    //================ Group 8: Turbulence Modeling (LES/RANS) ================
    PetscInt  les, rans;
    PetscInt  wallfunction, mixed, clark, dynamic_freq;
    PetscReal max_cs,Const_CS;
    PetscInt  testfilter_ik, testfilter_1d, i_homo_filter, j_homo_filter, k_homo_filter;
    PetscBool  averaging;
  
    //================ Group 9: Particle / DMSwarm Data & Settings ================
    PetscInt  np;
    PetscBool readFields;
    DM        dm_swarm;
    BoundingBox *bboxlist;
    PetscInt   ParticleInitialization;
    char particleRestartMode[16];

    //================ Group 10: Immersed Boundary & FSI Data Object Pointers ================
    IBMNodes  *ibm;
    IBMVNodes *ibmv;
    FSInfo    *fsi;
    PetscBool rstart_fsi;
    PetscInt  duplicate;

    //================ Group 11: Top-Level Managers, Custom Configuration,Logging/Monitoring ====
    UserMG    usermg;
    char      allowedFile[PETSC_MAX_PATH_LEN];
    PetscBool useCfg;
    char      **allowedFuncs;
    PetscInt  nAllowed;
    PetscInt  LoggingFrequency;
    PetscReal summationRHS;
    PetscReal MaxDiv;
    PetscInt  MaxDivFlatArg, MaxDivx,MaxDivy,MaxDivz;
    // Profiling 
    char      criticalFuncsFile[PETSC_MAX_PATH_LEN];
    PetscBool useCriticalFuncsCfg;
    char      **criticalFuncs;
    PetscInt  nCriticalFuncs;

    //================ Group 12: Post-Processing =================================================
    char      PostprocessingControlFile[PETSC_MAX_PATH_LEN];
    PostProcessParams *pps;
    
   //=============== Group 13: Miscellaneous =============================================
   PetscReal	r[101], tin[101], uinr[101][1001];

} SimCtx;


/**
 * @brief User-defined context containing data specific to a single computational
 *        grid level. This is the primary data structure passed to low-level
 *        numerical and solver functions.
 */
typedef struct UserCtx {

    // --- The Critical Link to Global Configuration ---
    SimCtx *simCtx;  ///< Back-pointer to the master simulation context.

    // --- Grid, Geometry & Parallelization (Per-Level) ---
    DM da, fda, fda2;
    DMDALocalInfo info;
    AO ao;
    PetscInt IM, JM, KM; // No.of grid points.
    PetscReal Max_X, Max_Y, Max_Z, Min_X, Min_Y, Min_Z;
    BoundingBox bbox;
    RankNeighbors neighbors;
    PetscInt GridOrientation, isc, jsc, ksc, _this;
    PetscReal rx, ry, rz;
    PetscInt cgrid;

    // --- Boundary Conditions & Interfacing ---
    BoundaryFaceConfig boundary_faces[6];
    PetscBool inletFaceDefined;
    BCFace    identifiedInletBCFace;
    BCS       Bcs;
    Vec       lFriction_Velocity; 
    PetscInt  bctype[6]; // Legacy BC setup
    PetscReal FluxIntpSum,FluxIntfcSum;

    // --- Primary Flow Fields (Global & Local Views) ---
    Vec Ucont, lUcont, Ucat, lUcat, P, lP, Phi, lPhi, Nvert, lNvert;

    // --- Time-Stepping & Solver Workspace Fields ---
    Vec Ucont_o, lUcont_o, Ucat_o, P_o, Nvert_o, lNvert_o;
    Vec Ucont_rm1, lUcont_rm1, Rhs, dUcont, pUcont;
    Vec CellFieldAtCorner, lCellFieldAtCorner;

    // --- Pressure-Poisson System ---
    Mat A, C; KSP ksp; MatNullSpace nullsp;
    PetscInt *KSKE;
    PetscBool multinullspace;
    Vec B,R;
    Mat MR, MP;
    PetscBool assignedA;
  
  
    // --- Grid Metrics (Global and Local) ---
  Vec Cent, lCent, Csi, Eta, Zet, Aj, lCsi, lEta, lZet, lAj, GridSpace,lGridSpace;
  Vec   Centx,Centy,Centz; // these are already local vectors.
  Vec 	ICsi, IEta, IZet, IAj, lICsi, lIEta, lIZet, lIAj;
  Vec 	JCsi, JEta, JZet, JAj, lJCsi, lJEta, lJZet, lJAj;
  Vec 	KCsi, KEta, KZet, KAj, lKCsi, lKEta, lKZet, lKAj;

  // --- Turbulence Modeling (LES/RANS) ---
  Vec Nu_t, lNu_t, CS, lCs, K_Omega, lK_Omega, K_Omega_o, lK_Omega_o, Distance;

  // --- Statistical Averaging ---
  Vec Ucat_sum, Ucat_cross_sum, Ucat_square_sum, P_sum;

  // --- Immersed Boundary Method (IBM) ---
  IBMNodes *ibm; IBMList *ibmlist;

  // --- Multigrid Hierarchy ---
  PetscInt thislevel, mglevels;
  UserCtx *user_f, *user_c;
  DM *da_f, *da_c;
  Vec *lNvert_c;

  // --- Particle System ---
  DM swarm;
  RankCellInfo *RankCellInfoMap;
  Vec ParticleCount, lParticleCount;
  Vec Psi, lPsi; //scalar dummy to demonstrate scatter.
  
  // --- Post-Processing ---
  DM  post_swarm;
  Vec P_nodal;
  Vec Ucat_nodal;
  Vec Qcrit;
  Vec Psi_nodal;
  
} UserCtx;

#ifdef __cplusplus
}
#endif

#endif // VARIABLES_H
