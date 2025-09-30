/**
 * @file setup.c  //  Setup code for running any simulation 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

 #include "setup.h"

#undef __FUNCT__
#define __FUNCT__ "CreateSimulationContext"

/**
 * @brief Allocates and populates the master SimulationContext object.
 *
 * This function serves as the single, authoritative entry point for all
 * simulation configuration. It merges the setup logic from both the legacy
 * FSI/IBM solver and the modern particle solver into a unified, robust
 * process.
 *
 * The function follows a strict sequence:
 * 1.  **Allocate Context & Set Defaults:** It first allocates the `SimulationContext`
 *     and populates every field with a sane, hardcoded default value. This
 *     ensures the simulation starts from a known, predictable state.
 * 2.  **Configure Logging System:** It configures the custom logging framework. It
 *     parses the `-func_config_file` option to load a list of function names
 *     allowed to produce log output. This configuration (the file path and the
 *     list of function names) is stored within the `SimulationContext` for
 *     later reference and cleanup.
 * 3.  **Parse All Options:** It performs a comprehensive pass of `PetscOptionsGet...`
 *     calls for every possible command-line flag, overriding the default
 *     values set in step 1.
 * 4.  **Log Summary:** After all options are parsed, it uses the now-active
 *     logging system to print a summary of the key simulation parameters.
 *
 * @param[in]  argc      Argument count passed from `main`.
 * @param[in]  argv      Argument vector passed from `main`.
 * @param[out] p_simCtx  On success, this will point to the newly created and
 *                       fully configured `SimulationContext` pointer. The caller
 *                       is responsible for eventually destroying this object by
 *                       calling `FinalizeSimulation()`.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode CreateSimulationContext(int argc, char **argv, SimCtx **p_simCtx)
{
    PetscErrorCode ierr;
    SimCtx *simCtx;
    char control_filename[PETSC_MAX_PATH_LEN] = ""; // Temporary placeholder for control file name.
    PetscBool control_flg; // Temporary placeholder for control file tag existence check flag.

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    // === 1. Allocate the Context Struct and Set ALL Defaults ==================
    ierr = PetscNew(p_simCtx); CHKERRQ(ierr);
    simCtx = *p_simCtx;

    // --- Group 1: Parallelism & MPI ---
    simCtx->rank = 0; simCtx->size = 1;

    // --- Group 2: Simulation Control, Time, and I/O ---
    simCtx->step = 0; simCtx->ti = 0.0; simCtx->StartStep = 0; simCtx->StepsToRun = 10;
    simCtx->tiout = 10; simCtx->StartTime = 0.0; simCtx->dt = 0.001;
    simCtx->OnlySetup = PETSC_FALSE;
    simCtx->logviewer = NULL; simCtx->OutputFreq = simCtx->tiout;
    strcpy(simCtx->eulerianSource,"solve");
    strcpy(simCtx->restart_dir,"results");
    strcpy(simCtx->output_dir,"results");
    strcpy(simCtx->log_dir,"logs");

    // --- Group 3: High-Level Physics & Model Selection Flags ---
    simCtx->immersed = 0; simCtx->movefsi = 0; simCtx->rotatefsi = 0;
    simCtx->sediment = 0; simCtx->rheology = 0; simCtx->invicid = 0;
    simCtx->TwoD = 0; simCtx->thin = 0; simCtx->moveframe = 0;
    simCtx->rotateframe = 0; simCtx->blank = 0;
    simCtx->dgf_x = 0; simCtx->dgf_y = 1; simCtx->dgf_z = 0;
    simCtx->dgf_ax = 1; simCtx->dgf_ay = 0; simCtx->dgf_az = 0;
    simCtx->st = 1.0;

    // --- Group 4: Specific Simulation Case Flags --- (DEPRICATED)
    simCtx->cop=0; simCtx->fish=0; simCtx->fish_c=0; simCtx->fishcyl=0;
    simCtx->eel=0; simCtx->pizza=0; simCtx->turbine=0; simCtx->Pipe=0;
    simCtx->wing=0; simCtx->hydro=0; simCtx->MHV=0; simCtx->LV=0;

    // --- Group 5: Solver & Numerics Parameters ---
    simCtx->implicit = 0; simCtx->implicit_type = 0; simCtx->imp_MAX_IT = 50;
    simCtx->imp_atol = 1e-7; simCtx->imp_rtol = 1e-4; simCtx->imp_stol = 1.e-8;
    simCtx->mglevels = 3; simCtx->mg_MAX_IT = 30; simCtx->mg_idx = 1;
    simCtx->mg_preItr = 1; simCtx->mg_poItr = 1;
    simCtx->poisson = 0; simCtx->poisson_tol = 5.e-9;
    simCtx->STRONG_COUPLING = 0;simCtx->central=0;
    simCtx->ren = 100.0; simCtx->cfl = 0.1; simCtx->vnn = 0.1;
    simCtx->cdisx = 0.0; simCtx->cdisy = 0.0; simCtx->cdisz = 0.0;
    simCtx->FieldInitialization = 0;
    simCtx->InitialConstantContra.x = 0.0;
    simCtx->InitialConstantContra.y = 0.0;
    simCtx->InitialConstantContra.z = 0.0;

    // --- Group 6: Physical & Geometric Parameters ---
    simCtx->NumberOfBodies = 1; simCtx->Flux_in = 1.0; simCtx->angle = 0.0;
    simCtx->max_angle = -54. * 3.1415926 / 180.;
    simCtx->CMx_c=0.0; simCtx->CMy_c=0.0; simCtx->CMz_c=0.0;

    // --- Group 7: Grid, Domain, and Boundary Condition Settings ---
    simCtx->block_number = 1; simCtx->inletprofile = 1;
    simCtx->grid1d = 0; simCtx->Ogrid = 0; simCtx->channelz = 0;
    simCtx->i_periodic = 0; simCtx->j_periodic = 0; simCtx->k_periodic = 0;
    simCtx->blkpbc = 10; simCtx->pseudo_periodic = 0;
    strcpy(simCtx->grid_file, "config/grid.run");
    simCtx->generate_grid = PETSC_FALSE;
    simCtx->da_procs_x = PETSC_DECIDE;
    simCtx->da_procs_y = PETSC_DECIDE;
    simCtx->da_procs_z = PETSC_DECIDE;
    simCtx->grid_rotation_angle  = 0.0;
    simCtx->Croty = 0.0; simCtx->Crotz = 0.0;
    simCtx->num_bcs_files = 1;
    ierr = PetscMalloc1(1, &simCtx->bcs_files); CHKERRQ(ierr);
    ierr = PetscStrallocpy("config/bcs.run", &simCtx->bcs_files[0]); CHKERRQ(ierr);
    simCtx->FluxInSum = 0.0; simCtx->FluxOutSum = 0.0; simCtx->Fluxsum = 0.0;
    simCtx->AreaInSum = 0.0; simCtx->AreaOutSum = 0.0;
    simCtx->U_bc = 0.0; simCtx->ccc = 0;
    simCtx->ratio = 0.0;
    
    
    // --- Group 8: Turbulence Modeling (LES/RANS) ---
    simCtx->averaging = PETSC_FALSE; simCtx->les = 0; simCtx->rans = 0;
    simCtx->wallfunction = 0; simCtx->mixed = 0; simCtx->clark = 0;
    simCtx->dynamic_freq = 1; simCtx->max_cs = 0.5;
    simCtx->testfilter_ik = 0; simCtx->testfilter_1d = 0;
    simCtx->i_homo_filter = 0; simCtx->j_homo_filter = 0; simCtx->k_homo_filter = 0;

    // --- Group 9: Particle / DMSwarm Data & Settings ---
    simCtx->np = 0; simCtx->readFields = PETSC_FALSE;
    simCtx->dm_swarm = NULL; simCtx->bboxlist = NULL;
    simCtx->ParticleInitialization = 0;
    strcpy(simCtx->particleRestartMode,"load");

    // --- Group 10: Immersed Boundary & FSI Data Object Pointers ---
    simCtx->ibm = NULL; simCtx->ibmv = NULL; simCtx->fsi = NULL;
    simCtx->rstart_fsi = PETSC_FALSE; simCtx->duplicate = 0;

    // --- Group 11: Logging and Custom Configuration ---
    strcpy(simCtx->allowedFile, "config/whitelist.run");
    simCtx->useCfg = PETSC_FALSE;
    simCtx->allowedFuncs = NULL;
    simCtx->nAllowed = 0;
    simCtx->LoggingFrequency = 10;
    simCtx->summationRHS = 0.0;
    simCtx->MaxDiv = 0.0;
    simCtx->MaxDivFlatArg = 0; simCtx->MaxDivx = 0; simCtx->MaxDivy = 0; simCtx->MaxDivz = 0;
    strcpy(simCtx->criticalFuncsFile, "config/profile.run");
    simCtx->useCriticalFuncsCfg =  PETSC_FALSE;
    simCtx->criticalFuncs = NULL;
    simCtx->nCriticalFuncs = 0;
    // --- Group 11: Post-Processing Information ---
    strcpy(simCtx->PostprocessingControlFile, "config/post.run");
    ierr = PetscNew(&simCtx->pps); CHKERRQ(ierr);

    // === 2. Get MPI Info and Handle Config File =============================
    // -- Group 1:  Parallelism & MPI Information
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &simCtx->rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &simCtx->size); CHKERRQ(ierr);
    
 // First, check if the -control_file argument was provided by the user/script.
    ierr = PetscOptionsGetString(NULL, NULL, "-control_file", control_filename, sizeof(control_filename), &control_flg); CHKERRQ(ierr);

    // If the flag is NOT present or the filename is empty, abort with a helpful error.
    if (!control_flg || strlen(control_filename) == 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "\n\n*** MANDATORY ARGUMENT MISSING ***\n"
                "The -control_file argument was not provided.\n"
                "This program must be launched with a configuration file.\n"
                "Example: mpiexec -n 4 ./picsolver -control_file /path/to/your/config.control\n"
                "This is typically handled automatically by the 'pic-flow' script.\n");
    }

    // At this point, we have a valid filename. Attempt to load it.
    LOG(GLOBAL, LOG_INFO, "Loading mandatory configuration from: %s\n", control_filename);
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, control_filename, PETSC_FALSE);
    if (ierr == PETSC_ERR_FILE_OPEN) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "The specified control file was not found or could not be opened: %s", control_filename);
    }
    CHKERRQ(ierr);

    // === 3. A Configure Logging System ========================================
    // This logic determines the logging configuration and STORES it in simCtx for
    // later reference and cleanup.
    ierr = PetscOptionsGetString(NULL, NULL, "-whitelist_config_file", simCtx->allowedFile, PETSC_MAX_PATH_LEN, &simCtx->useCfg); CHKERRQ(ierr);
    
    if (simCtx->useCfg) {
        ierr = LoadAllowedFunctionsFromFile(simCtx->allowedFile, &simCtx->allowedFuncs, &simCtx->nAllowed);
        if (ierr) {
            // Use direct PetscPrintf as logging system isn't fully active yet.
            PetscPrintf(PETSC_COMM_SELF, "[%s] WARNING: Failed to load allowed functions from '%s'. Falling back to default list.\n", __func__, simCtx->allowedFile);
            simCtx->useCfg = PETSC_FALSE; // Mark as failed.
            ierr = 0; // Clear the error to allow fallback.
        }
    }
    if (!simCtx->useCfg) {
        // Fallback to default logging functions if no file was used or if loading failed.
        simCtx->nAllowed = 2;
        ierr = PetscMalloc1(simCtx->nAllowed, &simCtx->allowedFuncs); CHKERRQ(ierr);
        ierr = PetscStrallocpy("main", &simCtx->allowedFuncs[0]); CHKERRQ(ierr);
        ierr = PetscStrallocpy("CreateSimulationContext", &simCtx->allowedFuncs[1]); CHKERRQ(ierr);
    }
    
    // Activate the configuration by passing it to the logging module's setup function.
    set_allowed_functions((const char**)simCtx->allowedFuncs, (size_t)simCtx->nAllowed);
    
    // Now that the logger is configured, we can use it.
    LOG_ALLOW(LOCAL, LOG_INFO, "Context created. Initializing on rank %d of %d.\n", simCtx->rank, simCtx->size);
    print_log_level(); // This will now correctly reflect the LOG_LEVEL environment variable.

    // === 3.B Configure Profiling System ========================================
    ierr = PetscOptionsGetString(NULL, NULL, "-profile_config_file", simCtx->criticalFuncsFile, PETSC_MAX_PATH_LEN, &simCtx->useCriticalFuncsCfg); CHKERRQ(ierr);
        if (simCtx->useCriticalFuncsCfg) {
        ierr = LoadAllowedFunctionsFromFile(simCtx->criticalFuncsFile, &simCtx->criticalFuncs, &simCtx->nCriticalFuncs);
        if (ierr) {
            PetscPrintf(PETSC_COMM_SELF, "[%s] WARNING: Failed to load critical profiling functions from '%s'. Falling back to default list.\n", __func__, simCtx->criticalFuncsFile);
            simCtx->useCriticalFuncsCfg = PETSC_FALSE;
            ierr = 0;
        }
    }
    if (!simCtx->useCriticalFuncsCfg) {
        // Fallback to a hardcoded default list if no file or loading failed
        simCtx->nCriticalFuncs = 4;
        ierr = PetscMalloc1(simCtx->nCriticalFuncs, &simCtx->criticalFuncs); CHKERRQ(ierr);
        ierr = PetscStrallocpy("Flow_Solver", &simCtx->criticalFuncs[0]); CHKERRQ(ierr);
        ierr = PetscStrallocpy("AdvanceSimulation", &simCtx->criticalFuncs[1]); CHKERRQ(ierr);
        ierr = PetscStrallocpy("LocateAllParticlesInGrid_TEST", &simCtx->criticalFuncs[2]); CHKERRQ(ierr);
        ierr = PetscStrallocpy("InterpolateAllFieldsToSwarm", &simCtx->criticalFuncs[3]); CHKERRQ(ierr);
    }

    // Initialize the profiling system with the current updated simulation context.
    ierr = ProfilingInitialize(simCtx); CHKERRQ(ierr);

    // === 4. Parse All Command Line Options ==================================
    LOG_ALLOW(GLOBAL, LOG_INFO, "Parsing command-line options...\n");

    // --- Group 2
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 2: Simulation Control,Time and I/O.\n");
    // Read the physical time to start from.
    // The default is already 0.0, so this will only be non-zero if the user provides it.
    ierr = PetscOptionsGetInt(NULL, NULL, "-start_step", &simCtx->StartStep, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-totalsteps", &simCtx->StepsToRun, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-only_setup", &simCtx->OnlySetup, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &simCtx->dt, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-tio", &simCtx->tiout, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-euler_field_source",simCtx->eulerianSource,64,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-output_dir",&simCtx->output_dir,sizeof(simCtx->output_dir),NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-restart_dir",&simCtx->restart_dir,sizeof(simCtx->restart_dir),NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-log_dir",&simCtx->log_dir,sizeof(simCtx->log_dir),NULL);CHKERRQ(ierr);

    simCtx->OutputFreq = simCtx->tiout; // backward compatibility related redundancy.
    if(strcmp(simCtx->eulerianSource,"solve")!= 0 && strcmp(simCtx->eulerianSource,"load") != 0){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Invalid value for -euler_field_source. Must be 'load' or 'solve'. You provided '%s'.",simCtx->eulerianSource);
    }

    //  --- Group 3
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 3: High-Level Physics & Model Selection Flags\n");    
    ierr = PetscOptionsGetInt(NULL, NULL, "-imm", &simCtx->immersed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-fsi", &simCtx->movefsi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rfsi", &simCtx->rotatefsi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-sediment", &simCtx->sediment, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rheology", &simCtx->rheology, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-inv", &simCtx->invicid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-TwoD", &simCtx->TwoD, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-thin", &simCtx->thin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mframe", &simCtx->moveframe, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rframe", &simCtx->rotateframe, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-blk", &simCtx->blank, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_z", &simCtx->dgf_z, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_y", &simCtx->dgf_y, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_x", &simCtx->dgf_x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_az", &simCtx->dgf_az, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_ay", &simCtx->dgf_ay, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_ax", &simCtx->dgf_ax, NULL); CHKERRQ(ierr);

    //  --- Group 4
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 4: Specific Simulation Case Flags \n");
    ierr = PetscOptionsGetInt(NULL, NULL, "-cop", &simCtx->cop, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-fish", &simCtx->fish, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pizza", &simCtx->pizza, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-turbine", &simCtx->turbine, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-fishcyl", &simCtx->fishcyl, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-eel", &simCtx->eel, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-cstart", &simCtx->fish_c, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-wing", &simCtx->wing, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mhv", &simCtx->MHV, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-hydro", &simCtx->hydro, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-lv", &simCtx->LV, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-Pipe", &simCtx->Pipe, NULL); CHKERRQ(ierr);

    //  --- Group 5
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 5: Solver & Numerics Parameters \n");
    ierr = PetscOptionsGetInt(NULL, NULL, "-imp", &simCtx->implicit, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-imp_type", &simCtx->implicit_type, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-imp_MAX_IT", &simCtx->imp_MAX_IT, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-imp_atol", &simCtx->imp_atol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-imp_rtol", &simCtx->imp_rtol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-imp_stol", &simCtx->imp_stol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-central", &simCtx->central, NULL); CHKERRQ(ierr);

    // --- Multigrid Options ---
    ierr = PetscOptionsGetInt(NULL, NULL, "-mg_level", &simCtx->mglevels, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mg_max_it", &simCtx->mg_MAX_IT, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mg_idx", &simCtx->mg_idx, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mg_pre_it", &simCtx->mg_preItr, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mg_post_it", &simCtx->mg_poItr, NULL); CHKERRQ(ierr);

    // --- Other Solver Options ---
    ierr = PetscOptionsGetInt(NULL, NULL, "-poisson", &simCtx->poisson, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-poisson_tol", &simCtx->poisson_tol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-str", &simCtx->STRONG_COUPLING, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-ren", &simCtx->ren, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-cfl", &simCtx->cfl, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-vnn", &simCtx->vnn, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-finit", &simCtx->FieldInitialization, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-ucont_x", &simCtx->InitialConstantContra.x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-ucont_y", &simCtx->InitialConstantContra.y, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-ucont_z", &simCtx->InitialConstantContra.z, NULL); CHKERRQ(ierr);
    // NOTE: cdisx,cdisy,cdisz haven't been parsed, add if necessary.

     //  --- Group 6
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 6: Physical & Geometric Parameters \n");   
    ierr = PetscOptionsGetInt(NULL, NULL, "-no_of_bodies", &simCtx->NumberOfBodies, NULL); CHKERRQ(ierr);
    // NOTE: angle is not parsed in the original code, it set programmatically. We will follow that.
    // NOTE: max_angle is calculated based on other flags (like MHV) in the legacy code.
    // We will defer that logic to a later setup stage and not parse them directly.
    // The Scaling Information is calculated here
    ierr = ParseScalingInformation(simCtx); CHKERRQ(ierr);


     //  --- Group 7
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 7: Grid, Domain, and Boundary Condition Settings \n");
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", &simCtx->block_number, NULL); CHKERRQ(ierr); // This is also a modern option
    ierr = PetscOptionsGetInt(NULL, NULL, "-inlet", &simCtx->inletprofile, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-Ogrid", &simCtx->Ogrid, NULL); CHKERRQ(ierr);
    // NOTE: channelz was not parsed, likely set programmatically. We will omit its parsing call.
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid1d", &simCtx->grid1d, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-grid", &simCtx->generate_grid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL, NULL, "-grid_file", simCtx->grid_file, PETSC_MAX_PATH_LEN, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-da_processors_x", &simCtx->da_procs_x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-da_processors_y", &simCtx->da_procs_y, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-da_processors_z", &simCtx->da_procs_z, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-i_periodic", &simCtx->i_periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-j_periodic", &simCtx->j_periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-k_periodic", &simCtx->k_periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pbc_domain", &simCtx->blkpbc, NULL); CHKERRQ(ierr);
    // NOTE: pseudo_periodic was not parsed. We will omit its parsing call. 
    ierr = PetscOptionsGetReal(NULL, NULL, "-grid_rotation_angle", &simCtx->grid_rotation_angle, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-Croty", &simCtx->Croty, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-Crotz", &simCtx->Crotz, NULL); CHKERRQ(ierr);
    PetscBool bcs_flg;
    char      file_list_str[PETSC_MAX_PATH_LEN * 10]; // Buffer for comma-separated list

    ierr = PetscOptionsGetString(NULL, NULL, "-bcs_files", file_list_str, sizeof(file_list_str), &bcs_flg); CHKERRQ(ierr);
     ierr = PetscOptionsGetReal(NULL, NULL, "-U_bc", &simCtx->U_bc, NULL); CHKERRQ(ierr);

    if (bcs_flg) {
      LOG_ALLOW(GLOBAL, LOG_DEBUG, "Found -bcs_files option, overriding default.\n");
    
      // A. Clean up the default memory we allocated in Phase 1.
      ierr = PetscFree(simCtx->bcs_files[0]); CHKERRQ(ierr);
      ierr = PetscFree(simCtx->bcs_files); CHKERRQ(ierr);
      simCtx->num_bcs_files = 0;
      simCtx->bcs_files = NULL;

      // B. Parse the user-provided comma-separated list.
      char *token;
      char *str_copy;
      ierr = PetscStrallocpy(file_list_str, &str_copy); CHKERRQ(ierr);

      // First pass: count the number of files.
      token = strtok(str_copy, ",");
      while (token) {
        simCtx->num_bcs_files++;
        token = strtok(NULL, ",");
      }
      ierr = PetscFree(str_copy); CHKERRQ(ierr);
    
      // Second pass: allocate memory and store the filenames.
      ierr = PetscMalloc1(simCtx->num_bcs_files, &simCtx->bcs_files); CHKERRQ(ierr);
      ierr = PetscStrallocpy(file_list_str, &str_copy); CHKERRQ(ierr);
      token = strtok(str_copy, ",");
      for (PetscInt i = 0; i < simCtx->num_bcs_files; i++) {
        ierr = PetscStrallocpy(token, &simCtx->bcs_files[i]); CHKERRQ(ierr);
        token = strtok(NULL, ",");
      }
      ierr = PetscFree(str_copy); CHKERRQ(ierr);
    }
    

     //  --- Group 8
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 8: Turbulence Modeling (LES/RANS) \n");
    ierr = PetscOptionsGetInt(NULL, NULL, "-les", &simCtx->les, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rans", &simCtx->rans, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-wallfunction", &simCtx->wallfunction, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mixed", &simCtx->mixed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-clark", &simCtx->clark, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dynamic_freq", &simCtx->dynamic_freq, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-max_cs", &simCtx->max_cs, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-testfilter_ik", &simCtx->testfilter_ik, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-testfilter_1d", &simCtx->testfilter_1d, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-i_homo_filter", &simCtx->i_homo_filter, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-j_homo_filter", &simCtx->j_homo_filter, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-k_homo_filter", &simCtx->k_homo_filter, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-averaging", &simCtx->averaging, NULL); CHKERRQ(ierr);

     //  --- Group 9
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 9:  Particle / DMSwarm Data & Settings \n");
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", &simCtx->np, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", &simCtx->readFields, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &simCtx->ParticleInitialization, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-particle_restart_mode",simCtx->particleRestartMode,sizeof(simCtx->particleRestartMode),NULL); CHKERRQ(ierr);
    // Validation for Particle Restart Mode
    if (strcmp(simCtx->particleRestartMode, "load") != 0 && strcmp(simCtx->particleRestartMode, "init") != 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Invalid value for -particle_restart_mode. Must be 'load' or 'init'. You provided '%s'.", simCtx->particleRestartMode);
    }
    // --- Group 10
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 10: Immersed Boundary & FSI Data Object Pointers \n");
    ierr = PetscOptionsGetBool(NULL, NULL, "-rs_fsi", &simCtx->rstart_fsi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-duplicate", &simCtx->duplicate, NULL); CHKERRQ(ierr);

    // --- Group 11
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 11: Top-Level Managers & Custom Configuration \n");
    ierr = PetscOptionsGetInt(NULL, NULL, "-logfreq", &simCtx->LoggingFrequency, NULL); CHKERRQ(ierr);

    if (simCtx->num_bcs_files != simCtx->block_number) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP, "Number of BC files (%d) does not match number of blocks (%d). Use -bcs_files \"file1.dat,file2.dat,...\".", simCtx->num_bcs_files, simCtx->block_number);
    }
    
    // --- Group 12
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Parsing Group 12: Post-Processing Information.\n");
    // This logic determines the Post Processing configuration and STORES it in simCtx for later reference and cleanup.
    ierr = PetscOptionsGetString(NULL,NULL,"-postprocessing_config_file",simCtx->PostprocessingControlFile,PETSC_MAX_PATH_LEN,NULL); CHKERRQ(ierr);
    ierr = PetscNew(&simCtx->pps); CHKERRQ(ierr);
    ierr = ParsePostProcessingSettings(simCtx);

    // === 5. Dependent Parameter Calculations ================================
    // Some parameters depend on others, so we calculate them here.
    simCtx->StartTime = (PetscReal)simCtx->StartStep*simCtx->dt;
    simCtx->ti = simCtx->StartTime;
    simCtx->step = simCtx->StartStep;

    // === 5. Log Summary and Finalize Setup ==================================
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "-- Console Output Functions [Total : %d] : --\n", simCtx->nAllowed);
    for (PetscInt i = 0; i < simCtx->nAllowed; ++i) {
      LOG_ALLOW(GLOBAL, LOG_DEBUG, "   [%2d] «%s»\n", i, simCtx->allowedFuncs[i]);
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Configuration complete. Key parameters:\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Run mode: %s\n", simCtx->OnlySetup ? "SETUP ONLY" : "Full Simulation");
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Time steps: %d (from %d to %d)\n", simCtx->StepsToRun, simCtx->StartStep, simCtx->StartStep + simCtx->StepsToRun);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Time step size (dt): %g\n", simCtx->dt);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Immersed Boundary: %s\n", simCtx->immersed ? "ENABLED" : "DISABLED");
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Particles: %d\n", simCtx->np);
    if (simCtx->StartStep > 0 && simCtx->np > 0) {
      LOG_ALLOW(GLOBAL, LOG_INFO, "    - Particle Restart Mode: %s\n", simCtx->particleRestartMode);
    }
    
    // --- Initialize PETSc's internal performance logging stage ---
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Finished CreateSimulationContext successfully on rank %d.\n", simCtx->rank);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AllocateContextHeirarchy"
/**
 * @brief Allocates the memory for the UserMG and UserCtx hierarchy.
 *
 * This function performs the foundational memory allocation for the solver's
 * data structures. It reads the number of multigrid levels and grid blocks
 * from the master SimulationContext, then allocates the arrays for MGCtx and
 * UserCtx.
 *
 * It performs three critical tasks:
 * 1.  Allocates the MGCtx array within the UserMG struct.
 * 2.  For each level, allocates the UserCtx array for all blocks.
 * 3.  Sets the `simCtx` back-pointer in every single UserCtx, linking it to the
 *     master configuration.
 *
 * @param simCtx The master SimulationContext, which contains configuration and
 *               will store the resulting allocated hierarchy in its `usermg` member.
 * @return PetscErrorCode
 */
static PetscErrorCode AllocateContextHierarchy(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    MGCtx          *mgctx;
    PetscInt       nblk = simCtx->block_number;
    PetscBool      found;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Allocating context hierarchy for %d levels and %d blocks...\n", simCtx->mglevels, nblk);

    // Store the number of levels in the UserMG struct itself
    usermg->mglevels = simCtx->mglevels;

    // --- 1. Allocate the array of MGCtx structs ---
    ierr = PetscMalloc(usermg->mglevels * sizeof(MGCtx), &usermg->mgctx); CHKERRQ(ierr);
    mgctx = usermg->mgctx;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Allocated MGCtx array of size %d.\n", simCtx->rank, usermg->mglevels);
    
    // --- 2. Parse semi-coarsening options (logic from MG_Initial) ---
    // These flags determine if a dimension is coarsened in the multigrid hierarchy.
    PetscInt *isc, *jsc, *ksc;
    ierr = PetscMalloc3(nblk, &isc, nblk, &jsc, nblk, &ksc); CHKERRQ(ierr);
    // Set defaults to FALSE (full coarsening)
    for (PetscInt i = 0; i < nblk; ++i) {
        isc[i] = 0; jsc[i] = 0; ksc[i] = 0;
    }

// Use a temporary variable for the 'count' argument to the parsing function.
    // This protects the original 'nblk' which is needed for the loop bounds.
    PetscInt n_opts_found = nblk;
    ierr = PetscOptionsGetIntArray(NULL, NULL, "-mg_i_semi", isc, &n_opts_found, &found); CHKERRQ(ierr);

    n_opts_found = nblk; // Reset the temp variable before the next call
    ierr = PetscOptionsGetIntArray(NULL, NULL, "-mg_j_semi", jsc, &n_opts_found, &found); CHKERRQ(ierr);

    n_opts_found = nblk; // Reset the temp variable before the next call
    ierr = PetscOptionsGetIntArray(NULL, NULL, "-mg_k_semi", ksc, &n_opts_found, &found); CHKERRQ(ierr);

    // --- 3. Loop over levels and blocks to allocate UserCtx arrays ---
    for (PetscInt level = 0; level < simCtx->mglevels; level++) {

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Setting up MG Level %d...\n", simCtx->rank, level);        
        // Allocate the array of UserCtx structs for this level
        ierr = PetscMalloc(nblk * sizeof(UserCtx), &mgctx[level].user); CHKERRQ(ierr);
        // It's good practice to zero out the memory to avoid uninitialized values
        ierr = PetscMemzero(mgctx[level].user, nblk * sizeof(UserCtx)); CHKERRQ(ierr);
        mgctx[level].thislevel = level;

        for (PetscInt bi = 0; bi < nblk; bi++) {
            UserCtx *currentUser = &mgctx[level].user[bi];
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   Initializing UserCtx for Level %d, Block %d.\n", simCtx->rank, level, bi);
	    
            // --- CRITICAL STEP: Set the back-pointer to the master context ---
            currentUser->simCtx = simCtx;

            // Initialize other per-context values
            currentUser->thislevel = level;
            currentUser->_this = bi; // Store the block index
            currentUser->mglevels = usermg->mglevels;

            // Assign semi-coarsening flags
            currentUser->isc = isc[bi];
            currentUser->jsc = jsc[bi];
            currentUser->ksc = ksc[bi];

            // Link to finer/coarser contexts for multigrid operations
            if (level > 0) {
                currentUser->user_c = &mgctx[level-1].user[bi];
		LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d:     -> Linked to coarser context (user_c).\n", simCtx->rank);
            }
            if (level < usermg->mglevels - 1) {
                currentUser->user_f = &mgctx[level+1].user[bi];
		LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d:     -> Linked to finer context (user_f).\n", simCtx->rank);
            }
        }
    }

    // Log a summary of the parsed flags on each rank.
    if (get_log_level() >= LOG_DEBUG && nblk > 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Final semi-coarsening configuration view:\n", simCtx->rank);
        for (PetscInt bi = 0; bi < nblk; ++bi) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   Block %d: i-semi=%d, j-semi=%d, k-semi=%d\n", simCtx->rank, bi, isc[bi], jsc[bi], ksc[bi]);
        }
    }

    // Clean up temporary arrays
    ierr = PetscFree3(isc, jsc, ksc); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Context hierarchy allocation complete.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetupSolverParameters"
static PetscErrorCode SetupSolverParameters(SimCtx *simCtx){
  
  PetscFunctionBeginUser;
  PROFILE_FUNCTION_BEGIN;

  LOG_ALLOW(GLOBAL,LOG_INFO, " -- Setting up solver parameters -- .\n");

  UserMG         *usermg = &simCtx->usermg;
  MGCtx          *mgctx = usermg->mgctx;
  PetscInt       nblk = simCtx->block_number;

  for (PetscInt level  = usermg->mglevels-1; level >=0; level--) {
    for (PetscInt bi = 0; bi < nblk; bi++) {
      UserCtx *user = &mgctx[level].user[bi];
      LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Setting up parameters for level %d, block %d\n", simCtx->rank, level, bi);

      user->assignedA = PETSC_FALSE;
      user->multinullspace = PETSC_FALSE;
    }
  }
  PROFILE_FUNCTION_END;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetupGridAndSolvers"
/**
 * @brief The main orchestrator for setting up all grid-related components.
 * (Implementation of the orchestrator function itself)
 */
PetscErrorCode SetupGridAndSolvers(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Starting Grid and Solvers Setup ---\n");

    // Phase 1: Allocate the UserMG and UserCtx hierarchy
    ierr = AllocateContextHierarchy(simCtx); CHKERRQ(ierr);

    // [ The next phases will be added here ]
    ierr = DefineAllGridDimensions(simCtx); CHKERRQ(ierr);
    ierr = InitializeAllGridDMs(simCtx); CHKERRQ(ierr);
    ierr = AssignAllGridCoordinates(simCtx);
    ierr = CreateAndInitializeAllVectors(simCtx); CHKERRQ(ierr);
    ierr = SetupSolverParameters(simCtx); CHKERRQ(ierr);
    ierr = CalculateAllGridMetrics(simCtx); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Grid and Solvers Setup Complete ---\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CreateAndInitializeAllVectors"
/**
 * @brief Creates and initializes all PETSc Vec objects for all fields.
 *
 * This function iterates through every UserCtx in the multigrid and multi-block
 * hierarchy. For each context, it creates the comprehensive set of global and
 * local PETSc Vecs required by the flow solver (e.g., Ucont, P, Nvert, metrics,
 * turbulence fields, etc.). Each vector is initialized to zero.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode
 */
PetscErrorCode CreateAndInitializeAllVectors(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    MGCtx          *mgctx = usermg->mgctx;
    PetscInt       nblk = simCtx->block_number;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Creating and initializing all simulation vectors...\n");
        
    for (PetscInt level  = usermg->mglevels-1; level >=0; level--) {
        for (PetscInt bi = 0; bi < nblk; bi++) {
            UserCtx *user = &mgctx[level].user[bi];

            if(!user->da || !user->fda) {
              SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "DMs not properly initialized in UserCtx before vector creation.");
            }

            LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Creating vectors for level %d, block %d\n", simCtx->rank, level, bi);

	    // --- Group A: Primary Flow Fields (Global and Local) ---
            // These are the core solution variables.
            ierr = DMCreateGlobalVector(user->fda, &user->Ucont); CHKERRQ(ierr); ierr = VecSet(user->Ucont, 0.0); CHKERRQ(ierr);
            ierr = DMCreateGlobalVector(user->fda, &user->Ucat);  CHKERRQ(ierr); ierr = VecSet(user->Ucat, 0.0); CHKERRQ(ierr);
            ierr = DMCreateGlobalVector(user->da,  &user->P);     CHKERRQ(ierr); ierr = VecSet(user->P, 0.0); CHKERRQ(ierr);
            ierr = DMCreateGlobalVector(user->da,  &user->Nvert); CHKERRQ(ierr); ierr = VecSet(user->Nvert, 0.0); CHKERRQ(ierr);

            ierr = DMCreateLocalVector(user->fda, &user->lUcont); CHKERRQ(ierr); ierr = VecSet(user->lUcont, 0.0); CHKERRQ(ierr);
            ierr = DMCreateLocalVector(user->fda, &user->lUcat);  CHKERRQ(ierr); ierr = VecSet(user->lUcat, 0.0); CHKERRQ(ierr);
            ierr = DMCreateLocalVector(user->da,  &user->lP);     CHKERRQ(ierr); ierr = VecSet(user->lP, 0.0); CHKERRQ(ierr);
            ierr = DMCreateLocalVector(user->da,  &user->lNvert); CHKERRQ(ierr); ierr = VecSet(user->lNvert, 0.0); CHKERRQ(ierr);

            // --- Group B: Time-Stepping & Workspace Fields (Finest Level Only) ---
            if (level == usermg->mglevels - 1) {
                ierr = VecDuplicate(user->Ucont, &user->Ucont_o);   CHKERRQ(ierr); ierr = VecSet(user->Ucont_o, 0.0); CHKERRQ(ierr);
                ierr = VecDuplicate(user->Ucont, &user->Ucont_rm1); CHKERRQ(ierr); ierr = VecSet(user->Ucont_rm1, 0.0); CHKERRQ(ierr);
                ierr = VecDuplicate(user->Ucat,  &user->Ucat_o);    CHKERRQ(ierr); ierr = VecSet(user->Ucat_o, 0.0); CHKERRQ(ierr);
                ierr = VecDuplicate(user->P,     &user->P_o);       CHKERRQ(ierr); ierr = VecSet(user->P_o, 0.0); CHKERRQ(ierr);
	              ierr = VecDuplicate(user->lUcont, &user->lUcont_o);   CHKERRQ(ierr); ierr = VecSet(user->lUcont_o, 0.0); CHKERRQ(ierr);
                ierr = VecDuplicate(user->lUcont, &user->lUcont_rm1); CHKERRQ(ierr); ierr = VecSet(user->lUcont_rm1, 0.0); CHKERRQ(ierr);
                ierr = DMCreateLocalVector(user->da, &user->lNvert_o); CHKERRQ(ierr); ierr = VecSet(user->lNvert_o, 0.0); CHKERRQ(ierr);
		            ierr = VecDuplicate(user->Nvert, &user->Nvert_o); CHKERRQ(ierr); ierr = VecSet(user->Nvert_o, 0.0); CHKERRQ(ierr);
		
            }

	    // --- Group C: Grid Metrics (Cell-Centered) ---
            ierr = DMCreateGlobalVector(user->fda, &user->Csi); CHKERRQ(ierr); ierr = VecSet(user->Csi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->Eta);         CHKERRQ(ierr); ierr = VecSet(user->Eta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->Zet);         CHKERRQ(ierr); ierr = VecSet(user->Zet, 0.0); CHKERRQ(ierr);
            ierr = DMCreateGlobalVector(user->da,  &user->Aj);  CHKERRQ(ierr); ierr = VecSet(user->Aj, 0.0); CHKERRQ(ierr);
	          ierr = VecDuplicate(user->Aj, &user->Phi);       CHKERRQ(ierr); ierr = VecSet(user->Phi, 0.0); CHKERRQ(ierr);

            ierr = DMCreateLocalVector(user->fda, &user->lCsi); CHKERRQ(ierr); ierr = VecSet(user->lCsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lEta);       CHKERRQ(ierr); ierr = VecSet(user->lEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lZet);       CHKERRQ(ierr); ierr = VecSet(user->lZet, 0.0); CHKERRQ(ierr);
            ierr = DMCreateLocalVector(user->da,  &user->lAj);  CHKERRQ(ierr); ierr = VecSet(user->lAj, 0.0); CHKERRQ(ierr);
	     ierr = VecDuplicate(user->lAj, &user->lPhi);       CHKERRQ(ierr); ierr = VecSet(user->lPhi, 0.0); CHKERRQ(ierr);
	    // --- Group D: Grid Metrics (Face-Centered) ---
            // Vector metrics are duplicated from Csi (DOF=3, fda-based)
            ierr = VecDuplicate(user->Csi, &user->ICsi); CHKERRQ(ierr); ierr = VecSet(user->ICsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->IEta); CHKERRQ(ierr); ierr = VecSet(user->IEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->IZet); CHKERRQ(ierr); ierr = VecSet(user->IZet, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->JCsi); CHKERRQ(ierr); ierr = VecSet(user->JCsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->JEta); CHKERRQ(ierr); ierr = VecSet(user->JEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->JZet); CHKERRQ(ierr); ierr = VecSet(user->JZet, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->KCsi); CHKERRQ(ierr); ierr = VecSet(user->KCsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->KEta); CHKERRQ(ierr); ierr = VecSet(user->KEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Csi, &user->KZet); CHKERRQ(ierr); ierr = VecSet(user->KZet, 0.0); CHKERRQ(ierr);
            // Scalar metrics are duplicated from Aj (DOF=1, da-based)
            ierr = VecDuplicate(user->Aj, &user->IAj); CHKERRQ(ierr); ierr = VecSet(user->IAj, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Aj, &user->JAj); CHKERRQ(ierr); ierr = VecSet(user->JAj, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Aj, &user->KAj); CHKERRQ(ierr); ierr = VecSet(user->KAj, 0.0); CHKERRQ(ierr);

            ierr = VecDuplicate(user->lCsi, &user->lICsi); CHKERRQ(ierr); ierr = VecSet(user->lICsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lIEta); CHKERRQ(ierr); ierr = VecSet(user->lIEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lIZet); CHKERRQ(ierr); ierr = VecSet(user->lIZet, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lJCsi); CHKERRQ(ierr); ierr = VecSet(user->lJCsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lJEta); CHKERRQ(ierr); ierr = VecSet(user->lJEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lJZet); CHKERRQ(ierr); ierr = VecSet(user->lJZet, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lKCsi); CHKERRQ(ierr); ierr = VecSet(user->lKCsi, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lKEta); CHKERRQ(ierr); ierr = VecSet(user->lKEta, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCsi, &user->lKZet); CHKERRQ(ierr); ierr = VecSet(user->lKZet, 0.0); CHKERRQ(ierr);

            ierr = VecDuplicate(user->lAj, &user->lIAj); CHKERRQ(ierr); ierr = VecSet(user->lIAj, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lAj, &user->lJAj); CHKERRQ(ierr); ierr = VecSet(user->lJAj, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lAj, &user->lKAj); CHKERRQ(ierr); ierr = VecSet(user->lKAj, 0.0); CHKERRQ(ierr);
	    
	    // --- Group E: Cell/Face Center Coordinates and Grid Spacing ---
            ierr = DMCreateGlobalVector(user->fda, &user->Cent); CHKERRQ(ierr); ierr = VecSet(user->Cent, 0.0); CHKERRQ(ierr);
            ierr = DMCreateLocalVector(user->fda, &user->lCent); CHKERRQ(ierr); ierr = VecSet(user->lCent, 0.0); CHKERRQ(ierr);
            
            ierr = VecDuplicate(user->Cent, &user->GridSpace); CHKERRQ(ierr); ierr = VecSet(user->GridSpace, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->lCent, &user->lGridSpace); CHKERRQ(ierr); ierr = VecSet(user->lGridSpace, 0.0); CHKERRQ(ierr);

            // Face-center coordinate vectors are GLOBAL to hold calculated values before scattering
            ierr = VecDuplicate(user->Cent, &user->Centx); CHKERRQ(ierr); ierr = VecSet(user->Centx, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Cent, &user->Centy); CHKERRQ(ierr); ierr = VecSet(user->Centy, 0.0); CHKERRQ(ierr);
            ierr = VecDuplicate(user->Cent, &user->Centz); CHKERRQ(ierr); ierr = VecSet(user->Centz, 0.0); CHKERRQ(ierr);

	    if(level == usermg->mglevels -1){
	    // --- Group F: Turbulence Models (Finest Level Only) ---
            if (simCtx->les || simCtx->rans) {
                ierr = DMCreateGlobalVector(user->da, &user->Nu_t); CHKERRQ(ierr); ierr = VecSet(user->Nu_t, 0.0); CHKERRQ(ierr);
                ierr = DMCreateLocalVector(user->da, &user->lNu_t); CHKERRQ(ierr); ierr = VecSet(user->lNu_t, 0.0); CHKERRQ(ierr);
		
	    // Add K_Omega, CS, etc. here as needed

            // Note: Add any other vectors from the legacy MG_Initial here as needed.
            // For example: Rhs, Forcing, turbulence Vecs (K_Omega, Nu_t)...
		
	        }
	    // --- Group G: Particle Methods 	
	    if(simCtx->np>0){
	      ierr = DMCreateGlobalVector(user->da,&user->ParticleCount); CHKERRQ(ierr); ierr = VecSet(user->ParticleCount,0.0); CHKERRQ(ierr);

	       ierr = DMCreateGlobalVector(user->da,&user->Psi); CHKERRQ(ierr); ierr = VecSet(user->Psi,0.0); CHKERRQ(ierr);

	      LOG_ALLOW(GLOBAL,LOG_DEBUG,"ParticleCount & Scalar(Psi) created for %d particles.\n",simCtx->np);
	      }
	    }
	    // --- Group H: Boundary Condition vectors needed by the legacy FormBCS ---
	    ierr = DMCreateGlobalVector(user->fda, &user->Bcs.Ubcs); CHKERRQ(ierr);
	    ierr = VecSet(user->Bcs.Ubcs, 0.0); CHKERRQ(ierr);
	    ierr = DMCreateGlobalVector(user->fda, &user->Bcs.Uch); CHKERRQ(ierr);
	    ierr = VecSet(user->Bcs.Uch, 0.0); CHKERRQ(ierr);
	    
      if(level == usermg->mglevels -1){
        if(simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR){
                LOG_ALLOW(LOCAL, LOG_DEBUG, "Post-processor mode detected. Allocating derived field vectors.\n");

                ierr = VecDuplicate(user->P, &user->P_nodal); CHKERRQ(ierr);
                ierr = VecSet(user->P_nodal, 0.0); CHKERRQ(ierr);

                ierr = VecDuplicate(user->Ucat, &user->Ucat_nodal); CHKERRQ(ierr);
                ierr = VecSet(user->Ucat_nodal, 0.0); CHKERRQ(ierr);
                
                ierr = VecDuplicate(user->P, &user->Qcrit); CHKERRQ(ierr);
                ierr = VecSet(user->Qcrit, 0.0); CHKERRQ(ierr);          
        }else{
                user->P_nodal = NULL;
                user->Ucat_nodal = NULL;
                user->Qcrit = NULL;
        }
	  }
	
  }
}

    LOG_ALLOW(GLOBAL, LOG_INFO, "All simulation vectors created and initialized.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateLocalGhosts"
/**
 * @brief Updates the local vector (including ghost points) from its corresponding global vector.
 *
 * This function identifies the correct global vector, local vector, and DM based on the
 * provided fieldName and performs the standard PETSc DMGlobalToLocalBegin/End sequence.
 * Includes optional debugging output (max norms before/after).
 *
 * @param user       The UserCtx structure containing the vectors and DMs.
 * @param fieldName  The name of the field to update ("Ucat", "Ucont", "P", "Nvert", etc.).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * @note This function assumes the global vector associated with fieldName has already
 *       been populated with the desired data (including any boundary conditions).
 */
PetscErrorCode UpdateLocalGhosts(UserCtx* user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Vec            globalVec = NULL;
    Vec            localVec = NULL;
    DM             dm = NULL; // The DM associated with this field pair

    PetscFunctionBeginUser; // Use User version for application code
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Starting ghost update for field '%s'.\n", rank, fieldName);

    // --- 1. Identify the correct Vectors and DM ---
    if (strcmp(fieldName, "Ucat") == 0) {
        globalVec = user->Ucat;
        localVec  = user->lUcat;
        dm        = user->fda;
    } else if (strcmp(fieldName, "Ucont") == 0) {
        globalVec = user->Ucont;
        localVec  = user->lUcont;
        dm        = user->fda;
    } else if (strcmp(fieldName, "P") == 0) {
        globalVec = user->P;
        localVec  = user->lP;
        dm        = user->da;
    } else if (strcmp(fieldName, "Csi") == 0) {
        globalVec = user->Csi;
        localVec  = user->lCsi;
        dm        = user->fda;
    } else if (strcmp(fieldName, "Eta") == 0) {
        globalVec = user->Eta;
        localVec  = user->lEta;
        dm        = user->fda;
    }  else if (strcmp(fieldName, "Zet") == 0) {
        globalVec = user->Zet;
        localVec  = user->lZet;
        dm        = user->fda;
    }else if (strcmp(fieldName, "Nvert") == 0) {
        globalVec = user->Nvert;
        localVec  = user->lNvert;
        dm        = user->da;
     // Add other fields as needed
    } else if (strcmp(fieldName, "Aj") == 0) {
        globalVec = user->Aj;
        localVec  = user->lAj;
        dm        = user->da;
    } else if (strcmp(fieldName, "Cent") == 0) {
        globalVec = user->Cent;
        localVec  = user->lCent;
        dm        = user->fda;
    }else if (strcmp(fieldName, "GridSpace") == 0) {
        globalVec = user->GridSpace;
        localVec  = user->lGridSpace;
        dm        = user->fda;
    }else if (strcmp(fieldName,"ICsi") == 0){
      globalVec = user->ICsi;
      localVec  = user->lICsi;
      dm        = user->fda;
    }else if (strcmp(fieldName,"IEta") == 0){
      globalVec = user->IEta;
      localVec  = user->lIEta;
      dm        = user->fda;
    }else if (strcmp(fieldName,"IZet") == 0){
      globalVec = user->IZet;
      localVec  = user->lIZet;
      dm        = user->fda;
    }else if (strcmp(fieldName,"JCsi") == 0){
      globalVec = user->JCsi;
      localVec  = user->lJCsi;
      dm        = user->fda;
    }else if (strcmp(fieldName,"JEta") == 0){
      globalVec = user->JEta;
      localVec  = user->lJEta;
      dm        = user->fda;
    }else if (strcmp(fieldName,"JZet") == 0){
      globalVec = user->JZet;
      localVec  = user->lJZet;
      dm        = user->fda;
    }else if (strcmp(fieldName,"KCsi") == 0){
      globalVec = user->KCsi;
      localVec  = user->lKCsi;
      dm        = user->fda;
    }else if (strcmp(fieldName,"KEta") == 0){
      globalVec = user->KEta;
      localVec  = user->lKEta;
      dm        = user->fda;
    }else if (strcmp(fieldName,"KZet") == 0){
      globalVec = user->KZet;
      localVec  = user->lKZet;
      dm        = user->fda;
    }else if (strcmp(fieldName,"IAj") == 0){
      globalVec = user->IAj;
      localVec  = user->lIAj;
      dm        = user->da;
    }else if (strcmp(fieldName,"JAj") == 0){
      globalVec = user->JAj;
      localVec  = user->lJAj;
      dm        = user->da;
    }else if (strcmp(fieldName,"KAj") == 0){
      globalVec = user->KAj;
      localVec  = user->lKAj;
      dm        = user->da;
    }else if (strcmp(fieldName,"Phi") == 0){
      globalVec = user->Phi;
      localVec  = user->lPhi;
      dm        = user->da;
    }else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Field '%s' not recognized for ghost update.", fieldName);
    }

    // --- 2. Check if components were found ---
    if (!globalVec) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Global vector for field '%s' is NULL.", fieldName);
    }
    if (!localVec) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Local vector for field '%s' is NULL.", fieldName);
    }
    if (!dm) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "DM for field '%s' is NULL.", fieldName);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Identified components for '%s': DM=%p, GlobalVec=%p, LocalVec=%p.\n",
              rank, fieldName, (void*)dm, (void*)globalVec, (void*)localVec);

    // --- 3. Optional Debugging: Norm Before Update ---
    // Use your logging convention check
    // if (get_log_level() >= LOG_LEVEL_DEBUG && is_function_allowed("UpdateLocalGhosts")) { // Example check
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){
        PetscReal norm_global_before;
        ierr = VecNorm(globalVec, NORM_INFINITY, &norm_global_before); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO,"Max norm '%s' (Global) BEFORE Ghost Update: %g\n", fieldName, norm_global_before);
        // Optional: Norm of local vector before update (might contain old ghost values)
        // PetscReal norm_local_before;
        // ierr = VecNorm(localVec, NORM_INFINITY, &norm_local_before); CHKERRQ(ierr);
        // LOG_ALLOW(GLOBAL, LOG_DEBUG,"Max norm '%s' (Local) BEFORE Ghost Update: %g\n", fieldName, norm_local_before);
    }

    // --- 4. Perform the Global-to-Local Transfer (Ghost Update) ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Calling DMGlobalToLocalBegin/End for '%s'.\n", rank, fieldName);
    ierr = DMGlobalToLocalBegin(dm, globalVec, INSERT_VALUES, localVec); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, globalVec, INSERT_VALUES, localVec); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Completed DMGlobalToLocalBegin/End for '%s'.\n", rank, fieldName);

    // --- 5. Optional Debugging: Norm After Update ---
    // Use your logging convention check
    // if (get_log_level() >= LOG_LEVEL_DEBUG && is_function_allowed("UpdateLocalGhosts")) { // Example check
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){ // Using your specific check
        PetscReal norm_local_after;
        ierr = VecNorm(localVec, NORM_INFINITY, &norm_local_after); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO,"Max norm '%s' (Local) AFTER Ghost Update: %g\n", fieldName, norm_local_after);

        // --- 6. Optional Debugging: Specific Point Checks (Example for Ucat on Rank 0/1) ---
        //    (Keep this conditional if it's only for specific debug scenarios)
        if (strcmp(fieldName, "Ucat") == 0) { // Only do detailed checks for Ucat for now
           PetscMPIInt rank_test;
           MPI_Comm_rank(PETSC_COMM_WORLD, &rank_test);

           // Get Local Info needed for indexing checks
           DMDALocalInfo info_check;
           ierr = DMDAGetLocalInfo(dm, &info_check); CHKERRQ(ierr); // Use the correct dm

           // Buffer for array pointer
           Cmpnts ***lUcat_arr_test = NULL;
           PetscErrorCode ierr_test = 0;

           LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Testing '%s' access immediately after ghost update...\n", rank_test, fieldName);
           ierr_test = DMDAVecGetArrayDOFRead(dm, localVec, &lUcat_arr_test); // Use correct dm and localVec

           if (ierr_test) {
               LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR %d getting '%s' array after ghost update!\n", rank_test, ierr_test, fieldName);
           } else if (!lUcat_arr_test) {
                LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR NULL pointer getting '%s' array after ghost update!\n", rank_test, fieldName);
           }
           else {
               // Check owned interior point (e.g., first interior point)
               PetscInt k_int = info_check.zs + (info_check.zm > 1 ? 1 : 0); // Global k index (at least zs+1 if possible)
               PetscInt j_int = info_check.ys + (info_check.ym > 1 ? 1 : 0); // Global j index
               PetscInt i_int = info_check.xs + (info_check.xm > 1 ? 1 : 0); // Global i index
                // Ensure indices are within global bounds if domain is very small
               //if (k_int >= info_check.mz-1) k_int = info_check.mz-2; if (k_int < 1) k_int = 1;
               //if (j_int >= info_check.my-1) j_int = info_check.my-2; if (j_int < 1) j_int = 1;
	       // if (i_int >= info_check.mx-1) i_int = info_check.mx-2; if (i_int < 1) i_int = 1;
	       // clamp k_int to [1 .. mz-2] 
	       if (k_int >= info_check.mz - 1) {
		 k_int = info_check.mz - 2;
	       }
	       if (k_int < 1) {
		 k_int = 1;
	       }

	       // clamp j_int to [1 .. my-2] 
	       if (j_int >= info_check.my - 1) {
		 j_int = info_check.my - 2;
	       }
	       if (j_int < 1) {
		 j_int = 1;
	       }

	       // clamp i_int to [1 .. mx-2]
	       if (i_int >= info_check.mx - 1) {
		 i_int = info_check.mx - 2;
	       }
	       if (i_int < 1) {
		 i_int = 1;
	       }

               // Only attempt read if indices are actually owned (relevant for multi-rank)
               if (k_int >= info_check.zs && k_int < info_check.zs + info_check.zm &&
                   j_int >= info_check.ys && j_int < info_check.ys + info_check.ym &&
                   i_int >= info_check.xs && i_int < info_check.xs + info_check.xm)
               {
                  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Attempting test read OWNED INTERIOR [%d][%d][%d] (Global)\n", rank_test, k_int, j_int, i_int);
                  Cmpnts test_val_owned_interior = lUcat_arr_test[k_int][j_int][i_int];
                  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: SUCCESS reading owned interior: x=%g\n", rank_test, test_val_owned_interior.x);
               } else {
                  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Skipping interior test read for non-owned index [%d][%d][%d].\n", rank_test, k_int, j_int, i_int);
               }


               // Check owned boundary point (e.g., first owned point)
               PetscInt k_bnd = info_check.zs; // Global k index
               PetscInt j_bnd = info_check.ys; // Global j index
               PetscInt i_bnd = info_check.xs; // Global i index
               LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Attempting test read OWNED BOUNDARY [%d][%d][%d] (Global)\n", rank_test, k_bnd, j_bnd, i_bnd);
               Cmpnts test_val_owned_boundary = lUcat_arr_test[k_bnd][j_bnd][i_bnd];
               LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: SUCCESS reading owned boundary: x=%g\n", rank_test, test_val_owned_boundary.x);


               // Check ghost point (e.g., one layer below in k, if applicable)
               if (info_check.zs > 0) { // Only if there's a rank below
                   PetscInt k_ghost = info_check.zs - 1;
                   PetscInt j_ghost = info_check.ys; // Use start of owned y, simple example
                   PetscInt i_ghost = info_check.xs; // Use start of owned x, simple example
                   LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Attempting test read GHOST [%d][%d][%d] (Global)\n", rank_test, k_ghost, j_ghost, i_ghost);
                   Cmpnts test_val_ghost = lUcat_arr_test[k_ghost][j_ghost][i_ghost];
                   LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: SUCCESS reading ghost: x=%g\n", rank_test, test_val_ghost.x);
               } else {
                   LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Skipping ghost test read (zs=0).\n", rank_test);
               }

               // Restore the array
               ierr_test = DMDAVecRestoreArrayDOFRead(dm, localVec, &lUcat_arr_test);
               if(ierr_test){ LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR %d restoring '%s' array after test read!\n", rank_test, ierr_test, fieldName); }
               LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Finished testing '%s' access.\n", rank_test, fieldName);
           }
        } // end if Ucat
    } // end debug logging check

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Completed ghost update for field '%s'.\n", rank, fieldName);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetupBoundaryConditions"
/**
 * @brief (Orchestrator) Sets up all boundary conditions for the simulation.
 */
PetscErrorCode SetupBoundaryConditions(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL,LOG_INFO, "--- Setting up Boundary Conditions ---\n");
    
    // --- Parse and Adapt for each block on the finest level ---
    LOG_ALLOW(GLOBAL,LOG_INFO,"Parsing BC configuration file and adapting to legacy system for finest grid.\n");
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(GLOBAL,LOG_DEBUG, "  -> Processing Block %d:\n", bi);

	// --- Generate the filename for the current block ---
	const char *current_bc_filename = simCtx->bcs_files[bi];
	LOG_ALLOW(GLOBAL,LOG_DEBUG,"  -> Processing Block %d using config file '%s'\n", bi, current_bc_filename);
       // This will populate user_finest[bi].boundary_faces
        ierr = ParseAllBoundaryConditions(&user_finest[bi],current_bc_filename); CHKERRQ(ierr);

        // Call the adapter to translate into the legacy format
        ierr = TranslateModernBCsToLegacy(&user_finest[bi]); CHKERRQ(ierr);

        // Call the function to calculate the center of the inlet face, which may be used to calculate Boundary values.
        ierr = CalculateInletCenter(&user_finest[bi]); CHKERRQ(ierr); 
    }

    LOG_ALLOW(GLOBAL,LOG_INFO, "--- Boundary Conditions setup complete ---\n");   


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Allocates a 3D array of PetscReal values using PetscCalloc.
 *
 * This function dynamically allocates memory for a 3D array of PetscReal values
 * with dimensions nz (layers) x ny (rows) x nx (columns). It uses PetscCalloc1
 * to ensure the memory is zero-initialized.
 *
 * The allocation is done in three steps:
 *  1. Allocate an array of nz pointers (one for each layer).
 *  2. Allocate a contiguous block for nz*ny row pointers and assign each layer’s row pointers.
 *  3. Allocate a contiguous block for all nz*ny*nx PetscReal values.
 *
 * This setup allows the array to be accessed as array[k][j][i], and the memory
 * for the data is contiguous, which improves cache efficiency.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers (z-direction).
 * @param[in]  ny    Number of rows (y-direction).
 * @param[in]  nx    Number of columns (x-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Allocate3DArrayScalar(PetscReal ****array, PetscInt nz, PetscInt ny, PetscInt nx)
{
  PetscErrorCode ierr;
  PetscReal      ***data;
  PetscReal      *dataContiguous;
  PetscInt       k, j;

  PetscFunctionBegin;
  /* Step 1: Allocate memory for an array of nz layer pointers (zero-initialized) */
  ierr = PetscCalloc1(nz, &data); CHKERRQ(ierr);

  /* Step 2: Allocate memory for all row pointers (nz * ny pointers) */
  ierr = PetscCalloc1(nz * ny, &data[0]); CHKERRQ(ierr);
  for (k = 1; k < nz; k++) {
    data[k] = data[0] + k * ny;
  }

  /* Step 3: Allocate one contiguous block for all data elements (nz*ny*nx) */
  ierr = PetscCalloc1(nz * ny * nx, &dataContiguous); CHKERRQ(ierr);

  /* Build the 3D pointer structure: each row pointer gets the correct segment of data */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      data[k][j] = dataContiguous + (k * ny + j) * nx;
      /* Memory is already zeroed by PetscCalloc1, so no manual initialization is needed */
    }
  }
  *array = data;
  PetscFunctionReturn(0);
}

/**
 * @brief Deallocates a 3D array of PetscReal values allocated by Allocate3DArrayScalar.
 *
 * This function frees the memory allocated for a 3D array of PetscReal values.
 * It assumes the memory was allocated using Allocate3DArrayScalar, which allocated
 * three separate memory blocks: one for the contiguous data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers (z-direction).
 * @param[in] ny    Number of rows (y-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Deallocate3DArrayScalar(PetscReal ***array, PetscInt nz, PetscInt ny)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
   if (!array || !array[0] || !array[0][0] ) { // Added more robust check
      LOG_ALLOW(GLOBAL, LOG_WARNING, "Deallocate3DArrayScalar called with potentially unallocated or NULL array.\n");
       if (array) {
           if (array[0]) { // Check if row pointers might exist
               // Cannot safely access array[0][0] if array[0] might be invalid/freed
               // Standard deallocation below assumes valid pointers.
                ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free row pointers if they exist
           }
            ierr = PetscFree(array); CHKERRQ(ierr); // Free layer pointers if they exist
       }
       PetscFunctionReturn(0);
   }

  // --- Standard Deallocation (assuming valid allocation) ---

  /* 1. Free the contiguous block of PetscReal values.
     The starting address was stored in array[0][0]. */
  ierr = PetscFree(array[0][0]); CHKERRQ(ierr); // Free the ACTUAL DATA

  /* 2. Free the contiguous block of row pointers.
     The starting address was stored in array[0]. */
  ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free the ROW POINTERS

  /* 3. Free the layer pointer array.
     The starting address is 'array' itself. */
  ierr = PetscFree(array); CHKERRQ(ierr); // Free the LAYER POINTERS

  PetscFunctionReturn(0);
}

/**
 * @brief Allocates a 3D array of Cmpnts structures using PetscCalloc.
 *
 * This function dynamically allocates memory for a 3D array of Cmpnts (vector) structures
 * with dimensions nz (layers) x ny (rows) x nx (columns). It uses PetscCalloc1 to ensure
 * that all allocated memory is zero-initialized.
 *
 * The allocation procedure is similar to Allocate3DArrayScalar:
 *  1. Allocate an array of nz pointers (one for each layer).
 *  2. Allocate a contiguous block for nz*ny row pointers.
 *  3. Allocate one contiguous block for nz*ny*nx Cmpnts structures.
 *
 * After allocation, the array can be accessed as array[k][j][i] and each element
 * (a Cmpnts structure) will have its x, y, and z fields initialized to 0.0.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers in the z-direction.
 * @param[in]  ny    Number of rows in the y-direction.
 * @param[in]  nx    Number of columns in the x-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Allocate3DArrayVector(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx)
{
  PetscErrorCode ierr;
  Cmpnts         ***data;
  Cmpnts         *dataContiguous;
  PetscInt       k, j;
  PetscMPIInt    rank;

  PetscFunctionBegin;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /* Step 1: Allocate memory for nz layer pointers (zeroed) */
  ierr = PetscCalloc1(nz, &data); CHKERRQ(ierr);

  LOG_ALLOW(LOCAL,LOG_DEBUG," [Rank %d] memory allocated for outermost layer (%d k-layer pointers).\n",rank,nz);
  
  /* Step 2: Allocate memory for all row pointers (nz * ny pointers) */
  ierr = PetscCalloc1(nz * ny, &data[0]); CHKERRQ(ierr);
  for (k = 1; k < nz; k++) {
    data[k] = data[0] + k * ny;
  }

  LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] memory allocated for %dx%d row pointers.\n",rank,nz,ny);
  
  /* Step 3: Allocate one contiguous block for nz*ny*nx Cmpnts structures (zeroed) */
  ierr = PetscCalloc1(nz * ny * nx, &dataContiguous); CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG,"[Rank %d] memory allocated for contigous block of %dx%dx%d Cmpnts structures).\n",rank,nz,ny,nx); 

  /* Build the 3D pointer structure for vector data */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      data[k][j] = dataContiguous + (k * ny + j) * nx;
      /* The PetscCalloc1 call has already initialized each Cmpnts to zero. */
    }
  }

  LOG_ALLOW(GLOBAL,LOG_DEBUG,"[Rank %d] 3D pointer structure for vector data created. \n",rank);
  
  *array = data;
  PetscFunctionReturn(0);
}

/**
 * @brief Deallocates a 3D array of Cmpnts structures allocated by Allocate3DArrayVector.
 *
 * This function frees the memory allocated for a 3D array of Cmpnts structures.
 * It assumes the memory was allocated using Allocate3DArrayVector, which created three
 * separate memory blocks: one for the contiguous vector data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers in the z-direction.
 * @param[in] ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
 PetscErrorCode Deallocate3DArrayVector(Cmpnts ***array, PetscInt nz, PetscInt ny)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // If array is NULL or hasn't been allocated properly, just return.
  if (!array || !array[0] || !array[0][0] ) {
      LOG_ALLOW(GLOBAL, LOG_WARNING, "Deallocate3DArrayVector called with potentially unallocated or NULL array.\n");
      // Attempt to free what might exist, but be cautious
      if (array) {
          if (array[0]) { // Check if row pointers were allocated
             // We don't have a direct pointer to the contiguous data block
             // saved separately in this allocation scheme. The allocation relies
             // on array[0][0] pointing to it. If array[0] was freed first,
             // accessing array[0][0] is unsafe.
             // The allocation scheme where the contiguous data block is not
             // stored separately makes safe deallocation tricky if freeing
             // happens out of order or if parts are NULL.

             // A SAFER ALLOCATION/DEALLOCATION would store the data pointer separately.
             // Given the current allocation scheme, the order MUST be:
             // 1. Free the data block (pointed to by array[0][0])
             // 2. Free the row pointer block (pointed to by array[0])
             // 3. Free the layer pointer block (pointed to by array)

             // Let's assume the allocation was successful and pointers are valid.
             // Get pointer to the contiguous data block *before* freeing row pointers
             Cmpnts *dataContiguous = array[0][0];
             ierr = PetscFree(dataContiguous); CHKERRQ(ierr); // Free data block

             // Now free the row pointers block
             ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free row pointers

          }
          // Finally, free the array of layer pointers
          ierr = PetscFree(array); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0); // Return gracefully if input was NULL initially
  }


  // --- Standard Deallocation (assuming valid allocation) ---

  /* 1. Free the contiguous block of Cmpnts structures.
     The starting address was stored in array[0][0] by Allocate3DArrayVector. */
  ierr = PetscFree(array[0][0]); CHKERRQ(ierr); // Free the ACTUAL DATA

  /* 2. Free the contiguous block of row pointers.
     The starting address was stored in array[0]. */
  ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free the ROW POINTERS

  /* 3. Free the layer pointer array.
     The starting address is 'array' itself. */
  ierr = PetscFree(array); CHKERRQ(ierr); // Free the LAYER POINTERS

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetOwnedCellRange"
/**
 * @brief Gets the global starting index of cells owned by this rank and the number of such cells.
 *
 * A cell's global index is considered the same as its origin node's global index.
 * This function assumes a node-centered DMDA where `info_nodes` provides all necessary
 * information:
 *  - `info_nodes->xs, ys, zs`: Global starting index of the first node owned by this rank (excluding ghosts).
 *  - `info_nodes->xm, ym, zm`: Number of nodes owned by this rank in each dimension (excluding ghosts).
 *  - `info_nodes->mx, my, mz`: Total number of global nodes in each dimension for the entire domain.
 *
 * A cell `C_k` (0-indexed) is defined by its origin node `N_k` and extends to node `N_{k+1}`.
 * Thus, the last node in the global domain cannot be an origin for a cell. The last possible
 * cell origin node index is `GlobalNodesInDim - 2`.
 *
 * @param[in] info_nodes Pointer to the DMDALocalInfo struct for the current rank.
 *                       This struct contains local ownership information (xs, xm, etc.)
 *                       and global domain dimensions (mx, my, mz for nodes).
 * @param[in] dim        The dimension for which to get the cell range (0 for i/x, 1 for j/y, 2 for k/z).
 * @param[out] xs_cell_global_out Pointer to store the global index of the first cell whose origin node
 *                                is owned by this rank. If the rank owns no valid cell origins in this
 *                                dimension, this will be the rank's starting node index, but
 *                                `xm_cell_local_out` will be 0.
 * @param[out] xm_cell_local_out  Pointer to store the number of cells for which this rank owns the
 *                                origin node AND that origin node is a valid cell origin within the
 *                                global domain.
 *
 * @return PetscErrorCode 0 on success, or an error code on failure.
 *
 * @note Example: If GlobalNodesInDim = 11 (nodes N0 to N10), there are 10 cells (C0 to C9).
 *       The last cell, C9, has its origin at node N9. So, N9 (index 9) is the last valid
 *       cell origin (GlobalNodesInDim - 2 = 11 - 2 = 9).
 *       If a rank owns nodes N8, N9, N10 (xs=8, xm=3):
 *         - First potential origin on rank = N8.
 *         - Last potential origin on rank (node that is not the last owned node) = N9.
 *         - Actual last origin this rank can form = min(N9, GlobalMaxOrigin=N9) = N9.
 *         - Number of cells = (N9 - N8 + 1) = 2 cells (C8, C9).
 *       If a rank owns only node N10 (xs=10, xm=1):
 *         - First potential origin on rank = N10.
 *         - Actual last origin rank can form = min(N9, GlobalMaxOrigin=N9) (since N10-1=N9).
 *         - first_potential_origin_on_rank (N10) > actual_last_origin_this_rank_can_form (N9) => 0 cells.
 */
PetscErrorCode GetOwnedCellRange(const DMDALocalInfo *info_nodes,
                                 PetscInt dim,
                                 PetscInt *xs_cell_global_out,
                                 PetscInt *xm_cell_local_out)
{
    PetscErrorCode ierr = 0; // Standard PETSc error code, not explicitly set here but good practice.
    PetscInt xs_node_global_rank;   // Global index of the first node owned by this rank in the specified dimension.
    PetscInt num_nodes_owned_rank;  // Number of nodes owned by this rank in this dimension (local count, excluding ghosts).
    PetscInt GlobalNodesInDim_from_info; // Total number of global nodes in this dimension, from DMDALocalInfo.

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- 1. Input Validation ---
    if (!info_nodes || !xs_cell_global_out || !xm_cell_local_out) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null pointer passed to GetOwnedCellRange.");
    }

    // --- 2. Extract Node Ownership and Global Dimension Information from DMDALocalInfo ---
    if (dim == 0) { // I-direction
        xs_node_global_rank = info_nodes->xs;       // Starting owned node index (global)
        num_nodes_owned_rank  = info_nodes->xm;     // Number of nodes owned by this rank (local count)
        GlobalNodesInDim_from_info = info_nodes->mx; // Total global nodes in this dimension
    } else if (dim == 1) { // J-direction
        xs_node_global_rank = info_nodes->ys;
        num_nodes_owned_rank  = info_nodes->ym;
        GlobalNodesInDim_from_info = info_nodes->my;
    } else if (dim == 2) { // K-direction
        xs_node_global_rank = info_nodes->zs;
        num_nodes_owned_rank  = info_nodes->zm;
        GlobalNodesInDim_from_info = info_nodes->mz;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d in GetOwnedCellRange. Must be 0, 1, or 2.", dim);
    }

    // --- 3. Handle Edge Cases for Global Domain Size ---
    // If the global domain has 0 or 1 node plane, no cells can be formed.
    if (GlobalNodesInDim_from_info <= 1) {
        *xs_cell_global_out = xs_node_global_rank; // Still report the rank's starting node
        *xm_cell_local_out = 0;                    // But 0 cells
        PetscFunctionReturn(0);
    }
    // Negative global dimension is an error (should be caught by DMDA setup, but defensive)
    if (GlobalNodesInDim_from_info < 0 ) {
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "GlobalNodesInDim %d from DMDALocalInfo must be non-negative for dimension %d.", GlobalNodesInDim_from_info, dim);
    }

    // --- 4. Determine Cell Ownership Based on Node Ownership ---
    // The first cell this rank *could* define has its origin at the first node this rank owns.
    *xs_cell_global_out = xs_node_global_rank;

    // If the rank owns no nodes in this dimension, it can't form any cell origins.
    if (num_nodes_owned_rank == 0) {
        *xm_cell_local_out = 0;
    } else {
        // Calculate the global index of the last possible node that can serve as a cell origin.
        // If GlobalNodesInDim = N (nodes 0 to N-1), cells are C_0 to C_{N-2}.
        // The origin of cell C_{N-2} is node N_{N-2}.
        // So, the last valid cell origin node index is (GlobalNodesInDim - 2).
        PetscInt last_possible_origin_global_idx = GlobalNodesInDim_from_info - 2;

        // Determine the range of nodes owned by this rank that could *potentially* be cell origins.
        // The first node owned by the rank is a potential origin.
        PetscInt first_potential_origin_on_rank = xs_node_global_rank;

        // A node can be an origin if there's at least one node after it to form the cell.
        // So, the last node owned by the rank that could *potentially* be an origin is
        // the second-to-last node it owns: (xs_node_global_rank + num_nodes_owned_rank - 1) - 1
        // which simplifies to: xs_node_global_rank + num_nodes_owned_rank - 2.
        PetscInt last_potential_origin_on_rank = xs_node_global_rank + num_nodes_owned_rank - 2;

        // The actual last origin this rank can provide is capped by the global domain limit.
        PetscInt actual_last_origin_this_rank_can_form = PetscMin(last_potential_origin_on_rank, last_possible_origin_global_idx);

        // If the first potential origin this rank owns is already beyond the actual last origin it can form,
        // then this rank forms no valid cell origins. This happens if:
        //  - num_nodes_owned_rank is 1 (so last_potential_origin_on_rank = first_potential_origin_on_rank - 1).
        //  - The rank only owns nodes at the very end of the global domain (e.g., only the last global node).
        if (first_potential_origin_on_rank > actual_last_origin_this_rank_can_form) {
            *xm_cell_local_out = 0;
        } else {
            // The number of cells is the count of valid origins this rank owns.
            // (Count = Last Index - First Index + 1)
            *xm_cell_local_out = actual_last_origin_this_rank_can_form - first_potential_origin_on_rank + 1;
        }
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "GetGhostedCellRange"

/**
 * @brief Gets the global cell range for a rank, including boundary cells.
 *
 * This function first calls GetOwnedCellRange to get the conservative range of
 * fully-contained cells. It then extends this range by applying the
 * "Lower-Rank-Owns-Boundary" principle. A rank claims ownership of the
 * boundary cells it shares with neighbors in the positive (+x, +y, +z)
 * directions.
 *
 * This results in a final cell range that is gap-free and suitable for building
 * the definitive particle ownership map.
 *
 * @param[in]  info_nodes       Pointer to the DMDALocalInfo struct.
 * @param[in]  neighbors        Pointer to the RankNeighbors struct containing neighbor info.
 * @param[in]  dim              The dimension (0 for i/x, 1 for j/y, 2 for k/z).
 * @param[out] xs_cell_global_out Pointer to store the final starting cell index.
 * @param[out] xm_cell_local_out  Pointer to store the final number of cells.
 *
 * @return PetscErrorCode 0 on success, or an error code on failure.
 */
PetscErrorCode GetGhostedCellRange(const DMDALocalInfo *info_nodes,
                                   const RankNeighbors *neighbors,
                                   PetscInt dim,
                                   PetscInt *xs_cell_global_out,
                                   PetscInt *xm_cell_local_out)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- 1. Get the base, conservative range from the original function ---
    ierr = GetOwnedCellRange(info_nodes, dim, xs_cell_global_out, xm_cell_local_out); CHKERRQ(ierr);

    // --- 2. Apply the "Lower-Rank-Owns-Boundary" correction ---
    // A rank owns the boundary if it has a neighbor in the positive direction.
    // We check if the neighbor's rank is valid (not MPI_PROC_NULL, which is < 0).
    if (dim == 0 && neighbors->rank_xp > -1) {
        (*xm_cell_local_out)++;
    } else if (dim == 1 && neighbors->rank_yp > -1) {
        (*xm_cell_local_out)++;
    } else if (dim == 2 && neighbors->rank_zp > -1) {
        (*xm_cell_local_out)++;
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeAndStoreNeighborRanks"
/**
 * @brief Computes and stores the Cartesian neighbor ranks for the DMDA decomposition.
 *
 * This function retrieves the neighbor information from the primary DMDA (user->da)
 * and stores the face neighbors (xm, xp, ym, yp, zm, zp) in the user->neighbors structure.
 * It assumes a standard PETSc ordering for the neighbors array returned by DMDAGetNeighbors.
 * If DMDAGetNeighbors returns a negative rank that is not MPI_PROC_NULL (which can happen
 * in some PETSc/MPI configurations for non-periodic boundaries if not fully standard),
 * this function will sanitize it to MPI_PROC_NULL to prevent issues.
 *
 * @param[in,out] user Pointer to the UserCtx structure where neighbor info will be stored.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeAndStoreNeighborRanks(UserCtx *user)
{
    PetscErrorCode    ierr;
    PetscMPIInt       rank;
    PetscMPIInt       size; // MPI communicator size
    const PetscMPIInt *neighbor_ranks_ptr; // Pointer to raw neighbor data from PETSc

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr); // Get MPI size for validation

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Computing DMDA neighbor ranks.\n", rank);

    if (!user || !user->da) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx or user->da is NULL in ComputeAndStoreNeighborRanks.");
    }

    // Get the neighbor information from the DMDA
    // neighbor_ranks_ptr will point to an internal PETSc array of 27 ranks.
    ierr = DMDAGetNeighbors(user->da, &neighbor_ranks_ptr); CHKERRQ(ierr);

    // Log the raw values from DMDAGetNeighbors for boundary-relevant directions for debugging
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d]Raw DMDAGetNeighbors: xm_raw=%d, xp_raw=%d, ym_raw=%d, yp_raw=%d, zm_raw=%d, zp_raw=%d. MPI_PROC_NULL is %d.\n",
                   rank,
                   neighbor_ranks_ptr[12], neighbor_ranks_ptr[14],
                   neighbor_ranks_ptr[10], neighbor_ranks_ptr[16],
                   neighbor_ranks_ptr[4],  neighbor_ranks_ptr[22],
                   (int)MPI_PROC_NULL);

    // PETSc standard indices for 3D face neighbors from the 27-point stencil:
    // Index = k_offset*9 + j_offset*3 + i_offset (where offsets -1,0,1 map to 0,1,2)
    // Center: (i_off=1, j_off=1, k_off=1) => 1*9 + 1*3 + 1 = 13
    // X-min:  (i_off=0, j_off=1, k_off=1) => 1*9 + 1*3 + 0 = 12
    // X-plus: (i_off=2, j_off=1, k_off=1) => 1*9 + 1*3 + 2 = 14
    // Y-min:  (i_off=1, j_off=0, k_off=1) => 1*9 + 0*3 + 1 = 10
    // Y-plus: (i_off=1, j_off=2, k_off=1) => 1*9 + 2*3 + 1 = 16
    // Z-min:  (i_off=1, j_off=1, k_off=0) => 0*9 + 1*3 + 1 = 4
    // Z-plus: (i_off=1, j_off=1, k_off=2) => 2*9 + 1*3 + 1 = 22

    if (neighbor_ranks_ptr[13] != rank) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Rank %d: DMDAGetNeighbors center index (13) is %d, expected current rank %d. Neighbor indexing might be non-standard or DMDA small.\n",
                  rank, neighbor_ranks_ptr[13], rank);
        // This warning is important. If the center isn't the current rank, the offsets are likely wrong.
        // However, PETSc should ensure this unless the DM is too small for a 3x3x3 stencil.
    }

    // Assign and sanitize each neighbor rank
    PetscMPIInt temp_neighbor;

    temp_neighbor = neighbor_ranks_ptr[12]; // xm
    if (temp_neighbor < 0  || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[Rank %d] Correcting invalid xm neighbor %d to MPI_PROC_NULL (%d).\n", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_xm = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_xm = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[14]; // xp
    if (temp_neighbor < 0 || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[Rank %d] Correcting invalid xp neighbor %d to MPI_PROC_NULL (%d).\n", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_xp = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_xp = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[10]; // ym
    if (temp_neighbor < 0 || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[Rank %d] Correcting invalid ym neighbor %d to MPI_PROC_NULL (%d).\n", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_ym = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_ym = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[16]; // yp
    if (temp_neighbor < 0 || temp_neighbor >= size) {
        // The log for index 16 was "zm" in your output, should be yp
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[Rank %d] Correcting invalid yp neighbor (raw index 16) %d to MPI_PROC_NULL (%d).\n", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_yp = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_yp = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[4]; // zm
    if (temp_neighbor < 0 || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[Rank %d] Correcting invalid zm neighbor %d to MPI_PROC_NULL (%d).\n", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_zm = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_zm = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[22]; // zp
    if (temp_neighbor < 0 || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[Rank %d] Correcting invalid zp neighbor %d to MPI_PROC_NULL (%d).\n", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_zp = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_zp = temp_neighbor;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d] Stored user->neighbors: xm=%d, xp=%d, ym=%d, yp=%d, zm=%d, zp=%d\n", rank,
              user->neighbors.rank_xm, user->neighbors.rank_xp,
              user->neighbors.rank_ym, user->neighbors.rank_yp,
              user->neighbors.rank_zm, user->neighbors.rank_zp);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); // Ensure logs are flushed

    // Note: neighbor_ranks_ptr memory is managed by PETSc, do not free it.
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetDMDAProcLayout"
/**
 * @brief Sets the processor layout for a given DMDA based on PETSc options.
 *
 * Reads the desired number of processors in x, y, and z directions using
 * PETSc options (e.g., -dm_processors_x, -dm_processors_y, -dm_processors_z).
 * If an option is not provided for a direction, PETSC_DECIDE is used for that direction.
 * Applies the layout using DMDASetNumProcs.
 *
 * Also stores the retrieved/decided values in user->procs_x/y/z if user context is provided.
 *
 * @param dm   The DMDA object to configure the layout for.
 * @param user Pointer to the UserCtx structure (optional, used to store layout values).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetDMDAProcLayout(DM dm, UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    size, rank;
    PetscInt       px = PETSC_DECIDE, py = PETSC_DECIDE, pz = PETSC_DECIDE;
    PetscBool      px_set = PETSC_FALSE, py_set = PETSC_FALSE, pz_set = PETSC_FALSE;
    SimCtx         *simCtx = user->simCtx;

    // Set no.of processors in direction 1
    if(simCtx->da_procs_x) {
      px_set = PETSC_TRUE;
      px = simCtx->da_procs_x;
    }
    // Set no.of processors in direction 2
    if(simCtx->da_procs_y) {
      py_set = PETSC_TRUE;
      py = simCtx->da_procs_y;
    }
    // Set no.of processors in direction 1
    if(simCtx->da_procs_z) {
      pz_set = PETSC_TRUE;
      pz = simCtx->da_procs_z;
    }

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm), &size); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm), &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Configuring DMDA processor layout for %d total processes.\n", rank, size);

    // --- Validate User Input (Optional but Recommended) ---
    // Check if specified processor counts multiply to the total MPI size
    if (px_set && py_set && pz_set) {
        if (px * py * pz != size) {
             SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_INCOMP,
                     "Specified processor layout %d x %d x %d = %d does not match MPI size %d",
                     px, py, pz, px * py * pz, size);
        }
         LOG_ALLOW(GLOBAL, LOG_INFO, "Using specified processor layout: %d x %d x %d\n", px, py, pz);
    } else if (px_set || py_set || pz_set) {
         // If only some are set, PETSC_DECIDE will be used for others
         LOG_ALLOW(GLOBAL, LOG_INFO, "Using partially specified processor layout: %d x %d x %d (PETSC_DECIDE for unspecified)\n", px, py, pz);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Using fully automatic processor layout (PETSC_DECIDE x PETSC_DECIDE x PETSC_DECIDE)\n");
    }
     // Additional checks: Ensure px, py, pz are positive if set
     if ((px_set && px <= 0) || (py_set && py <= 0) || (pz_set && pz <= 0)) {
         SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_OUTOFRANGE, "Specified processor counts must be positive.");
     }


    // --- Apply the layout to the DMDA ---
    ierr = DMDASetNumProcs(dm, px, py, pz); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank %d: DMDASetNumProcs called with px=%d, py=%d, pz=%d.\n", rank, px, py, pz);

    // --- Store the values in UserCtx (Optional) ---
    // Note: If PETSC_DECIDE was used, PETSc calculates the actual values during DMSetUp.
    // We store the *requested* values here. To get the *actual* values used,
    // you would need to call DMDAGetInfo after DMSetUp.
    /*
    if (user) {
        user->procs_x = px;
        user->procs_y = py;
        user->procs_z = pz;
    }
    */
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetupDomainRankInfo"
/**
 * @brief Sets up the full rank communication infrastructure for all blocks.
 *
 * This function orchestrates the setup of the parallel communication map. It uses the
 * existing helper functions to sequentially gather and broadcast the bounding box
 * information for each computational block. All ranks participate in each step,
 * assembling the final, unified multi-block list in parallel.
 *
 * @param[in,out] simCtx The master simulation context. After this call,
 *                       simCtx->bboxlist and the relevant user->neighbors fields
 *                       will be populated on all ranks.
 */
PetscErrorCode SetupDomainRankInfo(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscInt       nblk = simCtx->block_number;
    PetscInt       size = simCtx->size;
    BoundingBox    *final_bboxlist = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting full rank communication setup for %d block(s).\n", nblk);

    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    // --- Step 1: Compute neighbor ranks (unchanged) ---
    for (int bi = 0; bi < nblk; bi++) {
        ierr = ComputeAndStoreNeighborRanks(&user_finest[bi]); CHKERRQ(ierr);
    }
    LOG_ALLOW(GLOBAL, LOG_INFO, "Neighbor ranks computed and stored for all blocks.\n");

    // --- Step 2: Allocate the final, unified list on ALL ranks ---
    // Every rank will build this list in parallel.
    ierr = PetscMalloc1(size * nblk, &final_bboxlist); CHKERRQ(ierr);

    // --- Step 3: Loop through each block, gather then broadcast its bbox list ---
    for (int bi = 0; bi < nblk; bi++) {
        // This is a temporary pointer for the current block's list.
        BoundingBox *block_bboxlist = NULL;

        LOG_ALLOW(GLOBAL, LOG_INFO, "Processing bounding boxes for block %d...\n", bi);

        // A) GATHER: On rank 0, block_bboxlist is allocated and filled. On others, it's NULL.
        ierr = GatherAllBoundingBoxes(&user_finest[bi], &block_bboxlist); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "  -> Gather complete for block %d.", bi);

        // B) BROADCAST: On non-root ranks, block_bboxlist is allocated. Then, the data
        //    from rank 0 is broadcast to all ranks. After this call, ALL ranks have
        //    an identical, complete copy of the bounding boxes for the current block.
        ierr = BroadcastAllBoundingBoxes(&user_finest[bi], &block_bboxlist); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "  -> Broadcast complete for block %d.\n", bi);

        // C) ASSEMBLE: Every rank now copies the data for this block into the
        //    correct segment of its final, unified list.
        for (int r = 0; r < size; r++) {
            // The layout is [r0b0, r1b0, ..., r(size-1)b0, r0b1, r1b1, ...]
            final_bboxlist[bi * size + r] = block_bboxlist[r];
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "  -> Assembly into final list complete for block %d.\n", bi);

        // D) CLEANUP: Free the temporary list for this block on ALL ranks before the next iteration.
        //    Your helper functions use malloc, so we must use free.
        free(block_bboxlist);
    }

    // --- Step 4: Assign the final pointer and run the last setup step ---
    simCtx->bboxlist = final_bboxlist;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Final unified bboxlist created on all ranks and stored in SimCtx.\n");

    ierr = SetupDomainCellDecompositionMap(&user_finest[0]); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Domain Cell Composition set and broadcasted.\n");

    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankInfo: Completed successfully.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Contra2Cart"
/**
 * @brief Reconstructs Cartesian velocity (Ucat) at cell centers from contravariant
 *        velocity (Ucont) defined on cell faces.
 *
 * This function performs the transformation from a contravariant velocity representation
 * (which is natural on a curvilinear grid) to a Cartesian (x,y,z) representation.
 * For each interior computational cell owned by the rank, it performs the following:
 *
 * 1.  It averages the contravariant velocity components (U¹, U², U³) from the
 *     surrounding faces to get an estimate of the contravariant velocity at the cell center.
 * 2.  It averages the metric vectors (Csi, Eta, Zet) from the surrounding faces
 *     to get an estimate of the metric tensor at the cell center. This tensor forms
 *     the transformation matrix.
 * 3.  It solves the linear system `[MetricTensor] * [ucat] = [ucont]` for the
 *     Cartesian velocity vector `ucat = (u,v,w)` using Cramer's rule.
 * 4.  The computed Cartesian velocity is stored in the global `user->Ucat` vector.
 *
 * The function operates on local, ghosted versions of the input vectors (`user->lUcont`,
 * `user->lCsi`, etc.) to ensure stencils are valid across processor boundaries.
 *
 * @param[in,out] user      Pointer to the UserCtx structure. The function reads from
 *                          `user->lUcont`, `user->lCsi`, `user->lEta`, `user->lZet`, `user->lNvert`
 *                          and writes to the global `user->Ucat` vector.
 *
 * @return PetscErrorCode 0 on success.
 *
 * @note
 *  - This function should be called AFTER `user->lUcont` and all local metric vectors
 *    (`user->lCsi`, etc.) have been populated with up-to-date ghost values via `UpdateLocalGhosts`.
 *  - It only computes `Ucat` for interior cells (not on physical boundaries) and for
 *    cells not marked as solid/blanked by `user->lNvert`.
 *  - The caller is responsible for subsequently applying boundary conditions to `user->Ucat`
 *    and calling `UpdateLocalGhosts(user, "Ucat")` to populate `user->lUcat`.
 */
PetscErrorCode Contra2Cart(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Cmpnts       ***lcsi_arr, ***leta_arr, ***lzet_arr; // Local metric arrays
    Cmpnts       ***lucont_arr;                       // Local contravariant velocity array
    Cmpnts       ***gucat_arr;                        // Global Cartesian velocity array
    PetscReal    ***lnvert_arr;                       // Local Nvert array
    PetscReal    ***laj_arr;                          // Local Jacobian Determinant inverse array

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Starting Contravariant-to-Cartesian velocity transformation.\n");

    // --- 1. Get DMDA Info and Check for Valid Inputs ---
    // All inputs (lUcont, lCsi, etc.) and outputs (Ucat) are on DMs from the UserCtx.
    // We get local info from fda, which governs the layout of most arrays here.
    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);
    if (!user->lUcont || !user->lCsi || !user->lEta || !user->lZet || !user->lNvert || !user->Ucat) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Contra2Cart requires lUcont, lCsi/Eta/Zet, lNvert, and Ucat to be non-NULL.");
    }


    // --- 2. Get Read-Only Array Access to Local Input Vectors (with ghosts) ---
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcont, &lucont_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi,   &lcsi_arr);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta,   &leta_arr);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet,   &lzet_arr);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lNvert, &lnvert_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lAj, &laj_arr); CHKERRQ(ierr);

    // --- 3. Get Write-Only Array Access to the Global Output Vector ---
    // We compute for local owned cells and write into the global vector.
    // PETSc handles mapping the global indices to the correct local memory locations.
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &gucat_arr); CHKERRQ(ierr);


    // --- 4. Define Loop Bounds for INTERIOR Cells ---
    // We use adjusted bounds to avoid calculating Ucat on the physical domain boundaries,
    // as these are typically set explicitly by boundary condition functions.
    // The stencils use indices like i-1, j-1, k-1, so we must start loops at least at index 1.
    PetscInt i_start = (info.xs == 0) ? info.xs + 1 : info.xs;
    PetscInt i_end   = (info.xs + info.xm == info.mx) ? info.xs + info.xm - 1 : info.xs + info.xm;

    PetscInt j_start = (info.ys == 0) ? info.ys + 1 : info.ys;
    PetscInt j_end   = (info.ys + info.ym == info.my) ? info.ys + info.ym - 1 : info.ys + info.ym;

    PetscInt k_start = (info.zs == 0) ? info.zs + 1 : info.zs;
    PetscInt k_end   = (info.zs + info.zm == info.mz) ? info.zs + info.zm - 1 : info.zs + info.zm;

    // --- 5. Main Computation Loop ---
    // Loops over the GLOBAL indices of interior cells owned by this rank.
    for (PetscInt k_cell = k_start; k_cell < k_end; ++k_cell) {
        for (PetscInt j_cell = j_start; j_cell < j_end; ++j_cell) {
            for (PetscInt i_cell = i_start; i_cell < i_end; ++i_cell) {

                // Check if the cell is a fluid cell (not solid/blanked)
	      //    if (lnvert_arr[k_cell][j_cell][i_cell] > 0.1) continue; // Skip solid/blanked cells

                // Transformation matrix [mat] is the metric tensor at the cell center,
                // estimated by averaging metrics from adjacent faces.
                PetscReal mat[3][3];

		//	PetscReal aj_center = laj_arr[k_cell+1][j_cell+1][i_cell+1];
		
                mat[0][0] = 0.5 * (lcsi_arr[k_cell][j_cell][i_cell-1].x + lcsi_arr[k_cell][j_cell][i_cell].x); //* aj_center;
                mat[0][1] = 0.5 * (lcsi_arr[k_cell][j_cell][i_cell-1].y + lcsi_arr[k_cell][j_cell][i_cell].y); //* aj_center;
                mat[0][2] = 0.5 * (lcsi_arr[k_cell][j_cell][i_cell-1].z + lcsi_arr[k_cell][j_cell][i_cell].z); //* aj_center;

                mat[1][0] = 0.5 * (leta_arr[k_cell][j_cell-1][i_cell].x + leta_arr[k_cell][j_cell][i_cell].x); //* aj_center;
                mat[1][1] = 0.5 * (leta_arr[k_cell][j_cell-1][i_cell].y + leta_arr[k_cell][j_cell][i_cell].y); //* aj_center;
                mat[1][2] = 0.5 * (leta_arr[k_cell][j_cell-1][i_cell].z + leta_arr[k_cell][j_cell][i_cell].z); //* aj_center;

                mat[2][0] = 0.5 * (lzet_arr[k_cell-1][j_cell][i_cell].x + lzet_arr[k_cell][j_cell][i_cell].x); //* aj_center;
                mat[2][1] = 0.5 * (lzet_arr[k_cell-1][j_cell][i_cell].y + lzet_arr[k_cell][j_cell][i_cell].y); //* aj_center;
                mat[2][2] = 0.5 * (lzet_arr[k_cell-1][j_cell][i_cell].z + lzet_arr[k_cell][j_cell][i_cell].z); //* aj_center;
	
                // Contravariant velocity vector `q` at the cell center,
                // estimated by averaging face-based contravariant velocities.
                PetscReal q[3];
                q[0] = 0.5 * (lucont_arr[k_cell][j_cell][i_cell-1].x + lucont_arr[k_cell][j_cell][i_cell].x); // U¹ at cell center
                q[1] = 0.5 * (lucont_arr[k_cell][j_cell-1][i_cell].y + lucont_arr[k_cell][j_cell][i_cell].y); // U² at cell center
                q[2] = 0.5 * (lucont_arr[k_cell-1][j_cell][i_cell].z + lucont_arr[k_cell][j_cell][i_cell].z); // U³ at cell center

                // Solve the 3x3 system `mat * ucat = q` using Cramer's rule.
                PetscReal det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                                mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                                mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

		  if (PetscAbsReal(det) < 1.0e-18) {
		      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FLOP_COUNT, "Transformation matrix determinant is near zero at cell (%d,%d,%d) \n", i_cell, j_cell, k_cell);
		                }

                PetscReal det_inv = 1.0 / det;

                PetscReal det0 = q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                                 q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
                                 q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

                PetscReal det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                                  q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
                                  q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

                PetscReal det2 = q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
                                 q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
                                 q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

                // Store computed Cartesian velocity in the GLOBAL Ucat array at the
                // array index corresponding to the cell's origin node.
                gucat_arr[k_cell][j_cell][i_cell].x = det0 * det_inv;
                gucat_arr[k_cell][j_cell][i_cell].y = det1 * det_inv;
                gucat_arr[k_cell][j_cell][i_cell].z = det2 * det_inv;
            }
        }
    }

    // --- 6. Restore Array Access ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcont, &lucont_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi,   &lcsi_arr);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta,   &leta_arr);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet,   &lzet_arr);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lNvert, &lnvert_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lAj, &laj_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &gucat_arr); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed Contravariant-to-Cartesian velocity transformation. \n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetupDomainCellDecompositionMap"
/**
 * @brief Creates and distributes a map of the domain's cell decomposition to all ranks.
 * @ingroup DomainInfo
 *
 * This function is a critical part of the simulation setup. It determines the global
 * cell ownership for each MPI rank and makes this information available to all
 * other ranks. This "decomposition map" is essential for the robust "Walk and Handoff"
 * particle migration strategy, allowing any rank to quickly identify the owner of a
 * target cell.
 *
 * The process involves:
 * 1. Each rank gets its own node ownership information from the DMDA.
 * 2. It converts this node information into cell ownership ranges using the
 *    `GetOwnedCellRange` helper function.
 * 3. It participates in an `MPI_Allgather` collective operation to build a complete
 *    array (`user->RankCellInfoMap`) containing the ownership information for every rank.
 *
 * This function should be called once during initialization after the primary DMDA
 * (user->da) has been set up.
 *
 * @param[in,out] user Pointer to the UserCtx structure. The function will allocate and
 *                     populate `user->RankCellInfoMap` and set `user->num_ranks`.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *         Errors can occur if input pointers are NULL or if MPI communication fails.
 */
PetscErrorCode SetupDomainCellDecompositionMap(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  local_node_info;
    RankCellInfo   my_cell_info;
    PetscMPIInt    rank, size;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- 1. Input Validation and MPI Info ---
    if (!user) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx pointer is NULL in SetupDomainCellDecompositionMap.");
    }
    if (!user->da) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "user->da is not initialized in SetupDomainCellDecompositionMap.");
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Setting up domain cell decomposition map for %d ranks.\n", size);

    // --- 2. Determine Local Cell Ownership ---
    // Get the local node ownership information from the primary DMDA.
    ierr = DMDAGetLocalInfo(user->da, &local_node_info); CHKERRQ(ierr);

    // Use the robust helper function to convert node ownership to cell ownership.
    // A cell's index is defined by its origin node.
    
    ierr = GetGhostedCellRange(&local_node_info, &user->neighbors, 0, &my_cell_info.xs_cell, &my_cell_info.xm_cell); CHKERRQ(ierr);
    ierr = GetGhostedCellRange(&local_node_info, &user->neighbors, 1, &my_cell_info.ys_cell, &my_cell_info.ym_cell); CHKERRQ(ierr);
    ierr = GetGhostedCellRange(&local_node_info, &user->neighbors, 2, &my_cell_info.zs_cell, &my_cell_info.zm_cell); CHKERRQ(ierr);
    
    /*
    ierr = GetOwnedCellRange(&local_node_info, 0, &my_cell_info.xs_cell, &my_cell_info.xm_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&local_node_info, 1, &my_cell_info.ys_cell, &my_cell_info.ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&local_node_info, 2, &my_cell_info.zs_cell, &my_cell_info.zm_cell); CHKERRQ(ierr);
    */
    // Log the calculated local ownership for debugging purposes.
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Owns cells: i[%d, %d), j[%d, %d), k[%d, %d)\n",
              rank, my_cell_info.xs_cell, my_cell_info.xs_cell + my_cell_info.xm_cell,
              my_cell_info.ys_cell, my_cell_info.ys_cell + my_cell_info.ym_cell,
              my_cell_info.zs_cell, my_cell_info.zs_cell + my_cell_info.zm_cell);

    // --- 3. Allocate and Distribute the Global Map ---
    // Allocate memory for the global map that will hold information from all ranks.
    ierr = PetscMalloc1(size, &user->RankCellInfoMap); CHKERRQ(ierr);

    // Perform the collective communication to gather the `RankCellInfo` struct from every rank.
    // Each rank sends its `my_cell_info` and receives the complete array in `user->RankCellInfoMap`.
    // We use MPI_BYTE to ensure portability across different systems and struct padding.
    ierr = MPI_Allgather(&my_cell_info, sizeof(RankCellInfo), MPI_BYTE,
                         user->RankCellInfoMap, sizeof(RankCellInfo), MPI_BYTE,
                         PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Domain cell decomposition map created and distributed successfully.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BinarySearchInt64"
/**
 * @brief Performs a binary search for a key in a sorted array of PetscInt64.
 *
 * This is a standard binary search algorithm implemented as a PETSc-style helper function.
 * It efficiently determines if a given `key` exists within a `sorted` array.
 *
 * @param[in]  n      The number of elements in the array.
 * @param[in]  arr    A pointer to the sorted array of PetscInt64 values to be searched.
 * @param[in]  key    The PetscInt64 value to search for.
 * @param[out] found  A pointer to a PetscBool that will be set to PETSC_TRUE if the key
 *                    is found, and PETSC_FALSE otherwise.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *
 * @note The input array `arr` **must** be sorted in ascending order for the algorithm
 *       to work correctly.
 */
PetscErrorCode BinarySearchInt64(PetscInt n, const PetscInt64 arr[], PetscInt64 key, PetscBool *found)
{
    PetscInt low = 0, high = n - 1;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- 1. Input Validation ---
    if (!found) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer 'found' is NULL in PetscBinarySearchInt64.");
    }
    if (n > 0 && !arr) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input array 'arr' is NULL for n > 0.");
    }
    
    // Initialize output
    *found = PETSC_FALSE;

    // --- 2. Binary Search Algorithm ---
    while (low <= high) {
        // Use this form to prevent potential integer overflow on very large arrays
        PetscInt mid = low + (high - low) / 2;

        if (arr[mid] == key) {
            *found = PETSC_TRUE; // Key found!
            break;               // Exit the loop
        }
        
        if (arr[mid] < key) {
            low = mid + 1; // Search in the right half
        } else {
            high = mid - 1; // Search in the left half
        }
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


static PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)
{
  PetscInt nidx;
  DMDALocalInfo	info = user->info;

  PetscInt	mx = info.mx, my = info.my;
  
  AO ao;
  DMDAGetAO(user->da, &ao);
  nidx=i+j*mx+k*mx*my;
  
  AOApplicationToPetsc(ao,1,&nidx);
  
  return (nidx);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeDivergence"
PetscErrorCode ComputeDivergence(UserCtx *user )
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt ti = user->simCtx->step;
  
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj, ***nvert,***p;
  Cmpnts	***ucont;
  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscFunctionBeginUser;
  PROFILE_FUNCTION_BEGIN;

  DMDAVecGetArray(fda,user->lUcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);
   DMDAVecGetArray(da, user->P, &p);
  for (k=lzs; k<lze; k++) {
         for (j=lys; j<lye; j++){
	  for (i=lxs; i<lxe; i++) {
	    if (k==10 && j==10 && i==1){
	      LOG_ALLOW(LOCAL,LOG_INFO,"Pressure[10][10][1] =  %f | Pressure[10][10][0] =  %f  \n ",p[k][j][i],p[k][j][i-1]);
              }

	    if (k==10 && j==10 && i==mx-3)
     	      LOG_ALLOW(LOCAL,LOG_INFO,"Pressure[10][10][%d] =  %f | Pressure[10][10][%d] =  %f  \n ",mx-2,p[k][j][mx-2],mx-1,p[k][j][mx-1]);
   		}
	}
    }
   DMDAVecRestoreArray(da, user->P, &p);
 

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x +
		       ucont[k][j][i].y - ucont[k][j-1][i].y +
		       ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	    nvert[k][j+1][i] + nvert[k][j-1][i] +
	    nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	div[k][j][i] = maxdiv;

      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, Div, &div);
  PetscInt MaxFlatIndex;
  
  VecMax(Div, &MaxFlatIndex, &maxdiv);

  LOG_ALLOW(GLOBAL,LOG_INFO,"[Step %d]] The Maximum Divergence is %e at flat index %d.\n",ti,maxdiv,MaxFlatIndex); 

  user->simCtx->MaxDivFlatArg = MaxFlatIndex;
  user->simCtx->MaxDiv = maxdiv;
  
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (Gidx(i,j,k,user) == MaxFlatIndex) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"[Step %d] The Maximum Divergence(%e) is at location [%d][%d][%d]. \n", ti, maxdiv,k,j,i);
	  user->simCtx->MaxDivz = k;
	  user->simCtx->MaxDivy = j;
	  user->simCtx->MaxDivx = i;
	}
      }
    }
  }

  
 DMDAVecRestoreArray(da, user->lNvert, &nvert);
 DMDAVecRestoreArray(fda, user->lUcont, &ucont);
 DMDAVecRestoreArray(da, user->lAj, &aj);
 VecDestroy(&Div);
 
 PROFILE_FUNCTION_END;
 PetscFunctionReturn(0);
}
