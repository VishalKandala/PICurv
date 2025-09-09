 // --- Options from the modern (particle) code ---
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", &simCtx->np, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rstart", &simCtx->StartStep, &simCtx->rstart_flg); CHKERRQ(ierr);
    if(simCtx->rstart_flg) {
        ierr = PetscOptionsGetReal(NULL, NULL, "-ti", &simCtx->StartTime, NULL); CHKERRQ(ierr);
    }
    ierr = PetscOptionsGetInt(NULL,NULL, "-totalsteps", &simCtx->StepsToRun, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", &simCtx->block_number, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", &simCtx->readFields, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-only_setup", &simCtx->OnlySetup, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &simCtx->dt, NULL); CHKERRQ(ierr);

     // --- Options from the legacy (IBM/FSI) code ---
    ierr = PetscOptionsGetInt(NULL, NULL, "-tio", &simCtx->tiout, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-imm", &simCtx->immersed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-inv", &simCtx->invicid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-imp", &simCtx->implicit, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-imp_type", &simCtx->implicit_type, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-imp_MAX_IT", &simCtx->imp_MAX_IT, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-fsi", &simCtx->movefsi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rfsi", &simCtx->rotatefsi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-radi", &simCtx->radi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-inlet", &simCtx->inletprofile, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-str", &simCtx->STRONG_COUPLING, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-rs_fsi", &simCtx->rstart_fsi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-cop", &simCtx->cop, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-fish", &simCtx->fish, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pizza", &simCtx->pizza, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rheology", &simCtx->rheology, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-duplicate", &simCtx->duplicate, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-Pipe", &simCtx->Pipe, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-turbine", &simCtx->turbine, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-fishcyl", &simCtx->fishcyl, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-eel", &simCtx->eel, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-cstart", &simCtx->fish_c, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-wing", &simCtx->wing, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-sediment", &simCtx->sediment, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mhv", &simCtx->MHV, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-hydro", &simCtx->hydro, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-lv", &simCtx->LV, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-reg", &simCtx->regime, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-TwoD", &simCtx->TwoD, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-thin", &simCtx->thin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_z", &simCtx->dgf_z, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_y", &simCtx->dgf_y, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_x", &simCtx->dgf_x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_az", &simCtx->dgf_az, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_ay", &simCtx->dgf_ay, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dgf_ax", &simCtx->dgf_ax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-body", &simCtx->NumberOfBodies, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mframe", &simCtx->moveframe, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rframe", &simCtx->rotateframe, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-blk", &simCtx->blank, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-init1", &simCtx->InitialGuessOne, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-Ogrid", &simCtx->Ogrid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-x_c", &simCtx->CMx_c, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-y_c", &simCtx->CMy_c, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-z_c", &simCtx->CMz_c, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-imp_atol", &simCtx->imp_atol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-imp_rtol", &simCtx->imp_rtol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-imp_stol", &simCtx->imp_stol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-poisson_tol", &simCtx->poisson_tol, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-les", &simCtx->les, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rans", &simCtx->rans, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-wallfunction", &simCtx->wallfunction, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mixed", &simCtx->mixed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-clark", &simCtx->clark, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dynamic_freq", &simCtx->dynamic_freq, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-averaging", &simCtx->averaging, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid1d", &simCtx->grid1d, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-i_periodic", &simCtx->i_periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-j_periodic", &simCtx->j_periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-k_periodic", &simCtx->k_periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pbc_domain", &simCtx->blkpbc, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-testfilter_ik", &simCtx->testfilter_ik, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-testfilter_1d", &simCtx->testfilter_1d, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-poisson", &simCtx->poisson, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-i_homo_filter", &simCtx->i_homo_filter, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-j_homo_filter", &simCtx->j_homo_filter, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-k_homo_filter", &simCtx->k_homo_filter, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-max_cs", &simCtx->max_cs, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-St_exp", &simCtx->St_exp, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-wlngth", &simCtx->wavelength, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-cfl", &simCtx->cfl, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-vnn", &simCtx->vnn, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-chact_leng", &simCtx->chact_leng, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-mg_level", &simCtx->mglevels, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-grid", &simCtx->generate_grid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-block_number", &simCtx->block_number, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-grid_rotation_angle", &simCtx->grid_rotation_angle, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-Croty", &simCtx->Croty, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-Crotz", &simCtx->Crotz, NULL); CHKERRQ(ierr);

    // --- Array options that depend on nblk and generate_grid ---
    // Keep this block separate after the primary flags are parsed.
    PetscInt nblk = simCtx->block_number;
    PetscBool found;
    if (nblk > 0) {
      ierr = PetscMalloc1(nblk, &simCtx->cgrid); CHKERRQ(ierr);
      // Initialize cgrid defaults *before* parsing
      for (PetscInt i=0; i<nblk; i++) simCtx->cgrid[i] = 0;
      ierr = PetscOptionsGetIntArray(NULL, NULL, "-cgrid", simCtx->cgrid, &nblk, &found); CHKERRQ(ierr);
        
      if (simCtx->generate_grid) {
	ierr = PetscMalloc3(nblk, &simCtx->IM_arr, nblk, &simCtx->JM_arr, nblk, &simCtx->KM_arr); CHKERRQ(ierr);
	// It's good practice to set defaults for these arrays too.
	for (PetscInt i=0; i<nblk; i++) { simCtx->IM_arr[i]=10; simCtx->JM_arr[i]=10; simCtx->KM_arr[i]=10; }
	ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", simCtx->IM_arr, &nblk, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", simCtx->JM_arr, &nblk, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", simCtx->KM_arr, &nblk, &found); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-L_x", &simCtx->Lx, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-L_y", &simCtx->Ly, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-L_z", &simCtx->Lz, NULL); CHKERRQ(ierr);
      }
    }

////////////////////////////////////////////////////////////


PetscErrorCode Block_Interface_U(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	hostU;
  Cmpnts ***itfc,***ubcs, ucent;
  PetscReal *hostu,***nvert;

  VecScatter tolocalall;
  PetscErrorCode ierr;

  SimCtx *simCtx = user->simCtx;
  PetscInt block_number = simCtx->block_number;
    
  /* ucat is at cell centers while litfc is now on the cell corners! */
  
  for (bi=0; bi<block_number; bi++) {
    /* hostU is a parallel PETSc vector that will hold vector values 
       in the natural numbering, rather than in the PETSc parallel 
       numbering associated with the DA */
    ierr = DMDACreateNaturalVector(user[bi].fda, &hostU);
    

    // put the center node velocities in the natural ordering in nhostU 
    DMDAGlobalToNaturalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);
    DMDAGlobalToNaturalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);

    /* the velocities at cell centers are saved in user[bi].nhostU
       which is a sequential vector*/
   
    VecScatterCreateToAll(hostU, &tolocalall, &(user[bi].nhostU));

    VecScatterBegin(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		    SCATTER_FORWARD);
    VecScatterEnd(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		  SCATTER_FORWARD);
   
    VecScatterDestroy(&tolocalall);
    VecDestroy(&hostU);
   
  }

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt    lmx, lmy, lmz;
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes
       from the host nodes.
       itfc is the velocity at the cell corners
       hostU is the velocity at the cell centers of the host block */
    /**************************************************************************************************************************/
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
   
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];

      if ((x+y+z<-1.) &&
	  ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {
	itfc[ck][cj][ci].x = 0.;
	itfc[ck][cj][ci].y = 0.;
	itfc[ck][cj][ci].z = 0.;

      } else
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {

	VecGetArray(user[hb].nhostU, &hostu);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].z = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * z); //i+1,j+1,k+1

	VecRestoreArray(user[hb].nhostU, &hostu);
      }
    } // for itfcnumber
  
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMLocalToLocalBegin(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
			user[bi].lItfc);
    DMLocalToLocalEnd(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
		      user[bi].lItfc);
  } // for bi
  //  PetscBarrier(NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");

  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    Cmpnts ***ucont, ***kzet, ***jeta,***icsi;
   
    //    DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecGetArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecGetArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecGetArray(user[bi].fda, user[bi].lCsi, &icsi);

    /**************************************************************************************************************************/
    /* Create boundary condition for flux (ucont) and cell surface
       center velocity (ubcs) from the interpolated velocity (itfc).
       
       itfc is the velocity at the cell surface centers */
    /**************************************************************************************************************************/
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j-1][i  ].x) ;
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i].x*
			      icsi[k][j][i].x +
			      ubcs[k][j][i].y*
			      icsi[k][j][i].y +
			      ubcs[k][j][i].z*
			      icsi[k][j][i].z);

	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
/* 	  ubcs[k][j][i+1].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i+1].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i+1].z = itfc[k][j][i].z; */
	  ubcs[k][j][i+1].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j-1][i  ].x);
	  ubcs[k][j][i+1].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i+1].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i+1].x  *
			      icsi[k][j][i].x +
			      ubcs[k][j][i+1].y  *
			      icsi[k][j][i].y +
			      ubcs[k][j][i+1].z *
			      icsi[k][j][i].z);

	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j  ][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j  ][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j][i].x*
			       jeta[k][j][i].x +
			       ubcs[k][j][i].y *
			       jeta[k][j][i].y +
			       ubcs[k][j][i].z *
			       jeta[k][j][i].z);
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k][j+1][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j+1][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j+1][i].z = itfc[k][j][i].z; */
	  ubcs[k][j+1][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j  ][i-1].x);
	  ubcs[k][j+1][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j  ][i-1].y);
	  ubcs[k][j+1][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j+1][i].x*
			       jeta[k][j  ][i].x +
			       ubcs[k][j+1][i].y*
			       jeta[k][j  ][i].y +
			       ubcs[k][j+1][i].z*
			       jeta[k][j  ][i].z);
	 /*  if (i==1 && k==50)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==100)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==150)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==200)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==50)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==100)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==150)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==200)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
	}
      }
    }


    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j-1][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j-1][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k][j][i].x*
			       kzet[k][j][i].x +
			       ubcs[k][j][i].y*
			       kzet[k][j][i].y +
			       ubcs[k][j][i].z *
			       kzet[k][j][i].z);
	  
	}
      }
    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k+1][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k+1][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k+1][j][i].z = itfc[k][j][i].z; */
	  ubcs[k+1][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j-1][i-1].x);
	  ubcs[k+1][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j-1][i-1].y);
	  ubcs[k+1][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k+1][j][i].x*
			       kzet[k  ][j][i].x +
			       ubcs[k+1][j][i].y*
			       kzet[k  ][j][i].y +
			       ubcs[k+1][j][i].z *
			       kzet[k  ][j][i].z);

	}
      }
    }

    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    // This part is for blanking
    if(blank && bi==0){
           
      DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecGetArray(user[bi].fda, user[bi].lUcont, &ucont);
      for (sb=1; sb<block_number; sb++) {
	PetscInt ip, im, jp, jm, kp, km;
	PetscInt ii, jj, kk;
	  
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    for (i=lxs; i<lxe; i++) {
	      ip = (i<mx-2?(i+1):(i));
	      im = (i>1   ?(i-1):(i));
	      
	      jp = (j<my-2?(j+1):(j));
	      jm = (j>1   ?(j-1):(j));
	      
	      kp = (k<mz-2?(k+1):(k));
	      km = (k>1   ?(k-1):(k));
	      
	      if (((int)(nvert[k][j][i]+0.5) < sb*10) &&
		  ((int)(nvert[k][j][i]+0.5) > sb*10-3) ) {
		// flux in x direction
		kk=k;	jj=j;
		for (ii=im; ii<ip; ii++) {
		  ucent.x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
				    itfc[kk  ][jj-1][ii  ].x +
				    itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk-1][jj-1][ii  ].x) ;
		  ucent.y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
				    itfc[kk  ][jj-1][ii  ].y +
				    itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk-1][jj-1][ii  ].y);
		  ucent.z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
				    itfc[kk  ][jj-1][ii  ].z +
				    itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk-1][jj-1][ii  ].z);
		  ucont[kk][jj][ii].x = (ucent.x*
					 icsi[kk][jj][ii].x +
					 ucent.y*
					 icsi[kk][jj][ii].y +
					 ucent.z*
					 icsi[kk][jj][ii].z);
		}
		// flux in y direction
		kk=k;   ii=i;
		for (jj=jm; jj<jp; jj++) {
		  ucent.x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
				    itfc[kk  ][jj  ][ii-1].x +
				    itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk-1][jj  ][ii-1].x);
		  ucent.y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
				    itfc[kk  ][jj  ][ii-1].y +
				    itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk-1][jj  ][ii-1].y);
		  ucent.z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
				    itfc[kk  ][jj  ][ii-1].z +
				    itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk-1][jj  ][ii-1].z);
		  ucont[kk][jj][ii].y = ( ucent.x*
					  jeta[kk][jj][ii].x +
					  ucent.y *
					  jeta[kk][jj][ii].y +
					  ucent.z *
					  jeta[kk][jj][ii].z);
		}
		// flux in z direction
		jj=j;  ii=i;
		for (kk=km; kk<kp; kk++) {
		  ucent.x = 0.25 * (itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk  ][jj  ][ii-1].x +
				    itfc[kk  ][jj-1][ii  ].x +
				    itfc[kk  ][jj-1][ii-1].x);
		  ucent.y = 0.25 * (itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk  ][jj  ][ii-1].y +
				    itfc[kk  ][jj-1][ii  ].y +
				    itfc[kk  ][jj-1][ii-1].y);
		  ucent.z = 0.25 * (itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk  ][jj  ][ii-1].z +
				    itfc[kk  ][jj-1][ii  ].z +
				    itfc[kk  ][jj-1][ii-1].z);
		  ucont[kk][jj][ii].z = ( ucent.x*
					  kzet[kk][jj][ii].x +
					  ucent.y*
					  kzet[kk][jj][ii].y +
					  ucent.z *
					  kzet[kk][jj][ii].z);
		}
	      }// if (nvert)
	    }
	  }
	}
      } //for sb
      DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecRestoreArray(user[bi].fda, user[bi].lUcont, &ucont);
      //      PetscPrintf(PETSC_COMM_WORLD, "Local to global lUcont _U ");

      DMLocalToGlobalBegin(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);
      DMLocalToGlobalEnd(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);

      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    } // if blank
   
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lCsi, &icsi);
 
  } // bi

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostU);
    Contra2Cart(&(user[bi]));
    GhostNodeVelocity(&user[bi]);
  }

  return(0);      
}


//////////////////////////////////////////////////////////////////////////////

PetscErrorCode PoissonLHSNew(UserCtx *user)
{
  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create local variables to mirror the legacy globals for minimal code changes.
  const PetscInt les = simCtx->les;
  const PetscInt rans = simCtx->rans;
  // --- END CONTEXT ACQUISITION BLOCK ---
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt      IM=user->IM, JM=user->JM, KM=user->KM;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  PetscInt      gxs, gxe, gys, gye, gzs, gze;

  Vec		G11, G12, G13, G21, G22, G23, G31, G32, G33;
  PetscReal	***g11, ***g12, ***g13, ***g21, ***g22, ***g23;
  PetscReal	***g31, ***g32, ***g33;

  PetscReal	***nvert;
  PetscScalar	vol[19];
  PetscInt	idx[19], row;

  PetscInt	i, j, k, N;
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;
  //amir
  // just for turbulent channel flow

  
  if (!user->assignedA) {
    N = mx * my * mz;
    PetscInt M;
    VecGetLocalSize(user->Phi, &M);
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 19, PETSC_NULL, 19, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  MatZeroEntries(user->A);

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);

  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lIZet, &izet);

  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lJZet, &jzet);

  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(da, user->lNvert, &nvert);
 

  VecDuplicate(user->lAj, &G11);
  VecDuplicate(user->lAj, &G12);
  VecDuplicate(user->lAj, &G13);
  VecDuplicate(user->lAj, &G21);
  VecDuplicate(user->lAj, &G22);
  VecDuplicate(user->lAj, &G23);
  VecDuplicate(user->lAj, &G31);
  VecDuplicate(user->lAj, &G32);
  VecDuplicate(user->lAj, &G33);

  DMDAVecGetArray(da, G11, &g11);
  DMDAVecGetArray(da, G12, &g12);
  DMDAVecGetArray(da, G13, &g13);
  DMDAVecGetArray(da, G21, &g21);
  DMDAVecGetArray(da, G22, &g22);
  DMDAVecGetArray(da, G23, &g23);
  DMDAVecGetArray(da, G31, &g31);
  DMDAVecGetArray(da, G32, &g32);
  DMDAVecGetArray(da, G33, &g33);

  for (k=gzs; k<gze; k++) {
    for (j=gys; j<gye; j++) {
      for (i=gxs; i<gxe; i++) {
	if(i>-1 && j>-1 && k>-1 && i<IM+1 && j<JM+1 && k<KM+1){ //Mohsen April 2012
	g11[k][j][i] = (icsi[k][j][i].x * icsi[k][j][i].x +
			icsi[k][j][i].y * icsi[k][j][i].y +
			icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g12[k][j][i] = (ieta[k][j][i].x * icsi[k][j][i].x +
			ieta[k][j][i].y * icsi[k][j][i].y +
			ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
	g13[k][j][i] = (izet[k][j][i].x * icsi[k][j][i].x +
			izet[k][j][i].y * icsi[k][j][i].y +
			izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
					  			 
	g21[k][j][i] = (jcsi[k][j][i].x * jeta[k][j][i].x +
			jcsi[k][j][i].y * jeta[k][j][i].y +
			jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g22[k][j][i] = (jeta[k][j][i].x * jeta[k][j][i].x +
			jeta[k][j][i].y * jeta[k][j][i].y +
			jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
	g23[k][j][i] = (jzet[k][j][i].x * jeta[k][j][i].x +
			jzet[k][j][i].y * jeta[k][j][i].y +
			jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
					  			 
	g31[k][j][i] = (kcsi[k][j][i].x * kzet[k][j][i].x +
			kcsi[k][j][i].y * kzet[k][j][i].y +
			kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g32[k][j][i] = (keta[k][j][i].x * kzet[k][j][i].x +
			keta[k][j][i].y * kzet[k][j][i].y +
			keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	g33[k][j][i] = (kzet[k][j][i].x * kzet[k][j][i].x +
			kzet[k][j][i].y * kzet[k][j][i].y +
			kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
	}	
      }
    }
  }

  PetscInt m;
  PetscReal threshold;
  if (user->thislevel == 2) {
    threshold = 0.1;
  }
  else {
    threshold = 0.1;
  }
  //Mohsen April 2012 

  PetscInt x_str,x_end,y_str,y_end,z_str,z_end;

  if (user->bctype[0]==7) {
    x_end=mx-1;
    x_str=0;
  }
  if (user->bctype[2]==7) {
    y_end=my-1;
    y_str=0;
  }
  if (user->bctype[4]==7) {
    z_end=mz-1;
    z_str=0;
  }

  if (user->bctype[0]!=7){
    x_end=mx-2;
    x_str=1;
  }
  if (user->bctype[2]!=7){
    y_end=my-2;
    y_str=1;
  }
  if (user->bctype[4]!=7){
    z_end=mz-2;
    z_str=1;
  }
 
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	row = Gidx(i, j, k, user);
	if (i == 0 || i == mx-1 || j==0 || j==my-1 || k==0 || k==mz-1) {
	  vol[CP] = 1.;	idx[CP] = Gidx(i, j, k, user);
	  MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	}
	else {
	  if (nvert[k][j][i] > 0.1) { // i, j, k is not a fluid point
	    vol[CP] = 1.; idx[CP] = Gidx(i, j, k, user);
	    MatSetValues(user->A, 1, &row, 1, idx, vol, INSERT_VALUES);
	  }
	  else { // i, j, k is a fluid point
	    for (m=0; m<19; m++) {
	      vol[m] = 0.;
	    }
	    /* Contribution from i+1 - i */
	    if (nvert[k][j][i+1] < threshold && i != x_end) { // i+1, j, k is a fluid point
	      /* dpdc{i} = (p_{i+1} - p_{i}) * g11_{i} */
	      vol[CP] -= g11[k][j][i]; //i, j, k
	      vol[EP] += g11[k][j][i]; // i+1, j, k

	      /* dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1})
	       * 0.25 * g12[k][j][i] */
	      if ((j == my-2 && user->bctype[2]!=7) || nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] += g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] += g12[k][j][i] * 0.5; //i+1, j, k
		  vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
		  vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
		  vol[CP] += g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] += g12[k][j][i] * 0.5; //i+1, j, k
		  vol[SP] -= g12[k][j][i] * 0.5; //i, j-1, k
		  vol[SE] -= g12[k][j][i] * 0.5; //i+1, j-1, k
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5; //i, j+1, k
		  vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
		  vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5; //i, j+1, k
		  vol[NE] += g12[k][j][i] * 0.5; //i+1, j+1, k
		  vol[CP] -= g12[k][j][i] * 0.5; //i, j, k
		  vol[EP] -= g12[k][j][i] * 0.5; //i+1, j, k
		}
	      }
	      else {
		vol[NP] += g12[k][j][i] * 0.25; // i, j+1, k
		vol[NE] += g12[k][j][i] * 0.25; // i+1, j+1, k
		vol[SP] -= g12[k][j][i] * 0.25; // i, j-1, k
		vol[SE] -= g12[k][j][i] * 0.25; // i+1, j-1, k
	      }
	      
	      /* dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] += g13[k][j][i] * 0.5; // i, j, k
		  vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
		  vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
		  vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
		  vol[CP] += g13[k][j][i] * 0.5; // i, j, k
		  vol[EP] += g13[k][j][i] * 0.5; // i+1, j, k
		  vol[BP] -= g13[k][j][i] * 0.5; // i, j, k-1
		  vol[BE] -= g13[k][j][i] * 0.5; // i+1, j, k-1
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7) || nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
		  vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
		  vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
		  vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5;  // i, j, k+1
		  vol[TE] += g13[k][j][i] * 0.5; // i+1, j, k+1
		  vol[CP] -= g13[k][j][i] * 0.5;  // i, j, k
		  vol[EP] -= g13[k][j][i] * 0.5;  // i+1, j, k
		}
	      }
	      else {
		vol[TP] += g13[k][j][i] * 0.25; //i, j, k+1
		vol[TE] += g13[k][j][i] * 0.25; //i+1, j, k+1
		vol[BP] -= g13[k][j][i] * 0.25; //i, j, k-1
		vol[BE] -= g13[k][j][i] * 0.25; //i+1, j, k-1
	      }
	    }  // end of i+1 - i

	    /* Contribution from i - i-1 */
	    if (nvert[k][j][i-1] < threshold && i != x_str) { // i-1, j, k is a fluid point
	      /* -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i} */
	      vol[CP] -= g11[k][j][i-1];  //i, j, k
	      vol[WP] += g11[k][j][i-1];  //i-1, j, k

	      /* -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1})
	       * 0.25 * g12[k][j][i-1] */
	      if ((j == my-2 && user->bctype[2]!=7) || nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
		  vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
		  vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g12[k][j][i-1] * 0.5; // i-1, j, k
		  vol[SP] += g12[k][j][i-1] * 0.5; //i, j-1, k
		  vol[SW] += g12[k][j][i-1] * 0.5; // i-1, j-1, k
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; // i, j+1, k
		  vol[NW] -= g12[k][j][i-1] * 0.5; // i-1, j+1, k
		  vol[CP] += g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g12[k][j][i-1] * 0.5; // i-1, j, k
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; // i, j+1, k
		  vol[NW] -= g12[k][j][i-1] * 0.5; // i-1, j+1, k
		  vol[CP] += g12[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g12[k][j][i-1] * 0.5; // i-1, j, k
		}
	      }
	      else {
		vol[NP] -= g12[k][j][i-1] * 0.25; // i, j+1, k
		vol[NW] -= g12[k][j][i-1] * 0.25; //i-1, j+1, k
		vol[SP] += g12[k][j][i-1] * 0.25; // i, j-1, k
		vol[SW] += g12[k][j][i-1] * 0.25; // i-1, j-1, k
	      }

	      /* -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1})
	       * 0.25 * g13[k][j][i] */
	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
		  vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
		  vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] -= g13[k][j][i-1] * 0.5; // i-1, j, k
		  vol[BP] += g13[k][j][i-1] * 0.5; // i, j, k-1
		  vol[BW] += g13[k][j][i-1] * 0.5; // i-1, j, k-1
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7) || nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; // i, j, k+1
		  vol[TW] -= g13[k][j][i-1] * 0.5; //i-1, j, k+1
		  vol[CP] += g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g13[k][j][i-1] * 0.5; //i-1, j, k
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; // i, j, k+1
		  vol[TW] -= g13[k][j][i-1] * 0.5; //i-1, j, k+1
		  vol[CP] += g13[k][j][i-1] * 0.5; // i, j, k
		  vol[WP] += g13[k][j][i-1] * 0.5; //i-1, j, k
		}
	      }
	      else {
		vol[TP] -= g13[k][j][i-1] * 0.25;  // i, j, k+1
		vol[TW] -= g13[k][j][i-1] * 0.25; // i-1, j, k+1
		vol[BP] += g13[k][j][i-1] * 0.25;  // i, j, k-1
		vol[BW] += g13[k][j][i-1] * 0.25; // i-1, j, k-1
	      }
	    } // end of i - i-1

	    /* Contribution from j+1 - j */
	    if (nvert[k][j+1][i] < threshold && j != y_end) {
	      /* dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) *
		 0.25 */
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] += g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
		  vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
		  vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
		  vol[CP] += g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] += g21[k][j][i] * 0.5; // i, j+1, k
		  vol[WP] -= g21[k][j][i] * 0.5; // i-1, j, k
		  vol[NW] -= g21[k][j][i] * 0.5; // i-1, j+1, k
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7) || nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; // i+1, j, k
		  vol[NE] += g21[k][j][i] * 0.5; // i+1, j+1, k
		  vol[CP] -= g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] -= g21[k][j][i] * 0.5; // i, j+1, k
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 &&  nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; // i+1, j, k
		  vol[NE] += g21[k][j][i] * 0.5; // i+1, j+1, k
		  vol[CP] -= g21[k][j][i] * 0.5; // i, j, k
		  vol[NP] -= g21[k][j][i] * 0.5; // i, j+1, k
		}
	      }
	      else {
		vol[EP] += g21[k][j][i] * 0.25; //i+1, j, k
		vol[NE] += g21[k][j][i] * 0.25; //i+1, j+1, k
		vol[WP] -= g21[k][j][i] * 0.25; //i-1, j, k
		vol[NW] -= g21[k][j][i] * 0.25; //i-1, j+1, k
	      }
	     
	      /* dpde{j} = (p{j+1} - p{j}) * g22[k][j][i] */
	      vol[CP] -= g22[k][j][i];
	      vol[NP] += g22[k][j][i];

	      /* dpdz{j} = (p{j, k+1}+p{j+1,k+1} - p{j,k-1}-p{j+1,k-1}) *0.25*/
	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] += g23[k][j][i] * 0.5; //i,j,k
		  vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
		  vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
		  vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
		}
	      }
	      else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[CP] += g23[k][j][i] * 0.5; //i,j,k
		  vol[NP] += g23[k][j][i] * 0.5; //i, j+1, k
		  vol[BP] -= g23[k][j][i] * 0.5;//i, j, k-1
		  vol[BN] -= g23[k][j][i] * 0.5;//i, j+1, k-1
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; //i, j, k+1
		  vol[TN] += g23[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g23[k][j][i] * 0.5;//i, j, k
		  vol[NP] -= g23[k][j][i] * 0.5;//i, j+1, k
		}
	      }
	      else if ((k == 1 || k==mz-2 ) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; //i, j, k+1
		  vol[TN] += g23[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g23[k][j][i] * 0.5;//i, j, k
		  vol[NP] -= g23[k][j][i] * 0.5;//i, j+1, k
		}
	      }
	      else {
		vol[TP] += g23[k][j][i] * 0.25; // i, j, k+1
		vol[TN] += g23[k][j][i] * 0.25; // i, j+1, k+1
		vol[BP] -= g23[k][j][i] * 0.25; // i, j, k-1
		vol[BN] -= g23[k][j][i] * 0.25; // i, j+1, k-1
	      }
	    } // End of j+1 - j

	    /* Contribution j - j-1 */
	    if (nvert[k][j-1][i] < threshold && j != y_str) {
	      /* -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) *
		 0.25 * g21[k][j-1][i] */
	      if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
		  vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
		  vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
		  vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
		}
	      }
	      else  if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1) {
		  vol[CP] -= g21[k][j-1][i] * 0.5;// i, j, k
		  vol[SP] -= g21[k][j-1][i] * 0.5;// i, j-1, k
		  vol[WP] += g21[k][j-1][i] * 0.5;// i-1, j, k
		  vol[SW] += g21[k][j-1][i] * 0.5;// i-1, j-1, k
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5;//i+1, j, k
		  vol[SE] -= g21[k][j-1][i] * 0.5;//i+1, j-1, k
		  vol[CP] += g21[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g21[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5;//i+1, j, k
		  vol[SE] -= g21[k][j-1][i] * 0.5;//i+1, j-1, k
		  vol[CP] += g21[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g21[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[EP] -= g21[k][j-1][i] * 0.25;// i+1, j, k
		vol[SE] -= g21[k][j-1][i] * 0.25;// i+1, j-1, k
		vol[WP] += g21[k][j-1][i] * 0.25;// i-1, j, k
		vol[SW] += g21[k][j-1][i] * 0.25;// i-1, j-1, k
	      }
	      
	      /* -dpde{j-1} = -(p{j} - p{j-1}) * g22[k][j-1][i] */
	      vol[CP] -= g22[k][j-1][i];
	      vol[SP] += g22[k][j-1][i];

	      /* -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) *
		 0.25 * g23[k][j-1][i] */
	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
		  vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
		  vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 ) {
		  vol[CP] -= g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] -= g23[k][j-1][i] * 0.5;//i, j-1, k
		  vol[BP] += g23[k][j-1][i] * 0.5;//i, j, k-1
		  vol[BS] += g23[k][j-1][i] * 0.5;//i, j-1, k-1
		}
	      }
	     
	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5;//i, j, k+1
		  vol[TS] -= g23[k][j-1][i] * 0.5;//i, j-1, k+1
		  vol[CP] += g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g23[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5;//i, j, k+1
		  vol[TS] -= g23[k][j-1][i] * 0.5;//i, j-1, k+1
		  vol[CP] += g23[k][j-1][i] * 0.5;//i, j, k
		  vol[SP] += g23[k][j-1][i] * 0.5;//i, j-1, k
		}
	      }
	      else {
		vol[TP] -= g23[k][j-1][i] * 0.25;//i, j, k+1
		vol[TS] -= g23[k][j-1][i] * 0.25;//i, j-1, k+1
		vol[BP] += g23[k][j-1][i] * 0.25;//i, j, k-1
		vol[BS] += g23[k][j-1][i] * 0.25;//i, j-1, k-1
	      }
	    } // End of j - j-1

	    /* contribution from k+1 - k */
	    if (nvert[k+1][j][i] < threshold && k != z_end) {
	      /* dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) *
		 0.25 * g31[k][j][i] */
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] += g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
		  vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
		  vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
		  vol[CP] += g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] += g31[k][j][i] * 0.5;//i, j, k+1
		  vol[WP] -= g31[k][j][i] * 0.5;//i-1, j, k
		  vol[TW] -= g31[k][j][i] * 0.5;//i-1, j, k+1
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5;//i+1, j, k
		  vol[TE] += g31[k][j][i] * 0.5;//i+1, j, k+1
		  vol[CP] -= g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g31[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5;//i+1, j, k
		  vol[TE] += g31[k][j][i] * 0.5;//i+1, j, k+1
		  vol[CP] -= g31[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g31[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[EP] += g31[k][j][i] * 0.25;//i+1, j, k
		vol[TE] += g31[k][j][i] * 0.25;//i+1, j, k+1
		vol[WP] -= g31[k][j][i] * 0.25;//i-1, j, k
		vol[TW] -= g31[k][j][i] * 0.25;//i-1, j, k+1
	      }

	      /* dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 
		 0.25 * g32[k][j][i] */
	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] += g32[k][j][i] * 0.5;//i, j,k
		  vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
		  vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
		  vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[CP] += g32[k][j][i] * 0.5;//i, j,k
		  vol[TP] += g32[k][j][i] * 0.5;//i, j, k+1
		  vol[SP] -= g32[k][j][i] * 0.5;//i, j-1, k
		  vol[TS] -= g32[k][j][i] * 0.5;//i, j-1, k+1
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5;//i, j+1, k
		  vol[TN] += g32[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g32[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g32[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5;//i, j+1, k
		  vol[TN] += g32[k][j][i] * 0.5;//i, j+1, k+1
		  vol[CP] -= g32[k][j][i] * 0.5;//i, j, k
		  vol[TP] -= g32[k][j][i] * 0.5;//i, j, k+1
		}
	      }
	      else {
		vol[NP] += g32[k][j][i] * 0.25;//i, j+1, k
		vol[TN] += g32[k][j][i] * 0.25;//i, j+1, k+1
		vol[SP] -= g32[k][j][i] * 0.25;//i, j-1, k
		vol[TS] -= g32[k][j][i] * 0.25;//i, j-1, k+1
	      }

	      /* dpdz{k} = p{k+1} - p{k} */
	      vol[CP] -= g33[k][j][i]; //i, j, k
	      vol[TP] += g33[k][j][i]; //i, j, k+1
	    } // End of k+1 - k

	    /* Contribution from k - k-1 */
	    if (nvert[k-1][j][i] < threshold && k != z_str) {
	      /* -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) *
		 0.25 * g31[k-1][j][i] */
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
		  vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
		  vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1) {
		  vol[CP] -= g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] -= g31[k-1][j][i] * 0.5;//i, j, k-1
		  vol[WP] += g31[k-1][j][i] * 0.5;//i-1, j, k
		  vol[BW] += g31[k-1][j][i] * 0.5;//i-1, j, k-1
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5;//i+1, j, k
		  vol[BE] -= g31[k-1][j][i] * 0.5;//i+1, j, k-1
		  vol[CP] += g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g31[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5;//i+1, j, k
		  vol[BE] -= g31[k-1][j][i] * 0.5;//i+1, j, k-1
		  vol[CP] += g31[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g31[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[EP] -= g31[k-1][j][i] * 0.25;//i+1, j, k
		vol[BE] -= g31[k-1][j][i] * 0.25;//i+1, j, k-1
		vol[WP] += g31[k-1][j][i] * 0.25;//i-1, j, k
		vol[BW] += g31[k-1][j][i] * 0.25;//i-1, j, k-1
	      }
	      
	      /* -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) * 
		 0.25 * g32[k-1][j][i] */
	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
		  vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
		  vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
		  vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1) {
		  vol[CP] -= g32[k-1][j][i] * 0.5;//i, j,k
		  vol[BP] -= g32[k-1][j][i] * 0.5;//i, j, k-1
		  vol[SP] += g32[k-1][j][i] * 0.5;//i, j-1, k 
		  vol[BS] += g32[k-1][j][i] * 0.5;//i, j-1, k-1
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5;//i, j+1, k
		  vol[BN] -= g32[k-1][j][i] * 0.5;//i, j+1, k-1
		  vol[CP] += g32[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g32[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5;//i, j+1, k
		  vol[BN] -= g32[k-1][j][i] * 0.5;//i, j+1, k-1
		  vol[CP] += g32[k-1][j][i] * 0.5;//i, j, k
		  vol[BP] += g32[k-1][j][i] * 0.5;//i, j, k-1
		}
	      }
	      else {
		vol[NP] -= g32[k-1][j][i] * 0.25;//i, j+1, k
		vol[BN] -= g32[k-1][j][i] * 0.25;//i, j+1, k-1
		vol[SP] += g32[k-1][j][i] * 0.25;//i, j-1, k
		vol[BS] += g32[k-1][j][i] * 0.25;//i, j-1, k-1
	      }
	      
	      /* -dpdz{k-1} = -(p{k} - p{k-1}) * g33[k-1][j][i] */
	      vol[CP] -= g33[k-1][j][i]; // i, j, k
	      vol[BP] += g33[k-1][j][i]; //i, j, k-1
	    } // End of k - k-1
	    for (m=0; m<19; m++) {
	      vol[m] *= -aj[k][j][i];
	    }
	    //Mohsen April 2012
	    idx[CP] = Gidx(i  , j  , k  , user);
	    if (user->bctype[0]==7 && i==mx-2)
	      idx[EP] = Gidx(1, j  , k  , user);
	    else
	      idx[EP] = Gidx(i+1, j  , k  , user);
	    if (user->bctype[0]==7 && i==1)
	      idx[WP] = Gidx(mx-2, j  , k  , user);
	    else
	      idx[WP] = Gidx(i-1, j  , k  , user);
	    if (user->bctype[2]==7 && j==my-2)
	      idx[NP] = Gidx(i, 1  , k  , user);
	    else
	      idx[NP] = Gidx(i  , j+1, k  , user);
	    if (user->bctype[2]==7 && j==1)
	      idx[SP] = Gidx(i, my-2  , k  , user);
	    else
	      idx[SP] = Gidx(i  , j-1, k  , user);
	    
	    if (user->bctype[4]==7 && k==mz-2){
	      idx[TP] = Gidx(i  , j  , 1, user);
	    } else
	      idx[TP] = Gidx(i  , j  , k+1, user);
   	    if (user->bctype[4]==7 && k==1)
	      idx[BP] = Gidx(i  , j  , mz-2, user);
	    else
	      idx[BP] = Gidx(i  , j  , k-1, user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==mx-2 && j==my-2)
	      idx[NE] = Gidx(1, 1, k  , user);
	    else if (user->bctype[0]==7 && i==mx-2)
	      idx[NE] = Gidx(1, j+1, k  , user);
	    else if (user->bctype[2]==7 && j==my-2)
	      idx[NE] = Gidx(i+1, 1, k  , user);
	    else
	      idx[NE] = Gidx(i+1, j+1, k  , user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==mx-2 && j==1)
	      idx[SE] = Gidx(1, my-2, k  , user);
	    else if (user->bctype[0]==7 && i==mx-2)
	      idx[SE] = Gidx(1, j-1, k  , user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[SE] = Gidx(i+1, my-2, k  , user);
	    else
	      idx[SE] = Gidx(i+1, j-1, k  , user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==1 && j==my-2)
	      idx[NW] = Gidx(mx-2, 1, k  , user);
	    else if (user->bctype[0]==7 && i==1)
	      idx[NW] = Gidx(mx-2, j+1, k  , user);
	    else if (user->bctype[2]==7 && j==my-2)
	      idx[NW] = Gidx(i-1, 1, k  , user);
	    else
	      idx[NW] = Gidx(i-1, j+1, k  , user);
	    if (user->bctype[0]==7 && user->bctype[2]==7 && i==1 && j==1)
	      idx[SW] = Gidx(mx-2, my-2, k  , user);
	    else if (user->bctype[0]==7 && i==1)
	      idx[SW] = Gidx(mx-2, j-1, k  , user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[SW] = Gidx(i-1, my-2, k  , user);
	    else
	      idx[SW] = Gidx(i-1, j-1, k  , user);
	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==my-2 && k==mz-2)
	      idx[TN] = Gidx(i  , 1, 1, user);
	    else if (user->bctype[2]==7 && j==my-2)
	      idx[TN] = Gidx(i  , 1, k+1, user);
	    else if (user->bctype[4]==7 && k==mz-2)
	      idx[TN] = Gidx(i  , j+1, 1, user);
	    else
	      idx[TN] = Gidx(i  , j+1, k+1, user);
	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==my-2 && k==1)
	      idx[BN] = Gidx(i  , 1, mz-2 , user);
	    else if(user->bctype[2]==7 && j==my-2)
	      idx[BN] = Gidx(i  , 1, k-1, user);
	    else if (user->bctype[4]==7 && k==1)
	      idx[BN] = Gidx(i  , j+1, mz-2, user);
	    else
	      idx[BN] = Gidx(i  , j+1, k-1, user);
	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==1 && k==mz-2)
	      idx[TS] = Gidx(i  , my-2, 1, user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[TS] = Gidx(i  , my-2, k+1, user);
	    else if (user->bctype[4]==7 && k==mz-2)
	      idx[TS] = Gidx(i  , j-1, 1, user);
	    else
	      idx[TS] = Gidx(i  , j-1, k+1, user);
       	    if (user->bctype[2]==7 && user->bctype[4]==7 && j==1 && k==1)
	      idx[BS] = Gidx(i  , my-2, mz-2, user);
	    else if (user->bctype[2]==7 && j==1)
	      idx[BS] = Gidx(i  , my-2, k-1, user);
	    else if (user->bctype[4]==7 && k==1)
	      idx[BS] = Gidx(i  , j-1, mz-2, user);
	    else
	      idx[BS] = Gidx(i  , j-1, k-1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==mx-2 && k==mz-2)
	      idx[TE] = Gidx(1, j  , 1, user);
	    else if(user->bctype[0]==7 && i==mx-2)
	      idx[TE] = Gidx(1, j  , k+1, user);
	    else if(user->bctype[4]==7 && k==mz-2)
	      idx[TE] = Gidx(i+1, j  , 1, user);
	    else
	      idx[TE] = Gidx(i+1, j  , k+1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==mx-2 && k==1)
	      idx[BE] = Gidx(1, j  , mz-2, user);
	    else if(user->bctype[0]==7 && i==mx-2)
	      idx[BE] = Gidx(1, j  , k-1, user);
	    else if(user->bctype[4]==7 && k==1)
	      idx[BE] = Gidx(i+1, j  , mz-2, user);
	    else
	      idx[BE] = Gidx(i+1, j  , k-1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==1 && k==mz-2)
	      idx[TW] = Gidx(mx-2, j  , 1, user);
	    else if(user->bctype[0]==7 && i==1)
	      idx[TW] = Gidx(mx-2, j  , k+1, user);
	    else if (user->bctype[4]==7 && k==mz-2)
	      idx[TW] = Gidx(i-1, j  , 1, user);
	    else
	      idx[TW] = Gidx(i-1, j  , k+1, user);
	    if (user->bctype[0]==7 && user->bctype[4]==7 && i==1 && k==1)
	      idx[BW] = Gidx(mx-2, j  , mz-2, user);
	    else if (user->bctype[0]==7 && i==1)
	      idx[BW] = Gidx(mx-2, j  , k-1, user);
	    else if (user->bctype[4]==7 && k==1)
	      idx[BW] = Gidx(i-1, j  , mz-2, user);
	    else
	      idx[BW] = Gidx(i-1, j  , k-1, user);
	    MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);

	  } // End of fluid point
	} // End of interial points
      }
    }
  }
  
  
  
  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  DMDAVecRestoreArray(da, G11, &g11);
  DMDAVecRestoreArray(da, G12, &g12);
  DMDAVecRestoreArray(da, G13, &g13);
  DMDAVecRestoreArray(da, G21, &g21);
  DMDAVecRestoreArray(da, G22, &g22);
  DMDAVecRestoreArray(da, G23, &g23);
  DMDAVecRestoreArray(da, G31, &g31);
  DMDAVecRestoreArray(da, G32, &g32);
  DMDAVecRestoreArray(da, G33, &g33);
  
  VecDestroy(&G11);
  VecDestroy(&G12);
  VecDestroy(&G13);
  VecDestroy(&G21);
  VecDestroy(&G22);
  VecDestroy(&G23);
  VecDestroy(&G31);
  VecDestroy(&G32);
  VecDestroy(&G33);


  //  VecCopy(user->Phi, user->P);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);

  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);

  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);

  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);

  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  return 0;
}

/////////////////////////////////////////


PetscErrorCode UpdatePressure(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
 
  PetscInt	lxs, lys, lzs, lxe, lye, lze,gxs, gxe, gys, gye, gzs, gze;
  PetscInt	i, j, k;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  //amir
  Cmpnts ***cent,***ucont,***uch,***kzet;
  PetscReal ***p, ***phi, ***nvert, ***lp,***lphi,***kaj,dpdz ;

  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(da, user->Phi, &phi);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->Cent, &cent);
  

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
		p[k][j][i] += phi[k][j][i];
      }
    }
  }
  //amir    
	



  DMDAVecRestoreArray(da, user->Phi, &phi);
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->Cent, &cent);

  //Mohsen Aug 2012
  //Updating pressure at boundaries for Periodic BC's
  //Dec 2015
  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP,&lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    DMDAVecGetArray(da, user->P,&p);
    DMDAVecGetArray(da, user->Phi,&phi);
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if( k>0 && k<user->KM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k][j][i-2];
	      phi[k][j][i]=lphi[k][j][i-2];
	    }
	  }
	}
      }
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if( k>0 && k<user->KM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k][j][i+2];
	      phi[k][j][i]=lphi[k][j][i+2];	
	    }
	  }
	}
      }
    }
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMDAVecRestoreArray(da, user->P,&p);
    DMDAVecRestoreArray(da, user->Phi,&phi);

   /*  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*     DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi); */
/*     DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi); */
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    DMDAVecGetArray(da, user->P,&p);
    DMDAVecGetArray(da, user->Phi,&phi);
    
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if( k>0 && k<user->KM && i>0 && i<user->IM){
	      p[k][j][i]=lp[k][j-2][i];
	      phi[k][j][i]=lphi[k][j-2][i];
	    }
	  }
	}
      }
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if( k>0 && k<user->KM && i>0 && i<user->IM){
	      p[k][j][i]=lp[k][j+2][i];
	      phi[k][j][i]=lphi[k][j+2][i];
	    }
	  }
	}
      }
    }
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMDAVecRestoreArray(da, user->P,&p);
    DMDAVecRestoreArray(da, user->Phi,&phi);

  /*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*     DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi); */
/*     DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi); */

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    DMDAVecGetArray(da, user->P,&p);
    DMDAVecGetArray(da, user->Phi,&phi);

    if ((user->bctype[4]==7 || user->bctype[5]==7)){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if( i>0 && i<user->IM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k-2][j][i];
	      phi[k][j][i]=lphi[k-2][j][i];
	    }
	  }
	}
      }
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if( i>0 && i<user->IM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k+2][j][i];
	      phi[k][j][i]=lphi[k+2][j][i];
	    }
	  }
	}
      }
    } 
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMDAVecRestoreArray(da, user->P,&p);
    DMDAVecRestoreArray(da, user->Phi,&phi);

  /*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */
   
/*     DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi); */
/*     DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi); */

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);

    DMDAVecGetArray(da, user->P,&p);
    DMDAVecGetArray(da, user->Phi,&phi);
  
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    p[k][j][i]=lp[k][j][i-2];
	    phi[k][j][i]=lphi[k][j][i-2];
	  }
	}
      }
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    p[k][j][i]=lp[k][j][i+2];
	    phi[k][j][i]=lphi[k][j][i+2];
	  }
	}
      }
    }

     DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMDAVecRestoreArray(da, user->P,&p);
    DMDAVecRestoreArray(da, user->Phi,&phi);

   /*  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*     DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi); */
/*     DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi); */

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);
     
    DMDAVecGetArray(da, user->P,&p);
    DMDAVecGetArray(da, user->Phi,&phi);
 
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    p[k][j][i]=lp[k][j+2][i];
	    phi[k][j][i]=lphi[k][j+2][i];
	  }
	}
      }
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    p[k][j][i]=lp[k][j-2][i];
	    phi[k][j][i]=lphi[k][j-2][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMDAVecRestoreArray(da, user->P,&p);
    DMDAVecRestoreArray(da, user->Phi,&phi);

  /*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*     DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi); */
/*     DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi); */

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);
   
    DMDAVecGetArray(da, user->P,&p);
    DMDAVecGetArray(da, user->Phi,&phi);




    if ((user->bctype[4]==7 || user->bctype[5]==7)){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    p[k][j][i]=lp[k+2][j][i];
	    phi[k][j][i]=lphi[k+2][j][i];
	  }
	}
      }
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    p[k][j][i]=lp[k-2][j][i];
	    phi[k][j][i]=lphi[k-2][j][i];
	  }
	}
      }
    }
/*
       for (k=lzs; k<lze; k++) {
         for (j=lys; j<lye; j++){
	  for (i=lxs; i<lxe; i++) {
	    if (k==10 && j==10 && i==1){
     	     PetscPrintf(PETSC_COMM_SELF,  "%f  %f  \n ",p[k][j][i],p[k][j][i-1]);
              }

	    if (k==10 && j==10 && i==mx-3)
     	      PetscPrintf(PETSC_COMM_SELF,  "%f  %f  \n ",p[k][j][mx-2],p[k][j][mx-1]);
   		}
	}
    }
*/
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);

    DMDAVecRestoreArray(da, user->P,&p);
    DMDAVecRestoreArray(da, user->Phi,&phi);

    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);
   /*  DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */
/*     DMLocalToGlobalBegin(da, user->lPhi, INSERT_VALUES, user->Phi); */
/*     DMLocalToGlobalEnd(da, user->lPhi, INSERT_VALUES, user->Phi); */
  }


  return 0;
}

//////////////////////////////////////////////////////////////

PetscErrorCode Projection(UserCtx *user)
{

  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create local variables to mirror the legacy globals for minimal code changes.
  const PetscInt les = simCtx->les;
  const PetscInt rans = simCtx->rans;
  const PetscInt ti = simCtx->step; // Assuming simCtx->step is the new integer time counter
  const PetscInt TwoD = simCtx->TwoD;
  const PetscInt channelz = simCtx->channelz;
  const PetscInt fish = simCtx->fish;
  const PetscReal dt = simCtx->dt;
  const PetscReal st = simCtx->st;
  // --- END CONTEXT ACQUISITION BLOCK ---

  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***aj, ***iaj, ***jaj, ***kaj;
  PetscReal	***p;

  Cmpnts	***ucont;
  PetscReal	***nvert;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  
  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lIZet, &izet);

  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lJZet, &jzet);

  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->lIAj, &iaj);
  DMDAVecGetArray(da, user->lJAj, &jaj);
  DMDAVecGetArray(da, user->lKAj, &kaj);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lPhi, &p);
  
  DMDAVecGetArray(fda, user->Ucont, &ucont);


  PetscReal dpdc, dpde, dpdz;
  PetscInt i_end,j_end,k_end;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	// Mohsen Aug 2012
	
	if (user->bctype[0]==7) i_end=mx-1;
	else i_end=mx-2;
	
	if (i<i_end) {
	  
	  dpdc = p[k][j][i+1] - p[k][j][i];
	  
	  dpde = 0.;
	  dpdz = 0.;
	 
	  if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	    if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	      dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		      p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	    }
	  }
	  else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	    if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
	      dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		      p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	    }
	  }
	  else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	    if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	      dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		      p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	    }
	  }
	  else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	    if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	      dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		      p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	    }
	  }
	  else {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	  }

	  if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	    if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	      dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		      p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	    }
	  }
	  else  if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	    if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
	      dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		      p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	    }
	  }
	  else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	    if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	      dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		      p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	    }
	  }
	  else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	    if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	      dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		      p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	    }
	  }
	  else {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	  }

	  if (!(nvert[k][j][i] + nvert[k][j][i+1])) {
	    ucont[k][j][i].x -= 
	      (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		       icsi[k][j][i].y * icsi[k][j][i].y +
		       icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	       dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		       ieta[k][j][i].y * icsi[k][j][i].y +
		       ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] + 
	       dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		       izet[k][j][i].y * icsi[k][j][i].y +
		       izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i])
	      * dt * st / COEF_TIME_ACCURACY;
	 
	  }
	}

	if (user->bctype[2]==7) j_end=my-1;
	else j_end=my-2;
	
	if (j<j_end) {
	 
	  dpdc = 0.;
	  dpdz = 0.;

	  if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	    if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	      dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		      p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	    }
	  }
	  else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	    if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
	      dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		      p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	    }
	  }
	  else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	    if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	      dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		      p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	    }
	  }
	  else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	    if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	      dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		      p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	    }
	  }
	  else {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	  }

	  dpde = p[k][j+1][i] - p[k][j][i];

	  if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	    if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	      dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		      p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	    }
	  }
	  else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	    if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
	      dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		      p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	    }
	  }
	  else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	    if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	      dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		      p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	    }
	  }
	  else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	    if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	      dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		      p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	    }
	  }
	  else {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	  }

	  if (!(nvert[k][j][i] + nvert[k][j+1][i])) {
	    ucont[k][j][i].y -=
	      (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		       jcsi[k][j][i].y * jeta[k][j][i].y +
		       jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	       dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		       jeta[k][j][i].y * jeta[k][j][i].y +
		       jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	       dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		       jzet[k][j][i].y * jeta[k][j][i].y +
		       jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i])
	      * dt * st / COEF_TIME_ACCURACY;
	  }
	}
	//amir
	if (user->bctype[4]==7) k_end=mz-1;
	else k_end=mz-2;

	if (k < k_end) {
	 
	  dpdc = 0.;
	  dpde = 0.;
	  
	  if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	    if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	      dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		      p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	    }
	  }
	  else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	    if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
	      dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		      p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	    }
	  }
	  else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	    if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	      dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		      p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	    }
	  }
	  else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	    if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	      dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		      p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	    }
	  }
	  else {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	  }
	  
	  if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	    if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	      dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		      p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	    }
	  }
	  else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	    if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
	      dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		      p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	    }
	  }
	  else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	    if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	      dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		      p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	    }
	  }
	  else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	    if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	      dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		      p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	    }
	  }
	  else {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	  }

	  dpdz = p[k+1][j][i] - p[k][j][i];
	  if (!(nvert[k][j][i] + nvert[k+1][j][i])) {
	    
	    ucont[k][j][i].z -=
	      (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		       kcsi[k][j][i].y * kzet[k][j][i].y +
		       kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	       dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		       keta[k][j][i].y * kzet[k][j][i].y +
		       keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	       dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		       kzet[k][j][i].y * kzet[k][j][i].y +
		       kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i])
	      * dt * st / COEF_TIME_ACCURACY;
	   
	  }
	}
      }
    }
  }

  // Mohsen Sep 2012
  // Update velocity at boundary for periodic BC's//
  // i-direction

  if (user->bctype[0]==7 && xs==0){
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	i=xs;
	
	dpdc = p[k][j][i+1] - p[k][j][i];
	
	dpde = 0.;
	dpdz = 0.;
		
	if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		  p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	}
	
	if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else  if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		  p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	}
	
	
	
	if (!(nvert[k][j][i] + nvert[k][j][i+1])) {
	  ucont[k][j][i].x -=
	    (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		     icsi[k][j][i].y * icsi[k][j][i].y +
		     icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	     dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		     ieta[k][j][i].y * icsi[k][j][i].y +
		     ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	     dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		     izet[k][j][i].y * icsi[k][j][i].y +
		     izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i])
	    * dt * st / COEF_TIME_ACCURACY;
	  
	}
      }
    }
  }

 
/*     // j-direction */

  if (user->bctype[2]==7 && ys==0){
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	j=ys;
	
	dpdc = 0.;
	dpdz = 0.;
	if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		  p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	}
	
	dpde = p[k][j+1][i] - p[k][j][i];
	
	if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == mz-2 || k==1 ) && user->bctype[0]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		  p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	}
		
	if (!(nvert[k][j][i] + nvert[k][j+1][i])) {
	  ucont[k][j][i].y -=
	    (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		     jcsi[k][j][i].y * jeta[k][j][i].y +
		     jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	     dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		     jeta[k][j][i].y * jeta[k][j][i].y +
		     jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	     dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		     jzet[k][j][i].y * jeta[k][j][i].y +
		     jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i])
	    * dt * st / COEF_TIME_ACCURACY;
	}
      }
    }
  }
  //k+direction
  if (user->bctype[4]==7 && zs==0){
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	k=zs;
	dpdc = 0.;
	dpde = 0.;

	if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		  p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	}
	
	if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		  p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	}
	
	dpdz = p[k+1][j][i] - p[k][j][i];
	
	if (!(nvert[k][j][i] + nvert[k+1][j][i])) {
	  
	  ucont[k][j][i].z -=
	   (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
	  		       kcsi[k][j][i].y * kzet[k][j][i].y +
	  		       kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	  	       dpde * (keta[k][j][i].x * kzet[k][j][i].x +
	  		       keta[k][j][i].y * kzet[k][j][i].y +
	  		       keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	  	       dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
	  		       kzet[k][j][i].y * kzet[k][j][i].y +
	  		       kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i])
	  	      * dt * st / COEF_TIME_ACCURACY;
	  
	}
      }
    }
  }

  ////.................////
  //
  if (channelz){
  double lArea,AreaSum;
  double *FluxIn, *Fluxsum, *ratio;
    FluxIn=malloc((mz-1) *sizeof(double));	Fluxsum=malloc((mz-1) *sizeof(double));  ratio=malloc((mz-1) *sizeof(double));
  
    for (i=0;i<mz-1;i++) {
	 FluxIn[i]=0.0; ratio[i]=0.0; Fluxsum[i]=0.0;
	}
    lArea=0.0;
    if (zs==0) {
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k+1][j][i]<0.1){
	    
	    lArea += sqrt((zet[k+1][j][i].x) * (zet[k+1][j][i].x) +
			  (zet[k+1][j][i].y) * (zet[k+1][j][i].y) +
			  (zet[k+1][j][i].z) * (zet[k+1][j][i].z));
	  }
	}
      }
    }
    

    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); 
  
	for (k=zs; k<lze; k++){
	 for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	   if (nvert[k+1][j][i]<0.1){
	    FluxIn[k] += ucont[k][j][i].z;
	    
	  }
	}
      }
	}
 
 

     MPI_Allreduce(FluxIn,Fluxsum,mz-1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
    for (i=0;i<mz-1;i++) {
	 
	      ratio[i]=(-Fluxsum[i]+user->simCtx->FluxInSum)/AreaSum;
	      PetscPrintf(PETSC_COMM_WORLD, "Fluxsum %f   %d\n",ratio[i],mz-1); 
	    }
 
/*
     for (k=lzs; k<lze; k++) {
       for (j=lys; j<lye; j++) {
         for (i=lxs; i<lxe; i++) {
	   if (nvert[k][j][i]<0.1){
	  	   ucont[k][j][i].z+=ratio[k]* sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) + 
 					   (zet[k][j][i].y) * (zet[k][j][i].y) + 
 					   (zet[k][j][i].z) * (zet[k][j][i].z)); 
	}
      }
    } 
   }



    if (zs==0) {
	k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k+1][j][i]<0.1){
	   ucont[k][j][i].z+=ratio[k]* sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) + 
 					   (zet[k+1][j][i].y) * (zet[k+1][j][i].y) + 
 					   (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); 
	  }
	}
      }
    }  
  
    if (ze==mz) {
      k=mz-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k-1][j][i]<0.1){
 	       ucont[k-1][j][i].z+=ratio[k-1]* sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) + 
 					   (zet[k-1][j][i].y) * (zet[k-1][j][i].y) + 
 					   (zet[k-1][j][i].z) * (zet[k-1][j][i].z));  
	  
	   
					   
					   
	  }
	}
      }
    }
*/
free (ratio); free(Fluxsum);  free(FluxIn);
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (TwoD==1) */
/* 	  ucont[k][j][i].x =0.; */
/* 	else if (TwoD==2) */
/* 	  ucont[k][j][i].y =0.; */
/* 	else if (TwoD==3) */
/* 	  ucont[k][j][i].z =0.; */
/*       } */
/*     } */
/*   } */
  }

  DMDAVecRestoreArray(fda, user->Ucont, &ucont);

  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);

  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);

  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);

  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);

  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lIAj, &iaj);
  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  DMDAVecRestoreArray(da, user->lKAj, &kaj);
  DMDAVecRestoreArray(da, user->lPhi, &p);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
 
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);


 
  Contra2Cart(user);
//  FormBCS(user);
  GhostNodeVelocity(user);
  return(0);
}


//////////////////////////////////////////////////////



/**
 * @brief Initializes or updates the complete, consistent state of all Eulerian fields for a given timestep.
 *
 * This function is a high-level wrapper that orchestrates the entire process of preparing
 * the fluid fields for a single time step. It follows the standard procedure for a
 * curvilinear solver: first resolving contravariant velocities (`Ucont`) and then
 * converting them to Cartesian (`Ucat`).
 *
 * Its sequential operations are:
 * 1.  Update the INTERIOR of the domain:
 *     - For the initial step, it calls `SetInitialInteriorField` to generate values.
 *     - For subsequent steps, it calls the main fluid solver.
 *     - If restarting from a file, it reads the data, overwriting the whole field.
 *
 * 2.  Apply Boundary Conditions:
 *     - It then calls the modular `BoundarySystem_ExecuteStep` to enforce all configured
 *       boundary conditions on the domain edges.
 *
 * 3.  Convert to Cartesian and Finalize:
 *     - It calls `Contra2Cart` to compute `Ucat` from `Ucont`.
 *     - It calls `UpdateLocalGhosts` to ensure all parallel data is synchronized.
 *
 * @param user        Pointer to the UserCtx structure, containing all simulation data.
 * @param step        The current timestep number being processed.
 * @param StartStep   The initial timestep number of the simulation.
 * @param time        The current simulation time.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetEulerianFields(UserCtx *user)
{
  
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    SimCtx *simCtx = user->simCtx;
    PetscInt step = simCtx->step;
    PetscInt StartStep = simCtx->StartStep;
    PetscReal time = simCtx->ti;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preparing complete Eulerian state.\n", time, step);

    PetscReal unorm=0.0;
    // PetscReal ucont_max=0.0;
    //   PetscReal umin=0.0;

    // ==============================================================================
    // --- STEP 1: Update the INTERIOR of the domain based on the simulation phase ---
    // ==============================================================================

    if (step == StartStep && StartStep > 0) {
        // --- PATH 1: RESTART from file ---
        // This is the first time this function is called in a restarted run.
        LOG_ALLOW(GLOBAL, LOG_INFO, "RESTART condition: Reading all grid fields from file for step %d.\n", step);
        ierr = ReadSimulationFields(user, step); CHKERRQ(ierr); // Assumes this function reads Ucat, Ucont, etc.

        // After loading, we MUST update local ghosts to ensure consistency for any subsequent calculations.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for all fields after reading.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);

    } else {
        // --- PATH 2 & 3: FRESH START or TIME ADVANCEMENT ---
        // This block handles both generating initial fields and advancing the solver in time.

        if (step == 0) { // Condition is now simply step == 0 for a fresh start
            // --- PATH 2: Initial Field Setup (t=0) ---
            LOG_ALLOW(GLOBAL, LOG_INFO, "FRESH START: Generating INTERIOR fields for initial step 0.\n");
            ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);
        } else {
            // --- PATH 3: Advancing the simulation (step > 0) ---
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "TIME ADVANCE: Updating INTERIOR fields for step %d.\n", step);
            // This is the hook for the actual fluid dynamics solver.
            // ierr = YourNavierStokesSolver(user, user->dt); CHKERRQ(ierr);
	    // ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);
        }

        // The following logic is common to both fresh starts and time advancement,
        // but not to a file-based restart (which loads the final Ucat directly).

        // STEP 2: APPLY BOUNDARY CONDITIONS
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Executing LEGACY boundary condition system.\n", time, step);
        ierr = BoundarySystem_ExecuteStep_Legacy(user); CHKERRQ(ierr);

        // STEP 3: SYNCHRONIZE Ucont BEFORE CONVERSION
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucont.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

        // STEP 4: CONVERT CONTRAVARIANT TO CARTESIAN
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Converting Ucont to Ucat.\n", time, step);
        ierr = Contra2Cart(user); CHKERRQ(ierr);

        // STEP 5: Re-apply BCs and SYNCHRONIZE Ucat
        // It's often necessary to apply BCs again to ensure Ucat is correct at boundaries.
	//   ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucat.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    }

    ierr = VecNorm(user->Ucat,NORM_INFINITY,&unorm);
    //unorm = unorm/(3*(user->IM*user->JM*user->KM));
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Complete Eulerian state is now finalized and consistent. Max Ucat = %.6f \n", time, step,unorm);
    
    PetscFunctionReturn(0);
}


/////////////////////////////////////////////////////////


/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator,
 * which handles both location and migration in a single, robust, iterative process.
 *
 * Its sequential operations are:
 * 1. A single, comprehensive call to `LocateAllParticlesInGrid_TEST` to sort all particles
 *    to their correct owner ranks and find their initial host cells.
 * 2. If `user->ParticleInitialization == 0` (Surface Init), it re-initializes particles on the
 *    designated inlet surface, now that they are on the correct MPI ranks.
 * 3. A second call to `LocateAllParticlesInGrid_TEST` is needed after re-initialization to
 *    find the new, correct host cells for the surface-placed particles.
 * 4. Interpolates initial Eulerian fields to the settled particles.
 * 5. Scatters particle data to Eulerian fields (if applicable).
 * 6. Outputs initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscReal currentTime;
    PetscInt step;
    PetscInt OutputFreq;
    PetscInt StepsToRun;
    PetscInt StartStep;
    
    PetscFunctionBeginUser;

    SimCtx *simCtx = user->simCtx;

    currentTime = simCtx->ti;
    step = simCtx->step;
    OutputFreq = simCtx->OutputFreq;
    StepsToRun = simCtx->StepsToRun;
    StartStep = simCtx->StartStep;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures [TEST].\n", currentTime, step);

    // --- 1. Initial Particle Settlement (Location and Migration) ---
    // This single call replaces the old sequence of Locate -> Migrate. The new
    // orchestrator handles the iterative process internally until all particles are settled.
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Initial Settlement: Locating and migrating all particles to their correct ranks and cells.\n", currentTime, step);
    ierr = LocateAllParticlesInGrid_TEST(user,bboxlist); CHKERRQ(ierr);

    // --- 2. Re-initialize Particles on Inlet Surface (if applicable) ---
    if (simCtx->ParticleInitialization == 0 && user->inletFaceDefined) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Re-initializing particles on inlet surface now that they are on correct ranks.\n", currentTime, step);
        ierr = ReinitializeParticlesOnInletSurface(user, currentTime, step); CHKERRQ(ierr);

	LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Resetting statuses for post-reinitialization settlement.\n", currentTime, step);
        ierr = ResetAllParticleStatuses(user); CHKERRQ(ierr);
	
        // --- CRITICAL: After re-placing particles, we MUST locate them again. ---
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Post-Reinitialization Settlement: Finding host cells for newly placed inlet particles.\n", currentTime, step);
        ierr = LocateAllParticlesInGrid_TEST(user,bboxlist); CHKERRQ(ierr);
    }
    
    // --- 3. Finalize State for t=0 ---
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Interpolating initial fields to settled particles.\n", currentTime, step);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);
    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- 4. Initial Output ---
    if (OutputFreq > 0 || (StepsToRun == 0 && StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation data.\n", currentTime, step);
        ierr = LOG_PARTICLE_FIELDS(user, simCtx->LoggingFrequency); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "position", step, "dat"); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "velocity", step, "dat"); CHKERRQ(ierr);
	//  if (!readFields) {
        ierr = WriteSimulationFields(user); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Initializes or updates the complete, consistent state of all Eulerian fields for a given timestep.
 *
 * This function is a high-level wrapper that orchestrates the entire process of preparing
 * the fluid fields for a single time step. It follows the standard procedure for a
 * curvilinear solver: first resolving contravariant velocities (`Ucont`) and then
 * converting them to Cartesian (`Ucat`).
 *
 * Its sequential operations are:
 * 1.  Update the INTERIOR of the domain:
 *     - For the initial step, it calls `SetInitialInteriorField` to generate values.
 *     - For subsequent steps, it calls the main fluid solver.
 *     - If restarting from a file, it reads the data, overwriting the whole field.
 *
 * 2.  Apply Boundary Conditions:
 *     - It then calls the modular `BoundarySystem_ExecuteStep` to enforce all configured
 *       boundary conditions on the domain edges.
 *
 * 3.  Convert to Cartesian and Finalize:
 *     - It calls `Contra2Cart` to compute `Ucat` from `Ucont`.
 *     - It calls `UpdateLocalGhosts` to ensure all parallel data is synchronized.
 *
 * @param user        Pointer to the UserCtx structure, containing all simulation data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetEulerianFields(UserCtx *user);


/////////////////////////////////////

/**
 * @brief Gathers a PETSc vector onto rank 0 as a contiguous array of doubles.
 *
 * This function retrieves the local portions of the input vector \p inVec from
 * all MPI ranks via \c MPI_Gatherv and assembles them into a single array on rank 0.
 * The global size of the vector is stored in \p N, and a pointer to the newly
 * allocated array is returned in \p arrayOut (valid only on rank 0).
 *
 * @param[in]  inVec      The PETSc vector to gather.
 * @param[out] N          The global size of the vector (output).
 * @param[out] arrayOut   On rank 0, points to the newly allocated array of size \p N.
 *                        On other ranks, it is set to NULL.
 *
 * @return PetscErrorCode  Returns 0 on success, or a non-zero PETSc error code.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut)
{
    PetscErrorCode    ierr;
    MPI_Comm          comm;
    PetscMPIInt       rank, size;
    PetscInt          globalSize, localSize;
    const PetscScalar *localArr = NULL;

    PetscFunctionBeginUser;
    
    // Log entry into the function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "VecToArrayOnRank0 - Start gathering vector onto rank 0.\n");

    /* Get MPI comm, rank, size */
    ierr = PetscObjectGetComm((PetscObject)inVec, &comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);

    /* Get global size (for the entire Vec) */
    ierr = VecGetSize(inVec, &globalSize);CHKERRQ(ierr);
    *N = globalSize;

    /* Get local size (portion on this rank) */
    ierr = VecGetLocalSize(inVec, &localSize);CHKERRQ(ierr);

    // Log vector sizes and process info
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "VecToArrayOnRank0 - rank=%d of %d, globalSize=%d, localSize=%d.\n",
              rank, size, globalSize, localSize);

    /* Access the local array data */
    ierr = VecGetArrayRead(inVec, &localArr);CHKERRQ(ierr);

    /*
       We'll gather the local chunks via MPI_Gatherv:
       - First, gather all local sizes into recvcounts[] on rank 0.
       - Then set up a displacement array (displs[]) to place each chunk in the correct spot.
       - Finally, gather the actual data.
    */
    PetscMPIInt *recvcounts = NULL;
    PetscMPIInt *displs     = NULL;

    if (!rank) {
        ierr = PetscMalloc2(size, &recvcounts, size, &displs);CHKERRQ(ierr);
    }

    /* Convert localSize (PetscInt) to PetscMPIInt for MPI calls */
    PetscMPIInt localSizeMPI = (PetscMPIInt)localSize;

    /* Gather local sizes to rank 0 */
    ierr = MPI_Gather(&localSizeMPI, 1, MPI_INT,
                      recvcounts,      1, MPI_INT,
                      0, comm);CHKERRQ(ierr);

    /* On rank 0, build displacements and allocate the big array */
    if (!rank) {
        displs[0] = 0;
        for (PetscMPIInt i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        /* Allocate a buffer for the entire (global) array */
        ierr = PetscMalloc1(globalSize, arrayOut);CHKERRQ(ierr);
        if (!(*arrayOut)) SETERRQ(comm, PETSC_ERR_MEM, "Failed to allocate array on rank 0.");
    } else {
        /* On other ranks, we do not allocate anything */
        *arrayOut = NULL;
    }

    /* Gather the actual data (assuming real scalars => MPI_DOUBLE) */
    ierr = MPI_Gatherv((void *) localArr,         /* sendbuf on this rank */
                       localSizeMPI, MPI_DOUBLE,  /* how many, and type */
                       (rank == 0 ? *arrayOut : NULL),  /* recvbuf on rank 0 */
                       (rank == 0 ? recvcounts : NULL),
                       (rank == 0 ? displs    : NULL),
                       MPI_DOUBLE, 0, comm);CHKERRQ(ierr);

    /* Restore local array (cleanup) */
    ierr = VecRestoreArrayRead(inVec, &localArr);CHKERRQ(ierr);

    if (!rank) {
        ierr = PetscFree(recvcounts);
        ierr = PetscFree(displs);
    }

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_INFO, "VecToArrayOnRank0 - Successfully gathered data on rank 0.\n");

    PetscFunctionReturn(0);
}

////////////////////////////////////////


/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * This function first computes the local bounding box on each process by calling
 * `ComputeLocalBoundingBox`. It then uses an MPI gather operation to collect all
 * local bounding boxes on the root process (rank 0). On rank 0, it allocates an array
 * of `BoundingBox` structures to hold the gathered data and returns it via the
 * `allBBoxes` Pointer. On other ranks, `allBBoxes` is set to `NULL`.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 *                        This context must be properly initialized before calling this function.
 * @param[out] allBBoxes  Pointer to a Pointer where the array of gathered bounding boxes will be stored on rank 0.
 *                        On rank 0, this will point to the allocated array; on other ranks, it will be `NULL`.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure.
 */
/*
PetscErrorCode GatherAllBoundingBoxes(UserCtx *user, BoundingBox **allBBoxes)
{
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    BoundingBox *bboxArray = NULL;
    BoundingBox localBBox;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Entering the function. \n");

    // Validate input Pointers
    if (!user) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Input 'user' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!allBBoxes) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Output 'allBBoxes' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }

    // Get the rank and size of the MPI communicator
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (ierr != MPI_SUCCESS) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error getting MPI rank.\n");
        return ierr;
    }
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (ierr != MPI_SUCCESS) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error getting MPI size.\n");
        return ierr;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "MPI rank=%d, size=%d.\n", rank, size);

    // Compute the local bounding box on each process
    ierr = ComputeLocalBoundingBox(user, &localBBox);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error computing local bounding box.\n");
        return ierr;
    }
    
    PetscBarrier(PETSC_NULLPTR);

    // On rank 0, allocate memory for the array of bounding boxes
    if (rank == 0) {
        bboxArray = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!bboxArray) {
            LOG_ALLOW(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Memory allocation failed for bounding box array.\n");
            return PETSC_ERR_MEM;
        }
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAllBoundingBoxes: Allocated memory for bounding box array on rank 0.\n");
    }

    // Perform MPI_Gather to collect all local bounding boxes on rank 0
    // Corrected MPI_Gather call
ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
                  (rank == 0) ? bboxArray : NULL,  // Explicitly NULL on non-roots
                  sizeof(BoundingBox), MPI_BYTE,   // Recv count is ignored on non-roots
                  0, PETSC_COMM_WORLD); CHKERRMPI(ierr);
 
//   ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
// bboxArray, sizeof(BoundingBox), MPI_BYTE,
//                      0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Error during MPI_Gather operation.\n");
        if (rank == 0 && bboxArray) free(bboxArray); // Clean up if allocation was done
        return ierr;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Successfully gathered bounding boxes on rank 0.\n",rank);    
    
    // On rank 0, assign the gathered bounding boxes to the output Pointer
    if (rank == 0) {
        *allBBoxes = bboxArray;
    } else {
        *allBBoxes = NULL;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Exiting the function successfully.\n");
    return 0;
}
*/
/**
 * @brief Broadcasts the bounding box information collected on rank 0 to all other ranks.
 *
 * This function assumes that `GatherAllBoundingBoxes()` was previously called, so `bboxlist`
 * is allocated and populated on rank 0. All other ranks will allocate memory for `bboxlist`,
 * and this function will use MPI_Bcast to distribute the bounding box data to them.
 *
 * @param[in]     user      Pointer to the UserCtx structure. (Currently unused in this function, but kept for consistency.)
 * @param[in,out] bboxlist  Pointer to the array of BoundingBoxes. On rank 0, this should point to
 *                          a valid array of size 'size' (where size is the number of MPI ranks).
 *                          On non-root ranks, this function will allocate memory for `bboxlist`.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on MPI or PETSc-related errors.
 */
/*
PetscErrorCode BroadcastAllBoundingBoxes(UserCtx *user, BoundingBox **bboxlist) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;

    // Get MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // On non-root ranks, allocate memory for bboxlist before receiving the broadcast
    if (rank != 0) {
        *bboxlist = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!*bboxlist) SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_MEM, "Failed to allocate memory for bboxlist on non-root ranks.");
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Broadcasting bounding box information from rank 0.\n");

    // Broadcast bboxlist from rank 0 to all other ranks
    ierr = MPI_Bcast(*bboxlist, (PetscInt)(size * sizeof(BoundingBox)), MPI_BYTE, 0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_LIB, "MPI_Bcast failed for bboxlist.");
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Broadcasted bounding box information from rank 0.\n");    

    return 0;
}
*/

///////////////////////////


/**
 * @brief Writes data from a specific PETSc vector to a file.
 *
 * This function uses the field name to construct the file path and writes the data
 * from the provided PETSc vector to the corresponding file.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context.
 * @param[in] field_name Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[in] field_vec  PETSc vector containing the field data to write.
 * @param[in] ti         Time index for constructing the file name.
 * @param[in] ext        File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
/*
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    char filen[128];

    // Construct the file name
    snprintf(filen, sizeof(filen), "results/%s%05d_%d.%s", field_name, ti, user->_this, ext);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteFieldData - Attempting to write file: %s\n", filen);

    // Open the file for writing
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

    PetscReal vmin, vmax;
    ierr = VecMin(field_vec, NULL, &vmin); CHKERRQ(ierr);
    ierr = VecMax(field_vec, NULL, &vmax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"%s step %d  min=%.6e  max=%.6e\n",field_name, ti, (double)vmin, (double)vmax);

    
    // Write data from the vector
    ierr = VecView(field_vec, viewer); CHKERRQ(ierr);

    // Close the file viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);

    return 0;
}
*/

/*
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    char           filename[PETSC_MAX_PATH_LEN];
    PetscMPIInt    rank, size;
    MPI_Comm       comm;
    PetscInt       placeholder_int = 0;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)field_vec, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

    // Construct a filename that does NOT depend on the rank number.
    // This ensures a single, predictable output file.
    ierr = PetscSNPrintf(filename, sizeof(filename), "results/%s%05"PetscInt_FMT"_%d.%s", field_name, ti, placeholder_int, ext);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteFieldData - Preparing to write file: %s\n", filename);

    // Optional: Log min/max values. This is a collective operation.
    PetscReal vmin, vmax;
    ierr = VecMin(field_vec, NULL, &vmin); CHKERRQ(ierr);
    ierr = VecMax(field_vec, NULL, &vmax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"%s step %d  min=%.6e  max=%.6e\n", field_name, ti, (double)vmin, (double)vmax);

    if (size == 1) {
        // --- SERIAL CASE: Simple and direct write ---
        PetscViewer viewer;
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
        ierr = VecView(field_vec, viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    } else {
        // --- PARALLEL CASE: Gather on Rank 0 and write sequentially ---
        Vec         seq_vec = NULL; // A sequential vector on Rank 0
        VecScatter  scatter_ctx;

        // Create a scatter context to gather data from the parallel vector
        // onto a new sequential vector that will exist only on Rank 0.
        ierr = VecScatterCreateToZero(field_vec, &scatter_ctx, &seq_vec); CHKERRQ(ierr);

        // Perform the scatter (gather) operation.
        ierr = VecScatterBegin(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&scatter_ctx); CHKERRQ(ierr);

        // Now, only Rank 0 has the populated sequential vector and can write it.
        if (rank == 0) {
            PetscViewer viewer;
            ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
            ierr = VecView(seq_vec, viewer); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        }

        // All ranks must destroy the sequential vector, though it only "exists" on Rank 0.
        ierr = VecDestroy(&seq_vec); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);
    PetscFunctionReturn(0);
}
*/
/*
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    MPI_Comm       comm;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)field_vec, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

    if (size == 1) {
        // SERIAL CASE: This is simple and correct.
        PetscViewer viewer;
        char filename[PETSC_MAX_PATH_LEN];
        ierr = PetscSNPrintf(filename, sizeof(filename), "output/%s_step%d.%s", field_name, ti, ext); CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
        ierr = VecView(field_vec, viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    } else {
        // PARALLEL CASE: The robust "Gather on 0 and Write" pattern.
        Vec         seq_vec = NULL;
        VecScatter  scatter_ctx;
        IS          is_from, is_to;
        PetscInt    global_size;

        ierr = VecGetSize(field_vec, &global_size); CHKERRQ(ierr);
        ierr = ISCreateStride(comm, global_size, 0, 1, &is_from); CHKERRQ(ierr);
        ierr = ISCreateStride(PETSC_COMM_SELF, global_size, 0, 1, &is_to); CHKERRQ(ierr);

        // Rank 0 creates the destination sequential vector.
        if (rank == 0) {
            ierr = VecCreate(PETSC_COMM_SELF, &seq_vec); CHKERRQ(ierr);
            ierr = VecSetSizes(seq_vec, global_size, global_size); CHKERRQ(ierr);
            ierr = VecSetFromOptions(seq_vec); CHKERRQ(ierr);
        }

        // Create the general scatter context.
        ierr = VecScatterCreate(field_vec, is_from, seq_vec, is_to, &scatter_ctx); CHKERRQ(ierr);
        ierr = VecScatterBegin(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

        // Only Rank 0 writes the populated sequential vector.
        if (rank == 0) {
            PetscViewer viewer;
            char filename[PETSC_MAX_PATH_LEN];
            ierr = PetscSNPrintf(filename, sizeof(filename), "output/%s_step%d.%s", field_name, ti, ext); CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
            ierr = VecView(seq_vec, viewer); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        }

        // Clean up.
        ierr = ISDestroy(&is_from); CHKERRQ(ierr);
        ierr = ISDestroy(&is_to); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&scatter_ctx); CHKERRQ(ierr);
        if (rank == 0) {
            ierr = VecDestroy(&seq_vec); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);
    PetscFunctionReturn(0);
}
*/
/////////////////

/**
 * @brief Gathers a distributed DMSwarm field into a single C array on rank 0.
 *
 * This is a high-performance helper specifically for post-processing I/O. It
 * directly accesses local swarm data and uses MPI_Gatherv to collect it on rank 0,
 * avoiding the overhead of creating an intermediate PETSc Vec object.
 *
 * @param[in]  swarm             The DMSwarm object.
 * @param[in]  field_name        The name of the field to gather (e.g., "velocity").
 * @param[out] out_n_global      On rank 0, contains the total number of particles. 0 on other ranks.
 * @param[out] out_n_components  On rank 0, contains the number of components for the field. 0 on other ranks.
 * @param[out] out_data_rank0    On rank 0, a pointer to a newly allocated array with the gathered data.
 *                               The caller is responsible for freeing this memory.
 * @return PetscErrorCode
 */
PetscErrorCode SwarmFieldToArrayOnRank0(DM swarm, const char* field_name, PetscInt *out_n_global, PetscInt *out_n_components, PetscScalar **out_data_rank0)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    PetscInt       n_local, n_global, n_comp;
    PetscScalar    *local_data;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // Access the raw local data pointer and field info
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &n_global); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, &n_comp, NULL, (void**)&local_data); CHKERRQ(ierr);

    // --- MPI_Gatherv Implementation ---
    PetscInt *recvcounts = NULL, *displs = NULL;
    if (rank == 0) {
        ierr = PetscMalloc1(size, &recvcounts); CHKERRQ(ierr);
        ierr = PetscMalloc1(size, &displs); CHKERRQ(ierr);
        ierr = PetscMalloc1(n_global * n_comp, out_data_rank0); CHKERRQ(ierr);
    }
    PetscInt sendcount = n_local * n_comp;
    ierr = MPI_Gather(&sendcount, 1, MPIU_INT, recvcounts, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }
    }
    ierr = MPI_Gatherv(local_data, sendcount, MPIU_SCALAR,
                       (rank == 0) ? *out_data_rank0 : NULL, recvcounts, displs, MPIU_SCALAR,
                       0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    // --- End MPI_Gatherv ---

    ierr = DMSwarmRestoreField(swarm, field_name, &n_comp, NULL, (void**)&local_data); CHKERRQ(ierr);
    
    if (rank == 0) {
        *out_n_global = n_global;
        *out_n_components = n_comp;
        ierr = PetscFree(recvcounts); CHKERRQ(ierr);
        ierr = PetscFree(displs); CHKERRQ(ierr);
    } else {
        *out_n_global = 0; *out_n_components = 0; *out_data_rank0 = NULL;
    }
    PetscFunctionReturn(0);
}