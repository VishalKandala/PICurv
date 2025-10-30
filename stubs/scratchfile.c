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

///////////////////////////////////////

/**
 * @brief Identifies particles that have left the local MPI rank's domain and determines their target rank.
 *
 * This function iterates through all particles currently local to this MPI rank.
 * For each particle, it first checks if the particle is still within the rank's
 * primary bounding box (defined by `user->bbox`).
 *
 * If a particle is outside `user->bbox`:
 * 1. It preferentially checks the 6 immediate Cartesian neighbors (defined in `user->neighbors`).
 *    The bounding boxes for these neighbors are looked up in the global `bboxlist`.
 * 2. If the particle is not found in any of the immediate Cartesian neighbors, a fallback
 *    search is initiated. This fallback search iterates through *all other* MPI ranks
 *    (excluding the current rank and already checked immediate neighbors if optimized)
 *    using the `bboxlist` to find a rank whose bounding box contains the particle.
 *
 * Particles successfully assigned a `targetRank` are added to the `migrationList`.
 * If a particle leaves the local box but is not found in any other rank's bounding box
 * (even after the fallback), a warning is logged, as this particle might be lost from
 * the global computational domain (assuming `bboxlist` covers the entire domain).
 * Such "lost" particles should ideally be handled by a separate global boundary condition
 * check (e.g., `CheckAndRemoveOutOfBoundsParticles`) prior to calling this function.
 *
 * The `migrationList` is dynamically reallocated if its current capacity is exceeded.
 *
 * @param[in]      user           Pointer to the UserCtx structure. It must contain:
 *                                - `swarm`: The DMSwarm object.
 *                                - `bbox`: The BoundingBox of the current MPI rank.
 *                                - `neighbors`: A RankNeighbors struct with the ranks of the 6 Cartesian neighbors.
 *                                - `size`: The total number of MPI ranks in PETSC_COMM_WORLD (implicitly via MPI_Comm_size).
 * @param[in]      bboxlist       An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                This array must be up-to-date and available on all ranks.
 * @param[in,out]  migrationList  Pointer to an array of MigrationInfo structures. This array will be
 *                                populated with particles to be migrated. The function may reallocate
 *                                this array if more space is needed. The caller is responsible for
 *                                eventually freeing this memory if it's not NULL.
 * @param[out]     migrationCount Pointer to a PetscInt that will be set to the number of particles
 *                                identified for migration (i.e., the number of valid entries in `*migrationList`).
 * @param[in,out]  listCapacity   Pointer to a PetscInt representing the current allocated capacity of
 *                                `*migrationList` (in terms of number of MigrationInfo structs).
 *                                This function will update it if `*migrationList` is reallocated.
 *
 * @return PetscErrorCode 0 on success, non-zero on PETSc/MPI errors or if essential fields are missing.
 *
 * @note It is assumed that `user->neighbors` contains valid rank identifiers or `MPI_PROC_NULL`.
 * @note It is assumed that `bboxlist` is correctly populated and broadcast to all ranks before calling this.
 * @note For the fallback search to be effective, `bboxlist` should represent a tiling of the
 *       global domain. Particles not found in any box in `bboxlist` are effectively outside this
 *       tiled global domain.
 */
PetscErrorCode IdentifyMigratingParticles(UserCtx *user,
					  const BoundingBox *bboxlist,
					  MigrationInfo **migrationList,
					  PetscInt *migrationCount,
					  PetscInt *listCapacity)
{
  PetscErrorCode ierr;
  DM             swarm = user->swarm;
  PetscInt       nLocal, p;
  Cmpnts        *pos = NULL;
  PetscMPIInt    rank,size;
  BoundingBox    localBBox = user->bbox;
  RankNeighbors  neighbors = user->neighbors; // Use stored neighbors
  
  //PetscInt       currentMigrationCount = 0;
  //PetscInt       currentListCapacity = *listCapacity;
  //MigrationInfo *currentMigrationList = *migrationList;
  // Add PID pointer if logging PIDs
  PetscInt64    *pids = NULL;
  PetscFunctionBeginUser;

  // --- Input Validation and Initialization ---
  if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx 'user' is NULL.");
  if (!bboxlist) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Global 'bboxlist' is NULL.");
  if (!migrationList) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'migrationList' output pointer is NULL.");
  if (!migrationCount) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'migrationCount' output pointer is NULL.");
  if (!listCapacity) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'listCapacity' output pointer is NULL.");

  if (!swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx->swarm is NULL.");
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  if (nLocal == 0) {
    *migrationCount = 0;
    // Ensure output pointers are consistent even if no allocation happened
    // *migrationList = currentMigrationList;
    // *listCapacity = currentListCapacity;
    PetscFunctionReturn(0);
  }
  // Get read-only access to position
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr); // If logging PIDs
  if (!pos || !pids) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not access required DMSwarm fields.");
  /*
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "[IdentifyMigratingParticles - Rank %d INCOMING user->neighbors] xm=%d, xp=%d, ym=%d, yp=%d, zm=%d, zp=%d. MPI_PROC_NULL is %d. \n",
    rank, user->neighbors.rank_xm, user->neighbors.rank_xp,
    user->neighbors.rank_ym, user->neighbors.rank_yp,
    user->neighbors.rank_zm, user->neighbors.rank_zp, (int)MPI_PROC_NULL);
  */
  
  *migrationCount = 0;
  //currentMigrationCount = 0; // Reset count for this call
  for (p = 0; p < nLocal; p++) {
    // Check if particle is OUTSIDE the local bounding box
    if (!IsParticleInBox(&localBBox, &pos[p]))
      {
	PetscInt targetRank = MPI_PROC_NULL; // Target rank not yet found

	// Determine likely exit direction(s) to prioritize neighbor check
	PetscBool exit_xm = pos[p].x < localBBox.min_coords.x;
	PetscBool exit_xp = pos[p].x > localBBox.max_coords.x;
	PetscBool exit_ym = pos[p].y < localBBox.min_coords.y;
	PetscBool exit_yp = pos[p].y > localBBox.max_coords.y;
	PetscBool exit_zm = pos[p].z < localBBox.min_coords.z;
	PetscBool exit_zp = pos[p].z > localBBox.max_coords.z;

	// DEBUG ------------
	/*
	  if (rank == 1 && p < 5) { // Log for first few particles on Rank 1
	  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
	  "[Identify - Rank 1, p=%d] Particle Pos(%.2f,%.2f,%.2f). Exits(xm%d,xp%d,ym%d,yp%d,zm%d,zp%d). Neighbors(xm%d,xp%d,ym%d,yp%d,zm%d,zp%d)",
	  p, pos[p].x, pos[p].y, pos[p].z,
	  exit_xm, exit_xp, exit_ym, exit_yp, exit_zm, exit_zp,
	  neighbors.rank_xm, neighbors.rank_xp, neighbors.rank_ym, neighbors.rank_yp, neighbors.rank_zm, neighbors.rank_zp);
	  }
	*/
	// DEBUG -----------
        
	
	LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d[PID = %ld] at (%g,%g,%g) left local bbox. Checking neighbors...\n",rank, p,pids[p], pos[p].x, pos[p].y, pos[p].z);

	//1.  Check neighbors preferentially
	if (exit_xm && neighbors.rank_xm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xm], &pos[p])) targetRank = neighbors.rank_xm;
	else if (exit_xp && neighbors.rank_xp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xp], &pos[p])) targetRank = neighbors.rank_xp;
	else if (exit_ym && neighbors.rank_ym != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_ym], &pos[p])) targetRank = neighbors.rank_ym;
	else if (exit_yp && neighbors.rank_yp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_yp], &pos[p])) targetRank = neighbors.rank_yp;
	else if (exit_zm && neighbors.rank_zm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zm], &pos[p])) targetRank = neighbors.rank_zm;
	else if (exit_zp && neighbors.rank_zp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zp], &pos[p])) targetRank = neighbors.rank_zp;
	// Add checks for edge/corner neighbors if needed and if they were stored

	
	// 2.--- Fallback (if strict Cartesian neighbors aren't enough) ---
	
	if (targetRank == MPI_PROC_NULL){
	  /* ... loop through all ranks in bboxlist ... */
        
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %" PetscInt64_FMT " not in immediate neighbors. Fallback: checking all %d other ranks...",rank, pids[p],size);

	  for(PetscMPIInt r = 0; r < size; r++){

	    if(IsParticleInBox(&bboxlist[r], &pos[p])){

	      targetRank = r;

	      LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Particle %ld FOUND in bboxlist[%d] during fallback search. \n",rank, pids[p],r);
	      break;
	    }
	    
	  }
	  
	}
	// 3. ---- Add to migration list if a target rank was found.
	if (targetRank != MPI_PROC_NULL){

	  /* OLD LOGIC
          // Resize list if needed (using PetscRealloc for safety)
	  if (currentMigrationCount >= currentListCapacity) {
	    PetscInt OldCapacity = currentListCapacity;
	    
	      PetscInt newCapacity = (currentListCapacity == 0) ? 16 : currentListCapacity * 2;

	    ierr = PetscRealloc((size_t)newCapacity * sizeof(MigrationInfo),
				migrationList); CHKERRQ(ierr);

	    currentMigrationList = *migrationList;
	    
	    
	    currentListCapacity = newCapacity;

	    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reallocated migrationList capacity from %d to %d\n", rank,OldCapacity, newCapacity);
	   }
	    
	  // Add to migration list
	  currentMigrationList[currentMigrationCount].local_index = p;
	  currentMigrationList[currentMigrationCount].target_rank = targetRank;
	  // currentMigrationList[currentMigrationCount].pid = pids[p]; // If storing PID
	  currentMigrationCount++;
	  */

	  // NEW LOGIC (TEST)
          
	  ierr = AddToMigrationList(migrationList,    // Pointer to the list pointer
				    listCapacity,     // Pointer to the capacity variable
				    migrationCount,   // Pointer to the count variable
				    p,                // The particle's local index
				    targetRank);      // The destination rank
	  CHKERRQ(ierr);
	  
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d marked for migration to rank %d.\n", rank, p, targetRank);
	  
	}
	  else {
	  // Particle left local box but was not found in any *checked* neighbor box.
	  // Since CheckAndRemove should have run first, this might indicate a particle
	  // moved more than one cell width into a diagonal neighbor's domain not checked here,
	  // or there's an issue with BBox overlap/gaps.
          LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Particle %d[PID = %ld] at (%g,%g,%g) left local bbox but target neighbor rank not found! (May be lost or need wider neighbor check).\n",
                    rank, p, pids[p],pos[p].x, pos[p].y, pos[p].z);
          // Consider marking for removal here if this case should not happen:
          // Maybe add to a separate 'removalList' or use the marking technique
          // from the CheckAndRemove function if it's readily available.
	  }
      } // end if isOutsideLocal
  } // end particle loop
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr);
  //  *migrationList = currentMigrationList;
  // *migrationCount = currentMigrationCount;
  //  *listCapacity = currentListCapacity;
  
  //LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, currentMigrationCount);
  LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, *migrationCount);
  PetscFunctionReturn(0);
}


/**
 * @brief Performs one full cycle of particle migration: identify, set ranks, and migrate.
 *
 * This function encapsulates the three main steps of migrating particles between MPI ranks:
 * 1. Identify particles on the local rank that need to move based on their current
 *    positions and the domain decomposition (`bboxlist`).
 * 2. Determine the destination rank for each migrating particle.
 * 3. Perform the actual migration using PETSc's `DMSwarmMigrate`.
 * It also calculates and logs the global number of particles migrated.
 *
 * @param user Pointer to the UserCtx structure.
 * @param bboxlist Array of BoundingBox structures defining the spatial domain of each MPI rank.
 * @param migrationList_p Pointer to a pointer for the MigrationInfo array. This array will be
 *                        allocated/reallocated by `IdentifyMigratingParticles` if necessary.
 *                        The caller is responsible for freeing this list eventually.
 * @param migrationCount_p Pointer to store the number of particles identified for migration
 *                         on the local rank. This is reset to 0 after migration for the current cycle.
 * @param migrationListCapacity_p Pointer to store the current capacity of the `migrationList_p` array.
 * @param currentTime Current simulation time (used for logging).
 * @param step Current simulation step number (used for logging).
 * @param migrationCycleName A descriptive name for this migration cycle (e.g., "Preliminary Sort", "Main Loop")
 *                           for logging purposes.
 * @param[out] globalMigrationCount_out Pointer to store the total number of particles migrated
 *                                      across all MPI ranks during this cycle.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformSingleParticleMigrationCycle(UserCtx *user, const BoundingBox *bboxlist,
                                                   MigrationInfo **migrationList_p, PetscInt *migrationCount_p,
                                                   PetscInt *migrationListCapacity_p,
                                                   PetscReal currentTime, PetscInt step, const char *migrationCycleName,
                                                   PetscInt *globalMigrationCount_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank; // MPI rank of the current process

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

        /******************************************************************/
    /*                 START OF NEW TEST LOGIC                        */
    /******************************************************************/

    // --- TEST STEP 1: Take a snapshot of PIDs BEFORE migration ---
    PetscInt nlocal_before;
    ierr = DMSwarmGetLocalSize(user->swarm, &nlocal_before); CHKERRQ(ierr);
    
    PetscInt64 *pids_before_snapshot = NULL;
    // We need to get the PID field to pass to our snapshot function
    PetscInt64 *pid_field_for_snapshot;
    ierr = DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void**)&pid_field_for_snapshot); CHKERRQ(ierr);
    
    // Call our helper to create the sorted snapshot
    ierr = GetLocalPIDSnapshot(pid_field_for_snapshot, nlocal_before, &pids_before_snapshot); CHKERRQ(ierr);

    // Restore the field immediately
    ierr = DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void**)&pid_field_for_snapshot); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_INFO, "[TEST HARNESS - Rank %d] Created pre-migration PID snapshot with %d entries.\n", rank, nlocal_before);

    /******************************************************************/
    /*                 END OF NEW TEST LOGIC (PART 1)                 */
    /******************************************************************/

    // Step 1: Identify particles that need to migrate from the current rank
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] %s Migration: Identifying migrating particles.\n", currentTime, step, migrationCycleName);
    ierr = IdentifyMigratingParticles(user, bboxlist, migrationList_p, migrationCount_p, migrationListCapacity_p); CHKERRQ(ierr);

    // Ensure Identification is done for all ranks before sharing.
    ierr = PetscBarrier(NULL);
    
    // Step 2: Get the global count of migrating particles
    PetscInt localMigrationCount = *migrationCount_p; // Use a local variable for MPI_Allreduce
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d] Before MPI_Allreduce. localMigrationCount = %d.\n", rank, localMigrationCount);
    ierr = MPI_Allreduce(&localMigrationCount, globalMigrationCount_out, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d] After MPI_Allreduce. globalMigrationCount_out = %d.\n", rank, *globalMigrationCount_out);
     
    // Step 3: Perform migration if any particles are moving globally
    if (*globalMigrationCount_out > 0) {
        LOG_ALLOW(LOCAL, LOG_INFO, "[T=%.4f, Step=%d] %s Migration: Performing migration (%d particles globally, %d locally from rank %d).\n",
                  currentTime, step, migrationCycleName, *globalMigrationCount_out, localMigrationCount, rank);
        // Set the destination ranks for the locally identified migrating particles
        ierr = SetMigrationRanks(user, *migrationList_p, localMigrationCount); CHKERRQ(ierr);
        // Execute the migration
        ierr = PerformMigration(user); CHKERRQ(ierr);

	// Make sure all ranks finish Migration before proceding.
	ierr = PetscBarrier(NULL);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] %s Migration: No particles identified globally for migration.\n", currentTime, step, migrationCycleName);
    }

        /******************************************************************/
    /*                 START OF NEW TEST LOGIC (PART 2)               */
    /******************************************************************/
    
    // --- TEST STEP 2: Call the function we want to test ---
    // This happens AFTER the migration is complete.
    if (*globalMigrationCount_out > 0) {
        LOG_ALLOW(LOCAL, LOG_INFO, "[TEST HARNESS - Rank %d] Calling FlagNewcomersForLocation to verify identification.\n", rank);
        ierr = FlagNewcomersForLocation(user->swarm, nlocal_before, pids_before_snapshot); CHKERRQ(ierr);
    }

    // --- TEST STEP 3: Cleanup ---
    ierr = PetscFree(pids_before_snapshot); CHKERRQ(ierr);

    /******************************************************************/
    /*                 END OF NEW TEST LOGIC (PART 2)                 */
    /******************************************************************/
    
    // Reset local migration count for the next potential migration cycle by the caller.
    // The migrationList and its capacity persist and are managed by the caller.
    *migrationCount_p = 0;

    PetscFunctionReturn(0);
}


////////////////////////////////////////


/*
PetscErrorCode InflowFlux(UserCtx *user) 
{
  PetscInt     i, j, k, rank, fn;
  PetscReal    FluxIn, r, uin, uin_max, xc, yc, zc, d, H=4.1, Umax=1.5;
  PetscReal    lAreaIn, AreaSumIn;
  Vec          Coor,RFC;
  Cmpnts       ***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet, ***cent;  
  
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert, ***RR;

  // Get context from the user struct, similar to the original code
  SimCtx      *simCtx = user->simCtx;
  PetscInt    inletprofile = simCtx->inletprofile;
  PetscReal   CMx_c = simCtx->CMx_c;
  PetscReal   CMy_c = simCtx->CMy_c;
  PetscReal   CMz_c = simCtx->CMz_c;
  
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // This temporary vector is used for profiles needing pre-calculated radius
  VecDuplicate(user->lNvert, &RFC);
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMGetCoordinatesLocal(da, &Coor); 
  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, RFC, &RR);
  DMDAVecGetArray(fda, user->Cent, &cent);
  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  FluxIn = 0.0;
  lAreaIn = 0.0;
  uin = 0.0;


  // ====================== DEBUGGING BLOCK ======================
  if (rank == 0) {
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "--- InflowFlux Boundary Condition Setup ---\n");
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Face 0 (-X): %d\n", user->bctype[0]);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Face 1 (+X): %d\n", user->bctype[1]);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Face 2 (-Y): %d\n", user->bctype[2]);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Face 3 (+Y): %d\n", user->bctype[3]);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Face 4 (-Z): %d\n", user->bctype[4]);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Face 5 (+Z): %d\n", user->bctype[5]);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "(Note: INLET=%d, OUTLET=%d, WALL=%d, etc.)\n", INLET, OUTLET, WALL);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "-------------------------------------------\n");
  }
  PetscBarrier(NULL);
  
  // --- Pre-calculation Step ---
  // Some profiles require a single velocity value based on total area (pulsatile)
  // or a pre-calculated radius field (parabolic).
  
  // For pulsatile flow, we need the total area first to get a consistent velocity.
  if (inletprofile == 3) {
    lAreaIn = 0.0;
    for (fn=0; fn<6; fn++) {
      if (user->bctype[fn] == INLET) {
        // This logic is duplicated from the main loop below for area calculation only
        if (fn==0 && xs==0) { i=xs; for(k=lzs;k<lze;k++) for(j=lys;j<lye;j++) if(nvert[k][j][i+1]<0.1) lAreaIn+=sqrt(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z); }
        if (fn==1 && xe==mx) { i=mx-2; for(k=lzs;k<lze;k++) for(j=lys;j<lye;j++) if(nvert[k][j][i]<0.1)   lAreaIn+=sqrt(csi[k][j][i].x*csi[k][j][i].x+csi[k][j][i].y*csi[k][j][i].y+csi[k][j][i].z*csi[k][j][i].z); }
        if (fn==2 && ys==0) { j=ys; for(k=lzs;k<lze;k++) for(i=lxs;i<lxe;i++) if(nvert[k][j+1][i]<0.1) lAreaIn+=sqrt(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z); }
        if (fn==3 && ye==my) { j=my-2; for(k=lzs;k<lze;k++) for(i=lxs;i<lxe;i++) if(nvert[k][j][i]<0.1)   lAreaIn+=sqrt(eta[k][j][i].x*eta[k][j][i].x+eta[k][j][i].y*eta[k][j][i].y+eta[k][j][i].z*eta[k][j][i].z); }
        if (fn==4 && zs==0) { k=zs; for(j=lys;j<lye;j++) for(i=lxs;i<lxe;i++) if(nvert[k+1][j][i]<0.1) lAreaIn+=sqrt(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z); }
        if (fn==5 && ze==mz) { k=mz-2; for(j=lys;j<lye;j++) for(i=lxs;i<lxe;i++) if(nvert[k][j][i]<0.1)   lAreaIn+=sqrt(zet[k][j][i].x*zet[k][j][i].x+zet[k][j][i].y*zet[k][j][i].y+zet[k][j][i].z*zet[k][j][i].z); }
      }
    }
    MPI_Allreduce(&lAreaIn, &AreaSumIn, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    //  Flux_Waveform_Read(user);
    // uin = Pulsatile_Plug_Inlet_Flux(user, AreaSumIn); // uin = Q(t) / A_total
  } 
  // For fully-developed pipe flow, pre-calculate radius at the inlet face
  else if (inletprofile == 4) {
    for (fn=0; fn<6; fn++) {
      if (user->bctype[fn] == INLET && fn==4 && zs==0) { // Assuming profile 4 is only used on face 4 (-z)
        k=0;
        for(j=lys;j<lye;j++){
          for(i=lxs;i<lxe;i++){
            xc = cent[k+1][j][i].x - CMx_c;
            yc = cent[k+1][j][i].y - CMy_c;
            RR[k][j][i] = sqrt(xc*xc + yc*yc); // Store radius in temp array
          }
        }
      }
      // Add similar blocks for other faces if this profile can be used elsewhere
    }
  }
  
  // For other uniform profiles, get the max velocity now.
  if (inletprofile == 1 || inletprofile == 2 || inletprofile == 4 || inletprofile == 7) {
    PetscOptionsGetReal(NULL,NULL, "-uin", &uin_max, NULL); // Get Umax from options
  }

  // --- Main Application Loop ---
  lAreaIn = 0.0; // Reset for final summation
  for (fn=0; fn<6; fn++) {
   if (user->bctype[fn] == INLET) {
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Inlet detected at face: %d \n", fn);
    switch(fn){
    case 0: // -X face
       if (xs==0) {
       	i = xs;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    if (nvert[k][j][i+1]<0.1) {
              // Calculate uin for this point based on profile
              // Add profile-specific logic here if needed for this face
              if (inletprofile == 1) uin = 1.0;
              // ... other profiles
	      
	      d = sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
              lAreaIn += d;
              ucont[k][j][i].x = uin*d;
	      ubcs[k][j][i].x = uin*csi[k][j][i].x/d;
	      ubcs[k][j][i].y = uin*csi[k][j][i].y/d;
	      ubcs[k][j][i].z = uin*csi[k][j][i].z/d;
	      ucat[k][j][i+1] = ubcs[k][j][i];
	      FluxIn += ucont[k][j][i].x;
	    }
	  }
	}
      }
      break;
    case 1: // +X face
      if (xe==mx) {
	i = mx-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    if (nvert[k][j][i]<0.1) {
              // Calculate uin for this point based on profile
              if (inletprofile == 1) uin = 1.0;
              // ... other profiles

              d = sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
	      lAreaIn += d;
              ucont[k][j][i].x = -uin*d;
	      ubcs[k][j][i+1].x = -uin*csi[k][j][i].x/d;
	      ubcs[k][j][i+1].y = -uin*csi[k][j][i].y/d;
	      ubcs[k][j][i+1].z = -uin*csi[k][j][i].z/d;
	      ucat[k][j][i] = ubcs[k][j][i+1];
	      FluxIn += ucont[k][j][i].x;
	    }
	  }
	}
      }
      break;
    case 2: // -Y face
      if (ys==0) {
	j = ys;
        for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    if (nvert[k][j+1][i]<0.1) {
              // Calculate uin for this point based on profile
              if (inletprofile == 1) uin = 1.0;
              // ... other profiles

	      d = sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);
              lAreaIn += d;
              ucont[k][j][i].y = uin*d;
	      ubcs[k][j][i].x = uin*eta[k][j][i].x/d;
	      ubcs[k][j][i].y = uin*eta[k][j][i].y/d;
	      ubcs[k][j][i].z = uin*eta[k][j][i].z/d;
	      ucat[k][j+1][i] = ubcs[k][j][i];
	      FluxIn += ucont[k][j][i].y;
	    }
	  }
	}
      }
      break;
    case 3: // +Y face
      if (ye==my) {
	j = my-2;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    if (nvert[k][j][i]<0.1) {
              // Calculate uin for this point based on profile
              if (inletprofile == 1) uin = 1.0;
              // ... other profiles

 	      d = sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);   
	      lAreaIn += d;
              ucont[k][j][i].y = -uin*d;
	      ubcs[k][j+1][i].x = -uin*eta[k][j][i].x/d;
	      ubcs[k][j+1][i].y = -uin*eta[k][j][i].y/d;
	      ubcs[k][j+1][i].z = -uin*eta[k][j][i].z/d;
	      ucat[k][j][i] = ubcs[k][j+1][i];
	      FluxIn += ucont[k][j][i].y;
	    }
	  }
	}
      }
      break;
    case 4: // -Z face
      if (zs==0) {
	k = 0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if (nvert[k+1][j][i]<0.1) {
              // This is where most profiles from the original code were defined.
              // Calculate coordinates and radius for this specific point.
	      xc = (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25 - CMx_c;
	      yc = (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25 - CMy_c;
              // NOTE: original code used zc*zc + yc*yc for r. This depends on problem setup. Assuming xy-plane radius.
	      r = sqrt(xc*xc + yc*yc);
              
              // *** MERGED INLET PROFILE LOGIC ***
              if (inletprofile == 0 || inletprofile == 6 || inletprofile == 8) {
                uin = InletInterpolation(r, user);
              } else if (inletprofile == 1) {
                uin = 1.0;
              } else if (inletprofile == -1) {
                uin = -1.0;
              } else if (inletprofile == 2) { // 2D channel flow from target code
                uin = 4.0 * uin_max * yc * (H - yc) / (H * H);
              } else if (inletprofile == 3) {
                // uin already calculated for pulsatile plug flow
              } else if (inletprofile == 4) { // Fully-developed pipe flow
                r = RR[k][j][i]; // Use pre-calculated radius
                uin = uin_max * (1.0 - 4.0 * r * r); // Assumes max radius is 0.5
                if (r > 0.5) uin=0.; 
              } else if (inletprofile == 5) {
                uin = -InletInterpolation(r, user);
              } else if (inletprofile == 7) {
                uin = uin_max * (1.0 - 4.0 * yc * yc);
              } else if (inletprofile == 10) {
                double _y = (yc-1.5)*2., _x=xc-0.5;
                double A=1.;
                #ifndef M_PI 
                #define M_PI 3.14159265358979323846264338327950288
                #endif
                uin = 1 - _y*_y;
                int n;
                for(n=0; n<20; n++) {
                  uin -= 4.*pow(2./M_PI,3) * pow(-1., n) / pow(2*n+1.,3.) * cosh((2*n+1)*M_PI*_x/2.) * cos((2.*n+1)*M_PI*_y/2.) / cosh((2*n+1)*M_PI*A/2.);
                }
              } else if (inletprofile == 11) {
                uin = 0.185;
              } else {
                // LOG_ALLOW(LOCAL,LOG_DEBUG, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
                uin = 0.;
              }

	      d = sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
              lAreaIn += d;
              ucont[k][j][i].z = uin*d;
	      ubcs[k][j][i].x = uin*zet[k][j][i].x/d;
	      ubcs[k][j][i].y = uin*zet[k][j][i].y/d;
	      ubcs[k][j][i].z = uin*zet[k][j][i].z/d;
	      ucat[k+1][j][i] = ubcs[k][j][i];
	      FluxIn += ucont[k][j][i].z;
	    }
	  }
	}
      }

      LOG_LOOP_ALLOW_EXACT(LOCAL,LOG_DEBUG,i+j,10,"\n",FluxIn);
      break;
    case 5: // +Z face
      if (ze==mz) {	
	k = mz-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    if (nvert[k][j][i]<0.1) {
              // Add profile-specific logic here if needed for this face
              if (inletprofile == 1) uin = 1.0;
              // ... other profiles

	      d = sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
              lAreaIn += d;
              ucont[k][j][i].z = -uin*d;
	      ubcs[k+1][j][i].x = -uin*zet[k][j][i].x/d;
	      ubcs[k+1][j][i].y = -uin*zet[k][j][i].y/d;
	      ubcs[k+1][j][i].z = -uin*zet[k][j][i].z/d;
	      ucat[k][j][i] = ubcs[k+1][j][i];
	      FluxIn += ucont[k][j][i].z;
	    }
	  }
	}
      }

       LOG_LOOP_ALLOW_EXACT(LOCAL,LOG_DEBUG,i+j,10,"\n",FluxIn);
      break;
    }//end switch
   }// end inlet check
   else if(user->bctype[fn]==WALL) {
     LOG_ALLOW(GLOBAL,LOG_DEBUG,"Solid Wall detected at face: %d \n",fn);
   }
   else if(user->bctype[fn]==SYMMETRY) {
     LOG_ALLOW(GLOBAL,LOG_DEBUG,"Symmetry detected at face: %d \n",fn);
   }
   else if(user->bctype[fn]==OUTLET){
     LOG_ALLOW(GLOBAL,LOG_DEBUG,"Outlet detected at face: %d \n",fn);
   }
  }// end face counter 

  // Sum flux and area from all processes
  MPI_Allreduce(&FluxIn, &simCtx->FluxInSum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  if (inletprofile != 3) { // AreaSumIn already computed for pulsatile
     MPI_Allreduce(&lAreaIn, &AreaSumIn, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  }
  PetscBarrier(NULL);
  
  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Inflow Flux - Area:  %le - %le \n", simCtx->FluxInSum, AreaSumIn);    
  
  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, RFC, &RR);
  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  VecDestroy(&RFC);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  
  return 0; 
}
*/

//////////////////////////////

static PetscReal InletInterpolation(PetscReal r2, UserCtx *user)
{
  PetscInt i;
  PetscReal temp;
  PetscInt tstep, ts_p_cycle=2000;
  PetscInt inletprofile = user->simCtx->inletprofile;
  PetscInt ti = user->simCtx->step;
  PetscReal *r = user->simCtx->r;
  PetscReal *tin = user->simCtx->tin;
  PetscReal **uinr = user->simCtx->uinr;
  PetscReal Flux_in = user->simCtx->Flux_in;

  tstep = ti/2 - ((ti / ts_p_cycle) * ts_p_cycle);

  if (inletprofile == 8) {
    user->simCtx->Flux_in=sin(2*3.14159*ti/200.);
    return(Flux_in);
  }
  
  if (inletprofile==3 || inletprofile==6) 
    return(Flux_in);
  
  if (r2>1.) {
    temp = uinr[99][tstep];
    return(temp);
  }
  for (i=0; i<100; i++) {
    if (r2>= (r[i]) && r2< (r[i+1])) {
      temp = uinr[i][tstep] + (uinr[i+1][tstep] - uinr[i][tstep]) *
	(r2-r[i]) / (r[i+1]-r[i]);
    
      return (temp);
    }
  }
  return 0;   
}

////////////////////////////////


PetscErrorCode OutflowFlux(UserCtx *user) {
  
  PetscInt i, j, k;
  PetscReal FluxOut;
 
  Cmpnts	***ucont, ***ucat, ***zet;

  DM fda = user->fda;
  DMDALocalInfo	info = user->info;
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
  

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(fda, user->lZet,  &zet);
  
  FluxOut = 0;

  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Block %d calculating Outflow Flux .\n",user->_this);
  
  if (user->bctype[5] == 4 || user->bctype[5] == 0 || user->bctype[5] == 14) {    
    if (ze==mz) {
      k = mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }

  } else if (user->bctype[4] == 4 || user->bctype[4] == 0) {    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }
  }
  MPI_Allreduce(&FluxOut,&user->simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," Block %d FluxOutSum = %.6f .\n",user->_this,user->simCtx->FluxOutSum);

  DMDAVecRestoreArray(fda, user->Ucont, &ucont);

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  return 0;
}


def execute_command(command: list, run_dir: str, log_filename: str, monitor_cfg: dict):
    """!
    @brief Executes a command, streaming its output to the console and a log file.
    @param[in] command A list of strings representing the command and its arguments.
    @param[in] run_dir The directory in which to execute the command.
    @param[in] log_filename The name of the file to save the command's output to.
    @param[in] monitor_cfg The parsed monitor.yml dictionary, used to set LOG_LEVEL.
    @throws SystemExit if the command fails or is not found.
    """
    log_path = os.path.join(run_dir, "logs", log_filename)
    print(f"[INFO] Launching Command...\n  > {' '.join(command)}")
    print(f"       Log file: {os.path.relpath(log_path)}")
    print("-" * 60)
    run_env = os.environ.copy()
    verbosity = monitor_cfg.get('logging', {}).get('verbosity', 'ERROR').upper()
    run_env['LOG_LEVEL'] = verbosity
    print(f"[INFO] Setting LOG_LEVEL={verbosity} for C executable.")
    print("-" * 60)
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            cwd=run_dir, env=run_env, bufsize=1, universal_newlines=True, encoding='utf-8', errors='replace'
        )
        with open(log_path, "w") as log_file:
            for line in process.stdout:
                sys.stdout.write(line)
                log_file.write(line)
        process.wait()
        return_code = process.returncode
        print("-" * 60)
        if return_code == 0:
            print(f"[SUCCESS] Execution finished successfully.")
        else:
            print(f"[FATAL] Execution failed with exit code {return_code}. Check log: {os.path.relpath(log_path)}", file=sys.stderr)
            sys.exit(return_code)
    except FileNotFoundError:
        print(f"[FATAL] Command not found or is not executable: '{command[0]}'", file=sys.stderr)
        print("        Please check that the path is correct and the file has execute permissions.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[FATAL] An unexpected error occurred during execution: {e}", file=sys.stderr)
        sys.exit(1)

# ==============================================================================
# --- TEMPORARY DIAGNOSTIC FUNCTION ---
# ==============================================================================
def build_project(args):
    """!
    @brief Diagnostic function to inspect the Python process environment.
    """
    print("\n" + "="*20 + " PYTHON ENVIRONMENT DIAGNOSTICS " + "="*20)
    
    # Print the exact Python executable being used
    print(f"\n[INFO] Python executable path (from sys.executable):")
    print(f"  > {sys.executable}")

    # Print the PATH environment variable as Python sees it
    python_path = os.environ.get('PATH')
    print(f"\n[INFO] The PATH environment variable available to Python:")
    if python_path:
        # Print each path on a new line for clarity
        for p in python_path.split(':'):
            print(f"  > {p}")
    else:
        print("  > [FATAL] The PATH environment variable is not set at all!")

    # Use the definitive tool, shutil.which, to see what Python can find
    print("\n[INFO] Using shutil.which() to search the above PATH:")
    bash_location = shutil.which('bash')
    sh_location = shutil.which('sh')
    ls_location = shutil.which('ls')

    print(f"  > Location of 'bash': {bash_location}")
    print(f"  > Location of 'sh':   {sh_location}")
    print(f"  > Location of 'ls':   {ls_location}")

    # Check for the existence of /bin/bash directly
    print("\n[INFO] Checking existence of absolute paths using os.path.isfile():")
    bin_bash_exists = os.path.isfile('/bin/bash')
    bin_sh_exists = os.path.isfile('/bin/sh')
    print(f"  > Does '/bin/bash' exist? {bin_bash_exists}")
    print(f"  > Does '/bin/sh' exist?   {bin_sh_exists}")

    print("\n" + "="*25 + " DIAGNOSTICS COMPLETE " + "="*25)
    sys.exit(0) # Exit cleanly after printing diagnostics

////////////////////////////////

/**
 * @brief Interpolates a cell-centered field (scalar or vector) onto DMSwarm particles,
 *        converting the cell-center data to corner data first, then looping over particles.
 *
 * Steps:
 *   1) Check that the Vec has blockSize=1 or 3 (scalar vs. vector).
 *   2) Map the cell-centered Vec to a local array (fieldGlobal -> localPtr).
 *   3) Allocate a corner array (cornerPtr) via Allocate3DArray(...), sized (zm+1, ym+1, xm+1).
 *   4) Convert from cell-centers to corners via InterpolateFieldFromCenterToCorner(...).
 *   5) Restore the cell-centered array.
 *   6) Retrieve DMSwarm fields: "DMSwarm_CellID", "weight", and swarmOutFieldName.
 *   7) Loop over local particles, clamp i/j/k, skip or zero out if out of range, read (a1,a2,a3).
 *   8) Call InterpolateEulerFieldToSwarmForParticle(...) with cornerPtr to do final Interpolatetion.
 *   9) Restore swarm fields, free the corner array.
 *
 * @param[in]  user              User context with:
 *                                - user->da     (cell-centered DMDA),
 *                                - user->swarm  (DMSwarm).
 * @param[in]  fieldGlobal       Vec with blockSize=1 or 3, storing the cell-centered field.
 * @param[in]  fieldName         Human-readable field name for logging (e.g. "velocity").
 * @param[in]  swarmOutFieldName Name of the DMSwarm field where Interpolatetion results go.
 *
 * @return PetscErrorCode  0 on success, non-zero on error.
 */
/*
PetscErrorCode InterpolateEulerFieldToSwarm(
    UserCtx    *user,
    Vec         fieldGlobal,       //DMDA Vec containing cell‐center data 
    const char *fieldName,         // e.g., "Ucat" 
    const char *swarmOutFieldName) // Name of the output DMSwarm field 
{
  PetscErrorCode ierr;
  DM             fda    = user->fda;      // DM for cell‐center field data 
  //  DM             da     = user->da;        DM for grid information (local indices) 
  DM             swarm  = user->swarm;     DMSwarm for particles 
  PetscInt       bs;                      Block size: 1 (scalar) or 3 (vector) 
  DMDALocalInfo  info;                  // Local grid info 
void          *localPtr   = NULL;       // Pointer to cell‐center data from fda 
void          *cornerPtr  = NULL;       // Will hold the typed corner array 
  void          *swarmOut   = NULL;     // Pointer to the swarm output field 
  PetscInt    *cellIDs    = NULL;     // Particle cell indices from swarm 
  PetscReal     *weights    = NULL;     // Interpolation coefficients from swarm 
  PetscInt       nLocal;
  PetscMPIInt    rank;

  
  
  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  // (A) Check block size and get local domain info 
  ierr = VecGetBlockSize(fieldGlobal, &bs); CHKERRQ(ierr);
  if (bs != 1 && bs != 3) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
             "InterpolateEulerFieldToSwarm: blockSize must be 1 or 3, got %d.", (PetscInt)bs);
  }
  ierr = DMDAGetLocalInfo(fda, &info); CHKERRQ(ierr);

  PetscInt xs_node = info.xs;
  PetscInt ys_node = info.ys;
  PetscInt zs_node = info.zs;
  
  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Starting with field='%s', blockSize=%d, local domain: (%d x %d x %d)\n",
    fieldName, bs, info.xm, info.ym, info.zm);

  // (B) Map the cell-centered Vec to a local array using the DM attached to fieldGlobal 
  LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d] Dumping DM state (user->fda) BEFORE GetArray:\n", rank);

  if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) DMView(user->fda, PETSC_VIEWER_STDOUT_SELF);

  ierr = DMDAVecGetArrayRead(fda, fieldGlobal, &localPtr); CHKERRQ(ierr);

    // Corner domain is larger than cell center domain 
    PetscInt nz = info.zm; 
    PetscInt ny = info.ym;
    PetscInt nx = info.xm;
    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "[Rank %d]  Allocating corner array of size (%d x %d x %d), blockSize = %d\n",
	      rank,nx, ny, nz, bs);
    if (bs == 1) {
      // Declare a typed Pointer for scalar corners 
      PetscReal ***cornerScal = NULL;
      ierr = Allocate3DArray(&cornerScal, nz, ny, nx); CHKERRQ(ierr);
      if (!cornerScal) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateEulerFieldToSwarm: 3D Array Allocation for scalar corners failed.\n");
      }
      // Save typed Pointer into cornerPtr so later code can cast appropriately 
      cornerPtr = (void*) cornerScal;
      // (D) Convert cell-center data to corners for scalar field 
      ierr = InterpolateFieldFromCenterToCorner( (PetscReal ***) localPtr, cornerScal, user); CHKERRQ(ierr);
    } else {
      // For vector fields 
      Cmpnts ***cornerVec = NULL;
      ierr = Allocate3DArray(&cornerVec, nz, ny, nx); CHKERRQ(ierr);
      if (!cornerVec) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateEulerFieldToSwarm: 3D Array Allocation for vector corners failed.\n");
      }
      cornerPtr = (void*) cornerVec;

      // Comment out the actual call:
      ierr = InterpolateFieldFromCenterToCorner( (Cmpnts ***) localPtr, (Cmpnts ***)cornerPtr, user); CHKERRQ(ierr);
      */
      /*
      // --- DEBUG - BYPASS SECTION ---
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[%d] DEBUG: Bypassing InterpolateFieldFromCenterToCorner call.\n", rank);


      // Instead, populate cornerPtr with fixed values
      // Use local indices 0..nz-1, 0..ny-1, 0..nx-1 because cornerPtr is a local array
      if (bs == 3) {
	Cmpnts ***cornerVec_bypass = (Cmpnts ***)cornerPtr;
	for (PetscInt k_local = 0; k_local < nz; ++k_local) {
          for (PetscInt j_local = 0; j_local < ny; ++j_local) {
	    for (PetscInt i_local = 0; i_local < nx; ++i_local) {
	      cornerVec_bypass[k_local][j_local][i_local].x = 1.0; // Or 0.0, or rank number
	      cornerVec_bypass[k_local][j_local][i_local].y = 2.0;
	      cornerVec_bypass[k_local][j_local][i_local].z = 3.0;
	    }
          }
	}
	LOG_ALLOW(LOCAL,LOG_DEBUG, "[%d] DEBUG: Finished setting fixed values in cornerVec_bypass.\n", rank);
      
      }
      */
      /*
    LOG_ALLOW(LOCAL, LOG_INFO,
      "InterpolateEulerFieldToSwarm: Rank %d Completed center-to-corner Interpolatetion for field='%s'.\n",
	      rank,fieldName);
    }
  // (E) Restore the cell-centered array since we now have corner data 
  ierr = DMDAVecRestoreArrayRead(fda, fieldGlobal, &localPtr); CHKERRQ(ierr);

  // (F) Retrieve swarm fields: cell IDs, weights, and the output field 
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_CellID",  NULL, NULL, (void**)&cellIDs);  CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "weight",          NULL, NULL, (void**)&weights);  CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut);         CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Found %d local particles to process for field='%s'.\n",
    nLocal, fieldName);

  // (G) Loop over each local particle and perform the final Interpolatetion from corners 
  for (PetscInt p = 0; p < nLocal; p++) {
    PetscInt iCell_global = (PetscInt)cellIDs[3*p + 0];
    PetscInt jCell_global = (PetscInt)cellIDs[3*p + 1];
    PetscInt kCell_global = (PetscInt)cellIDs[3*p + 2];

// --- Convert GLOBAL cell index to LOCAL cell index ---
      // NOTE: This assumes the cell index corresponds directly to the
      //       "lower-left-front" node index for interpolation purposes.
      //       We need the local index relative to the cornerPtr array,
      //       which is indexed from 0 based on owned NODES.
    
      PetscInt iCell_local = iCell_global - xs_node;
      PetscInt jCell_local = jCell_global - ys_node;
      PetscInt kCell_local = kCell_global - zs_node;

      // --- Clamp LOCAL indices and check bounds for the BASE cell index ---
      // This check ensures the base index (i,j,k) for the interpolation
      // is within the bounds of the locally allocated cornerPtr array.
      // The interpolation function will handle the +1 offsets.
    // Boundary clamp: adjust indices to be within [0, mx), [0, my), [0, mz) 
    if (iCell_local >= user->info.mx) iCell_local = user->info.mx - 1;
    if (jCell_local >= user->info.my) jCell_local = user->info.my - 1;
    if (kCell_local >= user->info.mz) kCell_local = user->info.mz - 1;
    
    if (iCell_local < 0 ||
	jCell_local < 0 ||
	kCell_local < 0 ||
        iCell_local >= user->info.xm-1 ||
        jCell_local >= user->info.ym-1 ||
        kCell_local >= user->info.zm-1)
    {

      LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] Rank/Domain boundary reached(cell wise): (k,j,i) - (%d,%d,%d).\n",rank,kCell_local,jCell_local,iCell_local);
      
      // Out-of-range: set output to zero       
      if (bs == 1) {
        ((PetscReal*)swarmOut)[p] = 0.0;
      } else {
        ((PetscReal*)swarmOut)[3*p + 0] = 0.0;
        ((PetscReal*)swarmOut)[3*p + 1] = 0.0;
        ((PetscReal*)swarmOut)[3*p + 2] = 0.0;
      }
      continue;
    }
    
    
    // Retrieve Interpolatetion coefficients (a1, a2, a3) for this particle 
    PetscReal alpha1 = weights[3*p + 0];
    PetscReal alpha2 = weights[3*p + 1];
    PetscReal alpha3 = weights[3*p + 2];
      */
    /* (Optional) If your final Interpolatetion expects a corner offset, adjust here.
       For example, if the cell center corresponds to the average of corners at (iCell,jCell,kCell)
       and (iCell+1, jCell+1, kCell+1), you might add +1. For now, we pass indices as is.
    */
     /*  
    PetscInt iUse = iCell_local;  // + 1 if required
    PetscInt jUse = jCell_local;  // + 1 if required
    PetscInt kUse = kCell_local;  // + 1 if required
     */
    /* (H) Call the per-particle Interpolatetion function.
       This function will use your _Generic macro (TrilinearInterpolation or PiecewiseLinearInterpolation)
       on the corner data.
    */
     /* 
    ierr = InterpolateEulerFieldToSwarmForParticle(
              fieldName,    // e.g., "Ucat" 
              cornerPtr,    // typed Pointer: (PetscReal***) or (Cmpnts***) 
              iUse, jUse, kUse,
              alpha1, alpha2, alpha3,
              swarmOut,     // Pointer to swarm output array 
              p,            // particle index 
              bs);          // block size: 1 or 3 
    CHKERRQ(ierr);
  }

  // (I) Restore swarm fields 
  ierr = DMSwarmRestoreField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut);        CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "weight",          NULL, NULL, (void**)&weights); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID",  NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);

  // (J) Deallocate the corner array using the generic deallocation macro 
  if (bs == 1) {
    PetscReal ***cornerScal = (PetscReal ***) cornerPtr;
    ierr = Deallocate3DArray(cornerScal, info.zm, info.ym); CHKERRQ(ierr);
  } else {
    Cmpnts ***cornerVec = (Cmpnts ***) cornerPtr;
    ierr = Deallocate3DArray(cornerVec, info.zm, info.ym); CHKERRQ(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Rank %d Completed Interpolatetion of field='%s' for %d local particles.\n",
	    rank,fieldName, nLocal);

  PetscFunctionReturn(0);
}
*/

///////////////////

//////////////////////////////////

/**
 * @brief Resets the location-dependent state of a loaded swarm to force relocation.
 * @ingroup ParticleRestart
 *
 * This function is a critical part of the simulation restart procedure. It must be
 * called immediately after `ReadAllSwarmFields` has populated a swarm from restart
 * files. Its purpose is to invalidate the "location" state of the loaded particles,
 * ensuring that the `LocateAllParticlesInGrid_TEST` orchestrator performs a fresh,
 * comprehensive search for every particle based on its loaded position.
 *
 * It does this by performing two actions on every locally-owned particle:
 * 1.  It resets the `DMSwarm_CellID` field to a sentinel value of `(-1, -1, -1)`.
 *     This invalidates any cell index that might have been loaded or defaulted to 0.
 * 2.  It sets the `DMSwarm_location_status` field to `NEEDS_LOCATION`.
 *
 * This guarantees that the location logic will not mistakenly use a stale cell index
 * from a previous run and will instead use the robust "Guess -> Verify" strategy
 * appropriate for particles with unknown locations.
 *
 * @param[in,out] user Pointer to the UserCtx structure which contains the `DMSwarm` object
 *                     that has just been loaded with data from restart files.
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code if field access fails.
 */
PetscErrorCode PrepareLoadedSwarmForRelocation(UserCtx *user)
{
    PetscErrorCode ierr;
    DM             swarm;
    PetscInt       n_local;
    PetscInt      *cell_p;    // Pointer to the raw data for the CellID field
    PetscInt      *status_p;  // Pointer to the raw data for the location_status field
    PetscInt64    *PIDs;      // Pointer to the raw data for the Particle ID field.
    PetscMPIInt   rank,size;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);  CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);  CHKERRQ(ierr);

    SimCtx *simCtx = user->simCtx;
    
    // --- 1. Input Validation and Setup ---
    if (!user || !user->swarm) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx or its DMSwarm is NULL in PrepareLoadedSwarmForRelocation.");
    }
    swarm = user->swarm;

    // Get the number of particles on this MPI rank.
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);

    // If there are no local particles, there is nothing to do.
    if (n_local == 0) {
        PetscFunctionReturn(0);
    }

    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Preparing %d loaded particles for relocation by resetting their CellID and Status.\n", n_local);

    // --- 2. Get Writable Access to Swarm Fields ---
    // This provides direct pointers to the underlying data arrays for the fields.
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_p);   CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_p); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr); CHKERRQ(ierr);

    // --- 3. Determine Starting Global PID for this Rank ---
    PetscInt particles_per_rank_ideal = simCtx->np / size; // Assumes user->size is PETSC_COMM_WORLD size
    PetscInt remainder_particles = simCtx->np % size;
    PetscInt base_pid_for_rank = rank * particles_per_rank_ideal + PetscMin(rank, remainder_particles);
    // This calculation must match how particlesPerProcess was determined (e.g., in DistributeParticles).
    
    // --- 4. Loop Through All Local Particles and Reset State ---
    for (PetscInt p = 0; p < n_local; ++p) {

        
        // Reset the 3 components of the cell index vector.
        cell_p[3*p + 0] = -1;
        cell_p[3*p + 1] = -1;
        cell_p[3*p + 2] = -1;

        // Reset the status to ensure it will be processed by the location algorithm.
        status_p[p] = (ParticleLocationStatus)UNINITIALIZED;

    //set the PID for each particle
    PIDs[p] =  (PetscInt64)base_pid_for_rank + p;
    }

    // --- 4. Restore Fields ---
    // This returns control of the data arrays back to the DMSwarm. It is a mandatory
    // step to ensure data consistency and prevent memory issues.
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_p);   CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_p); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Successfully reset location state for all loaded particles.\n");

    PetscFunctionReturn(0);
}

/////////////////////////////////////////////////////////////////



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

    const PetscInt physical_nodes_in_dim = GlobalNodesInDim_from_info - 1; // Since last node is not really "physical".

    // --- 3. Handle Edge Cases for Global Domain Size ---
    // If the global domain has 0 or 1 node plane, no cells can be formed.
    if (physical_node_in_dim <= 1) {
        *xs_cell_global_out = xs_node_global_rank; // Still report the rank's starting node
        *xm_cell_local_out = 0;                    // But 0 cells
        PetscFunctionReturn(0);
    }
    // Negative global dimension is an error (should be caught by DMDA setup, but defensive)
    if (GlobalNodesInDim_from_info < 0 ) {
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "GlobalNodesInDim %d from DMDALocalInfo must be non-negative for dimension %d.", physical_node_in_dim, dim);
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
        // So, the last valid cell origin node index is (Physical_node_in_dim - 2).
        PetscInt last_possible_origin_global_idx = physical_node_in_dim - 2;

        // Determine the range of nodes owned by this rank that could *potentially* be cell origins.
        // The first node owned by the rank is a potential origin.
        PetscInt first_potential_origin_on_rank = xs_node_global_rank;

        // A node can be an origin if there's at least one node after it to form the cell.
        // So, the last node owned by the rank that could *potentially* be an origin is
        // the second-to-last node it owns: (xs_node_global_rank + num_nodes_owned_rank - 1) - 1
        // which simplifies to: xs_node_global_rank + num_nodes_owned_rank - 2.
        PetscInt last_potential_origin_on_rank = xs_node_global_rank + num_nodes_owned_rank - 1;

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

////////////////////////////////////


#undef __FUNCT__
#define __FUNCT__ "InterpolateFieldFromCenterToCorner_Scalar"
// -----------------------------------------------------------------------------
// Interpolation: Center -> Corner (scalar)
// -----------------------------------------------------------------------------
/**
 * @brief Interpolates a scalar field from cell centers to corner nodes.
 *
 * This function estimates the value of a scalar field at each grid node by averaging
 * the values from the cell centers of the cells surrounding that node (up to 8).
 * It handles physical boundaries by averaging only the available adjacent cells.
 *
 * Assumes input `centfield_arr` is from a ghosted local vector associated with `user->da` (DOF=1, s=2)
 * and output `field_arr` is from a ghosted local vector also associated with `user->da` (DOF=1, s=2).
 * Input array uses GLOBAL cell indices, output array uses GLOBAL node indices.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of scalar data at cell centers,
 *                            accessed via GLOBAL cell indices (k=0..KM-1, j=0..JM-1, i=0..IM-1).
 * @param[out] field_arr      Output: 3D array (ghosted) where interpolated node values are stored,
 *                            accessed via GLOBAL node indices (k=0..KM, j=0..JM, i=0..IM).
 * @param[in]  user           User context containing DMDA information (da).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar(
    PetscReal ***centfield_arr, /* Input: Ghosted local array based on da (cell indices) */
    PetscReal ***field_arr,     /* Output: Ghosted local array based on da (node indices) */
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
      "Rank %d starting interpolation.\n", rank);

    // Get local info based on the DMDA (da). This info primarily describes owned nodes.
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Define the range of NODES owned by this processor using GLOBAL indices
    PetscInt xs_node = info.xs;
    PetscInt xm_node = info.xm;
    PetscInt xe_node = xs_node + xm_node;
    PetscInt ys_node = info.ys;
    PetscInt ym_node = info.ym;
    PetscInt ye_node = ys_node + ym_node;
    PetscInt zs_node = info.zs;
    PetscInt zm_node = info.zm;
    PetscInt ze_node = zs_node + zm_node;

    // Get Global dimensions (number of cells IM, JM, KM)
    PetscInt IM = info.mx - 1;
    PetscInt JM = info.my - 1;
    PetscInt KM = info.mz - 1;


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank %d: Interpolating for owned NODES k=%d..%d, j=%d..%d, i=%d..%d\n",
              rank, zs_node, ze_node, ys_node, ye_node, xs_node, xe_node);

    // Loop over the GLOBAL indices of the NODES owned by this processor
    for (PetscInt k = zs_node; k < ze_node; k++) { // Global NODE index k (0..KM)
        for (PetscInt j = ys_node; j < ye_node; j++) { // Global NODE index j (0..JM)
            for (PetscInt i = xs_node; i < xe_node; i++) { // Global NODE index i (0..IM)

                PetscReal sum = 0.0;
                PetscInt count = 0;

                // Loop over the potential 8 CELL indices surrounding node (i,j,k)
                // Cell indices are (i & i-1), (j & j-1), (k & k-1)
                for (PetscInt dk = -1; dk <= 0; dk++) { // Relative cell k-index offset
                    for (PetscInt dj = -1; dj <= 0; dj++) { // Relative cell j-index offset
                        for (PetscInt di = -1; di <= 0; di++) { // Relative cell i-index offset

                            PetscInt ci = i + di; // Global CELL index of neighbor
                            PetscInt cj = j + dj; // Global CELL index of neighbor
                            PetscInt ck = k + dk; // Global CELL index of neighbor

                            // Check if this CELL index is within the valid global cell range (0..IM-1, etc.)
                            if (ci >= 0 && ci <= IM-1 &&
                                cj >= 0 && cj <= JM-1 &&
                                ck >= 0 && ck <= KM-1)
                            {
                                // Access the input 'centfield_arr' using GLOBAL cell indices.
                                // Relies on centfield_arr being from a ghosted local vector.
                                sum += centfield_arr[ck][cj][ci];
                                count++;
                                // Optional Debug Log
                                // LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, i*j*k, 10000,
                                // " Rank %d | Node(i,j,k)=%d,%d,%d | Using Cell(ci,cj,ck)=%d,%d,%d | Value=%.3f | Count=%d\n",
                                // rank, i, j, k, ci, cj, ck, centfield_arr[ck][cj][ci], count);
                            }
                        } // end di loop
                    } // end dj loop
                } // end dk loop

                // Calculate average and store in the output 'field_arr' at the NODE index (i,j,k)
                if (count > 0) {
                    // Store the result using the GLOBAL node index [k][j][i]
                    field_arr[k][j][i] = sum / (PetscReal)count;
                } else {
                    // This indicates an issue - a node should always be adjacent to at least one cell
                    LOG_ALLOW(GLOBAL, LOG_ERROR,
                              "Rank %d: Node (i=%d,j=%d,k=%d) had count=0 surrounding cells! Check logic/ghosting.\n", rank, i, j, k);
                    // Assign a default value or handle error
                    field_arr[k][j][i] = 0.0; // Defaulting to zero might hide issues
                }
            } // End loop i (nodes)
        } // End loop j (nodes)
    } // End loop k (nodes)

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
              "Rank %d completed interpolation.\n", rank);
    
    
    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InterpolateFieldFromCenterToCorner_Vector" 

// -----------------------------------------------------------------------------
// Interpolation: Center -> Corner (vector)
// -----------------------------------------------------------------------------

/**
 * @brief Interpolates a vector field from cell centers to corner nodes.
 *
 * This function estimates the value of a vector field at each grid node by averaging
 * the vector values from the cell centers of the cells surrounding that node (up to 8).
 * It handles physical boundaries by averaging only the available adjacent cells.
 *
 * Assumes input `centfield_arr` is from a ghosted local vector (e.g., representing ucat,
 * stored using node-indexing convention) and output `field_arr` is a ghosted local
 * vector associated with `user->fda` (DOF=3, s=2), accessed using global node indices.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of vector data conceptually at cell centers,
 *                            accessed via GLOBAL indices respecting the storage convention
 *                            (e.g., `ucat[k][j][i]` uses node index `i` but represents cell `C(i,j,k)` for interior).
 * @param[out] field_arr      Output: 3D array (ghosted) where interpolated node values are stored,
 *                            accessed via GLOBAL node indices (k=0..KM, j=0..JM, i=0..IM).
 * @param[in]  user           User context containing DMDA information (da and fda).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector( // Or original name if you prefer
    Cmpnts ***centfield_arr, /* Input: Ghosted local array (fda-based, cell data at node indices) */
    Cmpnts ***field_arr,     /* Output: local array (fda-based, true node data) */
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    DMDALocalInfo info_nodes;
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    // Node ownership range (GLOBAL indices)
    PetscInt xs_node = info_nodes.xs, xm_node = info_nodes.xm, xe_node = xs_node + xm_node;
    PetscInt ys_node = info_nodes.ys, ym_node = info_nodes.ym, ye_node = ys_node + ym_node;
    PetscInt zs_node = info_nodes.zs, zm_node = info_nodes.zm, ze_node = zs_node + zm_node;

    // Global grid dimensions (NODES) - Used for global cell index check
    PetscInt MX_node = user->IM; //info_nodes.mx;
    PetscInt MY_node = user->JM; //info_nodes.my;
    PetscInt MZ_node = user->KM; //info_nodes.mz;
    PetscInt IM = MX_node - 1; // Max cell index i
    PetscInt JM = MY_node - 1; // Max cell index j
    PetscInt KM = MZ_node - 1; // Max cell index k


    // Valid range for accessing the INPUT ghosted array (using NODE indices)
    PetscInt gxs_node = info_nodes.gxs, gxm_node = info_nodes.gxm, gxe_node = gxs_node + gxm_node;
    PetscInt gys_node = info_nodes.gys, gym_node = info_nodes.gym, gye_node = gys_node + gym_node;
    PetscInt gzs_node = info_nodes.gzs, gzm_node = info_nodes.gzm, gze_node = gzs_node + gzm_node;

    // Log only if this function is allowed by the list set in main()
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank %d: Interpolating for owned NODES k=%d..%d, j=%d..%d, i=%d..%d\n",
              rank, zs_node, ze_node-1, ys_node, ye_node-1, xs_node, xe_node-1);

    LOG_ALLOW(GLOBAL,LOG_DEBUG, "[%d] Dumping DM state (user->fda) BEFORE GetArray:\n", rank);

    
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){
    DMView(user->fda, PETSC_VIEWER_STDOUT_SELF);
    // Inside InterpolateFieldFromCenterToCorner_Vector, before the loops:
    
    PetscMPIInt    rank_check;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_check);
    if (rank_check == 1) { // Only for Rank 1
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] Attempting test read of OWNED INTERIOR centfield_arr[3][1][1]\n");
      Cmpnts test_val_owned_interior = centfield_arr[7][1][1];
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] SUCCESS reading owned interior: x=%f\n", test_val_owned_interior.x);
      
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] Attempting test read of OWNED BOUNDARY centfield_arr[3][0][0]\n");
      Cmpnts test_val_owned_boundary = centfield_arr[7][0][0];
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] SUCCESS reading owned boundary: x=%f\n", test_val_owned_boundary.x);

      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] Attempting test read of GHOST centfield_arr[2][0][0]\n");
      Cmpnts test_val_ghost = centfield_arr[7][0][0]; // This is the line that likely crashes
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] SUCCESS reading ghost: x=%f\n", test_val_ghost.x);
      
      }
    }

    ierr =  PetscBarrier(NULL);    
    
    // Proceed with the original loops...
    // Loop over the GLOBAL indices of the NODES owned by this processor (k, j, i)
    for (PetscInt k = zs_node; k < ze_node; k++) {
        for (PetscInt j = ys_node; j < ye_node; j++) {
            for (PetscInt i = xs_node; i < xe_node; i++) {
	      /*
              if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){
	      // --- TEMPORARY TEST --- WORKS ONLY WHEN IM,JM,KM=5 !!!
	      if (rank == 1 && k == 3 && j == 0 && i == 0) {
		LOG_ALLOW(GLOBAL,LOG_DEBUG, "[1] Test read inside loop for (3,0,0): Accessing centfield_arr[2][0][0]\n");
		Cmpnts test_val_loop = centfield_arr[2][0][0]; // The crashing access
		LOG_ALLOW(GLOBAL,LOG_DEBUG, "[1] Test read inside loop SUCCESS: x=%f\n", test_val_loop.x);
	        }
	      // --- END TEMPORARY TEST ---
	      }
	      */
	      
                Cmpnts sum = {0.0, 0.0, 0.0};
                PetscInt count = 0;
                PetscBool attempted_read = PETSC_FALSE; // Flag to track if read was tried

                // Loop over the 8 potential cells surrounding node N(k,j,i)
                for (PetscInt dk_offset = -1; dk_offset <= 0; dk_offset++) {
                    for (PetscInt dj_offset = -1; dj_offset <= 0; dj_offset++) {
                        for (PetscInt di_offset = -1; di_offset <= 0; di_offset++) {

                            // Calculate the NODE index where the relevant cell's data is stored
                            PetscInt node_idx_k = k + dk_offset;
                            PetscInt node_idx_j = j + dj_offset;
                            PetscInt node_idx_i = i + di_offset;

			    
                            // Check if this NODE index corresponds to a valid GLOBAL cell index
                            PetscInt cell_idx_k = node_idx_k; // Cell index is same as node index for storage
                            PetscInt cell_idx_j = node_idx_j;
                            PetscInt cell_idx_i = node_idx_i;
			    

			    
                            if (cell_idx_i >= 0 && cell_idx_i <= MX_node && // Cell index check
                                cell_idx_j >= 0 && cell_idx_j <= MY_node &&
                                cell_idx_k >= 0 && cell_idx_k <= MZ_node)
                            {
                            
			      /*
			    // Check if this NODE index is valid within the GLOBAL node domain (0..Mx-1)
                            // This implicitly checks if the corresponding cell is valid (0..Mx-2)
                            if (node_idx_i >= 0 && node_idx_i < MX_node -1 &&
                                node_idx_j >= 0 && node_idx_j < MY_node -1 &&
                                node_idx_k >= 0 && node_idx_k < MZ_node -1)
                            {
			    */
			      
			    // Check if the NODE index is within the accessible LOCAL+GHOST range of centfield_arr
                                if (node_idx_i >= gxs_node && node_idx_i < gxe_node &&
                                    node_idx_j >= gys_node && node_idx_j < gye_node &&
                                    node_idx_k >= gzs_node && node_idx_k < gze_node)
                                {
                                    // Log attempt just before read
				  LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"PRE-READ: Rank %d targeting Node(k,j,i)=%d,%d,%d. Reading input centfield_arr[%d][%d][%d] (for cell C(%d,%d,%d))\n",
					      rank,
					      k,j,i,
					      node_idx_k,node_idx_j,node_idx_i,
					      cell_idx_k, cell_idx_j, cell_idx_i);

				    attempted_read = PETSC_TRUE; // Mark that we are attempting a read
                                    
                                    // Convert GLOBAL neighbor index to LOCAL index for reading from the ghosted array
                                    PetscInt k_local_read = node_idx_k - gzs_node;
                                    PetscInt j_local_read = node_idx_j - gys_node;
                                    PetscInt i_local_read = node_idx_i - gxs_node;
                                    // ---> READ <---
                                    Cmpnts cell_val = centfield_arr[k_local_read][j_local_read][i_local_read];
//                                  Cmpnts cell_val = centfield_arr[node_idx_k][node_idx_j][node_idx_i];

                                    // Log success immediately after read
                                    LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"POST-READ: Rank %d successful read from [%d][%d][%d] -> (%.2f, %.2f, %.2f)\n",
					      rank,
					      node_idx_k,node_idx_j,node_idx_i,
					      cell_val.x, cell_val.y, cell_val.z);

                                    sum.x += cell_val.x;
                                    sum.y += cell_val.y;
                                    sum.z += cell_val.z;
                                    count++;

                                } else {
                                     LOG_ALLOW(GLOBAL, LOG_WARNING, /* ... Ghost range warning ... */);
                                }
                            } // end global cell check
                        } // end di_offset
                    } // end dj_offset
                } // end dk_offset

		// ---> Convert GLOBAL node indices (k,j,i) to LOCAL indices for field_arr <---
		// field_arr is dimensioned nx * ny * nz (local node dimensions)
		// Global indices (k,j,i) range from (zs,ys,xs) to (ze-1, ye-1, xe-1)
		// Local indices range from 0 to (zm-1, ym-1, xm-1)
		PetscInt k_local = k - zs_node; // Offset by the starting global index
		PetscInt j_local = j - ys_node;
        PetscInt i_local = i - xs_node;
                // Calculate average and write to output node (k,j,i)
                if (count > 0) {
		  LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"PRE-WRITE: Rank %d targeting Node(k,j,i)=%d,%d,%d. Writing avg (count=%d)\n", rank,k,j,i, count);

                     // --- Defensive check (optional but recommended) ---
                     if (k_local < 0 || k_local >= info_nodes.zm ||
                         j_local < 0 || j_local >= info_nodes.ym ||
                         i_local < 0 || i_local >= info_nodes.xm) {
		       SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Calculated local write index out of bounds!");
                     }
                     // --- End check ---

                     // ---> Write using LOCAL indices <---
                     field_arr[k_local][j_local][i_local].x = sum.x / (PetscReal)count;
                     field_arr[k_local][j_local][i_local].y = sum.y / (PetscReal)count;
                     field_arr[k_local][j_local][i_local].z = sum.z / (PetscReal)count;

                     LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"POST-WRITE: Rank %d successful write to field_arr[%d][%d][%d] (local)\n", rank, k_local,j_local,i_local); // Log local indices

		     

                } else {
                     LOG_ALLOW(GLOBAL, LOG_WARNING, // Use WARNING or ERROR
                               "Rank %d: Node (i=%d,j=%d,k=%d) had count=0 surrounding valid cells! Check logic/ghosting. Writing zero.\n", rank, i, j, k);
                     field_arr[k_local][j_local][i_local] = (Cmpnts){0.0, 0.0, 0.0};
                }
                 // Add a log entry even if count=0 or no read was attempted for this node
                 if (!attempted_read && count==0) {
                     LOG_ALLOW(LOCAL, LOG_DEBUG,"NODE-COMPLETE: Rank %d Node(k,j,i)=%d,%d,%d finished loops, no valid cells/reads attempted.\n", rank, k, j, i);
                 } else if (count == 0 && attempted_read) {
                     // This case should ideally not happen if the ghost region check is correct
                      LOG_ALLOW(LOCAL, LOG_ERROR,"NODE-COMPLETE: Rank %d Node(k,j,i)=%d,%d,%d finished loops, attempted reads but count=0!\n", rank, k, j, i);
                 } else {
                     // This is the normal completion case after writing
		   LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"NODE-COMPLETE: Rank %d Node(k,j,i)=%d,%d,%d finished loops and write.\n", rank, k, j, i);
                 }

            } // End loop node i
        } // End loop node j
    } // End loop node k

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, // Use INFO for completion
              "Rank %d completed interpolation function.\n", rank); // Changed message slightly
    
    PROFILE_FUNCTION_END;
              
    PetscFunctionReturn(0);
}


///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

///// WORKING BCS CODE: GROUND TRUTH /////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "InflowFlux"

/**
 * @brief Applies inlet boundary conditions based on the modern BC handling system.
 *
 * This function iterates through all 6 domain faces. For each face identified as an
 * INLET, it applies the velocity profile specified by its assigned handler and
 * parameters (e.g., 'constant_velocity' with vx,vy,vz or 'parabolic' with u_max).
 *
 * It calculates the contravariant flux (Ucont), Cartesian velocity on the face (Ubcs),
 * and the staggered Cartesian velocity (Ucat). It also computes the total incoming
 * flux and area across all MPI ranks.
 *
 * @param user The main UserCtx struct containing the BC configuration and PETSc objects.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InflowFlux(UserCtx *user) 
{
  PetscErrorCode ierr;
  PetscReal    lFluxIn = 0.0, lAreaIn = 0.0, AreaSumIn;
  Vec          lCoor;
  Cmpnts       ***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet, ***cent;  
  PetscReal    ***nvert;
  
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  // Get context for coordinate transformation if needed by a handler
  SimCtx      *simCtx = user->simCtx;
  PetscReal   CMx_c = simCtx->CMx_c;
  PetscReal   CMy_c = simCtx->CMy_c;
  PetscReal   CMz_c = simCtx->CMz_c;
  
  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
  if (xs==0) lxs = xs+1; if (ys==0) lys = ys+1; if (zs==0) lzs = zs+1;
  if (xe==mx) lxe = xe-1; if (ye==my) lye = ye-1; if (ze==mz) lze = ze-1;

  // --- Get PETSc arrays ---
  DMGetCoordinatesLocal(da,&lCoor);
  DMDAVecGetArray(fda, lCoor, &coor); // Use local coordinates
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(da,  user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->Cent, &cent);
  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  // --- Main Loop over all 6 faces ---
  for (PetscInt fn=0; fn<6; fn++) {
    BoundaryFaceConfig *face_config = &user->boundary_faces[fn];

    if (face_config->mathematical_type != INLET) {
        continue; // Skip non-inlet faces
    }
    
    // This processor only acts if it is on the boundary of the global domain
    PetscBool is_on_boundary = ( (fn==0 && xs==0) || (fn==1 && xe==mx) ||
                                 (fn==2 && ys==0) || (fn==3 && ye==my) ||
                                 (fn==4 && zs==0) || (fn==5 && ze==mz) );
    if (!is_on_boundary) continue;

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Applying INLET handler for face: %s \n", BCFaceToString(face_config->face_id));

    // --- Loop over the specific face geometry ---
    switch(face_config->face_id) {
      case BC_FACE_NEG_X: // -Xi
      case BC_FACE_POS_X: // +Xi
      {
        PetscReal sign = (face_config->face_id == BC_FACE_NEG_X) ? 1.0 : -1.0;
        PetscInt i = (face_config->face_id == BC_FACE_NEG_X) ? xs : mx - 2;
        for (PetscInt k=lzs; k<lze; k++) {
          for (PetscInt j=lys; j<lye; j++) {
            if ( (sign > 0 && nvert[k][j][i+1] > 0.1) || (sign < 0 && nvert[k][j][i] > 0.1) ) continue;

            PetscReal uin_this_point = 0.0;
            // --- Determine velocity based on the handler for this point ---
            if (face_config->handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY) {
                PetscBool found;
                ierr = GetBCParamReal(face_config->params, "vx", &uin_this_point, &found); CHKERRQ(ierr);
            } else if(face_config->handler_type== BC_HANDLER_INLET_PARABOLIC){
              PetscBool found;
              PetscReal umax,diameter=1.0;
              ierr = GetBCParamReal(face_config->params,"u_max",&umax,&found); CHKERRQ(ierr);
              ierr = GetBCParamReal(face_config->params,"diameter",&diameter,&found); CHKERRQ(ierr);

              // Radius is in the YZ-plane for an X-face inlet
              PetscReal yc = cent[k][j][i + (sign>0)].y - CMy_c;
              PetscReal zc = cent[k][j][i + (sign>0)].z - CMz_c;
              PetscReal r = sqrt(yc*yc + zc*zc);
              PetscReal r_norm = 2.0 * r / diameter;
              uin_this_point = umax * (1.0 - r_norm * r_norm);
              if (r_norm > 1.0) uin_this_point = 0.0; 
            }
            // Add other X-face handlers like 'else if (handler == ...)' here

            // --- Apply the calculated velocity ---
            PetscReal CellArea = sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
            lAreaIn += CellArea;
            ucont[k][j][i].x  = sign * uin_this_point * CellArea;
            lFluxIn          += ucont[k][j][i].x;
            ubcs[k][j][i + (sign < 0)].x = sign * uin_this_point * csi[k][j][i].x / CellArea;
            ubcs[k][j][i + (sign < 0)].y = sign * uin_this_point * csi[k][j][i].y / CellArea;
            ubcs[k][j][i + (sign < 0)].z = sign * uin_this_point * csi[k][j][i].z / CellArea;
            ucat[k][j][i + (sign > 0)] = ubcs[k][j][i + (sign < 0)];
          }
        }
      } break;

      case BC_FACE_NEG_Y: // -Eta
      case BC_FACE_POS_Y: // +Eta
      {
        PetscReal sign = (face_config->face_id == BC_FACE_NEG_Y) ? 1.0 : -1.0;
        PetscInt j = (face_config->face_id == BC_FACE_NEG_Y) ? ys : my - 2;
        for (PetscInt k=lzs; k<lze; k++) {
          for (PetscInt i=lxs; i<lxe; i++) {
            if ( (sign > 0 && nvert[k][j+1][i] > 0.1) || (sign < 0 && nvert[k][j][i] > 0.1) ) continue;

            PetscReal uin_this_point = 0.0;
            if (face_config->handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY) {
                PetscBool found;
                ierr = GetBCParamReal(face_config->params, "vy", &uin_this_point, &found); CHKERRQ(ierr);
            }else if(face_config->handler_type == BC_HANDLER_INLET_PARABOLIC){
              PetscBool found;
              PetscReal umax,diameter=1.0;
              ierr = GetBCParamReal(face_config->params,"umax",&umax,&found); CHKERRQ(ierr);
              ierr = GetBCParamReal(face_config->params,"diameter",&diameter,&found); CHKERRQ(ierr);
              
              // Radius is in the XZ-plane for a Y-face inlet
              PetscReal xc = cent[k][j + (sign>0)][i].x - CMx_c;
              PetscReal zc = cent[k][j + (sign>0)][i].z - CMz_c;
              PetscReal r = sqrt(xc*xc + zc*zc);
              PetscReal r_norm = 2.0 * r / diameter;
              uin_this_point = umax * (1.0 - r_norm * r_norm);
              if (r_norm > 1.0) uin_this_point = 0.0;

            }

            PetscReal CellArea = sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);
            lAreaIn += CellArea;
            ucont[k][j][i].y  = sign * uin_this_point * CellArea;
            lFluxIn          += ucont[k][j][i].y;
            ubcs[k][j + (sign < 0)][i].x = sign * uin_this_point * eta[k][j][i].x / CellArea;
            ubcs[k][j + (sign < 0)][i].y = sign * uin_this_point * eta[k][j][i].y / CellArea;
            ubcs[k][j + (sign < 0)][i].z = sign * uin_this_point * eta[k][j][i].z / CellArea;
            ucat[k][j + (sign > 0)][i] = ubcs[k][j + (sign < 0)][i];
          }
        }
      } break;

      case BC_FACE_NEG_Z: // -Zeta
      case BC_FACE_POS_Z: // +Zeta
      {
        PetscReal sign = (face_config->face_id == BC_FACE_NEG_Z) ? 1.0 : -1.0;
        PetscInt k = (face_config->face_id == BC_FACE_NEG_Z) ? zs : mz - 2;
        for (PetscInt j=lys; j<lye; j++) {
          for (PetscInt i=lxs; i<lxe; i++) {
            if ( (sign > 0 && nvert[k+1][j][i] > 0.1) || (sign < 0 && nvert[k][j][i] > 0.1) ) continue;

            PetscReal uin_this_point = 0.0;
            if (face_config->handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY) {
                PetscBool found;
                ierr = GetBCParamReal(face_config->params, "vz", &uin_this_point, &found); CHKERRQ(ierr);
            } else if (face_config->handler_type == BC_HANDLER_INLET_PARABOLIC) {
                PetscBool found;
                PetscReal umax, diameter=1.0; // Default diameter for r_norm = 2*r
                ierr = GetBCParamReal(face_config->params, "u_max", &umax, &found); CHKERRQ(ierr);
                ierr = GetBCParamReal(face_config->params, "diameter", &diameter, &found); CHKERRQ(ierr); // Optional
                
                // Radius in the XY-plane for a Z-face inlet.
                PetscReal xc = cent[k + (sign>0)][j][i].x - CMx_c;
                PetscReal yc = cent[k + (sign>0)][j][i].y - CMy_c;
                PetscReal r = sqrt(xc*xc + yc*yc);
                PetscReal r_norm = 2.0 * r / diameter; // Normalized radius
                uin_this_point = umax * (1.0 - r_norm * r_norm);
                if (r_norm > 1.0) uin_this_point = 0.0;
            }

            PetscReal CellArea = sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
            lAreaIn += CellArea;
            ucont[k][j][i].z  = sign * uin_this_point * CellArea;
            lFluxIn          += ucont[k][j][i].z;
            ubcs[k + (sign < 0)][j][i].x = sign * uin_this_point * zet[k][j][i].x / CellArea;
            ubcs[k + (sign < 0)][j][i].y = sign * uin_this_point * zet[k][j][i].y / CellArea;
            ubcs[k + (sign < 0)][j][i].z = sign * uin_this_point * zet[k][j][i].z / CellArea;
            ucat[k + (sign > 0)][j][i] = ubcs[k + (sign < 0)][j][i];
          }
        }
      } break;
    } // end switch(face_id)
  } // end for(fn)

  // --- Finalize: Sum flux and area from all processes ---
  ierr = MPI_Allreduce(&lFluxIn, &simCtx->FluxInSum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = MPI_Allreduce(&lAreaIn, &AreaSumIn, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
  
  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Inflow: Flux - Area:  %le - %le \n", simCtx->FluxInSum, AreaSumIn);    
  
  // --- Restore PETSc arrays ---
  DMDAVecRestoreArray(fda, lCoor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da,  user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  // --- Update local vectors for subsequent computations ---
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  
  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "OutflowFlux"

/**
 * @brief Calculates the total outgoing flux through all OUTLET faces for reporting.
 *
 * NOTE: In a mixed modern/legacy environment, this function is for DIAGNOSTICS ONLY.
 * It reads the contravariant velocities and calculates the total flux passing through
 * faces marked as OUTLET. It does NOT apply any boundary conditions itself, as that
 * is still the responsibility of the legacy monolith function.
 *
 * @param user The main UserCtx struct containing BC config and PETSc vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode OutflowFlux(UserCtx *user) {
  
    PetscErrorCode ierr;
    PetscReal      lFluxOut = 0.0;
    Cmpnts         ***ucont;

    DM             fda = user->fda;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    PetscFunctionBeginUser;
    
    PROFILE_FUNCTION_BEGIN;

    ierr = DMDAVecGetArrayRead(fda, user->Ucont, &ucont); CHKERRQ(ierr);
  
    // --- Loop over all 6 faces to find OUTLETS ---
    for (PetscInt fn = 0; fn < 6; fn++) {
        if (user->boundary_faces[fn].mathematical_type != OUTLET) {
            continue; 
        }

        PetscBool is_on_boundary = ( (fn==0 && xs==0) || (fn==1 && xe==mx) ||
                                     (fn==2 && ys==0) || (fn==3 && ye==my) ||
                                     (fn==4 && zs==0) || (fn==5 && ze==mz) );
        if (!is_on_boundary) continue;

        // --- Sum the flux for the appropriate face and component ---
        switch ((BCFace)fn) {
            case BC_FACE_NEG_X: case BC_FACE_POS_X: {
                PetscInt i = (fn == 0) ? xs : mx - 2;
                for (PetscInt k=info.zs; k<info.zs+info.zm; k++) for (PetscInt j=info.ys; j<info.ys+info.ym; j++) {
                    lFluxOut += ucont[k][j][i].x;
                }
            } break;

            case BC_FACE_NEG_Y: case BC_FACE_POS_Y: {
                PetscInt j = (fn == 2) ? ys : my - 2;
                for (PetscInt k=info.zs; k<info.zs+info.zm; k++) for (PetscInt i=info.xs; i<info.xs+info.xm; i++) {
                    lFluxOut += ucont[k][j][i].y;
                }
            } break;

            case BC_FACE_NEG_Z: case BC_FACE_POS_Z: {
                PetscInt k = (fn == 4) ? zs : mz - 2;
                for (PetscInt j=info.ys; j<info.ys+info.ym; j++) for (PetscInt i=info.xs; i<info.xs+info.xm; i++) {
                    lFluxOut += ucont[k][j][i].z;
                }
            } break;
        } // end switch
    } // end for loop

    ierr = DMDAVecRestoreArrayRead(fda, user->Ucont, &ucont); CHKERRQ(ierr);

    // --- Finalize: Sum and store the global total flux ---
    ierr = MPI_Allreduce(&lFluxOut, &user->simCtx->FluxOutSum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Reported Global FluxOutSum = %.6f\n", user->simCtx->FluxOutSum);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormBCS"

/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation/interface
   -1:  wallfunction
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   3:   slip wall/symmetry
   5:	Inlet
   4:	Outlet
   6:   farfield
   7:   periodic
   8:   Characteristic BC
   9:   Analytical Vortex
   10:  Oulet Junction
   11:  Annulus
   12:  Ogrid
   13:  Rheology
   14:  Outlet with Interface
   15:  No Gradient (Similar to Farfield)  
*/

PetscErrorCode FormBCS(UserCtx *user)
{
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  PetscReal	***nvert,***lnvert; //local working array

  PetscReal	***p,***lp;
  Cmpnts	***ucont, ***ubcs, ***ucat,***lucat, ***csi, ***eta, ***zet;
  Cmpnts	***cent,***centx,***centy,***centz,***coor;
  PetscScalar	FluxIn, FluxOut,ratio;
  PetscScalar   lArea, AreaSum;
 
  PetscScalar   FarFluxIn=0., FarFluxOut=0., FarFluxInSum, FarFluxOutSum;
  PetscScalar   FarAreaIn=0., FarAreaOut=0., FarAreaInSum, FarAreaOutSum;
  PetscScalar   FluxDiff, VelDiffIn, VelDiffOut;
  Cmpnts        V_frame;
 
  PetscReal Un, nx,ny,nz,A;

  SimCtx *simCtx = user->simCtx;  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx ) lxe = xe-1;
  if (ye==my ) lye = ye-1;
  if (ze==mz ) lze = ze-1;

  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  PetscInt ttemp;
  for (ttemp=0; ttemp<3; ttemp++) {
    DMDAVecGetArray(da, user->Nvert, &nvert); 
    DMDAVecGetArray(fda, user->lUcat,  &ucat);
    DMDAVecGetArray(fda, user->Ucont, &ucont);
/* ==================================================================================             */
/*   FAR-FIELD BC */
/* ==================================================================================             */
 

  // reset FAR FLUXES
  FarFluxIn = 0.; FarFluxOut=0.;
  FarAreaIn = 0.; FarAreaOut=0.;

  PetscReal lFlux_abs=0.0,FluxSum_abs=0.0,ratio=0.0;

    V_frame.x=0.;
    V_frame.y=0.;
    V_frame.z=0.;

 

  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x;
	  ubcs[k][j][i].y = ucat[k][j][i+1].y;
	  ubcs[k][j][i].z = ucat[k][j][i+1].z;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;
	  FarFluxIn += ucont[k][j][i].x;
	  lFlux_abs += fabs(ucont[k][j][i].x);
	  FarAreaIn += csi[k][j][i].x;
	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x;
	  ubcs[k][j][i].y = ucat[k][j][i-1].y;
	  ubcs[k][j][i].z = ucat[k][j][i-1].z;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	  FarFluxOut += ucont[k][j][i-1].x;
	  lFlux_abs  += fabs(ucont[k][j][i-1].x);
	  FarAreaOut += csi[k][j][i-1].x;
	}
      }
    }
  }

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j+1][i].x;
	  ubcs[k][j][i].y = ucat[k][j+1][i].y;
	  ubcs[k][j][i].z = ucat[k][j+1][i].z;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;
	  FarFluxIn += ucont[k][j][i].y;
	  lFlux_abs += fabs(ucont[k][j][i].y);
	  FarAreaIn += eta[k][j][i].y;
	}
      }
    }
  }
  
  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j-1][i].x;
	  ubcs[k][j][i].y = ucat[k][j-1][i].y;
	  ubcs[k][j][i].z = ucat[k][j-1][i].z;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
	  FarFluxOut += ucont[k][j-1][i].y;
	  lFlux_abs  += fabs(ucont[k][j-1][i].y);
	  FarAreaOut += eta[k][j-1][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==6 || user->bctype[4]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;
	  FarFluxIn += ucont[k][j][i].z;
	  lFlux_abs += fabs(ucont[k][j][i].z);
	  FarAreaIn += zet[k][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  FarFluxOut += ucont[k-1][j][i].z;
	  lFlux_abs  += fabs(ucont[k-1][j][i].z); 
	  FarAreaOut += zet[k-1][j][i].z;
	}
      }
    }
  }
  
  MPI_Allreduce(&FarFluxIn,&FarFluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarFluxOut,&FarFluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&lFlux_abs,&FluxSum_abs,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); 
  MPI_Allreduce(&FarAreaIn,&FarAreaInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarAreaOut,&FarAreaOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
  if (user->bctype[5]==6 || user->bctype[3]==6 || user->bctype[1]==6) {
  
    ratio=(FarFluxInSum - FarFluxOutSum)/FluxSum_abs;
    if (fabs(FluxSum_abs) <1.e-10) ratio = 0.;
    
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "/FluxSum_abs %le ratio %le \n", FluxSum_abs,ratio);
    
    FluxDiff = 0.5*(FarFluxInSum - FarFluxOutSum) ;
    VelDiffIn  = FluxDiff / FarAreaInSum ;
    
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = FluxDiff / FarAreaOutSum ;
  
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Far Flux Diff %d %le %le %le %le %le %le %le\n", simCtx->step, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  if (user->bctype[5]==10) {
    FluxDiff = simCtx->FluxInSum -( FarFluxOutSum -FarFluxInSum) ;
    VelDiffIn  = 1/3.*FluxDiff / (FarAreaInSum);// +  FarsimCtx->AreaOutSum);
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = 2./3.* FluxDiff / (FarAreaOutSum) ;//(FarAreaInSum +  FarsimCtx->AreaOutSum) ;
   
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Far Flux Diff %d %le %le %le %le %le %le %le\n", simCtx->step, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  
  // scale global mass conservation

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k-1][j][i].z  += ratio*fabs(ucont[k-1][j][i].z);
	  ubcs[k][j][i].z = ucont[k-1][j][i].z/zet[k-1][j][i].z;
	  //  ubcs[k][j][i].z = ucat[k-1][j][i].z + VelDiffOut ;//+ V_frame.z;
	  // ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {

	  ucont[k][j-1][i].y +=ratio*fabs(ucont[k][j-1][i].y);
	  ubcs[k][j][i].y = ucont[k][j-1][i].y /eta[k][j-1][i].y;
	  //	  ubcs[k][j][i].y = ucat[k][j-1][i].y + VelDiffOut;// + V_frame.y;
	  // ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;

	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucont[k][j][i-1].x +=ratio*fabs(ucont[k][j][i-1].x);
	  ubcs[k][j][i].x = ucont[k][j][i-1].x / csi[k][j][i-1].x ;
	  //  ubcs[k][j][i].x = ucat[k][j][i-1].x + VelDiffOut;// + V_frame.x;
	  // ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	}
      }
    }
  }


  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucont[k][j][i].x  -=ratio*fabs(ucont[k][j][i].x);
	  ubcs[k][j][i].x = ucont[k][j][i].x / csi[k][j][i].x;
	  // ubcs[k][j][i].x = ucat[k][j][i+1].x - VelDiffIn;// + V_frame.x;
	  // ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;
	}
      }
    }
  }
  

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].y -=ratio*fabs(ucont[k][j][i].y);
	  ubcs[k][j][i].y = ucont[k][j][i].y / eta[k][j][i].y;
	  //	  ubcs[k][j][i].y = ucat[k][j+1][i].y - VelDiffIn;// + V_frame.y;
	  // ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;
	}
      }
    }
  }
  
  
  if (user->bctype[4]==6 || user->bctype[5]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].z -=ratio*fabs(ucont[k][j][i].z);
	  ubcs[k][j][i].z =ucont[k][j][i].z / zet[k][j][i].z;
	  // ubcs[k][j][i].z = ucat[k+1][j][i].z - VelDiffIn;// + V_frame.z;
	  // ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;

	}
      }
    }
  }

//// Amir wall Ogrid
 
if (user->bctype[2]==1 || user->bctype[2]==-1)  {
    if (ys==0) {
    j= ys;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
        A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
               eta[k][j][i].y*eta[k][j][i].y +
               eta[k][j][i].x*eta[k][j][i].x);
        nx=eta[k][j][i].x/A;
        ny=eta[k][j][i].y/A;
        nz=eta[k][j][i].z/A;
        Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
        ubcs[k][j][i].x = 0.0;
        ubcs[k][j][i].y = 0.0;
        ubcs[k][j][i].z = 0.0;
        ucont[k][j][i].y = 0.;
      }
    }
    }
 }

/* ==================================================================================             */
/*   SOLID WALL BC (NO-SLIP / NO-PENETRATION) */
/* ==================================================================================             */

// NOTE: This block is added to explicitly handle bctype=1 (solid wall) for all faces.
// It ensures both no-slip (ubcs=0) and no-penetration (ucont_normal=0).
// ubcs is handled by the implicit-zero assumption, but ucont must be set explicitly.

// -X Face (i=0)
if (user->bctype[0]==1 || user->bctype[0]==-1)  {
    if (xs==0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
        for (j=lys; j<lye; j++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          ucont[k][j][i].x = 0.0; // Enforce no-penetration
        }
      }
    }
}

// +X Face (i=mx-1)
if (user->bctype[1]==1 || user->bctype[1]==-1)  {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
        for (j=lys; j<lye; j++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          // The relevant ucont is at the face, index i-1
          ucont[k][j][i-1].x = 0.0; // Enforce no-penetration
        }
      }
    }
}

// -Y Face (j=0)
if (user->bctype[2]==1 || user->bctype[2]==-1)  {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
        for (i=lxs; i<lxe; i++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          ucont[k][j][i].y = 0.0; // Enforce no-penetration
        }
      }
    }
}

// +Y Face (j=my-1)
if (user->bctype[3]==1 || user->bctype[3]==-1)  {
    if (ye==my) {
      j= ye-1;
      for (k=lzs; k<lze; k++) {
        for (i=lxs; i<lxe; i++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          // The relevant ucont is at the face, index j-1
          ucont[k][j-1][i].y = 0.0; // Enforce no-penetration
        }
      }
    }
}

/* Original "Amir wall Ogrid" block can now be removed or commented out
   as its functionality is included above.
if (user->bctype[2]==1 || user->bctype[2]==-1)  { ... }
*/
 
/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i].z*csi[k][j][i].z +
	       csi[k][j][i].y*csi[k][j][i].y +
	       csi[k][j][i].x*csi[k][j][i].x);
	nx=csi[k][j][i].x/A;
	ny=csi[k][j][i].y/A;
	nz=csi[k][j][i].z/A;
	Un=ucat[k][j][i+1].x*nx+ucat[k][j][i+1].y*ny+ucat[k][j][i+1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i+1].x-Un*nx;//-V_frame.x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i+1].z-Un*nz;
	ucont[k][j][i].x = 0.;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i-1].z*csi[k][j][i-1].z +
	       csi[k][j][i-1].y*csi[k][j][i-1].y +
	       csi[k][j][i-1].x*csi[k][j][i-1].x);
	nx=csi[k][j][i-1].x/A;
	ny=csi[k][j][i-1].y/A;
	nz=csi[k][j][i-1].z/A;
	Un=ucat[k][j][i-1].x*nx+ucat[k][j][i-1].y*ny+ucat[k][j][i-1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i-1].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i-1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i-1].z-Un*nz;
	ucont[k][j][i-1].x = 0.;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
	       eta[k][j][i].y*eta[k][j][i].y +
	       eta[k][j][i].x*eta[k][j][i].x);
	nx=eta[k][j][i].x/A;
	ny=eta[k][j][i].y/A;
	nz=eta[k][j][i].z/A;
	Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j+1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j+1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j+1][i].z-Un*nz;
	ucont[k][j][i].y = 0.;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j-1][i].z*eta[k][j-1][i].z +
	       eta[k][j-1][i].y*eta[k][j-1][i].y +
	       eta[k][j-1][i].x*eta[k][j-1][i].x);
	nx=eta[k][j-1][i].x/A;
	ny=eta[k][j-1][i].y/A;
	nz=eta[k][j-1][i].z/A;
	Un=ucat[k][j-1][i].x*nx+ucat[k][j-1][i].y*ny+ucat[k][j-1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j-1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j-1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j-1][i].z-Un*nz;
	ucont[k][j-1][i].y = 0.;
      }
    }
    }
  }
  

  if (user->bctype[4]==3) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k][j][i].z*zet[k][j][i].z +
	       zet[k][j][i].y*zet[k][j][i].y +
	       zet[k][j][i].x*zet[k][j][i].x);
	nx=zet[k][j][i].x/A;
	ny=zet[k][j][i].y/A;
	nz=zet[k][j][i].z/A;
	Un=ucat[k+1][j][i].x*nx+ucat[k+1][j][i].y*ny+ucat[k+1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k+1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k+1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k+1][j][i].z-Un*nz;
	ucont[k][j][i].z = 0.;
	}
      }
    }
  }

  if (user->bctype[5]==3) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z +
	       zet[k-1][j][i].y*zet[k-1][j][i].y +
	       zet[k-1][j][i].x*zet[k-1][j][i].x);
	nx=zet[k-1][j][i].x/A;
	ny=zet[k-1][j][i].y/A;
	nz=zet[k-1][j][i].z/A;
	Un=ucat[k-1][j][i].x*nx+ucat[k-1][j][i].y*ny+ucat[k-1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k-1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k-1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k-1][j][i].z-Un*nz;
	ucont[k-1][j][i].z = 0.;
	}
      }
    }
  }
 
/* ==================================================================================             */
/*     CHARACTERISTIC OUTLET BC :8 */
/* ==================================================================================             */

  if (user->bctype[5]==8) {
    if (ze == mz) {
      k = ze-2;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = simCtx->FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    // PetscGlobalSum(PETSC_COMM_WORLD,&FluxOut, &FluxOutSum);

    //ratio = FluxInSum / FluxOutSum;
    ratio = FluxIn / simCtx->FluxOutSum;
    if (fabs(simCtx->FluxOutSum) < 1.e-6) ratio = 1.;
    //if (fabs(FluxInSum) <1.e-6) ratio = 0.;
    if (fabs(FluxIn) <1.e-6) ratio = 0.;
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Char Ratio %d %le %le %le %le %d %d\n", simCtx->step, ratio, FluxIn, simCtx->FluxOutSum, FarFluxInSum,zs, ze);

    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  if (simCtx->step==0 || simCtx->step==1)
	    if (simCtx->inletprofile<0)
	      ubcs[k][j][i].z = -1.;
	    else if (user->bctype[4]==6)
	      ubcs[k][j][i].z = 0.;
	    else
	      ubcs[k][j][i].z = 1.;//ubcs[0][j][i].z;//-1.;//1.;
	  
	  else
	    ucont[k-1][j][i].z = ucont[k-1][j][i].z*ratio;
	  
	  ubcs[k][j][i].z = ucont[k-1][j][i].z / zet[k-1][j][i].z;
	}
      }
    }
  }
 
/* ==================================================================================             */
/*     OUTLETBC :4 */
/* ==================================================================================             */

  
  if (user->bctype[5]==OUTLET || user->bctype[5]==14 || user->bctype[5]==20) {
    lArea=0.;
    LOG_ALLOW(GLOBAL,LOG_VERBOSE,"+Zeta Outlet \n");
     LOG_ALLOW(GLOBAL,LOG_VERBOSE,"FluxOutSum before FormBCS applied = %.6f \n",simCtx->FluxOutSum);
    if (ze == mz) {
      //    k = ze-3;
      k=ze-1;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	  if ((nvert[k-1][j][i])<0.1) {
	  FluxOut += (ucat[k-1][j][i].x * (zet[k-1][j][i].x) +
		      ucat[k-1][j][i].y * (zet[k-1][j][i].y) +
		      ucat[k-1][j][i].z * (zet[k-1][j][i].z));

	  lArea += sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
			 (zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
			 (zet[k-1][j][i].z) * (zet[k-1][j][i].z));
	 
	  }
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   
    user->simCtx->AreaOutSum = AreaSum;

    LOG_ALLOW(GLOBAL,LOG_VERBOSE,"AreaOutSum = %.6f | FluxOutSum = %.6f \n",AreaSum,simCtx->FluxOutSum);
    
    if (simCtx->block_number>1 && user->bctype[5]==14) {
      simCtx->FluxOutSum += user->FluxIntfcSum;
      //      AreaSum    += user->AreaIntfcSum;
    }

    FluxIn = simCtx->FluxInSum + FarFluxInSum + user->FluxIntpSum;
    if (user->bctype[5]==20)
      ratio = (FluxIn / simCtx->FluxOutSum);
    else
      ratio = (FluxIn - simCtx->FluxOutSum) / AreaSum;
   
     LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Ratio for momentum correction = %.6f \n",ratio);
    
  /*   user->FluxOutSum += ratio*user->simCtx->AreaOutSum; */
    simCtx->FluxOutSum =0.0;
    FluxOut=0.0;
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if ((nvert[k-1][j][i])<0.1) {
	 
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;//+ratio;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z;// + ratio;//*n_z;

	    //  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	    if (user->bctype[5]==20)
	      ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				    ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				    ubcs[k][j][i].z * (zet[k-1][j][i].z ))*ratio;
	    
	    else{
	      ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				    ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				    ubcs[k][j][i].z * (zet[k-1][j][i].z ))
		+ ratio * sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
				(zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
				(zet[k-1][j][i].z) * (zet[k-1][j][i].z)); 

	      FluxOut += ucont[k-1][j][i].z;
	    }
	  }//if
	}
      }
    }
    
    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    LOG_ALLOW(GLOBAL,LOG_TRACE, "Timestep = %d | FluxInSum = %.6f | FlucOutSum = %.6f | FluxIntfcSum = %.6f | FluxIntpSum = %.6f \n", simCtx->step, simCtx->FluxInSum, simCtx->FluxOutSum, user->FluxIntfcSum,user->FluxIntpSum);

  } else if (user->bctype[5]==2) {
  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1;//sin(2*3.14*simCtx->step*simCtx->dt);//1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
  

  /*   OUTLET at k==0 */
  if (user->bctype[4]==OUTLET) {
    lArea=0.;
    if (zs == 0) {
      k = zs;
      //      k= zs + 1;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {


	  FluxOut += ucat[k+1][j][i].z * zet[k][j][i].z ;

	  lArea += zet[k][j][i].z;



	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = simCtx->FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  

    ratio = (simCtx->FluxInSum - simCtx->FluxOutSum) / AreaSum;
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Ratio b %d  %le %le %le %le %d %d\n", simCtx->step,ratio, simCtx->FluxInSum, simCtx->FluxOutSum, AreaSum,zs, ze);
    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
          ubcs[k][j][i].x = ucat[k+1][j][i].x;
          ubcs[k][j][i].y = ucat[k+1][j][i].y;
          ubcs[k][j][i].z = ucat[k+1][j][i].z;
          ucont[k][j][i].z = (ubcs[k][j][i].z+ratio) * zet[k][j][i].z;
	      }
      }
    }
  }


  
/* ==================================================================================             */
/*     Ogrid :77 */
/* ==================================================================================             */
  /* 
  if (user->bctype[3]=77 && Ogrid)
    {Cmpnts ***coor;
      lArea=0.;
      FluxOut=0.0;
      //    k = ze-3;
      
      Vec Coor; DMGetCoordinatesLocal(da, &Coor); 
      DMDAVecGetArray(fda,Coor,&coor);       
      if (ye==my) {
	j=my-2;
	for (k=lzs; k<lze; k++){
	  for (i=lxs; i<lxe; i++) {
	    
	      FluxOut += ucont[k][j][i].y;
	      lArea += sqrt( (eta[k-1][j][i].x) * (eta[k-1][j][i].x) +
			     (eta[k-1][j][i].y) * (eta[k-1][j][i].y) +
			     (eta[k-1][j][i].z) * (eta[k-1][j][i].z));
	 
	  }
	}
      }
      
      else {
	FluxOut = 0.;
      }
      
      MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      
      user->FluxOutSum = FluxOutSum;
      simCtx->AreaOutSum = AreaSum;
     export PL="$HOME/CSBL/Codes/fdf_project/pic_les"
 ratio=2*(FluxIn-FluxOutSum)/AreaSum;
       ratio=0.0;   
      if (ye==my){
	j=my-2;
	for (k=lzs; k<lze; k++){      
	  for (i=lxs; i<lxe; i++) {	
	    if ((nvert[k-1][j][i])<0.1 && coor[k][j][i].z >= 800.) {
	      ucont[k][j][i].y *= (1+ratio);
	      
	    }
	  }
	}
      }
            
      if (ye==my){
	j=my-2;
	for (k=lzs; k<lze; k++){      
	  for (i=lxs; i<lxe; i++) {	
	    
	    ubcs[k][j+1][i].z=ucat[k][j][i].z;
	    ubcs[k][j+1][i].y=ucat[k][j][i].y;
	    ubcs[k][j+1][i].x=ucat[k][j][i].x;
	  }
	}
      }
      
      //   LOG_ALLOW(GLOBAL,LOG_DEBUG, "  ratio %le ",ratio);

    }

  */


/* ==================================================================================             */
/*     Channelz */
/* ==================================================================================             */
 // Amir channel flow correction
  if (user->bctype[4]==7 && simCtx->channelz==1) {
 Vec Coor; DMGetCoordinatesLocal(da, &Coor); 
 DMDAVecGetArray(fda,Coor,&coor); 
    Cmpnts  ***uch;
    DMDAVecGetArray(fda, user->Bcs.Uch, &uch);

    lArea=0.0;
   // if (zs==0) {
   // k=0;
    FluxIn=0.0;
    
    double Fluxbcs=0.0, Fluxbcssum,ratiobcs;

   if (zs==0) {
    k=0;
     Fluxbcs=0.0;      
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i]<0.1){
	    Fluxbcs += ucont[k][j][i].z;

	    lArea +=  sqrt((zet[k][j][i].x) * (zet[k][j][i].x) +
			  (zet[k][j][i].y) * (zet[k][j][i].y) +
			  (zet[k][j][i].z) * (zet[k][j][i].z));
	  }
	}
      }
    }
    

    //int kk=(simCtx->step % (mz-2))+2;
 
    for (k=zs;k<lze;k++){      
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i]<0.1){
	    FluxIn += ucont[k][j][i].z /((mz)-1);
	  }
	}
      }
    }
   // else {
  //   FluxIn=0.0;
//  }
 
    MPI_Allreduce(&FluxIn,&simCtx->Fluxsum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&Fluxbcs,&Fluxbcssum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
    if (simCtx->step==simCtx->StartStep && simCtx->StartStep > 0 && simCtx->ccc==0) {
      simCtx->ccc=1;
      simCtx->FluxInSum=Fluxbcssum;
//	simCtx->FluxInSum=6.3908; 
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  FluxInSum %le .\n ",simCtx->FluxInSum);
    }
  
	simCtx->FluxInSum=6.35066; 
        ratio=(simCtx->FluxInSum-simCtx->Fluxsum)/AreaSum;
	simCtx->ratio=ratio;
        ratiobcs=(simCtx->FluxInSum-Fluxbcssum)/AreaSum;

    if (zs==0) {
	k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k+1][j][i]<0.1){
	    if (simCtx->fish) { 
			 ucont[k][j][i].z+=ratiobcs * /* (1.-exp(-500. * (1.-fabs(coor[k][j][i].y))))  */ sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) + 
 					   (zet[k+1][j][i].y) * (zet[k+1][j][i].y) + 
 					   (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); 
			}

	    uch[k][j][i].z=ratiobcs * /* (1.-exp(-500. * (1.-fabs(coor[k][j][i].y)))) */   sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) +
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
 	    if (simCtx->fish){
		   ucont[k-1][j][i].z+=ratiobcs * /*(1.-exp(-500. * (1.-fabs(coor[k][j][i].y))))  */    sqrt((zet[k-1][j][i].x) * (zet[k-1][j][i].x) + 
 					   (zet[k-1][j][i].y) * (zet[k-1][j][i].y) + 
 					   (zet[k-1][j][i].z) * (zet[k-1][j][i].z));  	
	  }
	    uch[k][j][i].z=ratiobcs *   /*(1.-exp(-500. * (1.-fabs(coor[k][j][i].y)))) */     sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) +
					   (zet[k+1][j][i].y) * (zet[k+1][j][i].y) +
					   (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); 
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->Bcs.Uch, &uch);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Ratio  %le %.15le %.15le  %.15le \n", ratio, ratiobcs, simCtx->FluxInSum,AreaSum);
    DMDAVecRestoreArray(fda,Coor,&coor); 

  ////.................////
   
   
  

    //just for check
/*         if (zs==0) { */
/*       k=10; */
/*       FluxIn=0.0; */
/*       for (j=lys; j<lye; j++) { */
/* 	for (i=lxs; i<lxe; i++) { */
/* 	  if (nvert[k+1][j][i]<0.1){ */
/* 	    FluxIn += ucont[k][j][i].z; */
/* 	    lArea += sqrt((zet[k+1][j][i].x) * (zet[k+1][j][i].x) + */
/* 			  (zet[k+1][j][i].y) * (zet[k+1][j][i].y) + */
/* 			  (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); */
/* 	  } */
/* 	} */
/*       } */
/* 	} */
/*     else { */
/*       FluxIn=0.0; */
/*     } */
/*     // LOG_ALLOW(PETSC_COMM_SELF, "  Fluxsum %le ",FluxIn); */

/*     MPI_Allreduce(&FluxIn,&Fluxsum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */



  }//channel



  /*  ==================================================================================== */
  /*     Cylinder O-grid */
  /*  ==================================================================================== */
  if (user->bctype[3]==12) {
  /* Designed to test O-grid for flow over a cylinder at jmax velocity is 1 (horizontal) 
   u_x = 1 at k==kmax */
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 1.;
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Annulus */
  /*  ==================================================================================== */
  /* designed to test periodic boundary condition for O-grid j=0 rotates */
  DMDAVecGetArray(fda, user->Cent, &cent);
  if (user->bctype[2]==11) {
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	
	/*   ubcs[k][j][i].x=0.0; */
	 
/* 	  ubcs[k][j][i].y = -cent[k][j+1][i].z/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y); */
	 
/* 	  ubcs[k][j][i].z =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y); */
	  ubcs[k][j][i].x = cent[k][j+1][i].y/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);;
	  ubcs[k][j][i].y =-cent[k][j+1][i].x/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  ubcs[k][j][i].z =0.0;
	  //  if(k==1)  LOG_ALLOW(PETSC_COMM_SELF, "@ i= %d j=%d k=%d ubcs.y is %le\n",i,j,k,ubcs[k][j][i].y);
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Rheology */
  /*  ==================================================================================== */
 
  if(simCtx->rheology && (user->bctype[2]==13 || user->bctype[3]==13 || user->bctype[4]==13 || user->bctype[5]==13)){
      LOG_ALLOW(GLOBAL,LOG_DEBUG, "moving plate velocity for rheology setup is %le \n",simCtx->U_bc);
  }
  
  if (user->bctype[2]==13){
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	   ubcs[k][j][i].x = 0.;
	  // ubcs[k][j][i].x = -simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = -simCtx->U_bc;
	  //ubcs[k][j][i].z =0.0;
	}
      }
    }
  }
  if (user->bctype[3]==13){
    if (ye==my){
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  // ubcs[k][j][i].x = simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = simCtx->U_bc;
	   //ubcs[k][j][i].z =0.0;
	}
      }
    }
  }

 if (user->bctype[4]==13){
    if (zs==0){
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =-simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  if (user->bctype[5]==13){
    if (ze==mz){
      k=ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

 
 
  Contra2Cart(user);
  UpdateLocalGhosts(user,"Ucat");
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */

  if (simCtx->wallfunction && user->bctype[2]==-1) {
  PetscReal ***friction_Velocity, ***aj;
  //Mohsen Dec 2015
  Vec Aj  =  user->lAj;
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(da, Aj,  &aj);
//  DMDAVecGetArray(da, user->Nvert, &nvert);
	DMDAVecGetArray(da, user->lFriction_Velocity, &friction_Velocity);
 
  // wall function for boundary
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {

	if( nvert[k][j][i]<1.1 &&  user->bctype[2]==-1 && j==1 )
	 {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc;
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  
	  
	  //Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  //if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	
         double AA=sqrt(eta[k][j][i].z*eta[k][j][i].z +
               eta[k][j][i].y*eta[k][j][i].y +
               eta[k][j][i].x*eta[k][j][i].x);
 	nj[0]=eta[k][j][i].x/AA;
        nj[1]=eta[k][j][i].y/AA;
        nj[2]=eta[k][j][i].z/AA;     
	  noslip (user, sc, sb, Ua, Uc, &ucat[k][j][i], nj[0], nj[1], nj[2]);
	wall_function_loglaw(user, 1.e-16, sc, sb, Ua, Uc, &ucat[k][j][i], &friction_Velocity[k][j][i], nj[0], nj[1], nj[2]);

	 // nvert[k][j][i]=1.;	/* set nvert to 1 to exclude from rhs */
	// if (k==1) 
	  // LOG_ALLOW(GLOBAL,LOG_DEBUG, " %d   %le   %le  %le   %le   %le   %le   %le   %le   %le\n",i, sb,aj[k][j][i],AA, nj[0], nj[1], nj[2],ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z);

	}
      }
  if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
        ubcs[k][j][i].x = 0.0;
        ubcs[k][j][i].y = 0.0;
        ubcs[k][j][i].z = 0.0;
        ucont[k][j][i].y = 0.;
      }
    }
    }

  DMDAVecRestoreArray(da, Aj,  &aj);
  DMDAVecRestoreArray(da, user->lFriction_Velocity, &friction_Velocity);
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
 // DMDAVecRestoreArray(da, user->Nvert, &nvert);
 // DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
 // DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
//   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); 
//   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); 

 
  }

/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Ucat, &ucat);
 
/* ==================================================================================             */
 
  // boundary conditions on ghost nodes
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }
 
  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
 /* ==================================================================================             */
  /*   Periodic BC *///Mohsen
/* ==================================================================================             */
  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    /* /\*   DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert); *\/ */
    /* /\*   DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert); *\/ */
  //Mohsen Dec 2015
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
   
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(j>0 && k>0 && j<user->JM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j][i-2];
	      p[k][j][i]=lp[k][j][i-2];
	      nvert[k][j][i]=lnvert[k][j][i-2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && k>0 && i<user->IM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j-2][i];
	      p[k][j][i]=lp[k][j-2][i];
	      nvert[k][j][i]=lnvert[k][j-2][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && j>0 && i<user->IM && j<user->JM){
	      ucat[k][j][i]=lucat[k-2][j][i];
	      nvert[k][j][i]=lnvert[k-2][j][i];
	      //amir
		p[k][j][i]=lp[k-2][j][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(j>0 && k>0 && j<user->JM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j][i+2];
	      p[k][j][i]=lp[k][j][i+2];
	      nvert[k][j][i]=lnvert[k][j][i+2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && k>0 && i<user->IM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j+2][i];
	      p[k][j][i]=lp[k][j+2][i];
	      nvert[k][j][i]=lnvert[k][j+2][i];
	    }
	  }
	}
      }
    }
  
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && j>0 && i<user->IM && j<user->JM){
	      ucat[k][j][i]=lucat[k+2][j][i];
      	      nvert[k][j][i]=lnvert[k+2][j][i]; 
	      //amir
		p[k][j][i]=lp[k+2][j][i];
	    }
	  }
	}
      }
    }

       


 
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);

  /*  /\*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); *\/ */
  /* /\*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); *\/ */

  /* /\*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); *\/ */
  /* /\*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); *\/ */

  /* /\*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); *\/ */
  /* /\*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); *\/ */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
}
 // 0 velocity on the corner point

  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->P, &p);
  
  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k+1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k][j+1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j+1][i]+p[k][j][i-1]);
      }
    }
  }
 
  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i-1]);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da, user->P, &p);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  //Mohsen Nov 2012
  //Velocity and Presurre at corners for Periodic BC's

  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){
  //i-direction

    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[0]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i-2];
	    p[k][j][i]=lp[k][j][i-2];
	    nvert[k][j][i]=lnvert[k][j][i-2];
	  }
	}
      }
    }
    if (user->bctype[1]==7){
      if (xe==mx){
	i=xe-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i].x=lucat[k][j][i+2].x;
	    p[k][j][i]=lp[k][j][i+2];
	    nvert[k][j][i]=lnvert[k][j][i+2];
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */
  
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    
    //j-direction
    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[2]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j-2][i];
	    p[k][j][i]=lp[k][j-2][i];
	    nvert[k][j][i]=lnvert[k][j-2][i];
	  }
	}
      }
    }
    
    if (user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j+2][i];
	    p[k][j][i]=lp[k][j+2][i];
	    nvert[k][j][i]=lnvert[k][j+2][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
    
/*   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    
    //k-direction
    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[4]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k-2][j][i];
	    nvert[k][j][i]=lnvert[k-2][j][i]; 
	    //amir   
	      p[k][j][i]=lp[k-2][j][i];
	  }
	}
      }
    }
    if (user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k+2][j][i];
	    nvert[k][j][i]=lnvert[k+2][j][i];
	      p[k][j][i]=lp[k+2][j][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
    
/*   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  }
  /* ==================================================================================             */
/*   Analytical Vortex BC */
/* ==================================================================================             */
 
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Cent, &cent); 
  DMDAVecGetArray(fda, user->Centx, &centx);
  DMDAVecGetArray(fda, user->Centy, &centy);
  DMDAVecGetArray(fda, user->Centz, &centz);

  if (user->bctype[0]==9) {
    if (xs==0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i+1].x)*sin(cent[k][j][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i+1].x)*cos(cent[k][j][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	 
	  ucont[k][j][i].x =-(cos(centx[k][j][i].x)*sin(centx[k][j][i].y)*csi[k][j][i].x)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i+1].x)*sin(cent[k][j+1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i+1].x)*cos(cent[k][j+1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i+1].x)*sin(cent[k][j-1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i+1].x)*cos(cent[k][j-1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }
  if (user->bctype[1]==9) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i-1].x)*sin(cent[k][j][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i-1].x)*cos(cent[k][j][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j][i-1].x =(-cos(centx[k][j][i-1].x)*sin(centx[k][j][i-1].y)*csi[k][j][i-1].x)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i-1].x)*sin(cent[k][j+1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i-1].x)*cos(cent[k][j+1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i-1].x)*sin(cent[k][j-1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i-1].x)*cos(cent[k][j-1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }

  if (user->bctype[2]==9) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i].x)*sin(cent[k][j+1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=sin(cent[k][j+1][i].x)*cos(cent[k][j+1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;

	  ucont[k][j][i].y=(sin(centy[k][j][i].x)*cos(centy[k][j][i].y)*eta[k][j][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
    }
  }
 
 
  if (user->bctype[3]==9) {
    if (ye==my) {
      j= ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i].x)*sin(cent[k][j-1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=sin(cent[k][j-1][i].x)*cos(cent[k][j-1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j-1][i].y=(sin(centy[k][j-1][i].x)*cos(centy[k][j-1][i].y)*eta[k][j-1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
    }
  }
  if (user->bctype[4]==9) {
    if (zs==0) {
      k= zs;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k+1][j][i].x;
	  ucat[k][j][i].y=ucat[k+1][j][i].y;
	  ucat[k][j][i].z=ucat[k+1][j][i].z;

	  ucont[k][j][i].z=0.0;
	}
      }
    }
  }
  if (user->bctype[5]==9) {
    if (ze==mz) {
      k= ze-1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k-1][j][i].x;
	  ucat[k][j][i].y=ucat[k-1][j][i].y;
	  ucat[k][j][i].z=ucat[k-1][j][i].z;

	  ucont[k-1][j][i].z=0.0;
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->Centx, &centx);
  DMDAVecRestoreArray(fda, user->Centy, &centy);
  DMDAVecRestoreArray(fda, user->Centz, &centz);

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(fda, user->Ucont,  &ucont);
 
  } // ttemp

 
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); 

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
  }

  /**
 * @brief Applies inlet boundary conditions based on the modern BC handling system.
 *
 * This function iterates through all 6 domain faces. For each face identified as an
 * INLET, it applies the velocity profile specified by its assigned handler and
 * parameters (e.g., 'constant_velocity' with vx,vy,vz or 'parabolic' with u_max).
 *
 * It calculates the contravariant flux (Ucont), Cartesian velocity on the face (Ubcs),
 * and the staggered Cartesian velocity (Ucat). It also computes the total incoming
 * flux and area across all MPI ranks.
 *
 * @param user The main UserCtx struct containing the BC configuration and PETSc objects.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InflowFlux(UserCtx *user); 

/**
 * @brief Calculates the total outgoing flux through all OUTLET faces for reporting.
 *
 * NOTE: In a mixed modern/legacy environment, this function is for DIAGNOSTICS ONLY.
 * It reads the contravariant velocities and calculates the total flux passing through
 * faces marked as OUTLET. It does NOT apply any boundary conditions itself, as that
 * is still the responsibility of the legacy FormBCS function.
 *
 * @param user The main UserCtx struct containing BC config and PETSc vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode OutflowFlux(UserCtx *user);
s
// TO BE FIXED
PetscErrorCode FormBCS(UserCtx *user);