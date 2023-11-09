/*
// Alireza Rashti
// August 2023
*/

/* various functions needed to make a physics object for BH-BH */


#include "bhbh_initialize.h"

/* decide how to initialize the new physics */
Physics_T *bhbh_initialize_new_physics(Physics_T *const old_phys)
{
  Physics_T *new_phys = 0;
  
  /* if already hit the stop */
  if (Pgeti(P_"STOP") == 1)
  {
    /* update things */
    new_phys = infer_new_physics(old_phys);
    /* save and dump */
    write_checkpoint(new_phys,Pgets(P_"my_directory"));
    free_physics(new_phys);
    
    printf(Pretty0"I'm done!  :)\n");
    return 0;
  }
  
  if (!old_phys)/* if empty, come up with a start off */
  {
    /* if we wanna use a particular checkpoint file, 
    // mostly for debug purposes, set:
    // set:
    // BHBH_start_off       = checkpoint_file
    // checkpoint_file_path = path_to_checkpoint_file */
    if (Pcmps(P_"start_off","checkpoint_file"))
    {
      /* modify output directories.
      // AD-HOC: put everything in top directory. */
      Psets(CHECKPOINT_SET_PARAM_ "top_directory", 
            Pgets("top_directory"));
      Psets(CHECKPOINT_SET_PARAM_ P_"my_directory", 
            Pgets("top_directory"));
      Psets(CHECKPOINT_SET_PARAM_ P_"Diagnostics", 
            Pgets("top_directory"));
      
      new_phys = bhbh_read_physics_from_checkpoint();
    }
    
    /* can we resume from a useful checkpoint file */
    else if (can_we_use_checkpoint(Pgets("top_directory")))
      new_phys = bhbh_read_physics_from_checkpoint();
      
    else 
      new_phys = guess_new_physics();
  }
  else/* use old physics, tune it and make new physics */
  {
    new_phys = infer_new_physics(old_phys);
  }
  
  return new_phys;
}


/* use old physics to infer the new physics */
static Physics_T *infer_new_physics(Physics_T *const old_bhbh)
{
  if (!old_bhbh) return 0;
  
  FUNC_TIC
  
  Physics_T *const bhbh    = init_physics(0,BHBH);/* the whole system */
  Physics_T *const bh1     = init_physics(bhbh,BH1);/* BH1 part */
  Physics_T *const bh2     = init_physics(bhbh,BH2);/* BH2 part */
  Physics_T *const old_bh1 = init_physics(old_bhbh,BH1);/* BH1 part */
  Physics_T *const old_bh2 = init_physics(old_bhbh,BH2);/* BH2 part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  old_bh1->grid_char = grid_char;
  old_bh1->igc       = Ibh1;
  old_bh2->grid_char = grid_char;
  old_bh2->igc       = Ibh2;
  
  /* update, adjust and tune */
  physics(old_bh1,BH_TUNE_SPIN);
  physics(old_bh2,BH_TUNE_SPIN);

  physics(old_bh1,BH_TUNE_RADIUS);
  physics(old_bh2,BH_TUNE_RADIUS);
  
  physics(old_bh1,BH_FIND_SURFACE);
  physics(old_bh2,BH_FIND_SURFACE);

  physics(old_bhbh,SYS_TUNE_P_ADM);
  
  /* new grid */
  create_new_grid(grid_char,bhbh);
  bh1->grid = bhbh->grid;
  bh2->grid = bhbh->grid;
  
  /* set and update parameters */
  update_params(bhbh);
  physics(bhbh,FREE_DATA_SET_PARAMS);
  physics(bhbh,ADM_SET_PARAMS);
  physics(bhbh,SYS_SET_PARAMS);
  physics(bhbh,STRESS_ENERGY_SET_PARAMS);
  physics(bhbh,BH_SET_PARAMS);
  physics(bhbh,OBSERVE_SET_PARAMS);  
  
  /* add fields */
  physics(bhbh,ADM_ADD_FIELDS);
  physics(bhbh,FREE_DATA_ADD_FIELDS);
  physics(bhbh,STRESS_ENERGY_ADD_FIELDS);
  physics(bhbh,SYS_ADD_FIELDS);
  physics(bhbh,OBSERVE_ADD_FIELDS);
  physics(bhbh,BH_ADD_FIELDS);
  
  /* populate fields */
  physics(bhbh,FREE_DATA_POPULATE);
  initialize_fields_using_previous_solve(bhbh,old_bhbh);
  
  /* move Jacobian if possible */
  move_jacobian(bhbh,old_bhbh);
  
  /* beta = B0+B1 */
  physics(bhbh,ADM_UPDATE_B1I);
  update_partial_derivatives(bhbh,".*","^dB0_U.+,^ddB0_U.+");
  physics(bhbh,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(bhbh,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  
  /* update AConf^{ij} */
  physics(bhbh,ADM_UPDATE_AConfIJ);

  /* update normal on AH */
  physics(bh1,BH_UPDATE_sConf);
  physics(bh2,BH_UPDATE_sConf);
  
  /* free */
  free_physics(bh1);
  free_physics(bh2);
  free_physics(old_bh1);
  free_physics(old_bh2);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return bhbh;
}  

/* use a known BH1 and BH2 solution to initialize the physics */
static Physics_T *guess_new_physics(void)
{
  FUNC_TIC
  
  Physics_T *const bhbh = init_physics(0,BHBH);/* the whole system */
  Physics_T *const bh1  = init_physics(bhbh,BH1);/* BH1 part */
  Physics_T *const bh2  = init_physics(bhbh,BH2);/* BH2 part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  
  /* set parameters */
  physics(bhbh,FREE_DATA_SET_PARAMS);
  physics(bhbh,ADM_SET_PARAMS);
  physics(bhbh,SYS_SET_PARAMS);
  physics(bhbh,STRESS_ENERGY_SET_PARAMS);
  physics(bhbh,OBSERVE_SET_PARAMS);
  physics(bhbh,BH_SET_PARAMS);
  
  /* create grid */
  bh1->grid_char = grid_char;
  bh1->igc       = Ibh1;
  bh2->grid_char = grid_char;
  bh2->igc       = Ibh2;

  physics(bh1,BH_START);
  physics(bh2,BH_START);
  
  physics(bh1,BH_FIND_SURFACE);
  physics(bh2,BH_FIND_SURFACE);
    
  create_new_grid(grid_char,bhbh);
  bh1->grid = bhbh->grid;
  bh2->grid = bhbh->grid;
  
  /* add fields */
  physics(bhbh,ADM_ADD_FIELDS);
  physics(bhbh,FREE_DATA_ADD_FIELDS);
  physics(bhbh,STRESS_ENERGY_ADD_FIELDS);
  physics(bhbh,SYS_ADD_FIELDS);
  physics(bhbh,OBSERVE_ADD_FIELDS);
  physics(bhbh,BH_ADD_FIELDS);
  
  /* populate fields */
  physics(bhbh,FREE_DATA_POPULATE);
  physics(bhbh,SYS_INITIALIZE_FIELDS);
  /* beta = B0+B1 */
  physics(bhbh,ADM_UPDATE_B1I);
  initial_B0I(bhbh,".*");
  update_partial_derivatives(bhbh,".*","^dB0_U.+,^ddB0_U.+");
  physics(bhbh,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(bhbh,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  
  /* update AConf^{ij} */
  physics(bhbh,ADM_UPDATE_AConfIJ);

  /* update normal on AH */
  physics(bh1,BH_UPDATE_sConf);
  physics(bh2,BH_UPDATE_sConf);

  /* free */
  free_physics(bh1);
  free_physics(bh2);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return bhbh;
}

/* based on grid character, make a new grid for the system. */
static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const bhbh)
{
  FUNC_TIC
  
  /* a new grid */
  Grid_T *const grid = elliptica_alloc_grid();
  const double bh1_box_len_ratio = 0.2;/* experimentally */
  const double bh2_box_len_ratio = 0.2;/* experimentally */
  
  /* grid for characters */
  grid_char->grid = grid;
  
  if (!Pcmps("grid_kind","SplitCubedSpherical(BH+BH)"))
    Error0(NO_OPTION);
  
  /* set "grid_BH1_central_box_length" and "grid_BH2_central_box_length"
  // automatically, if it is asked. 
  // NOTE: these params are set only for the very first time.
  // one must note that if mass of the object is being iterated, 
  // this auto option is not very ideal, since the final radius is not
  // known yet */
  if (Pcmps("grid_BH1_central_box_length","auto"))
    Psetd("grid_BH1_central_box_length",
          bh1_box_len_ratio*grid_char->params[Ibh1]->r_min);
          
  if (Pcmps("grid_BH2_central_box_length","auto"))
    Psetd("grid_BH2_central_box_length",
          bh2_box_len_ratio*grid_char->params[Ibh2]->r_min);
  
  /* separation */
  grid_char->S               = Pgetd("BHBH_separation");
  /* BH1 */
  grid_char->params[Ibh1]->l = Pgetd("grid_BH1_central_box_length");
  grid_char->params[Ibh1]->w = Pgetd("grid_BH1_central_box_length");
  grid_char->params[Ibh1]->h = Pgetd("grid_BH1_central_box_length");
  /* BH2 */
  grid_char->params[Ibh2]->l = Pgetd("grid_BH2_central_box_length");
  grid_char->params[Ibh2]->w = Pgetd("grid_BH2_central_box_length");
  grid_char->params[Ibh2]->h = Pgetd("grid_BH2_central_box_length");

  /* check central box length */
  if (grid_char->params[Ibh1]->l > grid_char->params[Ibh1]->r_min/2. ||
      grid_char->params[Ibh1]->w > grid_char->params[Ibh1]->r_min/2. ||
      grid_char->params[Ibh1]->h > grid_char->params[Ibh1]->r_min/2.)
    Error0("BH1 central box is too big!");
  
  /* check central box length */
  if (grid_char->params[Ibh2]->l > grid_char->params[Ibh2]->r_min/2. ||
      grid_char->params[Ibh2]->w > grid_char->params[Ibh2]->r_min/2. ||
      grid_char->params[Ibh2]->h > grid_char->params[Ibh2]->r_min/2.)
    Error0("BH2 central box is too big!");
  
  set_params_of_split_cubed_spherical_grid(grid_char);
    
  make_patches(grid);
  realize_interfaces(grid);
  
  bhbh->grid = grid;
  
  FUNC_TOC
}

/* update partial derivatives of the given field name regex match 
// at the given region */
static void update_partial_derivatives(Physics_T *const phys,
                                       const char *const region,
                                       const char *const regex)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  printf(Pretty0"%s\n",regex);
  fflush(stdout);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    partial_derivative_regex(patch,regex);
  }
  
  FUNC_TOC
}

/* initial B0^i in beta = B0+B1 
// also at outermost patches set psi = psi/r */
static void initial_B0I(Physics_T *const phys,
                       const char *const region)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    READ_v(beta_U0);
    READ_v(beta_U1);
    READ_v(beta_U2);
    READ_v(B1_U0);
    READ_v(B1_U1);
    READ_v(B1_U2);
    WRITE_v(psi);
    
    REALLOC_v_WRITE_v(B0_U0);
    REALLOC_v_WRITE_v(B0_U1);
    REALLOC_v_WRITE_v(B0_U2);
    
    FOR_ALL_ijk
    {
      double psim4 = pow(psi[ijk],-4.);
      
      B0_U0[ijk] = psim4*beta_U0[ijk]-B1_U0[ijk];
      B0_U1[ijk] = psim4*beta_U1[ijk]-B1_U1[ijk];
      B0_U2[ijk] = psim4*beta_U2[ijk]-B1_U2[ijk];
    }
    
    /* attenuate */
    if (IsItCovering(patch,"outermost"))
    {
      FOR_ALL_ijk
      {
        DEF_RELATIVE_x
        DEF_RELATIVE_y
        DEF_RELATIVE_z
        DEF_RELATIVE_r
        
        psi[ijk]   /= r;
        B0_U0[ijk] /= r;
        B0_U1[ijk] /= r;
        B0_U2[ijk] /= r;
      }
    }
  }
  
  FUNC_TOC
}

/* loading from checkpoint */
Physics_T *bhbh_read_physics_from_checkpoint(void)
{
  FUNC_TIC
  Physics_T *const bhbh = init_physics(0,BHBH);
  FILE *file = 0;
  
  /* first load grid and parameters */
  file = open_checkpoint_file_then_read_grid_and_params(bhbh);
  
  /* it already hit the stop */
  if (Pgeti(P_"STOP") == 1)
  {
    printf(Pretty0" All iterations have already done.\n");
    free_physics(bhbh);
    Fclose(file);
    
    FUNC_TOC
    return 0;
  }
  
  Physics_T *const bh1 = init_physics(bhbh,BH1);
  Physics_T *const bh2 = init_physics(bhbh,BH2);
  
  /* make the patches */
  make_patches(bhbh->grid);
  
  /* realizing the geometry */
  realize_interfaces(bhbh->grid);
  
  /* set parameters, it's important to add paramters 
  // since these call also reposible to set default functions. */
  physics(bhbh,FREE_DATA_SET_PARAMS);
  physics(bhbh,ADM_SET_PARAMS);
  physics(bhbh,SYS_SET_PARAMS);
  physics(bhbh,STRESS_ENERGY_SET_PARAMS);
  physics(bhbh,OBSERVE_SET_PARAMS);  
  physics(bhbh,BH_SET_PARAMS);
  
  /* now add fields */
  physics(bhbh,ADM_ADD_FIELDS);
  physics(bhbh,FREE_DATA_ADD_FIELDS);
  physics(bhbh,STRESS_ENERGY_ADD_FIELDS);
  physics(bhbh,SYS_ADD_FIELDS);
  physics(bhbh,OBSERVE_ADD_FIELDS);
  physics(bhbh,BH_ADD_FIELDS);
  
  /* populate free data fields */
  physics(bhbh,FREE_DATA_POPULATE);
  
  /* then read saved fields in param "checkpoint_save" */
  read_fields_from_checkpoint_file(bhbh,file);
  Fclose(file);
  
  /* beta = B0+B1 */
  physics(bhbh,ADM_UPDATE_B1I);
  update_partial_derivatives(bhbh,".*","^dB0_U.+,^ddB0_U.+");
  physics(bhbh,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(bhbh,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  
  /* update AConf^{ij} */
  physics(bhbh,ADM_UPDATE_AConfIJ);
  
  /* update normal on AH */
  physics(bh1,BH_UPDATE_sConf);
  physics(bh2,BH_UPDATE_sConf);
  
  free_physics(bh1);
  free_physics(bh2);

  FUNC_TOC
  return bhbh;
}

/* using copy or interpolation from old physics to 
// initialize fields for new physics */
static void initialize_fields_using_previous_solve
            (Physics_T *const new_phys,Physics_T *const old_phys)
{
  FUNC_TIC
  
  Physics_T *const old_bh1 = init_physics(old_phys,BH1);
  Physics_T *const new_bh1 = init_physics(new_phys,BH1);
  Physics_T *const old_bh2 = init_physics(old_phys,BH2);
  Physics_T *const new_bh2 = init_physics(new_phys,BH2);
  
  /* if resolution changed */
  if(Pgeti(P_"did_resolution_change?"))
  {
    interpolate_fields_from_old_grid_to_new_grid
      (old_phys->grid,new_phys->grid,"psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
  }
  else
  {
    const char *region1 = 0;
    const char *region2 = 0;
    if (new_phys->grid->kind == Grid_SplitCubedSpherical_BHBH)
    {
      /* since filling_box,outermost are fixed, only copy */
      region1 = "filling_box,outermost";
      region2 = "filling_box,outermost";
      interpolate_fields_from_old_grid_to_new_grid
        (mygrid(old_phys,region1),mygrid(new_phys,region2),
         "psi,alphaPsi,B0_U0,B0_U1,B0_U2",1);
      
      if (Pgeti("BH1_did_BH_surface_change?"))
      {
        /* if BH is empty, i.e., BH has not already been filled
        // we should fill the BH outerwise
        // "interpolate_fields_from_old_grid_to_new_grid" gives error.
        // note: obviously, the new grid has not been filled so this
        // only regards the old grid and so old_bh. */
        if (Pgeti("BH1_was_BH_filled?") == 0)
        {
          Psets("BH1_filler_fields","alphaPsi,psi,B0_U0,B0_U1,B0_U2");
          physics(old_bh1,BH_FILL);
          Pseti("BH1_was_BH_filled?",1);
        }
        region1 = "BH1,BH1_around";
        region2 = "BH1,BH1_around";
        interpolate_fields_from_old_grid_to_new_grid
          (mygrid(old_bh1,region1),mygrid(new_bh1,region2),
           "psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
      } // if (Pgeti("BH1_did_BH_surface_change?"))
      else
      {
        region1 = "BH1_around";
        region2 = "BH1_around";
        interpolate_fields_from_old_grid_to_new_grid
          (mygrid(old_bh1,region1),mygrid(new_bh1,region2),
           "psi,alphaPsi,B0_U0,B0_U1,B0_U2",1);
      }
      
      if (Pgeti("BH2_did_BH_surface_change?"))
      {
        /* if BH is empty, i.e., BH has not already been filled
        // we should fill the BH outerwise
        // "interpolate_fields_from_old_grid_to_new_grid" gives error.
        // note: obviously, the new grid has not been filled so this
        // only regards the old grid and so old_bh. */
        if (Pgeti("BH2_was_BH_filled?") == 0)
        {
          Psets("BH2_filler_fields","alphaPsi,psi,B0_U0,B0_U1,B0_U2");
          physics(old_bh2,BH_FILL);
          Pseti("BH2_was_BH_filled?",1);
        }
        region1 = "BH2,BH2_around";
        region2 = "BH2,BH2_around";
        interpolate_fields_from_old_grid_to_new_grid
          (mygrid(old_bh2,region1),mygrid(new_bh2,region2),
           "psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
      } // if (Pgeti("BH2_did_BH_surface_change?"))
      else
      {
        region1 = "BH2_around";
        region2 = "BH2_around";
        interpolate_fields_from_old_grid_to_new_grid
          (mygrid(old_bh2,region1),mygrid(new_bh2,region2),
           "psi,alphaPsi,B0_U0,B0_U1,B0_U2",1);
      }
         
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  
  free_physics(old_bh1);
  free_physics(new_bh1);
  free_physics(old_bh2);
  free_physics(new_bh2);
  
  FUNC_TOC
}

/* update some parameters for a new physics */
static void update_params(Physics_T *const phys)
{
  FUNC_TIC  
  UNUSED(phys);
  FUNC_TOC
}

/* move Jacobian if possible */
static void move_jacobian
            (Physics_T *const new_phys,Physics_T *const old_phys)
{
  FUNC_TIC
  
  if(Pgeti(P_"did_resolution_change?") || 
     !new_phys->grid                   || 
     !old_phys->grid)
  {
    FUNC_TOC
    return;
  }
  
  Physics_T *const old_bh1 = init_physics(old_phys,BH1);
  Physics_T *const new_bh1 = init_physics(new_phys,BH1);
  Physics_T *const old_bh2 = init_physics(old_phys,BH2);
  Physics_T *const new_bh2 = init_physics(new_phys,BH2);
  Grid_T *gnew = 0;
  Grid_T *gold = 0;
  const char *name1 = 0;
  const char *name2 = 0;
  
  if(new_phys->grid->kind == Grid_SplitCubedSpherical_BHBH)
  {
    /* move Jacobian of outermost */
    gnew = mygrid(new_phys,"outermost");
    gold = mygrid(old_phys,"outermost");
    
    FOR_ALL_p(gnew->np)
    {
      name1 = strchr(gnew->patch[p]->name,'_');
      for(Uint p2 = 0; p2 < gold->np; ++p2)
      {
        name2 = strchr(gold->patch[p2]->name,'_');
        if(!strcmp(name1,name2))
        {
          break;
        }
      }
    }
    
  }/* if(new_phys->grid->kind == Grid_SplitCubedSpherical_BHBH) */
  else
  {
    Error0(NO_OPTION);
  }
  
  free_physics(old_bh1);
  free_physics(new_bh1);
  free_physics(old_bh2);
  free_physics(new_bh2);

  FUNC_TOC
}

