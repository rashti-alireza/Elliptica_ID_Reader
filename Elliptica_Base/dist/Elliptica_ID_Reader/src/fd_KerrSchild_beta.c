/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2021 Alireza Rashti.
*/


#include "fd_header.h"
#include "fd_KerrSchild_header.h"


#undef x
#undef y
#undef z

#define KS_func_pass_args_sed KS_func_pass_args_macro

#define KS_set_args \
  struct Analytic_Func_Arg_S farg[1];\
  farg->x = x;\
  farg->y = y;\
  farg->z = z;\
  farg->X = fd_ks_X KS_func_pass_args_macro;\
  farg->Y = fd_ks_Y KS_func_pass_args_macro;\
  farg->Z = fd_ks_Z KS_func_pass_args_macro;\
  farg->R = fd_ks_R KS_func_pass_args_macro;\
  farg->dX_D0 = fd_ks_dX_D0 KS_func_pass_args_sed;\
  farg->dX_D1 = fd_ks_dX_D1 KS_func_pass_args_sed;\
  farg->dX_D2 = fd_ks_dX_D2 KS_func_pass_args_sed;\
  farg->dY_D0 = fd_ks_dY_D0 KS_func_pass_args_sed;\
  farg->dY_D1 = fd_ks_dY_D1 KS_func_pass_args_sed;\
  farg->dY_D2 = fd_ks_dY_D2 KS_func_pass_args_sed;\
  farg->dZ_D0 = fd_ks_dZ_D0 KS_func_pass_args_sed;\
  farg->dZ_D1 = fd_ks_dZ_D1 KS_func_pass_args_sed;\
  farg->dZ_D2 = fd_ks_dZ_D2 KS_func_pass_args_sed;


void fd_beta_KerrSchild_patch(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Beta);


void fd_beta_KerrSchild_patch(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Beta)
{

  /* declaring: */
  READ_v_STEM(ig_U1U1,ig)
  READ_v_STEM(ig_U1U2,ig)
  READ_v_STEM(ig_U0U0,ig)
  READ_v_STEM(ig_U0U1,ig)
  READ_v_STEM(ig_U0U2,ig)
  READ_v_STEM(ig_U2U2,ig)
  REALLOC_v_WRITE_v_STEM(beta_U1,Beta)
  REALLOC_v_WRITE_v_STEM(beta_U0,Beta)
  REALLOC_v_WRITE_v_STEM(beta_U2,Beta)


FOR_ALL_ijk
{
  double x,y,z;
  x = patch->node[ijk]->x[0]-BH_center_x;
  y = patch->node[ijk]->x[1]-BH_center_y;
  z = patch->node[ijk]->x[2]-BH_center_z;
  KS_set_args

  double betaD_D0 = fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_k0(x, y, z);
  double betaD_D1 = fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_k1(x, y, z);
  double betaD_D2 = fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_k2(x, y, z);
  double betaU_U1 = 
betaD_D0*ig_U0U1[ijk] + betaD_D1*ig_U1U1[ijk] + betaD_D2*
ig_U1U2[ijk];

  double betaU_U0 = 
betaD_D0*ig_U0U0[ijk] + betaD_D1*ig_U0U1[ijk] + betaD_D2*
ig_U0U2[ijk];

  double betaU_U2 = 
betaD_D0*ig_U0U2[ijk] + betaD_D1*ig_U1U2[ijk] + betaD_D2*
ig_U2U2[ijk];


  /* populating: */
  beta_U1[ijk] = betaU_U1;
  beta_U0[ijk] = betaU_U0;
  beta_U2[ijk] = betaU_U2;
}
}