/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "fd_header.h"
#include "fd_KerrSchild_header.h"


#define add_alloc_get(name) ADD_AND_ALLOC_FIELD(name);WRITE_v(name);


#define dfield_and_get_v(x) dField_di(x); READ_v(x);

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


void fd_Kij_trK_KerrSchild(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Chris,const char *const Kij,
 const char *const trK);

void fd_Kij_trK_KerrSchild(Patch_T *const patch,
 const double BH_center_x,const double BH_center_y,
 const double BH_center_z,const char *const ig,
 const char *const Chris,const char *const Kij,
 const char *const trK)

{

  /* declaring: */
  READ_v_STEM(ig_U1U1,ig)
  READ_v_STEM(ig_U1U2,ig)
  READ_v_STEM(ig_U0U0,ig)
  READ_v_STEM(ig_U0U1,ig)
  READ_v_STEM(ig_U0U2,ig)
  READ_v_STEM(ig_U2U2,ig)
  READ_v_STEM(Chris_U2D2D2,Chris)
  READ_v_STEM(Chris_U0D1D1,Chris)
  READ_v_STEM(Chris_U2D1D1,Chris)
  READ_v_STEM(Chris_U0D1D2,Chris)
  READ_v_STEM(Chris_U0D0D2,Chris)
  READ_v_STEM(Chris_U0D0D1,Chris)
  READ_v_STEM(Chris_U0D0D0,Chris)
  READ_v_STEM(Chris_U2D1D2,Chris)
  READ_v_STEM(Chris_U2D0D1,Chris)
  READ_v_STEM(Chris_U1D2D2,Chris)
  READ_v_STEM(Chris_U2D0D0,Chris)
  READ_v_STEM(Chris_U2D0D2,Chris)
  READ_v_STEM(Chris_U1D1D2,Chris)
  READ_v_STEM(Chris_U1D0D2,Chris)
  READ_v_STEM(Chris_U1D0D1,Chris)
  READ_v_STEM(Chris_U1D0D0,Chris)
  READ_v_STEM(Chris_U0D2D2,Chris)
  READ_v_STEM(Chris_U1D1D1,Chris)
  REALLOC_v_WRITE_v_STEM(Kij_D2D2,Kij)
  REALLOC_v_WRITE_v_STEM(Kij_D0D1,Kij)
  REALLOC_v_WRITE_v_STEM(Kij_D0D0,Kij)
  REALLOC_v_WRITE_v_STEM(Kij_D0D2,Kij)
  REALLOC_v_WRITE_v_STEM(Kij_D1D1,Kij)
  REALLOC_v_WRITE_v_STEM(Kij_D1D2,Kij)
  REALLOC_v_WRITE_v_STEM(trKij,trK)
  add_alloc_get(KS__beta_D2)
  add_alloc_get(KS__beta_D0)
  add_alloc_get(KS__beta_D1)
  ADD_FIELD(dKS__beta_D0D2)
  ADD_FIELD(dKS__beta_D0D1)
  ADD_FIELD(dKS__beta_D0D0)
  ADD_FIELD(dKS__beta_D1D2)
  ADD_FIELD(dKS__beta_D1D0)
  ADD_FIELD(dKS__beta_D1D1)
  ADD_FIELD(dKS__beta_D2D1)
  ADD_FIELD(dKS__beta_D2D0)
  ADD_FIELD(dKS__beta_D2D2)
  add_alloc_get(KS__alpha)


FOR_ALL_ijk
{
  double x,y,z;
  x = patch->node[ijk]->x[0]-BH_center_x;
  y = patch->node[ijk]->x[1]-BH_center_y;
  z = patch->node[ijk]->x[2]-BH_center_z;
  KS_set_args

  KS__beta_D0[ijk] = fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_k0(x, y, z);
  KS__beta_D1[ijk] = fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_k1(x, y, z);
  KS__beta_D2[ijk] = fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_k2(x, y, z);
  KS__alpha[ijk]   = 1./sqrt(1+fd_ks_c(x,y,z)*fd_ks_kt(x, y, z)*fd_ks_kt(x, y, z));
}
dfield_and_get_v(dKS__beta_D0D0);
dfield_and_get_v(dKS__beta_D0D1);
dfield_and_get_v(dKS__beta_D0D2);

dfield_and_get_v(dKS__beta_D1D0);
dfield_and_get_v(dKS__beta_D1D1);
dfield_and_get_v(dKS__beta_D1D2);

dfield_and_get_v(dKS__beta_D2D0);
dfield_and_get_v(dKS__beta_D2D1);
dfield_and_get_v(dKS__beta_D2D2);


FOR_ALL_ijk
{
  double DB_D0D2 = 
-Chris_U0D0D2[ijk]*KS__beta_D0[ijk] - Chris_U1D0D2[ijk]*
KS__beta_D1[ijk] - Chris_U2D0D2[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D0D2[ijk];

  double DB_D0D0 = 
-Chris_U0D0D0[ijk]*KS__beta_D0[ijk] - Chris_U1D0D0[ijk]*
KS__beta_D1[ijk] - Chris_U2D0D0[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D0D0[ijk];

  double DB_D0D1 = 
-Chris_U0D0D1[ijk]*KS__beta_D0[ijk] - Chris_U1D0D1[ijk]*
KS__beta_D1[ijk] - Chris_U2D0D1[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D0D1[ijk];

  double DB_D2D0 = 
-Chris_U0D0D2[ijk]*KS__beta_D0[ijk] - Chris_U1D0D2[ijk]*
KS__beta_D1[ijk] - Chris_U2D0D2[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D2D0[ijk];

  double DB_D2D1 = 
-Chris_U0D1D2[ijk]*KS__beta_D0[ijk] - Chris_U1D1D2[ijk]*
KS__beta_D1[ijk] - Chris_U2D1D2[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D2D1[ijk];

  double DB_D2D2 = 
-Chris_U0D2D2[ijk]*KS__beta_D0[ijk] - Chris_U1D2D2[ijk]*
KS__beta_D1[ijk] - Chris_U2D2D2[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D2D2[ijk];

  double DB_D1D2 = 
-Chris_U0D1D2[ijk]*KS__beta_D0[ijk] - Chris_U1D1D2[ijk]*
KS__beta_D1[ijk] - Chris_U2D1D2[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D1D2[ijk];

  double DB_D1D1 = 
-Chris_U0D1D1[ijk]*KS__beta_D0[ijk] - Chris_U1D1D1[ijk]*
KS__beta_D1[ijk] - Chris_U2D1D1[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D1D1[ijk];

  double DB_D1D0 = 
-Chris_U0D0D1[ijk]*KS__beta_D0[ijk] - Chris_U1D0D1[ijk]*
KS__beta_D1[ijk] - Chris_U2D0D1[ijk]*KS__beta_D2[ijk] + 
dKS__beta_D1D0[ijk];

  double ksKij_D2D2 = 
DB_D2D2/KS__alpha[ijk];

  double ksKij_D0D1 = 
(1.0/2.0)*(DB_D0D1 + DB_D1D0)/KS__alpha[ijk];

  double ksKij_D0D0 = 
DB_D0D0/KS__alpha[ijk];

  double ksKij_D1D2 = 
(1.0/2.0)*(DB_D1D2 + DB_D2D1)/KS__alpha[ijk];

  double ksKij_D0D2 = 
(1.0/2.0)*(DB_D0D2 + DB_D2D0)/KS__alpha[ijk];

  double ksKij_D1D1 = 
DB_D1D1/KS__alpha[ijk];

  double trk = 
ig_U0U0[ijk]*ksKij_D0D0 + 2.0*ig_U0U1[ijk]*ksKij_D0D1 + 2.0*
ig_U0U2[ijk]*ksKij_D0D2 + ig_U1U1[ijk]*ksKij_D1D1 + 2.0*ig_U1U2[ijk]*
ksKij_D1D2 + ig_U2U2[ijk]*ksKij_D2D2;


  /* populating: */
  Kij_D2D2[ijk] = ksKij_D2D2;
  Kij_D0D1[ijk] = ksKij_D0D1;
  Kij_D0D0[ijk] = ksKij_D0D0;
  Kij_D0D2[ijk] = ksKij_D0D2;
  Kij_D1D1[ijk] = ksKij_D1D1;
  Kij_D1D2[ijk] = ksKij_D1D2;
trKij[ijk] = trk;
}
remove_field_regex(patch,"^KS__");
remove_field_regex(patch,"^dKS__");
}
