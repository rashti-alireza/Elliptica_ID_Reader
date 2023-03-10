/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "dg_header.h"

#define STR_LEN (99)

#define read_ig(x) indx = strrchr(#x,'_'); sprintf(fname,"%s%s",ig,indx); \
    const double *const x = patch->fields[Ind(fname)]->v;

#define read_dg(x) indx = strrchr(#x,'_'); sprintf(fname,"%s%s",dg,indx); \
    const double *const x = patch->fields[Ind(fname)]->v;

#define write_Chris(x) indx = strrchr(#x,'_'); sprintf(fname,"%s%s",Chris,indx); \
    free_coeffs(patch->fields[Ind(fname)]); \
    double *const x = patch->fields[Ind(fname)]->v; assert(x);


void Christoffel_symbol_3d(Patch_T *const patch,const char *const ig,const char *const dg,const char *const Chris);
void Christoffel_symbol_3d(Patch_T *const patch,const char *const ig,const char *const dg,const char *const Chris)
{
char fname[STR_LEN] = {'\0'};
const char *indx = 0;

  /* declaring: */
  read_ig(ig_U1U1)
  read_ig(ig_U1U2)
  read_ig(ig_U0U0)
  read_ig(ig_U0U1)
  read_ig(ig_U0U2)
  read_ig(ig_U2U2)
  read_dg(dg_D1D2D0)
  read_dg(dg_D0D0D2)
  read_dg(dg_D0D0D1)
  read_dg(dg_D0D0D0)
  read_dg(dg_D1D2D1)
  read_dg(dg_D1D1D2)
  read_dg(dg_D0D2D1)
  read_dg(dg_D0D2D0)
  read_dg(dg_D0D2D2)
  read_dg(dg_D1D1D1)
  read_dg(dg_D0D1D2)
  read_dg(dg_D1D1D0)
  read_dg(dg_D0D1D0)
  read_dg(dg_D0D1D1)
  read_dg(dg_D2D2D2)
  read_dg(dg_D1D2D2)
  read_dg(dg_D2D2D1)
  read_dg(dg_D2D2D0)
  write_Chris(Chris_U2D2D2)
  write_Chris(Chris_U0D1D1)
  write_Chris(Chris_U2D1D1)
  write_Chris(Chris_U0D1D2)
  write_Chris(Chris_U0D0D2)
  write_Chris(Chris_U0D0D1)
  write_Chris(Chris_U0D0D0)
  write_Chris(Chris_U2D1D2)
  write_Chris(Chris_U2D0D1)
  write_Chris(Chris_U1D2D2)
  write_Chris(Chris_U2D0D0)
  write_Chris(Chris_U2D0D2)
  write_Chris(Chris_U1D1D2)
  write_Chris(Chris_U1D0D2)
  write_Chris(Chris_U1D0D1)
  write_Chris(Chris_U1D0D0)
  write_Chris(Chris_U0D2D2)
  write_Chris(Chris_U1D1D1)


FOR_ALL_ijk
{
  double GAMMA_U1D0D0 = 
0.5*dg_D0D0D0[ijk]*ig_U0U1[ijk] - 0.5*ig_U1U1[ijk]*(dg_D0D0D1[ijk] - 2*
dg_D0D1D0[ijk]) - 0.5*ig_U1U2[ijk]*(dg_D0D0D2[ijk] - 2*
dg_D0D2D0[ijk]);

  double GAMMA_U0D2D2 = 
0.5*dg_D2D2D2[ijk]*ig_U0U2[ijk] + 0.5*ig_U0U0[ijk]*(2*dg_D0D2D2[ijk] - 
dg_D2D2D0[ijk]) + 0.5*ig_U0U1[ijk]*(2*dg_D1D2D2[ijk] - 
dg_D2D2D1[ijk]);

  double GAMMA_U1D0D2 = 
0.5*dg_D0D0D2[ijk]*ig_U0U1[ijk] + 0.5*dg_D2D2D0[ijk]*ig_U1U2[ijk] + 0.5*
ig_U1U1[ijk]*(dg_D0D1D2[ijk] - dg_D0D2D1[ijk] + dg_D1D2D0[ijk]);

  double GAMMA_U2D2D2 = 
0.5*dg_D2D2D2[ijk]*ig_U2U2[ijk] + 0.5*ig_U0U2[ijk]*(2*dg_D0D2D2[ijk] - 
dg_D2D2D0[ijk]) + 0.5*ig_U1U2[ijk]*(2*dg_D1D2D2[ijk] - 
dg_D2D2D1[ijk]);

  double GAMMA_U2D1D1 = 
0.5*dg_D1D1D1[ijk]*ig_U1U2[ijk] + 0.5*ig_U0U2[ijk]*(2*dg_D0D1D1[ijk] - 
dg_D1D1D0[ijk]) - 0.5*ig_U2U2[ijk]*(dg_D1D1D2[ijk] - 2*
dg_D1D2D1[ijk]);

  double GAMMA_U1D0D1 = 
0.5*dg_D0D0D1[ijk]*ig_U0U1[ijk] + 0.5*dg_D1D1D0[ijk]*ig_U1U1[ijk] + 0.5*
ig_U1U2[ijk]*(-dg_D0D1D2[ijk] + dg_D0D2D1[ijk] + dg_D1D2D0[ijk]);

  double GAMMA_U1D2D2 = 
0.5*dg_D2D2D2[ijk]*ig_U1U2[ijk] + 0.5*ig_U0U1[ijk]*(2*dg_D0D2D2[ijk] - 
dg_D2D2D0[ijk]) + 0.5*ig_U1U1[ijk]*(2*dg_D1D2D2[ijk] - 
dg_D2D2D1[ijk]);

  double GAMMA_U2D1D2 = 
0.5*dg_D1D1D2[ijk]*ig_U1U2[ijk] + 0.5*dg_D2D2D1[ijk]*ig_U2U2[ijk] + 0.5*
ig_U0U2[ijk]*(dg_D0D1D2[ijk] + dg_D0D2D1[ijk] - dg_D1D2D0[ijk]);

  double GAMMA_U2D0D2 = 
0.5*dg_D0D0D2[ijk]*ig_U0U2[ijk] + 0.5*dg_D2D2D0[ijk]*ig_U2U2[ijk] + 0.5*
ig_U1U2[ijk]*(dg_D0D1D2[ijk] - dg_D0D2D1[ijk] + dg_D1D2D0[ijk]);

  double GAMMA_U2D0D1 = 
0.5*dg_D0D0D1[ijk]*ig_U0U2[ijk] + 0.5*dg_D1D1D0[ijk]*ig_U1U2[ijk] + 0.5*
ig_U2U2[ijk]*(-dg_D0D1D2[ijk] + dg_D0D2D1[ijk] + dg_D1D2D0[ijk]);

  double GAMMA_U2D0D0 = 
0.5*dg_D0D0D0[ijk]*ig_U0U2[ijk] - 0.5*ig_U1U2[ijk]*(dg_D0D0D1[ijk] - 2*
dg_D0D1D0[ijk]) - 0.5*ig_U2U2[ijk]*(dg_D0D0D2[ijk] - 2*
dg_D0D2D0[ijk]);

  double GAMMA_U0D0D2 = 
0.5*dg_D0D0D2[ijk]*ig_U0U0[ijk] + 0.5*dg_D2D2D0[ijk]*ig_U0U2[ijk] + 0.5*
ig_U0U1[ijk]*(dg_D0D1D2[ijk] - dg_D0D2D1[ijk] + dg_D1D2D0[ijk]);

  double GAMMA_U1D1D1 = 
0.5*dg_D1D1D1[ijk]*ig_U1U1[ijk] + 0.5*ig_U0U1[ijk]*(2*dg_D0D1D1[ijk] - 
dg_D1D1D0[ijk]) - 0.5*ig_U1U2[ijk]*(dg_D1D1D2[ijk] - 2*
dg_D1D2D1[ijk]);

  double GAMMA_U0D1D1 = 
0.5*dg_D1D1D1[ijk]*ig_U0U1[ijk] + 0.5*ig_U0U0[ijk]*(2*dg_D0D1D1[ijk] - 
dg_D1D1D0[ijk]) - 0.5*ig_U0U2[ijk]*(dg_D1D1D2[ijk] - 2*
dg_D1D2D1[ijk]);

  double GAMMA_U0D0D0 = 
0.5*dg_D0D0D0[ijk]*ig_U0U0[ijk] - 0.5*ig_U0U1[ijk]*(dg_D0D0D1[ijk] - 2*
dg_D0D1D0[ijk]) - 0.5*ig_U0U2[ijk]*(dg_D0D0D2[ijk] - 2*
dg_D0D2D0[ijk]);

  double GAMMA_U1D1D2 = 
0.5*dg_D1D1D2[ijk]*ig_U1U1[ijk] + 0.5*dg_D2D2D1[ijk]*ig_U1U2[ijk] + 0.5*
ig_U0U1[ijk]*(dg_D0D1D2[ijk] + dg_D0D2D1[ijk] - dg_D1D2D0[ijk]);

  double GAMMA_U0D0D1 = 
0.5*dg_D0D0D1[ijk]*ig_U0U0[ijk] + 0.5*dg_D1D1D0[ijk]*ig_U0U1[ijk] + 0.5*
ig_U0U2[ijk]*(-dg_D0D1D2[ijk] + dg_D0D2D1[ijk] + dg_D1D2D0[ijk]);

  double GAMMA_U0D1D2 = 
0.5*dg_D1D1D2[ijk]*ig_U0U1[ijk] + 0.5*dg_D2D2D1[ijk]*ig_U0U2[ijk] + 0.5*
ig_U0U0[ijk]*(dg_D0D1D2[ijk] + dg_D0D2D1[ijk] - dg_D1D2D0[ijk]);


  /* populating: */
  Chris_U2D2D2[ijk] = GAMMA_U2D2D2;
  Chris_U0D1D1[ijk] = GAMMA_U0D1D1;
  Chris_U2D1D1[ijk] = GAMMA_U2D1D1;
  Chris_U0D1D2[ijk] = GAMMA_U0D1D2;
  Chris_U0D0D2[ijk] = GAMMA_U0D0D2;
  Chris_U0D0D1[ijk] = GAMMA_U0D0D1;
  Chris_U0D0D0[ijk] = GAMMA_U0D0D0;
  Chris_U2D1D2[ijk] = GAMMA_U2D1D2;
  Chris_U2D0D1[ijk] = GAMMA_U2D0D1;
  Chris_U1D2D2[ijk] = GAMMA_U1D2D2;
  Chris_U2D0D0[ijk] = GAMMA_U2D0D0;
  Chris_U2D0D2[ijk] = GAMMA_U2D0D2;
  Chris_U1D1D2[ijk] = GAMMA_U1D1D2;
  Chris_U1D0D2[ijk] = GAMMA_U1D0D2;
  Chris_U1D0D1[ijk] = GAMMA_U1D0D1;
  Chris_U1D0D0[ijk] = GAMMA_U1D0D0;
  Chris_U0D2D2[ijk] = GAMMA_U0D2D2;
  Chris_U1D1D1[ijk] = GAMMA_U1D1D1;
}
}
