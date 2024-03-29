/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2023 Alireza Rashti.
*/


#include "obs_header.h"

void obs_BH_irreducible_mass_CS(Observe_T *const obs);
void obs_BH_irreducible_mass_CS(Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  Grid_T *const grid = mygrid(phys,Ftype("BH_around_IB"));
  double A_AH = 0;
  Uint p;

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Uint nn    = patch->nn;
    Uint ijk;

    ADD_FIELD(A_AH_integrand)

    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    READ_v(gConf_D0D0)
    READ_v(gConf_D0D1)
    READ_v(gConf_D0D2)
    READ_v(gConf_D1D1)
    READ_v(gConf_D1D2)
    READ_v(gConf_D2D2)
    READ_v(psi)
{
    REALLOC_v_WRITE_v(A_AH_integrand)
    for(ijk = 0; ijk < nn; ++ijk)
    {
    double psi4 =
pow(psi[ijk], 4);

      A_AH_integrand[ijk] = 1;
      /* metric */
      g00[ijk] = psi4*gConf_D0D0[ijk];
      g01[ijk] = psi4*gConf_D0D1[ijk];
      g02[ijk] = psi4*gConf_D0D2[ijk];
      g11[ijk] = psi4*gConf_D1D1[ijk];
      g12[ijk] = psi4*gConf_D1D2[ijk];
      g22[ijk] = psi4*gConf_D2D2[ijk];
    }
}

  DECLARE_FIELD(A_AH_integrand)
  Integration_T *I = init_integration();
  I->type = "Integral{f(x)dS},Spectral";
  I->Spectral->f = A_AH_integrand;
  I->g00 = g00;
  I->g01 = g01;
  I->g02 = g02;
  I->g11 = g11;
  I->g12 = g12;
  I->g22 = g22;
  I->Spectral->Z_surface = 1;
  I->Spectral->K         = 0;
  plan_integration(I);
  A_AH += execute_integration(I);

  free_integration(I);
  REMOVE_FIELD(A_AH_integrand)
  free(g00);
  free(g01);
  free(g02);
  free(g11);
  free(g12);
  free(g22);
  }

  obs->ret[0] = sqrt(A_AH/(16*M_PI));
  obs->ret[1] = A_AH;
}
