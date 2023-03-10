/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "obs_header.h"


double obs_ADM_mass_S2(Observe_T *const obs);
double obs_ADM_mass_S2(Observe_T *const obs)
{
  double adm_mass = 0;
  struct items_S **adm = obs->items;
  const Uint N = obs->Nitems;
  Uint p;

  for(p = 0; p < N; ++p)
  {
  Patch_T *patch = adm[p]->patch;


  /* declaring: */
  READ_v_UNUSED(psi)
  READ_v_UNUSED(dpsi_D0)
  READ_v_UNUSED(dpsi_D1)
  READ_v_UNUSED(dpsi_D2)
  READ_v_UNUSED(dgConf_D0D2D1)
  READ_v_UNUSED(dgConf_D0D2D2)
  READ_v_UNUSED(dgConf_D1D1D1)
  READ_v_UNUSED(dgConf_D1D2D1)
  READ_v_UNUSED(dgConf_D0D1D0)
  READ_v_UNUSED(dgConf_D0D1D1)
  READ_v_UNUSED(dgConf_D0D1D2)
  READ_v_UNUSED(dgConf_D0D2D0)
  READ_v_UNUSED(dgConf_D0D0D1)
  READ_v_UNUSED(dgConf_D0D0D0)
  READ_v_UNUSED(dgConf_D1D2D0)
  READ_v_UNUSED(dgConf_D0D0D2)
  READ_v_UNUSED(dgConf_D1D2D2)
  READ_v_UNUSED(dgConf_D1D1D0)
  READ_v_UNUSED(dgConf_D1D1D2)
  READ_v_UNUSED(dgConf_D2D2D1)
  READ_v_UNUSED(dgConf_D2D2D0)
  READ_v_UNUSED(dgConf_D2D2D2)
  READ_v_UNUSED(igConf_U2U2)
  READ_v_UNUSED(igConf_U1U2)
  READ_v_UNUSED(igConf_U1U1)
  READ_v_UNUSED(igConf_U0U2)
  READ_v_UNUSED(igConf_U0U0)
  READ_v_UNUSED(igConf_U0U1)


    Uint nn = patch->nn;
    Uint ijk;

    if (adm[p]->surface_integration_flg)
    {
      ADD_FIELD(ADM_mass_integrand_S)
      {
      const double *n_U0 = adm[p]->n_U0;
      const double *n_U1 = adm[p]->n_U1;
      const double *n_U2 = adm[p]->n_U2;
      REALLOC_v_WRITE_v(ADM_mass_integrand_S)
      for (ijk = 0; ijk < nn; ++ijk)
      {
      double dlnpsi_U2 = 
dpsi_D2[ijk]/psi[ijk];

      double dlnpsi_U1 = 
dpsi_D1[ijk]/psi[ijk];

      double dlnpsi_U0 = 
dpsi_D0[ijk]/psi[ijk];

      double M_s = 
n_U0[ijk]*(-8*dlnpsi_U0 + igConf_U0U1[ijk]*(dgConf_D0D0D1[ijk] -
dgConf_D0D1D0[ijk]) + igConf_U0U2[ijk]*(dgConf_D0D0D2[ijk] -
dgConf_D0D2D0[ijk]) + igConf_U1U1[ijk]*(dgConf_D0D1D1[ijk] -
dgConf_D1D1D0[ijk]) + igConf_U1U2[ijk]*(dgConf_D0D1D2[ijk] -
dgConf_D1D2D0[ijk]) + igConf_U1U2[ijk]*(dgConf_D0D2D1[ijk] -
dgConf_D1D2D0[ijk]) + igConf_U2U2[ijk]*(dgConf_D0D2D2[ijk] -
dgConf_D2D2D0[ijk])) - n_U1[ijk]*(8*dlnpsi_U1 + igConf_U0U0[ijk]*
(dgConf_D0D0D1[ijk] - dgConf_D0D1D0[ijk]) + igConf_U0U1[ijk]*
(dgConf_D0D1D1[ijk] - dgConf_D1D1D0[ijk]) - igConf_U0U2[ijk]*
(dgConf_D0D1D2[ijk] - dgConf_D0D2D1[ijk]) + igConf_U0U2[ijk]*
(dgConf_D0D2D1[ijk] - dgConf_D1D2D0[ijk]) - igConf_U1U2[ijk]*
(dgConf_D1D1D2[ijk] - dgConf_D1D2D1[ijk]) - igConf_U2U2[ijk]*
(dgConf_D1D2D2[ijk] - dgConf_D2D2D1[ijk])) - n_U2[ijk]*(8*dlnpsi_U2 +
igConf_U0U0[ijk]*(dgConf_D0D0D2[ijk] - dgConf_D0D2D0[ijk]) +
igConf_U0U1[ijk]*(dgConf_D0D1D2[ijk] - dgConf_D0D2D1[ijk]) +
igConf_U0U1[ijk]*(dgConf_D0D1D2[ijk] - dgConf_D1D2D0[ijk]) +
igConf_U0U2[ijk]*(dgConf_D0D2D2[ijk] - dgConf_D2D2D0[ijk]) +
igConf_U1U1[ijk]*(dgConf_D1D1D2[ijk] - dgConf_D1D2D1[ijk]) +
igConf_U1U2[ijk]*(dgConf_D1D2D2[ijk] - dgConf_D2D2D1[ijk]));

      ADM_mass_integrand_S[ijk] = M_s;
      }
      }
      DECLARE_FIELD(ADM_mass_integrand_S)
      Integration_T *I = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      I->Spectral->f = ADM_mass_integrand_S;
      I->g00 = adm[p]->g00;
      I->g01 = adm[p]->g01;
      I->g02 = adm[p]->g02;
      I->g11 = adm[p]->g11;
      I->g12 = adm[p]->g12;
      I->g22 = adm[p]->g22;
      I->Spectral->X_surface = adm[p]->X_surface;
      I->Spectral->Y_surface = adm[p]->Y_surface;
      I->Spectral->Z_surface = adm[p]->Z_surface;
      I->Spectral->I         = adm[p]->I;
      I->Spectral->J         = adm[p]->J;
      I->Spectral->K         = adm[p]->K;
      plan_integration(I);
      adm_mass += execute_integration(I);
      free_integration(I);
      REMOVE_FIELD(ADM_mass_integrand_S)
    }
    else
    {
      Error0("Wrong flag!\n");
    }

  }
  adm_mass /= (16*M_PI);
  return adm_mass;
}
