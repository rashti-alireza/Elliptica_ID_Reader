/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "obs_header.h"


#define add_and_get_field(name) \
  if (_Ind(#name) < 0) \
  {ADD_AND_ALLOC_FIELD(name);} \
  WRITE_v(name);


void obs_ADM_P_Stokes_SV_Ossokine(Observe_T *const obs);
void obs_ADM_P_Stokes_SV_Ossokine(Observe_T *const obs)
{
  struct items_S **adm = obs->items;
  const Uint N = obs->Nitems;
  const double CutOff = 1E8;
  Uint p;

  for(p = 0; p < N; ++p)
  {
  Patch_T *patch = adm[p]->patch;
  Uint nn = patch->nn;
  Uint ijk;

  /* declaring: */
  READ_v(AConfIJ_U0U1)
  READ_v(AConfIJ_U0U0)
  READ_v(AConfIJ_U0U2)
  READ_v(AConfIJ_U2U2)
  READ_v(AConfIJ_U1U1)
  READ_v(AConfIJ_U1U2)
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
  READ_v(gConf_D2D2)
  READ_v(igConf_U2U2)
  READ_v(igConf_U1U2)
  READ_v(igConf_U1U1)
  READ_v(igConf_U0U2)
  READ_v(igConf_U0U0)
  READ_v(igConf_U0U1)
  READ_v(psi)
  READ_v(dpsi_D0)
  READ_v(dpsi_D1)
  READ_v(dpsi_D2)
  READ_v(ChrisConf_U1D0D0)
  READ_v(ChrisConf_U0D2D2)
  READ_v(ChrisConf_U1D0D2)
  READ_v(ChrisConf_U2D2D2)
  READ_v(ChrisConf_U2D1D1)
  READ_v(ChrisConf_U2D0D0)
  READ_v(ChrisConf_U1D1D1)
  READ_v(ChrisConf_U0D1D1)
  READ_v(ChrisConf_U0D1D2)
  READ_v(ChrisConf_U1D1D2)
  READ_v(ChrisConf_U0D0D1)
  READ_v(ChrisConf_U0D0D0)
  READ_v(ChrisConf_U2D0D1)
  READ_v(ChrisConf_U0D0D2)
  READ_v(ChrisConf_U2D1D2)
  READ_v(ChrisConf_U1D2D2)
  READ_v(ChrisConf_U2D0D2)
  READ_v(ChrisConf_U1D0D1)
  READ_v(trK)
  add_and_get_field(obs__P_U1)
  add_and_get_field(obs__P_U0)
  add_and_get_field(obs__P_U2)
  add_and_get_field(obs__G_U0)
  add_and_get_field(obs__G_U1)
  add_and_get_field(obs__G_U2)



   if (adm[p]->surface_integration_flg)
   {
      const double *n_U0 = adm[p]->n_U0;
      const double *n_U1 = adm[p]->n_U1;
      const double *n_U2 = adm[p]->n_U2;
      for (ijk = 0; ijk < nn; ++ijk)
      {
      double psi4 = 
pow(psi[ijk], 4);

      double psi6 = 
pow(psi[ijk], 6);

      double P_U0U1 = 
AConfIJ_U0U1[ijk] - 0.66666666666666663*igConf_U0U1[ijk]*psi6*
trK[ijk];

      double P_U0U0 = 
AConfIJ_U0U0[ijk] - 0.66666666666666663*igConf_U0U0[ijk]*psi6*
trK[ijk];

      double P_U0U2 = 
AConfIJ_U0U2[ijk] - 0.66666666666666663*igConf_U0U2[ijk]*psi6*
trK[ijk];

      double P_U2U2 = 
AConfIJ_U2U2[ijk] - 0.66666666666666663*igConf_U2U2[ijk]*psi6*
trK[ijk];

      double P_U1U1 = 
AConfIJ_U1U1[ijk] - 0.66666666666666663*igConf_U1U1[ijk]*psi6*
trK[ijk];

      double P_U1U2 = 
AConfIJ_U1U2[ijk] - 0.66666666666666663*igConf_U1U2[ijk]*psi6*
trK[ijk];

      double Pn_U2 = 
psi4*(P_U0U2*gConf_D0D0[ijk]*n_U0[ijk] + P_U0U2*gConf_D0D1[ijk]*
n_U1[ijk] + P_U0U2*gConf_D0D2[ijk]*n_U2[ijk] + P_U1U2*gConf_D0D1[ijk]*
n_U0[ijk] + P_U1U2*gConf_D1D1[ijk]*n_U1[ijk] + P_U1U2*gConf_D1D2[ijk]*
n_U2[ijk] + P_U2U2*gConf_D0D2[ijk]*n_U0[ijk] + P_U2U2*gConf_D1D2[ijk]*
n_U1[ijk] + P_U2U2*gConf_D2D2[ijk]*n_U2[ijk]);

      double Pn_U1 = 
psi4*(P_U0U1*gConf_D0D0[ijk]*n_U0[ijk] + P_U0U1*gConf_D0D1[ijk]*
n_U1[ijk] + P_U0U1*gConf_D0D2[ijk]*n_U2[ijk] + P_U1U1*gConf_D0D1[ijk]*
n_U0[ijk] + P_U1U1*gConf_D1D1[ijk]*n_U1[ijk] + P_U1U1*gConf_D1D2[ijk]*
n_U2[ijk] + P_U1U2*gConf_D0D2[ijk]*n_U0[ijk] + P_U1U2*gConf_D1D2[ijk]*
n_U1[ijk] + P_U1U2*gConf_D2D2[ijk]*n_U2[ijk]);

      double Pn_U0 = 
psi4*(P_U0U0*gConf_D0D0[ijk]*n_U0[ijk] + P_U0U0*gConf_D0D1[ijk]*
n_U1[ijk] + P_U0U0*gConf_D0D2[ijk]*n_U2[ijk] + P_U0U1*gConf_D0D1[ijk]*
n_U0[ijk] + P_U0U1*gConf_D1D1[ijk]*n_U1[ijk] + P_U0U1*gConf_D1D2[ijk]*
n_U2[ijk] + P_U0U2*gConf_D0D2[ijk]*n_U0[ijk] + P_U0U2*gConf_D1D2[ijk]*
n_U1[ijk] + P_U0U2*gConf_D2D2[ijk]*n_U2[ijk]);


      /* populating: */
      obs__P_U1[ijk] = Pn_U1;
      obs__P_U0[ijk] = Pn_U0;
      obs__P_U2[ijk] = Pn_U2;
      }
   }
    else
    {
      for (ijk = 0; ijk < nn; ++ijk)
      {
      DEF_RELATIVE_x
      DEF_RELATIVE_y
      DEF_RELATIVE_z
      DEF_RELATIVE_r
      double att = r > CutOff ? 0:1;
      double psi6_ = 
pow(psi[ijk], 6);

      double P__U1U1 = 
AConfIJ_U1U1[ijk] - 0.66666666666666663*igConf_U1U1[ijk]*psi6_*
trK[ijk];

      double P__U1U2 = 
AConfIJ_U1U2[ijk] - 0.66666666666666663*igConf_U1U2[ijk]*psi6_*
trK[ijk];

      double P__U0U1 = 
AConfIJ_U0U1[ijk] - 0.66666666666666663*igConf_U0U1[ijk]*psi6_*
trK[ijk];

      double P__U0U0 = 
AConfIJ_U0U0[ijk] - 0.66666666666666663*igConf_U0U0[ijk]*psi6_*
trK[ijk];

      double P__U0U2 = 
AConfIJ_U0U2[ijk] - 0.66666666666666663*igConf_U0U2[ijk]*psi6_*
trK[ijk];

      double P__U2U2 = 
AConfIJ_U2U2[ijk] - 0.66666666666666663*igConf_U2U2[ijk]*psi6_*
trK[ijk];

      double G1_U1 = 
ChrisConf_U0D0D0[ijk]*P__U0U1 + ChrisConf_U0D0D1[ijk]*P__U1U1 + 
ChrisConf_U0D0D2[ijk]*P__U1U2 + ChrisConf_U1D0D0[ijk]*P__U0U0 + 3.0*
ChrisConf_U1D0D1[ijk]*P__U0U1 + 2.0*ChrisConf_U1D0D2[ijk]*P__U0U2 + 2.0*
ChrisConf_U1D1D1[ijk]*P__U1U1 + 3.0*ChrisConf_U1D1D2[ijk]*P__U1U2 + 
ChrisConf_U1D2D2[ijk]*P__U2U2 + ChrisConf_U2D0D2[ijk]*P__U0U1 + 
ChrisConf_U2D1D2[ijk]*P__U1U1 + ChrisConf_U2D2D2[ijk]*
P__U1U2;

      double G1_U0 = 
2.0*ChrisConf_U0D0D0[ijk]*P__U0U0 + 3.0*ChrisConf_U0D0D1[ijk]*P__U0U1 + 
3.0*ChrisConf_U0D0D2[ijk]*P__U0U2 + ChrisConf_U0D1D1[ijk]*P__U1U1 + 2.0*
ChrisConf_U0D1D2[ijk]*P__U1U2 + ChrisConf_U0D2D2[ijk]*P__U2U2 + 
ChrisConf_U1D0D1[ijk]*P__U0U0 + ChrisConf_U1D1D1[ijk]*P__U0U1 + 
ChrisConf_U1D1D2[ijk]*P__U0U2 + ChrisConf_U2D0D2[ijk]*P__U0U0 + 
ChrisConf_U2D1D2[ijk]*P__U0U1 + ChrisConf_U2D2D2[ijk]*
P__U0U2;

      double G1_U2 = 
ChrisConf_U0D0D0[ijk]*P__U0U2 + ChrisConf_U0D0D1[ijk]*P__U1U2 + 
ChrisConf_U0D0D2[ijk]*P__U2U2 + ChrisConf_U1D0D1[ijk]*P__U0U2 + 
ChrisConf_U1D1D1[ijk]*P__U1U2 + ChrisConf_U1D1D2[ijk]*P__U2U2 + 
ChrisConf_U2D0D0[ijk]*P__U0U0 + 2.0*ChrisConf_U2D0D1[ijk]*P__U0U1 + 3.0*
ChrisConf_U2D0D2[ijk]*P__U0U2 + ChrisConf_U2D1D1[ijk]*P__U1U1 + 3.0*
ChrisConf_U2D1D2[ijk]*P__U1U2 + 2.0*ChrisConf_U2D2D2[ijk]*
P__U2U2;

      double G2_U0 = 
-(2.0*P__U0U0*dpsi_D0[ijk]*gConf_D0D0[ijk]*igConf_U0U0[ijk] + 2.0*
P__U0U0*dpsi_D1[ijk]*gConf_D0D0[ijk]*igConf_U0U1[ijk] + 2.0*P__U0U0*
dpsi_D2[ijk]*gConf_D0D0[ijk]*igConf_U0U2[ijk] + 4.0*P__U0U1*
dpsi_D0[ijk]*gConf_D0D1[ijk]*igConf_U0U0[ijk] + 4.0*P__U0U1*
dpsi_D1[ijk]*gConf_D0D1[ijk]*igConf_U0U1[ijk] + 4.0*P__U0U1*
dpsi_D2[ijk]*gConf_D0D1[ijk]*igConf_U0U2[ijk] + 4.0*P__U0U2*
dpsi_D0[ijk]*gConf_D0D2[ijk]*igConf_U0U0[ijk] + 4.0*P__U0U2*
dpsi_D1[ijk]*gConf_D0D2[ijk]*igConf_U0U1[ijk] + 4.0*P__U0U2*
dpsi_D2[ijk]*gConf_D0D2[ijk]*igConf_U0U2[ijk] + 2.0*P__U1U1*
dpsi_D0[ijk]*gConf_D1D1[ijk]*igConf_U0U0[ijk] + 2.0*P__U1U1*
dpsi_D1[ijk]*gConf_D1D1[ijk]*igConf_U0U1[ijk] + 2.0*P__U1U1*
dpsi_D2[ijk]*gConf_D1D1[ijk]*igConf_U0U2[ijk] + 4.0*P__U1U2*
dpsi_D0[ijk]*gConf_D1D2[ijk]*igConf_U0U0[ijk] + 4.0*P__U1U2*
dpsi_D1[ijk]*gConf_D1D2[ijk]*igConf_U0U1[ijk] + 4.0*P__U1U2*
dpsi_D2[ijk]*gConf_D1D2[ijk]*igConf_U0U2[ijk] + 2.0*P__U2U2*
dpsi_D0[ijk]*gConf_D2D2[ijk]*igConf_U0U0[ijk] + 2.0*P__U2U2*
dpsi_D1[ijk]*gConf_D2D2[ijk]*igConf_U0U1[ijk] + 2.0*P__U2U2*
dpsi_D2[ijk]*gConf_D2D2[ijk]*igConf_U0U2[ijk])/psi[ijk];

      double G2_U1 = 
-(2.0*P__U0U0*dpsi_D0[ijk]*gConf_D0D0[ijk]*igConf_U0U1[ijk] + 2.0*
P__U0U0*dpsi_D1[ijk]*gConf_D0D0[ijk]*igConf_U1U1[ijk] + 2.0*P__U0U0*
dpsi_D2[ijk]*gConf_D0D0[ijk]*igConf_U1U2[ijk] + 4.0*P__U0U1*
dpsi_D0[ijk]*gConf_D0D1[ijk]*igConf_U0U1[ijk] + 4.0*P__U0U1*
dpsi_D1[ijk]*gConf_D0D1[ijk]*igConf_U1U1[ijk] + 4.0*P__U0U1*
dpsi_D2[ijk]*gConf_D0D1[ijk]*igConf_U1U2[ijk] + 4.0*P__U0U2*
dpsi_D0[ijk]*gConf_D0D2[ijk]*igConf_U0U1[ijk] + 4.0*P__U0U2*
dpsi_D1[ijk]*gConf_D0D2[ijk]*igConf_U1U1[ijk] + 4.0*P__U0U2*
dpsi_D2[ijk]*gConf_D0D2[ijk]*igConf_U1U2[ijk] + 2.0*P__U1U1*
dpsi_D0[ijk]*gConf_D1D1[ijk]*igConf_U0U1[ijk] + 2.0*P__U1U1*
dpsi_D1[ijk]*gConf_D1D1[ijk]*igConf_U1U1[ijk] + 2.0*P__U1U1*
dpsi_D2[ijk]*gConf_D1D1[ijk]*igConf_U1U2[ijk] + 4.0*P__U1U2*
dpsi_D0[ijk]*gConf_D1D2[ijk]*igConf_U0U1[ijk] + 4.0*P__U1U2*
dpsi_D1[ijk]*gConf_D1D2[ijk]*igConf_U1U1[ijk] + 4.0*P__U1U2*
dpsi_D2[ijk]*gConf_D1D2[ijk]*igConf_U1U2[ijk] + 2.0*P__U2U2*
dpsi_D0[ijk]*gConf_D2D2[ijk]*igConf_U0U1[ijk] + 2.0*P__U2U2*
dpsi_D1[ijk]*gConf_D2D2[ijk]*igConf_U1U1[ijk] + 2.0*P__U2U2*
dpsi_D2[ijk]*gConf_D2D2[ijk]*igConf_U1U2[ijk])/psi[ijk];

      double G2_U2 = 
-(2.0*P__U0U0*dpsi_D0[ijk]*gConf_D0D0[ijk]*igConf_U0U2[ijk] + 2.0*
P__U0U0*dpsi_D1[ijk]*gConf_D0D0[ijk]*igConf_U1U2[ijk] + 2.0*P__U0U0*
dpsi_D2[ijk]*gConf_D0D0[ijk]*igConf_U2U2[ijk] + 4.0*P__U0U1*
dpsi_D0[ijk]*gConf_D0D1[ijk]*igConf_U0U2[ijk] + 4.0*P__U0U1*
dpsi_D1[ijk]*gConf_D0D1[ijk]*igConf_U1U2[ijk] + 4.0*P__U0U1*
dpsi_D2[ijk]*gConf_D0D1[ijk]*igConf_U2U2[ijk] + 4.0*P__U0U2*
dpsi_D0[ijk]*gConf_D0D2[ijk]*igConf_U0U2[ijk] + 4.0*P__U0U2*
dpsi_D1[ijk]*gConf_D0D2[ijk]*igConf_U1U2[ijk] + 4.0*P__U0U2*
dpsi_D2[ijk]*gConf_D0D2[ijk]*igConf_U2U2[ijk] + 2.0*P__U1U1*
dpsi_D0[ijk]*gConf_D1D1[ijk]*igConf_U0U2[ijk] + 2.0*P__U1U1*
dpsi_D1[ijk]*gConf_D1D1[ijk]*igConf_U1U2[ijk] + 2.0*P__U1U1*
dpsi_D2[ijk]*gConf_D1D1[ijk]*igConf_U2U2[ijk] + 4.0*P__U1U2*
dpsi_D0[ijk]*gConf_D1D2[ijk]*igConf_U0U2[ijk] + 4.0*P__U1U2*
dpsi_D1[ijk]*gConf_D1D2[ijk]*igConf_U1U2[ijk] + 4.0*P__U1U2*
dpsi_D2[ijk]*gConf_D1D2[ijk]*igConf_U2U2[ijk] + 2.0*P__U2U2*
dpsi_D0[ijk]*gConf_D2D2[ijk]*igConf_U0U2[ijk] + 2.0*P__U2U2*
dpsi_D1[ijk]*gConf_D2D2[ijk]*igConf_U1U2[ijk] + 2.0*P__U2U2*
dpsi_D2[ijk]*gConf_D2D2[ijk]*igConf_U2U2[ijk])/psi[ijk];

      double G_U1 = 
att*(G1_U1 + G2_U1);

      double G_U0 = 
att*(G1_U0 + G2_U0);

      double G_U2 = 
att*(G1_U2 + G2_U2);


      /* populating: */
      obs__G_U0[ijk] = G_U0;
      obs__G_U1[ijk] = G_U1;
      obs__G_U2[ijk] = G_U2;
      }
    }


  }

  obs->ret[0] = obs_integral_SV
    (obs,"obs__P_U0","obs__G_U0",'+','-')/(8*M_PI);
  obs->ret[1] = obs_integral_SV
    (obs,"obs__P_U1","obs__G_U1",'+','-')/(8*M_PI);
  obs->ret[2] = obs_integral_SV
    (obs,"obs__P_U2","obs__G_U2",'+','-')/(8*M_PI);

  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = adm[p]->patch;
    remove_field_regex(patch,"^obs__P_U.$");
    remove_field_regex(patch,"^obs__G_U.$");
  }

}
