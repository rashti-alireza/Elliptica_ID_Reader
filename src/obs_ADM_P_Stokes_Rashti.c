/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "obs_header.h"


#define add_and_get_field(name) \
  if (_Ind(#name) < 0) \
  {ADD_AND_ALLOC_FIELD(name);} \
  WRITE_v(name);


void obs_ADM_P_Stokes_SV_Rashti(Observe_T *const obs);
void obs_ADM_P_Stokes_SV_Rashti(Observe_T *const obs)
{
  struct items_S **adm = obs->items;
  const Uint N = obs->Nitems;
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
  READ_v_UNUSED(dAConfIJ_U0U0D2)
  READ_v_UNUSED(dAConfIJ_U0U0D0)
  READ_v_UNUSED(dAConfIJ_U0U0D1)
  READ_v_UNUSED(dAConfIJ_U0U1D2)
  READ_v_UNUSED(dAConfIJ_U0U1D1)
  READ_v_UNUSED(dAConfIJ_U0U1D0)
  READ_v_UNUSED(dAConfIJ_U0U2D0)
  READ_v_UNUSED(dAConfIJ_U0U2D1)
  READ_v_UNUSED(dAConfIJ_U0U2D2)
  READ_v_UNUSED(dAConfIJ_U2U2D2)
  READ_v_UNUSED(dAConfIJ_U1U2D2)
  READ_v_UNUSED(dAConfIJ_U2U2D0)
  READ_v_UNUSED(dAConfIJ_U2U2D1)
  READ_v_UNUSED(dAConfIJ_U1U2D1)
  READ_v_UNUSED(dAConfIJ_U1U1D2)
  READ_v_UNUSED(dAConfIJ_U1U2D0)
  READ_v_UNUSED(dAConfIJ_U1U1D0)
  READ_v_UNUSED(dAConfIJ_U1U1D1)
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
  READ_v(gConf_D2D2)
  READ_v(dgConf_D0D2D1)
  READ_v(dgConf_D0D2D2)
  READ_v(dgConf_D1D1D1)
  READ_v(dgConf_D1D2D1)
  READ_v(dgConf_D0D1D0)
  READ_v(dgConf_D0D1D1)
  READ_v(dgConf_D0D1D2)
  READ_v(dgConf_D0D2D0)
  READ_v(dgConf_D0D0D1)
  READ_v(dgConf_D0D0D0)
  READ_v(dgConf_D1D2D0)
  READ_v(dgConf_D0D0D2)
  READ_v(dgConf_D1D2D2)
  READ_v(dgConf_D1D1D0)
  READ_v(dgConf_D1D1D2)
  READ_v(dgConf_D2D2D1)
  READ_v(dgConf_D2D2D0)
  READ_v(dgConf_D2D2D2)
  READ_v(psi)
  READ_v(dpsi_D0)
  READ_v(dpsi_D1)
  READ_v(dpsi_D2)
  READ_v_UNUSED(ChrisConf_U1D0D0)
  READ_v_UNUSED(ChrisConf_U0D2D2)
  READ_v_UNUSED(ChrisConf_U1D0D2)
  READ_v_UNUSED(ChrisConf_U2D2D2)
  READ_v_UNUSED(ChrisConf_U2D1D1)
  READ_v_UNUSED(ChrisConf_U2D0D0)
  READ_v_UNUSED(ChrisConf_U1D1D1)
  READ_v_UNUSED(ChrisConf_U0D1D1)
  READ_v_UNUSED(ChrisConf_U0D1D2)
  READ_v_UNUSED(ChrisConf_U1D1D2)
  READ_v_UNUSED(ChrisConf_U0D0D1)
  READ_v_UNUSED(ChrisConf_U0D0D0)
  READ_v_UNUSED(ChrisConf_U2D0D1)
  READ_v_UNUSED(ChrisConf_U0D0D2)
  READ_v_UNUSED(ChrisConf_U2D1D2)
  READ_v_UNUSED(ChrisConf_U1D2D2)
  READ_v_UNUSED(ChrisConf_U2D0D2)
  READ_v_UNUSED(ChrisConf_U1D0D1)
  READ_v(trK)
  READ_v(dtrK_D2)
  READ_v(dtrK_D1)
  READ_v(dtrK_D0)
  add_and_get_field(obs__DPi_D2)
  add_and_get_field(obs__DPi_D1)
  add_and_get_field(obs__DPi_D0)
  add_and_get_field(obs__Pin_D2)
  add_and_get_field(obs__Pin_D1)
  add_and_get_field(obs__Pin_D0)



   if (adm[p]->surface_integration_flg)
   {
      const double *n_U0 = adm[p]->n_U0;
      const double *n_U1 = adm[p]->n_U1;
      const double *n_U2 = adm[p]->n_U2;
      for (ijk = 0; ijk < nn; ++ijk)
      {
      double psi4 = 
pow(psi[ijk], 4);

      double psim6 = 
pow(psi[ijk], -6);

      double Pi_U2D2 = 
psim6*(AConfIJ_U0U2[ijk]*gConf_D0D2[ijk] + AConfIJ_U1U2[ijk]*
gConf_D1D2[ijk] + AConfIJ_U2U2[ijk]*gConf_D2D2[ijk]) - 
0.66666666666666663*trK[ijk];

      double Pi_U2D0 = 
psim6*(AConfIJ_U0U2[ijk]*gConf_D0D0[ijk] + AConfIJ_U1U2[ijk]*
gConf_D0D1[ijk] + AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]);

      double Pi_U2D1 = 
psim6*(AConfIJ_U0U2[ijk]*gConf_D0D1[ijk] + AConfIJ_U1U2[ijk]*
gConf_D1D1[ijk] + AConfIJ_U2U2[ijk]*gConf_D1D2[ijk]);

      double Pi_U0D0 = 
psim6*(AConfIJ_U0U0[ijk]*gConf_D0D0[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk] + AConfIJ_U0U2[ijk]*gConf_D0D2[ijk]) - 
0.66666666666666663*trK[ijk];

      double Pi_U0D1 = 
psim6*(AConfIJ_U0U0[ijk]*gConf_D0D1[ijk] + AConfIJ_U0U1[ijk]*
gConf_D1D1[ijk] + AConfIJ_U0U2[ijk]*gConf_D1D2[ijk]);

      double Pi_U1D1 = 
psim6*(AConfIJ_U0U1[ijk]*gConf_D0D1[ijk] + AConfIJ_U1U1[ijk]*
gConf_D1D1[ijk] + AConfIJ_U1U2[ijk]*gConf_D1D2[ijk]) - 
0.66666666666666663*trK[ijk];

      double Pi_U1D0 = 
psim6*(AConfIJ_U0U1[ijk]*gConf_D0D0[ijk] + AConfIJ_U1U1[ijk]*
gConf_D0D1[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]);

      double Pi_U0D2 = 
psim6*(AConfIJ_U0U0[ijk]*gConf_D0D2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D1D2[ijk] + AConfIJ_U0U2[ijk]*gConf_D2D2[ijk]);

      double Pi_U1D2 = 
psim6*(AConfIJ_U0U1[ijk]*gConf_D0D2[ijk] + AConfIJ_U1U1[ijk]*
gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*gConf_D2D2[ijk]);

      double Pin_D2 = 
psi4*(Pi_U0D2*gConf_D0D0[ijk]*n_U0[ijk] + Pi_U0D2*gConf_D0D1[ijk]*
n_U1[ijk] + Pi_U0D2*gConf_D0D2[ijk]*n_U2[ijk] + Pi_U1D2*gConf_D0D1[ijk]*
n_U0[ijk] + Pi_U1D2*gConf_D1D1[ijk]*n_U1[ijk] + Pi_U1D2*gConf_D1D2[ijk]*
n_U2[ijk] + Pi_U2D2*gConf_D0D2[ijk]*n_U0[ijk] + Pi_U2D2*gConf_D1D2[ijk]*
n_U1[ijk] + Pi_U2D2*gConf_D2D2[ijk]*n_U2[ijk]);

      double Pin_D0 = 
psi4*(Pi_U0D0*gConf_D0D0[ijk]*n_U0[ijk] + Pi_U0D0*gConf_D0D1[ijk]*
n_U1[ijk] + Pi_U0D0*gConf_D0D2[ijk]*n_U2[ijk] + Pi_U1D0*gConf_D0D1[ijk]*
n_U0[ijk] + Pi_U1D0*gConf_D1D1[ijk]*n_U1[ijk] + Pi_U1D0*gConf_D1D2[ijk]*
n_U2[ijk] + Pi_U2D0*gConf_D0D2[ijk]*n_U0[ijk] + Pi_U2D0*gConf_D1D2[ijk]*
n_U1[ijk] + Pi_U2D0*gConf_D2D2[ijk]*n_U2[ijk]);

      double Pin_D1 = 
psi4*(Pi_U0D1*gConf_D0D0[ijk]*n_U0[ijk] + Pi_U0D1*gConf_D0D1[ijk]*
n_U1[ijk] + Pi_U0D1*gConf_D0D2[ijk]*n_U2[ijk] + Pi_U1D1*gConf_D0D1[ijk]*
n_U0[ijk] + Pi_U1D1*gConf_D1D1[ijk]*n_U1[ijk] + Pi_U1D1*gConf_D1D2[ijk]*
n_U2[ijk] + Pi_U2D1*gConf_D0D2[ijk]*n_U0[ijk] + Pi_U2D1*gConf_D1D2[ijk]*
n_U1[ijk] + Pi_U2D1*gConf_D2D2[ijk]*n_U2[ijk]);


      /* populating: */
      obs__Pin_D2[ijk] = Pin_D2;
      obs__Pin_D1[ijk] = Pin_D1;
      obs__Pin_D0[ijk] = Pin_D0;
      }
   }
    else
    {
      for (ijk = 0; ijk < nn; ++ijk)
      {
      double psim6_ = 
pow(psi[ijk], -6);

      double Pi__U0D0 = 
psim6_*(AConfIJ_U0U0[ijk]*gConf_D0D0[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk] + AConfIJ_U0U2[ijk]*gConf_D0D2[ijk]) - 
0.66666666666666663*trK[ijk];

      double Pi__U0D1 = 
psim6_*(AConfIJ_U0U0[ijk]*gConf_D0D1[ijk] + AConfIJ_U0U1[ijk]*
gConf_D1D1[ijk] + AConfIJ_U0U2[ijk]*gConf_D1D2[ijk]);

      double Pi__U0D2 = 
psim6_*(AConfIJ_U0U0[ijk]*gConf_D0D2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D1D2[ijk] + AConfIJ_U0U2[ijk]*gConf_D2D2[ijk]);

      double Pi__U1D2 = 
psim6_*(AConfIJ_U0U1[ijk]*gConf_D0D2[ijk] + AConfIJ_U1U1[ijk]*
gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*gConf_D2D2[ijk]);

      double Pi__U1D1 = 
psim6_*(AConfIJ_U0U1[ijk]*gConf_D0D1[ijk] + AConfIJ_U1U1[ijk]*
gConf_D1D1[ijk] + AConfIJ_U1U2[ijk]*gConf_D1D2[ijk]) - 
0.66666666666666663*trK[ijk];

      double Pi__U2D2 = 
psim6_*(AConfIJ_U0U2[ijk]*gConf_D0D2[ijk] + AConfIJ_U1U2[ijk]*
gConf_D1D2[ijk] + AConfIJ_U2U2[ijk]*gConf_D2D2[ijk]) - 
0.66666666666666663*trK[ijk];

      double Pi__U1D0 = 
psim6_*(AConfIJ_U0U1[ijk]*gConf_D0D0[ijk] + AConfIJ_U1U1[ijk]*
gConf_D0D1[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]);

      double Pi__U2D0 = 
psim6_*(AConfIJ_U0U2[ijk]*gConf_D0D0[ijk] + AConfIJ_U1U2[ijk]*
gConf_D0D1[ijk] + AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]);

      double Pi__U2D1 = 
psim6_*(AConfIJ_U0U2[ijk]*gConf_D0D1[ijk] + AConfIJ_U1U2[ijk]*
gConf_D1D1[ijk] + AConfIJ_U2U2[ijk]*gConf_D1D2[ijk]);

      double dlnpsi_U2 = 
dpsi_D2[ijk]/psi[ijk];

      double dlnpsi_U1 = 
dpsi_D1[ijk]/psi[ijk];

      double dlnpsi_U0 = 
dpsi_D0[ijk]/psi[ijk];

      double dPi1_U1 = 
-6.0*psim6_*(AConfIJ_U0U0[ijk]*dlnpsi_U0*gConf_D0D1[ijk] + 
AConfIJ_U0U1[ijk]*dlnpsi_U0*gConf_D1D1[ijk] + AConfIJ_U0U1[ijk]*
dlnpsi_U1*gConf_D0D1[ijk] + AConfIJ_U0U2[ijk]*dlnpsi_U0*
gConf_D1D2[ijk] + AConfIJ_U0U2[ijk]*dlnpsi_U2*gConf_D0D1[ijk] + 
AConfIJ_U1U1[ijk]*dlnpsi_U1*gConf_D1D1[ijk] + AConfIJ_U1U2[ijk]*
dlnpsi_U1*gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*dlnpsi_U2*
gConf_D1D1[ijk] + AConfIJ_U2U2[ijk]*dlnpsi_U2*gConf_D1D2[ijk]);

      double dPi1_U0 = 
-6.0*psim6_*(AConfIJ_U0U0[ijk]*dlnpsi_U0*gConf_D0D0[ijk] + 
AConfIJ_U0U1[ijk]*dlnpsi_U0*gConf_D0D1[ijk] + AConfIJ_U0U1[ijk]*
dlnpsi_U1*gConf_D0D0[ijk] + AConfIJ_U0U2[ijk]*dlnpsi_U0*
gConf_D0D2[ijk] + AConfIJ_U0U2[ijk]*dlnpsi_U2*gConf_D0D0[ijk] + 
AConfIJ_U1U1[ijk]*dlnpsi_U1*gConf_D0D1[ijk] + AConfIJ_U1U2[ijk]*
dlnpsi_U1*gConf_D0D2[ijk] + AConfIJ_U1U2[ijk]*dlnpsi_U2*
gConf_D0D1[ijk] + AConfIJ_U2U2[ijk]*dlnpsi_U2*gConf_D0D2[ijk]);

      double dPi1_U2 = 
-6.0*psim6_*(AConfIJ_U0U0[ijk]*dlnpsi_U0*gConf_D0D2[ijk] + 
AConfIJ_U0U1[ijk]*dlnpsi_U0*gConf_D1D2[ijk] + AConfIJ_U0U1[ijk]*
dlnpsi_U1*gConf_D0D2[ijk] + AConfIJ_U0U2[ijk]*dlnpsi_U0*
gConf_D2D2[ijk] + AConfIJ_U0U2[ijk]*dlnpsi_U2*gConf_D0D2[ijk] + 
AConfIJ_U1U1[ijk]*dlnpsi_U1*gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*
dlnpsi_U1*gConf_D2D2[ijk] + AConfIJ_U1U2[ijk]*dlnpsi_U2*
gConf_D1D2[ijk] + AConfIJ_U2U2[ijk]*dlnpsi_U2*gConf_D2D2[ijk]);

      double dPi2_U0 = 
psim6_*(dAConfIJ_U0U0D0[ijk]*gConf_D0D0[ijk] + dAConfIJ_U0U1D0[ijk]*
gConf_D0D1[ijk] + dAConfIJ_U0U1D1[ijk]*gConf_D0D0[ijk] + 
dAConfIJ_U0U2D0[ijk]*gConf_D0D2[ijk] + dAConfIJ_U0U2D2[ijk]*
gConf_D0D0[ijk] + dAConfIJ_U1U1D1[ijk]*gConf_D0D1[ijk] + 
dAConfIJ_U1U2D1[ijk]*gConf_D0D2[ijk] + dAConfIJ_U1U2D2[ijk]*
gConf_D0D1[ijk] + dAConfIJ_U2U2D2[ijk]*gConf_D0D2[ijk]);

      double dPi2_U1 = 
psim6_*(dAConfIJ_U0U0D0[ijk]*gConf_D0D1[ijk] + dAConfIJ_U0U1D0[ijk]*
gConf_D1D1[ijk] + dAConfIJ_U0U1D1[ijk]*gConf_D0D1[ijk] + 
dAConfIJ_U0U2D0[ijk]*gConf_D1D2[ijk] + dAConfIJ_U0U2D2[ijk]*
gConf_D0D1[ijk] + dAConfIJ_U1U1D1[ijk]*gConf_D1D1[ijk] + 
dAConfIJ_U1U2D1[ijk]*gConf_D1D2[ijk] + dAConfIJ_U1U2D2[ijk]*
gConf_D1D1[ijk] + dAConfIJ_U2U2D2[ijk]*gConf_D1D2[ijk]);

      double dPi2_U2 = 
psim6_*(dAConfIJ_U0U0D0[ijk]*gConf_D0D2[ijk] + dAConfIJ_U0U1D0[ijk]*
gConf_D1D2[ijk] + dAConfIJ_U0U1D1[ijk]*gConf_D0D2[ijk] + 
dAConfIJ_U0U2D0[ijk]*gConf_D2D2[ijk] + dAConfIJ_U0U2D2[ijk]*
gConf_D0D2[ijk] + dAConfIJ_U1U1D1[ijk]*gConf_D1D2[ijk] + 
dAConfIJ_U1U2D1[ijk]*gConf_D2D2[ijk] + dAConfIJ_U1U2D2[ijk]*
gConf_D1D2[ijk] + dAConfIJ_U2U2D2[ijk]*gConf_D2D2[ijk]);

      double dPi3_U2 = 
-0.66666666666666663*dtrK_D2[ijk] + psim6_*(AConfIJ_U0U0[ijk]*
dgConf_D0D2D0[ijk] + AConfIJ_U0U1[ijk]*dgConf_D0D2D1[ijk] + 
AConfIJ_U0U1[ijk]*dgConf_D1D2D0[ijk] + AConfIJ_U0U2[ijk]*
dgConf_D0D2D2[ijk] + AConfIJ_U0U2[ijk]*dgConf_D2D2D0[ijk] + 
AConfIJ_U1U1[ijk]*dgConf_D1D2D1[ijk] + AConfIJ_U1U2[ijk]*
dgConf_D1D2D2[ijk] + AConfIJ_U1U2[ijk]*dgConf_D2D2D1[ijk] + 
AConfIJ_U2U2[ijk]*dgConf_D2D2D2[ijk]);

      double dPi3_U1 = 
-0.66666666666666663*dtrK_D1[ijk] + psim6_*(AConfIJ_U0U0[ijk]*
dgConf_D0D1D0[ijk] + AConfIJ_U0U1[ijk]*dgConf_D0D1D1[ijk] + 
AConfIJ_U0U1[ijk]*dgConf_D1D1D0[ijk] + AConfIJ_U0U2[ijk]*
dgConf_D0D1D2[ijk] + AConfIJ_U0U2[ijk]*dgConf_D1D2D0[ijk] + 
AConfIJ_U1U1[ijk]*dgConf_D1D1D1[ijk] + AConfIJ_U1U2[ijk]*
dgConf_D1D1D2[ijk] + AConfIJ_U1U2[ijk]*dgConf_D1D2D1[ijk] + 
AConfIJ_U2U2[ijk]*dgConf_D1D2D2[ijk]);

      double dPi3_U0 = 
-0.66666666666666663*dtrK_D0[ijk] + psim6_*(AConfIJ_U0U0[ijk]*
dgConf_D0D0D0[ijk] + AConfIJ_U0U1[ijk]*dgConf_D0D0D1[ijk] + 
AConfIJ_U0U1[ijk]*dgConf_D0D1D0[ijk] + AConfIJ_U0U2[ijk]*
dgConf_D0D0D2[ijk] + AConfIJ_U0U2[ijk]*dgConf_D0D2D0[ijk] + 
AConfIJ_U1U1[ijk]*dgConf_D0D1D1[ijk] + AConfIJ_U1U2[ijk]*
dgConf_D0D1D2[ijk] + AConfIJ_U1U2[ijk]*dgConf_D0D2D1[ijk] + 
AConfIJ_U2U2[ijk]*dgConf_D0D2D2[ijk]);

      double dPi_U1 = 
dPi1_U1 + dPi2_U1 + dPi3_U1;

      double dPi_U0 = 
dPi1_U0 + dPi2_U0 + dPi3_U0;

      double dPi_U2 = 
dPi1_U2 + dPi2_U2 + dPi3_U2;

      double G_U1 = 
ChrisConf_U0D0D1[ijk] + ChrisConf_U1D1D1[ijk] + ChrisConf_U2D1D2[ijk] + 
6*dlnpsi_U1;

      double G_U0 = 
ChrisConf_U0D0D0[ijk] + ChrisConf_U1D0D1[ijk] + ChrisConf_U2D0D2[ijk] + 
6*dlnpsi_U0;

      double G_U2 = 
ChrisConf_U0D0D2[ijk] + ChrisConf_U1D1D2[ijk] + ChrisConf_U2D2D2[ijk] + 
6*dlnpsi_U2;

      double DPi_D2 = 
G_U0*Pi__U0D2 + G_U1*Pi__U1D2 + G_U2*Pi__U2D2 + dPi_U2;

      double DPi_D0 = 
G_U0*Pi__U0D0 + G_U1*Pi__U1D0 + G_U2*Pi__U2D0 + dPi_U0;

      double DPi_D1 = 
G_U0*Pi__U0D1 + G_U1*Pi__U1D1 + G_U2*Pi__U2D1 + dPi_U1;


      /* populating: */
      obs__DPi_D2[ijk] = DPi_D2;
      obs__DPi_D1[ijk] = DPi_D1;
      obs__DPi_D0[ijk] = DPi_D0;
      }
    }


  }

  obs->ret[0] = obs_integral_SV
                 (obs,"obs__Pin_D0","obs__DPi_D0",'+','+')/(8*M_PI);
  obs->ret[1] = obs_integral_SV
                 (obs,"obs__Pin_D1","obs__DPi_D1",'+','+')/(8*M_PI);
  obs->ret[2] = obs_integral_SV
                 (obs,"obs__Pin_D2","obs__DPi_D2",'+','+')/(8*M_PI);

  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = adm[p]->patch;
    remove_field_regex(patch,"^obs__DPi_D.$");
    remove_field_regex(patch,"^obs__Pin_D.$");
  }

}