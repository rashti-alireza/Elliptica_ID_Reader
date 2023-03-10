/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "obs_header.h"



#define add_and_get_field(name) \
  if (_Ind(#name) < 0) \
  {ADD_AND_ALLOC_FIELD(name);} \
  WRITE_v(name);



double obs_Komar_mass(Observe_T *const obs);
double obs_Komar_mass(Observe_T *const obs)
{
  double Komar_mass = 0;
  struct items_S **const Komar = obs->items;
  const Uint N = obs->Nitems;
  Uint p;

  for(p = 0; p < N; ++p)
  {
    Patch_T *patch     = Komar[p]->patch;
    const double *n_U0 = Komar[p]->n_U0;
    const double *n_U1 = Komar[p]->n_U1;
    const double *n_U2 = Komar[p]->n_U2;


  /* declaring: */
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
  READ_v(gConf_D2D2)
  READ_v(AConfIJ_U0U1)
  READ_v(AConfIJ_U0U0)
  READ_v(AConfIJ_U0U2)
  READ_v(AConfIJ_U2U2)
  READ_v(AConfIJ_U1U1)
  READ_v(AConfIJ_U1U2)
  READ_v(alphaPsi)
  READ_v(dalphaPsi_D2)
  READ_v(dalphaPsi_D1)
  READ_v(dalphaPsi_D0)
  READ_v(beta_U1)
  READ_v(beta_U0)
  READ_v(beta_U2)
  READ_v(trK)
  READ_v(psi)
  READ_v(dpsi_D0)
  READ_v(dpsi_D1)
  READ_v(dpsi_D2)
  READ_v(EConf)
  READ_v(JConf_U0)
  READ_v(JConf_U1)
  READ_v(JConf_U2)
  READ_v(SConf)
  add_and_get_field(obs_komar_mass__s)
  add_and_get_field(obs_komar_mass__v)


  if (Komar[p]->surface_integration_flg)
  {

    FOR_ALL_ijk
    {
    double psim2 = 
pow(psi[ijk], -2);

    double psi2 = 
pow(psi[ijk], 2);

    double psi4 = 
pow(psi[ijk], 4);

    double dalpha_U1 = 
-alphaPsi[ijk]*dpsi_D1[ijk]/psi2 + dalphaPsi_D1[ijk]/
psi[ijk];

    double dalpha_U0 = 
-alphaPsi[ijk]*dpsi_D0[ijk]/psi2 + dalphaPsi_D0[ijk]/
psi[ijk];

    double dalpha_U2 = 
-alphaPsi[ijk]*dpsi_D2[ijk]/psi2 + dalphaPsi_D2[ijk]/
psi[ijk];

    double K_DD_D1D1 = 
0.33333333333333331*gConf_D1D1[ijk]*psi4*trK[ijk] + psim2*
(AConfIJ_U0U0[ijk]*pow(gConf_D0D1[ijk], 2) + 2.0*AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D1[ijk] + 2.0*AConfIJ_U0U2[ijk]*gConf_D0D1[ijk]*
gConf_D1D2[ijk] + AConfIJ_U1U1[ijk]*pow(gConf_D1D1[ijk], 2) + 2.0*
AConfIJ_U1U2[ijk]*gConf_D1D1[ijk]*gConf_D1D2[ijk] + AConfIJ_U2U2[ijk]*
pow(gConf_D1D2[ijk], 2));

    double K_DD_D1D2 = 
0.33333333333333331*gConf_D1D2[ijk]*psi4*trK[ijk] + psim2*
(AConfIJ_U0U0[ijk]*gConf_D0D1[ijk]*gConf_D0D2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D2[ijk] + AConfIJ_U0U1[ijk]*gConf_D0D2[ijk]*
gConf_D1D1[ijk] + AConfIJ_U0U2[ijk]*gConf_D0D1[ijk]*gConf_D2D2[ijk] + 
AConfIJ_U0U2[ijk]*gConf_D0D2[ijk]*gConf_D1D2[ijk] + AConfIJ_U1U1[ijk]*
gConf_D1D1[ijk]*gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*gConf_D1D1[ijk]*
gConf_D2D2[ijk] + AConfIJ_U1U2[ijk]*pow(gConf_D1D2[ijk], 2) + 
AConfIJ_U2U2[ijk]*gConf_D1D2[ijk]*gConf_D2D2[ijk]);

    double K_DD_D2D2 = 
0.33333333333333331*gConf_D2D2[ijk]*psi4*trK[ijk] + psim2*
(AConfIJ_U0U0[ijk]*pow(gConf_D0D2[ijk], 2) + 2.0*AConfIJ_U0U1[ijk]*
gConf_D0D2[ijk]*gConf_D1D2[ijk] + 2.0*AConfIJ_U0U2[ijk]*gConf_D0D2[ijk]*
gConf_D2D2[ijk] + AConfIJ_U1U1[ijk]*pow(gConf_D1D2[ijk], 2) + 2.0*
AConfIJ_U1U2[ijk]*gConf_D1D2[ijk]*gConf_D2D2[ijk] + AConfIJ_U2U2[ijk]*
pow(gConf_D2D2[ijk], 2));

    double K_DD_D0D0 = 
0.33333333333333331*gConf_D0D0[ijk]*psi4*trK[ijk] + psim2*
(AConfIJ_U0U0[ijk]*pow(gConf_D0D0[ijk], 2) + 2.0*AConfIJ_U0U1[ijk]*
gConf_D0D0[ijk]*gConf_D0D1[ijk] + 2.0*AConfIJ_U0U2[ijk]*gConf_D0D0[ijk]*
gConf_D0D2[ijk] + AConfIJ_U1U1[ijk]*pow(gConf_D0D1[ijk], 2) + 2.0*
AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*gConf_D0D2[ijk] + AConfIJ_U2U2[ijk]*
pow(gConf_D0D2[ijk], 2));

    double K_DD_D0D1 = 
0.33333333333333331*gConf_D0D1[ijk]*psi4*trK[ijk] + psim2*
(AConfIJ_U0U0[ijk]*gConf_D0D0[ijk]*gConf_D0D1[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D0[ijk]*gConf_D1D1[ijk] + AConfIJ_U0U1[ijk]*
pow(gConf_D0D1[ijk], 2) + AConfIJ_U0U2[ijk]*gConf_D0D0[ijk]*
gConf_D1D2[ijk] + AConfIJ_U0U2[ijk]*gConf_D0D1[ijk]*gConf_D0D2[ijk] + 
AConfIJ_U1U1[ijk]*gConf_D0D1[ijk]*gConf_D1D1[ijk] + AConfIJ_U1U2[ijk]*
gConf_D0D1[ijk]*gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]*
gConf_D1D1[ijk] + AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]*
gConf_D1D2[ijk]);

    double K_DD_D0D2 = 
0.33333333333333331*gConf_D0D2[ijk]*psi4*trK[ijk] + psim2*
(AConfIJ_U0U0[ijk]*gConf_D0D0[ijk]*gConf_D0D2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D0[ijk]*gConf_D1D2[ijk] + AConfIJ_U0U1[ijk]*gConf_D0D1[ijk]*
gConf_D0D2[ijk] + AConfIJ_U0U2[ijk]*gConf_D0D0[ijk]*gConf_D2D2[ijk] + 
AConfIJ_U0U2[ijk]*pow(gConf_D0D2[ijk], 2) + AConfIJ_U1U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D2[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*
gConf_D2D2[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]*gConf_D1D2[ijk] + 
AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]*gConf_D2D2[ijk]);

    double integrand_s = 
-n_U0[ijk]*(K_DD_D0D0*beta_U0[ijk] + K_DD_D0D1*beta_U1[ijk] + K_DD_D0D2*
beta_U2[ijk] - dalpha_U0) - n_U1[ijk]*(K_DD_D0D1*beta_U0[ijk] +
K_DD_D1D1*beta_U1[ijk] + K_DD_D1D2*beta_U2[ijk] - dalpha_U1) -
n_U2[ijk]*(K_DD_D0D2*beta_U0[ijk] + K_DD_D1D2*beta_U1[ijk] + K_DD_D2D2*
beta_U2[ijk] - dalpha_U2);


      obs_komar_mass__s[ijk] = integrand_s;
    }
    Integration_T *I = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->fields[Ind("obs_komar_mass__s")];
    I->g00 = Komar[p]->g00;
    I->g01 = Komar[p]->g01;
    I->g02 = Komar[p]->g02;
    I->g11 = Komar[p]->g11;
    I->g12 = Komar[p]->g12;
    I->g22 = Komar[p]->g22;
    I->Spectral->X_surface = Komar[p]->X_surface;
    I->Spectral->Y_surface = Komar[p]->Y_surface;
    I->Spectral->Z_surface = Komar[p]->Z_surface;
    I->Spectral->I         = Komar[p]->I;
    I->Spectral->J         = Komar[p]->J;
    I->Spectral->K         = Komar[p]->K;
    plan_integration(I);
    Komar_mass += execute_integration(I)/(4*M_PI);
    free_integration(I);
  }
  else
  {
  FOR_ALL_ijk
  {
  double psim6 = 
pow(psi[ijk], -6);

  double psi4_ = 
pow(psi[ijk], 4);

  double alpha = 
alphaPsi[ijk]/psi[ijk];

  double E = 
EConf[ijk]*psim6;

  double S = 
SConf[ijk]*psim6;

  double J_U2 = 
JConf_U2[ijk]*psim6;

  double J_U0 = 
JConf_U0[ijk]*psim6;

  double J_U1 = 
JConf_U1[ijk]*psim6;

  double integrand_v = 
alpha*(E + S) - 2.0*psi4_*(J_U0*beta_U0[ijk]*gConf_D0D0[ijk] + J_U0*
beta_U1[ijk]*gConf_D0D1[ijk] + J_U0*beta_U2[ijk]*gConf_D0D2[ijk] + J_U1*
beta_U0[ijk]*gConf_D0D1[ijk] + J_U1*beta_U1[ijk]*gConf_D1D1[ijk] + J_U1*
beta_U2[ijk]*gConf_D1D2[ijk] + J_U2*beta_U0[ijk]*gConf_D0D2[ijk] + J_U2*
beta_U1[ijk]*gConf_D1D2[ijk] + J_U2*beta_U2[ijk]*gConf_D2D2[ijk]);

   obs_komar_mass__v[ijk] = integrand_v;
  }
    Integration_T *I = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->fields[Ind("obs_komar_mass__v")];
    I->g00 = Komar[p]->g00;
    I->g01 = Komar[p]->g01;
    I->g02 = Komar[p]->g02;
    I->g11 = Komar[p]->g11;
    I->g12 = Komar[p]->g12;
    I->g22 = Komar[p]->g22;
    plan_integration(I);
    Komar_mass += execute_integration(I);
    free_integration(I);

  }

  REMOVE_FIELD(patch->fields[Ind("obs_komar_mass__s")]);
  REMOVE_FIELD(patch->fields[Ind("obs_komar_mass__v")]);
  }

  return Komar_mass;
}
