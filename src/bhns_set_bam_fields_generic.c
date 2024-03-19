/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2023 Alireza Rashti.
*/


#include "bhns_header.h"


#define add_alloc_get_field(name) ADD_FIELD(name) REALLOC_v_WRITE_v(name)


#define add_alloc_field(name) ADD_AND_ALLOC_FIELD(name) 


#define reav_v_if_exists(name) const double *name = 0; \
 if (_Ind(#name) >= 0) name = patch->fields[Ind(#name)]->v;


void bhns_set_bam_fields_generic(Grid_T *const grid);
void bhns_set_bam_fields_generic(Grid_T *const grid)
{
  Uint p;
  const Uint np = grid->np;

  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
  Patch_T *patch = grid->patch[p];

  add_alloc_get_field(bam_grhd_v_U0)
  add_alloc_get_field(bam_grhd_v_U1)
  add_alloc_get_field(bam_grhd_v_U2)
  add_alloc_get_field(bam_grhd_rho)
  add_alloc_get_field(bam_grhd_p)
  add_alloc_get_field(bam_grhd_epsl)
  add_alloc_get_field(bam_alpha)
  add_alloc_get_field(bam_beta_U0)
  add_alloc_get_field(bam_beta_U1)
  add_alloc_get_field(bam_beta_U2)
  add_alloc_get_field(bam_adm_g_D0D0)
  add_alloc_get_field(bam_adm_g_D0D1)
  add_alloc_get_field(bam_adm_g_D0D2)
  add_alloc_get_field(bam_adm_g_D1D1)
  add_alloc_get_field(bam_adm_g_D1D2)
  add_alloc_get_field(bam_adm_g_D2D2)
  add_alloc_get_field(bam_adm_Kij_D0D0)
  add_alloc_get_field(bam_adm_Kij_D0D1)
  add_alloc_get_field(bam_adm_Kij_D0D2)
  add_alloc_get_field(bam_adm_Kij_D1D1)
  add_alloc_get_field(bam_adm_Kij_D1D2)
  add_alloc_get_field(bam_adm_Kij_D2D2)
  reav_v_if_exists(gConf_D0D0)
  reav_v_if_exists(gConf_D0D1)
  reav_v_if_exists(gConf_D0D2)
  reav_v_if_exists(gConf_D1D1)
  reav_v_if_exists(gConf_D1D2)
  reav_v_if_exists(gConf_D2D2)
  reav_v_if_exists(igConf_U0U0)
  reav_v_if_exists(igConf_U0U1)
  reav_v_if_exists(igConf_U0U2)
  reav_v_if_exists(igConf_U1U1)
  reav_v_if_exists(igConf_U1U2)
  reav_v_if_exists(igConf_U2U2)
  reav_v_if_exists(adm_Kij_D0D0)
  reav_v_if_exists(adm_Kij_D0D1)
  reav_v_if_exists(adm_Kij_D0D2)
  reav_v_if_exists(adm_Kij_D1D1)
  reav_v_if_exists(adm_Kij_D1D2)
  reav_v_if_exists(adm_Kij_D2D2)
  reav_v_if_exists(beta_U0)
  reav_v_if_exists(beta_U1)
  reav_v_if_exists(beta_U2)
  reav_v_if_exists(psi)
  reav_v_if_exists(alphaPsi)
  reav_v_if_exists(enthalpy)
  reav_v_if_exists(W_U0)
  reav_v_if_exists(W_U1)
  reav_v_if_exists(W_U2)
  reav_v_if_exists(dphi_D0)
  reav_v_if_exists(dphi_D1)
  reav_v_if_exists(dphi_D2)
  reav_v_if_exists(u0)
  FOR_ALL_ijk
  {
  double psi4 =
pow(psi[ijk], 4);

   bam_alpha[ijk] = alphaPsi[ijk]/psi[ijk];
   bam_beta_U0[ijk] = beta_U0[ijk];
   bam_beta_U1[ijk] = beta_U1[ijk];
   bam_beta_U2[ijk] = beta_U2[ijk];
   double adm_g_D0D0 =
gConf_D0D0[ijk]*psi4;

   double adm_g_D0D1 =
gConf_D0D1[ijk]*psi4;

   double adm_g_D0D2 =
gConf_D0D2[ijk]*psi4;

   double adm_g_D1D1 =
gConf_D1D1[ijk]*psi4;

   double adm_g_D1D2 =
gConf_D1D2[ijk]*psi4;

   double adm_g_D2D2 =
gConf_D2D2[ijk]*psi4;


   bam_adm_g_D0D0[ijk] = adm_g_D0D0;
   bam_adm_g_D0D1[ijk] = adm_g_D0D1;
   bam_adm_g_D0D2[ijk] = adm_g_D0D2;
   bam_adm_g_D1D1[ijk] = adm_g_D1D1;
   bam_adm_g_D1D2[ijk] = adm_g_D1D2;
   bam_adm_g_D2D2[ijk] = adm_g_D2D2;
   double Kdd_D0D0 =
adm_Kij_D0D0[ijk];

   double Kdd_D0D1 =
adm_Kij_D0D1[ijk];

   double Kdd_D0D2 =
adm_Kij_D0D2[ijk];

   double Kdd_D1D1 =
adm_Kij_D1D1[ijk];

   double Kdd_D1D2 =
adm_Kij_D1D2[ijk];

   double Kdd_D2D2 =
adm_Kij_D2D2[ijk];


   bam_adm_Kij_D0D0[ijk] = Kdd_D0D0;
   bam_adm_Kij_D0D1[ijk] = Kdd_D0D1;
   bam_adm_Kij_D0D2[ijk] = Kdd_D0D2;
   bam_adm_Kij_D1D1[ijk] = Kdd_D1D1;
   bam_adm_Kij_D1D2[ijk] = Kdd_D1D2;
   bam_adm_Kij_D2D2[ijk] = Kdd_D2D2;
  }
  if (IsItCovering(patch,"NS"))
  {
  Physics_T *ns = init_physics(0,NS);
  EoS_T *eos    = init_EoS(ns);
  FOR_ALL_ijk
  {
  double psim4 =
pow(psi[ijk], -4);

  double grhd_v_U0 =
(W_U0[ijk] + psim4*(dphi_D0[ijk]*igConf_U0U0[ijk] + dphi_D1[ijk]*
igConf_U0U1[ijk] + dphi_D2[ijk]*igConf_U0U2[ijk]))/(bam_alpha[ijk]*
enthalpy[ijk]*u0[ijk]);

  double grhd_v_U1 =
(W_U1[ijk] + psim4*(dphi_D0[ijk]*igConf_U0U1[ijk] + dphi_D1[ijk]*
igConf_U1U1[ijk] + dphi_D2[ijk]*igConf_U1U2[ijk]))/(bam_alpha[ijk]*
enthalpy[ijk]*u0[ijk]);

  double grhd_v_U2 =
(W_U2[ijk] + psim4*(dphi_D0[ijk]*igConf_U0U2[ijk] + dphi_D1[ijk]*
igConf_U1U2[ijk] + dphi_D2[ijk]*igConf_U2U2[ijk]))/(bam_alpha[ijk]*
enthalpy[ijk]*u0[ijk]);


  bam_grhd_v_U0[ijk] = grhd_v_U0;
  bam_grhd_v_U1[ijk] = grhd_v_U1;
  bam_grhd_v_U2[ijk] = grhd_v_U2;
  eos->h = enthalpy[ijk];
  if(!isfinite(eos->h) || LSSEQL(eos->h,1.))
  {
    bam_grhd_rho[ijk]  = 0;
    bam_grhd_p[ijk]    = 0;
    bam_grhd_epsl[ijk] = 0;
  }
  else
  {
   bam_grhd_rho[ijk]  = eos->rest_mass_density(eos);
   bam_grhd_p[ijk]    = eos->pressure(eos);
   bam_grhd_epsl[ijk] = eos->specific_internal_energy(eos);
  }
  }
  free_physics(ns);
  free_EoS(eos);
  }
  }
}