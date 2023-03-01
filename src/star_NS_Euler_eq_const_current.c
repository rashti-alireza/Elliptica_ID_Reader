/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2021 Alireza Rashti.
*/


#include "star_header.h"
#include "physics_stress_energy_lib.h"


double star_NS_current_Euler_eq_const(Physics_T *const phys);
double star_NS_current_Euler_eq_const(Physics_T *const phys)
{
  AssureType(phys->type == NS)
  Grid_T *const grid = mygrid(phys,"NS");
  char cover[99] = {'\0'};

  sprintf(cover,"%s_%s",phys->spos,"central_box");
  FOR_ALL_p(grid->np)
  {
  Patch_T *patch = grid->patch[p];
  if (!IsItCovering(patch,cover))
    continue;

  Uint ijk = i_j_k_to_ijk(patch->n,patch->n[0]/2,patch->n[1]/2,patch->n[2]/2);

  Tij_NS_IF_XCTS_gConf_u0(patch);


  /* declaring: */
  READ_v(enthalpy)
  READ_v(igConf_U2U2)
  READ_v(igConf_U1U2)
  READ_v(igConf_U1U1)
  READ_v(igConf_U0U2)
  READ_v(igConf_U0U0)
  READ_v(igConf_U0U1)
  READ_v(W_U1)
  READ_v(W_U0)
  READ_v(W_U2)
  READ_v(dphi_D2)
  READ_v(dphi_D1)
  READ_v(dphi_D0)
  READ_v(beta_U1)
  READ_v(beta_U0)
  READ_v(beta_U2)
  READ_v(psi)
  READ_v(u0)


  double psim4 = 
pow(psi[ijk], -4);

  double dphiP = 
W_U0[ijk]*dphi_D0[ijk] + W_U1[ijk]*dphi_D1[ijk] + W_U2[ijk]*
dphi_D2[ijk] + psim4*(pow(dphi_D0[ijk], 2)*igConf_U0U0[ijk] + 2.0*
dphi_D0[ijk]*dphi_D1[ijk]*igConf_U0U1[ijk] + 2.0*dphi_D0[ijk]*
dphi_D2[ijk]*igConf_U0U2[ijk] + pow(dphi_D1[ijk], 2)*igConf_U1U1[ijk] +
2.0*dphi_D1[ijk]*dphi_D2[ijk]*igConf_U1U2[ijk] + pow(dphi_D2[ijk], 2)*
igConf_U2U2[ijk]);

  double Euler_C = 
(-dphiP - pow(enthalpy[ijk], 2) + enthalpy[ijk]*u0[ijk]*(beta_U0[ijk]*
dphi_D0[ijk] + beta_U1[ijk]*dphi_D1[ijk] + beta_U2[ijk]*dphi_D2[ijk]))/
(enthalpy[ijk]*u0[ijk]);

  return Euler_C;
  }
  Error0("Could not find NS central patch!");

  return DBL_MAX;
}
