/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2022 Alireza Rashti.
*/


#include "Tij_header.h"


void Tij_NS_IF_XCTS_gConf_psi6S(Patch_T *const patch,EoS_T *const eos)
{
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D0D2)
  READ_v(gConf_D1D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D2D2)
  READ_v(igConf_U0U0)
  READ_v(igConf_U0U1)
  READ_v(igConf_U0U2)
  READ_v(igConf_U1U1)
  READ_v(igConf_U1U2)
  READ_v(igConf_U2U2)
  READ_v(enthalpy)
  READ_v(W_U0)
  READ_v(W_U1)
  READ_v(W_U2)
  READ_v(dphi_D0)
  READ_v(dphi_D1)
  READ_v(dphi_D2)
  READ_v(psi)
  READ_v(rho0)
  REALLOC_v_WRITE_v(SConf)
  REALLOC_v_WRITE_v(SConfP)
  REALLOC_v_WRITE_v(SConfC)
FOR_ALL_ijk
{
  eos->h    = enthalpy[ijk];
  double p  = eos->pressure(eos);
  double e0 = eos->specific_internal_energy(eos);
  double psim4 =
pow(psi[ijk], -4);

  double psi4 =
pow(psi[ijk], 4);

  double psi6 =
pow(psi[ijk], 6);

  double P2 =
2.0*W_U0[ijk]*dphi_D0[ijk] + 2.0*W_U1[ijk]*dphi_D1[ijk] + 2.0*W_U2[ijk]*
dphi_D2[ijk] + psi4*(pow(W_U0[ijk], 2)*gConf_D0D0[ijk] + 2.0*W_U0[ijk]*
W_U1[ijk]*gConf_D0D1[ijk] + 2.0*W_U0[ijk]*W_U2[ijk]*gConf_D0D2[ijk] +
pow(W_U1[ijk], 2)*gConf_D1D1[ijk] + 2.0*W_U1[ijk]*W_U2[ijk]*
gConf_D1D2[ijk] + pow(W_U2[ijk], 2)*gConf_D2D2[ijk]) + psim4*
(pow(dphi_D0[ijk], 2)*igConf_U0U0[ijk] + 2.0*dphi_D0[ijk]*dphi_D1[ijk]*
igConf_U0U1[ijk] + 2.0*dphi_D0[ijk]*dphi_D2[ijk]*igConf_U0U2[ijk] +
pow(dphi_D1[ijk], 2)*igConf_U1U1[ijk] + 2.0*dphi_D1[ijk]*dphi_D2[ijk]*
igConf_U1U2[ijk] + pow(dphi_D2[ijk], 2)*igConf_U2U2[ijk]);

  double Sbar =
psi6*(P2*rho0[ijk]/enthalpy[ijk] + 3*p);

  SConf[ijk] = Sbar;
  double p_o_rho0 = enthalpy[ijk] - 1. - e0;
  double SbarP =
p*psi6*(P2/enthalpy[ijk] + 3*p_o_rho0);

  SConfP[ijk] = SbarP;
  SConfC[ijk] = p_o_rho0;
}
}