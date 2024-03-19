/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "adm_header.h"


void adm_add_3plus1_fields(Grid_T *const grid);
void adm_add_3plus1_fields(Grid_T *const grid)
{
 Uint p;
 FOR_ALL_PATCHES(p,grid)
 {
 Patch_T *patch = grid->patch[p];


  /* declaring: */
  ADD_FIELD(adm_g_D0D2)
  ADD_FIELD(adm_g_D0D1)
  ADD_FIELD(adm_g_D0D0)
  ADD_FIELD(adm_g_D1D2)
  ADD_FIELD(adm_g_D1D1)
  ADD_FIELD(adm_g_D2D2)
  ADD_FIELD(adm_Kij_D2D2)
  ADD_FIELD(adm_Kij_D0D2)
  ADD_FIELD(adm_Kij_D1D2)
  ADD_FIELD(adm_Kij_D0D0)
  ADD_FIELD(adm_Kij_D0D1)
  ADD_FIELD(adm_Kij_D1D1)
  ADD_FIELD(adm_KIJ_U0U2)
  ADD_FIELD(adm_KIJ_U0U0)
  ADD_FIELD(adm_KIJ_U0U1)
  ADD_FIELD(adm_KIJ_U1U2)
  ADD_FIELD(adm_KIJ_U1U1)
  ADD_FIELD(adm_KIJ_U2U2)
  ADD_FIELD(AConfIJ_U0U1)
  ADD_FIELD(AConfIJ_U0U0)
  ADD_FIELD(AConfIJ_U0U2)
  ADD_FIELD(AConfIJ_U2U2)
  ADD_FIELD(AConfIJ_U1U1)
  ADD_FIELD(AConfIJ_U1U2)
  ADD_FIELD(dAConfIJ_U0U0D2)
  ADD_FIELD(dAConfIJ_U0U0D0)
  ADD_FIELD(dAConfIJ_U0U0D1)
  ADD_FIELD(dAConfIJ_U0U1D2)
  ADD_FIELD(dAConfIJ_U0U1D1)
  ADD_FIELD(dAConfIJ_U0U1D0)
  ADD_FIELD(dAConfIJ_U0U2D0)
  ADD_FIELD(dAConfIJ_U0U2D1)
  ADD_FIELD(dAConfIJ_U0U2D2)
  ADD_FIELD(dAConfIJ_U2U2D2)
  ADD_FIELD(dAConfIJ_U1U2D2)
  ADD_FIELD(dAConfIJ_U2U2D0)
  ADD_FIELD(dAConfIJ_U2U2D1)
  ADD_FIELD(dAConfIJ_U1U2D1)
  ADD_FIELD(dAConfIJ_U1U1D2)
  ADD_FIELD(dAConfIJ_U1U2D0)
  ADD_FIELD(dAConfIJ_U1U1D0)
  ADD_FIELD(dAConfIJ_U1U1D1)
  ADD_FIELD(AConfIJ2)
  ADD_FIELD(ham1)
  ADD_FIELD(ham2)
  ADD_FIELD(mom1_U1)
  ADD_FIELD(mom1_U0)
  ADD_FIELD(mom1_U2)
  ADD_FIELD(mom2_U0)
  ADD_FIELD(mom2_U1)
  ADD_FIELD(mom2_U2)
  ADD_FIELD(beta_U1)
  ADD_FIELD(beta_U0)
  ADD_FIELD(beta_U2)
  ADD_FIELD(dbeta_U1D0)
  ADD_FIELD(dbeta_U2D2)
  ADD_FIELD(dbeta_U2D0)
  ADD_FIELD(dbeta_U2D1)
  ADD_FIELD(dbeta_U0D0)
  ADD_FIELD(dbeta_U1D1)
  ADD_FIELD(dbeta_U0D1)
  ADD_FIELD(dbeta_U0D2)
  ADD_FIELD(dbeta_U1D2)
  ADD_FIELD(ddbeta_U1D0D2)
  ADD_FIELD(ddbeta_U1D1D2)
  ADD_FIELD(ddbeta_U1D0D0)
  ADD_FIELD(ddbeta_U1D0D1)
  ADD_FIELD(ddbeta_U1D2D2)
  ADD_FIELD(ddbeta_U2D0D1)
  ADD_FIELD(ddbeta_U2D0D0)
  ADD_FIELD(ddbeta_U2D1D1)
  ADD_FIELD(ddbeta_U2D1D2)
  ADD_FIELD(ddbeta_U0D2D2)
  ADD_FIELD(ddbeta_U2D0D2)
  ADD_FIELD(ddbeta_U1D1D1)
  ADD_FIELD(ddbeta_U0D1D2)
  ADD_FIELD(ddbeta_U0D0D2)
  ADD_FIELD(ddbeta_U0D0D1)
  ADD_FIELD(ddbeta_U0D0D0)
  ADD_FIELD(ddbeta_U0D1D1)
  ADD_FIELD(ddbeta_U2D2D2)
  ADD_FIELD(B0_U2)
  ADD_FIELD(B0_U1)
  ADD_FIELD(B0_U0)
  ADD_FIELD(dB0_U2D2)
  ADD_FIELD(dB0_U2D0)
  ADD_FIELD(dB0_U0D0)
  ADD_FIELD(dB0_U2D1)
  ADD_FIELD(dB0_U0D1)
  ADD_FIELD(dB0_U0D2)
  ADD_FIELD(dB0_U1D2)
  ADD_FIELD(dB0_U1D1)
  ADD_FIELD(dB0_U1D0)
  ADD_FIELD(ddB0_U1D2D2)
  ADD_FIELD(ddB0_U1D0D0)
  ADD_FIELD(ddB0_U1D0D1)
  ADD_FIELD(ddB0_U1D0D2)
  ADD_FIELD(ddB0_U1D1D2)
  ADD_FIELD(ddB0_U2D0D0)
  ADD_FIELD(ddB0_U2D1D2)
  ADD_FIELD(ddB0_U0D2D2)
  ADD_FIELD(ddB0_U2D1D1)
  ADD_FIELD(ddB0_U0D0D1)
  ADD_FIELD(ddB0_U0D1D1)
  ADD_FIELD(ddB0_U0D1D2)
  ADD_FIELD(ddB0_U0D0D2)
  ADD_FIELD(ddB0_U1D1D1)
  ADD_FIELD(ddB0_U2D0D1)
  ADD_FIELD(ddB0_U2D2D2)
  ADD_FIELD(ddB0_U2D0D2)
  ADD_FIELD(ddB0_U0D0D0)
  ADD_FIELD(B1_U0)
  ADD_FIELD(B1_U1)
  ADD_FIELD(B1_U2)
  ADD_FIELD(dB1_U2D2)
  ADD_FIELD(dB1_U2D1)
  ADD_FIELD(dB1_U2D0)
  ADD_FIELD(dB1_U0D1)
  ADD_FIELD(dB1_U0D0)
  ADD_FIELD(dB1_U0D2)
  ADD_FIELD(dB1_U1D0)
  ADD_FIELD(dB1_U1D1)
  ADD_FIELD(dB1_U1D2)
  ADD_FIELD(ddB1_U2D1D1)
  ADD_FIELD(ddB1_U1D0D2)
  ADD_FIELD(ddB1_U0D2D2)
  ADD_FIELD(ddB1_U1D0D0)
  ADD_FIELD(ddB1_U2D1D2)
  ADD_FIELD(ddB1_U2D2D2)
  ADD_FIELD(ddB1_U1D0D1)
  ADD_FIELD(ddB1_U0D0D2)
  ADD_FIELD(ddB1_U2D0D1)
  ADD_FIELD(ddB1_U0D0D0)
  ADD_FIELD(ddB1_U0D0D1)
  ADD_FIELD(ddB1_U1D2D2)
  ADD_FIELD(ddB1_U2D0D0)
  ADD_FIELD(ddB1_U1D1D2)
  ADD_FIELD(ddB1_U0D1D2)
  ADD_FIELD(ddB1_U0D1D1)
  ADD_FIELD(ddB1_U1D1D1)
  ADD_FIELD(ddB1_U2D0D2)
  ADD_FIELD(alpha)
  ADD_FIELD(alphaPsi)
  ADD_FIELD(dalphaPsi_D2)
  ADD_FIELD(dalphaPsi_D1)
  ADD_FIELD(dalphaPsi_D0)
  ADD_FIELD(ddalphaPsi_D0D0)
  ADD_FIELD(ddalphaPsi_D0D1)
  ADD_FIELD(ddalphaPsi_D0D2)
  ADD_FIELD(ddalphaPsi_D1D1)
  ADD_FIELD(ddalphaPsi_D1D2)
  ADD_FIELD(ddalphaPsi_D2D2)
  ADD_FIELD(psi)
  ADD_FIELD(dpsi_D0)
  ADD_FIELD(dpsi_D1)
  ADD_FIELD(dpsi_D2)
  ADD_FIELD(ddpsi_D2D2)
  ADD_FIELD(ddpsi_D0D0)
  ADD_FIELD(ddpsi_D0D1)
  ADD_FIELD(ddpsi_D1D1)
  ADD_FIELD(ddpsi_D1D2)
  ADD_FIELD(ddpsi_D0D2)


 }
}