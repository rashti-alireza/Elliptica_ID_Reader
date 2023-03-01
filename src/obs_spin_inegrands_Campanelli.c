/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "obs_header.h"

#define declare_and_alloc_xi(name) \
  double *name = alloc_double(nn);

#define add_and_get_field(name) \
  if (_Ind(#name) < 0) \
  {ADD_AND_ALLOC_FIELD(name);} \
  WRITE_v(name);


void obs_populate_spin_integrands_Campanelli(Patch_T *const patch,const double xc[3],const double *const normal[3]);
void obs_populate_spin_integrands_Campanelli(Patch_T *const patch,const double xc[3],const double *const normal[3])
{
  const double x_c = xc[0];
  const double y_c = xc[1];
  const double z_c = xc[2];
  const Uint nn = patch->nn;
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
  READ_v(psi)
  add_and_get_field(SPIN_integrand_U0)
  add_and_get_field(SPIN_integrand_U1)
  add_and_get_field(SPIN_integrand_U2)
  declare_and_alloc_xi(xi_U2)
  declare_and_alloc_xi(xi_U0)
  declare_and_alloc_xi(xi_U1)



   for(ijk = 0; ijk < nn; ++ijk)
   {
   double x    = patch->node[ijk]->x[0];
   double y    = patch->node[ijk]->x[1];
   double z    = patch->node[ijk]->x[2];
   xi_U0[ijk] = x-x_c;
   xi_U1[ijk] = y-y_c;
   xi_U2[ijk] = z-z_c;
   }
   const double *const n_U0 = normal[0];
   const double *const n_U1 = normal[1];
   const double *const n_U2 = normal[2];
   for (ijk = 0; ijk < nn; ++ijk)
   {
   double psim2 = 
pow(psi[ijk], -2);

   double Pn_U2 = 
psim2*(AConfIJ_U0U0[ijk]*gConf_D0D0[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + 
AConfIJ_U0U0[ijk]*gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U1[ijk] + 
AConfIJ_U0U0[ijk]*pow(gConf_D0D2[ijk], 2)*n_U2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U0[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D2[ijk]*gConf_D1D1[ijk]*n_U1[ijk] + 2.0*AConfIJ_U0U1[ijk]*
gConf_D0D2[ijk]*gConf_D1D2[ijk]*n_U2[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D0[ijk]*gConf_D2D2[ijk]*n_U0[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D1[ijk]*gConf_D2D2[ijk]*n_U1[ijk] + AConfIJ_U0U2[ijk]*
pow(gConf_D0D2[ijk], 2)*n_U0[ijk] + AConfIJ_U0U2[ijk]*gConf_D0D2[ijk]*
gConf_D1D2[ijk]*n_U1[ijk] + 2.0*AConfIJ_U0U2[ijk]*gConf_D0D2[ijk]*
gConf_D2D2[ijk]*n_U2[ijk] + AConfIJ_U1U1[ijk]*gConf_D0D1[ijk]*
gConf_D1D2[ijk]*n_U0[ijk] + AConfIJ_U1U1[ijk]*gConf_D1D1[ijk]*
gConf_D1D2[ijk]*n_U1[ijk] + AConfIJ_U1U1[ijk]*pow(gConf_D1D2[ijk], 2)*
n_U2[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*gConf_D2D2[ijk]*
n_U0[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]*gConf_D1D2[ijk]*
n_U0[ijk] + AConfIJ_U1U2[ijk]*gConf_D1D1[ijk]*gConf_D2D2[ijk]*
n_U1[ijk] + AConfIJ_U1U2[ijk]*pow(gConf_D1D2[ijk], 2)*n_U1[ijk] + 2.0*
AConfIJ_U1U2[ijk]*gConf_D1D2[ijk]*gConf_D2D2[ijk]*n_U2[ijk] + 
AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]*gConf_D2D2[ijk]*n_U0[ijk] + 
AConfIJ_U2U2[ijk]*gConf_D1D2[ijk]*gConf_D2D2[ijk]*n_U1[ijk] + 
AConfIJ_U2U2[ijk]*pow(gConf_D2D2[ijk], 2)*n_U2[ijk]);

   double Pn_U1 = 
psim2*(AConfIJ_U0U0[ijk]*gConf_D0D0[ijk]*gConf_D0D1[ijk]*n_U0[ijk] + 
AConfIJ_U0U0[ijk]*pow(gConf_D0D1[ijk], 2)*n_U1[ijk] + AConfIJ_U0U0[ijk]*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D0[ijk]*gConf_D1D1[ijk]*n_U0[ijk] + AConfIJ_U0U1[ijk]*
pow(gConf_D0D1[ijk], 2)*n_U0[ijk] + 2.0*AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D1[ijk]*n_U1[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D2[ijk]*n_U2[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D2[ijk]*gConf_D1D1[ijk]*n_U2[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U0[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + 2.0*AConfIJ_U0U2[ijk]*
gConf_D0D1[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D1[ijk]*gConf_D2D2[ijk]*n_U2[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D2[ijk]*gConf_D1D2[ijk]*n_U2[ijk] + AConfIJ_U1U1[ijk]*
gConf_D0D1[ijk]*gConf_D1D1[ijk]*n_U0[ijk] + AConfIJ_U1U1[ijk]*
pow(gConf_D1D1[ijk], 2)*n_U1[ijk] + AConfIJ_U1U1[ijk]*gConf_D1D1[ijk]*
gConf_D1D2[ijk]*n_U2[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*
gConf_D1D2[ijk]*n_U0[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]*
gConf_D1D1[ijk]*n_U0[ijk] + 2.0*AConfIJ_U1U2[ijk]*gConf_D1D1[ijk]*
gConf_D1D2[ijk]*n_U1[ijk] + AConfIJ_U1U2[ijk]*gConf_D1D1[ijk]*
gConf_D2D2[ijk]*n_U2[ijk] + AConfIJ_U1U2[ijk]*pow(gConf_D1D2[ijk], 2)*
n_U2[ijk] + AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]*gConf_D1D2[ijk]*
n_U0[ijk] + AConfIJ_U2U2[ijk]*pow(gConf_D1D2[ijk], 2)*n_U1[ijk] + 
AConfIJ_U2U2[ijk]*gConf_D1D2[ijk]*gConf_D2D2[ijk]*n_U2[ijk]);

   double Pn_U0 = 
psim2*(AConfIJ_U0U0[ijk]*pow(gConf_D0D0[ijk], 2)*n_U0[ijk] + 
AConfIJ_U0U0[ijk]*gConf_D0D0[ijk]*gConf_D0D1[ijk]*n_U1[ijk] + 
AConfIJ_U0U0[ijk]*gConf_D0D0[ijk]*gConf_D0D2[ijk]*n_U2[ijk] + 2.0*
AConfIJ_U0U1[ijk]*gConf_D0D0[ijk]*gConf_D0D1[ijk]*n_U0[ijk] + 
AConfIJ_U0U1[ijk]*gConf_D0D0[ijk]*gConf_D1D1[ijk]*n_U1[ijk] + 
AConfIJ_U0U1[ijk]*gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U2[ijk] + 
AConfIJ_U0U1[ijk]*pow(gConf_D0D1[ijk], 2)*n_U1[ijk] + AConfIJ_U0U1[ijk]*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U2[ijk] + 2.0*AConfIJ_U0U2[ijk]*
gConf_D0D0[ijk]*gConf_D0D2[ijk]*n_U0[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D0[ijk]*gConf_D1D2[ijk]*n_U1[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D0[ijk]*gConf_D2D2[ijk]*n_U2[ijk] + AConfIJ_U0U2[ijk]*
gConf_D0D1[ijk]*gConf_D0D2[ijk]*n_U1[ijk] + AConfIJ_U0U2[ijk]*
pow(gConf_D0D2[ijk], 2)*n_U2[ijk] + AConfIJ_U1U1[ijk]*
pow(gConf_D0D1[ijk], 2)*n_U0[ijk] + AConfIJ_U1U1[ijk]*gConf_D0D1[ijk]*
gConf_D1D1[ijk]*n_U1[ijk] + AConfIJ_U1U1[ijk]*gConf_D0D1[ijk]*
gConf_D1D2[ijk]*n_U2[ijk] + 2.0*AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*
gConf_D0D2[ijk]*n_U0[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*
gConf_D1D2[ijk]*n_U1[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D1[ijk]*
gConf_D2D2[ijk]*n_U2[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]*
gConf_D1D1[ijk]*n_U1[ijk] + AConfIJ_U1U2[ijk]*gConf_D0D2[ijk]*
gConf_D1D2[ijk]*n_U2[ijk] + AConfIJ_U2U2[ijk]*pow(gConf_D0D2[ijk], 2)*
n_U0[ijk] + AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]*gConf_D1D2[ijk]*
n_U1[ijk] + AConfIJ_U2U2[ijk]*gConf_D0D2[ijk]*gConf_D2D2[ijk]*
n_U2[ijk]);

   double xiP_U2 = 
-Pn_U0*xi_U1[ijk] + Pn_U1*xi_U0[ijk];

   double xiP_U1 = 
Pn_U0*xi_U2[ijk] - Pn_U2*xi_U0[ijk];

   double xiP_U0 = 
-Pn_U1*xi_U2[ijk] + Pn_U2*xi_U1[ijk];


   /* populating: */
   SPIN_integrand_U0[ijk] = xiP_U0;
   SPIN_integrand_U1[ijk] = xiP_U1;
   SPIN_integrand_U2[ijk] = xiP_U2;
   }

   free(xi_U0);
   free(xi_U1);
   free(xi_U2);
}
