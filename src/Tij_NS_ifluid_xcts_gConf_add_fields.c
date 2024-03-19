/*
  These C codes generated by Cpi version 3.0
  Copyright (C) 2019-2023 Alireza Rashti.
*/


#include "Tij_header.h"


void Tij_NS_idealfluid_XCTS_gConf_add_fields(Grid_T *const grid)
{
 Uint p;
 FOR_ALL_PATCHES(p,grid)
 {
 Patch_T *patch = grid->patch[p];

  ADD_FIELD(enthalpy)
  ADD_FIELD(denthalpy_D0)
  ADD_FIELD(denthalpy_D1)
  ADD_FIELD(denthalpy_D2)
  ADD_FIELD(u0)
  ADD_FIELD(du0_D0)
  ADD_FIELD(du0_D1)
  ADD_FIELD(du0_D2)
  ADD_FIELD(rho0)
  ADD_FIELD(drho0_D0)
  ADD_FIELD(drho0_D1)
  ADD_FIELD(drho0_D2)
  ADD_FIELD(e0)
  ADD_AND_ALLOC_FIELD(EConf)
  ADD_AND_ALLOC_FIELD(JConf_U0)
  ADD_AND_ALLOC_FIELD(JConf_U1)
  ADD_AND_ALLOC_FIELD(JConf_U2)
  ADD_AND_ALLOC_FIELD(SConf)
  ADD_AND_ALLOC_FIELD(EConfP)
  ADD_AND_ALLOC_FIELD(JConfP_U0)
  ADD_AND_ALLOC_FIELD(JConfP_U1)
  ADD_AND_ALLOC_FIELD(JConfP_U2)
  ADD_AND_ALLOC_FIELD(SConfP)
  ADD_AND_ALLOC_FIELD(EConfC);WRITE_v(EConfC)
  ADD_AND_ALLOC_FIELD(JConfC);WRITE_v(JConfC)
  ADD_AND_ALLOC_FIELD(SConfC);WRITE_v(SConfC)

FOR_ALL_ijk
{
EConfC[ijk] = 1.;
SConfC[ijk] = 1.;
JConfC[ijk] = 1.;
}

 }
}
