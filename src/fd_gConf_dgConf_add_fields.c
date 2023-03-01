/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "fd_header.h"


void fd_add_fields_gConf_igConf_dgConf(Grid_T *const grid);
void fd_add_fields_gConf_igConf_dgConf(Grid_T *const grid)
{
 Uint p;
 FOR_ALL_PATCHES(p,grid)
 {
 Patch_T *patch = grid->patch[p];


  /* declaring: */
  ADD_AND_ALLOC_FIELD(gConf_D0D2)
  ADD_AND_ALLOC_FIELD(gConf_D0D0)
  ADD_AND_ALLOC_FIELD(gConf_D0D1)
  ADD_AND_ALLOC_FIELD(gConf_D1D2)
  ADD_AND_ALLOC_FIELD(gConf_D1D1)
  ADD_AND_ALLOC_FIELD(gConf_D2D2)
  ADD_AND_ALLOC_FIELD(igConf_U2U2)
  ADD_AND_ALLOC_FIELD(igConf_U1U2)
  ADD_AND_ALLOC_FIELD(igConf_U1U1)
  ADD_AND_ALLOC_FIELD(igConf_U0U2)
  ADD_AND_ALLOC_FIELD(igConf_U0U0)
  ADD_AND_ALLOC_FIELD(igConf_U0U1)
  ADD_AND_ALLOC_FIELD(dgConf_D0D2D1)
  ADD_AND_ALLOC_FIELD(dgConf_D0D2D2)
  ADD_AND_ALLOC_FIELD(dgConf_D1D1D1)
  ADD_AND_ALLOC_FIELD(dgConf_D1D2D1)
  ADD_AND_ALLOC_FIELD(dgConf_D0D1D0)
  ADD_AND_ALLOC_FIELD(dgConf_D0D1D1)
  ADD_AND_ALLOC_FIELD(dgConf_D0D1D2)
  ADD_AND_ALLOC_FIELD(dgConf_D0D2D0)
  ADD_AND_ALLOC_FIELD(dgConf_D0D0D1)
  ADD_AND_ALLOC_FIELD(dgConf_D0D0D0)
  ADD_AND_ALLOC_FIELD(dgConf_D1D2D0)
  ADD_AND_ALLOC_FIELD(dgConf_D0D0D2)
  ADD_AND_ALLOC_FIELD(dgConf_D1D2D2)
  ADD_AND_ALLOC_FIELD(dgConf_D1D1D0)
  ADD_AND_ALLOC_FIELD(dgConf_D1D1D2)
  ADD_AND_ALLOC_FIELD(dgConf_D2D2D1)
  ADD_AND_ALLOC_FIELD(dgConf_D2D2D0)
  ADD_AND_ALLOC_FIELD(dgConf_D2D2D2)


 }
}
