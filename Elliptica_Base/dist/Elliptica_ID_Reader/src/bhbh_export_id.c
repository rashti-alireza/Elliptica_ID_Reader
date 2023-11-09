/*
// Alireza Rashti
// August 2023
*/

/* exporting initial data for evolution codes */

#include "bhbh_export_id.h"


/* export ID for a general evolution code */
void bhbh_export_id_generic(void *vp)
{
  FUNC_TIC
  
  // sanity check
  assert(vp && "The input is null!");  
  
  Elliptica_ID_Reader_T *const idr = vp;
  Physics_T *bhbh = 0;
  ID_Reader_T *points = idr_init();
  double CM[3] = {0.};
  char fields_name[STR_LEN_MAX] = {'\0'};// elliptica field names
  char **sfield = 0;
  Uint f;
  
  /* adapt fields_notations for Elliptica and ensure given fields are expected.
  // NOTE: the following replaces the name with the same ORDER as the input. */
  assert(sprintf(fields_name,"%s",idr->ifields));
  
  /* metric fields */
  regex_replace(fields_name,"\\balpha\\b",EVO_"alpha",fields_name);
  
  regex_replace(fields_name,"\\bbetax\\b",EVO_"beta_U0",fields_name);
  regex_replace(fields_name,"\\bbetay\\b",EVO_"beta_U1",fields_name);
  regex_replace(fields_name,"\\bbetaz\\b",EVO_"beta_U2",fields_name);
  
  regex_replace(fields_name,"\\badm_gxx\\b",EVO_"adm_g_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_gxy\\b",EVO_"adm_g_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_gxz\\b",EVO_"adm_g_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_gyy\\b",EVO_"adm_g_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_gyz\\b",EVO_"adm_g_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_gzz\\b",EVO_"adm_g_D2D2",fields_name);
  
  regex_replace(fields_name,"\\badm_Kxx\\b",EVO_"adm_Kij_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_Kxy\\b",EVO_"adm_Kij_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_Kxz\\b",EVO_"adm_Kij_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_Kyy\\b",EVO_"adm_Kij_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_Kyz\\b",EVO_"adm_Kij_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_Kzz\\b",EVO_"adm_Kij_D2D2",fields_name);
  
  /* check if all fields are expected */
  sfield = read_separated_items_in_string(fields_name,',');
  f = 0;
  while(sfield[f])
  {
    if (!regex_search("^"EVO_,sfield[f]))
    {
      printf("field '%s' is Unexpected!\n",sfield[f]);
      Error1("Unexpected field! Add me!");
    }
    f++;
  }
  free_2d(sfield);
  
  /* set all parameters set in evo code. */
  // NOTE: adding both plain and modify version to make sure things kick in 
  // after loading from checkpoint file.
  for (Uint n = 0; n < idr->nparams; ++n)
  {
    char lv[STR_LEN_MAX];
    sprintf(lv,CHECKPOINT_SET_PARAM_"%s",idr->params_lv[n]);
    
    Psets(idr->params_lv[n],idr->params_rv[n]);
    Psets(lv,idr->params_rv[n]);
  }
  
  /* don't stop for checkpoint */
  Pseti(CHECKPOINT_SET_PARAM_ P_"STOP",0);
  /* read physics from checkpoint */
  Psets("checkpoint_file_path",idr->checkpoint_path);
  bhbh = bhbh_read_physics_from_checkpoint();
  points->grid = bhbh->grid;
  
  /* get (x,y,z) points from evo. NOTE: no allocation done for (x,y,z) */
  CM[0] = Pgetd(P_"x_CM");
  CM[1] = Pgetd(P_"y_CM");
  CM[2] = Pgetd(P_"z_CM");
  idr_find_XYZ_from_xyz(idr,points,CM);
  
  /* go from Omega x r to inertial coords sys asymptotically.
  // EVO needs this!? (this param should be set in the reader) */
  // Psets(CHECKPOINT_SET_PARAM_ "ADM_B1I_form","zero");
  physics(bhbh,ADM_UPDATE_Kij);/* need this for bhbh_set_evo_fields_generic */

  /* fill BH1 */
  Physics_T *const bh1  = init_physics(bhbh,BH1);
  Pseti("BH1_filler_verbose",1);/* make it verbose anyway. */
  /* fill these fields */
  Psets("BH1_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                           "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                           "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                           "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                           "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  physics(bh1,BH_FILL);
  free_physics(bh1);
 
  /* fill BH2 */
  Physics_T *const bh2  = init_physics(bhbh,BH2);
  Pseti("BH2_filler_verbose",1);/* make it verbose anyway. */
  /* fill these fields */
  Psets("BH2_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                           "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                           "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                           "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                           "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  physics(bh2,BH_FILL);
  free_physics(bh2);
  
  /* set bam fields based on initial data to be usable for evo */
  bhbh_set_evo_fields_generic(bhbh->grid);
  
  /* write into array */
  idr_interpolate_fields_and_save_in_array(idr,points,fields_name,idr->ifields);
  
  /* finishing up */
  // since no alocation done for (x,y,z):
  points->x = 0;
  points->y = 0;
  points->z = 0;
  idr_free(points);
  free_physics(bhbh);
  
  FUNC_TOC  
}

