/*
// Alireza Rashti
// February 2023
*/

#include "idr_main.h"

/* exporting Elliptica initial data for evolution codes */

/* field dictionary */
static const char *const Field_Dictionary[] =
{
"alpha",/* lapse: alpha */
"betax","betay","betaz",/* shift: beta^i */
"adm_gxx","adm_gxy","adm_gxz",/* symmetic ADM metric: g_ij */
"adm_gyy","adm_gyz","adm_gzz",/* symmetic ADM metric: g_ij */
"adm_Kxx","adm_Kxy","adm_Kxz",/* symmetic ADM extrinsic curvature: K_ij */
"adm_Kyy","adm_Kyz","adm_Kzz",/* symmetic ADM extrinsic curvature: K_ij */

/* matter part */
"grhd_rho",/* primitive rho (rest mass density, rho0 in Elliptica) */
"grhd_p",/* primitive p (pressure) */
"grhd_epsl",/* primitive epsilon (specific_internal_energy). 
            // note: total_energy_density = grhd_rho(1+grhd_epsl) */
"grhd_vx","grhd_vy","grhd_vz",/* primitive v, measured by an Eulerian observer, 
                              // v^i = u^i/(alpha u^0) + beta^i / alpha
                              // where u^{mu}=(u^0,u^i) is the 4-velocity of the fluid */
                         
  0/* --> detemine the last pointer */
};



/* tutorial: how to use ID reader in an evolution code
// ---------------------------------------------------

const char *checkpnt_path  = "full/path/to/elliptica/checkpoint/file"
const int Npnts = 16*16*16; // for all x,y,z coords

// initialize
Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path,"generic");

// the list of fields to be interpolated. should be comma separated
idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy";
idr->npoints  = Npnts;
// here we convert a 3D index to 1D index, e.g., (i,j,k) -> ijk.
idr->x_coords = a pointer to double type 1D(ijk) array of Cartesian x coord values;
idr->y_coords = a pointer to double type 1D(ijk) array of Cartesian y coord values;
idr->z_coords = a pointer to double type 1D(ijk) array of Cartesian z coord values;

// set parameters for elliptica. These won't change Elliptica's checkpoint file.
// example:
idr->set_param("BH_filler_method","ChebTn_Ylm_perfect_s2",idr);
idr->set_param("ADM_B1I_form","zero",idr);

// get double parameters for elliptica, for example:
idr->get_param_dbl("BHNS_angular_velocity",idr);

// interpolate
elliptica_id_reader_interpolate(idr);

// now in an evolution code one can get interpolated values as follows:
int ijk = 0;
for(i,j,k)
{
    double g_xx  = idr->field[idr->indx("adm_gxx")][ijk];
    double betax = idr->field[idr->indx("betax")][ijk];
    double alpha = idr->field[idr->indx("alpha")][ijk];
    
    ijk++;
}

// free
elliptica_id_reader_free(idr);

// note: for optimization inside the for(i,j,k) loop, 
// one can save the field indices beforehand like:

const int iell_alpha     = idr->indx("alpha");
const int iell_betax     = idr->indx("betax");
const int iell_betay     = idr->indx("betay");
const int iell_betaz     = idr->indx("betaz");

const int iell_adm_gxx   = idr->indx("adm_gxx");
const int iell_adm_gxy   = idr->indx("adm_gxy");
const int iell_adm_gxz   = idr->indx("adm_gxz");
const int iell_adm_gyy   = idr->indx("adm_gyy");
const int iell_adm_gyz   = idr->indx("adm_gyz");
const int iell_adm_gzz   = idr->indx("adm_gzz");

const int iell_adm_Kxx   = idr->indx("adm_Kxx");
const int iell_adm_Kxy   = idr->indx("adm_Kxy");
const int iell_adm_Kxz   = idr->indx("adm_Kxz");
const int iell_adm_Kyy   = idr->indx("adm_Kyy");
const int iell_adm_Kyz   = idr->indx("adm_Kyz");
const int iell_adm_Kzz   = idr->indx("adm_Kzz");

const int iell_grhd_rho  = idr->indx("grhd_rho");
const int iell_grhd_epsl = idr->indx("grhd_epsl");
const int iell_grhd_p    = idr->indx("grhd_p");
const int iell_grhd_vx   = idr->indx("grhd_vx");
const int iell_grhd_vy   = idr->indx("grhd_vy");
const int iell_grhd_vz   = idr->indx("grhd_vz");

----------------------------------------------------------------------
if you want a thread safe interpolation
----------------------------------------------------------------------

// initialize (not thread safe)
Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path,"generic_MT_safe");

// this is needed only for BHNS system
idr->set_param("BH_filler_method","ChebTn_Ylm_perfect_s2",idr);

// this is needed for both BHNS and BNS systems
idr->set_param("ADM_B1I_form","zero",idr);

// preparation and setting some interpolation settings (not thread safe)
elliptica_id_reader_interpolate(idr);

// now in an evolution code one can get interpolated values as follows:
for(i,j,k)
{
    // set x,y,z coord of the evolution code.
    double x = x_coord_of_evo_code(...);
    double y = y_coord_of_evo_code(...);
    double z = z_coord_of_evo_code(...);
    
    // populate point wise. (THREAD SAFE)
    double g_xx  = idr->fieldx(idr,"adm_gxx",x,y,z);
    double betax = idr->fieldx(idr,"betax",x,y,z);
    double alpha = idr->fieldx(idr,"alpha",x,y,z);
    
    ijk++;
}

// free (thread safe)
elliptica_id_reader_free(idr);

*/

// for adding project
int Initial_Data_Reader(void *vp)
{
  UNUSED(vp);
  return EXIT_SUCCESS;
}

/* -> create a struct for initial data reader. 
// checkpnt: full path to the checkpoint file.
// ---------
// note checkpoint file has all info about the system.
//
// option: for any special treatment/request from elliptica reader.
// -------
// note: it is case insensitive.
// options include:
// generic: a standard (default) method for exporting ID, what is done for BAM. 
// generic_MT_safe: a standard (default) method for exporting ID, what is done for BAM 
//                  and multi thread safe.
//
// --> return value: alloc and return a pointer to workspace for ID reader struct */
Elliptica_ID_Reader_T *elliptica_id_reader_init (
  const char *const checkpnt/* path/to/elliptica/checkpoint/file */,
  const char *const option/* the option for elliptica reader */
  )
{
// NOTE: in this function one should not add or use any global variables of Elliptica
// as they are not initialized yet, e.g., timer, parameters, etc. 
// ONLY set things in idr struct.

  Uint nf;
  
  // some sanity checks
  if (!checkpnt)
  {
    Error1("No checkpoint file path.");
  }
  if (!option)
  {
    Error1("No option specified.");
  }
  
  Elliptica_ID_Reader_T *idr = calloc(1,sizeof(*idr));
  IsNull(idr);
  FILE *file = 0;
  Parameter_T *par = 0;
  
  // set path
  idr->checkpoint_path = dup_s(checkpnt);
  
  // set index finder
  idr->indx = find_field_index;
  
  // settig param functions
  idr->set_param     = set_param_from_evo;
  idr->get_param_dbl = get_param_double_from_checkpoint;
  
  // read checkpoint file
  file = Fopen(idr->checkpoint_path,"r");

  // which ID system
  par = parameter_query_from_checkpoint("Project",file);
  if (!par)
  {
    Errors("I could not find the parameter name '%s'!","Project");
  }
  
  idr->system = dup_s(par->rv);
  // add a project parameter
  idr->set_param("Project",idr->system,idr);
  // which option
  idr->option = dup_s(option);
  
  // alloc field. NOTE: is allocs index for all field anyway
  for (nf = 0; Field_Dictionary[nf]; nf++);
  idr->nfield = nf;
  idr->field  = calloc(nf,sizeof(*idr->field));
  IsNull(idr->field);

  free_given_parameter(par);
  Fclose(file);

  return idr;
}

/* set parameter for reader from the evo code side.
// e.g., set_param_from_evo("bh_filler", "on",idr); --> bh_filler = on. */
static void set_param_from_evo(
          const char *const lv/* e.g., force_balance */, 
          const char *const rv/* e.g., on */,
          Elliptica_ID_Reader_T *const idr)
{
  // checks
  if(!lv && !rv)
  {
    Error1("A param with empty content.");
  }
  
  Uint nparams = idr->nparams;
  
  idr->params_lv = realloc(idr->params_lv,(nparams+1)*(sizeof(*idr->params_lv)));
  IsNull(idr->params_lv);
  idr->params_lv[nparams] = dup_s(lv);

  idr->params_rv = realloc(idr->params_rv,(nparams+1)*(sizeof(*idr->params_rv)));
  IsNull(idr->params_rv);
  idr->params_rv[nparams] = dup_s(rv);
  
  idr->nparams++;
}

/* get a double parameter from elliptica's checkpoint file. 
// note: it opens the checkpint file and takes a query.
// e.g.,  double Omega = get_param_double_from_evo("BHNS_angular_velocity",idr); */
static double get_param_double_from_checkpoint( 
          const char *const lv/* e.g., "BHNS_angular_velocity" */,
          Elliptica_ID_Reader_T *const idr)
{
  double ret = DBL_MAX;
  // check
  if(!lv) return ret;
  
  FILE *file = 0;
  Parameter_T *par = 0;

  // read checkpoint file
  file = Fopen(idr->checkpoint_path,"r");
  par  = parameter_query_from_checkpoint(lv,file);
  if (!par)
  {
    Fclose(file);
    Errors("I could not find the parameter name '%s'!",lv);
  }
  ret = strtod(par->rv,0);
  
  // free
  free_given_parameter(par);
  Fclose(file);
  
  return ret;
}

/* given the field name, it returns the index number of the field.
// gives error if it couldn't find the field. 
// NOTE: it's case sensitive
// NOTE: we are assuming the order of field are as Field_Dictionary.
// ->return: field index. */
static Uint find_field_index(const char *const fname)
{
  Uint indx = UINT_MAX;
  Uint nf   = 0;
  
  while (Field_Dictionary[nf])
  {
    if (!strcmp(fname,Field_Dictionary[nf]))
    {
      indx = nf;
      break;
    }
    nf++;
  }
  
  if (indx == UINT_MAX)
  {
    Errors("could not find the field '%s'",fname);
  }
  
  return indx;
}

/* call the pertinent exporter to interpolate fields. */
/* ->return: success */
int elliptica_id_reader_interpolate(Elliptica_ID_Reader_T *const idr)
{
  // NOTE: since we add params, using time, and etc we need this.
  init_global_variables(".");

  FUNC_TIC

  if (strcmp_i(idr->system,"BH_NS_binary_initial_data") &&
      strcmp_i(idr->option,"generic"))
  {
    SANITY_CHECK
    Psets(P_"BHNS_export_id","generic");
    BH_NS_Binary_Initial_Data(idr);
  }
  else if (strcmp_i(idr->system,"BH_NS_binary_initial_data") &&
           strcmp_i(idr->option,"generic_MT_safe"))
  {
    // no sanity check needed! if ifield is not set we set them automatically.
    Psets(P_"BHNS_export_id","generic_MT_safe");
    idr->fieldx = idr_interpolate_field_thread_safe;
    BH_NS_Binary_Initial_Data(idr);
  }
  else if (strcmp_i(idr->system,"NS_NS_binary_initial_data") &&
           strcmp_i(idr->option,"generic"))
  {
    SANITY_CHECK
    Psets(P_"NSNS_export_id","generic");
    NS_NS_Binary_Initial_Data(idr);
  }
  else if (strcmp_i(idr->system,"BH_BH_binary_initial_data") &&
           strcmp_i(idr->option,"generic"))
  {
    SANITY_CHECK
    Psets(P_"BHBH_export_id","generic");
    BH_BH_Binary_Initial_Data(idr);
  }
  else if (strcmp_i(idr->system,"Single_NS_initial_data") &&
           strcmp_i(idr->option,"generic"))
  {
    SANITY_CHECK
    Psets(P_"SNS_export_id","generic");
    Single_NS_Initial_Data(idr);
  }
  else
  {
    Error1(NO_OPTION);
  }
  
  FUNC_TOC
 
  return EXIT_SUCCESS;
}

/* ->return success. free everything */
int elliptica_id_reader_free(Elliptica_ID_Reader_T *idr)
{
  Uint n;

  /* free parameter data base */  
  free_parameter_db();
  /* free grid data base */
  free_grid_db();
  
  if (!idr)
    return EXIT_SUCCESS;
  
  Free(idr->checkpoint_path);
  Free(idr->system);
  Free(idr->option);
  free_2d(idr->id_field_names);
  
  // free fields
  for (n = 0; n < idr->nfield; n++)
  {
    Free(idr->field[n]);
  }
  Free(idr->field);
  
  // free params
  for (n = 0; n < idr->nparams; n++)
  {
    Free(idr->params_lv[n]);
    Free(idr->params_rv[n]);
    
  }
  Free(idr->params_lv);
  Free(idr->params_rv);
  
  Free(idr);
  return EXIT_SUCCESS;
}

