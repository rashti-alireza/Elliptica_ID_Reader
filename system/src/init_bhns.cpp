/* (c) 2023 Pedro Espino & Alireza Rashti */

#include "util.hpp"

// Auxiliary variables for interepolating Beta equilibrium table

void Elliptica_BHNS_initialize(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT n_tab_betaloc;
  CCTK_REAL Y_e_tabloc[MAX_NTAB], log_rho0_tab_betaloc[MAX_NTAB];
  CCTK_INT n_nearest_betaloc;

  if (verbose)
  {
    CCTK_INFO("Entering Elliptica_BHNS_initialize");
  }

  if (init_real)
  {
    CCTK_INFO("(with realistic EOS)");
    eid_load_beta_equilloc(
      beta_file, log_rho0_tab_betaloc, Y_e_tabloc, &n_tab_betaloc);
  }

  // Other quantities in terms of Cactus units
  CCTK_INT keyerr = 0, anyerr = 0;

  //Get EOS_Omni handle
  if (!(*init_eos_key = EOS_Omni_GetHandle(eos_table)))
    CCTK_WARN(0, "Cannot get initial eos handle, aborting...");
  CCTK_VInfo(
    CCTK_THORNSTRING, "Elliptica will use the %s equation of state.",
    eos_table);
  CCTK_VInfo(
    CCTK_THORNSTRING, "Elliptica will use the %d eos handle",
    *init_eos_key);

  //Set up local coordinate arrays
  CCTK_INFO("Setting up coordinates");
  CCTK_INT const N_points = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  CCTK_REAL *xx = (CCTK_REAL *)calloc(N_points, sizeof(*xx));
  assert(xx);
  CCTK_REAL *yy = (CCTK_REAL *)calloc(N_points, sizeof(*yy));
  assert(yy);
  CCTK_REAL *zz = (CCTK_REAL *)calloc(N_points, sizeof(*zz));
  assert(zz);

#pragma omp parallel for
  for (CCTK_INT i = 0; i < N_points; ++i)
  {
    xx[i] = x[i];
    yy[i] = y[i];
    zz[i] = z[i];
  }

  // --------------------------------------------------------------
  //   CHECKING FILE NAME EXISTENCE
  // --------------------------------------------------------------
  FILE *file;
  if ((file = fopen(Elliptica_id_file, "r")) != NULL)
    fclose(file);
  else
  {
    CCTK_VError(
      __LINE__, __FILE__, CCTK_THORNSTRING,
      "File \"%s\" does not exist. ABORTING", Elliptica_id_file);
  }
  // when using EOS, check for EOS file.
  //if (strlen(eos_table_filepath) > 0) {
  //  if (setenv("LORENE_TABULATED_EOS_PATH", eos_table_filepath, 1)) {
  //    CCTK_ERROR("Unable to set environment variable LORENE_TABULATED_EOS_PATH");

  //  }
  //}

  CCTK_VInfo(CCTK_THORNSTRING, "Reading from file \"%s\"", Elliptica_id_file);

  try
  {
    CCTK_VInfo(CCTK_THORNSTRING, "Calling Elliptica_ID_Reader_T");
    Elliptica_ID_Reader_T *idr =
      elliptica_id_reader_init(Elliptica_id_file, Elliptica_bhns_option);
    CCTK_REAL K     = poly_K;     // make sure ths is in polytropic units
    CCTK_REAL Gamma = poly_gamma; // make sure ths is in polytropic units
    idr->ifields =
      "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,grhd_rho,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
    idr->npoints  = N_points;
    idr->x_coords = xx;
    idr->y_coords = yy;
    idr->z_coords = zz;
    idr->set_param("ADM_B1I_form", "zero", idr);
    CCTK_REAL Omega = 0.0;
    if (ID_type == "BHNS")
    {
      Omega = idr->get_param_dbl("BHNS_angular_velocity", idr);
      idr->set_param("BH_filler_method", BH_filler_method, idr);
    }
    else
    {
      Omega = idr->get_param_dbl("NSNS_angular_velocity", idr);
    }

    elliptica_id_reader_interpolate(idr);

#pragma omp parallel for
    for (CCTK_INT i = 0; i < N_points; ++i)
    {

      if (CCTK_EQUALS(initial_lapse, "Elliptica"))
      {
        alp[i] = idr->field[idr->indx("alpha")][i];
      }

      //TODO: this is modified by a negative sign from LORENE ID. Is that needed here?
      if (CCTK_EQUALS(initial_shift, "Elliptica"))
      {
        betax[i] = idr->field[idr->indx("betax")][i];
        betay[i] = idr->field[idr->indx("betay")][i];
        betaz[i] = idr->field[idr->indx("betaz")][i];
      }

      if (CCTK_EQUALS(initial_data, "Elliptica"))
      {
        gxx[i] = idr->field[idr->indx("adm_gxx")][i];
        gxy[i] = idr->field[idr->indx("adm_gxy")][i];
        gxz[i] = idr->field[idr->indx("adm_gxz")][i];
        gyy[i] = idr->field[idr->indx("adm_gyy")][i];
        gyz[i] = idr->field[idr->indx("adm_gyz")][i];
        gzz[i] = idr->field[idr->indx("adm_gzz")][i];

        kxx[i] = idr->field[idr->indx("adm_Kxx")][i];
        kxy[i] = idr->field[idr->indx("adm_Kxy")][i];
        kxz[i] = idr->field[idr->indx("adm_Kxz")][i];
        kyy[i] = idr->field[idr->indx("adm_Kyy")][i];
        kyz[i] = idr->field[idr->indx("adm_Kyz")][i];
        kzz[i] = idr->field[idr->indx("adm_Kzz")][i];
      }

      if (CCTK_EQUALS(initial_data, "Elliptica"))
      {
        rho[i] = idr->field[idr->indx("grhd_rho")][i];
        //if using a realistic EOS, set the beta equilibrium conditions
        if (init_real)
        {
          if (rho[i] >= 1e-7)
          {
            CCTK_REAL yeres = eid_interploc(
              log_rho0_tab_betaloc, Y_e_tabloc, n_tab_betaloc, log10(rho[i]),
              &n_nearest_betaloc);
            if (yeres <= 0.036)
            {
              yeres = 0.036;
            }
            Y_e[i]         = yeres;
            temperature[i] = 0.1;
          }
          else
          {
            Y_e[i] = 0.25;
          }
          temperature[i] = 0.1;
        }
        if (!recalculate_eps)
        { //we don't know the temperature, so assume epsilon from ID is correct.
          eps[i] = idr->field[idr->indx("grhd_epsl")][i];
        }
        // Pressure from EOS_Omni call
        if (
          CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::temperature") > 0
          && CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::Y_e") > 0)
        {
          EOS_Omni_press(
            *init_eos_key, recalculate_eps, eos_precision, 1, &(rho[i]),
            &(eps[i]), &(temperature[i]), &(Y_e[i]), &(press[i]), &keyerr,
            &anyerr);
        }
        else
        {
          EOS_Omni_press(
            *init_eos_key, recalculate_eps, eos_precision, 1, &(rho[i]),
            &(eps[i]), NULL, NULL, &(press[i]), &keyerr, &anyerr);
        }

        vel[i]                = idr->field[idr->indx("grhd_vx")][i];
        vel[i + N_points]     = idr->field[idr->indx("grhd_vy")][i];
        vel[i + 2 * N_points] = idr->field[idr->indx("grhd_vz")][i];

        // Especially the velocity is set to strange values outside of the
        // matter region, so take care of this in the following way
        if (rho[i] < 1.e-15)
        {
          rho[i]                = 1.e-15;
          vel[i]                = 0.0;
          vel[i + N_points]     = 0.0;
          vel[i + 2 * N_points] = 0.0;
          eps[i]                = K * pow(rho[i], Gamma - 1.) / (Gamma - 1.);
          press[i]              = K * pow(rho[i], Gamma);
        }
      }

    } // for i

    {
      // These initial data assume a helical Killing vector field

      if (CCTK_EQUALS(initial_lapse, "Elliptica"))
      {
        if (CCTK_EQUALS(initial_dtlapse, "Elliptica"))
        {
          CCTK_INFO("Calculating time derivatives of lapse");
          eid_set_dt_from_domega(CCTK_PASS_CTOC, alp, dtalp, Omega);
        }
        else if (
          CCTK_EQUALS(initial_dtlapse, "none")
          or CCTK_EQUALS(initial_dtlapse, "zero"))
        {
          // do nothing
        }
        else
        {
          CCTK_WARN(CCTK_WARN_ABORT, "internal error setting dtlapse");
        }
      }

      if (CCTK_EQUALS(initial_shift, "Elliptica"))
      {
        if (CCTK_EQUALS(initial_dtshift, "Elliptica"))
        {
          CCTK_INFO("Calculating time derivatives of shift");
          eid_set_dt_from_domega(CCTK_PASS_CTOC, betax, dtbetax, Omega);
          eid_set_dt_from_domega(CCTK_PASS_CTOC, betay, dtbetay, Omega);
          eid_set_dt_from_domega(CCTK_PASS_CTOC, betaz, dtbetaz, Omega);
        }
        else if (
          CCTK_EQUALS(initial_dtshift, "none")
          or CCTK_EQUALS(initial_dtshift, "zero"))
        {
          // do nothing
        }
        else
        {
          CCTK_WARN(CCTK_WARN_ABORT, "internal error setting dtshift");
        }
      }
    }

    // free
    elliptica_id_reader_free(idr);
    CCTK_INFO("Done.");
  }
  catch (ios::failure e)
  {
    CCTK_VWarn(
      CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
      "Could not read initial data from file '%s': %s", Elliptica_id_file,
      e.what());
  }

  if (verbose)
  {
    CCTK_INFO("Exiting Elliptica_BHNS_initialize");
  }

  if (xx)
    free(xx);
  if (yy)
    free(yy);
  if (zz)
    free(zz);
}
