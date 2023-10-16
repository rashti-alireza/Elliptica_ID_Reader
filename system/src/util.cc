#include "util.hpp"

/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/***************************************************************************/
void eid_huntloc(CCTK_REAL xx[], CCTK_INT n, CCTK_REAL x, CCTK_INT *jlo)
{
  CCTK_INT jm, jhi, inc, ascnd;

  ascnd = (xx[n] > xx[1]);
  if (*jlo <= 0 || *jlo > n)
  {
    *jlo = 0;
    jhi  = n + 1;
  }
  else
  {
    inc = 1;
    if (x >= xx[*jlo] == ascnd)
    {
      if (*jlo == n)
        return;
      jhi = (*jlo) + 1;
      while (x >= xx[jhi] == ascnd)
      {
        *jlo = jhi;
        inc += inc;
        jhi = (*jlo) + inc;
        if (jhi > n)
        {
          jhi = n + 1;
          break;
        }
      }
    }
    else
    {
      if (*jlo == 1)
      {
        *jlo = 0;
        return;
      }
      jhi = (*jlo);
      *jlo -= 1;
      while (x < xx[*jlo] == ascnd)
      {
        jhi = (*jlo);
        inc += inc;
        *jlo = jhi - inc;
        if (*jlo < 1)
        {
          *jlo = 0;
          break;
        }
      }
    }
  }
  while (jhi - (*jlo) != 1)
  {
    jm = (jhi + (*jlo)) >> 1;
    if (x > xx[jm] == ascnd)
      *jlo = jm;
    else
      jhi = jm;
  }
}
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */
/*************************************************************************/
CCTK_REAL eid_interploc(
  CCTK_REAL xp[],
  CCTK_REAL yp[],
  CCTK_INT np,
  CCTK_REAL xb,
  CCTK_INT *n_nearest_pt)
{
  CCTK_INT k, /* index of 1st point */
    m = 4;    /* degree of interpolation */

  CCTK_REAL y; /* intermediate value */

  eid_huntloc(xp, np, xb, n_nearest_pt);

  k = IMIN(IMAX((*n_nearest_pt) - (m - 1) / 2, 1), np + 1 - m);

  if (xb == xp[k] || xb == xp[k + 1] || xb == xp[k + 2] || xb == xp[k + 3])
    xb += DBL_EPSILON;

  y = (xb - xp[k + 1]) * (xb - xp[k + 2]) * (xb - xp[k + 3]) * yp[k]
      / ((xp[k] - xp[k + 1]) * (xp[k] - xp[k + 2]) * (xp[k] - xp[k + 3]))

    + (xb - xp[k]) * (xb - xp[k + 2]) * (xb - xp[k + 3]) * yp[k + 1]
      / ((xp[k + 1] - xp[k]) * (xp[k + 1] - xp[k + 2]) * (xp[k + 1] - xp[k + 3]))

    + (xb - xp[k]) * (xb - xp[k + 1]) * (xb - xp[k + 3]) * yp[k + 2]
      / ((xp[k + 2] - xp[k]) * (xp[k + 2] - xp[k + 1]) * (xp[k + 2] - xp[k + 3]))

    + (xb - xp[k]) * (xb - xp[k + 1]) * (xb - xp[k + 2]) * yp[k + 3]
      / ((xp[k + 3] - xp[k]) * (xp[k + 3] - xp[k + 1]) * (xp[k + 3] - xp[k + 2]));

  return (y);
}
/*************************************************************************/
/* Load Beta equil file.                                                        */
/*************************************************************************/
void eid_load_beta_equilloc(
  const char beta_equil_file[],
  CCTK_REAL log_rho0_table[MAX_NTAB],
  CCTK_REAL Y_e_table[MAX_NTAB],
  CCTK_INT *n_tab_betaloc)
{
  CCTK_INT i; /* counter */

  CCTK_REAL rho0, /* density */
    ye;           /* electron fraction */

  FILE *f_beta; /* pointer to beta_equil_file */

  /* OPEN FILE TO READ */

  if ((f_beta = fopen(beta_equil_file, "r")) == NULL)
  {
    CCTK_VERROR("cannot open beta-equil. file:  %s\n", beta_equil_file);
  }

  /* READ NUMBER OF TABULATED POINTS */

  fscanf(f_beta, "%d\n", n_tab_betaloc);

  /* READ ENERGY DENSITY, P, H, N0 AND CONVERT TO CACTUS UNITS */

  for (i = 1; i <= (*n_tab_betaloc); i++)
  {
    fscanf(f_beta, "%lf %lf \n", &rho0, &ye);
    log_rho0_table[i] = log10(rho0); /* multiply by C^2 to get energy density */
    if (ye <= 0.036)
    {
      Y_e_table[i] = 0.036;
    }
    else
    {
      Y_e_table[i] = ye;
    }
  }
}

void eid_set_dt_from_domega(
  CCTK_ARGUMENTS,
  CCTK_REAL const * const var,
  CCTK_REAL * const dtvar,
  CCTK_REAL const &omega)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INT const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<CCTK_REAL> dxvar(npoints), dyvar(npoints);

  Diff_gv(cctkGH, 0, var, &dxvar[0], -1);
  Diff_gv(cctkGH, 1, var, &dyvar[0], -1);

#pragma omp parallel for
  for (CCTK_INT i = 0; i < npoints; ++i)
  {
    CCTK_REAL const ephix    = +y[i];
    CCTK_REAL const ephiy    = -x[i];
    CCTK_REAL const dphi_var = ephix * dxvar[i] + ephiy * dyvar[i];
    dtvar[i]                 = omega * dphi_var;
  }
}
