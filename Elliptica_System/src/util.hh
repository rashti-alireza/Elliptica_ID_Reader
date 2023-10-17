#ifndef EID_UTIL_HPP_
#define EID_UTIL_HPP_

#include <cassert>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cstdlib>
#include <stdio.h>
#include <vector>

#define MAX_NTAB   (16001)
#define IMAX(a, b) (a > b ? a : b)
#define IMIN(a, b) (a < b ? a : b)
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#ifndef DBL_EPSILON
#define DBL_EPSILON (1e-15)
#endif

using namespace std;

void eid_huntloc(CCTK_REAL xx[], CCTK_INT n, CCTK_REAL x, CCTK_INT *jlo);
CCTK_REAL eid_interploc(
  CCTK_REAL xp[],
  CCTK_REAL yp[],
  CCTK_INT np,
  CCTK_REAL xb,
  CCTK_INT *n_nearest_pt);

void eid_load_beta_equilloc(
  const char beta_equil_file[],
  CCTK_REAL log_rho0_table[MAX_NTAB],
  CCTK_REAL Y_e_table[MAX_NTAB],
  CCTK_INT *n_tab_betaloc);

void eid_set_dt_from_domega(
  CCTK_ARGUMENTS,
  CCTK_REAL const * const var,
  CCTK_REAL * const dtvar,
  CCTK_REAL const &omega);

#endif // EID_UTIL_HPP_
