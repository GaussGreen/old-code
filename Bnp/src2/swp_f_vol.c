/* -------------------------------------------------------------------------

   FILE NAME	: swp_f_vol.c

   PURPOSE		: compute volatilities via a pointer to a function
   VolFunc passed to the SrtExternalFunctions static ( see srt_f_external_fct.c
   )
   ------------------------------------------------------------------------- */

#include "swp_h_all.h"
#include "swp_h_external_fct.h"
#include "swp_h_vol.h"

/* ---------------------------------------------------------------------------
        Calculate volatility   , using the VolFunc that has been attached to the
        library (see srt_f_external_fct.c)
  --------------------------------------------------------------------------- */

Err swp_f_vol(char *vol_id, Ddate start, Ddate end, double strike,
              double *volatility, double *power) {
  Err err;
  VolFuncType vol_func;

  /* Gets volatility function attached to the library */
  err = swp_f_GetVolFunc(&vol_func);
  if (err)
    return err;

  /* Computes the volatility for the end date*/
  err = vol_func(vol_id, start, end, strike, volatility, power);
  if (err)
    return err;

  return NULL;

} /* END srt_f_vol */

/* ---------------------------------------------------------------------------
        Calculate volatility   , using the SABRVolFunc that has been attached to
  the library (see srt_f_external_fct.c)
  ---------------------------------------------------------------------------
        case 0:"ATMLOGNORMAL";
        case 1:"ATMNORMAL";
        case 2:"ATMSIGMABETA";
        case 3:"LOGNORMAL";
        case 4:"NORMAL";
        case 5:"ALPHA";
        case 6:"BETA";
        case 7:"RHO";
        case 8: 1 if SABR or SABRAf or any SmileModel
        case 9:"ZETA";

*/

Err swp_f_SABRvol(char *vol_id, Ddate start, Ddate end, double strike,
                  double *volatility, double *power, int component) {
  Err err;
  SABRVolFuncType sabrvol_func;

  /* Gets volatility function attached to the library */
  err = swp_f_GetSABRVolFunc(&sabrvol_func);
  if (err)
    return err;

  /* Computes the volatility for the end date*/
  /*	if (component == 9)
          {
                  err = sabrvol_func(vol_id  , start  , end  , strike  , &dtemp
     , power  , (int)8); if (err) return err;

                  if (dtemp == 1)
                          return "not a SABR vol ";
                  else
                  {
                          *volatility = 0.6;
                          return NULL;
                  }
          }
          else
          {*/
  err = sabrvol_func(vol_id, start, end, strike, volatility, power, component);
  if (err)
    return err;
  //	}

  return NULL;

} /* END srt_f_SABRvol */

/* ---------------------------------------------------------------------------
        Calculate non standard volatility   , using Pat Hagan's trick
  --------------------------------------------------------------------------- */
Err swp_f_transformvol(char *yc_id, char *vol_id, char *freq_std,
                       char *basis_std, char *refRateCode, Ddate start,
                       Ddate end, double strike, char *freq, char *basis,
                       double *volatility, double *power) {
  Err err;
  VolFuncType vol_func;
  double R_std, R, Ravg, K;
  double a, b;
  SrtCompounding ifreq, ifreq_std;
  SrtBasisCode ibasis, ibasis_std;

  /* Gets volatility function attached to the library */
  err = swp_f_GetVolFunc(&vol_func);
  if (err)
    return err;

  err = interp_basis(basis_std, &ibasis_std);
  err = interp_basis(basis, &ibasis);
  err = interp_compounding(freq_std, &ifreq_std);
  err = interp_compounding(freq, &ifreq);

  a = coverage((long)(start), (long)(end), ibasis) /
      coverage((long)(start), (long)(end), ibasis_std);
  b = (ifreq - ifreq_std) / (2.0 * ifreq * ifreq_std);

  err = swp_f_ForwardRate((long)(start), (long)(end), freq_std, basis_std,
                          yc_id, refRateCode, &R_std);
  err = swp_f_ForwardRate((long)(start), (long)(end), freq, basis, yc_id,
                          refRateCode, &R);

  K = R_std + a * (strike - R) + a * b * (strike * strike - R * R);
  Ravg = 0.5 * (R + strike);

  /* Computes the standard volatility at the modified strike*/
  err = vol_func(vol_id, start, end, K, volatility, power);
  if (err)
    return err;

  /* Computes the non standard volatility*/
  if (*power == 1) {
    *volatility = *volatility / (1 + b * Ravg);
  } else if (*power == 0) {
    *volatility = *volatility / (a + 2 * a * b * Ravg);
  } else {
    err = "Beta vol type not supported in transform vol";
    return err;
  }

  return NULL;

} /* END srt_f_transformvol */

/* ---------------------------------------------------------------------------
        Calculate non standard volatility   , using Pat Hagan's trick
  --------------------------------------------------------------------------- */
Err swp_f_transformvol_func(
    char *yc_id, char *vol_id,
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    char *freq_std, char *basis_std, char *refRateCode, Ddate start, Ddate end,
    double strike, char *freq, char *basis, double *volatility, double *power) {
  Err err;
  double R_std, R, Ravg, K;
  double a, b;
  SrtCompounding ifreq, ifreq_std;
  SrtBasisCode ibasis, ibasis_std;

  /* Gets volatility function attached to the library */

  err = interp_basis(basis_std, &ibasis_std);
  err = interp_basis(basis, &ibasis);
  err = interp_compounding(freq_std, &ifreq_std);
  err = interp_compounding(freq, &ifreq);

  a = coverage((long)(start), (long)(end), ibasis) /
      coverage((long)(start), (long)(end), ibasis_std);
  b = (ifreq - ifreq_std) / (2.0 * ifreq * ifreq_std);

  err = swp_f_ForwardRate((long)(start), (long)(end), freq_std, basis_std,
                          yc_id, refRateCode, &R_std);
  err = swp_f_ForwardRate((long)(start), (long)(end), freq, basis, yc_id,
                          refRateCode, &R);

  K = R_std + a * (strike - R) + a * b * (strike * strike - R * R);
  Ravg = 0.5 * (R + strike);

  /* Computes the standard volatility at the modified strike*/
  err = get_cash_vol(vol_id, start, end, K, 0, refRateCode, volatility, power);
  if (err)
    return err;

  /* Computes the non standard volatility*/
  if (*power == 1) {
    *volatility = *volatility / (1 + b * Ravg);
  } else if (*power == 0) {
    *volatility = *volatility / (a + 2 * a * b * Ravg);
  } else {
    err = "Beta vol type not supported in transform vol";
    return err;
  }

  return NULL;

} /* END srt_f_transformvol */

Err swp_f_IsSmileVol(
    char *vol_id,
    double *yesno) // return 1 if Model is SABR or SABRAFor any smile model.
{
  double dontcare;

  return swp_f_SABRvol(vol_id, 0, 0, 0, yesno, &dontcare, 8);
}

/* ======================================================================= */
