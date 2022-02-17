/* ===========================================================================

   FILENAME:     swp_f_vol_interpol.cxx

   PURPOSE:      Provide functions to interpolate on a volatility matrix

  =============================================================================
*/
#include "math.h"
#include "swp_h_all.h"
#include "swp_h_vol_interpol.h"

/* ========================================================================
This function interpolates on a matrix that looks like:
         Swaption Expiry as a Date
    S    0	01/01/97	01/01/98	01/01/99
    w    91    17%...
    a   182
    p   365
        547
    M   730
    a  1095
    t  1825
        where:
                - array[0][i] is the ith expiry date of the swaptions
                - array[j][0] is the jth underlying swap maturity
                - array[j][i] is a [i] into a [j] length swaption
   ======================================================================== */

Err vol_interp_from_frd(double **array, int nrow, int ncol, double exp_date,
                        double daysfromstart, int method, double today,
                        double *vol) {
  int i, j;
  double vola, volb;
  int exit = 0;

  if (exp_date <= 0 || daysfromstart <= 0)
    return serror("Vol dates too small");

  i = 1;
  while (array[0][i] < exp_date) {
    i++;
    if (i >= ncol)
      break;
  }
  j = 1;
  while (array[j][0] < daysfromstart) {
    j++;
    if (j >= nrow)
      break;
  }

  if (i == 1) {
    if (j == 1)
      *vol = array[1][1];
    else if (j >= nrow)
      *vol = array[nrow - 1][1];
    else
      *vol = array[j - 1][1] + (daysfromstart - array[j - 1][0]) /
                                   (array[j][0] - array[j - 1][0]) *
                                   (array[j][1] - array[j - 1][1]);
  } else if (i >= ncol) {
    if (j == 1)
      *vol = array[1][ncol - 1];
    else if (j >= nrow)
      *vol = array[nrow - 1][ncol - 1];
    else
      *vol = array[j - 1][ncol - 1] +
             (daysfromstart - array[j - 1][0]) /
                 (array[j][0] - array[j - 1][0]) *
                 (array[j][ncol - 1] - array[j - 1][ncol - 1]);
  } else {
    if (j == 1) {
      vola = array[1][i - 1];
      volb = array[1][i];
    } else if (j >= nrow) {
      vola = array[nrow - 1][i - 1];
      volb = array[nrow - 1][i];
    } else {
      vola = array[j - 1][i - 1] + (daysfromstart - array[j - 1][0]) /
                                       (array[j][0] - array[j - 1][0]) *
                                       (array[j][i - 1] - array[j - 1][i - 1]);
      volb = array[j - 1][i] + (daysfromstart - array[j - 1][0]) /
                                   (array[j][0] - array[j - 1][0]) *
                                   (array[j][i] - array[j - 1][i]);
    }
    if (method == 0) /* interpolate linearly */
      *vol = vola + (exp_date - array[0][i - 1]) /
                        (array[0][i] - array[0][i - 1]) * (volb - vola);
    else /* interpolate via (sigma^2 * t) */
      *vol = sqrt(
          (vola * vola * (array[0][i - 1] - today) * (array[0][i] - exp_date) +
           volb * volb * (array[0][i] - today) * (exp_date - array[0][i - 1])) /
          ((array[0][i] - array[0][i - 1]) * (exp_date - today)));
  }

  return NULL;
}

/* ========================================================================= */

/* ---------------------------------------------------------------------- */

/* ========================================================================
This function was originally in acc_addin from New York
This function interpolates on a matrix that looks like:
                Underlying Swap Maturity
        0		  1	 2	 3	 4	 5	 10
    O   01/01/97	17%...
    p   01/06/97
    t   01/01/98
        01/01/99
    E   01/01/00
    x   01/01/01
    p     ...
    i
    r
    y
        where:
         - exp_dates_vec[i] is the ith expiry date of the swaptions(date)
         - und_mats_vec[j]  is the jth underlying swap maturity (years)
         - array[i][j] is a [i] into a [j] length swaption
   ======================================================================== */

double vol_interpol_function(double tgt_exp_date, double tgt_und_mat,
                             double *exp_dates_vec, long l_exp_dates,
                             double *und_mats_vec, long l_und_mats,
                             double **mkt_vol_2darray, long m_exp_dates,
                             long m_und_mats) {
  /* note rows index "exp_dates_vec"        , columns index "und_mats_vec" */
  long i, j;
  double tgt_vol_by_date_l, tgt_vol_by_date_u, interpolated_vol;

  /*** NOW PROCEED WITH REQUISITE VOLATILITY INTERPOLATION ***/

  /***	find times on either side of tgt_und_mat & tgt_exp_date ***/

  for (j = 0; und_mats_vec[j] < tgt_und_mat && j < l_und_mats; j++)
    ;
  for (i = 0; exp_dates_vec[i] < tgt_exp_date && i < l_exp_dates; i++)
    ;

  /**	boundary conditions **/
  /* boundary corners */
  if (i == 0 && j == 0) {
    interpolated_vol = mkt_vol_2darray[0][0];
  }

  if (i == l_exp_dates && j == l_und_mats) {
    interpolated_vol = mkt_vol_2darray[l_exp_dates - 1][l_und_mats - 1];
  }

  if (i == 0 && j == l_und_mats) {
    interpolated_vol = mkt_vol_2darray[0][l_und_mats - 1];
  }

  if (i == l_exp_dates && j == 0) {
    interpolated_vol = mkt_vol_2darray[l_exp_dates - 1][0];
  }

  /* boundary top row        , i.e.interpolating on und_mats_vec only */
  if (i == 0 && j < l_und_mats && j > 0) {
    interpolated_vol = mkt_vol_2darray[0][j - 1] +
                       (mkt_vol_2darray[0][j] - mkt_vol_2darray[0][j - 1]) *
                           (tgt_und_mat - und_mats_vec[j - 1]) /
                           (und_mats_vec[j] - und_mats_vec[j - 1]);
  }

  /* boundary bottom row        , i.e.interpolating on und_mats_vec only */
  if (i == l_exp_dates && j < l_und_mats && j > 0) {
    interpolated_vol = mkt_vol_2darray[l_exp_dates - 1][j - 1] +
                       (mkt_vol_2darray[l_exp_dates - 1][j] -
                        mkt_vol_2darray[l_exp_dates - 1][j - 1]) *
                           (tgt_und_mat - und_mats_vec[j - 1]) /
                           (und_mats_vec[j] - und_mats_vec[j - 1]);
  }

  /* boundary leftmost column        , i.e. interpolating on exp_dates_vec
   * only*/
  if (j == 0 && i < l_exp_dates && i > 0) {
    interpolated_vol = mkt_vol_2darray[i - 1][0] +
                       (mkt_vol_2darray[i][0] - mkt_vol_2darray[i - 1][0]) *
                           (tgt_exp_date - exp_dates_vec[i - 1]) /
                           (exp_dates_vec[i] - exp_dates_vec[i - 1]);
  }

  /* boundary rightmost column i.e. interpolating on exp_dates_vec only*/
  if (j == l_und_mats && i < l_exp_dates && i > 0) {
    interpolated_vol = mkt_vol_2darray[i - 1][l_und_mats - 1] +
                       (mkt_vol_2darray[i][l_und_mats - 1] -
                        mkt_vol_2darray[i - 1][l_und_mats - 1]) *
                           (tgt_exp_date - exp_dates_vec[i - 1]) /
                           (exp_dates_vec[i] - exp_dates_vec[i - 1]);
  }

  /**     standard        , interior solution        , note interpolating on
  2Dim array_grid using mat_date vector as root index vector        , i.e. first
  interpolating across exp_dates_vec and then across und_mats_vec
  **/
  if (i > 0 && i < l_exp_dates && j > 0 && j < l_und_mats) {
    tgt_vol_by_date_l =
        mkt_vol_2darray[i - 1][j - 1] +
        (mkt_vol_2darray[i][j - 1] - mkt_vol_2darray[i - 1][j - 1]) *
            (tgt_exp_date - exp_dates_vec[i - 1]) /
            (exp_dates_vec[i] - exp_dates_vec[i - 1]);
    tgt_vol_by_date_u = mkt_vol_2darray[i - 1][j] +
                        (mkt_vol_2darray[i][j] - mkt_vol_2darray[i - 1][j]) *
                            (tgt_exp_date - exp_dates_vec[i - 1]) /
                            (exp_dates_vec[i] - exp_dates_vec[i - 1]);

    interpolated_vol =
        tgt_vol_by_date_l + (tgt_vol_by_date_u - tgt_vol_by_date_l) *
                                (tgt_und_mat - und_mats_vec[j - 1]) /
                                (und_mats_vec[j] - und_mats_vec[j - 1]);
  }
  return (interpolated_vol);
}
