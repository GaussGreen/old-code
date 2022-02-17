/*******************************************************************************
**
**	srt_f_optprtdbl.c
**
*******************************************************************************/

/* ==========================================================================
   include files
   ========================================================================== */

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/* ==========================================================================
   part_bar_fn(...)
        - Gaussian Method
   Here  , mu and sig are local variables (not cumulative ones)
   ========================================================================== */

double part_dbar_fn(double xx, double up, double dwn, double mu1, double mu2,
                    double sig1, double sig2, double mat1, double mat2,
                    int nb_term) {
  double result, coeff1, coeff2, cor;
  double x1_n, x2_n, bdash;
  double sig1_2, sig2_2, sig1_rt, sig2_rt, x, y;
  double n1, n2, n3, n4;
  int i;

  sig1_2 = sig1 * sig1;
  sig2_2 = sig2 * sig2;

  sig1_rt = sig1 * sqrt(mat1);
  sig2_rt = sig2 * sqrt(mat2 - mat1);

  bdash = sqrt(sig1_2 * mat1 + sig2_2 * (mat2 - mat1));
  cor = sig1 * sqrt(mat1) / bdash;
  /* Check for Numerical Roundings that could create too big a correlation */
  if (cor > 1.0)
    cor = 1.0;
  else if (cor < -1.0)
    cor = -1.0;

  result = 0.0;

  for (i = -nb_term; i <= nb_term; i++) {
    x1_n = 2 * i * (up - dwn);
    x2_n = 2 * up - x1_n;

    coeff1 = exp(mu1 * x1_n / sig1_2);

    x = (-dwn + x1_n + mu1 * mat1) / sig1_rt;
    y = (x1_n - xx + mu1 * mat1 + mu2 * (mat2 - mat1)) / bdash;
    n1 = bivar(x, y, cor);

    x = (-up + x1_n + mu1 * mat1) / sig1_rt;
    n2 = bivar(x, y, cor);

    coeff2 = exp(mu1 * x2_n / sig1_2);

    x = (-dwn + x2_n + mu1 * mat1) / sig1_rt;
    y = (x2_n - xx + mu1 * mat1 + mu2 * (mat2 - mat1)) / bdash;
    n3 = bivar(x, y, cor);

    x = (-up + x2_n + mu1 * mat1) / sig1_rt;
    n4 = bivar(x, y, cor);

    result += coeff1 * (n1 - n2) - coeff2 * (n3 - n4);
  }

  return (result);
}

/* ==========================================================================
   srt_f_optprtdbl(...)
        - Gaussian Method for a partial double barrier
   ========================================================================== */

double srt_f_optprtdbl(double fwd1, double fwd2, double spot, double strike,
                       double b_do, double b_up, double vol1, double vol2,
                       double mat1, double mat2, double disc,
                       SrtCallPutType call_put, int nb_term,
                       SrtGreekType greek) {
  double up, dwn, k, mu1, mu2, bound_dwn, bound_up, premium, shift, answer;

  double drift1, drift2, price_k, price_s, sig1, sig2, sig1_2, sig2_2;

  int cp;

  /* If the maturity 1 is less than 0 we have an european */
  if (mat1 <= 0) {
    premium =
        srt_f_optblksch(fwd2, strike, vol2, mat2, disc, call_put, PREMIUM);
  } else if ((spot <= b_do) || (spot >= b_up)) {
    premium = 0.0;
  } else {
    cp = (call_put == SRT_CALL) ? 1 : -1;

    up = log(b_up / spot);
    dwn = log(b_do / spot);

    k = log(strike / spot);

    drift1 = log(fwd1 / spot) / mat1;
    sig1 = vol1;

    if (mat1 < mat2) {
      drift2 = log(fwd2 / fwd1) / (mat2 - mat1);
      sig2 = sqrt((vol2 * vol2 * mat2 - vol1 * vol1 * mat1) / (mat2 - mat1));
    } else {
      drift2 = drift1;
      sig2 = sig1;
    }

    sig1_2 = sig1 * sig1;
    sig2_2 = sig2 * sig2;

    mu1 = drift1 - sig1_2 / 2;
    mu2 = drift2 - sig2_2 / 2;

    if (call_put == SRT_CALL) {
      bound_up = up; /*  call  */
      bound_dwn = dwn;
    } else if (call_put == SRT_PUT) {
      bound_up = (-dwn); /*  put  */
      bound_dwn = (-up);
      k = (-k);
      mu1 = (-mu1);
      mu2 = (-mu2);
      sig1_2 = (-sig1_2);
      sig2_2 = (-sig2_2);
    }

    price_s = part_dbar_fn(k, bound_up, bound_dwn, mu1 + sig1_2, mu2 + sig2_2,
                           sig1, sig2, mat1, mat2, nb_term);

    price_k = part_dbar_fn(k, bound_up, bound_dwn, mu1, mu2, sig1, sig2, mat1,
                           mat2, nb_term);

    premium = fwd2 * price_s - strike * price_k;

    premium *= cp * disc;
  }

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTA_FWD1: /*** DELTA FWD at T1 ***/
    shift = fwd1 / 10000;
    answer =
        (srt_f_optprtdbl(fwd1 + shift, fwd2, spot, strike, b_do, b_up, vol1,
                         vol2, mat1, mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA_FWD2: /*** DELTA FWD at T2 ***/
    shift = fwd2 / 10000;
    answer =
        (srt_f_optprtdbl(fwd1, fwd2 + shift, spot, strike, b_do, b_up, vol1,
                         vol2, mat1, mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA: /*** DELTA SPOT + FWD1 + FWD2 ***/
    shift = spot / 10000;
    answer =
        (srt_f_optprtdbl(fwd1 * (1 + shift / spot), fwd2 * (1 + shift / spot),
                         spot + shift, strike, b_do, b_up, vol1, vol2, mat1,
                         mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case GAMMA: /*** GAMMA ***/
    shift = spot / 10000;
    answer =
        srt_f_optprtdbl(fwd1 * (1 + shift / spot), fwd2 * (1 + shift / spot),
                        spot + shift, strike, b_do, b_up, vol1, vol2, mat1,
                        mat2, disc, call_put, nb_term, PREMIUM);
    answer =
        srt_f_optprtdbl(fwd1 * (1 - shift / spot), fwd2 * (1 - shift / spot),
                        spot - shift, strike, b_do, b_up, vol1, vol2, mat1,
                        mat2, disc, call_put, nb_term, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA1: /*** VEGA at T1***/
    shift = GVOPT.vol_add;
    answer =
        (srt_f_optprtdbl(fwd1, fwd2, spot, strike, b_do, b_up, vol1 + shift,
                         vol2, mat1, mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case VEGA2: /*** VEGA at T2***/
    shift = GVOPT.vol_add;
    answer = (srt_f_optprtdbl(fwd1, fwd2, spot, strike, b_do, b_up, vol1,
                              vol2 + shift, mat1, mat2, disc, call_put, nb_term,
                              PREMIUM) -
              premium) /
             shift;
    break;

  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer = srt_f_optprtdbl(fwd1, fwd2, spot, strike, b_do, b_up, vol1, vol2,
                             mat1 - shift, mat2 - shift,
                             disc * exp(-shift * log(disc) / mat2), call_put,
                             nb_term, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);
}
