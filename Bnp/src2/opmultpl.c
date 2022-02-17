/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/******************************************************************************
        computes the following probabilities:
                Prob ( x(T) > k_x ; y(T) > k_y ): 	cp_x = 1; cp_y = 1
                Prob ( x(T) > k_x ; y(T) < k_y ): 	cp_x = 1; cp_y = -1
                Prob ( x(T) < k_x ; y(T) > k_y ): 	cp_x = -1; cp_y = 1
                Prob ( x(T) < k_x ; y(T) < k_y ): 	cp_x = -1; cp_y = -1
******************************************************************************/

static double proba_joint_x_y(double mu_x, double mu_y, double sig_x,
                              double sig_y, double mat, double k_x, double k_y,
                              double corr,
                              double cp_x, /*1 if call:X>K  -1 if put:X<K */
                              double cp_y  /*1 if call:Y>K  -1 if put:Y<K */
) {
  double prob;
  double sigx_sqrt = sig_x * sqrt(mat);
  double sigy_sqrt = sig_y * sqrt(mat);

  prob = bivar(cp_x * (-k_x + mu_x) / sigx_sqrt,
               cp_y * (-k_y + mu_y) / sigy_sqrt, cp_x * cp_y * corr);

  return (prob);
}

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optmultpl(...)
 *
 * PURPOSE      	: Calculates the premium for a Multiple Option
 *
 * DESCRIPTION  	: At maturity  , the holder receives the price difference
 *between some asset S1 and a strike K1  , times the price difference between a
 *2nd asset S2 and another strike K2  , times a gearing.
 *
 * CALLS		: bivar() in gen_math.c
 *		  srt_f_optblksch()
 *
 * PARAMETERS   	: fwd1     	- forward price of 1st underlying
 *              	: strike1   	- strike price for 1st underlying
 *              	: fwd2        	- forward price of 2nd underlying
 *              	: strike2       - forward price of 2nd underlying
 *              	: vol1          - annual volatility of 1st underlying
 *              	: vol2          - annual volatility of 2nd underlying
 *              	: corr          - correlation between 1st and 2nd
 *underlyings : mat           - option maturity  , in years : disc          -
 *discount factor : gearing       - gearing : call_put1    	- option on 1st
 *underlying (0: call  , 1: put) : call_put2    	- option on 2nd
 *underlying (0: call  , 1: put)
 *
 * RETURNS      	: premium       - option premium
 *
 *******************************************************************************/

double srt_f_optmultpl(double fwd1, double strike1, double fwd2, double strike2,
                       double vol1, double vol2, double corr, double mat,
                       double disc, SrtCallPutType call_put1,
                       SrtCallPutType call_put2, SrtGreekType greek) {
  double cp_x;
  double cp_y;
  double mu_x, mu_x_x;
  double mu_y, mu_y_y;
  double k_x;
  double k_y;

  double n = 0.0;
  double n_x = 0.0;
  double n_y = 0.0;
  double n_x_y = 0.0;
  double cov;

  double premium = 0.0;
  double shift, shift1, shift2;
  double answer;

  cp_x = (call_put1 == SRT_CALL) ? 1 : -1;
  cp_y = (call_put2 == SRT_CALL) ? 1 : -1;

  if (disc == 0) {
    premium = 0;
  } else if ((vol1 == 0.0) || (vol2 == 0.0) || (mat == 0)) {
    premium =
        srt_f_optblksch(fwd2, strike2, vol2, mat, disc, call_put2, PREMIUM) *
        srt_f_optblksch(fwd1, strike1, vol1, mat, disc, call_put1, PREMIUM) /
        disc;
  } else {
    mu_x = -0.5 * vol1 * vol1 * mat;
    mu_y = -0.5 * vol2 * vol2 * mat;
    k_x = log(strike1 / fwd1);
    k_y = log(strike2 / fwd2);

    mu_x_x = mu_x + vol1 * vol1 * mat;
    mu_y_y = mu_y + vol2 * vol2 * mat;

    cov = corr * vol1 * vol2 * mat;

    n_x_y = proba_joint_x_y(mu_x_x + cov, mu_y_y + cov, vol1, vol2, mat, k_x,
                            k_y, corr, cp_x, cp_y);

    n_x = proba_joint_x_y(mu_x_x, mu_y + cov, vol1, vol2, mat, k_x, k_y, corr,
                          cp_x, cp_y);

    n_y = proba_joint_x_y(mu_x + cov, mu_y_y, vol1, vol2, mat, k_x, k_y, corr,
                          cp_x, cp_y);

    n = proba_joint_x_y(mu_x, mu_y, vol1, vol2, mat, k_x, k_y, corr, cp_x,
                        cp_y);

    premium = disc * cp_x * cp_y *
              (exp(cov) * fwd1 * fwd2 * n_x_y - fwd1 * strike2 * n_x -
               fwd2 * strike1 * n_y + strike1 * strike2 * n);
  }

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTAX: /*** DELTA FWD 1***/
    shift = fwd1 / 10000;
    answer = (srt_f_optmultpl(fwd1 + shift, strike1, fwd2, strike2, vol1, vol2,
                              corr, mat, disc, call_put1, call_put2, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTAY: /*** DELTA FWD 2***/
    shift = fwd2 / 10000;
    answer = (srt_f_optmultpl(fwd1, strike1, fwd2 + shift, strike2, vol1, vol2,
                              corr, mat, disc, call_put1, call_put2, PREMIUM) -
              premium) /
             shift;
    break;

  case GAMMAX: /*** GAMMA FWD 1***/
    shift = fwd1 / 10000;
    answer = srt_f_optmultpl(fwd1 + shift, strike1, fwd2, strike2, vol1, vol2,
                             corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer += srt_f_optmultpl(fwd1 - shift, strike1, fwd2, strike2, vol1, vol2,
                              corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case GAMMAY: /*** GAMMA FWD 2***/
    shift = fwd2 / 10000;
    answer = srt_f_optmultpl(fwd1, strike1, fwd2 + shift, strike2, vol1, vol2,
                             corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer += srt_f_optmultpl(fwd1, strike1, fwd2 - shift, strike2, vol1, vol2,
                              corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case GAMMAXY: /*** CROSS GAMMA ***/
    shift1 = fwd1 / 10000;
    shift2 = fwd2 / 10000;
    answer =
        srt_f_optmultpl(fwd1 + shift1, strike1, fwd2 + shift2, strike2, vol1,
                        vol2, corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer +=
        srt_f_optmultpl(fwd1 - shift1, strike1, fwd2 - shift2, strike2, vol1,
                        vol2, corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer -=
        srt_f_optmultpl(fwd1 - shift1, strike1, fwd2 + shift2, strike2, vol1,
                        vol2, corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer -=
        srt_f_optmultpl(fwd1 + shift1, strike1, fwd2 - shift2, strike2, vol1,
                        vol2, corr, mat, disc, call_put1, call_put2, PREMIUM);
    answer /= 4 * shift1 * shift2;
    break;

  case VEGAX: /*** VEGA 1***/
    shift = GVOPT.vol_add;
    answer = (srt_f_optmultpl(fwd1, strike1, fwd2, strike2, vol1 + shift, vol2,
                              corr, mat, disc, call_put1, call_put2, PREMIUM) -
              premium) /
             shift;
    break;

  case VEGAY: /*** VEGA 2***/
    shift = GVOPT.vol_add;
    answer = (srt_f_optmultpl(fwd1, strike1, fwd2, strike2, vol1, vol2 + shift,
                              corr, mat, disc, call_put1, call_put2, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer = srt_f_optmultpl(fwd1, strike1, fwd2, strike2, vol1, vol2, corr,
                             mat - shift, disc * exp(-shift * log(disc) / mat),
                             call_put1, call_put2, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
  }

  return (answer);

} /* END srt_f_optmultpl(...) */

/******************************************************************************/
