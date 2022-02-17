/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optsladig(...)
 *
 * PURPOSE      	: European pricing for a slalom digital
 *
 * DESCRIPTION  	: Pays $1 at T2 if ( X(T1)<K  &  X(T2) > K )
 *
 * CALLS		: gauss(...) in gen_math.cxx
 *		  norm(...)  in gen_math.cxx
 *		  srt_f_intnrm_i1        ,i2        ,j1        ,j2  in
 *gen_math.cxx
 *
 * PARAMETERS
 *	INPUT	: fwd1		- forward underlying price
 *              	: fwd2	        -
 *		: barrier      	- strike price
 *              	: vol1         	- annual volatility
 *              	: vol2         	- annual volatility
 *		: mat1 		- date of slalom gate
 *		: mat2 		- maturity
 *              	: disc        	- discount factor to expiry
 *
 * RETURNS      	: premium       - option premium
 *
 *******************************************************************************/

double srt_f_opteursla(double fwd1, double fwd2, double barrier, double vol1,
                       double vol2, double mat1, double mat2, double disc,
                       SrtBarrierType below_above, SrtGreekType greek) {
  double mu;
  double mu1;
  double mu2;
  double s1;
  double s2;
  double a0;
  double b0;
  double a1;
  double b1;
  double integral_part3;
  double prob_max_gt_k;
  double bar_prime;
  double premium, answer, shift;

  if ((mat1 <= 0) || (vol1 == 0)) {
    if (below_above == SRT_DOWN) {
      if (fwd1 >= barrier)
        premium = srt_f_opteurdig(fwd2, barrier, vol2, mat2, disc, below_above,
                                  PREMIUM);
      else
        premium = 0.0;
    } else if (below_above == SRT_UP) {
      if (fwd1 <= barrier)
        premium = srt_f_opteurdig(fwd2, barrier, vol2, mat2, disc, below_above,
                                  PREMIUM);
      else
        premium = 0.0;
    }
  } else if ((mat2 - mat1) <= 0.0) {
    premium = 0.0;
  } else {
    s1 = vol1;
    s2 = sqrt((vol2 * vol2 * mat2 - vol1 * vol1 * mat1) / (mat2 - mat1));

    mu = log(fwd2 / fwd1) / (mat2 - mat1);
    mu1 = -s1 * s1 / 2;
    mu2 = mu - s2 * s2 / 2;
    bar_prime = log(barrier / fwd1);

    /*** calculate third part of horrendous integral ***/

    a0 = mu2 * (mat2 - mat1);
    b0 = s2 * sqrt(mat2 - mat1);
    a1 = bar_prime - mu1 * mat1;
    b1 = s1 * sqrt(mat1);

    integral_part3 = srt_f_intnrm_i2(a0, b0, a1, b1);

    /*** sum up separate parts of horrendous integral ***/

    prob_max_gt_k = integral_part3;

    premium = disc * prob_max_gt_k;
  }

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTA_FWD1: /*** DELTA FWD ***/
    shift = fwd1 / 10000;
    answer = (srt_f_opteursla(fwd1 + shift, fwd2, barrier, vol1, vol2, mat1,
                              mat2, disc, below_above, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA_FWD2: /*** DELTA FWD 2 ***/
    shift = fwd2 / 10000;
    answer = (srt_f_opteursla(fwd1, fwd2 + shift, barrier, vol1, vol2, mat1,
                              mat2, disc, below_above, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA: /*** DELTA SPOT + FWD ***/
    shift = fwd1 / 10000;
    answer = (srt_f_opteursla(fwd1 * (1 + shift / fwd1),
                              fwd2 * (1 + shift / fwd1), barrier, vol1, vol2,
                              mat1, mat2, disc, below_above, PREMIUM) -
              premium) /
             shift;
    break;

  case GAMMA: /*** GAMMA ***/
    shift = fwd1 / 10000;
    answer = srt_f_opteursla(fwd1 * (1 + shift / fwd1),
                             fwd2 * (1 + shift / fwd1), barrier, vol1, vol2,
                             mat1, mat2, disc, below_above, PREMIUM);
    answer += srt_f_opteursla(fwd1 * (1 - shift / fwd1),
                              fwd2 * (1 - shift / fwd1), barrier, vol1, vol2,
                              mat1, mat2, disc, below_above, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA1: /*** VEGA ***/
    shift = GVOPT.vol_add;
    answer = (srt_f_opteursla(fwd1, fwd2, barrier, vol1 + shift, vol2, mat1,
                              mat2, disc, below_above, PREMIUM) -
              premium) /
             shift;
    break;

  case VEGA2: /*** VEGA ***/
    shift = GVOPT.vol_add;
    answer = (srt_f_opteursla(fwd1, fwd2, barrier, vol1, vol2 + shift, mat1,
                              mat2, disc, below_above, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer = srt_f_opteursla(
                 fwd1, fwd2, barrier, vol1, vol2, mat1 - shift, mat2 - shift,
                 disc * exp(-shift * log(disc) / mat2), below_above, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);

} /* END srt_f_opteursla() */

/******************************************************************************/
