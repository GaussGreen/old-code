/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optmovbar(...)
 *
 * PURPOSE      	: pricing a moving barrier
 *
 * DESCRIPTION  	: calculate the price for an option with the spot
 *staying above (or below) bar1 from T=0 to mat1 then above (or below) bar2 from
 *mat1 to mat2.
 *
 * CALLS			: SELF
 *
 * PARAMETERS   	: spot     	- spot price
 *              	: forward1  - forward price of underlying at mat1
 *              	: forward2  - forward price of underlying at mat2
 *              	: barrier1  - barrier level from T=0 to mat1
 *              	: barrier2  - barrier level from T=mat1 to mat2
 *              	: vol1      - volatility of underlying from T=0 to mat1
 *              	: vol2      - volatility of underlying from T=0 to mat2
 *              	: mat1		- maturity of 1st extinguishable period
 *              	: mat2        - maturity of 2nd extinguishable period
 *				: disc_fact	- discount factor
 *              	: greek     - ??
 *
 * RETURNS      	: premium   - ??
 *
 *******************************************************************************/

/* Integrand Function */
static double prob_inf_bar2(double x, va_list argptr) {
  double spot, mat1, mat2, bar1, bar2, mu01, mu12, sigma01, sigma12;
  double answer;

  spot = va_arg(argptr, double);
  mat1 = va_arg(argptr, double);
  mat2 = va_arg(argptr, double);
  bar1 = va_arg(argptr, double);
  bar2 = va_arg(argptr, double);
  mu01 = va_arg(argptr, double);
  mu12 = va_arg(argptr, double);
  sigma01 = va_arg(argptr, double);
  sigma12 = va_arg(argptr, double);

  answer =
      pow(exp(1), (-mu01 + pow(sigma01, 2)) * mat1) *
      (bar2 * norm((2 * mu12 * mat1 - pow(sigma12, 2) * mat1 - 2 * mu12 * mat2 +
                    pow(sigma12, 2) * mat2 + 2 * log(bar2) - 2 * log(x)) /
                   (2 * sigma12 * sqrt(-mat1 + mat2))) -
       pow(bar2 / x, 2 * mu12 / pow(sigma12, 2)) * x *
           norm((2 * mu12 * mat1 - pow(sigma12, 2) * mat1 - 2 * mu12 * mat2 +
                 pow(sigma12, 2) * mat2 - 2 * log(bar2) + 2 * log(x)) /
                (2 * sigma12 * sqrt(-mat1 + mat2)))) *
      (-(pow(bar1, 2 * mu01 / pow(sigma01, 2)) * pow(spot, 3) *
         gauss((2 * mu01 * mat1 - 3 * pow(sigma01, 2) * mat1 + 4 * log(bar1) -
                2 * log(spot) - 2 * log(x)) /
               (2 * sigma01 * sqrt(mat1)))) +
       pow(bar1, 3) * pow(spot, 2 * mu01 / pow(sigma01, 2)) *
           gauss((2 * mu01 * mat1 - 3 * pow(sigma01, 2) * mat1 + 2 * log(spot) -
                  2 * log(x)) /
                 (2 * sigma01 * sqrt(mat1)))) /
      (pow(bar1, 3) * bar2 * sigma01 *
       pow(spot, (2 * mu01 + pow(sigma01, 2)) / pow(sigma01, 2)) * sqrt(mat1));

  return answer;
}

static double prob_sup_bar2(double x, va_list argptr) {
  double spot, mat1, mat2, bar1, bar2, mu01, mu12, sigma01, sigma12;
  double answer;

  spot = va_arg(argptr, double);
  mat1 = va_arg(argptr, double);
  mat2 = va_arg(argptr, double);
  bar1 = va_arg(argptr, double);
  bar2 = va_arg(argptr, double);
  mu01 = va_arg(argptr, double);
  mu12 = va_arg(argptr, double);
  sigma01 = va_arg(argptr, double);
  sigma12 = va_arg(argptr, double);

  answer =
      pow(exp(1), (-mu01 + pow(sigma01, 2)) * mat1) *
      (pow(bar2 * x, 2 * mu12 / pow(sigma12, 2)) / x *
           norm((-2 * mu12 * mat1 + pow(sigma12, 2) * mat1 + 2 * mu12 * mat2 -
                 pow(sigma12, 2) * mat2 + 2 * log(bar2) + 2 * log(x)) /
                (2 * sigma12 * sqrt(-mat1 + mat2))) -
       bar2 *
           norm((-2 * mu12 * mat1 + pow(sigma12, 2) * mat1 + 2 * mu12 * mat2 -
                 pow(sigma12, 2) * mat2 - 2 * log(bar2) - 2 * log(x)) /
                (2 * sigma12 * sqrt(-mat1 + mat2)))) *
      (-(pow(bar1, 2 * mu01 / pow(sigma01, 2)) * pow(spot, 3) *
         gauss((2 * mu01 * mat1 - 3 * pow(sigma01, 2) * mat1 + 4 * log(bar1) -
                2 * log(spot) + 2 * log(x)) /
               (2 * sigma01 * sqrt(mat1)))) +
       pow(bar1, 3) * pow(spot, 2 * mu01 / pow(sigma01, 2)) *
           gauss((2 * mu01 * mat1 - 3 * pow(sigma01, 2) * mat1 + 2 * log(spot) +
                  2 * log(x)) /
                 (2 * sigma01 * sqrt(mat1)))) /
      (pow(bar1, 3) * bar2 * sigma01 *
       pow(spot, (2 * mu01 + pow(sigma01, 2)) / pow(sigma01, 2)) * sqrt(mat1));

  answer *= (-1 / x / x);

  return answer;
}

/* PRICING function */
double srt_f_optmovbar(double spot, double fwd1, double fwd2, double bar1,
                       double bar2, double sigma1, double sigma2, double mat1,
                       double mat2, double disc_fact,
                       SrtBarrierType below_above, SrtGreekType greek) {
  double shift, premium, answer;
  double mu01, mu12, sigma01, sigma12;
  SrtMinmaxType min_max;

  if (mat2 == 0) {
    /* we don't care about mat1 */
    if (((below_above == SRT_DOWN) && (spot < bar2) && (spot < bar1)) ||
        ((below_above == SRT_UP) && (spot > bar2) && (spot > bar1))) {
      premium = 1;
    } else
      premium = 0;
  } else if (mat1 == 0) {
    if (below_above == SRT_UP) {
      min_max = SRT_MAX;
    } else {
      min_max = SRT_MIN;
    }

    /* we have an american digital */
    premium = srt_f_optamedig(fwd2, spot, bar2, sigma2, mat2, disc_fact,
                              below_above, min_max, greek);
  } else {
    mu01 = log(fwd1 / spot) / mat1;
    mu12 = log(fwd2 / fwd1) / (mat2 - mat1);
    sigma01 = sigma1;
    sigma12 =
        sqrt((sigma2 * sigma2 * mat2 - sigma1 * sigma1 * mat1) / (mat2 - mat1));

    if (below_above == SRT_DOWN) {
      premium = sm_qsimp_list(&prob_inf_bar2, EPS, bar1 - EPS, spot, mat1, mat2,
                              bar1, bar2, mu01, mu12, sigma01, sigma12);
    } else if (below_above == SRT_UP) {
      premium = sm_qsimp_list(&prob_sup_bar2, EPS, 1 / (bar1 + EPS), spot, mat1,
                              mat2, bar1, bar2, mu01, mu12, sigma01, sigma12);
    }

    premium *= disc_fact;
  }

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTA_FWD1: /*** DELTA FWD ***/
    shift = fwd1 / 10000;
    answer =
        (srt_f_optmovbar(spot, fwd1 + shift, fwd2, bar1, bar2, sigma1, sigma2,
                         mat1, mat2, disc_fact, below_above, PREMIUM) -
         premium) /
        shift;
    break;
  case DELTA_FWD2: /*** DELTA FWD ***/
    shift = fwd2 / 10000;
    answer =
        (srt_f_optmovbar(spot, fwd1, fwd2 + shift, bar1, bar2, sigma1, sigma2,
                         mat1, mat2, disc_fact, below_above, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA: /*** DELTA FWD1 + FWD2 ***/
    shift = fwd1 / 10000;
    answer =
        (srt_f_optmovbar(spot, fwd1 * (1 + shift / fwd1),
                         fwd2 * (1 + shift / fwd1), bar1, bar2, sigma1, sigma2,
                         mat1, mat2, disc_fact, below_above, PREMIUM) -
         premium) /
        shift;
    break;
  case GAMMA: /*** GAMMA ***/
    shift = fwd1 / 10000;
    answer = srt_f_optmovbar(
        spot, fwd1 * (1 + shift / fwd1), fwd2 * (1 + shift / fwd1), bar1, bar2,
        sigma1, sigma2, mat1, mat2, disc_fact, below_above, PREMIUM);
    answer += srt_f_optmovbar(
        spot, fwd1 * (1 - shift / fwd1), fwd2 * (1 - shift / fwd1), bar1, bar2,
        sigma1, sigma2, mat1, mat2, disc_fact, below_above, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;
  case VEGA1: /*** VEGA AT mat1***/
    shift = GVOPT.vol_add;
    answer =
        (srt_f_optmovbar(spot, fwd1, fwd2, bar1, bar2, sigma1 + shift, sigma2,
                         mat1, mat2, disc_fact, below_above, PREMIUM) -
         premium) /
        shift;
    break;
  case VEGA2: /*** VEGA AT mat2***/
    shift = GVOPT.vol_add;
    answer =
        (srt_f_optmovbar(spot, fwd1, fwd2, bar1, bar2, sigma1, sigma2 + shift,
                         mat1, mat2, disc_fact, below_above, PREMIUM) -
         premium) /
        shift;
    break;
  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer = srt_f_optmovbar(spot, fwd1, fwd2, bar1, bar2, sigma1, sigma2,
                             mat1 - shift, mat2 - shift,
                             disc_fact * exp(shift * log(disc_fact) / mat2),
                             below_above, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);

} /* srt_f_optmovbar() */
/******************************************************************************/
