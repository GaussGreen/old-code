/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optscdout(14)
 *
 * PURPOSE      	: SCUD Option
 *
 * DESCRIPTION  	: Barrier on a secondary underlying
 *
 * CALLS		: bivar()
 *		: srt_f_optscdout()
 *
 * PARAMETERS   	: fwdx     	- forward price of 1st underlying x
 *              	: fwdy          - forward price of 2nd underlying y
 *              	: spoty         - spot price of 2nd underlying y
 *              	: strike        - strike price
 *              	: barrier       - barrier level
 *              	: sigx          - ?? of 1st underlying x
 *              	: sigy          - ?? of 2nd underlying y
 *              	: rho           - ??
 *              	: mat           - maturity of option        , in years
 *              	: disc          - discount factor
 *              	: call_put      - type of option: 0 call        , 1 put
 *              	: down_up     - type of scud: ??
 *              	: greek        	- ??
 *
 * RETURNS      	: ??          	- ??
 *
 *******************************************************************************/

double srt_f_optscdout(double fwdx, double fwdy, double spoty, double strike,
                       double barrier, double sigx, double sigy, double rho,
                       double mat1, double mat2, double disc,
                       SrtCallPutType call_put, SrtBarrierType down_up,
                       SrtGreekType greek) {
  double mux;
  double muy;

  double l;
  double a;
  double alpha;
  double beta;

  double d1x;
  double d2x;
  double d1y;
  double d2y;

  double premium;
  double result;

  double shift;
  double shiftx;
  double shifty;

  int cp;
  int du;

  if ((down_up == SRT_DOWN) && (spoty < barrier)) {
    premium = 0;
  } else if ((down_up == SRT_UP) && (spoty > barrier)) {
    premium = 0;
  } else if (mat1 <= 0) {
    premium =
        srt_f_optblksch(fwdx, strike, sigx, mat2, disc, call_put, PREMIUM);
  } else {
    mux = -(sigx * sigx / 2);
    muy = log(fwdy / spoty) / mat1 - (sigy * sigy / 2);

    l = log(barrier / spoty);
    a = exp(2 * muy * l / (sigy * sigy));
    alpha = 2 * rho * l * sigx / sigy;
    beta = (2 * l);

    d1x = (log(fwdx / strike) + mux * mat2 + sigx * sigx * mat2) /
          (sigx * sqrt(mat2));
    d2x = d1x - sigx * sqrt(mat2);

    d1y = (log(spoty / barrier) + muy * mat1 + rho * sigy * sigx * mat1) /
          (sigy * sqrt(mat1));
    d2y = d1y - rho * sigx * sqrt(mat1);

    cp = (call_put == SRT_CALL) ? 1 : -1;
    du = (down_up == SRT_DOWN) ? 1 : -1;

    premium = fwdx * bivar(cp * d1x, du * d1y, cp * du * rho);
    premium -= fwdx * a * exp(alpha) *
               bivar(cp * (d1x + (alpha / (sigx * sqrt(mat1)))),
                     du * (d1y + (beta / (sigy * sqrt(mat1)))), cp * du * rho);
    premium -= strike * bivar(cp * d2x, du * d2y, cp * du * rho);
    premium += strike * a *
               bivar(cp * (d2x + (alpha / (sigx * sqrt(mat1)))),
                     du * (d2y + (beta / (sigy * sqrt(mat1)))), cp * du * rho);
    premium *= disc * cp;
  }

  switch (greek) {
  case PREMIUM:
    return (premium);
    break;

  case DELTAX:
    shift = fwdx / 10000;
    result =
        (srt_f_optscdout(fwdx + shift, fwdy, spoty, strike, barrier, sigx, sigy,
                         rho, mat1, mat2, disc, call_put, down_up, PREMIUM) -
         premium) /
        shift;
    return (result);
    break;

  case DELTAY:
    shift = fwdy / 10000;
    result = (srt_f_optscdout(fwdx, fwdy + shift, spoty + shift, strike,
                              barrier, sigx, sigy, rho, mat1, mat2, disc,
                              call_put, down_up, PREMIUM) -
              premium) /
             shift;
    return (result);
    break;

  case GAMMAX:
    shift = fwdx / 1000;
    result =
        srt_f_optscdout(fwdx + shift, fwdy, spoty, strike, barrier, sigx, sigy,
                        rho, mat1, mat2, disc, call_put, down_up, PREMIUM);
    result +=
        srt_f_optscdout(fwdx - shift, fwdy, spoty, strike, barrier, sigx, sigy,
                        rho, mat1, mat2, disc, call_put, down_up, PREMIUM);
    result -= 2 * premium;
    result /= shift * shift;
    return (result);
    break;

  case GAMMAY:
    shift = fwdy / 1000;
    result = srt_f_optscdout(fwdx, fwdy + shift, spoty + shift, strike, barrier,
                             sigx, sigy, rho, mat1, mat2, disc, call_put,
                             down_up, PREMIUM);
    result += srt_f_optscdout(fwdx, fwdy - shift, spoty - shift, strike,
                              barrier, sigx, sigy, rho, mat1, mat2, disc,
                              call_put, down_up, PREMIUM);
    result -= 2 * premium;
    result /= shift * shift;
    return (result);
    break;

  case GAMMAXY:
    shifty = fwdy / 1000;
    shiftx = fwdx / 1000;
    result = srt_f_optscdout(
        fwdx + shiftx, fwdy + shifty, spoty * (1 + shifty / fwdy), strike,
        barrier, sigx, sigy, rho, mat1, mat2, disc, call_put, down_up, PREMIUM);
    result += srt_f_optscdout(
        fwdx - shiftx, fwdy - shifty, spoty * (1 - shifty / fwdy), strike,
        barrier, sigx, sigy, rho, mat1, mat2, disc, call_put, down_up, PREMIUM);
    result -= srt_f_optscdout(
        fwdx - shiftx, fwdy + shifty, spoty * (1 + shifty / fwdy), strike,
        barrier, sigx, sigy, rho, mat1, mat2, disc, call_put, down_up, PREMIUM);
    result -= srt_f_optscdout(
        fwdx + shiftx, fwdy - shifty, spoty * (1 - shifty / fwdy), strike,
        barrier, sigx, sigy, rho, mat1, mat2, disc, call_put, down_up, PREMIUM);
    result /= 4 * shiftx * shifty;
    return (result);
    break;

  case VEGAX:
    shift = GVOPT.vol_add;
    result =
        (srt_f_optscdout(fwdx, fwdy, spoty, strike, barrier, sigx + shift, sigy,
                         rho, mat1, mat2, disc, call_put, down_up, PREMIUM) -
         premium) /
        shift;
    return (result);
    break;

  case VEGAY:
    shift = GVOPT.vol_add;
    result =
        (srt_f_optscdout(fwdx, fwdy, spoty, strike, barrier, sigx, sigy + shift,
                         rho, mat1, mat2, disc, call_put, down_up, PREMIUM) -
         premium) /
        shift;
    return (result);
    break;

  case THETA:
    shift = YEARS_IN_DAY;
    result = srt_f_optscdout(fwdx, fwdy, spoty, strike, barrier, sigx, sigy,
                             rho, mat1 - shift, mat2 - shift,
                             disc * exp(-shift * log(disc) / mat2), call_put,
                             down_up, PREMIUM) -
             premium;
    return (result);
    break;

  default:
    return (UNKNOWN_GREEK);
    break;
  }

} /* END srt_f_optscdout() */

/******************************************************************************/
