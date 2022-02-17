

/******************************************************************************/
/*******************************************************************************
*
* FUNCTION     	: srt_f_optblkschbetaquick(...)
*                 srt_f_optblkschbetastochquick(...)
*
* PURPOSE      	: Premium of a European call option for a single risk factor
*				  using dF=a*F^{beta}dW
*                       da = v a dZ
*                       <dZ        ,dW> = rho dt
* DESCRIPTION  	: Transforms the Beta vol in the equivalent BS vol        , then
compute the required greeks accordingly
*
* CALLS		:     srt_f_optbetavoltoblkvol
                          srt_f_optblksch
*
* PARAMETERS
*	INPUT	    : fwd_price	    - forward underlying price
*              	: strike      	- strike price
*              	: vol         	- annual volatility
*              	: mat         	- initial time        , in years
*				: beta          - dF = a*F^{beta}dW
*              	: call_put      - type of option: 0 call        , 1 put
*				: greek		- info wanted (premium       ,
greeks...)
*
*
* RETURNS      	: premium       - option premium
*
*******************************************************************************/

/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/* ------------------------------------------------------------------------------
 */
double srt_f_optblkschbetaquick(double forward, double strike, double betavol,
                                double mat, double beta, double disc,
                                SrtCallPutType call_put, SrtGreekType greek) {
  Err err = NULL;
  double bsvol;
  double ans;
  double shift;
  double premium;
  /* Simple checks */
  if (mat < 0.)
    return 0.0;

  /* Transforms the BetaVol into its equivalent BsVol */
  err = srt_f_optbetavoltoblkvol(forward, strike, betavol, mat, beta, &bsvol);

  /* The premium is computed using the lognormal BS formula */
  premium =
      srt_f_optblksch(forward, strike, bsvol, mat, disc, call_put, PREMIUM);
  switch (greek) {
  case PREMIUM:
    ans = premium;
    break;

  case DELTA:
    shift = GVOPT.delta_shift;
    /*	shift = forward / 10000; */
    ans = (srt_f_optblkschbetaquick(forward + shift, strike, betavol, mat, beta,
                                    disc, call_put, PREMIUM) -
           premium) /
          shift;
    break;

  case GAMMA:
    /*shift = forward / 10000; */
    shift = GVOPT.delta_shift;

    ans = srt_f_optblkschbetaquick(forward + shift, strike, betavol, mat, beta,
                                   disc, call_put, PREMIUM);
    ans += srt_f_optblkschbetaquick(forward - shift, strike, betavol, mat, beta,
                                    disc, call_put, PREMIUM);
    ans -= 2 * premium;
    ans /= shift * shift;
    break;

  case VEGA:
    /* this does a multiplicative shift and multiplies the result by 10 (as
     * aksed for by traders) */
    shift = GVOPT.sabrvol_shift;
    /*  shift = betavol / 100; */
    ans = (srt_f_optblkschbetaquick(forward, strike, betavol * shift, mat, beta,
                                    disc, call_put, PREMIUM) -
           premium) *
          10;
    break;

  case THETA:
    /*	shift = YEARS_IN_DAY; */
    shift = GVOPT.theta_shift;
    ans = (srt_f_optblkschbetaquick(forward, strike, betavol, mat - shift, beta,
                                    disc * exp(-shift * log(disc) / mat),
                                    call_put, PREMIUM) -
           premium) /
          shift;
    break;

  default:
    ans = UNKNOWN_GREEK;
  }

  return ans;
}

/* ------------------------------------------------------------------------------
 */
double srt_f_optblkschbetastochquick(double forward, double strike,
                                     double maturity, double betavol,
                                     double alpha, double beta, double rho,
                                     double disc, SrtDiffusionType typeinput,
                                     SrtDiffusionType log_norm,
                                     SrtCallPutType call_put,
                                     SrtGreekType greek) {
  Err err = NULL;
  double volequ;
  double ans;
  double shift;
  double premium;

  /* Simple checks */
  if (maturity < 0.)
    return 0.0;

  /* Transforms the Vol Input into its equivalent BsVol */
  err = srt_f_optsarbvol(forward, strike, maturity, betavol, alpha, beta, rho,
                         typeinput, log_norm, &volequ);
  if (err)
    return -1;

  /* The premium is computed using the BS formula */

  /* Lognormal diffusion */
  if (log_norm == SRT_LOGNORMAL) {
    premium = srt_f_optblksch(forward, strike, volequ, maturity, disc, call_put,
                              PREMIUM);
  }
  /* Normal diffusion */
  else {
    premium = srt_f_optblknrm(forward, strike, volequ, maturity, disc, call_put,
                              PREMIUM);
  }

  switch (greek) {
  case PREMIUM:
    ans = premium;
    break;

  case DELTA:
    shift = GVOPT.delta_shift;
    /*		shift = forward / 10000; */
    ans = (srt_f_optblkschbetastochquick(
               forward + shift, strike, maturity, betavol, alpha, beta, rho,
               disc, typeinput, log_norm, call_put, PREMIUM) -
           premium) /
          shift;
    break;

  case GAMMA:
    shift = GVOPT.delta_shift;
    /*		shift = forward / 10000; */
    ans = srt_f_optblkschbetastochquick(forward + shift, strike, maturity,
                                        betavol, alpha, beta, rho, disc,
                                        typeinput, log_norm, call_put, PREMIUM);
    ans += srt_f_optblkschbetastochquick(
        forward - shift, strike, maturity, betavol, alpha, beta, rho, disc,
        typeinput, log_norm, call_put, PREMIUM);
    ans -= 2 * premium;
    ans /= shift * shift;
    break;

  case VEGA:
    shift = GVOPT.vol_add;
    /*		shift = betavol / 100; */
    ans = (srt_f_optblkschbetastochquick(
               forward, strike, maturity, betavol + shift, alpha, beta, rho,
               disc, typeinput, log_norm, call_put, PREMIUM) -
           premium) /
          shift;
    break;

  case THETA:
    shift = GVOPT.theta_shift;
    /*		shift = YEARS_IN_DAY; */
    ans = (srt_f_optblkschbetastochquick(
               forward, strike, maturity - shift, betavol, alpha, beta, rho,
               disc * exp(-shift * log(disc) / maturity), typeinput, log_norm,
               call_put, PREMIUM) -
           premium) /
          shift;
    break;

  case RHO:
    shift = GVOPT.rho_shift;
    /*		shift = rho / 100; */
    if (rho + shift > 1.0) {
      shift = -GVOPT.rho_shift;
    }
    ans = (srt_f_optblkschbetastochquick(
               forward, strike, maturity, betavol, alpha, beta, rho + shift,
               disc, typeinput, log_norm, call_put, PREMIUM) -
           premium) /
          shift;
    break;

  case ALPHA:
    shift = GVOPT.alpha_shift;
    /*		shift = alpha / 100; */
    ans = (srt_f_optblkschbetastochquick(
               forward, strike, maturity, betavol, alpha + shift, beta, rho,
               disc, typeinput, log_norm, call_put, PREMIUM) -
           premium) /
          shift;
    break;

  case BETA:
    shift = GVOPT.beta_shift;
    /*		shift = beta / 100; */
    if (beta + shift > 1.0) {
      shift = -GVOPT.beta_shift;
      /*		shift = -beta / 100;  */
    }
    ans = (srt_f_optblkschbetastochquick(
               forward, strike, maturity, betavol, alpha, beta + shift, rho,
               disc, typeinput, log_norm, call_put, PREMIUM) -
           premium) /
          shift;
    break;

  default:
    ans = UNKNOWN_GREEK;
  }

  return ans;
}
