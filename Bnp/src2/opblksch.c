/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optblksch(...)
 *
 * PURPOSE      	: European option pricing for a single risk factor
 *
 * DESCRIPTION  	: Black-Scholes: calculates the option premium and
 *		  theoretical delta  , gamma  , vega and theta for an
 *		  underlying which is lognormal
 *
 * CALLS		: gauss(...) in gen_math.c
 *		  norm(...)  in gen_math.c
 *
 * PARAMETERS
 *	INPUT	: fwd_price	- forward underlying price
 *              	: strike      	- strike price
 *              	: vol         	- annual volatility
 *              	: mat         	- maturity  , in years
 *              	: disc        	- discount factor to expiry
 *              	: call_put      - type of option: 0 call  , 1 put
 *		: greek		- info wanted (premium  , greeks...)
 *
 * RETURNS      	: premium       - option premium
 *
 *******************************************************************************/

double srt_f_optblksch(double fwd_price, double strike, double vol, double mat,
                       double disc, SrtCallPutType call_put,
                       SrtGreekType greek) {
  double stdev;
  double d1;
  double d2;
  double answer;
  double nd1;
  double nd2;
  double cp;

  /* Check on negative maturity */
  if (mat < 0)
    return 0;

  if (call_put == SRT_STRADDLE) {
    answer =
        srt_f_optblksch(fwd_price, strike, vol, mat, disc, SRT_CALL, greek) +
        srt_f_optblksch(fwd_price, strike, vol, mat, disc, SRT_PUT, greek);
    return answer;
  }

  /* Set the sign of cp depending on call or put */
  cp = (call_put == SRT_CALL) ? 1 : -1;

  /* Set-up d1 and d2 (see Hull  , p. 224) to avoid divide-by-zeros */
  /* Added fwd_price >0 */
  if ((vol != 0.0) && (mat != 0.0) && (strike > 0.0) && (fwd_price > 0.0)) {
    stdev = vol * sqrt(mat);
    d1 = (log(fwd_price / strike) + (vol * vol / 2) * mat) / stdev;
    d2 = d1 - stdev;
    nd1 = norm(cp * d1);
    nd2 = norm(cp * d2);
  } else {
    /* Allow intrinsic value of underlying to be returned if required     */

    if (cp * (fwd_price - strike) > 0.0) {
      d1 = INFINITY; /* Equivalent of +oo */
      d2 = INFINITY; /* Equivalent of +oo */
      nd1 = 1;
      nd2 = 1;
    } else {
      if (fwd_price == strike) {
        d1 = 0;
        d2 = 0;
        nd1 = 0.5;
        nd2 = 0.5;
      } else {
        d1 = -INFINITY; /* Equivalent of -oo */
        d2 = -INFINITY; /* Equivalent of -oo */
        nd1 = 0;
        nd2 = 0;
      }
    }
  }

  /* ==========================================================================
   */
  /* Evaluate option-type dependent greeks (spot-deltas  , theta) and premium */
  /* ==========================================================================
   */

  if ((mat == 0.0) && (greek != PREMIUM))
    answer = 0.0;
  else {
    switch (greek) {
    case DELTA:     /*** DELTA UNDERLYING ***/
    case DELTA_FWD: /*** DELTA_FWD UNDERLYING ***/
      answer = cp * nd1 * disc;
      break;

    case GAMMA:     /*** GAMMA UNDERLYING ***/
    case GAMMA_FWD: /*** GAMMA_FWD UNDERLYING ***/
      answer = disc / (fwd_price * vol * sqrt(mat)) * gauss(d1);
      break;

    case VEGA: /*** VEGA ***/
      answer = 0.01 * disc * fwd_price * sqrt(mat) * gauss(d1);
      break;

    case THETA: /*** THETA ***/
      answer = -(1 / 365.00) *
               (0.5 * fwd_price * disc * vol / sqrt(mat) * gauss(d1) -
                cp * log(disc) * disc * strike * nd2 / mat);
      break;

    case PREMIUM: /*** PREMIUM ***/
      answer = disc * cp * (fwd_price * nd1 - strike * nd2);
      break;

    case RHO: /*** RHO ***/
      answer = disc * cp * strike * mat * nd2;
      break;

    case RHORHO:
      answer = -disc * cp * strike * mat * mat * nd2 +
               cp * disc * mat * strike * gauss(d2) * (sqrt(mat) / vol);
      break;

    case RHODELTA:
      answer =
          -cp * mat * disc * nd1 + cp * disc * gauss(d1) * (sqrt(mat) / vol);
      break;

    case VANNA:
      answer = -disc * d2 * gauss(d1) / vol;
      break;

    case VOLGA:
      answer = disc * fwd_price * sqrt(mat) * d2 * d1 * gauss(d1) / vol;
      break;

    default:
      answer = UNKNOWN_GREEK;
      break;
    }
  }

  return (answer);

} /* END srt_f_optblksch() */

/* ---------------------------------------------------------------------------------
 */

/* A function that checks the Black-Scholes inputs */
Err srt_f_optchkblksch(double fwd, double strike, double vol, double mat) {
  if (mat < 0)
    return serror("Negative Maturity in BlackScholes inputs");

  if (fwd < 0)
    return serror("Negative Forward in BlackScholes inputs");

  if (strike < 0)
    return serror("Negative Strike in BlackScholes inputs");

  if (vol < 0)
    return serror("Negative Volatility in BlackScholes inputs");

  return NULL;
}

/* ---------------------------------------------------------------------------
 */

/* A Higher level interface function to call Black Scholes */
Err OptBlkSch(double fwd_price, double strike, double vol, double mat,
              double disc, char *call_put_str, char *greek_str,
              double *result) {
  Err err = NULL;
  SrtCallPutType call_put_type;
  SrtGreekType greek_type = PREMIUM;

  /* Transforms the call_put string into a type */
  err = interp_call_put(call_put_str, &call_put_type);
  if (err)
    return err;

  /* Transforms the greek string into a type */
  if (strcmp(greek_str, "") != 0) {
    err = interp_greeks(greek_str, &greek_type);
    if (err)
      return err;
  }

  /* Check the inputs */
  err = srt_f_optchkblksch(fwd_price, strike, vol, mat);
  if (err)
    return err;

  /* Compute the requested greek */
  *result = srt_f_optblksch(fwd_price, strike, vol, mat, disc, call_put_type,
                            greek_type);

  /* Return a success message */
  return NULL;
}

double srt_f_optblksch_accurate(double fwd_price, double strike, double vol,
                                double mat, double disc,
                                SrtCallPutType call_put, SrtGreekType greek) {
  double stdev;
  double d1;
  double d2;
  double answer;
  double nd1;
  double nd2;
  double cp;

  /* Check on negative maturity */
  if (mat < 0)
    return 0;

  if (call_put == SRT_STRADDLE) {
    answer = srt_f_optblksch_accurate(fwd_price, strike, vol, mat, disc,
                                      SRT_CALL, greek) +
             srt_f_optblksch(fwd_price, strike, vol, mat, disc, SRT_PUT, greek);
    return answer;
  }

  /* Set the sign of cp depending on call or put */
  cp = (call_put == SRT_CALL) ? 1 : -1;

  /* Set-up d1 and d2 (see Hull  , p. 224) to avoid divide-by-zeros */
  if ((vol != 0.0) && (mat != 0.0) && (strike > 0)) {
    stdev = vol * sqrt(mat);
    d1 = (log(fwd_price / strike) + (vol * vol / 2) * mat) / stdev;
    d2 = d1 - stdev;
    nd1 = norm_accurate(cp * d1);
    nd2 = norm_accurate(cp * d2);
  } else {
    /* Allow intrinsic value of underlying to be returned if required     */

    if (cp * (fwd_price - strike) > 0.0) {
      d1 = INFINITY; /* Equivalent of +oo */
      d2 = INFINITY; /* Equivalent of +oo */
      nd1 = 1;
      nd2 = 1;
    } else {
      if (fwd_price == strike) {
        d1 = 0;
        d2 = 0;
        nd1 = 0.5;
        nd2 = 0.5;
      } else {
        d1 = -INFINITY; /* Equivalent of -oo */
        d2 = -INFINITY; /* Equivalent of -oo */
        nd1 = 0;
        nd2 = 0;
      }
    }
  }

  /* ==========================================================================
   */
  /* Evaluate option-type dependent greeks (spot-deltas  , theta) and premium */
  /* ==========================================================================
   */

  if ((mat == 0.0) && (greek != PREMIUM))
    answer = 0.0;
  else {
    switch (greek) {
    case DELTA:     /*** DELTA UNDERLYING ***/
    case DELTA_FWD: /*** DELTA_FWD UNDERLYING ***/
      answer = cp * nd1 * disc;
      break;

    case GAMMA:     /*** GAMMA UNDERLYING ***/
    case GAMMA_FWD: /*** GAMMA_FWD UNDERLYING ***/
      answer = disc / (fwd_price * vol * sqrt(mat)) * gauss(d1);
      break;

    case VEGA: /*** VEGA ***/
      answer = 0.01 * disc * fwd_price * sqrt(mat) * gauss(d1);
      break;

    case THETA: /*** THETA ***/
      answer = -(1 / 365.00) *
               (0.5 * fwd_price * disc * vol / sqrt(mat) * gauss(d1) -
                cp * log(disc) * disc * strike * nd2 / mat);
      break;

    case PREMIUM: /*** PREMIUM ***/
      answer = disc * cp * (fwd_price * nd1 - strike * nd2);
      break;

    case RHO: /*** RHO ***/
      answer = disc * cp * strike * mat * nd2;
      break;

    case RHORHO:
      answer = -disc * cp * strike * mat * mat * nd2 +
               cp * disc * mat * strike * gauss(d2) * (sqrt(mat) / vol);
      break;

    case RHODELTA:
      answer =
          -cp * mat * disc * nd1 + cp * disc * gauss(d1) * (sqrt(mat) / vol);
      break;

    case VANNA:
      answer = -disc * d2 * gauss(d1) / vol;
      break;

    case VOLGA:
      answer = disc * fwd_price * sqrt(mat) * d2 * d1 * gauss(d1) / vol;
      break;

    default:
      answer = UNKNOWN_GREEK;
      break;
    }
  }

  return (answer);

} /* END srt_f_optblksch_accurate() */

/* A Higher level interface function to call Black Scholes */
Err OptBlkSchAccurate(double fwd_price, double strike, double vol, double mat,
                      double disc, char *call_put_str, char *greek_str,
                      double *result) {
  Err err = NULL;
  SrtCallPutType call_put_type;
  SrtGreekType greek_type = PREMIUM;

  /* Transforms the call_put string into a type */
  err = interp_call_put(call_put_str, &call_put_type);
  if (err)
    return err;

  /* Transforms the greek string into a type */
  if (strcmp(greek_str, "") != 0) {
    err = interp_greeks(greek_str, &greek_type);
    if (err)
      return err;
  }

  /* Check the inputs */
  err = srt_f_optchkblksch(fwd_price, strike, vol, mat);
  if (err)
    return err;

  /* Compute the requested greek */
  *result = srt_f_optblksch_accurate(fwd_price, strike, vol, mat, disc,
                                     call_put_type, greek_type);

  /* Return a success message */
  return NULL;
}
