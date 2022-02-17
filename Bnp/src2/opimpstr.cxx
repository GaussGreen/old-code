/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optimpstr()
 *
 * PURPOSE      	: Calculates the strike implied from a premium
 *		  	computed using the Black-Scholes formula
 *
 * DESCRIPTION  	:
 *
 * CALLS		: srt_f_optblksch()
 *
 * PARAMETERS
 *	INPUT	: premium	- premium on European option
 *              	: fwd_price    	- forward price of underlying price
 *              	: vol           - annual volatility
 *              	: mat           - maturity        , in years
 *              	: disc          - discount factor
 *              	: call_put      - type of option: 0 call        , 1 put
 *
 * RETURNS      	: strike_guess  - guess at implied strike
 *
 *******************************************************************************/
static double optvalstrike(double fwd, double strike, double vol, double mat,
                           double disc_fact, SrtCallPutType call_put,
                           SrtGreekType greek, SrtDiffusionType log_or_norm) {
  if (log_or_norm == SRT_LOGNORMAL) {
    return srt_f_optblksch(fwd, strike, vol, mat, disc_fact, call_put, greek);
  } else {
    return srt_f_optblknrm(fwd, strike, vol, mat, disc_fact, call_put, greek);
  }
}

double srt_f_optimpstr(double premium, double fwd_price, double vol, double mat,
                       double disc, SrtCallPutType call_put,
                       SrtDiffusionType log_norm) {
  int i;

  double strike_guess;
  double tmp_strike;
  double strike_init;
  double prem_shift;
  double prem_new;
  double deriv;

  strike_init = fwd_price;
  strike_guess = strike_init;

  /* ==========================================================================
   */
  /* XXX     		       		   			      	      */
  /* ==========================================================================
   */

  for (i = 0; i <= MAX_ITER; i++) {
    /* ------------------------------------------------------------------ */
    /* XXX  	      */
    /* ------------------------------------------------------------------ */

    prem_new = optvalstrike(fwd_price, strike_guess, vol, mat, disc, call_put,
                            PREMIUM, log_norm);

    if (fabs(prem_new - premium) < PREM_TOL) {
      break;
    }

    /* ------------------------------------------------------------------ */
    /* XXX  	      */
    /* ------------------------------------------------------------------ */

    prem_shift = optvalstrike(fwd_price,
                              strike_guess + (double)(STRIKE_SHIFT * fwd_price),
                              vol, mat, disc, call_put, PREMIUM, log_norm);

    deriv = (prem_shift - prem_new) / (double)(STRIKE_SHIFT * fwd_price);

    /* ------------------------------------------------------------------ */
    /* XXX  	      */
    /* ------------------------------------------------------------------ */

    if (deriv != 0) {
      /* ---------------------------------------------------------- */
      /* XXX  	      */
      /* ---------------------------------------------------------- */

      tmp_strike = strike_guess - ((prem_new - premium) / deriv);

      if (tmp_strike > (strike_guess * 10))
        strike_guess = strike_guess * 2.0;
      else if (tmp_strike < (strike_guess / 10))
        strike_guess = strike_guess / 2.0;
      else
        strike_guess = tmp_strike;
    } else {
      /* ---------------------------------------------------------- */
      /* XXX  	      */
      /* ---------------------------------------------------------- */

      if (prem_new < premium)
        strike_guess = strike_guess * 2.0;
      else
        strike_guess = strike_guess / 1.5;
    }
  }

  /* ==========================================================================
   */
  /* Return implied strike guess       	      				      */
  /* ==========================================================================
   */

  return (strike_guess);

} /* END srt_f_optimpstr() */

/******************************************************************************/
