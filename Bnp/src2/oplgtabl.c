/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optlgtabl
 *
 * PURPOSE      	: lightable barrier options
 *
 * DESCRIPTION  	: closed form on black-scholes
 *
 * CALLS			: srt_f_optblksch()
 *				  srt_f_optexting   ()
 *
 * PARAMETERS   	: fwd	 	- forward price of underlying
 *              	: spot     	- spot price of underlying
 *              	: strike    	- strike price
 *              	: barrier   	- barrier ??
 *              	: vol       	- volatility of underlying
 *              	: mat       	- maturity of underlying
 *              	: disc      	- discount factor
 *              	: call_put      - type of option: 0 call  , 1 put
 *              	: greek      	- Greek asked
 *
 * RETURNS      	: premium   	- premium_or_delta
 *
 *******************************************************************************/

double srt_f_optlgtabl(double fwd, double spot, double strike, double barrier,
                       double vol, double mat, double disc,
                       SrtCallPutType call_put, SrtBarrierType down_up,
                       SrtGreekType greek) {
  double premium, answer, shift;

  /* An extinguishable and a lightable give an european */
  premium = srt_f_optblksch(fwd, strike, vol, mat, disc, call_put, PREMIUM) -
            srt_f_optexting(fwd, spot, strike, barrier, vol, mat, disc,
                            call_put, down_up, PREMIUM);

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTA_FWD: /*** DELTA FWD ***/
    shift = fwd / 10000;
    answer = (srt_f_optlgtabl(fwd + shift, spot, strike, barrier, vol, mat,
                              disc, call_put, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA: /*** DELTA SPOT + FWD ***/
    shift = spot / 10000;
    answer =
        (srt_f_optlgtabl(fwd * (1 + shift / spot), spot + shift, strike,
                         barrier, vol, mat, disc, call_put, down_up, PREMIUM) -
         premium) /
        shift;
    break;

  case GAMMA: /*** GAMMA ***/
    shift = spot / 10000;
    answer =
        srt_f_optlgtabl(fwd * (1 + shift / spot), spot + shift, strike, barrier,
                        vol, mat, disc, call_put, down_up, PREMIUM);
    answer +=
        srt_f_optlgtabl(fwd * (1 - shift / spot), spot - shift, strike, barrier,
                        vol, mat, disc, call_put, down_up, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA: /*** VEGA ***/
    shift = GVOPT.vol_add;
    answer = (srt_f_optlgtabl(fwd, spot, strike, barrier, vol + shift, mat,
                              disc, call_put, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer = srt_f_optlgtabl(fwd, spot, strike, barrier, vol, mat - shift,
                             disc * exp(-shift * log(disc) / mat), call_put,
                             down_up, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);

} /* srt_f_optlgtabl() */

/******************************************************************************/
