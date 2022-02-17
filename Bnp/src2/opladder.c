/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optladder()
 *
 * PURPOSE     	: ladder options for stocks  , bonds & swaps
 *
 * DESCRIPTION  	: XX
 *
 * CALLS		: srt_f_optrebate()
 *		: srt_f_optonetou()
 *		: srt_f_optblksch()
 *
 * PARAMETERS  	: fwd_price 	- forward price of underlying
 *             	: spot     	- spot price of underlying
 *             	: strike   	- strike price
 *             	: barrier   	- barrier ??
 *             	: vol         	- volatility of underlying
 *             	: mat        	- maturity of option
 *             	: disc       	- discount factor
 *             	: call_put  	- type of option: 0 call  , 1 put
 *             	: grad_best     - type of ladder: best or graduating
 *             	: greek   - premium (0 ONLY !!!)
 *
 * RETURNS      	: premium    	- price of the option
 *
 *******************************************************************************/

double srt_f_optladder(double fwd_price, double spot, double strike,
                       double hist_ext, double rungs[], double vol, double mat,
                       double disc, SrtCallPutType call_put,
                       SrtLadderType grad_best, SrtGreekType greek) {
  /* BEGIN srt_f_optladder() */

  int i;
  int num_rungs;
  SrtBarrierType down_up;
  int first_rung;
  int cp, gb;
  double premium = 0.0;
  double answer = 0.0;
  double shift;
  double du = 0.0;
  double a = 0.0;

  cp = (call_put == SRT_CALL) ? 1 : -1;
  gb = (grad_best == SRT_GRADUATED) ? 0 : 1;

  /***
          deals with historical extreme
          finds whether the option is increasing
          or decreasing  , depending on the order of rungs
  ***/

  if (call_put == SRT_CALL && (spot > hist_ext))
    hist_ext = spot;
  else if (call_put == SRT_PUT && (spot < hist_ext))
    hist_ext = spot;

  for (i = 1; rungs[i] > 0 && i <= MAX_RUNG; i++)
    ; /* set number of rungs */
  num_rungs = i - 1;

  if (rungs[1] <= strike && call_put == SRT_CALL)
    return (0.0);
  else if (rungs[1] >= strike && call_put == SRT_PUT)
    return (0.0);

  if (num_rungs == 1) {
    down_up = (call_put == SRT_CALL) ? SRT_UP : SRT_DOWN;
  } else {
    for (i = 2; i <= num_rungs; i++) /* find if rungs are in increasing */
    {                                /* or decreasing order */
      if ((rungs[i] - rungs[i - 1]) > 0)
        a++;
      else if ((rungs[i] - rungs[i - 1]) < 0)
        a--;
    }

    if (a == (double)(num_rungs - 1))
      down_up = SRT_UP; /*if increasing  , then up barrier*/
    else if (a == (double)-(num_rungs - 1))
      down_up = SRT_DOWN; /*if decreasing  , then down barrier*/
    else
      return (0.0); /*if rungs are not ordered*/

    if (((down_up == SRT_DOWN) &&
         (call_put == SRT_CALL)) /*  if decreasing  , then put*/
        || ((down_up == SRT_UP) &&
            (call_put == SRT_PUT))) /*  if increasing  , then call*/
      return (0.0);
  }

  du = (down_up == SRT_DOWN) ? 1 : -1;
  first_rung = 1;

  for (i = 1; i <= num_rungs; i++) {    /* look for the rung */
    if (du * hist_ext >= du * rungs[i]) /* reached by the spot */
      first_rung++;                     /* given historical extrema */
  }

  if (first_rung > num_rungs) {
    premium = disc * cp * (rungs[num_rungs] - strike) +
              gb * srt_f_optblksch(fwd_price, rungs[num_rungs], vol, mat, disc,
                                   call_put, PREMIUM);
  } else {
    if (first_rung == 1) {
      premium = gb * srt_f_optonetou(fwd_price, spot, strike, rungs[first_rung],
                                     vol, mat, disc, call_put, PREMIUM) +
                (1 - gb) * srt_f_optrebate(fwd_price, spot, rungs[first_rung],
                                           cp * (rungs[first_rung] - strike),
                                           vol, mat, disc, call_put, PREMIUM);

      for (i = first_rung; i <= num_rungs - 1; i++) {
        premium +=
            gb * srt_f_optonetou(fwd_price, spot, rungs[i], rungs[i + 1], vol,
                                 mat, disc, call_put, PREMIUM) +
            (1 - gb) * srt_f_optrebate(fwd_price, spot, rungs[i + 1],
                                       cp * (rungs[i + 1] - rungs[i]), vol, mat,
                                       disc, call_put, PREMIUM);
      }

      premium += gb * srt_f_optblksch(fwd_price, rungs[num_rungs], vol, mat,
                                      disc, call_put, PREMIUM);

    } else {
      premium = disc * cp * (rungs[first_rung - 1] - strike);
      for (i = first_rung - 1; i <= num_rungs - 1; i++) {
        premium +=
            gb * srt_f_optonetou(fwd_price, spot, rungs[i], rungs[i + 1], vol,
                                 mat, disc, call_put, PREMIUM) +
            (1 - gb) * srt_f_optrebate(fwd_price, spot, rungs[i + 1],
                                       cp * (rungs[i + 1] - rungs[i]), vol, mat,
                                       disc, call_put, PREMIUM);
      }

      premium += gb * srt_f_optblksch(fwd_price, rungs[num_rungs], vol, mat,
                                      disc, call_put, PREMIUM);
    }
  }

  if (greek == PREMIUM) {
    answer = premium;
  } else if (greek == DELTA_FWD) { /* DELTA FWD */
    shift = fwd_price / 1000;
    answer = (srt_f_optladder(fwd_price + shift, spot, strike, hist_ext, rungs,
                              vol, mat, disc, call_put, grad_best, PREMIUM) -
              premium) /
             shift;
  } else if (greek == DELTA) { /* DELTA */
    shift = spot / 1000;
    answer = (srt_f_optladder(fwd_price * (1 + shift / spot), spot + shift,
                              strike, hist_ext, rungs, vol, mat, disc, call_put,
                              grad_best, PREMIUM) -
              premium) /
             shift;
  } else if (greek == GAMMA_FWD) { /*   GAMMA FWD  */
    shift = fwd_price / 1000;
    answer = (srt_f_optladder(fwd_price + shift, spot, strike, hist_ext, rungs,
                              vol, mat, disc, call_put, grad_best, PREMIUM) +
              srt_f_optladder(fwd_price - shift, spot, strike, hist_ext, rungs,
                              vol, mat, disc, call_put, grad_best, PREMIUM) -
              2 * premium) /
             (shift * shift);
  } else if (greek == GAMMA) { /*   GAMMA */
    shift = spot / 1000;
    answer =
        (srt_f_optladder(fwd_price + shift, spot + shift, strike, hist_ext,
                         rungs, vol, mat, disc, call_put, grad_best, PREMIUM) +
         srt_f_optladder(fwd_price - shift, spot - shift, strike, hist_ext,
                         rungs, vol, mat, disc, call_put, grad_best, PREMIUM) -
         2 * premium) /
        (shift * shift);
  } else if (greek == VEGA) { /*   VEGA    */
    shift = GVOPT.vol_add;
    answer =
        (srt_f_optladder(fwd_price, spot, strike, hist_ext, rungs, vol + shift,
                         mat, disc, call_put, grad_best, PREMIUM) -
         premium);
  } else if (greek == THETA) { /*   THETA   */
    shift = YEARS_IN_DAY;
    answer = (srt_f_optladder(fwd_price, spot, strike, hist_ext, rungs, vol,
                              mat - shift, disc * exp(-shift * log(disc) / mat),
                              call_put, grad_best, PREMIUM) -
              premium);
  } else
    answer = UNKNOWN_GREEK;

  return (answer);

} /* END srt_f_optladder() */

/******************************************************************************/
