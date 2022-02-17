/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optonetou()
 *
 * PURPOSE     	: one touch options for stocks  , bonds & swaps
 *
 * DESCRIPTION  	: XX
 *
 * CALLS		: norm() in gen_math.c
 *
 * PARAMETERS  	: fwd_price 	- forward price of underlying
 *             	: spot     	- spot price of underlying
 *             	: strike   	- strike price
 *             	: barrier   	- barrier of the "touch"
 *             	: vol         	- volatility of underlying
 *             	: mat        	- maturity of option
 *             	: disc       	- discount factor
 *             	: call_put      - type of option: 0 call  , 1 put
 *             	: greek   - premium (0 ONLY !!!)
 *
 * RETURNS      	: premium    	- price of the option
 *
 *******************************************************************************/

double srt_f_optonetou(double fwd_price, double spot, double strike,
                       double barrier, double vol, double mat, double disc,
                       SrtCallPutType call_put, SrtGreekType greek) {

  /**** ONE_TOUCH payoff : (barrier-strike) if barrier reached	****/
  /**** 		      	 max(spot-strike  ,0) 	otherwise  	****/
  /**** 			(in the case of the call)		****/
  /**** In other words  , the payoff can be locked at a certain 	****/
  /**** level  , but is also capped at that same level.		****/

  double premium, chi, answer, shift;
  double cp;
  double x1, x2, y1, y2, xx1, xx2, yy1, yy2;
  double ndx1, ndx2, ndy1, ndy2, ndxx1, ndxx2, ndyy1, ndyy2;

  cp = (call_put == SRT_CALL) ? 1 : -1;

  if (mat <= 0.0) {
    if (spot >= barrier)
      premium = cp * (barrier - strike);
    else {
      if (cp * (spot - strike) > 0)
        premium = cp * (spot - strike);
      else
        premium = 0.0;
    }
  } else if (cp * (spot - barrier) >= 0) {
    premium = cp * disc * (barrier - strike);
  } else if (vol == 0) {
    if ((cp * (fwd_price - barrier) < 0) && (cp * (fwd_price - strike) > 0)) {
      premium = cp * (fwd_price - strike) * disc;
    } else {
      if (cp * (fwd_price - barrier) > 0)
        premium = cp * (barrier - strike) * disc;
      else
        premium = 0.0;
    }
  } else {
    chi = (-log(spot / fwd_price) / (vol * vol * mat)) + 0.5;

    x1 = (log(fwd_price / strike) + vol * vol / 2 * mat) / (vol * sqrt(mat));
    x2 = x1 - vol * sqrt(mat);
    ndx1 = norm(x1);
    ndx2 = norm(x2);

    y1 = (log((barrier * barrier * fwd_price) / (spot * spot * strike)) +
          vol * vol / 2 * mat) /
         (vol * sqrt(mat));
    y2 = y1 - vol * sqrt(mat);
    ndy1 = norm(y1);
    ndy2 = norm(y2);

    xx1 = (log(fwd_price / barrier) + vol * vol / 2 * mat) / (vol * sqrt(mat));
    xx2 = xx1 - vol * sqrt(mat);
    ndxx1 = norm(xx1);
    ndxx2 = norm(xx2);

    yy1 = (log((barrier * fwd_price) / (spot * spot)) + vol * vol / 2 * mat) /
          (vol * sqrt(mat));
    yy2 = yy1 - vol * sqrt(mat);
    ndyy1 = norm(yy1);
    ndyy2 = norm(yy2);

    if (call_put == SRT_CALL) /*  if call */
    {

      if (spot <= barrier && strike <= barrier) {
        /*     up & out     */
        premium = fwd_price * (ndx1 - ndxx1) - strike * (ndx2 - ndxx2) -
                  fwd_price * pow(spot / barrier, -2 * chi) * (ndy1 - ndyy1) +
                  strike * pow(spot / barrier, -2 * chi + 2) *
                      (ndy2 - ndyy2) /******MODIFIED***/
                  +
                  (barrier - strike) *
                      (ndxx2 + pow(spot / barrier, -2 * chi + 2) * (1 - ndyy2));
        premium *= disc;

      } else if (spot > barrier && strike <= barrier)
        premium = (barrier - strike) * disc;
      else
        premium = 0.0;
    } else if (call_put == SRT_PUT) /*  if put  */
    {
      if (spot >= barrier && strike >= barrier) {
        /*     down  & out     */
        premium = -fwd_price * (1 - ndx1) + strike * (1 - ndx2) +
                  fwd_price * (1 - ndxx1) - strike * (1 - ndxx2) -
                  fwd_price * pow(spot / barrier, -2 * chi) * (ndy1 - ndyy1) +
                  strike * pow(spot / barrier, -2 * chi + 2) * (ndy2 - ndyy2) +
                  (strike - barrier) *
                      (1 - ndxx2 + pow(spot / barrier, -2 * chi + 2) * ndyy2);
        premium *= disc;
      } else if (spot < barrier && strike >= barrier)
        premium = (strike - barrier) * disc;
      else
        premium = 0.0;
    }
  }

  if (greek == PREMIUM) {
    answer = premium;
  } else if (greek == DELTA_FWD) { /*   DELTA  FWD UP*/
    shift = fwd_price / 1000;
    answer = (srt_f_optonetou(fwd_price + shift, spot, strike, barrier, vol,
                              mat, disc, call_put, PREMIUM) -
              premium) /
             shift;
  } else if (greek == DELTA) { /*   DELTA  SPOT DOWN*/
    shift = spot / 1000;
    answer = (premium - srt_f_optonetou(fwd_price - shift, spot - shift, strike,
                                        barrier, vol, mat, disc, call_put,
                                        PREMIUM)) /
             shift;
  } else if (greek == GAMMA_FWD) { /*   GAMMA  FWD*/
    shift = fwd_price / 1000;
    answer = (srt_f_optonetou(fwd_price + shift, spot, strike, barrier, vol,
                              mat, disc, call_put, PREMIUM) +
              srt_f_optonetou(fwd_price - shift, spot, strike, barrier, vol,
                              mat, disc, call_put, PREMIUM) -
              2 * premium) /
             (shift * shift);
  } else if (greek == GAMMA) { /*   GAMMA  */
    shift = spot / 1000;
    answer = (srt_f_optonetou(fwd_price + shift, spot + shift, strike, barrier,
                              vol, mat, disc, call_put, PREMIUM) +
              srt_f_optonetou(fwd_price - shift, spot - shift, strike, barrier,
                              vol, mat, disc, call_put, PREMIUM) -
              2 * premium) /
             (shift * shift);
  } else if (greek == VEGA) { /*   VEGA    */
    shift = GVOPT.vol_add;
    answer = (srt_f_optonetou(fwd_price, spot, strike, barrier, vol + shift,
                              mat, disc, call_put, PREMIUM) -
              premium);
  } else if (greek == THETA) { /*   THETA   */
    shift = YEARS_IN_DAY;
    answer = (srt_f_optonetou(fwd_price, spot, strike, barrier, vol,
                              mat - shift, disc * exp(-shift * log(disc) / mat),
                              call_put, PREMIUM) -
              premium);
  } else
    answer = UNKNOWN_GREEK;

  return (answer);

} /* srt_f_optonetou() */

/******************************************************************************/
