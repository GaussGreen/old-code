/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optrebate()
 *
 * PURPOSE     	: rebate options for stocks, bonds & swaps
 *
 * DESCRIPTION  	: XX
 *
 * CALLS		: norm() in gen_math.c
 *
 * PARAMETERS  	: fwd_price 	- forward price of underlying
 *             	: spot     	- spot price of underlying
 *             	: barrier   	- barrier of the option
 *             	: rebate   	- rebate value when barrier reached
 *             	: vol         	- volatility of underlying
 *             	: mat        	- maturity of option
 *             	: disc       	- discount factor
 *             	: call_put      - type of option: 0 call, 1 put
 *             	: greek   - premium (0 ONLY!!!)
 *
 * RETURNS      	: premium    	- price of the option
 *
 *******************************************************************************/

double srt_f_optrebate(
    double         fwd_price,
    double         spot,
    double         barrier,
    double         rebate,
    double         vol,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    SrtGreekType   greek)
{ /** BEGIN srt_f_optrebate() **/

    /****  REBATE payoff : (barrier - strike) if barrier reached	****/
    /**** 		       0  otherwise  (in the case of the call) 	****/
    /**** In other words, this product is similar to a 		****/
    /**** path-dependent digital.					****/

    double premium, chi, answer, shift, cp;
    double xx1, xx2, yy1, yy2;
    double ndxx1, ndxx2, ndyy1, ndyy2;
    double caux, paux;
    double nprimexx1, nprimeyy1;

    cp = (call_put == SRT_CALL) ? 1 : -1;

    if (vol != 0 && mat != 0)
    {
        chi = (-log(spot / fwd_price) / (vol * vol * mat)) + 0.5;

        xx1       = (log(fwd_price / barrier) + vol * vol / 2 * mat) / (vol * sqrt(mat));
        xx2       = xx1 - vol * sqrt(mat);
        ndxx1     = norm(xx1);
        ndxx2     = norm(xx2);
        nprimexx1 = exp(-(xx1 * xx1) / 2) / sqrt(2 * SRT_PI);

        yy1 =
            (log((barrier * fwd_price) / (spot * spot)) + vol * vol / 2 * mat) / (vol * sqrt(mat));
        yy2       = yy1 - vol * sqrt(mat);
        ndyy1     = norm(yy1);
        ndyy2     = norm(yy2);
        nprimeyy1 = exp(-(yy1 * yy1) / 2) / sqrt(2 * SRT_PI);

        caux = disc * rebate *
               (nprimexx1 * (barrier / fwd_price + pow(spot / fwd_price, -2 * chi + 2)) /
                    (vol * sqrt(mat) * fwd_price) +
                ndxx2 * pow(spot / barrier, -2 * chi + 2) *
                    ((-2 * chi + 2) / spot +
                     2 * vol / vol / mat * (1 / spot - 1 / fwd_price) * log(spot / barrier)));

        paux = disc * rebate *
               (-nprimexx1 * (barrier / fwd_price + pow(spot / fwd_price, -2 * chi + 2)) /
                    (vol * sqrt(mat) * fwd_price) -
                ndxx2 * pow(spot / barrier, -2 * chi + 2) *
                    ((-2 * chi + 2) / spot +
                     2 * vol / vol / mat * (1 / spot - 1 / fwd_price) * log(spot / barrier)));
    };

    if (cp * (spot - barrier) > 0)
    {
        premium = rebate;
    }
    else if (mat == 0.0)
    {
        premium = 0.0;
    }
    else if (vol == 0)
    {
        if (cp * (fwd_price - barrier) > 0)
            premium = rebate * disc;
        else
            premium = 0.0;
    }
    else if (call_put == SRT_CALL) /*  if call */
    {
        if (spot <= barrier)
        {
            premium = ndxx2 + pow(spot / barrier, -2 * chi + 2) * (1 - ndyy2);
            premium *= disc * rebate;
        }
        else
        {
            premium = disc * rebate;
        };
    }
    else if (call_put == SRT_PUT) /*  if put */
    {
        if (spot >= barrier)
        {
            premium = 1 - ndxx2 + pow(spot / barrier, -2 * chi + 2) * ndyy2;
            premium *= disc * rebate;
        }
        else
        {
            premium = disc * rebate;
        }
    };

    if (greek == PREMIUM)
    {
        answer = premium;
    }
    else if (greek == DELTA_FWD)
    { /* DELTA FWD */
        shift  = fwd_price / 10000;
        answer = (srt_f_optrebate(
                      fwd_price + shift, spot, barrier, rebate, vol, mat, disc, call_put, PREMIUM) -
                  premium) /
                 shift;
    }
    else if (greek == DELTA)
    { /* DELTA */
        shift  = spot / 10000;
        answer = (srt_f_optrebate(
                      fwd_price + shift,
                      spot + shift,
                      barrier,
                      rebate,
                      vol,
                      mat,
                      disc,
                      call_put,
                      PREMIUM) -
                  premium) /
                 shift;
    }
    else if (greek == GAMMA_FWD)
    { /* GAMMA FWD */
        shift  = fwd_price / 10000;
        answer = (srt_f_optrebate(
                      fwd_price + shift, spot, barrier, rebate, vol, mat, disc, call_put, PREMIUM) +
                  srt_f_optrebate(
                      fwd_price - shift, spot, barrier, rebate, vol, mat, disc, call_put, PREMIUM) -
                  2 * premium) /
                 (shift * shift);
    }
    else if (greek == GAMMA)
    { /* GAMMA */
        shift  = spot / 10000;
        answer = (srt_f_optrebate(
                      fwd_price + shift,
                      spot + shift,
                      barrier,
                      rebate,
                      vol,
                      mat,
                      disc,
                      call_put,
                      PREMIUM) +
                  srt_f_optrebate(
                      fwd_price - shift,
                      spot - shift,
                      barrier,
                      rebate,
                      vol,
                      mat,
                      disc,
                      call_put,
                      PREMIUM) -
                  2 * premium) /
                 (shift * shift);
    }
    else if (greek == VEGA)
    { /* VEGA */
        shift = 0.01;
        answer =
            (srt_f_optrebate(
                 fwd_price, spot, barrier, rebate, vol + shift, mat, disc, call_put, PREMIUM) -
             premium);
    }
    else if (greek == THETA)
    { /* THETA */
        shift = 1 / 365;
        answer =
            (srt_f_optrebate(
                 fwd_price,
                 spot,
                 barrier,
                 rebate,
                 vol,
                 mat - shift,
                 disc * exp(-shift * log(disc) / mat),
                 call_put,
                 PREMIUM) -
             premium);
    }

    return (answer);

} /* END srt_f_optrebate() */

/******************************************************************************/
