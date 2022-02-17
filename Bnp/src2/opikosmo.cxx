/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: iko_ext_fct_binomial3(14)
 *
 * PURPOSE      	: ??
 *
 * DESCRIPTION  	: iko + iko ext with premium smoothed
 *
 * CALLS		: comb()
 *		: smooth()
 *		: srt_f_optexting()
 *
 * PARAMETERS   	: fwd_price2    - forward price of 2nd underlying
 *              	: fwd_price1  	- forward price of 1st underlying
 *              	: spot_price   	- ??
 *              	: barrier   	- barrier level
 *              	: vol1          - volatility of 1st underlying
 *              	: vol2          - volatility of 2nd underlying
 *              	: mat1          - maturity of 1st underlying
 *              	: mat2          - maturity of 2nd underlying
 *              	: disc          - discount factor
 *              	: call_put 	- type of option: 0 call        , 1 put
 *              	: num_rungs     - number of rungs
 *              	: rungs[] 	- rungs ??
 *              	: strikes[]   	- strike prices
 *              	: step_num 	- number of steps for ??
 *
 * RETURNS      	: premium        - ??
 *
 *******************************************************************************/

/*** iko + ikoext with premium smoothed ***/

static double iko_ext_fct_binomial3(double fwd_price2, double fwd_price1,
                                    double spot_price, double barrier,
                                    double vol1, double vol2, double mat1,
                                    double mat2, double disc,
                                    SrtCallPutType call_put, int num_rungs,
                                    double rungs[], double strikes[],
                                    int step_num) {
  double disc1; /* discount factor for first time period  */
  double disc2; /* discount factor for second time period */

  double p, d, u, r; /* see hull p.75 			  */
  double cp;         /* call or put   */
  double prob;       /* number of ways to have a given highest value */
  double spot_end;   /* spot price at end of first time period */
  double strike_end; /* adjusted strike price after first time period */
  double premium[MAX_STEP + 3]; /* price for increasing strike option */
  double price[69];
  double opt;      /* black scholes value */
  double max_spot; /* max spot price during first time period */
  double rate1;    /* risk free rate at time 0 --- T1 */
  double rate2;    /* risk free rate at time  T1 --- T2 */
  double pts;

  int i;
  int j;
  int k;
  int q;
  int min;
  int n;
  double dd;
  double answer;

  cp = (call_put == SRT_CALL) ? 1 : -1;

  if (cp * (spot_price - barrier) <= 0.0)
    answer = 0.0;
  else {
    spot_price -= 13 * 0.05;

    for (k = 1; k <= 25; k++) {
      spot_price += 0.05;

      /* compute using a binary discretization */

      rate1 = log(fwd_price1 / spot_price) / mat1;
      rate2 = log(fwd_price2 / fwd_price1) / (mat2 - mat1);
      disc1 = exp(-rate1 * mat1);
      disc2 = exp(-rate2 * (mat2 - mat1));

      u = exp(vol1 * sqrt(mat1 / step_num));
      d = 1 / u;
      p = (exp(rate1 * mat1 / step_num) - d) / (u - d);
      r = exp(rate1 * mat1 / step_num);

      for (i = 0; i <= step_num; i++) {
        premium[i] = 0.0;

        /***
                iterate over different spot prices at end of first period
                spot_end is spot at end of period
                max_spot is highest the max spot value could have gone ---
                (this is at least the spot) while still ending up at
                spot_end.
        ***/

        spot_end = spot_price * pow(u, step_num - 2 * i);
        min = (step_num > (2 * i) ? i : step_num - i);
        max_spot = spot_price * pow(u, step_num - i);

        for (j = 0; j <= min; j++) {
          /***
                  iterate over different max prices during the first period
          ***/
          strike_end = strikes[1];

          for (q = 1; q <= num_rungs; q++) {
            /***
            set value of strike at end of first period
            ***/

            if ((cp * max_spot) > (cp * rungs[q])) {
              strike_end = ((cp * strikes[q]) > (cp * strikes[1]) ? strikes[q]
                                                                  : strikes[1]);
            }
          }

          dd = (double)j;
          prob = comb(dd, step_num, p, i) - comb(dd - 1, step_num, p, i);

          opt = srt_f_optexting(spot_end / disc2, spot_end, strike_end, barrier,
                                vol2, mat2 - mat1, disc, call_put, SRT_DOWN,
                                PREMIUM);

          premium[i] += prob * opt;
          max_spot /= u;
        }

        dd = (double)i;
        premium[i] /= comb(dd, step_num, p, i);
      }

      /* binomial tree with payoffs at t=T1 set at premium[step_num] obtained
       * above */

      for (i = step_num - 1; i >= 0; i--) {
        for (j = 0; j <= i; j++) {
          spot_end = spot_price * pow(u, i - j * 2);

          if (spot_end > barrier)
            premium[j] = ((p * premium[j]) + ((1 - p) * premium[j + 1]));
          else
            premium[j] = 0;
        }
      }

      price[k] = premium[0];

    } /*end of loop for smoothing */

    pts = 13.0;
    n = 10;

    answer = smooth(price, n, pts);

  } /* end of loop for test spot_barrier */

  return (answer);

} /* END iko_ext_fct_binomial3(14) */

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optikosmo()
 *
 * PURPOSE      	: ??
 *
 * DESCRIPTION  	: ??
 *
 * CALLS			: iko_ext_fct_binomial3()
 *
 * PARAMETERS
 *	INPUT		: fwd_price2  	- forward price of 2nd underlying
 *              	: fwd_price1  	- forward price of 1st underlying
 *              	: spot_price   	- ??
 *              	: barrier   	- barrier level
 *              	: vol1          - volatility of 1st underlying
 *              	: vol2          - volatility of 2nd underlying
 *              	: mat1          - maturity of 1st underlying
 *              	: mat2          - maturity of 2nd underlying
 *              	: disc          - discount factor
 *              	: call_put 	- type of option: 0 call        , 1 put
 *              	: num_rungs     - number of rungs
 *              	: rungs[]  	- rungs ??
 *              	: strikes[]   	- strike prices
 *              	: step_num 	- number of steps for the tree
 *
 * RETURNS      	: premium     	- ??
 *
 *******************************************************************************/

double srt_f_optikosmo(double fwd_price2, double fwd_price1, double spot_price,
                       double barrier, double hist_extreme, double vol1,
                       double vol2, double mat1, double mat2, double disc,
                       SrtCallPutType call_put, int num_rungs, double rungs[],
                       double strikes[], int step_num, SrtGreekType greek) {
  int i;
  int hist_rung;
  double premium;
  double answer;
  double shift;
  double cp = (call_put == SRT_CALL) ? 1 : -1;

  /***
          deal with historical extreme
  ***/

  /***
          first find what strike must have been
          reset to by past movement of the spot.
  ***/

  hist_rung = 1;

  for (i = 2; i <= num_rungs; i++) {
    if ((cp * hist_extreme) >= (cp * rungs[i]))
      hist_rung++;
  }

  /***
          hist_rung is the index of the bottom
          of the interval which corresponds to
          the new strike value.

          now change strikes corresponding to
          rungs below the historical extreme to
          reset (due to historical movement) strike value.
  ***/

  for (i = 1; i < hist_rung; i++)
    strikes[i] = strikes[hist_rung];

  premium = iko_ext_fct_binomial3(fwd_price1, fwd_price2, spot_price, barrier,
                                  vol1, vol2, mat1, mat2, disc, call_put,
                                  num_rungs, rungs, strikes, step_num);

  switch (greek) {
  case PREMIUM:
    answer = premium;
    break;

  case DELTA_FWD1:
    shift = fwd_price1 / 10000;
    answer =
        (srt_f_optikosmo(fwd_price2, fwd_price1 + shift, spot_price, barrier,
                         hist_extreme, vol1, vol2, mat1, mat2, disc, call_put,
                         num_rungs, rungs, strikes, step_num, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA_FWD2:
    shift = fwd_price2 / 10000;
    answer =
        (srt_f_optikosmo(fwd_price2 + shift, fwd_price1, spot_price, barrier,
                         hist_extreme, vol1, vol2, mat1, mat2, disc, call_put,
                         num_rungs, rungs, strikes, step_num, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA:
    shift = spot_price / 10000;
    answer = (srt_f_optikosmo(fwd_price1 * (1 + shift / spot_price),
                              fwd_price2 * (1 + shift / spot_price),
                              spot_price + shift, barrier, hist_extreme, vol1,
                              vol2, mat1, mat2, disc, call_put, num_rungs,
                              rungs, strikes, step_num, PREMIUM) -
              premium) /
             shift;
    break;

  case GAMMA:
    shift = spot_price / 10000;
    answer = srt_f_optikosmo(fwd_price1 * (1 + shift / spot_price),
                             fwd_price2 * (1 + shift / spot_price),
                             spot_price + shift, barrier, hist_extreme, vol1,
                             vol2, mat1, mat2, disc, call_put, num_rungs, rungs,
                             strikes, step_num, PREMIUM);
    answer += srt_f_optikosmo(fwd_price1 * (1 - shift / spot_price),
                              fwd_price2 * (1 - shift / spot_price),
                              spot_price - shift, barrier, hist_extreme, vol1,
                              vol2, mat1, mat2, disc, call_put, num_rungs,
                              rungs, strikes, step_num, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA1:
    shift = GVOPT.vol_add;
    answer = (srt_f_optikosmo(fwd_price2, fwd_price1, spot_price, barrier,
                              hist_extreme, vol1 + shift, vol2, mat1, mat2,
                              disc, call_put, num_rungs, rungs, strikes,
                              step_num, PREMIUM) -
              premium) /
             shift;
    break;

  case VEGA2:
    shift = GVOPT.vol_add;
    answer = (srt_f_optikosmo(fwd_price2, fwd_price1, spot_price, barrier,
                              hist_extreme, vol1, vol2 + shift, mat1, mat2,
                              disc, call_put, num_rungs, rungs, strikes,
                              step_num, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA:
    shift = YEARS_IN_DAY;
    answer =
        (srt_f_optikosmo(fwd_price2, fwd_price1, spot_price, barrier,
                         hist_extreme, vol1, vol2, mat1 - shift, mat2 - shift,
                         disc * exp(-shift * log(disc) / mat2), call_put,
                         num_rungs, rungs, strikes, step_num, PREMIUM) -
         premium) /
        shift;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);

} /* srt_f_optikosmo() */

/******************************************************************************/
