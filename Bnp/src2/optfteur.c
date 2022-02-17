/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_opttfteur()
 *
 * PURPOSE      	: Calculates the premium and the greeks for a European
 *		  	spread option with strike
 *
 * DESCRIPTION  	: Pricing is accomplished using a two-factor extension
 *		  	of binomial tree
 *
 * CALLS		: <none>
 *
 * PARAMETERS
 *	INPUT	: fwd1        	- forward price of 1st underlying
 *              	: fwd2          - forward price of 2nd underlying
 *              	: strike        - strike price
 *              	: a             - gearing for 1st underlying
 *              	: b             - gearing for 2nd underlying
 *              	: vol1      	- annual volatility of 1st underlying
 *              	: vol2          - annual volatility of 2nd underlying
 *              	: rho12      	- correlation between underlyings
 *              	: mat           - maturity  , in years
 *              	: disc_factor   - discount factor up to maturity
 *              	: call_put     	- option type: 0 call  , 1 put
 *              	: step          - number of steps in tree (< 58)
 *		: greek	- info wanted: premium  , greeks...
 *
 * RETURNS      	: ??		- ??
 *
 *******************************************************************************/

double srt_f_opttfteur(double fwd1, double fwd2, double strike, double a,
                       double b, double vol1, double vol2, double rho12,
                       double mat, double disc_factor, SrtCallPutType call_put,
                       int step, SrtGreekType greek) {
  struct greek_struct greeks;

  double rr;  /* (1+r)^(mat/n)              */
  double h;   /* maturity/step	      */
  double u1;  /* Up move size asset 1       */
  double u2;  /* Up move size asset 2       */
  double d1;  /* Down move size asset 1     */
  double d2;  /* Down move size asset 2     */
  double p1;  /* Probability of u1  , u2  , u3  */
  double p2;  /* Probability of u1  , u2  , d3  */
  double p3;  /* Probability of u1  , d2  , u3  */
  double p4;  /* Probability of u1  , d2  , d3  */
  double pp1; /* Probability of u1          */
  double pp2; /* Probability of u2          */

  double asset_1;  /* [MAX_STEP]                 */
                   /* Temporary asset 1 	      */
  double *prem_i0; /* Temporary premium [i]      */
  double *prem_i1; /* Temporary premium [i+1]    */
  double asset1[100];
  double asset2[100];

  double intrinsic;
  double act_premium;
  double premium[60][60];
  double answer;
  int i, j;
  int width_matrix;
  int cp;

  /* ==========================================================================
     Step must be within [3  , MAX_STEP-3] & mat > 0
     Size of the matrix [0......MAX_STEP-1]
     ==========================================================================
   */

  if (step > (MAX_STEP - 3) || step < 3) {
    return (0);
  }

  if (mat <= 0) {
    return (0);
  }

  h = mat / step; /* bug if step = 0 or <= 0    */

  u1 = exp(vol1 * sqrt(h)); /* bug if h <= 0	      */
  u2 = exp(vol2 * sqrt(h));

  d1 = 1 / u1;
  d2 = 1 / u2;

  rr = 1;

  if (u1 == d1) {
    pp1 = 0.5;
  } else {
    /* bug if vol1 = 0 */
    pp1 = (rr - d1) / (u1 - d1);
  }

  if (u2 == d2) {
    pp2 = 0.5;
  } else {
    pp2 = (rr - d2) / (u2 - d2);
  } /* bug if vol2 = 0 */

  /* ==========================================================================
   */
  /* Evaluate Probabilities */
  /* ==========================================================================
   */

  p1 = pp1 * pp2 + rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));
  p2 = pp1 * (1 - pp2) - rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));
  p3 = (1 - pp1) * pp2 - rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));
  p4 = (1 - pp1) * (1 - pp2) + rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));

  /* ==========================================================================
   */
  /* Initialisation of matrix	      	      				      */
  /* ==========================================================================
   */

  cp = (call_put == SRT_CALL) ? 1 : -1;

  for (i = 0; i <= step; i++) {
    asset1[i] = fwd1 * pow(u1, step - i) * pow(d1, i);
  }

  for (j = 0; j <= step; j++) {
    asset2[j] = fwd2 * pow(u2, step - j) * pow(d2, j);
  }

  for (i = 0; i <= step; i++) {
    asset_1 = asset1[i];
    prem_i0 = premium[i];

    for (j = 0; j <= step; j++) {
      intrinsic = ((a * asset_1) - (b * asset2[j]) - strike) * cp;

      if (intrinsic > 0)
        prem_i0[j] = intrinsic; /* = premium[i][j] */
      else
        prem_i0[j] = 0;
    }
  }

  /* ==========================================================================
   */
  /* Dynamic calculation in case of european      	 		      */
  /* ==========================================================================
   */

  for (width_matrix = step - 1; width_matrix >= 1; width_matrix--) {
    for (i = 0; i <= width_matrix; i++) {
      prem_i0 = premium[i];
      prem_i1 = premium[i + 1];

      for (j = 0; j <= width_matrix; j++) {
        premium[i][j] = p1 * prem_i0[j] +
                        /* = premium[i][j]     */
                        p2 * prem_i0[j + 1] +
                        /* = premium[i][j+1]   */
                        p3 * prem_i1[j] +
                        /* = premium[i+1][j]   */
                        p4 * prem_i1[j + 1];
        /* = premium[i+1][j+1] */
      }
    }
  }

  /* ==========================================================================
   */
  /* Discounting after the premium */
  /* ==========================================================================
   */

  for (i = 0; i <= 1; i++) {
    for (j = 0; j <= 1; j++) {
      premium[i][j] = disc_factor * premium[i][j];
    }
  }
  act_premium = (p1 * premium[0][0]) + (p2 * premium[0][1]) +
                (p3 * premium[1][0]) + (p4 * premium[1][1]);

  /* ==========================================================================
   */
  /* Calculation of delta */
  /* ==========================================================================
   */

  if (u1 != d1 && fwd1 != 0) {
    greeks.fwd_delta_underlying =
        (p1 + p3) * (premium[0][0] - premium[1][0]) / ((u1 - d1) * fwd1) +
        (p2 + p4) * (premium[0][1] - premium[1][1]) / ((u1 - d1) * fwd1);
  } else {
    greeks.fwd_delta_underlying = 0;
  }
  greeks.fwd_gamma = 0;
  greeks.vega = 0;
  greeks.theta = 0;

  /* ==========================================================================
   */
  /* Calculation of premium   	   	      				      */
  /* ==========================================================================
   */

  switch (greek) {
  case PREMIUM:
    answer = act_premium;
    break;
  case DELTA_FWD: /* fwd delta underlying  */
    answer = greeks.fwd_delta_underlying;
    break;
  case DELTA: /* spot delta underlying */
    answer = greeks.fwd_delta_underlying;
    break;
  case GAMMA: /* fwd gamma             */
    answer = greeks.fwd_gamma;
    break;
  case VEGA: /* vega */
    answer = greeks.vega;
    break;
  case THETA: /* theta                 */
    answer = greeks.theta;
    break;
  }

  return (answer);

} /* srt_f_opttfteur() */

/******************************************************************************/
