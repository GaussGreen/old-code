/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_opttftame()
 *
 * PURPOSE      	: Calculates the premium and the greeks for a American
 *		  	spread option with strike
 *
 * DESCRIPTION  	: Pricing is accomplished using a two-factor extension
 *		  	of binomial tree
 *
 * CALLS		: <none>
 *
 * PARAMETERS
 *	INPUT	: fwd1        	- forward price of 1st underlying
 *              	: fwd2         	- forward price of 2nd underlying
 *              	: strike        - strike price
 *              	: a            	- gearing for 1st underlying
 *              	: b           	- gearing for 2nd underlying
 *              	: vol1         	- annual volatility of 1st underlying
 *              	: vol2         	- annual volatility of 2nd underlying
 *              	: rho12         - correlation between underlyings
 *              	: mat           - maturity  , in years
 *              	: disc_factor   - discount factor up to maturity
 *              	: call_put      - option type: 0 call  , 1 put
 *              	: step          - number of steps in tree (< 58)
 *		: greek	_ info wanted (premium  , greeks)
 *
 * RETURNS      	: ??            - ??
 *
 *******************************************************************************/

double srt_f_opttftame(double fwd1, double fwd2, double spot1, double spot2,
                       double strike, double a, double b, double vol1,
                       double vol2, double rho12, double mat,
                       double disc_factor, SrtCallPutType call_put, int step,
                       SrtGreekType greek) {
  double rr1;             /* (1 + r) ^ (mat / n)       */
  double rr2;             /* (1 + r) ^ (mat / n)       */
  double small_disc_fact; /* DF fir one step 	     */
  double h;               /* maturity / step           */
  double u1;              /* Up move size asset 1      */
  double u2;              /* Up move size asset 2      */
  double d1;              /* Down move size asset 1    */
  double d2;              /* Down move size asset 2    */
  double p1;              /* Probability of u1  , u2  , u3 */
  double p2;              /* Probability of u1  , u2  , d3 */
  double p3;              /* Probability of u1  , d2  , u3 */
  double p4;              /* Probability of u1  , d2  , d3 */
  double pp1;             /* Probability of u1         */
  double pp2;             /* Probability of u2         */
  double asset_1;         /* [MAX_STEP]              */
                          /* Temporary asset 1       */
                          /* Temporary asset 2       */
  double *prem_i0;        /* Temporary premium [i  ] */
  double *prem_i1;        /* Temporary premium [i+1] */
  double *asset1;
  double *asset2;

  double intrinsic;
  double act_premium;
  double answer;
  double **premium;
  int i, j;
  int width_matrix;
  int cp;
  Greekstruct greeks;

  /* step must be within [3  ,MAX_STEP] & mat > 0 */
  /* size of the matrix [0......MAX_STEP]       */

  if (mat <= 0) {
    return (0);
  }

  h = mat / step; /* bug if step = 0 or <= 0 */
  small_disc_fact = pow(disc_factor, 1 / step);

  u1 = exp(vol1 * sqrt(h)); /* bug if h <= 0 */
  u2 = exp(vol2 * sqrt(h));

  d1 = 1 / u1;
  d2 = 1 / u2;

  rr1 = (1 / mat) * log(fwd1 / spot1);
  rr2 = (1 / mat) * log(fwd2 / spot2);

  if (u1 == d1) {
    pp1 = 0.5;
  } else {
    /* bug if vol1 = 0 */
    pp1 = (exp(rr1 * h) - d1) / (u1 - d1);
  }

  if (u2 == d2) {
    pp2 = 0.5;
  } else {
    pp2 = (exp(rr2 * h) - d2) / (u2 - d2);
  } /* bug if vol2 = 0 */

  p1 = pp1 * pp2 + rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));
  p2 = pp1 * (1 - pp2) - rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));
  p3 = (1 - pp1) * pp2 - rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));
  p4 = (1 - pp1) * (1 - pp2) + rho12 * u1 * u2 / ((u1 + 1) * (u2 + 1));

  /**
  for (i=0;i<step;i++)
           for (j=0;j<step;j++)
                  premium[i][j]=0;
  **/

  /* ==========================================================================
   */
  /* Initialisation of Matrix */
  /* ==========================================================================
   */

  cp = (call_put == SRT_CALL) ? 1 : -1;

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Add Two Steps
   * !!!!!!!!!!!!!!!!!!!!!!!!!*/

  step += 2;

  /* Allocation for the memory */
  if ((asset1 = dvector(0, step)) == NULL)
    return MEMORY_ERR;
  if ((asset2 = dvector(0, step)) == NULL) {
    free_dvector(asset1, 0, step);
    return MEMORY_ERR;
  }
  if ((premium = dmatrix(0, step, 0, step)) == NULL) {
    free_dvector(asset1, 0, step);
    free_dvector(asset2, 0, step);
    return MEMORY_ERR;
  }

  /* Computation */
  asset1[0] = spot1 * pow(u1, step);
  asset2[0] = spot2 * pow(u2, step);

  for (i = 1; i <= step; i++) {
    asset1[i] = asset1[i - 1] * d1 * d1;
    asset2[i] = asset2[i - 1] * d2 * d2;
  }

  /* for(j=0; j<=step; j++) {asset2[j] = fwd2 * pow(u2  ,step-j)*pow(d2  ,j); }
   */

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
  /* Dynamic Calculation in case of American
   */
  /* ==========================================================================
   */

  for (width_matrix = step - 1; width_matrix >= 2; width_matrix--) {
    for (j = 0; j <= width_matrix; j++) {
      asset2[j] *= d2;
    }

    for (i = 0; i <= width_matrix; i++) {
      prem_i0 = premium[i];
      prem_i1 = premium[i + 1];
      asset1[i] *= d1;
      asset_1 = asset1[i];

      for (j = 0; j <= width_matrix; j++) {
        premium[i][j] = p1 * prem_i0[j] +
                        /* = premium[i][j]  */
                        p2 * prem_i0[j + 1] +
                        /* = premium[i][j+1] */
                        p3 * prem_i1[j] +
                        /* = premium[i+1][j] */
                        p4 * prem_i1[j + 1];
        /* = premium[i+1][j+1] */

        premium[i][j] *= small_disc_fact;

        intrinsic = ((a * asset_1) - (b * asset2[j]) - strike) * cp;

        if (intrinsic > premium[i][j]) {
          premium[i][j] = intrinsic;
        }
      }
    }
  }

  act_premium = premium[1][1];

  /* ==========================================================================
   */
  /* Calculation of delta   	   	      				      */
  /* ==========================================================================
   */

  /*
  ( premium(uu  , ud  , ud) - premium(dd  , ud  , ud) ) / ((uu * fwd) - (dd *
  fwd)) )
  */

  if ((u1 != d1) && (fwd1 != 0)) {
    greeks.fwd_delta_underlying =
        (p1 + p3) * (premium[0][0] - premium[1][0]) / ((u1 - d1) * fwd1) +
        (p2 + p4) * (premium[0][1] - premium[1][1]) / ((u1 - d1) * fwd1);
    greeks.spot_delta_underlying = greeks.fwd_delta_underlying / disc_factor;

  } else {
    greeks.fwd_delta_underlying = 0;
  }

  greeks.fwd_gamma = 0;
  greeks.vega = 0;
  greeks.theta = 0;

  /* ==========================================================================
   */
  /* Calculation of gamma */
  /* ==========================================================================
   */

  switch (greek) {
  case PREMIUM: /* premium               */
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
  case VEGA: /* vega                  */
    answer = greeks.vega;
    break;
  case THETA: /* theta                 */
    answer = greeks.theta;
    break;
  default:
    answer = UNKNOWN_GREEK;
    break;
  }
  /* ==========================================================================
   */
  /* Calculation of premium   	   	      				      */
  /* ==========================================================================
   */

  /* Free the memory */
  free_dvector(asset1, 0, step);
  free_dvector(asset2, 0, step);
  free_dmatrix(premium, 0, step, 0, step);

  return (answer);

} /* srt_f_opttftame() */

/******************************************************************************/
