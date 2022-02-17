/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include <OPFNCTNS.H>
#include <math.h"
#include <utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optfirsec ( )
 *
 * PURPOSE      	: First Second Digital        , with payment at maturity
 *
 * DESCRIPTION  	: Pays off one at maturity if the underlying hits
 *		  barrier A first then barrier B . Barrier A can be up
 * 		  or down. Uses formula developed by Dimitri Spoliansky.
 *		  ' Probability densities & pricing for double
 *		  barrier options' 1994
 *
 * PARAMETERS
 *	INPUT		: fwd			- forward underlying price
 *              	: spot      	- spot price
 *              	: fir_bar      	- first barrier
 *              	: sec_bar      	- second barrier
 *              	: vol         	- volatility
 *              	: mat         	- maturity        , in years
 *              	: disc        	- discount fact to expiry
 *              	: dufirst       - 0 : hit down barrier first
 *                               - 1 : hit up barrier first
 *		: nterms	- number of terms in expansion
 *
 * RETURNS      	: premium       - option premium
 *
 *******************************************************************************/

double srt_f_optfirsec(double fwd, double spot, double fir_bar, double sec_bar,
                       double vol, double mat, double disc, int nterms,
                       SrtGreekType greek) {
  double rate;     /* r  					*/
  double drift;    /* mu 					*/
  double a_d, b_u; /* normalised barrier 	*/

  double vol_squared;  /* is (vol * vol)    	*/
  double vol_root_mat; /* is (vol * sqrt(mat))	*/
  double temp;
  double fact1, fact2, fact3, fact4;
  double d1, d2, d3, d4;
  double n1, n2, n3, n4;

  double bar_d, bar_u;
  int dufirst;

  double premium = 0.0, answer = 0.0, shift = 0.0;
  double probability = 0;
  int p;

  if ((mat <= 0.0) || (vol == 0.0)) {
    premium = 0.0;
  } else {
    bar_d = (fir_bar <= sec_bar) ? fir_bar : sec_bar;
    bar_u = (fir_bar >= sec_bar) ? fir_bar : sec_bar;
    dufirst = (fir_bar == bar_u) ? 1 : 0;

    vol_squared = (vol * vol);

    rate = log(fwd / spot) / mat;
    drift = rate - (vol_squared / 2);
    a_d = log(bar_d / spot);
    b_u = log(bar_u / spot);

    if (dufirst == 1) /* barrier up to be hit first */
    {
      drift = (-drift);

      temp = a_d;
      a_d = (-b_u);
      b_u = (-temp);
    }

    vol_root_mat = vol * sqrt(mat);

    for (p = 1; p < nterms + 1; p++) {
      fact1 = exp(2.0 * drift * p * (b_u - a_d) / vol_squared);
      d1 = (b_u - 2.0 * p * (b_u - a_d) - drift * mat) / vol_root_mat;
      n1 = fact1 * norm(d1);

      fact2 = exp(2.0 * drift * (b_u + p * (b_u - a_d)) / vol_squared);
      d2 = (-b_u - 2.0 * p * (b_u - a_d) - drift * mat) / vol_root_mat;
      n2 = fact2 * norm(d2);

      fact3 = exp(2.0 * drift * (b_u - p * (b_u - a_d)) / vol_squared);
      d3 = (b_u - 2.0 * p * (b_u - a_d) + drift * mat) / vol_root_mat;
      n3 = fact3 * norm(d3);

      fact4 = 1.0 / fact1;
      d4 = (-b_u - 2.0 * p * (b_u - a_d) + drift * mat) / vol_root_mat;
      n4 = fact4 * norm(d4);

      probability += n1 - n2 + n3 - n4;
    }

    premium = probability * disc;
  }

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTA_FWD: /*** DELTA FWD ***/
    shift = fwd / 10000;
    answer = (srt_f_optfirsec(fwd + shift, spot, fir_bar, sec_bar, vol, mat,
                              disc, nterms, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA: /*** DELTA SPOT ***/
    shift = spot / 10000;
    answer = (srt_f_optfirsec(fwd * (1 + shift / spot), spot + shift, fir_bar,
                              sec_bar, vol, mat, disc, nterms, PREMIUM) -
              premium) /
             shift;
    break;

  case GAMMA: /*** GAMMA ***/
    shift = spot / 10000;
    answer = srt_f_optfirsec(fwd * (1 + shift / spot), spot + shift, fir_bar,
                             sec_bar, vol, mat, disc, nterms, PREMIUM);
    answer += srt_f_optfirsec(fwd * (1 - shift / spot), spot - shift, fir_bar,
                              sec_bar, vol, mat, disc, nterms, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA: /*** VEGA ***/
    shift = GVOPT.vol_add;
    answer = (srt_f_optfirsec(fwd, spot, fir_bar, sec_bar, vol + shift, mat,
                              disc, nterms, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer =
        srt_f_optfirsec(fwd, spot, fir_bar, sec_bar, vol, mat - shift,
                        disc * exp(-shift * log(disc) / mat), nterms, PREMIUM) -
        premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);
}

/* ==========================================================================
   END srt_f_optfirsec
   ========================================================================== */
