/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optblkd12(6)
 *
 * PURPOSE      	: Calculates d1 and d2 for the Black-Scholes option
 *		  	pricing for a single risk factor
 *
 * DESCRIPTION  	: ??
 *
 * CALL		: <none>
 *
 * PARAMETERS
 *	INPUT	: fwd_price    		- forward price of underlying
 *              	: strike             	- strike price
 *              	: vol                   - annual volatility of
 *underlying OUTPUT  : d1                    - Black-Scholes d1 : d2 -
 *Black-Scholes d2
 *
 * RETURNS      	: void
 *
 *******************************************************************************/

void srt_f_optblkd12(double fwd_price, double strike, double vol, double mat,
                     double *d1, double *d2) {
  /* ==========================================================================
   */
  /* XXX      	      				      */
  /* ==========================================================================
   */

  if (vol != 0 && mat != 0) {
    *d1 = (log(fwd_price / strike) + (vol * vol / 2) * mat) / (vol * sqrt(mat));
    *d2 = *d1 - (vol * sqrt(mat));
  } else {
    /* ------------------------------------------------------------------ */
    /* XXX  	      */
    /* ------------------------------------------------------------------ */

    if ((fwd_price / strike) > 1) {
      *d1 = INFINITY; /* Equivalent of +oo */
      *d2 = INFINITY; /* Equivalent of +oo */
    } else {
      if ((fwd_price / strike) == 1) {
        /* -------------------------------------------------- */
        /* XXX  	      */
        /* -------------------------------------------------- */

        *d1 = 0;
        *d2 = 0;
      } else {
        /* -------------------------------------------------- */
        /* XXX  	      */
        /* -------------------------------------------------- */

        *d1 = -INFINITY; /* Equivalent of +oo */
        *d2 = -INFINITY; /* Equivalent of +oo */
      }
    }
  }

} /* srt_f_optblkd12() */

/******************************************************************************/
