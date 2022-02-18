/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_opteurdig( 7)
 *
 * PURPOSE      	: European Digital, with payment at maturity
 *
 * DESCRIPTION  	: ??
 *
 * PARAMETERS      ???
 *	INPUT	: fwd_price	- forward underlying price
 *              	: strike      	- strike price
 *              	: vol         	- annual volatility
 *              	: mat         	- maturity, in years
 *              	: disc        	- discount factor to expiry
 *              	: type        	- type of option: 0 call, 1 put
 *	OUTPUT	: greeks        - greeks for given option pricing
 *
 * RETURNS      	: premium       - option premium
 *
 *******************************************************************************/

double srt_f_opteurdig(
    double         fwd,
    double         barrier,
    double         vol,
    double         mat,
    double         disc,
    SrtBarrierType below_above,
    SrtGreekType   greek)
{
    double vol_squared;  /* is (vol * vol)     	      */
    double vol_root_mat; /* is (vol * sqrt( mat)       */

    double d2; /* d2 as Black-Scholes	      */

    double premium;
    double answer;
    double shift;

    if (mat < 0)
    {
        premium = 0.0;
    }
    else if ((vol == 0) || (mat == 0))
    {
        if (below_above == SRT_UP)
        {
            if (fwd >= barrier)
                premium = disc;
            else
                premium = 0.0;
        }
        else if (below_above == SRT_DOWN)
        {
            if (fwd > barrier)
                premium = 0.0;
            else
                premium = disc;
        }
    }
    else if (barrier != 0)
    {
        /* ------------------------------------------------------------------
           Expected  case: no divide-by-zeros
           ------------------------------------------------------------------ */

        vol_squared  = (vol * vol);
        vol_root_mat = vol * sqrt(mat);

        d2 = (log(fwd / barrier) - ((vol_squared / 2) * mat)) / vol_root_mat;

        if (below_above == SRT_DOWN)
        {
            /* ----------------------------------------------------------
               BELOW: probability of ending below the barrier
               ---------------------------------------------------------- */

            premium = disc * norm(-d2);
        }
        else if (below_above == SRT_UP)
        {
            /* ----------------------------------------------------------
               ABOVE:probability of ending above the barrier
               ---------------------------------------------------------- */

            premium = disc * norm(d2);
        }
    }
    else /* barrier = 0 */
    {
        if (below_above == SRT_DOWN)
        {
            /* ----------------------------------------------------------
               BELOW: probability of ending below the barrier
               ---------------------------------------------------------- */

            premium = 0;
        }
        else if (below_above == SRT_UP)
        {
            /* ----------------------------------------------------------
               ABOVE:probability of ending above the barrier
               ---------------------------------------------------------- */

            premium = disc;
        }
    }

    switch (greek)
    {
    case PREMIUM: /*** PREMIUM ***/
        answer = premium;
        break;

    case DELTA_FWD: /*** DELTA FWD ***/
        shift  = fwd / 10000;
        answer = (srt_f_opteurdig(fwd + shift, barrier, vol, mat, disc, below_above, PREMIUM) -
                  premium) /
                 shift;
        break;

    case GAMMA: /*** GAMMA ***/
        shift  = fwd / 10000;
        answer = srt_f_opteurdig(fwd + shift, barrier, vol, mat, disc, below_above, PREMIUM);
        answer += srt_f_opteurdig(fwd - shift, barrier, vol, mat, disc, below_above, PREMIUM);
        answer -= 2 * premium;
        answer /= shift * shift;
        break;

    case VEGA: /*** VEGA ***/
        shift  = GVOPT.vol_add;
        answer = (srt_f_opteurdig(fwd, barrier, vol + shift, mat, disc, below_above, PREMIUM) -
                  premium) /
                 shift;
        break;

    case THETA: /*** THETA  ***/
        shift  = YEARS_IN_DAY;
        answer = srt_f_opteurdig(
                     fwd,
                     barrier,
                     vol,
                     mat - shift,
                     disc * exp(-shift * log(disc) / mat),
                     below_above,
                     PREMIUM) -
                 premium;
        break;

    default:
        answer = UNKNOWN_GREEK;
        break;
    }

    return (answer);

} /* END srt_f_opteurdig() */
