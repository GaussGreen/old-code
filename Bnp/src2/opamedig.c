/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optamedig( 7)
 *
 * PURPOSE      	: American Digital, with payment at maturity
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

double srt_f_optamedig(
    double         fwd,
    double         spot,
    double         barrier,
    double         vol,
    double         mat,
    double         disc,
    SrtBarrierType below_above,
    SrtMinmaxType  min_max,
    SrtGreekType   greek)
{
    double rate;  /* r		  			      */
    double drift; /* mu 					      */
    double x_bar; /* X bar: normalised barrier  */

    double vol_squared;  /* is (vol * vol)     	      */
    double vol_root_mat; /* is (vol * sqrt( mat)       */
    double factor;

    double first_term;
    double second_term;

    double prob_floor_above = 0;
    double prob_floor_below = 0;
    double prob_cap_above   = 0;
    double prob_cap_below   = 0;

    double probability = 0.0;
    double premium     = 0.0;
    double answer      = 0.0;
    double shift       = 0.0;

    /* ==========================================================================
       Lets do it!
       ========================================================================== */

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
            else if (fwd < spot)
                premium = 0.0;
        }
        else if (below_above == SRT_DOWN)
        {
            if (fwd <= barrier)
                premium = disc;
            else if (fwd > spot)
                premium = 0.0;
        }
    }
    else if ((spot != 0) && (barrier != 0))
    {
        /* ------------------------------------------------------------------
           Expected  case: no divide-by-zeros
           ------------------------------------------------------------------ */
        vol_squared = (vol * vol);

        rate  = log(fwd / spot) / mat;
        drift = rate - (vol_squared / 2);
        x_bar = log(barrier / spot);

        vol_root_mat = vol * sqrt(mat);
        factor       = exp((2 * drift * x_bar) / vol_squared);

        if (min_max == SRT_MIN)
        {
            /* ----------------------------------------------------------
               Digital Floor: boundary condition B <= S
                                                       0
               ---------------------------------------------------------- */

            if ((below_above == SRT_DOWN) && (spot <= barrier))
                probability = 1.0;
            else if ((below_above == SRT_UP) && (spot <= barrier))
                probability = 0.0;
            else
            {
                first_term = norm((-x_bar + (drift * mat)) / vol_root_mat);

                second_term = factor * norm((x_bar + (drift * mat)) / vol_root_mat);

                prob_floor_above = first_term - second_term;
                prob_floor_below = 1 - prob_floor_above;

                if (below_above == SRT_DOWN)
                {
                    /* ------------------------------------------
                       Prob[ Min S  <= Barrier ]
                                  T
                       ------------------------------------------ */
                    probability = prob_floor_below;
                }
                else if (below_above == SRT_UP)
                {
                    /* ------------------------------------------
                       Prob[ Min S  >= Barrier ]
                                  T
                    ------------------------------------------ */

                    probability = prob_floor_above;
                }
            }
        }
        else if (min_max == SRT_MAX)
        {
            /* ----------------------------------------------------------
               Digital Cap: boundary condition is B >= S
                                                        0
               ---------------------------------------------------------- */

            if ((below_above == SRT_UP) && (spot >= barrier))
                probability = 1.0;
            else if ((below_above == SRT_DOWN) && (spot >= barrier))
                probability = 0.0;
            else
            {
                first_term = norm((x_bar - (drift * mat)) / vol_root_mat);

                second_term = factor * norm((-x_bar - (drift * mat)) / vol_root_mat);

                prob_cap_below = first_term - second_term;
                prob_cap_above = 1 - prob_cap_below;

                if (below_above == SRT_DOWN)
                {
                    /* ------------------------------------------
                       Prob[ Max S  <= Barrier ]
                                  T
                       ------------------------------------------ */

                    probability = prob_cap_below;
                }
                else if (below_above == SRT_UP)
                {
                    /* ------------------------------------------
                       Prob[ Max S  >= Barrier ]
                                  T
                       ------------------------------------------ */
                    probability = prob_cap_above;
                }
            }
        }

        premium = probability * disc;
    }
    else if (spot != 0)
    {
        if (below_above == SRT_UP)
            probability = 1.0;
        else if (below_above == SRT_DOWN)
            probability = 0.0;

        premium = probability * disc;
    }

    switch (greek)
    {
    case PREMIUM: /*** PREMIUM ***/
        answer = premium;
        break;

    case DELTA_FWD: /*** DELTA FWD ***/
        shift  = fwd / 10000;
        answer = (srt_f_optamedig(
                      fwd + shift, spot, barrier, vol, mat, disc, below_above, min_max, PREMIUM) -
                  premium) /
                 shift;
        break;

    case DELTA: /*** DELTA SPOT + FWD ***/
        shift  = spot / 10000;
        answer = (srt_f_optamedig(
                      fwd * (1 + shift / spot),
                      spot + shift,
                      barrier,
                      vol,
                      mat,
                      disc,
                      below_above,
                      min_max,
                      PREMIUM) -
                  premium) /
                 shift;
        break;

    case GAMMA: /*** GAMMA ***/
        shift  = spot / 10000;
        answer = srt_f_optamedig(
            fwd * (1 + shift / spot),
            spot + shift,
            barrier,
            vol,
            mat,
            disc,
            below_above,
            min_max,
            PREMIUM);
        answer += srt_f_optamedig(
            fwd * (1 - shift / spot),
            spot - shift,
            barrier,
            vol,
            mat,
            disc,
            below_above,
            min_max,
            PREMIUM);
        answer -= 2 * premium;
        answer /= shift * shift;
        break;

    case VEGA: /*** VEGA ***/
        shift  = GVOPT.vol_add;
        answer = (srt_f_optamedig(
                      fwd, spot, barrier, vol + shift, mat, disc, below_above, min_max, PREMIUM) -
                  premium) /
                 shift;
        break;

    case THETA: /*** THETA  ***/
        shift  = YEARS_IN_DAY;
        answer = srt_f_optamedig(
                     fwd,
                     spot,
                     barrier,
                     vol,
                     mat - shift,
                     disc * exp(-shift * log(disc) / mat),
                     below_above,
                     min_max,
                     PREMIUM) -
                 premium;
        break;

    default:
        answer = UNKNOWN_GREEK;
        break;
    }

    return (answer);

} /*  END srt_f_optamedig() */
