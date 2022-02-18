/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optratcht()
 *
 * PURPOSE      	: price of a ratchet option
 * 		  this is an option where strike is revised on initial period
 *		  when the spot hits pre-set barriers levels
 *
 * DESCRIPTION  	: this option can be decomposed in a sum of extinguish and
 *		  lightable options with barrier and stirke being
 *		  equal to each rung defined  and the associated strike
 * CALLS		: srt_f_optblksch()
 *		: srt_f_optprtext()
 *		: srt_f_optprtlgt()
 *
 * PARAMETERS   	: fwd1     	- forward price of underlying at T1
 *              	: fwd2        	- foward price of underlying at T2
 *              	: spot         	- spot price
 *              	: hist_ext     	- historical extreme
 *              	: vol1          - volatility of underlying upto T1
 *              	: vol2          - volatility of underlying upto T2
 *              	: mat1        	- maturity of 1st underlying
 *              	: mat2        	- maturity of 2nd underlying
 *              	: disc        	- discount function
 *              	: call_put      - call_put of option:
 *              	: rungs[]     	- rungs: levels where strike changes
 *              	: strikes[]    	- strike prices
 *		: greek   - premium or greeks
 *
 * RETURNS      	: ??          	- ??
 *
 *******************************************************************************/

double srt_f_optratcht(
    double         fwd1,
    double         fwd2,
    double         spot,
    double         hist_ext,
    double         vol1,
    double         vol2,
    double         mat1,
    double         mat2,
    double         disc,
    SrtCallPutType call_put,
    SrtGreekType   greek,
    double         rungs[],
    double         strikes[])
{
    int i;
    int num_strikes;
    int num_rungs;
    int first_rung;

    SrtBarrierType down_up;

    double premium = 0.0;
    double a       = 0.0;
    double answer;
    double shift;

    /***
            deals with historical extreme
            finds whether the option is increasing
            or decreasing, depending on the order of rungs
    ***/

    for (i = 1; strikes[i] > 0 && i <= MAX_RUNG; i++)
        ; /* set no. of strikes */
    num_strikes = i - 1;
    for (i = 2; rungs[i] > 0 && i <= MAX_RUNG; i++)
        ;              /* set no. of rungs   */
    num_rungs = i - 1; /* in fact, there is  */
                       /*	one rung less */
    if (num_rungs != num_strikes)
        return (0.0);

    if (num_rungs >= 3)
    {
        /* find if rungs are in increasing */
        /* or decreasing order */
        for (i = 3; i <= num_rungs; i++)
        {
            if ((rungs[i] - rungs[i - 1]) > 0)
                a++;
            else if ((rungs[i] - rungs[i - 1]) < 0)
                a--;
        }

        if (a == (double)(num_rungs - 2))
            down_up = SRT_UP; /* if increasing, then up barrier   */
        else if (a == (double)-(num_rungs - 2))
            down_up = SRT_DOWN; /* if decreasing, then down barrier */
        else
            return (0.0); /* if rungs are not ordered         */
    }
    else if (num_rungs == 2)
    { /* if only one rung                 */
        if ((rungs[2] - spot) > 0)
            down_up = SRT_UP; /* compare to spot                  */
        else if ((rungs[2] - spot) < 0)
            down_up = SRT_DOWN;
        else
            return (0.0);
    }

    a          = (down_up == SRT_UP) ? 1 : -1;
    first_rung = 2;

    /* look for the rung reached by the spot
       given historical extremum */
    for (i = 2; i <= num_rungs; i++)
    {
        if ((a * hist_ext) >= (a * rungs[i]))
            first_rung++;
    }

    if (mat1 <= 0)
    {
        if (first_rung > num_rungs)
            first_rung = num_rungs;

        premium = srt_f_optblksch(fwd2, strikes[first_rung], vol2, mat2, disc, call_put, PREMIUM);
    }
    else if (first_rung > num_rungs)
    {
        premium = srt_f_optblksch(fwd2, strikes[num_strikes], vol2, mat2, disc, call_put, PREMIUM);
    }
    else
    {
        premium = srt_f_optprtext(
            fwd2,
            fwd1,
            spot,
            strikes[first_rung - 1],
            rungs[first_rung],
            vol1,
            vol2,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);

        for (i = first_rung; i <= num_rungs; i++)
        {
            premium += srt_f_optprtlgt(
                fwd2,
                fwd1,
                spot,
                strikes[i],
                rungs[i],
                vol1,
                vol2,
                mat1,
                mat2,
                disc,
                call_put,
                down_up,
                PREMIUM);
            if (i <= num_rungs - 1)
                premium -= srt_f_optprtlgt(
                    fwd2,
                    fwd1,
                    spot,
                    strikes[i],
                    rungs[i + 1],
                    vol1,
                    vol2,
                    mat1,
                    mat2,
                    disc,
                    call_put,
                    down_up,
                    PREMIUM);
        }
    }

    switch (greek)
    {
    case PREMIUM: /*** PREMIUM ***/
        answer = premium;
        break;

    case DELTA_FWD1: /*** DELTA FWD 2***/
        shift  = fwd1 / 10000;
        answer = (srt_f_optratcht(
                      fwd1 + shift,
                      fwd2,
                      spot,
                      hist_ext,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM,
                      rungs,
                      strikes) -
                  premium) /
                 shift;
        break;

    case DELTA_FWD2: /*** DELTA FWD 2***/
        shift  = fwd2 / 10000;
        answer = (srt_f_optratcht(
                      fwd1,
                      fwd2 + shift,
                      spot,
                      hist_ext,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM,
                      rungs,
                      strikes) -
                  premium) /
                 shift;
        break;

    case DELTA: /*** DELTA SPOT + FWD1 + FWD 2 ***/
        shift  = spot / 10000;
        answer = (srt_f_optratcht(
                      fwd1 * (1 + shift / spot),
                      fwd2 * (1 + shift / spot),
                      spot + shift,
                      hist_ext,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM,
                      rungs,
                      strikes) -
                  premium) /
                 shift;
        break;

    case GAMMA: /*** GAMMA ***/
        shift  = spot / 10000;
        answer = srt_f_optratcht(
            fwd1 * (1 + shift / spot),
            fwd2 * (1 + shift / spot),
            spot + shift,
            hist_ext,
            vol1,
            vol2,
            mat1,
            mat2,
            disc,
            call_put,
            PREMIUM,
            rungs,
            strikes);
        answer += srt_f_optratcht(
            fwd1 * (1 - shift / spot),
            fwd2 * (1 - shift / spot),
            spot + shift,
            hist_ext,
            vol1,
            vol2,
            mat1,
            mat2,
            disc,
            call_put,
            PREMIUM,
            rungs,
            strikes);
        answer -= 2 * premium;
        answer /= shift * shift;
        break;

    case VEGA1: /*** VEGA1 ***/
        shift  = GVOPT.vol_add;
        answer = (srt_f_optratcht(
                      fwd1,
                      fwd2,
                      spot,
                      hist_ext,
                      vol1 + shift,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM,
                      rungs,
                      strikes) -
                  premium) /
                 shift;
        break;

    case VEGA2: /*** VEGA2 ***/
        shift  = GVOPT.vol_add;
        answer = (srt_f_optratcht(
                      fwd1,
                      fwd2,
                      spot,
                      hist_ext,
                      vol1,
                      vol2 + shift,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM,
                      rungs,
                      strikes) -
                  premium) /
                 shift;
        break;

    case THETA: /*** THETA  ***/
        shift  = YEARS_IN_DAY;
        answer = srt_f_optratcht(
                     fwd1,
                     fwd2,
                     spot,
                     hist_ext,
                     vol1,
                     vol2,
                     mat1 - shift,
                     mat2 - shift,
                     disc * exp(-shift * log(disc) / mat2),
                     call_put,
                     PREMIUM,
                     rungs,
                     strikes) -
                 premium;
        break;

    default:
        answer = UNKNOWN_GREEK;
        break;
    }

    return (answer);
} /* END srt_f_optratcht() */

/******************************************************************************/
