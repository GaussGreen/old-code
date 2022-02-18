/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optpaloba()
 *
 * PURPOSE     	: partial lookback options for stocks, bonds & swaps
 *
 * DESCRIPTION  	: XX
 *
 * CALLS		: norm() in gen_math.c
 *               : srt_f_optlokbak
 *               : srt_f_optblksch
 *
 * PARAMETERS  	: fwd_price1 	- forward price1 of underlying
 *             	: fwd_price2 	- forward price2 of underlying
 *             	: spot     	- spot price of underlying
 *             	: hist_extr  	- historical extrema for hedging purposes
 *             	: vol1         	- volatility of underlying for subperiod
 *             	: vol2         	- volatility of underlying for global period
 *             	: mat1        	- subperiod duration
 *             	: mat2        	- maturity of the option
 *             	: disc       	- discount factor at maturity
 *             	: call_put      - type of option: 0 call, 1 put
 *             	: greek   - premium (0 ONLY !!!)
 *
 * RETURNS      	: premium    	- price of the option
 *
 *******************************************************************************/

double srt_f_optpaloba(
    double         fwd_price1,
    double         fwd_price2,
    double         spot,
    double         hist_extr,
    double         vol1,
    double         vol2,
    double         mat1,
    double         mat2,
    double         disc,
    SrtCallPutType call_put,
    SrtGreekType   greek)

{
    /* BEGIN srt_f_optpaloba */

    double l_ext;               /** = log(MIN..MAX /SPOT)	**/
    double mu01, mu12, mu02;    /** mu = log(fwd/spot)/mat  	**/
    double var01, var12, var02; /** = vol*vol*mat		**/
    double corr12;
    double d_ext01, d_ext02, d_12; /** used in B&S call_put formula**/
    double nu;                     /** = 2*mu01/(vol1*vol1)	**/
    double vol;

    double premium1, premium2, premium3, premium4, premium5;
    double premium;
    double answer, shift;

    /** 	   	call_put 	0=CALL	 	1=PUT
            -----------------------------------------
            sub=[0,T1]		S-min		Max-S

    **/

    /**IN ORDER TO PREVENT ERRORS IN EXTREMA **/
    if ((call_put == SRT_CALL) && (hist_extr > spot))
        hist_extr = spot;
    if ((call_put == SRT_PUT) && (hist_extr < spot))
        hist_extr = spot;

    /** IN ORDER TO PREVENT ERRORS IN CALCULATION DUE TO BIVAR OR NORM **/
    if (fabs(fwd_price1 - spot) <= 0.00001)
        fwd_price1 = spot;

    l_ext = log(hist_extr / spot);

    var01 = vol1 * vol1 * mat1;
    var02 = vol2 * vol2 * mat2;
    var12 = var02 - var01;

    /** VARIABLES DEFINED IF VOLATILTIES ARE != 0 **/
    if ((vol1 != 0) && (vol2 != 0))
    {
        corr12 = sqrt(var01 / var02);
        mu01   = log(fwd_price1 / spot) / mat1;
        if (fabs(mu01) < 0.000001)
            mu01 = 0.0;
        d_ext01 = (-l_ext + mu01 * mat1 + 0.5 * var01) / sqrt(var01);
        mu02    = log(fwd_price2 / spot) / mat2;
        d_ext02 = (-l_ext + mu02 * mat2 + 0.5 * var02) / sqrt(var02);
        mu12    = log(fwd_price2 / fwd_price1) / (mat2 - mat1);
        d_12    = (mu12 * (mat2 - mat1) + 0.5 * var12) / sqrt(var12);
        nu      = 2 * mu01 / (vol1 * vol1);
    };

    /** IN ORDER TO ALLOW PERIOD=SUB_PERIOD   **/
    if (mat2 == mat1)
    {
        if ((fwd_price1 == fwd_price2) && (vol1 == vol2))
        {
            premium =
                srt_f_optlokbak(fwd_price1, spot, hist_extr, vol1, mat1, disc, call_put, PREMIUM);
        }
        else
            premium = 0.0;
    }
    else
        /** REAL BEGINNING OF CALCULATIONS **/
        if (call_put == SRT_CALL)
    {
        /** CALL  S - min S(s)  **/
        if (mat1 == 0.0)
        {
            premium = srt_f_optblksch(fwd_price2, hist_extr, vol2, mat2, disc, call_put, PREMIUM);
        }
        else if (vol1 == 0.0)
        {
            if (fwd_price1 < hist_extr)
                hist_extr = fwd_price1;

            if (vol2 != 0.0)
            {
                vol     = sqrt(var12 / (mat2 - mat1));
                premium = srt_f_optblksch(
                    fwd_price2, hist_extr, vol, mat2 - mat1, disc, call_put, PREMIUM);
            }
            else
            {
                if (fwd_price2 > hist_extr)
                {
                    premium = fwd_price2 - hist_extr;
                    premium *= disc;
                }
                else
                    premium = 0.0;
            }
        }
        else
        {
            premium1 = fwd_price2 * bivar(d_ext01, d_ext02, corr12);
            premium1 -= hist_extr * bivar(d_ext01 - sqrt(var01), d_ext02 - sqrt(var02), corr12);

            if (mu01 != 0.0)
            {
                premium2 = fwd_price2 * norm(d_12);
                premium2 -= (1 + 1 / nu) * fwd_price1 * norm(d_12 - sqrt(var12));
                premium2 *= norm(-d_ext01);

                premium3 = bivar(
                    d_ext01 - sqrt(var01) + 2 * l_ext / sqrt(var01),
                    d_ext02 - sqrt(var02) + 2 * l_ext / sqrt(var02),
                    corr12);
                premium3 *= 1 / nu * spot * exp(l_ext * nu);

                premium4 = bivar(
                    -(d_12 - sqrt(var12)) + nu * sqrt(var12),
                    d_ext01 * sqrt(var01 / var02) - (d_12 - sqrt(var12)) * sqrt(var12 / var02) +
                        nu * var12 / sqrt(var02),
                    sqrt(var12 / var02));
                premium4 -= norm(-(d_12 - sqrt(var12)) + nu * sqrt(var12));
                premium4 *= 1 / nu * fwd_price1 *
                            exp(nu * (-mu12 * (mat2 - mat1) + 0.5 * var12 * (nu + 1)));

                premium5 = 0.0;
            }
            else
            {
                premium2 = fwd_price2 * norm(d_12);
                premium2 -= fwd_price1 * sqrt(var12) * gauss(d_12 - sqrt(var12));
                premium2 *= norm(-d_ext01);

                premium3 = bivar(-d_ext01, d_ext02 - sqrt(var02) + 2 * l_ext / sqrt(var02), corr12);
                premium3 *= -spot - fwd_price1 * (d_ext01 * sqrt(var01));

                premium4 = norm(-(d_12 - sqrt(var12)));
                premium4 -= bivar(
                    -(d_12 - sqrt(var12)),
                    d_ext01 * sqrt(var01 / var02) - (d_12 - sqrt(var12)) * sqrt(var12 / var02),
                    sqrt(var12 / var02));
                premium4 *= fwd_price1 * (1 + d_12 * sqrt(var12) - var12);

                premium5 = sqrt(var01) * gauss(-d_ext01) * norm(d_12 - sqrt(var12));
                premium5 += sqrt(var02) *
                            gauss(
                                (-d_ext01 * sqrt(var01) + (d_12 - sqrt(var12)) * sqrt(var12)) /
                                sqrt(var02)) *
                            norm(
                                ((-d_12 + sqrt(var12)) * sqrt(var01) - d_ext01 * sqrt(var12)) /
                                sqrt(var02));
                premium5 *= fwd_price1;
            }

            premium = disc * (premium1 + premium2 + premium3 + premium4 + premium5);
        }
    }
    else if (call_put == SRT_PUT)
    {
        /** PUT  Max S(s) - S  **/
        if (mat1 == 0.0)
        {
            premium = srt_f_optblksch(fwd_price2, hist_extr, vol2, mat2, disc, call_put, PREMIUM);
        }
        else if (vol1 == 0.0)
        {
            if (fwd_price1 > hist_extr)
                hist_extr = fwd_price1;

            if (vol2 != 0.0)
            {
                vol     = sqrt(var12 / (mat2 - mat1));
                premium = srt_f_optblksch(
                    fwd_price2, hist_extr, vol, mat2 - mat1, disc, call_put, PREMIUM);
            }
            else
            {
                if (fwd_price2 < hist_extr)
                {
                    premium = hist_extr - fwd_price2;
                    premium *= disc;
                }
                else
                    premium = 0.0;
            }
        }
        else
        {
            premium1 = hist_extr * bivar(-d_ext01 + sqrt(var01), -d_ext02 + sqrt(var02), corr12);
            premium1 -= fwd_price2 * bivar(-d_ext01, -d_ext02, corr12);

            if (mu01 != 0.0)
            {
                premium2 = (1 + 1 / nu) * fwd_price1 * norm(-d_12 + sqrt(var12));
                premium2 -= fwd_price2 * norm(-d_12);
                premium2 *= norm(d_ext01);

                premium3 = -bivar(
                    -d_ext01 + sqrt(var01) - 2 * l_ext / sqrt(var01),
                    -d_ext02 + sqrt(var02) - 2 * l_ext / sqrt(var02),
                    corr12);
                premium3 *= 1 / nu * spot * exp(l_ext * nu);

                premium4 = -bivar(
                    (d_12 - sqrt(var12)) - nu * sqrt(var12),
                    -d_ext01 * sqrt(var01 / var02) + (d_12 - sqrt(var12)) * sqrt(var12 / var02) -
                        nu * var12 / sqrt(var02),
                    sqrt(var12 / var02));
                premium4 += norm(d_12 - sqrt(var12) - nu * sqrt(var12));
                premium4 *= 1 / nu * fwd_price1 *
                            exp(nu * (-mu12 * (mat2 - mat1) + 0.5 * var12 * (nu + 1)));

                premium5 = 0.0;
            }
            else
            {
                premium2 = -fwd_price1 * sqrt(var12) * gauss(-d_12 + sqrt(var12));
                premium2 -= fwd_price2 * norm(-d_12);
                premium2 *= norm(d_ext01);

                premium3 = bivar(d_ext01, -d_ext02 + sqrt(var02) - 2 * l_ext / sqrt(var02), corr12);
                premium3 *= spot + fwd_price1 * (d_ext01 * sqrt(var01));

                premium4 = bivar(
                    d_12 - sqrt(var12),
                    -d_ext01 * sqrt(var01 / var02) + (d_12 - sqrt(var12)) * sqrt(var12 / var02),
                    sqrt(var12 / var02));
                premium4 -= norm(d_12 - sqrt(var12));
                premium4 *= fwd_price1 * (1 + d_12 * sqrt(var12) - var12);

                premium5 = sqrt(var01) * gauss(d_ext01) * norm(-d_12 + sqrt(var12));
                premium5 +=
                    sqrt(var02) *
                    gauss(
                        (d_ext01 * sqrt(var01) - (d_12 - sqrt(var12)) * sqrt(var12)) /
                        sqrt(var02)) *
                    norm(
                        ((d_12 - sqrt(var12)) * sqrt(var01) + d_ext01 * sqrt(var12)) / sqrt(var02));
                premium5 *= fwd_price1;
            }

            premium = disc * (premium1 + premium2 + premium3 + premium4 + premium5);
        }
    }

    if (greek == PREMIUM)
    {
        answer = premium;
    }
    else if (greek == DELTA_FWD1)
    { /* DELTA FWD1 UP */
        shift  = fwd_price1 / 10000;
        answer = (srt_f_optpaloba(
                      fwd_price1 + shift,
                      fwd_price2,
                      spot,
                      hist_extr,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM) -
                  premium) /
                 shift;
    }
    else if (greek == DELTA_FWD2)
    { /* DELTA FWD2 UP */
        shift  = fwd_price2 / 10000;
        answer = (srt_f_optpaloba(
                      fwd_price1,
                      fwd_price2 + shift,
                      spot,
                      hist_extr,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM) -
                  premium) /
                 shift;
    }
    else if (greek == DELTA)
    { /* DELTA */
        shift  = spot / 10000;
        answer = (srt_f_optpaloba(
                      fwd_price1 + shift,
                      fwd_price2 + shift,
                      spot + shift,
                      hist_extr,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM) -
                  premium) /
                 shift;
    }
    /*
            else if ( greek == GAMMA_FWD1 )
            {
                    shift = fwd_price1 / 10000;
                    answer = ( srt_f_optpaloba(	fwd_price1+shift,
                                                                            fwd_price2,
                                                                            spot,
                                                                            hist_extr,
                                                                            vol1,
                                                                            vol2,
                                                                            mat1,
                                                                            mat2,
                                                                            disc,
                                                                            call_put,
                                                                            PREMIUM)
                                    + srt_f_optpaloba(	fwd_price1-shift,
                                                                            fwd_price2,
                                                                            spot,
                                                                            hist_extr,
                                                                            vol1,
                                                                            vol2,
                                                                            mat1,
                                                                            mat2,
                                                                            disc,
                                                                            call_put,
                                                                            PREMIUM)
                                    - 2*premium) / (shift*shift);
            }
            else if ( greek == GAMMA_FWD2 )
            {
                    shift = fwd_price2 / 10000;
                    answer = ( srt_f_optpaloba(	fwd_price1,
                                                                            fwd_price2+shift,
                                                                            spot,
                                                                            hist_extr,
                                                                            vol1,
                                                                            vol2,
                                                                            mat1,
                                                                            mat2,
                                                                            disc,
                                                                            call_put,
                                                                            PREMIUM)
                                    + srt_f_optpaloba(	fwd_price1,
                                                                            fwd_price2-shift,
                                                                            spot,
                                                                            hist_extr,
                                                                            vol1,
                                                                            vol2,
                                                                            mat1,
                                                                            mat2,
                                                                            disc,
                                                                            call_put,
                                                                            PREMIUM)
                                    - 2*premium) / (shift*shift);
            }
    */
    else if (greek == GAMMA)
    { /*   GAMMA SPOT  */
        shift  = spot / 10000;
        answer = (srt_f_optpaloba(
                      fwd_price1 * (1 + shift / spot),
                      fwd_price2 * (1 + shift / spot),
                      spot + shift,
                      hist_extr,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM) +
                  srt_f_optpaloba(
                      fwd_price1 * (1 - shift / spot),
                      fwd_price2 * (1 - shift / spot),
                      spot - shift,
                      hist_extr,
                      vol1,
                      vol2,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      PREMIUM) -
                  2 * premium) /
                 (shift * shift);
    }
    else if (greek == VEGA1)
    { /*   VEGA 1   */
        shift = GVOPT.vol_add;
        answer =
            (srt_f_optpaloba(
                 fwd_price1,
                 fwd_price2,
                 spot,
                 hist_extr,
                 vol1 + shift,
                 vol2,
                 mat1,
                 mat2,
                 disc,
                 call_put,
                 PREMIUM) -
             premium);
    }
    else if (greek == VEGA2)
    { /*   VEGA 2   */
        shift = GVOPT.vol_add;
        answer =
            (srt_f_optpaloba(
                 fwd_price1,
                 fwd_price2,
                 spot,
                 hist_extr,
                 vol1,
                 vol2 + shift,
                 mat1,
                 mat2,
                 disc,
                 call_put,
                 PREMIUM) -
             premium);
    }
    else if (greek == THETA)
    { /*   THETA   */
        shift = YEARS_IN_DAY;
        answer =
            (srt_f_optpaloba(
                 fwd_price1,
                 fwd_price2,
                 spot,
                 hist_extr,
                 vol1,
                 vol2,
                 mat1 - shift,
                 mat2 - shift,
                 disc * exp(-shift * log(disc) / mat2),
                 call_put,
                 PREMIUM) -
             premium);
    }
    else
        answer = UNKNOWN_GREEK;

    return (answer);

} /* END srt_f_optpaloba() */

/******************************************************************************/
