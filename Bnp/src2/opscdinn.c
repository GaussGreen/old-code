/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optscdinn( 14 )
 *
 * PURPOSE      	: SCUD
 *
 * DESCRIPTION  	: Barrier option on secondary underlying
 *
 * CALLS		: srt_f_optblksch()
 *		: scud_out_function()
 *		: SELF
 *
 * PARAMETERS   	: fwdx          - forward price of 1st underlying x
 *              	: fwdy          - forward price of 2nd underlying y
 *              	: spoty         - spot price of 2nd underlying y
 *              	: strike        - strike price
 *              	: barrier       - barrier level (on 2nd underlying y)
 *              	: sigx          - ?? of 1st underlying x
 *              	: sigy          - ?? of 2nd underlying y
 *              	: rho           - ??
 *              	: mat           - maturity of option, in years
 *              	: disc          - discount factor
 *              	: call_put    	- type of option: 0 call, 1 put
 *              	: down_up     - type of scud: ??
 *              	: greek         	- ??
 *
 * RETURNS      	: ??            - ??
 *
 *******************************************************************************/

double srt_f_optscdinn(
    double         fwdx,
    double         fwdy,
    double         spoty,
    double         strike,
    double         barrier,
    double         sigx,
    double         sigy,
    double         rho,
    double         mat1,
    double         mat2,
    double         disc,
    SrtCallPutType call_put,
    SrtBarrierType down_up,
    SrtGreekType   greek)
{
    double premium;
    double result;

    double shift;
    double shiftx;
    double shifty;

    premium = srt_f_optblksch(fwdx, strike, sigx, mat2, disc, call_put, PREMIUM);
    premium -= srt_f_optscdout(
        fwdx,
        fwdy,
        spoty,
        strike,
        barrier,
        sigx,
        sigy,
        rho,
        mat1,
        mat2,
        disc,
        call_put,
        down_up,
        PREMIUM);

    switch (greek)
    {
    case PREMIUM:
        return (premium);
        break;

    case DELTAX:
        shift  = fwdx / 10000;
        result = (srt_f_optscdinn(
                      fwdx + shift,
                      fwdy,
                      spoty,
                      strike,
                      barrier,
                      sigx,
                      sigy,
                      rho,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      down_up,
                      PREMIUM) -
                  premium) /
                 shift;
        return (result);
        break;

    case DELTAY:
        shift  = fwdy / 10000;
        result = (srt_f_optscdinn(
                      fwdx,
                      fwdy + shift,
                      spoty + shift,
                      strike,
                      barrier,
                      sigx,
                      sigy,
                      rho,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      down_up,
                      PREMIUM) -
                  premium) /
                 shift;
        return (result);
        break;

    case GAMMAX:
        shift  = fwdx / 1000;
        result = srt_f_optscdinn(
            fwdx + shift,
            fwdy,
            spoty,
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result += srt_f_optscdinn(
            fwdx - shift,
            fwdy,
            spoty,
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result -= 2 * premium;
        result /= shift * shift;
        return (result);
        break;

    case GAMMAY:
        shift  = fwdy / 1000;
        result = srt_f_optscdinn(
            fwdx,
            fwdy + shift,
            spoty * (1 + shift / fwdy),
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result += srt_f_optscdinn(
            fwdx,
            fwdy - shift,
            spoty * (1 - shift / fwdy),
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result -= 2 * premium;
        result /= shift * shift;
        return (result);
        break;

    case GAMMAXY:
        shifty = fwdy / 1000;
        shiftx = fwdx / 1000;
        result = srt_f_optscdinn(
            fwdx + shiftx,
            fwdy + shifty,
            spoty * (1 + shifty / fwdy),
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result += srt_f_optscdinn(
            fwdx - shiftx,
            fwdy - shifty,
            spoty * (1 - shifty / fwdy),
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result -= srt_f_optscdinn(
            fwdx + shiftx,
            fwdy - shifty,
            spoty * (1 - shifty / fwdy),
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result -= srt_f_optscdinn(
            fwdx - shiftx,
            fwdy + shifty,
            spoty * (1 + shifty / fwdy),
            strike,
            barrier,
            sigx,
            sigy,
            rho,
            mat1,
            mat2,
            disc,
            call_put,
            down_up,
            PREMIUM);
        result /= 4 * shiftx * shifty;
        return (result);
        break;

    case VEGAX:
        shift  = GVOPT.vol_add;
        result = (srt_f_optscdinn(
                      fwdx,
                      fwdy,
                      spoty,
                      strike,
                      barrier,
                      sigx + shift,
                      sigy,
                      rho,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      down_up,
                      PREMIUM) -
                  premium) /
                 shift;
        return (result);
        break;

    case VEGAY:
        shift  = GVOPT.vol_add;
        result = (srt_f_optscdinn(
                      fwdx,
                      fwdy,
                      spoty,
                      strike,
                      barrier,
                      sigx,
                      sigy + shift,
                      rho,
                      mat1,
                      mat2,
                      disc,
                      call_put,
                      down_up,
                      PREMIUM) -
                  premium) /
                 shift;
        return (result);
        break;

    case THETA:
        shift  = YEARS_IN_DAY;
        result = srt_f_optscdinn(
                     fwdx,
                     fwdy,
                     spoty,
                     strike,
                     barrier,
                     sigx,
                     sigy,
                     rho,
                     mat1 - shift,
                     mat2 - shift,
                     disc * exp(-shift * log(disc) / mat2),
                     call_put,
                     down_up,
                     PREMIUM) -
                 premium;
        return (result);
        break;

    default:
        return (UNKNOWN_GREEK);
        break;
    }

} /* END srt_f_optscdout() */

/******************************************************************************/
