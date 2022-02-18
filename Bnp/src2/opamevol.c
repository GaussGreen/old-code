/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optamevol(...)
 *
 * PURPOSE      	: calculates volatility of american option implied by  premium
 *
 * DESCRIPTION  	: iteration
 *
 * CALLS		: srt_f_optamevol(...)
 *
 * PARAMETERS   	: premium   	- option premium
 * 		: spot     	- spot price
 *              	: forward     	- forward price of underlying
 *              	: strike   	- strike price
 *              	: disc_fact    	- discount factor
 *              	: opt_mat 	- maturity of option
 *              	: call_put 	- type of option: 0 call, 1 put
 *              	: step_num 	- number of steps for tree calculation
 *
 * RETURNS      	: volatility
 *
 *******************************************************************************/

static double         RTSAFE_PREMIUM;
static double         RTSAFE_SPOT;
static double         RTSAFE_FORWARD;
static double         RTSAFE_STRIKE;
static double         RTSAFE_DISCFACT;
static double         RTSAFE_OPTMAT;
static SrtCallPutType RTSAFE_CALLPUT;
static int            RTSAFE_STEPNUM;

/*--------------------------------------------------------------------*/

static Err rtsafe_opamevol(double guess_vol, double* diff_prem, double* deriv)
{
    *diff_prem = srt_f_optameopt(
                     RTSAFE_SPOT,
                     RTSAFE_FORWARD,
                     RTSAFE_STRIKE,
                     RTSAFE_DISCFACT,
                     guess_vol,
                     RTSAFE_OPTMAT,
                     RTSAFE_CALLPUT,
                     RTSAFE_STEPNUM,
                     PREMIUM) -
                 RTSAFE_PREMIUM;

    *deriv = srt_f_optameopt(
        RTSAFE_SPOT,
        RTSAFE_FORWARD,
        RTSAFE_STRIKE,
        RTSAFE_DISCFACT,
        guess_vol,
        RTSAFE_OPTMAT,
        RTSAFE_CALLPUT,
        RTSAFE_STEPNUM,
        VEGA);
    return NULL;
}

/*--------------------------------------------------------------------*/

Err srt_f_optamevol(
    double         premium,
    double         spot,
    double         forward,
    double         strike,
    double         disc_fact,
    double         opt_mat,
    SrtCallPutType call_put,
    int            step_num,
    double*        implied_vol)
{
    double intrinsic;
    SrtErr err = NULL;

    /* Compute option intrinsic value */
    intrinsic = srt_f_optameopt(
        spot, forward, strike, disc_fact, NULL_VOL, opt_mat, call_put, step_num, PREMIUM);

    /* Checks target premium is above intrinsic value */
    if (intrinsic > premium)
    {
        *implied_vol = 0.0;
        return serror("Intrinsic higher than option premium");
    }

    /* Renormalise the option price by the spot value */
    if (spot != 0.00)
    {
        forward /= spot;
        strike /= spot;
        premium /= spot;
        spot = 1.0;
    }

    RTSAFE_SPOT     = spot;
    RTSAFE_FORWARD  = forward;
    RTSAFE_STRIKE   = strike;
    RTSAFE_DISCFACT = disc_fact;
    RTSAFE_OPTMAT   = opt_mat;
    RTSAFE_CALLPUT  = call_put;
    RTSAFE_PREMIUM  = premium;
    RTSAFE_STEPNUM  = step_num;

    /* Start Newton like iterations */
    err = rtsafe(
        rtsafe_opamevol,
        1e-11,     /* inf of the solution */
        2.0,       /* max of the solution */
        0.0000001, /* accuracy */
        100,       /* number of iterations */
        implied_vol);

    return err;

} /* END srt_f_optamevol(...) */

/******************************************************************************/
