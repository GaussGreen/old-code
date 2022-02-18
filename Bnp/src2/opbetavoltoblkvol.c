

/******************************************************************************/
/*******************************************************************************
*
* FUNCTION     	: srt_f_optbetavoltoblkvol(...)
*                 srt_f_optbetastochvoltoblkvol(...)
*
* PURPOSE      	: A Quick vol transformation to go from a BS vol to a BS beta vol
*
* DESCRIPTION  	: The volatility transformation is based on a small noise expansion
                  of the density for both BS model and its equivalent BS beta
                                  (see Pat Hagan's paper on the subject)
*
*
* PARAMETERS
*	INPUT	    : fwd_price	    - forward underlying price
*              	: strike      	- strike price
*              	: betavol       - annual beta volatility
*              	: mat         	- initial time, in years
*				: beta          - dF = a*F^{beta}dW
*							      da = v a dZ
*                                 <dZ,dW> = rho dt
* RETURNS      	: bsvol         - equivalent BS vol
*
*******************************************************************************/

/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/* ------------------------------------------------------------------------------ */

Err srt_f_optbetavoltoblkvol(
    double forward, double strike, double betavol, double mat, double beta, double* bsvol)
{
    double middle;
    double distance;
    double space_correction;
    double time_correction;
    double variance;
    Err    err = NULL;

    /* A few security checks */
    if ((mat <= 0.0) || (betavol == 0.0))
    {
        *bsvol = 0.0;
        return NULL;
    }

    if ((beta > 1.0) || (beta < 0.0))
        return serror("Beta has to be between 0.0 and 1.0");

    /* The point around which the expansion is made : 0.5*(forward + strike) */
    middle = 0.5 * (forward + strike);

    /* The distance between the forward and the strike : (forward - strike ) */
    distance = forward - strike;

    /* The equivalent lognormal variance of the distribution */
    variance = betavol * betavol * pow(middle, 2.0 * (beta - 1.0)) * mat;

    /* The leading order term for corection in space (diffusion term) */
    space_correction =
        (1.0 - beta) * (2.0 + beta) / 24.0 * (distance / middle) * (distance / middle);

    /* The leading order term for corection due to time */
    time_correction = (1.0 - beta) * (1.0 - beta) / 24.0 * variance;

    /* The equivalent BlackScholes BEta vol (constant elasticity) */
    *bsvol = betavol * pow(middle, beta - 1.0) * (1.0 + space_correction + time_correction);

    /* Return a success message */
    return NULL;
}

Err srt_f_optbetastochvoltoblkvol(
    double  forward,
    double  strike,
    double  betavol,
    double  vovol,
    double  rho,
    double  mat,
    double  beta,
    double* bsvol)
{
    double zz;
    double xx;
    double first_term;
    double second_term;
    double third_term;
    Err    err = NULL;

    /* A few security checks */
    if ((mat <= 0.0) || (betavol == 0.0))
    {
        *bsvol = 0.0;
        return NULL;
    }

    if ((beta > 1.0) || (beta < 0.0))
        return serror("Beta has to be between 0.0 and 1.0");

    if ((rho > 1.0) || (rho < -1.0))
        return serror("Rho has to be between -1 and 1");

    if (strike < 0.0)
        return serror("Strike has to be positive");

    /* The new variable coming from Optic expansion is zz and xx */
    zz = vovol / betavol * log(forward / strike) * pow(forward * strike, (1.0 - beta) / 2.);

    xx = log((sqrt(1.0 - 2.0 * rho * zz + zz * zz) + zz - rho) / (1.0 - rho));

    /* First term */
    first_term = pow(forward * strike, (1.0 - beta) / 2.0) *
                 (1.0 + 1.0 / 24.0 * pow((1.0 - beta) * log(forward / strike), 2.0));

    /* Second Term */
    second_term = pow(1.0 - beta, 2.0) / 24.0;
    second_term += 0.25 * rho * beta * vovol / betavol * pow(forward * strike, (1.0 - beta) / 2.0);
    second_term += (2.0 - 3.0 * rho * rho) / 24.0 * vovol * vovol / (betavol * betavol) *
                   pow(forward * strike, (1.0 - beta));
    second_term *= betavol * betavol * mat / pow(forward * strike, (1.0 - beta));
    second_term += 1.0;

    /* Third Term */
    if (fabs(zz) < 0.00005)
    {
        third_term = betavol / (1 + 0.5 * zz * rho + 1. / 6. * zz * zz * (3. * rho * rho - 1.));
    }
    else
        third_term = betavol * (zz / xx);

    /* The equivalent BlackScholes betastoch vol (constant elasticity and stoch vol) */
    *bsvol = third_term * second_term / first_term;

    /* Return a success message */
    return NULL;
}
/* ------------------------------------------------------------------------------- */
