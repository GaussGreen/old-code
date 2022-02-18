/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

#define INTEG_LIMIT 49.0         /* "infinity" for integration     */
#define NORMAL_LIMIT 6.0         /* "infinity" for arg.of normal dis. */
#define TINY_EPSILON 1.e-10      /* To prevent division by zero */
#define DEFAULT_SPRNUM_EPS 1.e-8 /* For Simpson */

static double NI_a;
static double NI_b;
static double NI_k;
static double NI_rho;
static double NI_sigX;
static double NI_sigY;
static double NI_sigZ;
static double NI_muX;
static double NI_muY;
static double NI_limitX;

/*************************************************************

        integration of AX-BY-K with AX-BY-K>0
        where X and Y are correlated log normal assets and
        K is strike

Essential functionality: For any a > 0, K > 0, and any b:

Writing x=rho*y+ksi*z, and integrating a first time on z, we get:

Premium(call) = disc * integral{Gaussian(u)*max{0,[S*N(d1)-Q*N(d2)]}}du
where:
    d1(u) = [logS/Q +.5*sigx*sigx*(1-rho*rho)*t]/[sigx*sqrt({1-rho*rho}*t)]
    d2(u) = d1 - sigx*sqrt[(1-rho*rho)*t]

     S(u) = a*FwdX*exp(-.5*rho*rho*sigx*sigx*t) *  exp(+rho*sigx*sqrt{t}*u)
     Q(u) = K + b*FwdY*exp(-.5*sigy*sigy*t) * exp(+sigy*sqrt{t}*u)

Note: d1 and d2 should be replaced by +infinity if Q <= 0

Premium(put) = K*disc - a*X + b*Y + Premium(call)

Numerical accuracy: convergence criterion 1/10^6
   integration is limited to best of O(h^2) or O(h/sigx*sqrt(1-rho^2)*t)^4
**************************************************************************/

static double int_in(double var)
{
    double kprime, d1, d2, bs, integrd;

    kprime = NI_a * exp((NI_sigX - NI_rho * NI_sigY) * var + NI_muX) -
             NI_k * exp(-NI_rho * NI_sigY * var);

    if ((NI_b > 0) || ((NI_k > 0) && (NI_a < 0)) ||
        ((NI_k > 0) && (NI_a > 0) && (var < NI_limitX)) ||
        ((NI_k < 0) && (NI_a < 0) && (var > NI_limitX)))
    {
        if (fabs(kprime) < TINY_EPSILON)
        {
            kprime = (NI_b > 0 ? TINY_EPSILON : -TINY_EPSILON);
        }

        // careful => no catch of (kprime / NI_b <0)   (J.B)

        d1 = log(kprime / NI_b) - NI_muY;

        if (NI_sigZ < TINY_EPSILON)
        {
            if (d1 >= 0)
            {
                d1 = NORMAL_LIMIT;
            }
            else
            {
                d2 = NORMAL_LIMIT;
            }
        }
        else
        {
            d1 /= NI_sigZ;
        }

        d2 = d1 - NI_sigZ;
    }
    else
    {
        d1 = d2 = -NORMAL_LIMIT;
    }

    if (NI_b > 0)
    {
        bs = kprime * norm(d1) - NI_b * exp(NI_muY + NI_sigZ * NI_sigZ / 2.0) * norm(d2);
    }
    else
    {
        bs = kprime * norm(-d1) - NI_b * exp(NI_muY + NI_sigZ * NI_sigZ / 2.0) * norm(-d2);
    }

    integrd = exp(-0.5 * var * var) * exp(NI_rho * NI_sigY * var) *
              bs; /* Here, we forget the sqrt(2*PI)*/

    return (integrd);

} /* END int_in(...) */

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optsprnum(...)
 *
 * PURPOSE      	: ??
 *
 * DESCRIPTION  	: called by int_doub_addin (in option_tools.c) <--
 * 			@spread_numeric(...)
 *
 * CALLS		: int_in(...)
 *		: sm_qsimp(...)
 *		: SELF
 *
 * PARAMETERS   	: fwdx    	- forward price of 1st underlying
 *              	: fwdy     	- forward price of 2nd underlying
 *              	: sigx        	- vol of 1st und
 *              	: sigy        	- vol of 2nd und
 *              	: rho        	- correlation
 *              	: disc        	- discount factor
 *              	: a          	- gearing of 1st und
 *              	: b          	- gearing of 2nd und
 *              	: k          	- strike
 *              	: t          	- maturity
 *              	: call_put      - option type
 *              	: greek         - greek wanted
 *
 * RETURNS      	: ??          	- ??
 *
 *******************************************************************************/

double srt_f_optsprnum(
    double         fwdx,
    double         fwdy,
    double         sigx,
    double         sigy,
    double         rho,
    double         disc,
    double         a,
    double         b,
    double         k,
    double         t,
    SrtCallPutType call_put,
    SrtGreekType   greek)
{
    double result;
    double answer;
    double price;
    double shift, shiftx, shifty;

    double glimithigh, glimitlow;

    /* set globals for int_in  */
    NI_sigX = sigx * sqrt(t);
    NI_sigY = sigy * sqrt(t);
    NI_sigZ = sqrt(1.0 - rho * rho) * NI_sigY;
    NI_muX  = log(fwdx) - sigx * sigx * t / 2.0;
    NI_muY  = log(fwdy) - sigy * sigy * t / 2.0;
    NI_rho  = rho;
    NI_a    = a;
    NI_b    = b;
    NI_k    = k;
    if (k * a > 0.0)
    {
        NI_limitX = (log(k / a) - NI_muX) / NI_sigX;
    }
    else
    {
        NI_limitX = 0.0;
    }

    /* sets integration limits */
    glimitlow  = (fabs(NI_limitX) > INTEG_LIMIT ? -INTEG_LIMIT - fabs(NI_limitX) : -INTEG_LIMIT);
    glimithigh = (fabs(NI_limitX) > INTEG_LIMIT ? +INTEG_LIMIT + fabs(NI_limitX) : +INTEG_LIMIT);

    /* Different limit cases depending on a, b, k */
    if ((k * a > 0.0) && (b > 0))
    {
        if (k > 0)
        {
            glimitlow = NI_limitX;
        }
        else
        {
            glimithigh = NI_limitX;
        }
    }

    /* Different cases depending on a,b, and k */
    if ((b > 0) && (k > 0) && (a < 0))
    {
        result = 0;
    }
    else
    {
        result = sm_qsimp(int_in, glimitlow, glimithigh, DEFAULT_SPRNUM_EPS);
        result *= INV_SQRT_TWO_PI;
    }

    /* call_put parity  */
    if (call_put == SRT_PUT)
        result = result + (k - (a * fwdx) + (b * fwdy));

    price  = result * disc;
    result = price;
    switch (greek)
    {
    case PREMIUM:
        return (price);
        break;

    case DELTAX:
        shift  = fwdx / 10000;
        answer = (srt_f_optsprnum(
                      fwdx + shift, fwdy, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case DELTAY:
        shift  = fwdy / 10000;
        answer = (srt_f_optsprnum(
                      fwdx, fwdy + shift, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case GAMMAX:
        shift  = fwdx / 1000;
        answer = srt_f_optsprnum(
            fwdx + shift, fwdy, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer += srt_f_optsprnum(
            fwdx - shift, fwdy, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer -= 2 * result;
        answer /= shift * shift;
        return (answer);
        break;

    case GAMMAY:
        shift  = fwdy / 1000;
        answer = srt_f_optsprnum(
            fwdx, fwdy + shift, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer += srt_f_optsprnum(
            fwdx, fwdy - shift, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer -= 2 * result;
        answer /= shift * shift * disc * disc;
        return (answer);
        break;

    case GAMMAXY:
        shiftx = fwdx / 1000;
        shifty = fwdy / 1000;
        answer = srt_f_optsprnum(
            fwdx + shiftx, fwdy + shifty, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer += srt_f_optsprnum(
            fwdx - shiftx, fwdy - shifty, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer -= srt_f_optsprnum(
            fwdx + shiftx, fwdy - shifty, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer -= srt_f_optsprnum(
            fwdx - shiftx, fwdy + shifty, sigx, sigy, rho, disc, a, b, k, t, call_put, PREMIUM);
        answer /= 4 * shiftx * shifty;
        return (answer);
        break;

    case VEGAX:
        shift  = GVOPT.vol_add;
        answer = (srt_f_optsprnum(
                      fwdx, fwdy, sigx + shift, sigy, rho, disc, a, b, k, t, call_put, PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case VEGAY:
        shift  = GVOPT.vol_add;
        answer = (srt_f_optsprnum(
                      fwdx, fwdy, sigx, sigy + shift, rho, disc, a, b, k, t, call_put, PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case THETA:
        shift  = YEARS_IN_DAY;
        answer = srt_f_optsprnum(
                     fwdx,
                     fwdy,
                     sigx,
                     sigy,
                     rho,
                     disc * exp(-shift * log(disc) / t),
                     a,
                     b,
                     k,
                     t - shift,
                     call_put,
                     PREMIUM) -
                 result;
        return (answer);
        break;

    default:
        return UNKNOWN_GREEK;
    }

} /* END srt_f_optsprnum(...) */

/******************************************************************************/
