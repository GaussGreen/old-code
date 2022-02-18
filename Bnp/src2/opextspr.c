/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

#define CONVERR 1.e-6            /* convergence criterion for integral  */
#define GAUSS_LIMIT 6.0          /* "infinity" for arg. of Gaussian     */
#define NORMAL_LIMIT 5.5         /* "infinity" for arg.of normal dis. */
#define TINY_EPSILON 1.e-20      /* To prevent division by zero */
#define DEFAULT_SPRNUM_EPS 1.e-4 /* For Simpson */

static double NI_scoef;
static double NI_qcoef;
static double NI_aeff; /* Note that 1/sqrt(2pi) has been     */
static double NI_apeff;
static double NI_beff; /* absorbed into aeff, beff, and keff */
static double NI_bpeff;
static double NI_keff;
static double NI_stddev;

static double        NI_nlimithigh;
static double        NI_nlimitlow;
static SrtMinmaxType NI_min_max;

/*************************************************************

        integration of aX+a' - min(bY+b',K) with aX+a' - min(bY+b',K) > 0
        where X and Y are correlated log normal assets and
        K is strike

Essential functionality: For any a,b > 0, K > 0, and any a',b':

Writing x=rho*y+ksi*z, and integrating a first time on z, we get:

Premium(call) = disc * integral{Gaussian(u)*max{0,[S*N(d1)-Q*N(d2)]}}du
where:
    d1(u) = [logS/Q +.5*sigx*sigx*(1-rho*rho)*t]/[sigx*sqrt({1-rho*rho}*t)]
    d2(u) = d1 - sigx*sqrt[(1-rho*rho)*t]

     S(u) = a*FwdX*exp(-.5*rho*rho*sigx*sigx*t) *  exp(+rho*sigx*sqrt{t}*u)
     Q(u) = K + b*FwdY*exp(-.5*sigy*sigy*t) * exp(+sigy*sqrt{t}*u)

Note: d1 and d2 should be replaced by +infinity if Q <= 0

Premium(put) = MIN/MAX(K*disc,b*Y+b') - (a*X + a') + Premium(call)

Numerical accuracy: convergence criterion 1/10^6
   integration is limited to best of O(h^2) or O(h/sigx*sqrt(1-rho^2)*t)^4
**************************************************************************/

static double int_in(double var)
{
    double s, q, d1, d2, bs, integrd;

    s = NI_aeff * exp(NI_scoef * var) + NI_apeff;

    if (NI_min_max == SRT_MIN)
        q = DMIN(NI_keff, NI_beff * exp(NI_qcoef * var) + NI_bpeff); /* s and q                   */
    else
        q = DMAX(NI_keff, NI_beff * exp(NI_qcoef * var) + NI_bpeff); /* s and q                   */

    if ((s / q) >= NI_nlimithigh) /* d1 and d2 are + infinity */
    {
        integrd = (s - q) * exp(-.5 * var * var);
        return (integrd);
    }
    else if ((s / q) <= NI_nlimitlow) /* d1 and d2 are - infinity */
    {
        integrd = 0.;
        return (integrd);
    }
    else /* d1 and d2 are order 1    */
    {
        if (NI_stddev > TINY_EPSILON)
        {
            d1 = log(s / q) / NI_stddev + 0.5 * NI_stddev;
            d2 = d1 - NI_stddev;
            bs = s * norm(d1) - q * norm(d2);
        }
        else
        {
            if (s > q)
                bs = s - q;
            else
                bs = 0.0;
        }

        integrd = exp(-0.5 * var * var) * bs; /* Here, we forget the sqrt(2*PI)*/
    }

    return (integrd);

} /* END int_in(...) */

/******************************************************************************/

/*******************************************************************************
*
* FUNCTION     	: srt_f_optextspr(...)
*
* PURPOSE      	: ??
*
* DESCRIPTION  	:
*
* CALLS		: int_in(...)
*			: sm_qsimp(...)
*			: SELF
*
* PARAMETERS   	: fwdx    	- forward price of 1st underlying
*              	: fwdy     	- forward price of 2nd underlying
*              	: sigx        	- vol of 1st und
*              	: sigy        	- vol of 2nd und
*              	: rho        	- correlation
*              	: disc        	- discount factor
*              	: a          	- gearing of 1st und  a*x + a'
*              	: a'          	-
*              	: b          	- gearing of 2nd und  b*x + b'
*              	: b'          	-
*              	: k          	- strike
*              	: t          	- maturity
*				: min_max       -
*              	: call_put      - option type
*              	: greek         - greek wanted


*
* RETURNS      	:           	- premium(or other greek desired)
*
*******************************************************************************/

double srt_f_optextspr(
    double         fwdx,
    double         fwdy,
    double         sigx,
    double         sigy,
    double         rho,
    double         disc,
    double         a,
    double         ap,
    double         b,
    double         bp,
    double         k,
    double         t,
    SrtMinmaxType  min_max,
    SrtCallPutType call_put,
    SrtGreekType   greek)
{
    double result;
    double answer;
    double price;
    double shift, shiftx, shifty;

    double glimithigh, glimitlow;

    double (*f)(double);

    /* set globals for int_in  */
    NI_scoef   = rho * sigx * sqrt(t);
    NI_qcoef   = sigy * sqrt(t);
    NI_aeff    = a * fwdx * exp(-0.5 * NI_scoef * NI_scoef);
    NI_apeff   = ap;
    NI_beff    = b * fwdy * exp(-0.5 * NI_qcoef * NI_qcoef);
    NI_bpeff   = bp;
    NI_keff    = k;
    NI_stddev  = sigx * sqrt(t * (1 - rho * rho));
    NI_min_max = min_max;

    /* sets integration limits */
    NI_nlimithigh = exp((NORMAL_LIMIT + 0.5 * NI_stddev) * NI_stddev);
    NI_nlimitlow  = 1.0 / NI_nlimithigh;
    glimithigh    = (NI_qcoef > NI_scoef) ? (GAUSS_LIMIT + NI_qcoef) : (GAUSS_LIMIT + NI_scoef);
    glimitlow     = (NI_scoef > 0.) ? (-GAUSS_LIMIT) : (-GAUSS_LIMIT + NI_scoef);

    f      = int_in;
    result = sm_qsimp(f, glimitlow, glimithigh, DEFAULT_SPRNUM_EPS);
    result *= INV_SQRT_TWO_PI;

    /* call_put parity  */
    if (call_put == SRT_PUT)
        if (NI_min_max == SRT_MIN)
            result = result + (-(a * fwdx + ap) + DMIN(k, (b * fwdy + bp)));
        else
            result = result + (-(a * fwdx + ap) + DMAX(k, (b * fwdy + bp)));

    price = result * disc;

    switch (greek)
    {
    case PREMIUM:
        return (price);
        break;

    case DELTAX:
        shift  = fwdx / 10000;
        answer = (srt_f_optextspr(
                      fwdx + shift,
                      fwdy,
                      sigx,
                      sigy,
                      rho,
                      disc,
                      a,
                      ap,
                      b,
                      bp,
                      k,
                      t,
                      min_max,
                      call_put,
                      PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case DELTAY:
        shift  = fwdy / 10000;
        answer = (srt_f_optextspr(
                      fwdx,
                      fwdy + shift,
                      sigx,
                      sigy,
                      rho,
                      disc,
                      a,
                      ap,
                      b,
                      bp,
                      k,
                      t,
                      min_max,
                      call_put,
                      PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case GAMMAX:
        shift  = fwdx / 1000;
        answer = srt_f_optextspr(
            fwdx + shift,
            fwdx,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer += srt_f_optextspr(
            fwdx - shift,
            fwdy,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer -= 2 * result;
        answer /= shift * shift;
        return (answer);
        break;

    case GAMMAY:
        shift  = fwdy / 1000;
        answer = srt_f_optextspr(
            fwdx,
            fwdy + shift,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer += srt_f_optextspr(
            fwdx,
            fwdy - shift,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer -= 2 * result;
        answer /= shift * shift * disc * disc;
        return (answer);
        break;

    case GAMMAXY:
        shiftx = fwdx / 1000;
        shifty = fwdy / 1000;
        answer = srt_f_optextspr(
            fwdx + shiftx,
            fwdy + shifty,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer += srt_f_optextspr(
            fwdx - shiftx,
            fwdy - shifty,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer -= srt_f_optextspr(
            fwdx + shiftx,
            fwdy - shifty,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer -= srt_f_optextspr(
            fwdx - shiftx,
            fwdy + shifty,
            sigx,
            sigy,
            rho,
            disc,
            a,
            ap,
            b,
            bp,
            k,
            t,
            min_max,
            call_put,
            PREMIUM);
        answer /= 4 * shiftx * shifty;
        return (answer);
        break;

    case VEGAX:
        shift  = GVOPT.vol_add;
        answer = (srt_f_optextspr(
                      fwdx,
                      fwdy,
                      sigx + shift,
                      sigy,
                      rho,
                      disc,
                      a,
                      ap,
                      b,
                      bp,
                      k,
                      t,
                      min_max,
                      call_put,
                      PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case VEGAY:
        shift  = GVOPT.vol_add;
        answer = (srt_f_optextspr(
                      fwdx,
                      fwdy,
                      sigx,
                      sigy + shift,
                      rho,
                      disc,
                      a,
                      ap,
                      b,
                      bp,
                      k,
                      t,
                      min_max,
                      call_put,
                      PREMIUM) -
                  result) /
                 shift;
        return (answer);
        break;

    case THETA:
        shift  = YEARS_IN_DAY;
        answer = srt_f_optextspr(
                     fwdx,
                     fwdy,
                     sigx,
                     sigy,
                     rho,
                     disc * exp(-shift * log(disc) / t),
                     a,
                     ap,
                     b,
                     bp,
                     k,
                     t - shift,
                     min_max,
                     call_put,
                     PREMIUM) -
                 result;
        return (answer);
        break;

    default:
        return UNKNOWN_GREEK;
    }

} /* END srt_f_optsprnum(...) */