/*-----------------------------------------------------------------------
C SOURCE FILE:  binorm.c

CREATED BY:     10/14/94  By Neil Yang
MODIFIED BY:    9/2/99 By Charles Irwin/Peter Taylor

PURPOSE:        Generate Bivariate Normal Cumulative

$Header$
---------------------------------------------------------------------- 
Proprietary software, whether developed for Morgan by in-house
staff or external contractors, must not be copied or removed from Morgan
premises without the approval of a Senior Vice President and Audit.

This proprietary software has been developed strictly for Morgan's own
internal use.  Any use or misuse, intentional or otherwise which contradicts
or places this policy in jeopardy is strictly forbidden.  Also, any actions or 
inactions resulting in a breach of this goal will be dealt with harshly.

Do Not Distribute this to anyone outside the Analytics Library
Group without authorization from its management. 

Copyright 1997 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */

#include <math.h>
#include <float.h>
#include "error2.h"
#include "assert.h"
#include "proba_utils.h"

#ifndef MAX
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)>(b) ? (b) : (a))
#endif

#define IS_ALMOST_ZERO(x) (fabs(x)<3e-15)?1:0

#ifndef PI
#define PI	   3.14159265358979323846264338328      /* pi */
#endif


static double BiNormalCumAux (double a, double b, double rho);
static double BNC1 (double a, double b, double rho);
static double Phi (double a, double b, double rho);

/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --
Function: BiNormalCum

Created by: 10/14/94 Neil Yang

Modified by: 9/19/97 David Hait
Modified by: 9/2/99 Charles Irwin to use an approximation

Old Description: This function calculates the value for cumulative
bivariate normal distribution by performing a numeric quadrature using
transformation of variables.  See Analytics Library Technical Memo for
mathematical background.

New Description: The new method is based on an approximation that is
well known and fast (Gauss quadrature). The weights and abscissae were
computed to 15 significant digits using Mathematica. The accuracy for
a 10 point quadrature in two dimensions is approximately 10^(-12).
--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- */

int BiNormalCum(
   double a,            /* (I) normal point for first variable */
   double b,            /* (I) normal point for second variable */
   double rho,          /* (I) correlation of the two variables */
   double *biNormalCum) /* (O) probability */
{
    static char routine[]="BiNormalCum";
    int status = FAILURE;

    if (rho < -1. || 
        rho > 1.)
    {
        DR_Error("%s: abs(rho) (%f) >= 1.\n", routine, fabs(rho));
        goto done; /* failure */
    }

    *biNormalCum = BiNormalCumAux(a, b, rho);

    /*
     * Perhaps the routine returned something just outside the range 
     * [0,1]. This may be due to instability in GtoNormalCum.
     * We (silently) correct it here.
     */

    if ((*biNormalCum >= -2.0) &&
        (*biNormalCum <= 3.0))
    {
        *biNormalCum = MIN(*biNormalCum, 1.0);
        *biNormalCum = MAX(*biNormalCum, 0.0);
    }
    else
    {
        /*
         * We do a sanity check for numbers grossly out of range.
         * The routine should be stable enough that this case never 
         * actually occurs.
         */
        DR_Error("%s: Probability (%g) outide the range [0,1] returned "
                  "from integration.\n",
                  routine, *biNormalCum);
        goto done;
    }

    status = SUCCESS;

 done:
    if (status==FAILURE)
    {
        DR_Error ("%s: Failed on inputs a=%g, b=%g, rho=%g.\n",
                   routine, a, b, rho);
    }
    return (status);
}

/*
** FUNCTION: BiNormalCum
** AUTHOR:   Peter Taylor (27 July 1999)
**
** This is an algorithm extracted from a book. This is based on
** Drezner's method.
*/
static double BiNormalCumAux (double a, double b, double rho)
{
    double BNC; /* To be returned */

    if ((IS_ALMOST_ZERO(a)) &&
        (IS_ALMOST_ZERO(b)))
    {
        /*
         * An easy case where a and b are both zero.
         */
        return (0.25 + asin(rho)*(1./(2.*PI)));
    }

    /*
     * We trap the case where rho is either +1 or -1. In this case,
     * the bivariate normal degenerates.
     */
    if (IS_ALMOST_ZERO(fabs(rho) -  1.0))
    {
        /*
        ** Degenerate case - perfectly correlated.
        */
        if (rho > 0.0)
        {
            /*
            ** Perfectly correlated
            */
            BNC = NormalCum (MIN (a,b));
        }
        else
        {
            /*
            ** Perfectly negatively correlated
            */
            if (a < -b)
            {
                BNC = 0.0;
            }
            else
            {
                BNC = NormalCum(a) + NormalCum(b) - 1.0;
            }
        }
        return BNC;
    }

    if (a*b*rho <= 0.0)
    {
        BNC = BNC1 (a,b,rho);
    }
    else
    {
        /*
         * a*b*rho > 0 in this case.
         */
        double rho1;
        double rho2;
        double den2;
        double den;

        if (fabs(a) <= fabs(b))
        {
            double z = a/b;

            den2 = (z-1)*(z-1)+2*(1-rho)*z;
            if ((den2 <= 0) ||
                (IS_ALMOST_ZERO(den2)))
            {
                /*
                 * This case implies that z=1 or z=-1, which implies that
                 * rho=1 or rho=-1.
                 *
                 * This case is theoretically impossible, but we do a reasonable
                 * thing anyway.
                 */
                if (rho > 0.0)
                {
                    rho = 1.0;
                }
                else
                {
                    rho = -1.0;
                }
                return (BiNormalCumAux(a, b, rho));
            }

            den = sqrt(den2);

            rho1 = (rho*z - 1)/den;
            rho2 = (rho   - z)/den;

            if (((a > 0.0) && (b < 0.0)) ||
                ((a < 0.0) && (b > 0.0)))
            {
                /*
                 * a and b are of opposite sign.
                 */
                rho1 = -rho1;
            }
        }
        else
        {
            double z = b/a;

            den2 = (z-1)*(z-1)+2*(1-rho)*z;
            if ((den2 <= 0) ||
                (IS_ALMOST_ZERO(den2)))
            {
                /*
                 * This case implies that z=1 or z=-1, which implies that
                 * rho=1 or rho=-1.
                 *
                 * This case is theoreticall impossible, but we do a reasonable
                 * thing anyway.
                 */
                if (rho > 0.0)
                {
                    rho = 1.0;
                }
                else
                {
                    rho = -1.0;
                }
                return (BiNormalCumAux(a, b, rho));
            }

            den = sqrt(den2);

            rho1 = (rho   - z)/den;
            rho2 = (rho*z - 1)/den;

            if (((a > 0.0) && (b < 0.0)) ||
                ((a < 0.0) && (b > 0.0)))
            {
                /*
                 * a and b are of opposite sign.
                 */
                rho2 = -rho2;
            }
        }

        /*
         * Extra safety check to make absolutely sure rho1 and rho2 are in the 
         * correct range [-1, 1].
         */

        rho1 = MIN(rho1,  1.0);
        rho1 = MAX(rho1, -1.0);
        rho2 = MIN(rho2,  1.0);
        rho2 = MAX(rho2, -1.0);

        if (((a > 0.0) && (b > 0.0)) ||
            ((a < 0.0) && (b < 0.0)))
        {
            /*
             * a and b are of the same sign.
             */
            
            BNC = BNC1 (a,0.0,rho1) + BNC1 (b,0.0,rho2);
        }
        else
        {
            /*
             * a and b are of different sign.
             */
            BNC = BNC1 (a,0.0,rho1) + BNC1 (b,0.0,rho2) - 0.5;
        }
    }
    return BNC;
}

/*
 * This routine performs a 2D numerical integration using a 10 point Gauss 
 * quadrature with weighting function Exp(-x^2) over the interval [0,Infinity].
 *
 * Previous methods integrated over [-Infinity,Infinity] and used the Hermite
 * polynomials's roots as the abscissae, but this leads to an instability at the 
 * origin. It is better to derive quadrature results for the half-infinite 
 * interval.
 *
 * The coefficients were generated by computing the orthogonal polynomials with
 * respect to the inner product <f,g>=Integral(0,Infinity,Exp(-x^2)*f(x)*g(x)).
 *
 * These were computed algebraically exact with Mathematica. Then, the roots of
 * this polynomial were found to 100 decimal places. After this, the weights
 * were easily determined (also to 100 decimal places).
 *
 * The resulting weights and abscissae were checked against the polynomials
 * x^i for i=0, 1, 2, ..., 19.
 *
 * The difference between the theoretical value of the integral and the 
 * quadrature value was on the order of 10^(-90).
 *
 * The values listed here were rounded to exactly 15 digits in scientific 
 * notation.
 *
 * The array ww[i][j] holds the result of multiplying weight i against weight
 * j (We precomputed these in Mathematica for speed).
 *
 * For details, see Z. Drezner, "Computation of the Bivariate Normal Integral"
 * Math. of Computation, Vol 32, Number 141, pp. 277-279 (1978).
 */

static double Phi (double a1, double b1, double rho)
{
    int    i;
    int    j;
    double sum;
    double BNC;

    static double x[10] = 
    {
        3.87385243256994e-2,
        1.98233304012949e-1,
        4.65201111814507e-1,
        8.16861885591907e-1,
        1.23454132402774e+0,
        1.70679814968865e+0,
        2.22994008892444e+0,
        2.80910374689825e+0,
        3.46387241949537e+0,
        4.25536180636561e+0,
    };                   /* Roots of the 10th degree orthogonal polynomial */

    static double ww[10][10] = 
    {
         {
             9.71251592540161e-3,
             2.05656611704430e-2,
             2.48402225754016e-2,
             1.95807584549757e-2,
             9.57910843855530e-3,
             2.66331292395503e-3,
             3.74956200622751e-4,
             2.25572193453660e-5,
             4.28242844331612e-7,
             1.22967113105880e-9,
         },
         {
             2.05656611704430e-2,
             4.35465354832845e-2,
             5.25976590213905e-2,
             4.14610639445273e-2,
             2.02831789389435e-2,
             5.63940297298984e-3,
             7.93946927345222e-4,
             4.77635386718788e-5,
             9.06778151287977e-7,
             2.60375375717955e-9,
         },
         {
             2.48402225754016e-2,
             5.25976590213905e-2,
             6.35300536271685e-2,
             5.00787233660737e-2,
             2.44990265668761e-2,
             6.81154979071513e-3,
             9.58968361136659e-4,
             5.76911640119526e-5,
             1.09525149314805e-6,
             3.14494254883637e-9,
         },
         {
             1.95807584549757e-2,
             4.14610639445273e-2,
             5.00787233660737e-2,
             3.94754669765184e-2,
             1.93118044788806e-2,
             5.36932834442915e-3,
             7.55924299324728e-4,
             4.54761121433375e-5,
             8.63351963521471e-7,
             2.47905831828260e-9,
         },
         {
             9.57910843855530e-3,
             2.02831789389435e-2,
             2.44990265668761e-2,
             1.93118044788806e-2,
             9.44753338706181e-3,
             2.62673065355268e-3,
             3.69805942462381e-4,
             2.22473818155004e-5,
             4.22360660759302e-7,
             1.21278082822666e-9,
         },
         {
             2.66331292395503e-3,
             5.63940297298984e-3,
             6.81154979071513e-3,
             5.36932834442915e-3,
             2.62673065355268e-3,
             7.30319083684035e-4,
             1.02818435790040e-4,
             6.18551714843314e-6,
             1.17430407389781e-7,
             3.37193682946561e-10,
         },
         {
             3.74956200622751e-4,
             7.93946927345222e-4,
             9.58968361136659e-4,
             7.55924299324728e-4,
             3.69805942462381e-4,
             1.02818435790040e-4,
             1.44753587500177e-5,
             8.70831958198587e-7,
             1.65325144471073e-8,
             4.74720266981928e-11,
         },
         {
             2.25572193453660e-5,
             4.77635386718788e-5,
             5.76911640119526e-5,
             4.54761121433375e-5,
             2.22473818155004e-5,
             6.18551714843314e-6,
             8.70831958198587e-7,
             5.23889122553910e-8,
             9.94589645655809e-10,
             2.85589868155718e-12,
         },
         {
             4.28242844331612e-7,
             9.06778151287977e-7,
             1.09525149314805e-6,
             8.63351963521471e-7,
             4.22360660759302e-7,
             1.17430407389781e-7,
             1.65325144471073e-8,
             9.94589645655809e-10,
             1.88820214174986e-11,
             5.42184812670272e-14,
         },
         {
             1.22967113105880e-9,
             2.60375375717955e-9,
             3.14494254883637e-9,
             2.47905831828260e-9,
             1.21278082822666e-9,
             3.37193682946561e-10,
             4.74720266981928e-11,
             2.85589868155718e-12,
             5.42184812670272e-14,
             1.55684799095647e-16,
         },
    };                  /* Array of weights precomputed for 2D integration */

    sum = 0.0;

    for (i = 0; i < 10; ++i)
    {
        for (j = 0; j < 10; ++j)
        {
            double sumand;

            /*
             * Since Phi is called with all NEGATIVE arguments, the exponent 
             * is negative (check that each term is negative). Thus it is 
             * impossible that the exponent can cause an overflow problem. 
             * The worst possible case is that the exponent is -Inf, in which 
             * case exp(-Inf)=0.
             */

            sumand = ww[i][j] * exp (a1 * (x[i]+x[i]-a1) +
                                     b1 * (x[j]+x[j]-b1) +
                                     2.0 * rho * (x[i]-a1) * (x[j]-b1));
            sum += sumand;
        }
    }

    BNC = sqrt (1.0 - rho * rho) / PI * sum;

    return BNC;
}

static double BNC1 (double a, double b, double rho)
{
    /*
    ** Only use if a * b * rho is zero or negative.
    */

    double BNC; /* To be returned */

    double a1;
    double b1;
    double c = sqrt (2.0 * (1.0 - rho * rho));

    /*
     * NOTE: We must make a check here, because this can be called with a
     * rho (very close to 1) from BiNormalCum that makes this c zero.
     */

    if (IS_ALMOST_ZERO (c))
    {
        /*
        ** Degenerate case - perfectly correlated.
        */
        if (rho > 0.0)
        {
            /*
            ** Perfectly correlated
            */
            BNC = NormalCum (MIN (a,b));
        }
        else
        {
            /*
            ** Perfectly negatively correlated
            */
            if (a < -b)
            {
                BNC = 0.0;
            }
            else
            {
                BNC = NormalCum (a) + NormalCum(b) - 1.0;
            }
        }
        return BNC;
    }

    /*
     * This is now a "safe" division, although the answer might be +Inf or -Inf.
     */
    a1 = a/c;
    b1 = b/c;

    /*
     * Note: It is acceptable for a1 or b1 to be either +Inf or -Inf.
     * The Phi routine handles this case correctly, because it gets called
     * with NEGATIVE arguments.
     */
    if (a <= 0.0 && b <= 0.0 && rho <= 0.0)
    {
        BNC = Phi (a1,b1,rho);
    }
    else if (a <= 0.0 && b >= 0.0 && rho >= 0.0)
    {
        BNC = NormalCum(a) - Phi (a1, -b1, -rho);
    }
    else if (a >= 0.0 && b <= 0.0 && rho >= 0.0)
    {
        BNC = NormalCum(b) - Phi (-a1, b1, -rho);
    }
    else if (a >= 0.0 && b >= 0.0 && rho <= 0.0)
    {
        BNC = NormalCum(a) + NormalCum(b) - 1 + Phi (-a1, -b1, rho);
    }
    else
    {
        /*
         * This case should not occur under any circumstance.
         */
        assert(0);
        BNC = -999;
    }

    return BNC;
}
