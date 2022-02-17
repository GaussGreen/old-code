/******************************************************************************
 * Module:	Q3
 * Submodule:
 * File: bas.c	
 * Function:	
 * Author:	Interest Rates DR
 * Revision:	$Header$
 *****************************************************************************/
#include <math.h>


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "q3.h"


/* Coefficients for cumulative normal */
#define     h0   0.2316419
#define     h1   0.31938153
#define     h2  -0.356563782
#define     h3   1.781477937
#define     h4  -1.821255978
#define     h5   1.330274429

/* Coefficients for inverse cumulative normal */
#define   a0     2.50662823884
#define   a1   -18.61500062529
#define   a2    41.39119773534
#define   a3   -25.44106049637
#define   b0    -8.47351093090
#define   b1    23.08336743743
#define   b2   -21.06224101826
#define   b3     3.13082909833
#define   c0     0.3374754822726147
#define   c1     0.9761690190917186
#define   c2     0.1607979714918209
#define   c3     0.0276438810333863
#define   c4     0.0038405729373609
#define   c5     0.0003951896511919
#define   c6     0.0000321767881768
#define   c7     0.0000002888167364
#define   c8     0.0000003960315187

/* coefficients for double precision cumulative normal -- courtesy of ALIB */

#define SQRT2   1.414213562373095049     /* sqrt(2) */
#define SQRTPI  1.772453850905516027     /* sqrt(pi) */

/* Coefficients in expression of erf(x) for -0.46875<=x<=0.46875 */
#define P10 3209.377589138469472562    /* Numerator */
#define P11 377.4852376853020208137
#define P12 113.8641541510501556495
#define P13 3.161123743870565596947
#define P14 0.1857777061846031526730
#define Q10 2844.236833439170622273   /* Denominator */
#define Q11 1282.616526077372275645
#define Q12 244.0246379344441733056
#define Q13 23.60129095234412093499
#define Q14 1.0

/* Coefficients in expression of erfc(x) for 0.46875<=x<=4.0 */
#define P20 1230.33935479799725272  /* Numerator */
#define P21 2051.07837782607146532
#define P22 1712.04761263407058314
#define P23 881.952221241769090411
#define P24 298.635138197400131132
#define P25 66.1191906371416294775
#define P26 8.88314979438837594118
#define P27 0.564188496988670089180
#define P28 2.15311535474403846343e-8
#define Q20 1230.33935480374942043  /* Denominator */
#define Q21 3439.36767414372163696
#define Q22 4362.61909014324715820
#define Q23 3290.79923573345962678
#define Q24 1621.38957456669018874
#define Q25 537.181101862009857509
#define Q26 117.693950891312499305
#define Q27 15.7449261107098347253
#define Q28 1.0

/* Coefficients in expression of erfc(x) for x>= 4.0 */
#define P30 -6.58749161529837803157E-4    /* Numerator */
#define P31 -1.60837851487422766278E-2
#define P32 -1.25781726111229246204E-1
#define P33 -3.60344899949804439429E-1
#define P34 -3.05326634961232344035E-1
#define P35 -1.63153871373020978498E-2
#define Q30  2.33520497626869185443E-3    /* Denominator */
#define Q31  6.05183413124413191178E-2
#define Q32  5.27905102951428412248E-1
#define Q33  1.87295284992346047209
#define Q34  2.56852019228982242072
#define Q35  1.0



/*f----------------------------------------------------------------------------
 * Normal density function
 */

double  NormDens (double x) 
{
    return (INVSQRT2PI * exp (-0.5 * x * x));
}


/*f----------------------------------------------------------------------------
 * Double precision cumulative normal function
 */
double NormCum (double x)
{
   return (x>0.0) ? 1. - 0.5 * ExpCErrFcn(0,0,x) : 0.5 * ExpCErrFcn(0,0,x);
}


/*f----------------------------------------------------------------------------
 * Inverse cumulative normal distribution
 *
 * Based on Risk Magazine
 */

double NormCumInv (double prob)
{
    double  t;
    double  x;
    double  r;
   
    t = (prob < 0.5) ? (1. - prob) : prob;
   
    x = t - 0.5;
    if (fabs (x) < 0.42) {
        r = x * x;
        r = x * (((a3 * r + a2) * r + a1) * r + a0) 
            / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.);
        return (prob < 0.5) ? -r : r;
    } else {
        r = t;
        if (x > 0.) r = 1. - t;
        r = log (- log (r));
        r = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r 
            * (c6 + r * (c7 + r * c8)))))));
        if (x < 0.) r = -r;
        return (prob < 0.5) ? -r : r;
   }
} /* NormCumInv */


/*f----------------------------------------------------------------------------
 * Exp(a^2) times difference of cumulative normal distribution
 * Computes exp(b+a^2/2)[N(x2-a)-N(x1-a)] (via erfc)
 */

double  BSQInt (double  a, double b, double x1, double x2)
{
    double erfc1, erfc2; 
    double erfc = 0; 
    int    sign = 0;

    erfc1 = ExpCErrFcn (a,b,x1); 
    erfc2 = ExpCErrFcn (a,b,x2);

    /* process x1 */
    if (x1-a > 0) {
        erfc += erfc1; 
        sign -= 1;
    } else {
        erfc -= erfc1;
    }

    /* process x2 */
    if (x2-a > 0) {
      erfc -= erfc2;
      sign += 1;
    } else {
      erfc += erfc2;
    }

    /* NormCum has an extra coefficient of 1/2 */
    erfc *= 0.5;

    /* put back value at infinity as in NormCum */
    if (sign) erfc += sign * exp (b + 0.5 * a * a);
 
    return erfc;
} /* ExpA2DiffNormCum */


/*f----------------------------------------------------------------------------
 * Single Q Black Pricer
 *
 * Price of a Call/Put with Vega using Q version of Black&Scholes.
 */

int BSQPricer(
    double  Y,                  /* Fwd yield                     */
    double  K,                  /* Strike                        */
    double  T,                  /* Option expiration in years    */
    double  s,                  /* Annualized volatility         */
    double  Q,                  /* Q weight                      */
    long    I,                  /* Instrument                    */
    double *P)                  /* Price & Vega                  */
{

    double  C;	                /* Call price                    */
    double  V=0;                /* Vega                          */
    double  d;	                /* d in N(d) in Black & Scholes  */
    double  st;	                /* Sigma * sqrt(T) * Q           */
    double  Y1, K1;             /* Modified fwd and strike       */
    double  r;                  /* Adjusted strike               */
    long    vOn;                /* Vega calc on or off           */

    if (Y < Q3_MIN_FWD) return FAILURE;

    /* compute vega? */
    vOn = (Q3_PAY_TYPE(I) == Q3_VEGA);

    /* for zero volatility return intrinsic */
    if (s * sqrt(T) < Q3_MIN_VOL) {
        C = MAX(Y-K, 0.);
        V = 0.;

        switch(Q3_COP_TYPE(I)) {
            case Q3_CALL: 
                P[0] = C; 
                if (vOn) P[1] = V; 
                break;
            case Q3_PUT:            
                P[0] = C -(Y-K); 
                if (vOn) P[1] = V; 
                break;
	    default:
	        return FAILURE;
        }

        return SUCCESS;
    }

    if (fabs(Q) > Q3_Q_SHIFT) {
        st  = s * sqrt(T) * Q;
        r   = (K/Y-1.) * Q + 1.;
        /* if K is not in the range of the distribution */
        if (r < TINY) {
	    C = MAX(Y-K, 0.);
            V = 0.;
	} else {
            d   = - log(r) / st;
            st /= 2.;

            Y1  = Y / Q;
            K1  = K - Y + Y1;

            C   = Y1 * NormCum(d+st) - K1 * NormCum(d-st);        
            if (vOn) V = Y * exp(-(d*d+st*st)/2.) * sqrt(r*T) * INVSQRT2PI;
	}
    } else {
        st  = s * sqrt(T);
        d   = - (K/Y-1.)/st;

        C   = Y * st * exp(-d*d/2) * INVSQRT2PI + (Y-K) * NormCum(d);
        if (vOn) V = Y * exp(-(d*d)/2.) * sqrt(T) * INVSQRT2PI;
    }

    switch(Q3_COP_TYPE(I)) {
        case Q3_CALL: 
            P[0] = C; 
            if (vOn) P[1] = V; 
            break;
        case Q3_PUT:            
            P[0] = C- (Y-K); 
            if (vOn) P[1] = V; 
            break;
	default:
	    return FAILURE;
    }

    return SUCCESS;

}  /* BSQPricer */


/*f----------------------------------------------------------------------------
 *      Implied vol of a call or put using single Q Black-Scholes.
 */
#define VOLDELTA 1E-4
#define PREQERR  1E-8

int     BSQImpVol   (
    double yield,                     /* (I) Fwd yield              */
    double strike,                    /* (I) Strike                 */
    double expiry,                    /* (I) Option expiration      */
    double Q,                         /* (I) Q parameter            */
    double price,                     /* (I) Price of option        */
    long   optType,                   /* (I) Call or put            */
    double volGuess,                  /* (I) Initial vol guess      */
    double *impVol)                   /* (O) Implied BS vol         */
{
    int    found, count;
    double Pr, Pr2, dPr, scale, iVol;

    static char      routine[]="BSQImpVol";
    int              status = FAILURE;

        
    scale = yield * sqrt(expiry / 6.28);
    iVol  = volGuess;
    count = 0;
    found = FALSE;

    do {
        if (BSQPricer (yield,
                       strike,
                       expiry,
                       iVol,
                       Q,
                       optType,
                       &Pr) == FAILURE) goto RETURN;

        if (BSQPricer (yield,
                       strike,
                       expiry,
                       iVol+VOLDELTA,
                       Q,
                       optType,
                       &Pr2) == FAILURE) goto RETURN;

        if (fabs(Pr - price) > scale * PREQERR) {
            dPr = (Pr2 - Pr) / VOLDELTA;
            if (fabs(dPr) < Q3_MQ_RESN) {
                Q3ErrMsg("%s: numerical derivative is 0.\n", routine);
                goto RETURN;
            }

            iVol -= (Pr - price) / dPr; 
        } else {
            found = TRUE;
        }
        count++;
    } while ((found != TRUE) && (count < 10));
                  
    /* failure if exceed 10 iterations */
    if (found == TRUE) {
        *impVol = iVol;
        status = SUCCESS;
    } else {
        Q3ErrMsg ("%s: exceeded maximum nb iterations.\n", routine);        
        *impVol = -999;
    }

  RETURN:

    return status;

} /* BSQImpVol */


/*f----------------------------------------------------------------------------
 * Cumulative error function weighted by exponential factor. Computes
 * exp(b+a^2/2)erfc(x-a)
 *
 * Alib comment:  The routine has a relative accuracy no worse than 1.0E-14, 
 * where relative accuracy is defined as (computed - truth)/truth, and truth
 * comes from a continued fraction calculation.  This is essentially 
 * accurate to the next to last decimal digit of machine accuracy on the Sun.
 */

double ExpCErrFcn (double a, double b, double x)
{
    double numerator;            /* numerator of polynomial in expression */
    double denominator;          /* denominator of polynomial in expression */
    double y,y2;                 /* y = abs(x)/sqrt(2), y2 = y*y */
    double erfc;                 /* return value */

    double W = b + 0.5 * a * a;

    y  = fabs(x-a) / SQRT2;
    y2 = y * y;

    if (y < 0.46875) {
        numerator   = P10 + y2*(P11 + y2*(P12 + y2*(P13 +y2*P14)));
        denominator = Q10 + y2*(Q11 + y2*(Q12 + y2*(Q13 +y2*Q14)));
        erfc = exp(W) * (1 - y * numerator / denominator);
    }
    else if (y < 4.0) {
        numerator   = P20 + y*(P21 + y*(P22 + y*(P23 +
                            y*(P24 + y*(P25 + y*(P26 + y*(P27 + y*P28)))))));
        denominator = Q20 + y*(Q21 + y*(Q22 + y*(Q23 +
                            y*(Q24 + y*(Q25 + y*(Q26 + y*(Q27 + y*Q28)))))));
        erfc = exp(-y2 + W) * numerator / denominator;
    }
    else /* (y > 4.0) */ {
        double z2 = 1/y2; 
        numerator   = P30 + z2*(P31 + z2*(P32 + z2*(P33 + z2*(P34 +z2*P35))));
        denominator = Q30 + z2*(Q31 + z2*(Q32 + z2*(Q33 + z2*(Q34 +z2*Q35))));
        erfc = (exp(-y2 + W)/y) * 
               (1.0 / SQRTPI + numerator / (denominator * y2)); 
    }
    
    return erfc;
} /*ExpCerrFcn */
