/****************************************************************************/
/*      Black & Scholes formulas.                                           */
/****************************************************************************/
/*      BAS.c                                                               */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esl_bas.h"
#include "esl_error.h"



/*****  NormalH  ************************************************************/
/**
*      Normal cumulative distribution accrding to J.Hull
*/
double  NormalH (double  x)
{
    /* Coefficients for cumulative normal */
    const double h0=  0.2316419;
    const double h1=  0.31938153;
    const double h2= -0.356563782;
    const double h3=  1.781477937;
    const double h4= -1.821255978;
    const double h5=  1.330274429;

    double  k;
    double  y;
    double  norm;

    if (x > 0) y = x; else y = -x;

    k = 1. / (1. + h0 * y);

    norm  = exp (-0.5 * y * y) * INVSQRT2PI;
    norm *= k * (h1 + k * (h2 + k * (h3 + k * (h4 + k * h5))));
    if (x > 0) norm  = 1. - norm;

    return (norm);
}


/*****  NormalHP  ************************************************************/
/*                                                                           */
/*      Normal cumulative distribution according to W. J. Cody               */
/*      This function evaluates near-minimax approximations from             */
/*      "Rational Chebyshev approximations for the error function",          */
/*      Math. Comp., 1969, PP. 631-638.                                      */
/*                                                                           */
/*      The original algorithm uses rational functions that theoretically    */
/*      approximate erf(x) to at least 18 significant decimal digits.        */
/*      The result is then mapped into the cumulative normal distribution    */
/*      N(u):=(erfc(-u/sqrt(2))/2;                                           */
/*                                                                           */
double  NormalHP (double u)
{

    const double SQRT2   = 1.414213562373095049;     /* sqrt(2)  */
    const double SQRTPI  = 1.772453850905516027;     /* sqrt(pi) */

    const double a0 = 1.161110663653770e-002;
    const double a1 = 3.951404679838207e-001;
    const double a2 = 2.846603853776254e+001;
    const double a3 = 1.887426188426510e+002;
    const double a4 = 3.209377589138469e+003;

    const double b0 = 1.767766952966369e-001;
    const double b1 = 8.344316438579620e+000;
    const double b2 = 1.725514762600375e+002;
    const double b3 = 1.813893686502485e+003;
    const double b4 = 8.044716608901563e+003;

    const double c0 = 2.15311535474403846e-8;
    const double c1 = 5.64188496988670089e-1;
    const double c2 = 8.88314979438837594e00;
    const double c3 = 6.61191906371416295e01;
    const double c4 = 2.98635138197400131e02;
    const double c5 = 8.81952221241769090e02;
    const double c6 = 1.71204761263407058e03;
    const double c7 = 2.05107837782607147e03;
    const double c8 = 1.23033935479799725e03;

    const double d0 = 1.00000000000000000e00;
    const double d1 = 1.57449261107098347e01;
    const double d2 = 1.17693950891312499e02;
    const double d3 = 5.37181101862009858e02;
    const double d4 = 1.62138957456669019e03;
    const double d5 = 3.29079923573345963e03;
    const double d6 = 4.36261909014324716e03;
    const double d7 = 3.43936767414372164e03;
    const double d8 = 1.23033935480374942e03;

    const double p0 = 1.63153871373020978e-2;
    const double p1 = 3.05326634961232344e-1;
    const double p2 = 3.60344899949804439e-1;
    const double p3 = 1.25781726111229246e-1;
    const double p4 = 1.60837851487422766e-2;
    const double p5 = 6.58749161529837803e-4;

    const double q0 = 1.00000000000000000e00;
    const double q1 = 2.56852019228982242e00;
    const double q2 = 1.87295284992346047e00;
    const double q3 = 5.27905102951428412e-1;
    const double q4 = 6.05183413124413191e-2;
    const double q5 = 2.33520497626869185e-3;

    register double y, z;


    y = fabs(u);
    if (y <= 0.46875 * SQRT2) 
    {
        /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
        z = y * y;
        y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4)
              / ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);

        return (0.5 + y);
    }
 
    z = exp(- y*y/2) / 2;
    if (y <= 4.0 * SQRT2) 
    {
        /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
        y = y / SQRT2;
        y = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8)
          / ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

        y = z * y;
    } 
    else 
    {
        /* evaluate erfc() for |u| > sqrt(2)*4.0 */
        z = z * SQRT2 / y;
        y = 2 / (y * y);
        y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5)
              / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
        
        y = z * (1./SQRTPI - y);
    }

    return (u < 0.0 ? y : 1 - y);

}


/*****  Normal_InvHP  ********************************************************/
/*                                                                           */
/*      Inverse of normal cumulative distribution                            */
/*      This function calculates the normal deviate Z corresponding to a     */
/*      given lower tail area of p. Z is accurate to about 1 part in 10**16. */
/*                                                                           */
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, PP 477-484    */
/*                                                                           */
double  Normal_InvHP (double p)
{
    const double a0 = 3.3871328727963666080e0000;
    const double a1 = 1.3314166789178437745e0002;
    const double a2 = 1.9715909503065514427e0003;
    const double a3 = 1.3731693765509461125e0004;
    const double a4 = 4.5921953931549871457e0004;
    const double a5 = 6.7265770927008700853e0004;
    const double a6 = 3.3430575583588128105e0004;
    const double a7 = 2.5090809287301226727e0003;

    const double b1 = 4.2313330701600911252e0001;
    const double b2 = 6.8718700749205790830e0002;
    const double b3 = 5.3941960214247511077e0003;
    const double b4 = 2.1213794301586595867e0004;
    const double b5 = 3.9307895800092710610e0004;
    const double b6 = 2.8729085735721942674e0004;
    const double b7 = 5.2264952788528545610e0003;

    const double c0 = 1.42343711074968357734e000;
    const double c1 = 4.63033784615654529590e000;
    const double c2 = 5.76949722146069140550e000;
    const double c3 = 3.64784832476320460504e000;
    const double c4 = 1.27045825245236838258e000;
    const double c5 = 2.41780725177450611770e-01;
    const double c6 = 2.27238449892691845833e-02;
    const double c7 = 7.74545014278341407640e-04;

    const double d1 = 2.05319162663775882187e000;
    const double d2 = 1.67638483018380384940e000;
    const double d3 = 6.89767334985100004550e-01;
    const double d4 = 1.48103976427480074590e-01;
    const double d5 = 1.51986665636164571966e-02;
    const double d6 = 5.47593808499534494600e-04;
    const double d7 = 1.05075007164441684324e-09;

    const double e0 = 6.65790464350110377720e000;
    const double e1 = 5.46378491116411436990e000;
    const double e2 = 1.78482653991729133580e000;
    const double e3 = 2.96560571828504891230e-01;
    const double e4 = 2.65321895265761230930e-02;
    const double e5 = 1.24266094738807843860e-03;
    const double e6 = 2.71155556874348757815e-05;
    const double e7 = 2.01033439929228813265e-07;

    const double f1 = 5.99832206555887937690e-01;
    const double f2 = 1.36929880922735805310e-01;
    const double f3 = 1.48753612908506148525e-02;
    const double f4 = 7.86869131145613259100e-04;
    const double f5 = 1.84631831751005468180e-05;
    const double f6 = 1.42151175831644588870e-07;
    const double f7 = 2.04426310338993978564e-15;

    const double const1 = 0.180625;
    const double const2 = 1.6;

    const double split1 = 0.425;
    const double split2 = 5.0;

    register double out, q, r;


    if (p >  1.0 || p <  0.0) { return (-999); }           /* not accepted input values */

    q = p - 0.5;
    if (p >= 1.0 || p <= 0.0) { return ( (q > 0.0) ? 9.0 : -9.0); }  /* for p=1 and p=0 */  

    if (fabs(q) <= split1)
    {
        r = const1 - q*q;
        out = ( q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0 )
                  / (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + 1.0) );

        return (out);
    }
    else
    {
        r = (q < 0.0) ? p : 1.0 - p;
        r = sqrt(-log(r));
        if (r <= split2)
        {
            r -= const2;
            out = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0 )
                / (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + 1.0);
        }
        else
        {
            r -= split2;
            out = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0 )
                / (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + 1.0);           
        }
    }

    return (q < 0.0 ? -out : out);
}


/*****  Normal_InvH  ********************************************************/
/**
*      Normal cumulative inverse distribution Risk Magazine
*/
double Normal_InvH (double prob)
{

    /* Coefficients for inverse cumulative normal */
    const double   a0 =   2.50662823884;
    const double   a1 = -18.61500062529;
    const double   a2 =  41.39119773534;
    const double   a3 = -25.44106049637;
    const double   b0 =  -8.47351093090;
    const double   b1 =  23.08336743743;
    const double   b2 = -21.06224101826;
    const double   b3 =   3.13082909833;
    const double   c0 =   0.3374754822726147;
    const double   c1 =   0.9761690190917186;
    const double   c2 =   0.1607979714918209;
    const double   c3 =   0.0276438810333863;
    const double   c4 =   0.0038405729373609;
    const double   c5 =   0.0003951896511919;
    const double   c6 =   0.0000321767881768;
    const double   c7 =   0.0000002888167364;
    const double   c8 =   0.0000003960315187;

    double  t;
    double  x;
    double  r;
   
    t = (prob < 0.5) ? (1. - prob) : prob;
   
    x = t - 0.5;
    if (fabs (x) < 0.42)
    {
        r = x * x;
        r = x * (((a3 * r + a2) * r + a1) * r + a0) 
            / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.);
        return (prob < 0.5) ? -r : r;
    }
    else
    {
        r = t;
        if (x > 0.) r = 1. - t;
        r = log (- log (r));
        r = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r 
                                            * (c6 + r * (c7 + r * c8)))))));
        if (x < 0.) r = -r;
        return (prob < 0.5) ? -r : r;
   }
}



/*****  NDensity  ************************************************************/
/**
*       Normal density.
*/
double  NDensity (double x)
{
    return (exp (-0.5 * x * x) * INVSQRT2PI);
}

/*****  pdf  ****************************************************************/
/**
*      Density function of the normal distribution.
*/
double  pdf (   double  x)
{
        return (exp (-.5 * x * x) / sqrt(2. * PI));

}  /* pdf */



/*****  BS_Density  **********************************************************/
/**
*      Density function in a Black & Scholes environment.
*/
double  BS_Density (    double S,  /**< (I) Price of the underlying       */
                        double S0, /**< (I) Initial price or expectation  */
                        double T,  /**< (I) Option expiration in years    */
                        double s   /**< (I) Annualized volatility         */
        )
{
        double
                nd,                         /* Density                       */
                d,                          /* d in N(d) in Black & Scholes  */
                st;                         /* Sigma * sqrt(T)               */


        if (S <= TINY * S0)
                return (0.);

        st  = s * sqrt (T);

        d  = log (S/S0) / st + .5 * st;
        nd = pdf (d);

        return (nd);

}  /* BS_Density */

/*****  Binary_BS  **********************************************************/
/**
       Binary price using Black&Scholes.
*/
double  Binary_BS ( double  S,    /**< Price of the underlying         */
                    double  K,    /**< Strike                          */
                    double  T,    /**< Binary expiration in years      */
                    double  r,    /**< Risk free rate                  */
                    double  s,    /**< Annualized volatility           */
                    char    CoP)  /**< Call or put binary ('C' or 'P') */
{

    double  D;  /* Delta */
    double  d;  /* d in N(d) in Black & Scholes formula */


    if (fabs(S) < TINY)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d - d / 2.;
    D  = NormalH (d);
        
    if (CoP == 'C')
        return (D);
    else
        return (1. - D);

}  /* Binary_BS */



/*****  Call_BS  ************************************************************/
/**
       Price of a Call using Black&Scholes.
*/
double  Call_BS (   double S,  /**< Price of the underlying    */
                    double K,  /**< Strike                     */
                    double T,  /**< Option expiration in years */
                    double r,  /**< Risk free rate             */
                    double s)  /**< Annualized volatility      */
{

    double  C;     /* Call price */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */


    if (S <= 0.)
        return (0.);
    if (K <= 0.0)
        return S;

    if (K <= 0.0)
      return S;             

    K  *= exp (-r * T);
    st  = s * sqrt (T);

    if (st <= 0.0)
        return MAX(S-K, 0.0);

    d   = log (S / K) / st;
    st /= 2.;

    C   = S * NormalH (d + st) - K * NormalH (d - st);

    return (C);

}  /* Call_BS */



/*****  Put_BS  *************************************************************/
/**
       Price of a Put using Black&Scholes.
*/
double  Put_BS (double S,  /**< Price of the underlying    */
                double K,  /**< Strike                     */
                double T,  /**< Option expiration in years */
                double r,  /**< Risk free rate             */
                double s)  /**< Annualized volatility      */
{

    double  P;     /* Put price */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */


    K *= exp (-r * T);

    if(S <= 0.)
        return (K);

    

    st  = s * sqrt (T);

    if (st <= 0.0)
        return MAX(K-S, 0.0);

    d   = log (K / S) / st;
    st /= 2.;
    P   = K * NormalH (d + st) - S * NormalH (d - st);

    return (P);

}  /* Put_BS */



/*****  D_Call_BS  **********************************************************/
/**
       Delta of a Call using Black&Scholes.
*/
double  D_Call_BS ( double S,   /**< Price of the underlying    */
                    double K,   /**< Strike                     */
                    double T,   /**< Option expiration in years */
                    double r,   /**< Risk free rate             */
                    double s)   /**< Annualized volatility      */
{

    double  D;   /* Delta */
    double  d;   /* d in N(d) in Black & Scholes formula */


    if (fabs(S) < TINY)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d + d / 2.;
    D  = NormalH (d);

    return (D);

}  /* D_Call_BS */


/*****  D_Put_BS  ***********************************************************/
/**
       Delta of a Put using Black&Scholes.
*/
double  D_Put_BS (  double S,  /**< Price of the underlying    */
                    double K,  /**< Strike                     */
                    double T,  /**< Option expiration in years */
                    double r,  /**< Risk free rate             */
                    double s)  /**< Annualized volatility      */
{
    return (D_Call_BS (S, K, T, r, s) - 1.);

}  /* D_Put_BS */



/*****  T_Call_BS  **********************************************************/
/**
       Theta of a Call using Black&Scholes.
*/
double  T_Call_BS ( double S, /**< Price of the underlying    */
                    double K, /**< Strike                     */
                    double T, /**< Option expiration in years */
                    double r, /**< Risk free rate             */
                    double s) /**< Annualized volatility      */
{

    double  Th;     /* Theta */
    double  d;      /* d in N(d) in Black & Scholes formula */
    double  st;     /* .5 * Sigma * sqrt(T) */


    if (fabs(S) < TINY)
        return (0.);

    K *= exp (-r * T);
    st  = s * sqrt (T);
    d   = log (S / K) / st;
    st /= 2.;

    Th  = -S * s / (2. * sqrt(T)) * exp (-(d + st) * (d + st) / 2.) 
        / sqrt(2. * PI) - K * r * NormalH (d - st);

    return (Th);

}  /* T_Call_BS */



/*****  T_Put_BS  ***********************************************************/
/**
       Theta of a Put using Black&Scholes.
*/
double  T_Put_BS (  double S,   /**< Price of the underlying    */
                    double K,   /**< Strike                     */
                    double T,   /**< Option expiration in years */
                    double r,   /**< Risk free rate             */
                    double s)   /**< Annualized volatility      */
{
    
    double  Th;     /* Theta */

    Th = T_Call_BS (S, K, T, r, s) + r * K * exp (-r * T);

    return (Th);

}  /* T_Put_BS */



/*****  R_Call_BS  **********************************************************/
/**
       Rho of a Call using Black&Scholes.
*/
double  R_Call_BS ( double S, /**< Price of the underlying    */
                    double K, /**< Strike                     */
                    double T, /**< Option expiration in years */
                    double r, /**< Risk free rate             */
                    double s) /**< Annualized volatility      */
{

    double  R;    /* rho */
    double  d;    /* d in N(d) in Black & Scholes formula */


    K *= exp (-r * T);
    if (fabs(S) < TINY)
        return (T * K);

    d  = s * sqrt (T);
    d  = log (S / K) / d - d /2.;
    R  = T * K * NormalH (d);

    return (R);

}  /* R_Call_BS */



/*****  R_Put_BS  ***********************************************************/
/**
       Rho of a Put using Black&Scholes.
*/
double  R_Put_BS (  double  S,  /**< Price of the underlying    */
                    double  K,  /**< Strike                     */
                    double  T,  /**< Option expiration in years */
                    double  r,  /**< Risk free rate             */
                    double  s)  /**< Annualized volatility      */
{
    double  R;        /* rho */

    R = R_Call_BS (S, K, T, r, s) - T * K * exp (-r * T);

    return (R);

}  /* R_Put_BS */


/*****  Call_BSQ  ***********************************************************/
/**
*      Price of a Call using Q version of Black&Scholes.
*/
double  Call_BSQ (double Y,  /**< (I) Fwd yield                     */
                  double K,  /**< (I) Strike                        */
                  double T,  /**< (I) Option expiration in years    */
                  double s,  /**< (I) Annualized volatility         */
                  double Q)  /**< (I) Q weight                      */
{

    double  C;                             /* Call price                    */
    double  d;                             /* d in N(d) in Black & Scholes  */
    double  st;                            /* Sigma * sqrt(T) * Q           */
    double  Y1, K1;                        /* Modified fwd and strike       */


    if (fabs(Y) < TINY) return (0.);

    if (fabs(Q) > QCUTOFF)
    {
        st  = s * sqrt (T) * Q;
        d   = - log ((K / Y - 1.) * Q + 1.) / st;
        st /= 2.;

        Y1  = Y / Q;
        K1  = K - Y + Y1;

        C   = Y1 * NormalH (d + st) - K1 * NormalH (d - st);
    }
    else
    {
        st  = s * sqrt (T);
        d   = - (K / Y - 1.) / st;

        C   = Y * st * exp (- d*d / 2) * INVSQRT2PI + (Y - K) * NormalH (d);
    }

    if(s<TINY)
    {
    C = MAX(Y-K,0);
    }

    return (C);

}  /* Call_BSQ */



/*****  Put_BSQ  ************************************************************/
/**
*      Price of a Put using Black&Scholes.
*/
double  Put_BSQ (double Y,   /**< (I) Fwd yield                     */
                 double K,   /**< (I) Strike                        */
                 double T,   /**< (I) Option expiration in years    */
                 double s,   /**< (I) Annualized volatility         */
                 double Q)   /**< (I) Q weight                      */
{

    double  P;                             /* Put price                     */
    double  d;                             /* d in N(d) in Black & Scholes  */
    double  st;                            /* Sigma * sqrt(T) * Q           */
    double  Y1, K1;                        /* Modified fwd and strike       */


    if (fabs(Y) < TINY) return (K);

    if (fabs(Q) > QCUTOFF)
    {
        st  = s * sqrt (T) * Q;
        d   = - log ((K / Y - 1.) * Q + 1.) / st;
        st /= 2.;

        Y1  = Y / Q;
        K1  = K - Y + Y1;

        P   = K1 * NormalH (- d + st) - Y1 * NormalH (- d - st);
    }
    else
    {
        st  = s * sqrt (T);
        d   = (K / Y - 1.) / st;

        P   = Y * st * exp (- d*d / 2) * INVSQRT2PI + (K - Y) * NormalH (d);
    }

    if(s<TINY)
    {
        P = MAX(K-Y,0);
    }

    return (P);

}  /* Put_BSQ */



/*****  Option_BSQ  ********************************************************/
/*
*      Price of a call or put using Q version of Black&Scholes.
*/
double  Option_BSQ (double Y,                /* Fwd yield                     */
                    double K,                /* Strike                        */
                    double T,                /* Option expiration in years    */
                    double s,                /* Annualized volatility         */
                    char   CoP,              /* Call or put                   */
                    double Q)                /* Q weight                      */
{
    if (CoP == 'C')
    {
        if (s < TINY)
        {
            return (Y - K);
        }
        else
        {
            return Call_BSQ(Y,K,T,s,Q);
        }
    }
    else
    {
        if (s < TINY)
        {
            return (K - Y);
        }
        else
        {
            return Put_BSQ(Y,K,T,s,Q);
        }
    }
}


/*****  ImpVol_BSQ  ********************************************************/
/*
*      Implied vol of a call or put using Q version of Black&Scholes.
*/
double  ImpVol_BSQ (double Y,                /* (I) Fwd yield              */
                    double K,                /* (I) Strike                 */
                    double T,                /* (I) Option expiration      */
                    double P,                /* (I) Price of option        */
                    char   CoP,              /* (I) Call or put            */
                    double Q,                /* (I) Q                      */
                    double VolGuess)         /* (I) Initial vol guess      */
{
    int
           found,
           count;

    const double VOLRATIO = 0.01;
    const double PREQERR  = 0.00000001;
    const double BAD_IMPL = -9.99;

    double 
           Pr, Pr2, dPr, scale,
           iVol;


    scale = Y * sqrt(T / 6.28);
    iVol  = VolGuess;
    count = 0;
    found = FALSE;

    if (((CoP == 'C') && (P > Y)) || 
        ((CoP == 'P') && (P > K)))
    {
        return (BAD_IMPL);
    }
    if (((CoP == 'C') && (P <= Y-K)) || 
        ((CoP == 'P') && (P <= K-Y)))
    {
        return (0.);
    }

    do 
    {
        Pr  = Option_BSQ (Y,K,T,iVol,CoP,Q);
        if (Pr < 0.0) 
        {
            return (BAD_IMPL);
        }

        Pr2 = Option_BSQ (Y,K,T,iVol+iVol*VOLRATIO,CoP,Q);
        if (Pr2 < 0.0)
        {
            return (BAD_IMPL);
        }

        if ((fabs(Pr - P) > scale * PREQERR))
        {
            dPr = (Pr2 - Pr) / (iVol*VOLRATIO);
            if (fabs(dPr) < TINY)
            {
                DR_Error("ImpVol_BSQ: numerical derivative is 0!");
                return (BAD_IMPL);
            }

            iVol -= (Pr - P) / dPr;
        }
        else 
        {
            found = TRUE;
        }
        count++ ;
    }
    while ( (found != TRUE) && (count < 20) );
                  
    if (found != TRUE) 
    {
        DR_Error("ImpVol_BSQ: exceeded maimum nb iterations!");        
        iVol = BAD_IMPL;
    }

    return (iVol);
} /* ImpVol_BSQ */



/*****  Vega_BSQ  ***********************************************************/
/**
*      Vega of put/call Black&Scholes.
*/
double  Vega_BSQ (double Y,  /**< (I) Fwd yield                     */
                  double K,  /**< (I) Strike                        */
                  double T,  /**< (I) Option expiration in years    */
                  double s,  /**< (I) Annualized volatility         */
                  double Q)  /**< (I) Q weight                      */
{

    double  V;                             /* Vega                          */
    double  d;                             /* d in N(d) in Black & Scholes  */
    double  sqrtT;                         /* sqrt of T                     */
    double  st;                            /* Sigma * sqrt(T) * Q           */


    if (fabs(Y) < TINY) return (0);

    if (fabs(Q) > QCUTOFF)
    {
        sqrtT = sqrt (T);
        st    = s * sqrtT * Q;
        d     = - log ((K / Y - 1.) * Q + 1.) / st;
        st   /= 2.;
        d    += st;

        V     = sqrtT * Y * exp(- d * d / 2.) * INVSQRT2PI;
    }
    else
    {
        sqrtT = sqrt (T);
        st    = s * sqrtT;
        d     = - (K / Y - 1.) / st;

        V     = sqrtT * Y * exp (- d * d / 2.) * INVSQRT2PI;
    }

    return (V);

}  /* Vega_BSQ */


/*****  Gamma_BS  **********************************************************/
/**
*      Gamma of an option using Black&Scholes.
*      Remember that the Gamma of a Call is equal to the Gamma of a Put.
*/
double  Gamma_BS (  double        S,   /**< (I) Price of the underlying    */
                    double        K,   /**< (I) Strike                     */
                    double        T,   /**< (I) Option expiration in years */
                    double        r,   /**< (I) Risk free rate             */
                    double        s)   /**< (I) Annualized volatility      */
{

    double  G;    /* Gamma */
    double  d;    /* d in N(d) in Black & Scholes formula */


    if (fabs(S) < TINY)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d + d / 2.;

    G  = exp (-d*d / 2.) / (sqrt (2. * PI) * S * s * sqrt (T));

    return (G);

}  /* Gamma_BS */


/*****  Vega_BS  ************************************************************/
/**
*      Vega of an option using Black&Scholes.
*      Remember that the Vega of a Call is equal to the Vega of a Put.
*/
double  Vega_BS (   double         S,   /**< (I) Price of the underlying    */
                    double         K,   /**< (I) Strike                     */
                    double         T,   /**< (I) Option expiration in years */
                    double         r,   /**< (I) Risk free rate             */
                    double         s)   /**< (I) Annualized volatility      */
{

    double  V;    /* Vega */
    double  d;    /* d in N(d) in Black & Scholes formula */


    if (fabs(S) < TINY)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d + d / 2.;
    V  = S * sqrt (T) * exp (-d*d / 2.) / sqrt (2. * PI);

    return (V);

}  /* Vega_BS */


/*****  CCEquation_BS2Q  ****************************************************/
/*
*       Evalute calibration equation for calibration constant.
*/
double  CCEquation_BS2Q (double    S,      /**< (I) Annualized volatility  */
                         double    QLeft,  /**< (I) Q left                 */
                         double    QRight, /**< (I) Q right                */
                         double    FwdSh,  /**< (I) Fwd shift              */
                         double    CC)     /**< (I) Constant               */
{
    double
            CCEq,
            CCIS,
            SQR,
            SQL; 

    SQR  = S * QRight;
    SQL  = S * QLeft;
    CCIS = CC / S;
    CCEq = FwdSh;

    if (fabs (QRight) > QCUTOFF)
    {
        CCEq -= (exp (CCIS * SQR + 0.5 * SQR * SQR) * NormalH (+ CCIS + SQR)
                 - NormalH (+ CCIS) ) / QRight;
    }
    else
    {
        CCEq -= (CCIS * NormalH(+ CCIS) + NDensity(CCIS)) * S;  
    }

    if (fabs (QLeft) > QCUTOFF)
    {
        CCEq -= (exp (CCIS * SQL + 0.5 * SQL * SQL) * NormalH (- CCIS - SQL)
                 - NormalH (- CCIS) ) / QLeft;
    }
    else
    {
        CCEq -= (CCIS * NormalH(- CCIS) - NDensity(CCIS)) * S;  
    }

    return CCEq;
}



/*****  ConvexityC_BS2Q  ****************************************************/
/**
*      Calibrate constant in 2q formulas. 
*/
int     ConvexityC_BS2Q (double    S,      /**< (I) Annualized volatility  */
                         double    QLeft,  /**< (I) Q left                 */
                         double    QRight, /**< (I) Q right                */
                         double    FwdSh,  /**< (I) Fwd shift              */
                         double    *CC)    /**< (O) Constant               */
{
    int
           count = 0,
           found,
           status = FAILURE;      /* Error status = FAILURE initially */

    double CCDELTA = 0.0001;
    double CCEQERR = 0.00000001;  /* JBL 30/1/03. Fix CET calibration problem.*/ 

    double 
           QM,     
           CC0,
           CCEq,
           CCEq2,
           dCC;

    QM = (QLeft + QRight) * 0.5;
    if (fabs (QM) > QCUTOFF)
    {
        CC0 = log (FwdSh * QM + 1.) / QM - 0.5 * QM * S * S;
    }
    else
    {
        CC0 = FwdSh;
    }

    found = FALSE;
    do 
    {
        CCEq  = CCEquation_BS2Q (S, QLeft, QRight, FwdSh, CC0);
        CCEq2 = CCEquation_BS2Q (S, QLeft, QRight, FwdSh, CC0 + CCDELTA);
        
        if (fabs(CCEq) > CCEQERR)
        {

            dCC   = (CCEq2 - CCEq) / CCDELTA;
            if (fabs(dCC) < TINY) goto RETURN;

            CC0  -= CCEq / dCC;
        }
        else
        {
            found = TRUE;
        }
        count++ ;
    }
    while ( (found != TRUE) && (count < 10) );

    if (found == TRUE)
    {
        *CC = CC0;
        status = SUCCESS;
    }
    else
    {
    *CC= 1.0e10;
    }

    RETURN:

    return (status);

}  /* ConvexityC_BS2Q */



/*****  Option_BS2Q  ********************************************************/
/**
*      Price of a call or put using 2Q version of Black&Scholes.
*      Note: vol is defined by Y = Y/(1+fsh) * F((1+fsh)/(1+q*fsh)*vol)
*/
double  Option_BS2Q (double Y,      /**< (I) Fwd yield              */
                     double K,      /**< (I) Strike                 */
                     double T,      /**< (I) Option expiration      */
                     double S,      /**< (I) Annualized volatility  */
                     char   CoP,    /**< (I) Call or put            */
                     double QLeft,  /**< (I) Q left                 */
                     double QRight, /**< (I) Q right                */
                     double FwdSh)  /**< (I) Fwd shift              */
{
    double
            QM,                            /* Q for fwd shift correction    */
            Q,                             /* Q used in part of density     */
            CC,                            /* Calibration const             */
            YCorr,                         /* Yield correction              */
            C, P,                          /* Call, put price               */
            d,                             /* d in N(d) in Black & Scholes  */
            St,                            /* Sigma * sqrt(T) * Q           */
            Y1, K1;                        /* Modified fwd and strike       */

    double  Ymin = -BIG;
    double  Ymax =  BIG;

    /* Non-singular FwdSh already checked in Param_Check */
    YCorr = Y / (1. + FwdSh);

    /* Set Ymin/max if distribution is bounded */
    if (QLeft  > 0.) Ymin = YCorr * (1. - 1./QLeft);
    if (QRight < 0.) Ymax = YCorr * (1. - 1./QRight);

    /* Deal with degenerate cases first. If strike is outside bounds */
    /* we define call and put as follows to avoid ill-defined log(K) */

    if (fabs(Y) < TINY)
    {
        C = MAX(-K, 0.);
        P = MAX( K, 0.);
        goto RETURN;
    }
    else if (K < Ymin + TINY)
    {
        P = 0.;
        C = Y-K;
        goto RETURN;
    }
    else if (K > Ymax - TINY)
    {
        C = 0.;
        P = K-Y;
        goto RETURN;
    }

    QM  = (QLeft + QRight) / 2;
    St  = S * sqrt (T) * (1 + FwdSh) / (1 + QM * FwdSh);

    if (ConvexityC_BS2Q (St, QLeft, QRight, FwdSh, &CC) == FAILURE)
    {
        DR_Error ("Could not find calib const for Call_BS2Q !");
        C = P = BS2Q_FAILURE_VALUE;
        goto RETURN;
    }

    Q = (K > YCorr) ? QRight : QLeft;
    
    if (fabs(Q) > QCUTOFF)
    {
        d   = (CC - log ((K / YCorr - 1.) * Q + 1.) / Q) / St;

        Y1  = exp (CC * Q + 0.5 * St * St * Q * Q);
        Y1 *= YCorr / Q;
        K1  = K - YCorr * (1. - 1. / Q);

        if (K > YCorr)
        {
            C = + Y1 * NormalH (+ d + St * Q) - K1 * NormalH (+ d);
            P = C - Y + K;
        }
        else
        {
            P = - Y1 * NormalH (- d - St * Q) + K1 * NormalH (- d);
            C = Y - K + P;
        }
    }
    else
    {
        d   = (CC - (K / YCorr - 1.)) / St;

        K1  = K - YCorr * (1. + CC);

        if (K > YCorr)
        {
            C = + YCorr * St * NDensity(d) - K1 * NormalH (+ d);  
            P = C - Y + K;
        }
        else
        {
            P = + YCorr * St * NDensity(d) + K1 * NormalH (- d);  
            C = Y - K + P;
        }
    }

    RETURN:

    return ((CoP == 'C') ? C : P);

}  /* Opt_BS2Q */


/*****  ImpVol_BS2QD  ********************************************************/
/**
        Derive the Implied Vol from a 2Q option price using the dichotomy method.

        Performs the same function as ImpVol_BS2Q, except it is more robust than
        the Newton-Raphson used there.
 */

int     ImpVol_BS2QD (double Y,       /**< (I) Fwd yield          */
                     double K,        /**< (I) Strike             */
                     double T,        /**< (I) Option expiration  */
                     double P,        /**< (I) Price of option    */
                     char   CoP,      /**< (I) Call or put        */
                     double QLeft,    /**< (I) Q left             */
                     double QRight,   /**< (I) Q right            */
                     double FwdSh,    /**< (I) Fwd shift          */
                     double VolGuess, /**< (I) Initial vol guess  */
                     double *ImpVol)  /**< (I) Initial vol guess  */
{
    int    i,
           found,
           count,
           status;

    const double VOLDELTA = 0.001;
    const double PREQERR  = 0.00001;

    double  scale,
            rP,
            lP,
            mP,
            res,
           iVol;

    status = FAILURE;
    scale = Y * sqrt(T / 6.28);
    iVol  = VolGuess;
    count = 0;
    found = FALSE;


    /*Middle point*/
    mP = VolGuess;
    res = Option_BS2Q (Y,K,T,VolGuess,CoP,QLeft,QRight,FwdSh);

    if(fabs(P - res)<PREQERR){
        *ImpVol = (VolGuess);
        status = SUCCESS;
        goto RETURN;
    }


    /*Find  right or left point */
    if(res<P){
        lP = VolGuess;
        rP = VolGuess;
        while(res<P){
            rP += VOLDELTA;
            res = Option_BS2Q (Y,K,T,rP,CoP,QLeft,QRight,FwdSh);
            if(rP >2){
                 goto RETURN;
            }
        }

    }
    else{
        
        lP = VolGuess;
        rP = VolGuess;
        while(res>=P){
            lP -= VOLDELTA;
            res = Option_BS2Q (Y,K,T,lP,CoP,QLeft,QRight,FwdSh);
            if(lP <0){
                 goto RETURN;
            }
        }

    }

    mP = (rP + lP)/(double)(2);
    res = Option_BS2Q (Y,K,T,mP,CoP,QLeft,QRight,FwdSh);
    i = 0;

    while(fabs(res-P)>PREQERR){
    
        if( (res> P)){
            /*Move left*/
            rP = mP;
            mP = (rP + lP)/(double)(2);
        }
        else{
        /*Move right*/
            lP = mP;
            mP = (rP + lP)/(double)(2);
        }

            res = Option_BS2Q (Y,K,T,mP,CoP,QLeft,QRight,FwdSh);
            i++;
            if(i >1000000){
                goto RETURN;
            }
    }

    *ImpVol = mP;

    status = SUCCESS;

RETURN:

    return status;
}            



/*****  ImpVol_BS2Q  ********************************************************/
/**
*      Implied vol of a call or put using 2Q version of Black&Scholes.
*/
double  ImpVol_BS2Q (double Y,         /**< (I) Fwd yield              */
                     double K,         /**< (I) Strike                 */
                     double T,         /**< (I) Option expiration      */
                     double P,         /**< (I) Price of option        */
                     char   CoP,       /**< (I) Call or put            */
                     double QLeft,     /**< (I) Q left                 */
                     double QRight,    /**< (I) Q right                */
                     double FwdSh,     /**< (I) Fwd shift              */
                     double VolGuess)  /**< (I) Initial vol guess      */
{
    int
           found,
           count;

    const double VOLDELTA = 0.0001;
    const double PREQERR  = 0.00000001;

    double 
           Pr, Pr2, dPr, scale,
           iVol;


    scale = Y * sqrt(T / 6.28);
    iVol  = VolGuess;
    count = 0;
    found = FALSE;

    do 
    {
        Pr  = Option_BS2Q (Y,K,T,iVol,CoP,QLeft,QRight,FwdSh);
        if (Pr < 0.0) 
        {
            DR_Error("ImpVol_BS2Q: Option_BS2Q price is < 0!");
            return (BS2Q_FAILURE_VALUE);
        }

        Pr2 = Option_BS2Q (Y,K,T,iVol+VOLDELTA,CoP,QLeft,QRight,FwdSh);
        if (Pr2 < 0.0)
        {
            DR_Error("ImpVol_BS2Q: tweaked Option_BS2Q price is < 0!");
            return (BS2Q_FAILURE_VALUE);
        }

        if ((fabs(Pr - P) > scale * PREQERR))
        {
            dPr = (Pr2 - Pr) / VOLDELTA;
            if (fabs(dPr) < TINY)
            {
                DR_Error("ImpVol_BS2Q: numerical derivative is 0!");
                return (BS2Q_FAILURE_VALUE);
            }

            iVol -= (Pr - P) / dPr;
        }
        else 
        {
            found = TRUE;
        }
        count++ ;
    }
    while ( (found != TRUE) && (count < 10) );
                  
    /* failure if exceed 10 iterations */
    if (found != TRUE) 
    {
        DR_Error("ImpVol_BS2Q: exceeded maimum nb iterations!");        
        iVol = BS2Q_FAILURE_VALUE;
    }

    return (iVol);


}


                  
/*****  NormalAS  *************************************************************/
/**
*       Cumulative function of the normal distribution.
*       See Handbook of Mathematical Functions 
*       by Milton Abramowitz and Irene A. Stegun
*       26.2.19 page 932                       
*/

double  NormalAS (double  x)
{
    const double d1 =  .0498673470; 
    const double d2 =  .0211410061; 
    const double d3 =  .0032776263;
    const double d4 =  .0000380036;
    const double d5 =  .0000488906;
    const double d6 =  .0000053830;

        double  y;
        double  z;
                
        int     i;        


        if (x >= BIG)          /*  x=+Œ => N(x)=1 */
                return (1.);
        else if (x <= -BIG)
                return (0.);

        y = fabs (x);
        y = 1. + y * (d1 + y * (d2 + y * (d3 + y * (d4 + y * (d5 + y * d6)))));
       


        for (i = 1, z = y; i < 16; i++)      /* For efficiency */
                z *= y;


        if (x >= 0.)
                return (1. - .5 / z);
        else
                return (.5 / z);

}  /* NormalAS */


/*****  KO_ExpAS  **************************************************************/
/**
*   Get the conditional knock-out probability, expectation and variance.
*/
void    KO_ExpAS (
                double  *P,  /**< Output: probability of not knocking-out   */
                double  *ES, /**< Output: expected value cond on not K-out  */
                double  *ES2,/**< Output: expected sqr value cond on not KO */
                double  S,   /**< Underlying expectation                    */
                double  S0,  /**< Initial price                             */
                double  H,   /**< Barrier                                   */
                char    UoD, /**< Up & out or down & out ('U' or 'D')       */
                double  T,   /**< Option expiration in years                */
                double  s)   /**< Annualized volatility                     */
{
        double
                Nd1, Nd2, Nd3, d, x,    /* Intermediate doubles              */
                st,                     /* Sigma * sqrt(T)                   */
                mut;                    /* Lognormal drift of the underlying */
        int
                Phi;                    /* +1 for down & out,-1 for up & out */


        Phi = (UoD == 'D') - (UoD == 'U');

        if ((H <= ERROR) && (UoD == 'D'))      /* If barrier=0 we have no KO */ 
        {                                      /* probability for down & out */
                *P   = 1.;
                *ES  = S;
                *ES2 = S * S * exp (s * s * T);

        }  /* if */
        
        if (Phi * S0 < Phi * H)   /* Already knocked-out */
        {        
                *P   = 0.;
                *ES  = 0.;
                *ES2 = 0.;
                
        }  /* if */

                
        st  = s * sqrt (T);
        mut = log (S / S0) - .5 * st * st;

        d   = Phi * (log (H/S0) + mut) / st;
        Nd1 = -NormalAS (d);
        d   += Phi * st;
        Nd2 = -NormalAS (d);
        d   += Phi * st;
        Nd3 = -NormalAS (d);
        
        x   = pow (H/S0, 2.*mut/st/st);
        Nd1 *= x;
        x   *= H*H/S0/S0;
        Nd2 *= x;
        x   *= H*H/S0/S0;
        Nd3 *= x;
        
        d   = Phi * (log (S0/H) + mut) / st;
        Nd1 += NormalAS (d);
        d   += Phi * st;
        Nd2 += NormalAS (d);
        d   += Phi * st;
        Nd3 += NormalAS (d);
        
        Nd2 *= S;
        Nd3 *= S * S * exp (st * st);
        
        *P   = Nd1;
        *ES  = Nd2;
        *ES2 = Nd3;
        
        return;
        
}  /* KO_ExpAS */



/*****  Binary_BSAS  *********************************************************/
/**
*      Binary price using Black&Scholes.
*/
double  Binary_BSAS (
                   double       S,   /**< Price of the underlying         */
                   double       K,   /**< Strike                          */
                   double       T,   /**< Binary expiration in years      */
                   double       r,   /**< Risk free rate                  */
                   double       s,   /**< Annualized volatility           */
                   char         CoP) /**< Call or put binary ('C' or 'P') */    
{
        double
                D,            /* Delta                                */
                d;            /* d in N(d) in Black & Scholes formula */


        if (fabs(S) < TINY) return (0.);

        K *= exp (-r * T);
        d  = s * sqrt (T);
        d  = log (S / K) / d - d / 2.;
        D  = NormalAS (d);
        
        if (CoP == 'C')
                return (D);
        else
                return (1. - D);

}  /* Binary_BSAS */



/*****  Call_BSAS ************************************************************/
/**
*      Price of a Call using Black&Scholes.
*/
double  Call_BSAS (
         double S, /**< Price of the underlying    */
                 double K, /**< Strike                     */
                 double T, /**< Option expiration in years */
                 double r, /**< Risk free rate             */
                 double s) /**< Annualized volatility      */
{
        double
                C,                /*  Call price                             */
                d,                /*  d in N(d) in Black & Scholes formula   */
                st;               /*  Sigma * sqrt(T)                        */


        if (fabs(S) < TINY) return (0.);

        K  *= exp (-r * T);
        st  = s * sqrt (T);
        d   = log (S / K) / st;
        st /= 2.;

        C   = S * NormalAS (d + st) - K * NormalAS (d - st);

        return (C);

}  /* Call_BSAS */



/*****  Put_BSAS  ************************************************************/
/**
*      Price of a Put using Black&Scholes.
*/
double  Put_BSAS (
        double S, /**< Price of the underlying    */
                double K, /**< Strike                     */
                double T, /**< Option expiration in years */
                double r, /**< Risk free rate             */
                double s) /**< Annualized volatility      */
{
        double
                P,                 /*  Put price                             */
                d,                 /*  d in N(d) in Black & Scholes formula  */
                st;                /*  Sigma * sqrt(T)                       */

        K *= exp (-r * T);

        if(fabs(S) < TINY)
                return (K);

        st  = s * sqrt (T);
        d   = log (K / S) / st;
        st /= 2.;
        P   = K * NormalAS (d + st) - S * NormalAS (d - st);

        return (P);

}  /* Put_BSAS */



/*****  D_Call_BSAS  *********************************************************/
/**
*      Delta of a Call using Black&Scholes.
*/
double  D_Call_BSAS (
           double  S, /**< Price of the underlying    */
                   double  K, /**< Strike                     */
                   double  T, /**< Option expiration in years */
                   double  r, /**< Risk free rate             */
                   double  s) /**< Annualized volatility      */
{
        double
                D,                  /* Delta                                  */
                d;                  /* d in N(d) in Black & Scholes formula   */


        if (fabs(S) < TINY)
                return (0.);

        K *= exp (-r * T);
        d  = s * sqrt (T);
        d  = log (S / K) / d + d / 2.;
        D  = NormalAS (d);

        return (D);

}  /* D_Call_BSAS */



/*****  D_Put_BS  ************************************************************/
/**
*      Delta of a Put using Black&Scholes.
*/
double  D_Put_BSAS (
          double  S, /**< Price of the underlying    */
                  double  K, /**< Strike                     */
                  double  T, /**< Option expiration in years */
                  double  r, /**< Risk free rate             */
                  double  s) /**< Annualized volatility      */
{
        return (D_Call_BSAS (S, K, T, r, s) - 1.);

}  /* D_Put_BS */


/*****  T_Call_BSAS  *********************************************************/
/**
*      Theta of a Call using Black&Scholes.
*/
double  T_Call_BSAS (
           double S, /**< Price of the underlying    */
                   double K, /**< Strike                     */
                   double T, /**< Option expiration in years */
                   double r, /**< Risk free rate             */
                   double s) /**< Annualized volatility      */
{
        double
                Th,  /* Theta                                */
                d,   /* d in N(d) in Black & Scholes formula */
                st;  /* .5 * Sigma * sqrt(T)                 */


        if (fabs(S) < TINY)
                return (0.);

        K *= exp (-r * T);
        st  = s * sqrt (T);
        d   = log (S / K) / st;
        st /= 2.;

        Th  = -S * s / (2. * sqrt(T)) * exp (-(d + st) * (d + st) / 2.) 
                       / sqrt(2. * PI) - K * r * NormalAS (d - st);

        return (Th);

}  /* T_Call_BSAS */



/*****  T_Put_BSAS  **********************************************************/
/**
*      Theta of a Put using Black&Scholes.
*/
double  T_Put_BSAS (
          double S, /**< Price of the underlying    */
                  double K, /**< Strike                     */
                  double T, /**< Option expiration in years */
                  double r, /**< Risk free rate             */
                  double s) /**< Annualized volatility      */
{
        double
                Th;  /* Theta */

        Th = T_Call_BSAS (S, K, T, r, s) + r * K * exp (-r * T);

        return (Th);

}  /* T_Put_BSAS */



/*****  R_Call_BSAS  *********************************************************/
/**
*      Rho of a Call using Black&Scholes.
*/
double  R_Call_BSAS (
           double S, /**< Price of the underlying    */
                   double K, /**< Strike                     */
                   double T, /**< Option expiration in years */
                   double r, /**< Risk free rate             */
                   double s) /**< Annualized volatility      */
{
        double
                R,      /* rho */
                d;      /* d in N(d) in Black & Scholes formula */


        K *= exp (-r * T);
        if (fabs(S) < TINY)
                return (T * K);

        d  = s * sqrt (T);
        d  = log (S / K) / d - d /2.;
        R  = T * K * NormalAS (d);

        return (R);

}  /* R_Call_BSAS */



/*****  R_Put_BSAS  **********************************************************/
/**
*      Rho of a Put using Black&Scholes.
*/
double  R_Put_BSAS (
          double S,  /**< Price of the underlying    */
                  double K,  /**< Strike                     */
                  double T,  /**< Option expiration in years */
                  double r,  /**< Risk free rate             */
                  double s)  /**< Annualized volatility      */
{
        double
                R; /* rho */

        R = R_Call_BSAS (S, K, T, r, s) - T * K * exp (-r * T);

        return (R);

}  /* R_Put_BSAS */



/*****  Option_Normal  *******************************************************/
/**
       Price of an option using Normal Model
*/
double  Option_Normal ( double  S,  /**< Price of the underlying         */
                        double  K,  /**< Strike                          */
                        char    CoP,/**< Call or put                     */
                        double  T,  /**< Option expiration in years      */
                        double  r,  /**< (I) Risk free rate              */
                        double  s)  /**< Annualized normal volatility    */
{
    double  C;     /* call price */
    double  P;     /* Put price  */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */

    st  = s * sqrt (T);
    d   = (K - S) / st;

    P   = (K - S) * NormalH (d) + st * NDensity (d);
    C   = S - K + P;
    
    if (CoP == 'C')
        return (C);
    else
        return (P);

}  /* Option_Normal */



/*****  Vega_Normal  *********************************************************/
/**
*      Vega of an option using Normal Model
*      
*/
double  Vega_Normal (double         S,   /**< (I) Price of the underlying    */
                    double         K,   /**< (I) Strike                     */
                    double         T,   /**< (I) Option expiration in years */
                    double         r,   /**< (I) Risk free rate             */
                    double         s)   /**< (I) Annualized volatility      */
{

    double  V;    /* Vega */
    double  st;
    double  d;    

    st  = s * sqrt (T);
    d = (K - S) / st;
    V  = sqrt (T) * NDensity (d);

    return (V);

}  /* Vega_Normal */


/*****  Option_E2Q  ********************************************************/
/**
*      Price of a call or put using E2Q version of Black&Scholes.
*      
*/
double  Option_E2Q  (
            double          Y,              /**< (I) Fwd yield              */
            double          K,              /**< (I) Strike                 */
            double          T,              /**< (I) Option expiration      */
            double          S,              /**< (I) Annualized volatility  */
            char            CoP,            /**< (I) Call or put            */
            double          QLeft,          /**< (I) QLeft                  */
            double          QRight,         /**< (I) QRight                 */
            double          a,              /**< (I) A mapping parameter    */
            double          b,              /**< (I) B mapping parameter    */
            double          FwdShift)       /**< (I) forward shift          */
{
    double  price;

    price = Option_BS2Q (a+b*Y, K+a+(b-1)*Y, T, S, CoP, QLeft, QRight, FwdShift);
    
    return (price);

}  /* Option_E2Q */



/*****  ImpVol_E2Q  ********************************************************/
/**
*      Implied vol of a call or put using E2Q version of Black&Scholes.
*/
double  ImpVol_E2Q  (
            double          Y,              /**< (I) Fwd yield              */
            double          K,              /**< (I) Strike                 */
            double          T,              /**< (I) Option expiration      */
            double          P,              /**< (I) Price of option        */
            char            CoP,            /**< (I) Call or put            */
            double          QLeft,          /**< (I) QLeft                  */
            double          QRight,         /**< (I) QRight                 */
            double          a,              /**< (I) A mapping parameter    */
            double          b,              /**< (I) B mapping parameter    */
            double          FwdSh,          /**< (I) Forward shift          */
            double          VolGuess)       /**< (I) Initial vol guess      */
{
    double iVol;

    iVol = ImpVol_BS2Q(a+b*Y, K+a+(b-1)*Y, T, P, CoP, QLeft, QRight, FwdSh, VolGuess);

    return (iVol);

} /* ImpVol_E2Q */



/*****  ImpMarketVol   ***************************************************************/
/*
*       
*/
double ImpMarketVol (
                double  pY,             /* (I) Fwd par yield                */
                double  K,              /* (I) strike                       */
                double  price,          /* (I) atm price                    */
                char    CoP,            /* (I) call or put                  */
                int     VolTypeFlag,    /* (I) Normal or Lognormal          */
                double  T,              /* (I) Time to expiry               */
                double  FwdSh,          /* (I) Forward shift                */
                double  VolGuess)
{
    double Vol;

    if (VolTypeFlag == 1)
    {
        /* use implied volE2q with q = 1, a=0, b =1 for lognormal case */
        Vol = ImpVol_E2Q(pY,K,T,price,CoP, 1., 1., 0., 1., FwdSh, VolGuess);
    }
    else
    {
        /* use implied volE2q with q =0, a=1, b =0 for normal case */
        Vol = ImpVol_E2Q (pY,K,T,price,CoP, 0., 0., 1., 0.,FwdSh, VolGuess);
    }

    return (Vol);
} /* ImpMarketVol */


/* Market to X-space vol adjustment. The single and     */
/* two q cases are treated differently for consistency  */
/* with old model.                                      */
int Conv_MarketVol_To_XSpace (  
      double         pY,          /* (I) Fwd par yield                */
      double const   Vol,         /* (I) Market Vol                   */
      long const     VolDate,     /* (I) Volatility date              */
      int            VolTypeFlag, /* (I) Normal or Lognormal          */
      double*        XVol,        /* (O) Vol in X space               */
      double          T,          /* (I) Time to expiry               */
      double          FwdSh,      /* (I) Forward shift                */
      double          QLeft,      /* (I) QLeft                        */
      double          QRight,     /* (I) QRight                       */
      double          a,          /* (I) A mapping parameter          */
      double          b)          /* (I) B mapping parameter          */
{
    
    double tempXVol;
    double  atmPr;             /* Market price of atm option       */
    int     status = FAILURE;  /* Error status = FAILURE initially */

    int    IS_Classic;         /* Check if classic 2q model        */

    IS_Classic = (IS_EQUAL(a, 0.0) && (IS_EQUAL(b, 1.0)));


 
    /* if VolTypeFlag = 0 use normal model for ATM price    */
    /* otherwise use lognormal model*/
    if (VolTypeFlag == 0)    atmPr = sqrt(T) * Vol / sqrt (2. * PI);
    else
    {
         if (IS_Classic) 
             atmPr= Option_BS2Q (pY,pY,T,Vol,'C',1.,1.,0.);    
         else
             atmPr = pY * (2 * NormalH (.5 * sqrt(T) * Vol) - 1.);
    }
    

    
    if ((fabs(QLeft - QRight) < TINY) && (fabs(FwdSh) < TINY) && 
        (VolTypeFlag == 1) && IS_Classic)
    {
        if (fabs(QLeft) > QCUTOFF)                                             
            tempXVol = Normal_InvH (QLeft * (NormalH (.5 * sqrt(T) * Vol) -.5) + .5) / (.5 * QLeft);
        else
            tempXVol = (2. * NormalH (.5 * sqrt(T) * Vol) - 1.) * sqrt (2. * PI);

        *XVol = SQUARE(tempXVol);
    }
    else
    {
        tempXVol = ImpVol_E2Q (pY, pY, T, atmPr, 'C', QLeft, QRight, a, b, FwdSh,Vol);
        if (tempXVol < 0.0)
        {
            DR_Error ("Conv MarketVol To XSpace: problem in E2Q implied vol at %ld ", VolDate);
            goto RETURN;
        }
    
        *XVol = T * SQUARE (tempXVol);
    }
    
    status = SUCCESS;

    RETURN:

    return (status);

}/* Conv_MarketVol_To_XSpace */



/* Convert from q measure back to Market */
int Conv_XSpaceVol_To_Market (
       double       pY,          /* (I) Fwd par yield            */
       double*      Vol,         /* (O) Market Vol               */  
       int          VolTypeFlag, /* (I) Normal or Lognormal      */        
       double       XVol,        /* (I) Vol in X space           */
       double       T,           /* (I) Time to expiry           */
       double       FwdSh,       /* (I) Forward shift            */
       double       QLeft,       /* (I) QLeft                    */
       double       QRight,      /* (I) QRight                   */
       double       a,           /* (I) A mapping parameter      */
       double       b)           /* (I) B mapping parameter      */
{
    double atmPr;                      /* Model price of atm option        */
    int     status = FAILURE;          /* Error status = FAILURE initially */
    int    IS_Classic;                 /* Check if classic 2q model        */

    IS_Classic = (IS_EQUAL(a, 0.0) && (IS_EQUAL(b, 1.0)));
    

    /* compute model atm price*/
    if ((fabs(QLeft - QRight) < TINY) && (fabs(FwdSh) < TINY) && 
         (VolTypeFlag == 1) && IS_Classic)
    {
        XVol *=   sqrt (T);
        if (fabs(QLeft) > QCUTOFF)                                             
            XVol = 2. * Normal_InvH ((NormalH (.5 * QLeft * XVol) - .5) / QLeft + .5);
        else
            XVol = 2. * Normal_InvH (.5 * (1. + XVol / sqrt(2.*PI)));

        *Vol = XVol / sqrt(T);
    }
    else
    {
        atmPr = Option_E2Q (pY,pY,T,XVol,'C',QLeft, QRight, a, b, FwdSh);
        if (atmPr < 0.0)
        {
            DR_Error("Conv_XSpaceVol_To_Market: problem in E2Q price.");
            goto RETURN;
        }

        /* compute corresponding atm vol*/
        *Vol = ImpMarketVol (pY,pY, atmPr,'C',VolTypeFlag,T, 0.,XVol); 
    }

    status = SUCCESS;

    RETURN:

    return (status);

}/* Conv_XSpaceVol_To_Market */


