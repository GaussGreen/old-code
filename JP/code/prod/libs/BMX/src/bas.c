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
#include "bmx123head.h"


#define     QCUTOFF         1E-4           /* Normal model for |q|<QCUTOFF  */

#ifndef    BAD_IMPL
#define    BAD_IMPL  -9.99
#endif

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



/*****  NormalH  ************************************************************/
/*
*      Normal cumulative distribution accrding to J.Hull
*/
double  NormalH (double  x)
{
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



/*****  Normal_InvH  ********************************************************/
/*
*      Normal cumulative inverse distribution Risk Magazine
*/
double Normal_InvH (double prob)
{
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
/*
*       Normal density.
*/
double  NDensity (double x)
{
    return (exp (-0.5 * x * x) * INVSQRT2PI);
}



/*****  Call_BSQ  ***********************************************************/
/*
*      Price of a Call using Q version of Black&Scholes.
*/
double  Call_BSQ (double Y,                /* Fwd yield                     */
                  double K,                /* Strike                        */
                  double T,                /* Option expiration in years    */
                  double s,                /* Annualized volatility         */
                  double Q)                /* Q weight                      */
{

    double  C;	                           /* Call price                    */
    double  d;	                           /* d in N(d) in Black & Scholes  */
    double  st;	                           /* Sigma * sqrt(T) * Q           */
    double  Y1, K1;                        /* Modified fwd and strike       */


    if (Y == 0.) return (0.);

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

    return (C);

}  /* Call_BSQ */



/*****  Put_BSQ  ************************************************************/
/*
*      Price of a Put using Black&Scholes.
*/
double  Put_BSQ (double Y,                 /* Fwd yield                     */
                 double K,                 /* Strike                        */
                 double T,                 /* Option expiration in years    */
                 double s,                 /* Annualized volatility         */
                 double Q)                 /* Q weight                      */
{

    double  P;	                           /* Put price                     */
    double  d;	                           /* d in N(d) in Black & Scholes  */
    double  st;	                           /* Sigma * sqrt(T) * Q           */
    double  Y1, K1;                        /* Modified fwd and strike       */


    if (Y == 0.) return (K);

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

    return (P);

}  /* Put_BSQ */



/*****  Vega_BSQ  ***********************************************************/
/*
*      Vega of put/call Black&Scholes.
*/
double  Vega_BSQ (double Y,                /* Fwd yield                     */
                  double K,                /* Strike                        */
                  double T,                /* Option expiration in years    */
                  double s,                /* Annualized volatility         */
                  double Q)                /* Q weight                      */
{

    double  V;	                           /* Vega                          */
    double  d;                             /* d in N(d) in Black & Scholes  */
    double  sqrtT;                         /* sqrt of T                     */
    double  st;	                           /* Sigma * sqrt(T) * Q           */


    if (Y == 0.) return (0);

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



/*****  Binary_BS  **********************************************************/
/*
*      Binary price using Black&Scholes.
*/
double  Binary_BS (	double  S,    /* Price of the underlying         */
                    double  K,    /* Strike                          */
                    double  T,    /* Binary expiration in years      */
                    double  r,    /* Risk free rate                  */
                    double  s,    /* Annualized volatility           */
                    char    CoP)  /* Call or put binary ('C' or 'P') */
{

    double  D;  /* Delta */
    double  d;  /* d in N(d) in Black & Scholes formula */


    if (S == 0.)
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
/*
*      Price of a Call using Black&Scholes.
*/
double  Call_BS (	double S,  /* Price of the underlying    */
                    double K,  /* Strike                     */
                    double T,  /* Option expiration in years */
                    double r,  /* Risk free rate             */
                    double s)  /* Annualized volatility      */
{

    double  C;     /* Call price */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */


    if (S == 0.)
        return (0.);

    K  *= exp (-r * T);
    st  = s * sqrt (T);
    d   = log (S / K) / st;
    st /= 2.;

    C   = S * NormalH (d + st) - K * NormalH (d - st);

    return (C);

}  /* Call_BS */



/*****  Put_BS  *************************************************************/
/*
*      Price of a Put using Black&Scholes.
*/
double  Put_BS (double S,  /* Price of the underlying    */
                double K,  /* Strike                     */
                double T,  /* Option expiration in years */
                double r,  /* Risk free rate             */
                double s)  /* Annualized volatility      */
{

    double  P;     /* Put price */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */


    K *= exp (-r * T);

    if(S == 0.)
        return (K);

    st  = s * sqrt (T);
    d   = log (K / S) / st;
    st /= 2.;
    P   = K * NormalH (d + st) - S * NormalH (d - st);

    return (P);

}  /* Put_BS */



/*****  D_Call_BS  **********************************************************/
/*
*      Delta of a Call using Black&Scholes.
*/
double  D_Call_BS (	double       S,   /* Price of the underlying    */
                    double       K,   /* Strike                     */
                    double       T,   /* Option expiration in years */
                    double       r,   /* Risk free rate             */
                    double       s)   /* Annualized volatility      */
{

    double  D;   /* Delta */
    double  d;   /* d in N(d) in Black & Scholes formula */


    if (S == 0.)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d + d / 2.;
    D  = NormalH (d);

    return (D);

}  /* D_Call_BS */



/*****  D_Put_BS  ***********************************************************/
/*
*      Delta of a Put using Black&Scholes.
*/
double  D_Put_BS (  double        S,  /* Price of the underlying    */
                    double        K,  /* Strike                     */
                    double        T,  /* Option expiration in years */
                    double        r,  /* Risk free rate             */
                    double        s)  /* Annualized volatility      */
{
    return (D_Call_BS (S, K, T, r, s) - 1.);

}  /* D_Put_BS */



/*****  Gamma_BS  **********************************************************/
/*
*      Gamma of an option using Black&Scholes.
*      Remember that the Gamma of a Call is equal to the Gamma of a Put.
*/
double  Gamma_BS (  double        S,   /* Price of the underlying    */
                    double        K,   /* Strike                     */
                    double        T,   /* Option expiration in years */
                    double        r,   /* Risk free rate             */
                    double        s)   /* Annualized volatility      */
{

    double  G;    /* Gamma */
    double  d;    /* d in N(d) in Black & Scholes formula */


    if (S == 0.)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d + d / 2.;

    G  = exp (-d*d / 2.) / (sqrt (2. * PI) * S * s * sqrt (T));

    return (G);

}  /* Gamma_BS */



/*****  Vega_BS  ************************************************************/
/*
*      Vega of an option using Black&Scholes.
*      Remember that the Vega of a Call is equal to the Vega of a Put.
*/
double  Vega_BS (	double         S,   /* Price of the underlying    */
                    double         K,   /* Strike                     */
                    double         T,   /* Option expiration in years */
                    double         r,   /* Risk free rate             */
                    double         s)   /* Annualized volatility      */
{

    double  V;    /* Vega */
    double  d;    /* d in N(d) in Black & Scholes formula */


    if (S == 0.)
        return (0.);

    K *= exp (-r * T);
    d  = s * sqrt (T);
    d  = log (S / K) / d + d / 2.;
    V  = S * sqrt (T) * exp (-d*d / 2.) / sqrt (2. * PI);

    return (V);

}  /* Vega_BS */



/*****  T_Call_BS  **********************************************************/
/*
*      Theta of a Call using Black&Scholes.
*/
double  T_Call_BS (	double       S,     /* Price of the underlying    */
                    double       K,     /* Strike                     */
                    double       T,     /* Option expiration in years */
                    double       r,     /* Risk free rate             */
                    double       s)     /* Annualized volatility      */
{

    double  Th;     /* Theta */
    double  d;      /* d in N(d) in Black & Scholes formula */
    double  st;     /* .5 * Sigma * sqrt(T) */


    if (S == 0.)
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
/*
*      Theta of a Put using Black&Scholes.
*/
double  T_Put_BS (  double        S,   /* Price of the underlying    */
                    double        K,   /* Strike                     */
                    double        T,   /* Option expiration in years */
                    double        r,   /* Risk free rate             */
                    double        s)   /* Annualized volatility      */
{
    
    double  Th;     /* Theta */

    Th = T_Call_BS (S, K, T, r, s) + r * K * exp (-r * T);

    return (Th);

}  /* T_Put_BS */



/*****  R_Call_BS  **********************************************************/
/*
*      Rho of a Call using Black&Scholes.
*/
double  R_Call_BS (	double       S,     /* Price of the underlying    */
                    double       K,     /* Strike                     */
                    double       T,     /* Option expiration in years */
                    double       r,     /* Risk free rate             */
                    double       s)     /* Annualized volatility      */
{

    double  R;    /* rho */
    double  d;    /* d in N(d) in Black & Scholes formula */


    K *= exp (-r * T);
    if (S == 0.)
        return (T * K);

    d  = s * sqrt (T);
    d  = log (S / K) / d - d /2.;
    R  = T * K * NormalH (d);

    return (R);

}  /* R_Call_BS */



/*****  R_Put_BS  ***********************************************************/
/*
*      Rho of a Put using Black&Scholes.
*/
double  R_Put_BS (  double        S,      /* Price of the underlying    */
                    double        K,      /* Strike                     */
                    double        T,      /* Option expiration in years */
                    double        r,      /* Risk free rate             */
                    double        s)      /* Annualized volatility      */
{

    double  R;        /* rho */


    R = R_Call_BS (S, K, T, r, s) - T * K * exp (-r * T);

    return (R);

}  /* R_Put_BS */



/*****  CCEquation_BS2Q  ****************************************************/
/*
*       Evalute calibration equation for calibration constant.
*/
double  CCEquation_BS2Q (double    S,         /* (I) Annualized volatility  */
                         double    QLeft,     /* (I) Q left                 */
                         double    QRight,    /* (I) Q right                */
                         double    FwdSh,     /* (I) Fwd shift              */
                         double    CC)        /* (I) Constant               */
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
/*
*      Calibrate constant in 2q formulas. 
*/
#define CCDELTA  0.0001
#define CCEQERR  0.000001

int     ConvexityC_BS2Q (double    S,         /* (I) Annualized volatility  */
                         double    QLeft,     /* (I) Q left                 */
                         double    QRight,    /* (I) Q right                */
                         double    FwdSh,     /* (I) Fwd shift              */
                         double    *CC)       /* (O) Constant               */
{
    int
           count = 0,
           found,
           status = FAILURE;      /* Error status = FAILURE initially */

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
    while ( (found != TRUE) && (count < 20) );

    if (found == TRUE)
    {
        *CC = CC0;
        status = SUCCESS;
    }

    RETURN:

    return (status);

}  /* ConvexityC_BS2Q */



/*****  Option_BS2Q  ********************************************************/
/*
*      Price of a call or put using 2Q version of Black&Scholes.
*      Note: vol is defined by Y = Y/(1+fsh) * F((1+fsh)/(1+q*fsh)*vol)
*/
double  Option_BS2Q (double Y,                /* (I) Fwd yield              */
                     double K,                /* (I) Strike                 */
                     double T,                /* (I) Option expiration      */
                     double S,                /* (I) Annualized volatility  */
                     char   CoP,              /* (I) Call or put            */
                     double QLeft,            /* (I) Q left                 */
                     double QRight,           /* (I) Q right                */
                     double FwdSh)            /* (I) Fwd shift              */
{
    double
            QM,                            /* Q for fwd shift correction    */
            Q,                             /* Q used in part of density     */
            CC,                            /* Calibration const             */
            YCorr,                         /* Yield correction              */
            C, P,	                       /* Call, put price               */
            d,                             /* d in N(d) in Black & Scholes  */
            St,	                           /* Sigma * sqrt(T) * Q           */
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
        C = P = BAD_IMPL;
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



/*****  ImpVol_BS2Q  ********************************************************/
/*
*      Implied vol of a call or put using 2Q version of Black&Scholes.
*/
#define VOLDELTA  0.0001
#define PREQERR  0.00000001

double  ImpVol_BS2Q (double Y,                /* (I) Fwd yield              */
                     double K,                /* (I) Strike                 */
                     double T,                /* (I) Option expiration      */
                     double P,                /* (I) Price of option        */
                     char   CoP,              /* (I) Call or put            */
                     double QLeft,            /* (I) Q left                 */
                     double QRight,           /* (I) Q right                */
                     double FwdSh,            /* (I) Fwd shift              */
                     double VolGuess)         /* (I) Initial vol guess      */
{
    int
           found,
           count;

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
            return (BAD_IMPL);
        }

        Pr2 = Option_BS2Q (Y,K,T,iVol+VOLDELTA,CoP,QLeft,QRight,FwdSh);
        if (Pr2 < 0.0)
        {
            DR_Error("ImpVol_BS2Q: tweaked Option_BS2Q price is < 0!");
            return (BAD_IMPL);
        }

        if ((fabs(Pr - P) > scale * PREQERR))
        {
            dPr = (Pr2 - Pr) / VOLDELTA;
            if (fabs(dPr) < TINY)
            {
                DR_Error("ImpVol_BS2Q: numerical derivative is 0!");
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
                  
    /* failure if exceed 20 iterations */
    if (found != TRUE) 
    {
        DR_Error("ImpVol_BS2Q: exceeded maimum nb iterations!");        
        iVol = BAD_IMPL;
    }

    return (iVol);
}
                  

                  
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

    return 0;
}


/*****  ImpVol_BSQ  ********************************************************/
/*
*      Implied vol of a call or put using Q version of Black&Scholes.
*/
#define VOLRATIO  0.01
#define PREQERR  0.00000001

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
}
