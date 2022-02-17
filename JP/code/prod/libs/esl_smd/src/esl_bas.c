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
double  BS_Density (	double S,  /**< (I) Price of the underlying       */
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

    K  *= exp (-r * T);
    st  = s * sqrt (T);
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

    double  C;	                           /* Call price                    */
    double  d;	                           /* d in N(d) in Black & Scholes  */
    double  st;	                           /* Sigma * sqrt(T) * Q           */
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

    double  P;	                           /* Put price                     */
    double  d;	                           /* d in N(d) in Black & Scholes  */
    double  st;	                           /* Sigma * sqrt(T) * Q           */
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

    double  V;	                           /* Vega                          */
    double  d;                             /* d in N(d) in Black & Scholes  */
    double  sqrtT;                         /* sqrt of T                     */
    double  st;	                           /* Sigma * sqrt(T) * Q           */


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


        if (x >= BIG)          /*  x=+î => N(x)=1 */
                return (1.);
        else if (x <= -BIG)
                return (0.);

        y = fabs (x);
        y = 1. + y * (d1 + y * (d2 + y * (d3 + y * (d4 + y * (d5 + y * d6)))));
       


        for (i = 1, z = y; i < 16; i++)	     /* For efficiency */
                z *= y;


        if (x >= 0.)
                return (1. - .5 / z);
        else
                return (.5 / z);

}  /* NormalAS */


/*****  KO_ExpAS  **************************************************************/
/**
*	Get the conditional knock-out probability, expectation and variance.
*/
void	KO_ExpAS (
		double 	*P,  /**< Output: probability of not knocking-out   */
                double 	*ES, /**< Output: expected value cond on not K-out  */
                double	*ES2,/**< Output: expected sqr value cond on not KO */
                double 	S,   /**< Underlying expectation                    */
                double 	S0,  /**< Initial price                             */
                double 	H,   /**< Barrier                                   */
                char   	UoD, /**< Up & out or down & out ('U' or 'D')       */
                double 	T,   /**< Option expiration in years                */
                double 	s)   /**< Annualized volatility                     */
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


        if (fabs(S) < TINY)
                return (0.);

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


        if (fabs(S) < TINY)
                return (0.);

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



