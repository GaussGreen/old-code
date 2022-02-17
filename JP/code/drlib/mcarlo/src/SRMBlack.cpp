
//----------------------------------------------------------------------------
//
//   Group       : QR cross asset
//
//   Filename    : SRMBlack.hpp
//
//   Description : 2Q Black-Scholes utility functions
//
//   Date        : 2 May 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMBlack.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SRMConstants.hpp"

DRLIB_BEGIN_NAMESPACE

/* Coefficients for cumulative normal */
#define     INVSQRT2PI      0.3989422804   /* 1/sqrt(2*pi)                  */
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

/*****  Normal_InvH  ********************************************************/
/*
*      Normal cumulative inverse distribution Risk Magazine
*/
double SRMBlack::Normal_InvH(double prob)
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
*       Normal density. From bas::NDensity
*/
double SRMBlack::NDensity(double x)
{
    return (exp (-0.5 * x * x) * INVSQRT2PI);
}


/*  Normal cumulative distribution accrding to J.Hull
    This is here to make the numbers match. From bas::NormalH */
double SRMBlack::NormalH(double  x){
    double  y;
    if (x > 0){
        y = x;
    } else {
        y = -x;
    }
    double k = 1. / (1. + h0 * y);

    double norm  = exp (-0.5 * y * y) * INVSQRT2PI;
    norm *= k * (h1 + k * (h2 + k * (h3 + k * (h4 + k * h5))));
    if (x > 0) norm  = 1. - norm;

    return (norm);
}

/*****  Call_BS  ************************************************************/
/*
*      Price of a Call using Black&Scholes. From bas::Call_BS
*/
double SRMBlack::Call_BS(
			    double S,  /* Price of the underlying    */
                double K,  /* Strike                     */
				double T,  /* Option expiration in years */
				double r,  /* Risk free rate             */
				double s)  /* Annualized volatility      */
{

    double  C;     /* Call price */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */


    K  *= exp (-r * T);
    st  = s * sqrt (T);

    if (K <= 0.)
    {
        C = - K + S;
        return C;
    }

    if ((st == 0.) || (S == 0.))
    {
        C = Maths::max(S - K, 0.);
        return C;
    }

    d   = log (S / K) / st;
    C   = S * NormalH (d + 0.5 * st) - K * NormalH (d - 0.5 * st);

    return (C);

}  /* Call_BS */


/*****  Put_BS  *************************************************************/
/*
*      Price of a Put using Black&Scholes. From bas::Put_BS
*/
double SRMBlack::Put_BS(
			   double S,  /* Price of the underlying    */
               double K,  /* Strike                     */
			   double T,  /* Option expiration in years */
			   double r,  /* Risk free rate             */
			   double s)  /* Annualized volatility      */
{
    double  P;     /* Put price */
    double  d;     /* d in N(d) in Black & Scholes formula */
    double  st;    /* Sigma * sqrt(T) */


    K  *= exp (-r * T);
    st  = s * sqrt (T);
    
    if (K <= 0.)
        return (0.);

    if ((st == 0.) || (S == 0.))
    {
        P = Maths::max(K - S, 0.);
        return P;
    }

    d   = log (K / S) / st;
    P   = K * NormalH (d + 0.5 * st) - S * NormalH (d - 0.5 * st);

    return (P);

}  /* Put_BS */

/*****  CCEquation_BS2Q  ****************************************************
 *       Evalute calibration equation for calibration constant.
 *       From bas::CCEquation_BS2Q
 */
double SRMBlack::CCEquation_BS2Q(
    double S,         /* (I) Annualized volatility  */
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

    if (fabs (QRight) >SRMConstants::QCUTOFF)
    {
        CCEq -= (exp (CCIS * SQR + 0.5 * SQR * SQR) * NormalH (+ CCIS + SQR)
                 - NormalH (+ CCIS) ) / QRight;
    }
    else
    {
        CCEq -= (CCIS * NormalH(+ CCIS) + NDensity(CCIS)) * S;  
    }

    if (fabs (QLeft) >SRMConstants::QCUTOFF)
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
*      Calibrate constant in 2q formulas. From bas::ConvexityC_BS2Q
*/
#define CCDELTA  0.0001
#define CCEQERR  0.000001
double SRMBlack::ConvexityC_BS2Q(
					   double    S,     /* (I) Annualized volatility  */
                       double    QLeft, /* (I) Q left                 */
                       double    QRight,/* (I) Q right                */
                       double    FwdSh) /* (I) Fwd shift              */
{
    int count = 0;

    double 
        QM,     
        CC0,
        CCEq,
        CCEq2;

    QM = (QLeft + QRight) * 0.5;
    if (fabs (QM) >SRMConstants::QCUTOFF)
    {
        CC0 = log (FwdSh * QM + 1.) / QM - 0.5 * QM * S * S;
    }
    else
    {
        CC0 = FwdSh;
    }

    bool found = false;
    do 
    {
        CCEq  = CCEquation_BS2Q(S, QLeft, QRight, FwdSh, CC0);
        CCEq2 = CCEquation_BS2Q(S, QLeft, QRight, FwdSh, CC0 + CCDELTA);
        
        if (fabs(CCEq) > CCEQERR)
        {
            CC0  -= CCEq * CCDELTA / (CCEq2 - CCEq);
        }
        else
        {
            found = true;
        }
        count++ ;
    }
    while (!found && (count < 10));

    if (count < 10) {
        return CC0;
    } else {
        throw ModelException("ConvexityC_BS2Q", "Failed to calibrate constant"
                             " in 2q formulas");
    }
}  /* ConvexityC_BS2Q */

/*****  Option_BS2Q  ********************************************************/
/*
*      Price of a call or put using 2Q version of Black&Scholes.
*      Note: vol is defined by Y = Y/(1+fsh) * F((1+fsh)/(1+q*fsh)*vol)
*      From bas::Option_BS2Q
*/
double SRMBlack::Option_BS2Q(
    double Y,                /* (I) Fwd yield              */
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
            C, P,                          /* Call, put price               */
            d,                             /* d in N(d) in Black & Scholes  */
            St,                            /* Sigma * sqrt(T) * Q           */
            Y1, K1;                        /* Modified fwd and strike       */

    /* lognormal case */
    if (Maths::areEqualWithinTol(QLeft, 1.0, SRMConstants::SRM_TINY) &&
        Maths::areEqualWithinTol(QRight, 1.0, SRMConstants::SRM_TINY) &&
        Maths::areEqualWithinTol(FwdSh, 0.0, SRMConstants::SRM_TINY))
    {
        return (CoP == 'C'? Call_BS(Y, K, T, 0.0, S): Put_BS(Y, K, T, 0.0, S));
    }

    if (Y == 0.) return (0.);

    QM  = (QLeft + QRight) / 2;
    St  = S * sqrt (T) * (1 + FwdSh) / (1 + QM * FwdSh);

    YCorr = Y / (1. + FwdSh);

    CC = ConvexityC_BS2Q (St, QLeft, QRight, FwdSh);

    Q = (K > YCorr) ? QRight : QLeft;
    
    if (fabs(Q) >SRMConstants::QCUTOFF)
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

    return ((CoP == 'C') ? C : P);
}  /* Opt_BS2Q */


/*****  ImpVol_BS2Q  ********************************************************/
/*
*      Price of a call or put using 2Q version of Black&Scholes.
*      From bas::ImpVol_BS2Q
*/
#define VOLDELTA  0.0001
#define PREQERR  0.000001

double SRMBlack::ImpVol_BS2Q(
    double Y,                /* (I) Fwd yield              */
    double K,                /* (I) Strike                 */
    double T,                /* (I) Option expiration      */
    double P,                /* (I) Price of option        */
    char   CoP,              /* (I) Call or put            */
    double QLeft,            /* (I) Q left                 */
    double QRight,           /* (I) Q right                */
    double FwdSh,            /* (I) Fwd shift              */
    double VolGuess){         /* (I) Initial vol guess      */
    double scale = Y * sqrt(T / 6.28);
    double iVol  = VolGuess;
    int count = 0;
    bool found = false;

    do {
        double Pr  = Option_BS2Q (Y,K,T,iVol,CoP,QLeft,QRight,FwdSh);
        double Pr2 = Option_BS2Q (Y,K,T,iVol+VOLDELTA,CoP,QLeft,QRight,FwdSh);

        if ((fabs(Pr - P) > scale * PREQERR)) {
            double dPr = (Pr2 - Pr) / VOLDELTA;
            if (fabs(dPr) < scale * PREQERR) {
                return (-1.0e5);
            }
            iVol -= (Pr - P) / dPr;
        } else {
            found = true;
        }
        count++ ;
    }
    while (!found && (count < 10));
                  
    /* failure if exceed 10 iterations */
    if (count >= 10) iVol = -1.0e5;

    return (iVol);
}
                  
DRLIB_END_NAMESPACE
