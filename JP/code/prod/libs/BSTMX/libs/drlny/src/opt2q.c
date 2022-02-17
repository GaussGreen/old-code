/****************************************************************
 * Module:	DRL
 * Submodule:	OPTIO
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include <math.h>	
#include <string.h>
#include <float.h>	

#include "drlmath.h"		/* DrlCumNorm */

#include "drloptio.h"		/* prototype consistency */


#undef	ISFLAG
#define	ISFLAG(b)		(!strncmp(what,b,strlen(b)))




#define     INVSQRT2PI      0.3989422804   /* 1/sqrt(2*pi)                  */
#define     QCUTOFF         1E-4           /* Normal model for |q|<QCUTOFF  */



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
static	double  NormalH (double  x)
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
static	double Normal_InvH (double prob)
{
    double  t;
    double  x;
    double  r;
   
    if (prob < 0e0 ||
	prob > 1e0)
	GtoErrMsg("Normal_InvH failed.  Input probability (%lf) is invalid.\n",
		  prob);

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
static	double  NDensity (double x)
{
    return (exp (-0.5 * x * x) * INVSQRT2PI);
}



/*****  CCEquation_BS2Q  ****************************************************/
/*
*       Evalute calibration equation for calibration constant.
*/
static	double  CCEquation_BS2Q(
	double    S,         /* (I) Annualized volatility  */
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

DLL_EXPORT(int)
ConvexityC_BS2Q (
	double    S,         /* (I) Annualized volatility  */
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
           CCEq2;




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
            CC0  -= CCEq * CCDELTA / (CCEq2 - CCEq);
        }
        else
        {
            found = TRUE;
        }
        count++ ;
    }
    while (!found && (count < 10));

    if (count < 10)
    {
        *CC = CC0;
        status = SUCCESS;
    }
    else
    {
        *CC = 1.0e10;
    }

    return (status);

}  /* ConvexityC_BS2Q */



/*f--------------------------------------------------------------
 * Options formulae : 2Q pricing
 *
 * <br><br>
 * Note: vol is defined by Y = Y/(1+fsh) * F((1+fsh)/(1+q*fsh)*vol)
 */

DLL_EXPORT(int)
DrlOptBS2Q(
	double Y,                /* (I) Fwd yield              */
	double K,                /* (I) Strike                 */
	double T,                /* (I) Option expiration      */
	double S,                /* (I) Annualized volatility  */
	char   CoP,              /* (I) Call or put            */
	double QLeft,            /* (I) Q left                 */
	double QRight,           /* (I) Q right                */
	double FwdSh,            /* (I) Fwd shift              */
	double *value)		/* (O) option value */
{
static	char	routine[] ="DrlOptBS2Q";
    double
            QM,                            /* Q for fwd shift correction    */
            Q,                             /* Q used in part of density     */
            CC,                            /* Calibration const             */
            YCorr,                         /* Yield correction              */
            C, P,	                       /* Call, put price               */
            d,                             /* d in N(d) in Black & Scholes  */
            St,	                           /* Sigma * sqrt(T) * Q           */
            Y1, K1;                        /* Modified fwd and strike       */


    if (Y == 0.) return (0);

    QM  = (QLeft + QRight) / 2;
    St  = S * sqrt (T) * (1 + FwdSh) / (1 + QM * FwdSh);

    YCorr = Y / (1. + FwdSh);

    if (ConvexityC_BS2Q (St, QLeft, QRight, FwdSh, &CC) == FAILURE)
    {
        GtoErrMsg("Could not find calib const for Call_BS2Q !");
        C = P = -1.0e10;
	return (FAILURE);
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

    FREE_MEM_AND_RETURN:

	*value = ((CoP == 'C') ? C : P);
	return (SUCCESS);

}




#define VOLDELTA  0.0001
#define PREQERR  0.000001

/*f--------------------------------------------------------------
 * Options formulae : implied 2Q volatility
 *
 * <br><br>
 */

DLL_EXPORT(int)
DrlOptBS2QImplVol(
	double Y,                /* (I) Fwd yield              */
	double K,                /* (I) Strike                 */
	double T,                /* (I) Option expiration      */
	double P,                /* (I) Price of option        */
	char   CoP,              /* (I) Call or put            */
	double QLeft,            /* (I) Q left                 */
	double QRight,           /* (I) Q right                */
	double FwdSh,            /* (I) Fwd shift              */
	double VolGuess,         /* (I) Initial vol guess      */
	double *implVol)	/* (O) implied volatility */
{
static	char	routine[] = "DrlOptBS2QImplVol";
    int
           found,
           count;

    double 
           Pr, Pr2, dPr,
           iVol;


    iVol  = VolGuess;
    count = 0;
    found = FALSE;

    do 
    {
	IF_FAILED_DONE( DrlOptBS2Q(
		Y,K,T,iVol,CoP,QLeft,QRight,FwdSh, &Pr));
	IF_FAILED_DONE( DrlOptBS2Q(
		Y,K,T,iVol+VOLDELTA,CoP,QLeft,QRight,FwdSh, &Pr2));

        dPr = (Pr2 - Pr) / VOLDELTA;

        if ((fabs(Pr - P) > PREQERR))
        {
            iVol -= (Pr - P) / dPr;
        }
        else 
        {
            found = TRUE;
        }
        count++ ;
    }
    while (!found && (count < 10));

	if (!found) {
		GtoErrMsg("%s: failed.\n", routine);
		goto done;
	}


	*implVol = iVol;
	return (SUCCESS);
done:
	GtoErrMsg("%s: failed.\n", routine);
	return (FAILURE);
}
                  


/*f--------------------------------------------------------------
 * Options formulae : convert ATM to 2Q volatility.
 *
 * <br><br>
 */

DLL_EXPORT(int)
DrlATMTo2QVol(
	double texp,	/* (I) time to expiration */
	double fwd,	/* (I) exp. value of underlying */
	double vol,	/* (I) volatility */
	KVolType vType,	/* (I) LOGVOL, NORMVOL */
	double strike,	/* (I) strike */

	double QLeft,	/* (I) Q left */
	double QRight,	/* (I) Q right */
	double FwdSh,	/* (I) Fwd shift */

	char *callPut,	/* (I) "C" for call, "P" for put */
	double *retVal)	/* (O) return value */
{
static	char	routine[] = "DrlATMTo2QVol";
	int	status = FAILURE;

	double	T = texp,
		atmPr,
		y;


	/*
	 * A FEW SPECIAL CASES AND SANITARY CHECKS
	 */


	/*
	 * In the case where fwd <= 0, only qL=qR=0 and volType=0 are allowed
	 * in the current 2q mapping, i.e. normal backbone and normal 
	 * distribution.  The benchmark vol input is bp vol.
	 */
	if (IS_ALMOST_ZERO(fwd) ||
	    fwd < 0e0)
	{
	    if(!IS_ALMOST_ZERO(QLeft) ||
	       !IS_ALMOST_ZERO(QRight))
	    {
		GtoErrMsg("%s: in the case of fwd<=0 (at texp =%lf), both "
			"qL(=%lf) and qR(=%lf) must be 0(normal).\n",
			routine,
			texp, QLeft, QRight);
		goto done;
	    }
	
	    if(vType != NORMVOL)
	    {
		GtoErrMsg("%s: in the case of fwd<=0 (at texp =%lf), only "
			"normal vol is allowed.\n",
			routine,
			texp);
		goto done;
	    }

	    *retVal = vol;
	    return (SUCCESS);
	}


	/* 
	 * Trivial cases: 
	 * 1. lognormal smile and % vol input
	 */
	if ((IS_ALMOST_ZERO(QLeft  - 1e0)) &&
            (IS_ALMOST_ZERO(QRight - 1e0)) &&
	    (vType == LOGVOL)) 
	{
	    *retVal = fwd * vol;
	    return (SUCCESS);
	}
		
	/*
	 * 2. normal smile and bp vol input
	 */
        if ((IS_ALMOST_ZERO(QLeft  - 0e0)) &&
            (IS_ALMOST_ZERO(QRight - 0e0)) &&
	    (vType == NORMVOL))
	{
	    *retVal = vol;
	    return (SUCCESS);
		
	}

	 
	/* 
	 * Currently only lognormal or normal bechmark vol inputs
	 * allowed. Compute the ATM price.
	 */
	if (vType == LOGVOL)
		atmPr = fwd * (2. * NormalH(.5 * sqrt(T) * vol) - 1.);
	else if (vType == NORMVOL)
		atmPr = sqrt(T) * vol / sqrt(2. * PI);
	else
	{
		GtoErrMsg("%s: invalid vol type (vType=%lf). "
		  "only lognormal (%d) or normal vol (%d) allowed.\n",
		  routine,
		  vType,
		  LOGVOL,
		  NORMVOL);
		goto done;
	}


        /* Lognormal to X-space vol adjustment. The single and */
        /* two q cases are treated differently for consistency */
        /* with old model.                                     */

	/* "y" is the bp vol 				       */
        if (IS_ALMOST_ZERO(QLeft - QRight) && IS_ALMOST_ZERO(FwdSh))
        {
            if (fabs(QLeft) > QCUTOFF) 
	    {	
		double prob;
		prob = .5 * (QLeft * atmPr / fwd + 1.);

		if (prob < 0e0 ||
		    prob > 1e0)
		{
		    GtoErrMsg("%s: can not solve for 2q vol.  The "
			      "cumulative probability (%lf) is invalid.\n",
			      routine, prob);

	    	    GtoErrMsg("QLeft = QRight = %lf, exp = %lf, fwd=%lf, "
			      "atmPr=%lf.\n", 
			      QLeft, T, fwd, atmPr);

		    goto done;
		}

                y = fwd * Normal_InvH (.5 * (QLeft * atmPr / fwd + 1.)) / (.5 * QLeft);
	    }
            else
                y = atmPr * sqrt (2. * PI);
  
            y /= sqrt(T);
        }
        else
        {
	   	IF_FAILED_DONE ( DrlOptBS2QImplVol(
			fwd, fwd, 
			T,atmPr, callPut[0], QLeft,QRight,FwdSh,vol,
            		&y));
		y *= fwd;
  
        }


	*retVal = y;
	return (SUCCESS);
done:
	GtoErrMsg("%s: failed.\n", routine);
	return (FAILURE);
}

