/*********************************************************************************
 * OPTBSQ.C
 * basic BS and BSQ option routines
 *
 ********************************************************************************/
#include <crxflow/include/optbsq_x.h>
#include <alib/rtbrent.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "crxmacros.h"

#define  INVSQRT2PI             0.398942280401433   /* 1/sqrt(2*pi) */

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
 *      BS price formula
 *      no discounting!
 */
int BSPricer(
    double     *result,                   /* (O) Price, Vega, or Delta          */
    KOptType   optType,                   /* (I) Option type                    */
	double     Y,                         /* (I) Fwd yield                      */
	double     K,                         /* (I) Strike                         */
	double     T,                         /* (I) Option expiration in years     */
	double     s,                         /* (I) Annualized volatility          */
    KOptResult optResult)                 /* (I) Return result type             */
{
    static char      routine[]="BSPricer";
    int              status    = FAILURE;
	double  C;                            /* Call price                    */
	double  d1, d2;                       /* d1, d2 in BS formula          */
	double  dtmp;
    double  price, delta;
    double  exProb=0.0;                     /* Ex prob = N(d2) for C, 1-N(d2) P */
    double  nd1, nd2;                     /* N(d1) and N(d2)               */

    if (Y < OPT_MIN_FWD)
    {
        DR_Error("%s: Fwd Yield %g too small.\n", routine, Y);
        goto RETURN;
    }
	
    /* for zero volatility return intrinsic */
    if (s * sqrt(T) < OPT_MIN_VOL) 
    {
        C = MAX(Y-K, 0.);

        switch(optType) 
        {
        case OPT_CALL: 
            price = C;
            delta = (Y > K)?1.0:0.0;
            break;
        case OPT_PUT:            
            price = C -(Y-K);
            delta = (K > Y)?-1.0:0.0;
            break;
        default:
            DR_Error("%s: Unknown option type.\n", routine);
            goto RETURN;
        }
    } else {
        dtmp = s*sqrt(T);
        d1 = (log(Y/K) + (s*s/2)*T)/dtmp;
        d2 = d1 - dtmp;
        
        nd1 = NormCum(d1);
		nd2 = NormCum(d2);
        C = Y * nd1 - K * nd2;
        
        switch(optType) 
        {
		case OPT_CALL: 
            price = C;
            delta = nd1;
			exProb = nd2;
            break;
		case OPT_PUT:            
            price = C -(Y-K);
            delta = nd1 - 1.0;
			exProb = 1 - nd2;
            break;
		default:
            DR_Error("%s: Unknown option type.\n", routine);
            goto RETURN;
        }
    }

    switch(optResult){
    case OPT_PRICE:
        *result = price;
        break;
    case OPT_DELTA:
        *result = delta;
        break;
    case OPT_GAMMA:
        if (s * sqrt(T) < OPT_MIN_VOL)
        {
            *result = 0.0;
            if(fabs(Y-K) < TINY) *result = 1.0/TINY;
        } else {
            dtmp = s*sqrt(T);
            d1 = (log(Y/K) + (s*s/2)*T)/dtmp;

            *result = exp(-d1 * d1/2.0) * INVSQRT2PI /dtmp/Y;
        }
        break;
	case OPT_EXERCISE_PROBABILITY:
		*result = exProb;
		break;
    default:
        DR_Error("%s: Unknown option result.\n", routine);
        goto RETURN;
    }
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        DR_Error("%s: Failed.\n", routine);
    }
    
    return status;
}

#define VOLDELTA 1E-4
#define PREQERR  1E-8

/*f----------------------------------------------------------------------------
 *      Implied vol of a call or put using Black-Scholes.
 */
int     BSImpVol   (
    double     *impVol,                   /* (O) Implied BS vol                 */
    double     yield,                     /* (I) Fwd yield                      */
    double     strike,                    /* (I) Strike                         */
    double     expiry,                    /* (I) Option expiration              */
    double     price,                     /* (I) Price of option                */
    KOptType   optType,                   /* (I) Option type                    */
    double     volGuess)                  /* (I) Initial vol guess              */
{
    int    found, count;
    double Pr, Pr2, dPr, scale, iVol;

    static char      routine[]="BSImpVol";
    int              status = FAILURE;
        
    scale = yield * sqrt(expiry / 6.28);
    iVol  = volGuess;
    count = 0;
    found = FALSE;

    do {
        if (BSPricer (&Pr,
                      optType,
                      yield,
                      strike,
                      expiry,
                      iVol,
                      OPT_PRICE) == FAILURE) goto RETURN;

        if (BSPricer (&Pr2,
                      optType,
                      yield,
                      strike,
                      expiry,
                      iVol+VOLDELTA,
                      OPT_PRICE) == FAILURE) goto RETURN;

        if (fabs(Pr - price) > scale * PREQERR) 
        {
            dPr = (Pr2 - Pr) / VOLDELTA;
            if (fabs(dPr) < OPT_MQ_RESN) 
            {
                DR_Error("%s: numerical derivative is 0.\n", routine);
                goto RETURN;
            }

            iVol -= (Pr - price) / dPr; 
        } 
        else 
        {
            found = TRUE;
        }
        count++;
    } while ((found != TRUE) && (count < 10));
                  
    /* failure if exceed 10 iterations */
    if (found == TRUE) 
    {
        *impVol = iVol;
        status = SUCCESS;
    } 
    else 
    {
        DR_Error ("%s: exceeded maximum nb iterations.\n", routine);        
        *impVol = -999;
    }

 RETURN:

    return status;

} /* BSImpVol */


double  BSQIntCR (double x1, double x2,double a,double b,double q,double x0,double C,double sig)
{
  double sx2 = x2/sig;
  double sx1 = x1/sig;
  double dtmp;
  

  if(IS_ZERO(q))
  {

	dtmp = C*((a - b*x0)*(NormCum(sx2) - NormCum(sx1))+
			  b * sig * (exp(-0.5*sx1*sx1) - exp(-0.5*sx2*sx2))*INVSQRT2PI);
	
	//		printf("here C %f a %f b %f x0 %f sx1 %f sx2 %f,dtmp %f,sig %f\n",C,a,b,x0,sx1,sx2,dtmp,sig);
	return dtmp;

  } else {

	return C*((a-b/q) * (NormCum(sx2) - NormCum(sx1)) 
	  + (b/q) * exp(0.5*q*q*sig*sig - q*x0)*(NormCum(sx2 - q*sig) - NormCum(sx1 - q*sig)));

  }
}


/******************************************************************************
 * If fabs(q) is too mall, replace it with 1E-7
 *****************************************************************************/
#define CRX_Q_CUTOFF 1E-7

/******************************************************************************
 * Crx2QCalibrate
 * Given the Q-Smile details, calculates the 2-Q 'A' proportionality constant
 * and the 2Q total volatility sigmaQ
 *****************************************************************************/
int Crx2QCalibrate(
    double  forward,        /**<The forward that you want to match           */
    double  atmOptionPrice, /**<The ATM option price you want to match       */
    double  qL,             /**<Left q parameter                             */
    double  qR,             /**<Right q parameter                            */
    double  muQ,            /**<Q forward shift parameter                    */
    double* A,              /**<(O) Calibrated Q A proprtionality parameter  */
    double* sigmaQ          /**<(O) Calibrated Q total volatility parameter  */
                   ) {

    int status = FAILURE;
    const static char* function = "Crx2QCalibrate";
    double a, s, slast; // temp values for A and sigmaQ
    double f, p; // temp values for forward and option price
    int maxIter = 30; // maximum number of iterations
    int it; // current iteration
    double errBound = 1E-4; // max percentage error
    double errSize; // current percentage error

    /*=========================================================================
     * INPUT ARGUMENT CHECKS
     *=======================================================================*/
    if (forward<DBL_EPSILON) {
        DR_Error("%s failed - forward too small or negative %lf", 
            function, forward);
        goto RETURN;
    } else if (atmOptionPrice<DBL_EPSILON) {
        DR_Error("%s failed - ATM option price too small or negative %lf", 
            function, atmOptionPrice);
        goto RETURN;
    }

    /*=========================================================================
     * TO AVOID SCALING AND ROUNDING ERRORS, SET FORWARD TO 1 AND RESCALE
     * OPTION PRICE - UNSCALE AGAIN BEFORE RETURNING
     *=======================================================================*/
    atmOptionPrice /= forward;

    /*=========================================================================
     * INITIALIZATION - FIRST VOL GUESS IS BS IMPLED TOTAL VOL
     *=======================================================================*/
    if (FAILURE==BSImpVol(&s,1.0,1.0,1.0,atmOptionPrice,OPT_CALL,0.2)) {
        DR_Error("%s failed - couldn't compute initial total vol guess "
            "value from forward=%lf, atmOptionPrice=%lf ", 
            function, 1.0, atmOptionPrice);
        goto RETURN;
    }
    if (FAILURE==Crx2QForward(1.0,muQ,s,qL,qR,&f)) {
       DR_Error("%s failed - couldn't compute initial forward "
            "value from a=%lf, s=%lf ", 
            function, 1.0, s);
        goto RETURN;
    }
    a = 1.0/f;
    if (FAILURE==Crx2QOptionPrice(GtoOPTION_CALL,1.0,a,muQ,s,qL,qR,&p)) {
        DR_Error("%s failed - couldn't compute initial option price "
            "from a=%lf, s=%lf ", 
            function, a, s);
        goto RETURN;
    }
    
    /*=========================================================================
     * ONE DIMENSIONAL NUMERICAL NEWTON ITERATION
     * USE ABSOLUTE PERCENTAGE DIFFERENCE AS CONVERGENCE TEST
     *=======================================================================*/
    it = 0;
    errSize = fabs(p/atmOptionPrice-1);
    slast = 0.0;
    while (it<maxIter && errSize>errBound) {
        double dpds;
        double dsdp;
        double stweak;

        /* Compute derivative */
        stweak = fabs(s-slast)/100.0;
        if (stweak<1E-9) stweak = 1E-5;
        // recompute tweaked a
        if (FAILURE==Crx2QForward(a,muQ,s+stweak,qL,qR,&f)) {
            DR_Error("%s failed - failed to compute dfds at iteration %d: "
                "a=%lf, s=%lf, f=%lf, p=%lf, stweak=%lf", 
                function, a, s, f, p, stweak);
            goto RETURN;
        }
        a *= 1.0/f;
        f = 1.0;
        // compute tweaked option price
        if (FAILURE==Crx2QOptionPrice(GtoOPTION_CALL,1.0,a,muQ,s+stweak,qL,qR,&dpds)) {
            DR_Error("%s failed - failed to compute dpda at iteration %d: "
                "a=%lf, s=%lf, f=%lf, p=%lf, stweak=%lf", 
                function, a, s, f, p, stweak);
            goto RETURN;
        }
        dpds = (dpds-p)/stweak;
        if (fabs(dpds)<DBL_EPSILON) {
            DR_Error("%s failed - zero derivative at iteration %d: "
                "a=%lf, s=%lf, f=%lf, p=%lf, dpds=%lf", 
                function, a, s, f, p, dpds);
            goto RETURN;
        }
        dsdp = 1/dpds;


        /* Next iteration value */
        slast = s;
        s += (atmOptionPrice-p)*dsdp;
        if (s<DBL_EPSILON) s=slast/2; // keep in reasonable range if step below zero

        // recalculate forward scale parameter so forward still 1
        if (FAILURE==Crx2QForward(a,muQ,s,qL,qR,&f)) {
            DR_Error("%s failed - couldn't compute forward at iteration %d"
                "value from a=%lf, s=%lf ", 
                function, it, a, s);
            goto RETURN;
        }
        a *= 1.0/f;
        f = 1.0;
        if (FAILURE==Crx2QOptionPrice(GtoOPTION_CALL,1.0,a,muQ,s,qL,qR,&p)) {
            DR_Error("%s failed - couldn't compute option price at iteration %d"
                "from a=%lf, s=%lf ", 
                function, it, a, s);
            goto RETURN;
        }
        
        errSize = fabs(p/atmOptionPrice-1);
        it++;
    }
    if (it>=maxIter || errSize>errBound) {
        DR_Error("%s failed - Newton's method failed to converge "
            "with iterations=%d, current pct error=%lf ", 
            function, it, errSize);
        goto RETURN;
    }

    *A = a*forward; // rescale to actual forward
    *sigmaQ = s;
    status = SUCCESS;
RETURN:
    return status;
}

/******************************************************************************
 * Crx2QForward
 * Calculates the forward price under a 2-Q model
 * The underlying price, P, is assumed to be a function of a normal R Var, X:
 * P = A(1 + (exp(X qL)-1))/qL) if P<=F
 *   = A(1 + (exp(X qL)-1))/qR) if P>=F
 * where X ~ Normal(F, s^2)
 * Returns FAILURE if anything goes wrong, otherwise SUCCESS
 *****************************************************************************/
int Crx2QForward(
    double A,       /**<(I) Proportionality factor */
    double muQ,       /**<(I) Expectation of driver variable, a.k.a. forward shift */
    double sigQ,       /**<(I) Total volatility (s^2 = sigma^2 T) of driver variable */
    double qL,      /**<(I) Left side Q                                      */
    double qR,      /**<(I) Right side Q                                     */
    double* forward /**<(O) Calculated forward                               */
    ) 

{
    int status = FAILURE;
    const static char* function = "Crx2QForward";
    double alpha;

    /*=========================================================================
     * ARGUMENT CHECKS
     *=======================================================================*/
    if (sigQ<0)
    {
        DR_Error("%s failed - negative driver volatility %lf", function, sigQ);
        goto RETURN;
    }
    
    /*=========================================================================
     * ZERO VOLATILITY CASE
     *=======================================================================*/
    if (sigQ<DBL_EPSILON) 
    {
        // effectively zero volatility
        if (muQ>=0) 
        {
            // right side q at expectation
            if (fabs(qR)<CRX_Q_CUTOFF)
            {   
                *forward = A*(1.0 + muQ);
            }
            else 
            {
                *forward = A*(1 + (exp(qR*muQ)-1)/qR);
            }
        }
        else 
        {
            // left side q at expectation
            if (fabs(qL)<CRX_Q_CUTOFF)
            {   
                *forward = A*(1.0 + muQ);
            }
            else 
            {
                *forward = A*(1 + (exp(qR*muQ)-1)/qR);
            }
        }
        status = SUCCESS;
        goto RETURN;
    }

    alpha = muQ/sigQ;  // std normal space equiv of X=0    
    *forward = 1.0;
    /*=========================================================================
     * LEFT Q CONTRIBUTION
     *=======================================================================*/
    if (fabs(qL)<CRX_Q_CUTOFF) 
    {
        // normal distribution
        *forward += muQ*GtoNormalCum(-alpha) + sigQ*exp(-alpha*alpha/2)*INVSQRT2PI;
    }
    else 
    {
        // non-normal q-distribution
        double sL = sigQ*qL;
        *forward += 
            (-GtoNormalCum(-alpha) 
            + exp(sL*(alpha+0.5*sL)) * GtoNormalCum(-alpha - sL)
            )
            /qL;
    }
    /*=========================================================================
     * RIGHT Q CONTRIBUTION
     *=======================================================================*/
    if (fabs(qR)<CRX_Q_CUTOFF) 
    {
        // normal distribution
        *forward += muQ*GtoNormalCum(alpha) - sigQ*exp(-alpha*alpha/2)*INVSQRT2PI;
    }
    else 
    {
        // non-normal q-distribution
        double sR = sigQ*qR;
        *forward += 
            (-GtoNormalCum(alpha) 
            + exp(sR*(alpha+0.5*sR)) * GtoNormalCum(alpha + sR)
            )
            /qR;
    }


    *forward *= A;

    status = SUCCESS;
RETURN:
    return status;

} /* Crx2QForward */


/******************************************************************************
 * Crx2QOptionPrice
 * Calculates the price of call and put options under a 2-Q model
 * The underlying price, P, is assumed to be a function of a normal RV, X:
 * P = A(1 + (exp(X qL)-1))/qL) if P<=A
 *   = A(1 + (exp(X qL)-1))/qR) if P>=A
 * where X ~ Normal(F, s^2 T)
 * Returns FAILURE if anything goes wrong, otherwise SUCCESS
 *****************************************************************************/
int Crx2QOptionPrice(
    int optType, /**<GtoOPTION_CALL or GtoOPTION_PUT */
    double   K, /**<(I) Strike                         */
    double   A, /**<(I) Proportionality Factor                       */
    double   F, /**<(I) Expectation of driver variable, X (Forward Shift)*/
    double   s, /**<(I) Total volatility of driver, X (sigma sqrt(T)        */
    double   qL, /**<(I) Left side Q */
    double   qR, /**<(I) Right side Q*/
    double*  price /**<(O) Option price*/
    ) 
{
    int status = FAILURE;
    const static char* function = "Crx2QOptionPrice";

    /*=========================================================================
     * ARGUMENT CHECKS
     *=======================================================================*/
    if (s<0)
    {
        DR_Error("%s failed - negative total volatility %lf", function, s);
        goto RETURN;
    }

    /*=========================================================================
     * ZERO VOL - JUST INTRINSIC VALUE
     *=======================================================================*/
    if (s<DBL_EPSILON)
    {
        double fwd;
        if (FAILURE==Crx2QForward(A,F,s,qL,qR,&fwd))
        {
            DR_Error("%s failed - bad 0 vol forward for F=%lf, A=%lf",
                function,F,A);
            goto RETURN;
        }

        switch (optType)
        {
        case GtoOPTION_CALL:
            *price = (K>=fwd ? 0.0 : fwd-K);
            break;
        case GtoOPTION_PUT:
        default:
            *price = (K<=fwd ? 0.0 : K-fwd);
            break;
        }
        status = SUCCESS;
        goto RETURN;
    }

    /*=========================================================================
     * USE CALL-PUT PARITY TO GUARANTEE ONLY 1-SIDE INTEGRATION
     *=======================================================================*/
    if (K>=A)
    {
        /*=====================================================================
         * RIGHT-SIDE OPTION
         *===================================================================*/
        // evaluate call price
        if (fabs(qR<DBL_EPSILON)) 
        {
            /* NORMAL CASE */
            double d1 = (1 - K/A + F)/s;
            *price = A*(GtoNormalCum(d1)*(F + 1 - K/A) + s*GtoNormalDen(d1));
        }
        else 
        {   
            /* Q!=0 CASE */
            double qSigma = qR*s;
            double d1 = -log( (K/A-1)*qR + 1 )/qSigma + qSigma + F/s;
            double d2 = d1 - qSigma;
            double nd1 = GtoNormalCum(d1);
            double nd2 = GtoNormalCum(d2);
            *price = A*
                (GtoNormalCum(d1)*exp(qR*(F + qR*s*s/2))/qR 
                + GtoNormalCum(d2)*(1 - 1/qR - K/A));
        }
        // if put, use call-put parity
        if (optType==GtoOPTION_PUT)
        {
            double fwd;
            if (FAILURE==Crx2QForward(A,F,s,qL,qR,&fwd))
            {
                DR_Error("%s failed - bad fwd for F=%lf, A=%lf, s=%lf",
                    function,F,A,s);
                goto RETURN;
            }
            *price += (K-fwd);
        }
    }
    else if (K<A)
    {
        /*=====================================================================
         * LEFT-SIDE OPTION
         *===================================================================*/
        // evaluate put price
        if (fabs(qR<DBL_EPSILON)) 
        {
            double d1 = (1 - K/A + F)/s;
            *price = -A*(GtoNormalCum(-d1)*(F + 1 - K/A) - s*GtoNormalDen(-d1));
        }
        else 
        {   
            double qSigma = qL*s;
            double d1 = -log((K/A-1)*qL + 1)/qSigma + qSigma + F/s;
            double d2 = d1 - qSigma;
            *price = A*
                (-GtoNormalCum(-d1)*exp(qL*(F + qL*s*s/2))/qL 
                - GtoNormalCum(-d2)*(1 - 1/qL - K/A));
        }
        // if put, use call-put parity
        if (optType==GtoOPTION_CALL)
        {
            double fwd;
            if (FAILURE==Crx2QForward(A,F,s,qL,qR,&fwd))
            {
                DR_Error("%s failed - bad fwd for F=%lf, A=%lf, s=%lf",
                    function,F,A,s);
                goto RETURN;
            }
            *price += (fwd-K);
        }
    }

    status = SUCCESS;
RETURN:
    return status;

} /* Crx2QOptionPrice */

/* VALUE FUNCTION FOR ALIB GtoRootFindBrent */
int Crx2QVolFinder(double x,  /* evaluate at this point */
                   void* data,/* other arguments */
                   double*y   /* result */
                   )
{
    const static char* function = "Crx2QVolFinder";
    int status = FAILURE;
    
    double* vals = (double*)data;
    double target = vals[0];
    int optType = (vals[1]==0.0 ? GtoOPTION_CALL : GtoOPTION_PUT);
    double px;
    
    if (FAILURE==Crx2QOptionPrice(
        optType,
        vals[4] /*K*/,
        //vals[5] /*T*/, // THIS IS NOT USED - TOTAL VOLATILITY INSTEAD!
        vals[3] /*A*/,
        vals[2] /*F*/,
        x,
        vals[6] /*qL*/,
        vals[7] /*qR*/,
        &px))
    {
        // error will have been logged in Crx2QOptionPrice
        goto RETURN;
    }

    *y = target - px;
    status = SUCCESS;
RETURN:
    return status;
}

/******************************************************************************
 * Crx2QImpliedVol
 * Calculates the implied vol of call and put options under a 2-Q model
 * The underlying price, P, is assumed to be a function of a normal RV, X:
 * P = A(1 + (exp(X qL)-1))/qL) if P<=A
 *   = A(1 + (exp(X qL)-1))/qR) if P>=A
 * where X ~ Normal(F, s^2 T)
 * This routine calculates s so that the price matches that given
 * Returns FAILURE if anything goes wrong, otherwise SUCCESS
 *****************************************************************************/
int Crx2QImpliedVol(
    int optType,    /**<GtoOPTION_CALL or GtoOPTION_PUT                      */
    double   K,    /**<(I) Strike                                            */
    double   A,     /**<(I) Proportional forward                                    */
    double   F,     /**<(I) Forward Shift                */
    double   price, /**<(I) Option price to match                            */
    double   qL,    /**<(I) Left side Q (0 is normal)                        */
    double   qR,    /**<(I) Right side Q (0 is normal)                       */
    double*  vol    /**<(O) Resulting 2-Q total volatility                         */
    ) 
{
    
    static const char* function = "Crx2QImpliedVol";
    int status = FAILURE;
    double data[9];
    data[0] = price;
    data[1] = (optType==GtoOPTION_CALL ? 0.0 : 1.0);
    data[2] = F;
    data[3] = A;
    data[4] = K;
    //data[5] = T;
    data[6] = qL;
    data[7] = qR;

    if (price<=DBL_EPSILON)
    {
        DR_Error("%s failed - non positive price %lf", function, price);
        goto RETURN;
    }
    else if (FAILURE==
        GtoRootFindBrent(Crx2QVolFinder, 
            (void*)data, 
            DBL_EPSILON /*Min vol*/,
            100 /*Max vol*/, 
            30 /* max interations */,
            0.1 /* initial guess */,
            0.001 /* Initial step */,
            0.0 /* no knowledge of initial derivative (vega) */,
            0.0005 /* 1/20% accuracy desired */,
            1E-7 /* 1/1000 of a b.p. of the premium is f-tolerance */,
            vol)) 
    {
        DR_Error("%s failed - GtoRootFindBrent failed for " 
            "Px=%lf, F=%lf, A=%lf, K=%lf, qL=%lf, qR=%lf",
            function, price, F, A, K, qL, qR);
    }

    status = SUCCESS;
 RETURN:
    return status;

} /* Crx2QImpliedVol */



