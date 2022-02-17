/******************************************************************************
 * Module:      CRXQ
 * File:        crxq.c       
 * Author:      Credit QRD (JC Porras orig + modifications for 
 *              bivariates by Charles Morcom)
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "crxq.h"

#include "q3.h" // for BSQInt and a few other functions
#include "crxerror.h" // for DR_Error

/* simple function to calculate (exp(qi(x-xi))-1)/qi */
double qxq(double q, double x, double x0)
{
  if(fabs(q) < TINY) {
    return x-x0;
  }
  else {
    return (exp(q*(x-x0))-1.0)/q;
  }
}

/* BSQ macro: Q Black-Scholes on an interval BSQInt from q3.h */
#define  BSQ_INT(m,s,q,x,y) \
    (BSQInt((q)*(s),(q)*(m),(((x)-(m))/(s)),(((y)-(m))/(s))))

/* MQ pricer for ___/|___ or ___|\___ payoff -> see end of file */
static int MQSharkFin(double,double,double,double,double *,double *);

/*f----------------------------------------------------------------------------
 * Single Q Black Pricer
 *
 * Price of a Call/Put with Vega using Q version of Black&Scholes.
 */

int CRXQBSQPricer(
    double  Y,                  /* Fwd yield                     */
    double  K,                  /* Strike                        */
    double  T,                  /* Option expiration in years    */
    double  s,                  /* Annualized volatility         */
    double  Q,                  /* Q weight                      */
    long    I,                  /* Instrument                    */
    double *P)                  /* Price & Vega                  */
{

    static char routine[] = "CRXQBSQPricer";
    double  C;                  /* Call price                    */
    double  V = 0.;             /* Vega                          */
    double  d;                  /* d in N(d) in Black & Scholes  */
    double  st;                 /* Sigma * sqrt(T) * Q           */
    double  Y1, K1;             /* Modified fwd and strike       */
    double  r;                  /* Adjusted strike               */
    long    vOn;                /* Vega calc on or off           */


    if (Y < CRXQ_MIN_FWD) 
    {
        DR_Error("%s failed: forward %lf less than minimum permitted %lf.\n", 
            routine, Y, CRXQ_MIN_FWD);
        return FAILURE;
    }

    /* compute vega? */
    vOn = (CRXQ_PAY_TYPE(I) == CRXQ_VEGA);

    /* for zero volatility return intrinsic */
    if (s * sqrt(T) < CRXQ_MIN_VOL) 
    {
        C = MAX(Y-K, 0.);
        V = 0.;

        switch(CRXQ_COP_TYPE(I)) 
        {
        case CRXQ_CALL: 
            P[0] = C; 
            if (vOn) P[1] = V; 
            break;
        case CRXQ_PUT:            
            P[0] = C -(Y-K); 
            if (vOn) P[1] = V; 
            break;
        default:
            DR_Error("%s failed: illegal option type %d.\n", 
                routine, CRXQ_COP_TYPE(I));
            return FAILURE;
        }

        return SUCCESS;
    }

    if (fabs(Q) > CRXQ_Q_SHIFT) 
    {
        st  = s * sqrt(T) * Q;
        r   = (K/Y-1.) * Q + 1.;
        /* if K is not in the range of the distribution */
        if (r < TINY) 
        {
            C = MAX(Y-K, 0.);
            V = 0.;
        } 
        else 
        {
            d   = - log(r) / st;
            st /= 2.;

            Y1  = Y / Q;
            K1  = K - Y + Y1;

            C   = Y1 * NormCum(d+st) - K1 * NormCum(d-st);        
            if (vOn) V = Y * exp(-(d*d+st*st)/2.) * sqrt(r*T) * INVSQRT2PI;
        }
    } 
    else 
    {
        st  = s * sqrt(T);
        d   = - (K/Y-1.)/st;

        C   = Y * st * exp(-d*d/2) * INVSQRT2PI + (Y-K) * NormCum(d);
        if (vOn) V = Y * exp(-(d*d)/2.) * sqrt(T) * INVSQRT2PI;
    }

    switch(CRXQ_COP_TYPE(I)) 
    {
    case CRXQ_CALL: 
        P[0] = C; 
        if (vOn) P[1] = V; 
        break;
    case CRXQ_PUT:            
        P[0] = C- (Y-K); 
        if (vOn) P[1] = V; 
        break;
    default:
        DR_Error("%s failed: X illegal option type %d.\n", 
                routine, CRXQ_COP_TYPE(I));
        return FAILURE;
    }

    return SUCCESS;

}  /* CRXQBSQPricer */



/*f----------------------------------------------------------------------------
 * MultiQ Black pricer     
 *
 * Uses Q parameters and the normal mean/stdev to compute call/put price
 */
int CRXQPricer(
    const CRXQDATA     *mq,           /* (I) MQ data  */
    long              optType,       /* (I) option type */
    double            strike,        /* (I) option strike */
    double            *price)        /* (O) option price & vega */
{
    static char      routine[]="CRXQPricer";
    int              status = FAILURE;

    double fwdRate = mq->fwdRate;
    double muMQ    = mq->muMQ;
    double sigMQ   = mq->sigMQ;

    long   nq;          /* number of intervals                    */
    double strkMult;    /* strike/fwdRate                         */
    double strkMultMQ;  /* strike/fwdRate after affine adjustment */
    const double *q;          /* pointer to q parameters                */
    const double *k;          /* pointer to strike ratios               */
    const double *d;          /* pointer to strike slopes               */
    double *x;          /* pointer to gaussian space boundaries   */

    /* affine adjustment constants */
    double K = mq->K;
    double C = mq->C;

    double temp, xtemp, d2q;
    long   i, idx;
    int    sign;

    /* check option type */
    if (optType != CRXQ_CALL && optType != CRXQ_PUT) goto RETURN;

    /* check forward */
    if (fwdRate < CRXQ_MIN_FWD) goto RETURN;

    /* for zero volatility, price intrinsic */
    if (sigMQ < CRXQ_MIN_VOL) 
    {
        *price = MAX(CRXQ_COP_TO_COP(optType) * (fwdRate - strike), 0);
        status = SUCCESS;
        goto RETURN;
    }

    /* branch into low/high strikes (left/right); subtle      */
    /* choice needed in calibration: for ATM option, use left */
    /* intervals for put and high intervals for call          */
    strkMult = strike/fwdRate;

    /* apply affine adjustment */
    strkMultMQ = (K-1)/(K*C) + strkMult/ (K*C);

    if (strkMultMQ < 1.) 
    {
        nq = mq->nbQL;
        q = mq->qL;
        k = mq->kL;
        d = mq->dL;
        x = (double*)mq->xL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        k = mq->kR;
        d = mq->dR;
        x = (double*)mq->xR;
        sign = 1;
    }

    /* Gaussian integration cutoff for the last value at "infinity" */
    x[nq] = muMQ + sign * (CRXQ_CUM_NUM_STDEV) * sigMQ;
      
    /* find strike location */
    idx = 0;
    while (idx < nq - 1 && sign * (strkMultMQ - k[idx+1]) >= 0.) idx++;

    /* find normal point corresponding to strike */
    d2q = d[idx]/q[idx];
    temp = 1 + (strkMultMQ - k[idx])/d2q;
    xtemp = (temp > 0 ? x[idx] + log(temp)/q[idx] : x[nq]);

    /* key step(!): price put for low strikes, call for high strikes */

    /* last interval */
    temp = exp(-q[idx] * x[idx]);
    *price = 
        strkMultMQ         * BSQ_INT(muMQ, sigMQ, 0,      x[nq],    xtemp) -
        (k[idx] - d2q)     * BSQ_INT(muMQ, sigMQ, 0,      x[idx+1], xtemp) -
        d2q * temp         * BSQ_INT(muMQ, sigMQ, q[idx], x[idx+1], xtemp);

    /* other intervals */
    for (i=idx+1; i<nq; i++) 
    {
        d2q = d[i]/q[i];
        temp = exp(-q[i] * x[i]);
        *price -= 
            ((k[i] - d2q)  * BSQ_INT(muMQ, sigMQ, 0,      x[i+1],   x[i]) +
             d2q * temp     * BSQ_INT(muMQ, sigMQ, q[i],   x[i+1],   x[i]));
    }

    /* apply affine adjustment and rescale by forward */
    *price *= K * C * fwdRate;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }
    
    return status;
} /* CRXQPricer */

/*f----------------------------------------------------------------------------
 * CRXQLevPricer
 *
 * MultiQ Black pricer with leverage     
 *
 * Uses Q parameters and the normal mean/stdev to compute call/put price
 */
int CRXQLevPricer(
    const CRXQDATA *mq,       /* (I) MQ data             */
    long    optType,  /* (I) option type         */
    double  strike,   /* (I) option strike       */
    double  leverage, /* (I) leverage            */
    double *price     /* (O) option price & vega */
    )
{
    static char routine[] = "CRXQLevPricer";
    int         status    = FAILURE;

    double lPrice, lStrike;
    long   lOptType;

    /* handle -ve leverage flips call <-> put */
    if (leverage < 0.)
    {
        if (optType == CRXQ_CALL)
        {
            lOptType = CRXQ_PUT;
        }
        else if(optType == CRXQ_PUT)
        {
            lOptType = CRXQ_CALL;
        }
        else
        {
            goto RETURN;
        }
    }
    else
    {
        lOptType = optType;
    }

    /* adjust strike */
    lStrike = strike / leverage;

    if (CRXQPricer(
        mq,
        lOptType,
        lStrike,
        &lPrice) == FAILURE) goto RETURN;

    *price = fabs(leverage) * lPrice;
        
    status = SUCCESS;
    
 RETURN:

    if (status == FAILURE)
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQLevPricer */


/*f----------------------------------------------------------------------------
 * CRXQBinPricer
 *
 * MultiQ binary pricer     
 *
 * Uses Q parameters and the normal mean/stdev to compute binary price
 */
int CRXQBinPricer(
    const CRXQDATA  *mq,      /* (I) MQ data  */
    long     optType, /* (I) call/put */
    double   strike,  /* (I) option strike */
    double  *price)   /* (O) option price & vega */
{
    static char     routine[] ="CRXQBinPricer";
    double          fwdRate   = mq->fwdRate;
    double          muMQ      = mq->muMQ;
    double          sigMQ     = mq->sigMQ;
    int             status    = FAILURE;

    double          strkMultMQ, strkMult, temp, xtemp, d2q;
    const double          *q, *k, *d;
    double* x;
    long            nq, idx;         
    int             sign;

    /* affine adjustment constants */
    double C = mq->C;
    double K = mq->K;

    /* check option type */
    if (optType != CRXQ_CALL && optType != CRXQ_PUT) 
    {
        DR_Error("%s failed: illegal option type - must be CRXQ_CALL or CRXQ_PUT\n", routine);
        goto RETURN;
    }

    /* check forward */
    if (fwdRate < CRXQ_MIN_FWD) 
    {
        DR_Error("%s failed: forward rate (%lf) too low\n", routine, fwdRate);
        goto RETURN;
    }

    /* branch into low/high strikes (left/right); subtle
     * choice needed in calibration: for ATM option, use left
     * intervals for put and high intervals for call */
    strkMult = strike/fwdRate;

    /* Invert affine map */
    strkMultMQ = (K-1)/(K*C) + strkMult/ (K*C);;

    if (strkMultMQ < 1.0) 
    {
        nq = mq->nbQL;
        q = mq->qL;
        k = mq->kL;
        x = (double*)mq->xL;
        d = mq->dL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        k = mq->kR;
        x = (double*)mq->xR;
        d = mq->dR;
        sign = 1;
    }

    /* Gaussian integration cutoff for the last value at "infinity" */
    x[nq] = muMQ + sign * (CRXQ_CUM_NUM_STDEV) * sigMQ;
      
    /* find strike location */
    idx = 0;
    while (idx < nq - 1 && sign * (strkMultMQ - k[idx+1]) > 0.) idx++;
    
    /* find normal point corresponding to strike */
    d2q = d[idx]/q[idx];
    temp = 1 + (strkMultMQ - k[idx])/d2q;
    xtemp = (temp > 0 ? x[idx] + log(temp)/q[idx] : x[nq]);

    /* for zero volatility, price intrinsic */
    if (sigMQ > CRXQ_MIN_VOL) 
    {
        *price = NormCum((xtemp - muMQ) / sigMQ);
    }
    else
    {
        *price = (strkMult >= 1.0 ? 1. : 0.);
    }

    if (optType == CRXQ_CALL) *price = 1. - *price;

    status = SUCCESS;

 RETURN:
    
    return status;

} /* CRXQBinPricer */

int CRXQCalibrateAffineAdjustments(CRXQDATA* mq) {

    static char   routine[] = "CRXQCalibrateAffineAdjustments";
    int           status    = FAILURE;
    double        expecH = 0.0; // expected value of H(x)
    int i; // loop counter
    double        atmCallPrice;  // the calculated ATM call price
    double d2q;

    mq->K = 1.0;
    mq->C = 1.0;

    /* Special case for zero vol: must make sure that H(muMQ) maps to the forward 
       Set volC=1.0, since degenerate case. */
    if (mq->sigMQ<TINY)
    {
        mq->C = 1.0;
        mq->K = 1.0;
        mq->optATM = 0.0;
        if (CRXQMap(mq, mq->muMQ, &expecH)==FAILURE)
        {
            DR_Error("%s failed: unable to calculate mapping of muMQ for zero vol calibration.\n", routine);
            goto RETURN;
        }
        mq->C = 1/expecH;
        status = SUCCESS;
        goto RETURN;
    }

    /* FIND THE EXPECTATION OF THE MAPPING FUNCTION, H(x) */
    /* Left side first */
    for (i=0; i<mq->nbQL; i++) {
        d2q = mq->dL[i]/mq->qL[i];
        expecH += (mq->kL[i] - d2q) * fabs(BSQ_INT(mq->muMQ, mq->sigMQ , 0.0, mq->xL[i], mq->xL[i+1]));
        expecH += d2q*exp(-mq->qL[i]*mq->xL[i]) * fabs(BSQ_INT(mq->muMQ, mq->sigMQ , mq->qL[i], mq->xL[i], mq->xL[i+1]));
    }
    /* Now right side */
    for (i=0; i<mq->nbQR; i++) {
        d2q = mq->dR[i]/mq->qR[i];
        expecH += (mq->kR[i] - d2q) * fabs(BSQ_INT(mq->muMQ, mq->sigMQ , 0.0, mq->xR[i], mq->xR[i+1]));
        expecH += d2q*exp(-mq->qR[i]*mq->xR[i]) * fabs(BSQ_INT(mq->muMQ, mq->sigMQ , mq->qR[i], mq->xR[i], mq->xR[i+1]));
    }
    mq->C = 1/expecH;

    /* Now, since E[Y] = fwdRate, (K (C expecH-1) + 1) = 1. Also, if fwdX is the point
       in gaussian X-space corresonding to the forward rate,
       fwdRate = fwdRate (K (C H(fwdX)-1) + 1)
       fwdX = H_inv(expecH)
       Since the price of the call is E[(Y-fwdRate)+] = fwdRate K C E[(H(x) - H(fwdX)+]
       we can price the call assuming that volC=1 and then rescale to get the correct call price
       */
    mq->K=1.0;
    if (CRXQPricer(mq, CRXQ_CALL, mq->fwdRate, &atmCallPrice)==FAILURE) 
    {
        DR_Error("%s failed: unable to price ATM option with forward rate %lf.\n", routine, mq->fwdRate);
        goto RETURN;
    }
    /* The target price is the simple BS ATM option price */
    if (CRXQBSQPricer(mq->fwdRate, mq->fwdRate, mq->optExpy, mq->sigATM, 1.0, CRXQ_CALL, &(mq->optATM))==FAILURE) 
    {
        DR_Error("%s failed: unable to price ATM option with forward rate=%lf, expiry=%lf, vol=%lf.\n", 
            routine, mq->fwdRate, mq->optExpy, mq->sigATM);
        goto RETURN;
    }
    //DR_Error("I have Q-call price=%lf and BS call price=%lf.\n", atmCallPrice, mq->optATM);
    //goto RETURN;
    /* Now set K so that the ATM option price matches, and the expectation is the
       forward rate */
    mq->K = (mq->optATM)/atmCallPrice;

    status = SUCCESS;
    RETURN:
  if (status == FAILURE) 
  {
    DR_Error("%s: failed\n", routine);
  }

  return status;
}


/* Calibrate a CRXQ data structure from an input set of variables */
int CRXQCalibrateCRXQDATA(double optionExpiry, 
                     double forwardRate, 
                     double volATM, 
                     int numberOfQs,
                     const double *inputQ,
                     const double *inputDelta,    // middle value to both L and R, same 
                     CRXQDATA* mq
                     ) {

    static char      routine[] = "CRXQCalibrateCRXQDATA";
    int              status    = FAILURE;
    double           sigTotal  = volATM * sqrt(optionExpiry);
    int i; // loop counter

    /* SOME SIMPLE CHECKS FIRST */
    if (optionExpiry<TINY)
    {   
        DR_Error("%s failed: non-positive time to expiry.\n", routine);
        goto RETURN;
    } else if (volATM<TINY) 
    {
        DR_Error("%s failed: non-positive ATM vol.\n", routine);
        goto RETURN;
    } else if (numberOfQs<=0 || numberOfQs%2!=0) 
    {
        DR_Error("%s failed: You must provide an even number of q-parameters.\n", routine);
        goto RETURN;
    }


    /* need to check that the order of deltas is correct and no deltas are beyond the cutoff */
    for (i=0; i < numberOfQs-2; i++)
    {
        if (inputDelta[i+1]-inputDelta[i] < TINY)
        {
          DR_Error("%s failed: Deltas %lf [%d] and  %lf [%d] are too close together for no. qs %d.\n", 
              routine, inputDelta[i], i, inputDelta[i], (i+1), numberOfQs);
          goto RETURN;
        }
    }
    /* note: add checks for delta cutoff and more resonable spacing */

    /* set-up */
    /* misc set-up */
    mq->nbQL = numberOfQs/2;
    mq->nbQR = mq->nbQL;
    mq->fwdRate   = forwardRate;
    mq->sigATM    = volATM;
    mq->optExpy   = optionExpiry;

      // these are chosen. The calibration comes from the volC and fwdC params
      mq->sigMQ     = sigTotal;
      mq->muMQ      = 0.0; // this is by choice: 

    /* LEFT SIDE */
    for(i = 0;i < mq->nbQL; i++)
          /* Reset all the x values to match the strike boundaries, deltaL.   */
    {
        mq->qL[i] = inputQ[mq->nbQL - i - 1];
        /* If the q's are zero, set them to a manageable small number */
        if (fabs(mq->qL[i])<CRXQ_Q_SHIFT) mq->qL[i] = CRXQ_Q_SHIFT;
        mq->xL[i] = sigTotal * NormCumInv(inputDelta[mq->nbQL - i - 1]);
     }
    mq->xL[mq->nbQL] = - sigTotal * CRXQ_CUM_NUM_STDEV; // cut off beyond this

  /* Build k and d coefficients using freedom to set first values 
   * so that H(0)=1 and H'(0)=1 and then
   * by requiring continuity and differentiability */
  mq->kL[0] = 1.0;
  mq->dL[0] = 1.0;
  for(i = 1;i < mq->nbQL;i++)
  {
    mq->dL[i] = mq->dL[i-1]*exp(mq->qL[i-1] * (mq->xL[i]-mq->xL[i-1]));
    mq->kL[i] = mq->kL[i-1] + mq->dL[i-1] * qxq(mq->qL[i-1], mq->xL[i], mq->xL[i-1]);
  }  
   
  /* RIGHT SIDE - SIMILAR TO LEFT, SO NO COMMENTS */
  for(i = 0;i < mq->nbQR;i++)
  {
    mq->qR[i] = inputQ[mq->nbQR + i];
    /* If the q's are zero, set them to a manageable small number */
    if (fabs(mq->qR[i])<CRXQ_Q_SHIFT) mq->qR[i] = CRXQ_Q_SHIFT;
    mq->xR[i] = sigTotal * NormCumInv(inputDelta[mq->nbQR -1 + i]);
  }
  mq->xR[mq->nbQR] = sigTotal * CRXQ_CUM_NUM_STDEV;
  mq->kR[0] = 1.0;
  mq->dR[0] = 1.0;
  for(i = 1;i < mq->nbQR;i++)
  {
    mq->dR[i] = mq->dR[i-1]*exp(mq->qR[i-1] * (mq->xR[i]-mq->xR[i-1]));
    mq->kR[i] = mq->kR[i-1] + mq->dR[i-1] * qxq(mq->qR[i-1], mq->xR[i], mq->xR[i-1]);
  } 
  
  /* Now you must calibrate the affine adjustment parameters so that the expected value
   * is the forward, and the ATM option price matches BS */
  if (CRXQCalibrateAffineAdjustments(mq)==FAILURE) 
  {
      DR_Error("%s failed: unable to calculate volC and fwdC parameters in Q-calibration.", 
          routine);
      goto RETURN;
  }
 
  status = SUCCESS;
    

    RETURN:
  if (status == FAILURE) 
  {
    DR_Error("%s: Failed\n", routine);
  }

  return status;
}



/*f----------------------------------------------------------------------------
 * MultiQ: mapping function
 *
 * Evaluates mapping function at normal point x given CRXQDATA data
 */
int CRXQMap(
    const CRXQDATA      *mq,      /* (I) MQ data       */
    double             xval,    /* (I) normal point  */
    double            *yield)   /* (O) mapped yield  */
{
    static char      routine[]="CRMQMap";
    int              status = FAILURE;

    /* distribution data */

    double fwdRate = mq->fwdRate;

    long   nq; /* number of intervals      */
    const double *q; /* pointer to q parameters  */
	const double *d; /* pointer to d parameters  */
    const double *k; /* pointer to strike ratios */
	const double *x; /* pointer to gaussian space strike boundaries */
    
    long   idx;
    int    sign;

    /* check forward */
    if (fwdRate < CRXQ_MIN_FWD) 
    {
        DR_Error("%s: fwdRate = %lf is too small.\n", routine, fwdRate);
        goto RETURN;
    }
	
    /* for zero volatility, return forward no matter what the point is */
    if (mq->sigMQ < CRXQ_MIN_VOL) 
    {
        *yield = fwdRate;
        status = SUCCESS;
        goto RETURN;
    }    

    /* branch into low/high strikes (left/right) */
    if (xval <= mq->xL[0])
	{
		nq = mq->nbQL;
		q = mq->qL;
		d = mq->dL;
		k = mq->kL;
		x = mq->xL;
		sign = -1;
	} 
	else 
	{
		nq = mq->nbQR;
		q = mq->qR;
		d = mq->dR;
		k = mq->kR;
		x = mq->xR;
		sign = 1;
	}
	
    /* find x location */
    idx = 0;
    while (idx < nq - 1 && sign * (xval - x[idx+1]) >= 0.) idx++;
	
    /* evaluate mapping function */

	*yield = k[idx] + d[idx] * (exp(q[idx]*(xval-x[idx])) - 1)/q[idx];

	*yield = mq->fwdRate*(mq->K*(mq->C * *yield - 1.0) + 1.0);
	
    status = SUCCESS;
 RETURN:
    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }
    
    return status;
} /* CRXQMap */

/*-----------------------------------------------------------------------------
 * CRXQGridPricer
 *
 * Uses numerical grid integration to price a payoff on a multi-q distribution.
 *
 */
int CRXQGridPricer(
    CRXQPAYOFF  *optPayoff, /* (I) option payoff structure */
    CRXQFPAYOFF *payFunc,   /* (I) option payoff function  */
    CRXQDATA  *mq,        /* (I) measure data            */
    double  *premium    /* (O) fwd premium/rate        */
    )
{
    double grid[CRXQ_1D_INT_PTS];
    double yield[CRXQ_1D_INT_PTS];
    double dens[CRXQ_1D_INT_PTS];
    double point[2];
    double densNorm;
    double fwdRate;

    int status = FAILURE;

    /* check for degenerate case */
    if (mq->sigMQ > CRXQ_MIN_VOL)
    {
        /* calculate MQ yield density */
        if (CRXQDens(            
            mq,
            -CRXQ_1D_NUM_STDEV,   /* start */
            CRXQ_1D_NUM_STDEV,    /* end   */
            CRXQ_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;

        /* integrate density * payoff (from q3.h) */
        if (CRXQSimpsonPricer1D(
            optPayoff,
            payFunc,
            grid,
            yield,
            dens,
            densNorm,
            premium) == FAILURE) goto RETURN;
    }
    else
    {   
        fwdRate = mq->fwdRate;
        if( mq->calcFwd == TRUE)
        {
            if(CRXQMap(
               mq,
               mq->muMQ,
               &fwdRate) == FAILURE) goto RETURN;
        }
        point[0] = fwdRate;
        point[1] = mq->muMQ;
        if ((*payFunc)(
            optPayoff,
            point,
            0., 
            premium) == FAILURE) goto RETURN;
    }
    
    status = SUCCESS;

 RETURN:

    return status;

} /* CRXQGridPricer */


/*f----------------------------------------------------------------------------
 * CRXQDens
 *
 * Generates density for multi-q distribution from underlying measure.
 *
 */
int CRXQDens(
    const CRXQDATA *mq,         /* (I) measure */
    double  start,      /* (I) start in std normal space */
    double  end,        /* (I) end in std normal space   */
    long    numGridPts, /* (I) num pts to generate       */
    double *grid,       /* (O) normal grid points        */
    double *yields,     /* (O) yield points              */
    double *dens,       /* (O) density                   */
    double *normC       /* (O) density normalization     */
    )
{
    static char routine[] = "CRXQDens";
    int         status    = FAILURE;

    double step, xval;
    long   i;

    step  = (end - start)/(numGridPts-1);
    xval  = start - step;

    /* calculate mapping and density */
    for (i = 0; i < numGridPts; i++)
    {
        xval   += step;
        grid[i] = mq->sigMQ * xval + mq->muMQ;
        if (CRXQMap(
            mq,
            grid[i],
            &(yields[i])) == FAILURE) goto RETURN;
             
        dens[i] = NormDens(xval);
    }

    /* normalize distribution */
    *normC = Q3SimpsIntegral(
        step,
        numGridPts,
        dens);
    if (*normC < TINY) 
    {
        DR_Error("%s: Non-positive total mass.\n", routine);
        goto RETURN;
    }

    status = SUCCESS;

 RETURN:

    return status;
} /* CRXQDens */

/* Does one Newton iteration on f(q) = [exp(q.deltaX)-1]/q for target value targetF. 
    f'(q) = (deltaX - 1/q)f(q) */
int qSolveNewton1Iteration(double* q1, double q0, double targetF, double deltaX) {
    if (fabs(q0)<TINY) 
    {
        // for small q, Taylor series is f(q) ~ deltaX + 0.5 q deltaX^2
        // f'(q) ~ 0.5 deltaX^2
        *q1 = q0 + 2*(targetF - deltaX - 0.5*q0*deltaX*deltaX)/(deltaX*deltaX);
        return SUCCESS;
    } 
    else 
    {
        // f'(q) = deltaX exp(q.deltaX)/q - (exp(q.deltaX) - 1)/q^2
        double expVal = exp(q0*deltaX);
        double fVal = (expVal-1)/q0;
        double fPrimeVal = (deltaX*expVal - fVal)/q0;
        if (fabs(fPrimeVal)<TINY) {
            return FAILURE;
        }
        else 
        {
            *q1 = q0 + (targetF - fPrimeVal)/fPrimeVal;
            return SUCCESS;
        }
    }
}


/******************************************************************************
 * CRXQApproximateFAWithMQ
 * Calibrate a CRXQDATA multi-q measure to approximate a forward-adjusted 
 * measure (used for faster bivariate pricing)
 * TODO: Add more details about how this is done
 *       Change so that uses full number of q intervals - CRXQ_NBQ
 *****************************************************************************/ 
int CRXQApproximateFAWithMQ(
    const CRXQFADATA* fa,   /**<(I) Given FA measure to be approximated      */
    int               numQ, /**<(I) No. of qs on each side of fwd to fit     */
    CRXQDATA*         pa    /**<(O) Resulting multi-q measure                */    
    )             
{

    static char routine[] = "CRXQApproximateFAWithMQ";
    int         status    = FAILURE;

    CRXQPAYOFF pf; 
    /* Numerical FA density */
    double grid[CRXQ_1D_INT_PTS];
    double yield[CRXQ_1D_INT_PTS];
    double dens[CRXQ_1D_INT_PTS];
	double cdfDens[CRXQ_1D_INT_PTS];
    double densNorm;

    int ptrK;
    double thisK;
    double thisY, lastY; // adjacent yield values for CDF interp
    double thisP, lastP; // adjacent P-values for CDF interp
    double deltaX;
    double fQTarget;
    double qNext;
    double qLast;
    double* k = 0;
    double* x = 0;
    double* d = 0;
    double* q = 0;
    double* cdf = 0;
    double cdfL[CRXQ_NBQ];
    double cdfR[CRXQ_NBQ];
    int nbIT; // number of iterations in Newton's method for q
    int    i, side, nbQ;
    double strikeStepL;
    double strikeStepR;

    if (numQ<0 || numQ>CRXQ_NBQ) 
    {
        DR_Error("%s failed: illegal numQ, %d, should be positive and <=%d.\n", 
            routine, numQ, CRXQ_NBQ);
    }

    /* calculate the FA yield density                              */
    /* note: use CRXQFADens & Q3SimpsonPricer1D separately instead */
    /* of CRXQFAPricer for efficiency from reuse of density        */
    if (CRXQFADens(            
        fa,
        -CRXQ_1D_NUM_STDEV,   /* start */
        CRXQ_1D_NUM_STDEV,    /* end   */
        CRXQ_1D_INT_PTS,
        &(grid[0]),
        &(yield[0]),
        &(dens[0]),
        &densNorm) == FAILURE) 
    {
        DR_Error("%s failed: unable to compute numerical FA density.\n", 
            routine);
        goto RETURN;
    }

	/* calculate fwd and set in new MQ */
    if (CRXQSimpsonPricer1D(
        &pf,
        CRXQPay1D_Yield,
        grid,
        yield,
        dens,
        densNorm,
        &(pa->fwdRate)) == FAILURE) goto RETURN;

	if (pa->fwdRate < TINY) 
    {
        DR_Error("%s failed: Adjusted forward: %f < TINY.\n", 
            routine, pa->fwdRate);        
        goto RETURN;
    }
    // intervals between adjacent k points equal spaced over yield range,
    // kR[0]=kL[0] = 1.0;
    strikeStepL = (1.0-yield[0]/pa->fwdRate)/CRXQ_NBQ;
    strikeStepR = (yield[CRXQ_1D_INT_PTS]/pa->fwdRate - 1.0)/CRXQ_NBQ;
    
    /* Calculate the ATM call price and implied BS vol for the FA measure */
    pf.cop = 1; // call
    pf.strike = pa->fwdRate;
    if (CRXQSimpsonPricer1D(
        &pf,
        CRXQPay1D_Vnl,
        grid,
        yield,
        dens,
        densNorm,
        &(pa->optATM)) == FAILURE) 
    {
        DR_Error("%s failed: unable to price ATM call under FA measure", routine);
        goto RETURN;
    }
    /* Compute the BS implied vol for the ATM option (BSQImpVol from Q3/bas.c) */
    if (FAILURE==BSQImpVol(
        pa->fwdRate,
        pa->fwdRate, 
        fa->mq->optExpy, 
        1.0, pa->optATM, 
        Q3_CALL, 
        fa->mq->sigATM, 
        &(pa->sigATM)))
    {
        DR_Error("%s failed: unable to calculate implied vol for ATM option with price %lf.\n", 
            routine, pa->optATM);
        goto RETURN;
    }

	/* Compute the (unnormalized) cumulative probabilities that you will be matching 
       and recalculate the normalization constant. */
	cdfDens[0] = 0.0;
    for (i=1; i<CRXQ_1D_INT_PTS; i++) {
		cdfDens[i] = cdfDens[i-1] + (yield[i]-yield[i-1]) * 0.5 * (dens[i] + dens[i-1]);
	}
    /* Used a slightly different integration method from CRXQFADens,
     * so you need to make sure that you match the computed CDF */
    densNorm = cdfDens[CRXQ_1D_INT_PTS-1];

    /* Copy other parameters that you know will stay
       the same in the new measure */
    pa->optExpy = fa->mq->optExpy;
    pa->sigMQ = pa->sigATM;
    pa->muMQ = 0.0;
	pa->nbQL = CRXQ_NBQ; // use full number of Qs on each side to
	pa->nbQR = CRXQ_NBQ; // make approximation as accurate as possible

	/* Initialize strike points */
    pa->kL[0] = 1.0;
    pa->kR[0] = 1.0;
    pa->dL[0] = 1.0;
    pa->dR[0] = 1.0;
    pa->qR[CRXQ_NBQ - 1] = 0.0; // normal in upper tail
    pa->qL[CRXQ_NBQ - 1] = 1.0; // lognormal in lower tail
	for (i=1; i<CRXQ_NBQ; i++) {
		pa->kL[i] = pa->kL[i-1] - strikeStepL;
		pa->kR[i] = pa->kR[i-1] + strikeStepR;
	}

	// no affine adjustments until calibrate at last stage
	pa->K = 1.0;
	pa->C = 1.0;

    /* skip further calibration if vol too small - all you can do 
       is match the forward */
    if (pa->sigMQ < CRXQ_MIN_VOL_CALIB)
    {
        pa->sigMQ = 0.;
        pa->muMQ = 0.;
        status = SUCCESS;
        goto RETURN;
    }

    /* Find all the strike boundary cumulative probabilities */
    // LEFT SIDE
    ptrK = pa->nbQL-1;
    thisK = pa->kL[ptrK];
    lastY = 0.0;
    lastP = 0.0;
    i = 0;
    while (i<CRXQ_1D_INT_PTS && ptrK>=0) {
        thisY = yield[i];
        thisP = cdfDens[i];
        if (thisY>=pa->fwdRate*thisK) {
            // you have passed kL[ptrK] - interpolate cdfL
            if (fabs(lastP)<TINY || fabs(thisY-lastY)<TINY) 
            {
                // lastP too small to bother with, first step
                // or y-steps too close together
                cdfL[ptrK] = thisP;
            } else 
            {
                cdfL[ptrK] = (thisP*(pa->fwdRate*thisK - lastY) + lastP*(thisY - pa->fwdRate*thisK))/(thisY-lastY);
            }
            ptrK--;
        }
        lastY = thisY;
        lastP = thisP;
        i++;
    }
    // RIGHT SIDE
    ptrK = 1;
    cdfR[0] = cdfL[0];
    thisK = pa->kL[ptrK];
    while (i<CRXQ_1D_INT_PTS && ptrK<pa->nbQR) {
        thisK = yield[i];
        thisP = cdfDens[i];
        if (thisY>=pa->fwdRate*thisK) {
            // you have passed kL[ptrK] - interpolate cdfL
            if (fabs(lastP)<TINY || fabs(thisY-lastY)<TINY) 
            {
                // lastP too small to bother with,
                // or y-steps too close together
                cdfL[ptrK] = thisP;
            } else 
            {
                cdfL[ptrK] = (thisP*(pa->fwdRate*thisK - lastY) + lastP*(thisY - pa->fwdRate*thisK))/(thisY-lastY);
            }
            ptrK++;
        }
        lastY = thisY;
        lastP = thisP;
        i++;
    }

    /* Now calibrate the x, q, d values to match the CDF */
	for (side = 0; side<2; side++) {
		if (side==0) {
			/* LEFT SIDE */
			nbQ = fa->mq->nbQL;
			x = fa->mq->xL;
			q = fa->mq->qL;
			d = fa->mq->dL;
            k = fa->mq->kL;
            cdf = cdfL;
		} else {
			/* RIGHT SIDE */
            nbQ = fa->mq->nbQR;
			x = fa->mq->xR;
			q = fa->mq->qR;
			d = fa->mq->dR;
            k = fa->mq->kR;
            cdf = cdfR;
		}

        for (i=0; i<nbQ; i++) {

            /* x[i] is determined by the CDF:
               H(x[i]) = k[i], and so
               CDF[k[i].fwdRate] = Phi((x[i]-muMQ)/sigMQ) */
            x[i] = NormCumInv(cdf[i]/densNorm)*pa->sigMQ + pa->muMQ;
            if (fabs(d[i])<TINY)
            {
                DR_Error("%s failed: zero d-value at i=%d on %c side.", routine, i, (side==0 ? 'L' : 'R'));
                goto RETURN;
            }
            /* Note that d[0] is set already */
            if (i>0) {
                // calibrate d, q to match x and k:
                deltaX = x[i] - x[i-1];
                fQTarget = (k[i] - k[i-1])/d[i];
                nbIT = 0;
                qNext = q[i-1]; // start with the old q as your first guess
                qLast = qNext+1; // guarantee 1 iteration at least
                while (fabs(qNext-qLast)>=1E-7 /*TODO: WHAT TOLERANCE SHOULD I USE?*/)
                {
                    nbIT++;
                    qLast = qNext;
                    if(FAILURE==qSolveNewton1Iteration(&qNext, qLast, fQTarget, deltaX))
                    {
                        DR_Error("%s failed: zero derivative in Newton's method solving for q at %lf after %d iterations with target %lf and deltaX %lf", 
                            routine, qLast, nbIT, fQTarget, deltaX);
                        goto RETURN;
                    }
                    else if (nbIT>30)
                    {
                        DR_Error("%s failed: no convergence after 30 iterations with last q %lf after %d iterations with target %lf and deltaX %lf",
                            routine, qLast, nbIT, fQTarget, deltaX);
                    }
                }
                /* Note that this sets q[i-1]: q[nbQ-1] will stay the same as the
                   tail-q in the original distribution, as it was initialized that way above */
                q[i-1] = qNext;
                
                /* Differentiability of H(x) requires d[i] = d[i-1].exp(q[i-1].deltaX) */
                d[i] = d[i-1] * exp(q[i-1]*deltaX);
                    
            }   
        }
	}

    /* Finally, calibrate C and K to match the ATM call price and the forward. */
    if (FAILURE==CRXQCalibrateAffineAdjustments(pa)) 
    {
        DR_Error("%s failed: unable to calibrate affine parameters C, K to match forward and ATM call price.");
        goto RETURN;
    }
    
    status = SUCCESS;
 RETURN:
    return status;

} /* CRXQApproximateFAWithMQ */



