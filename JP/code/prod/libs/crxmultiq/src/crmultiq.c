/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File: multiq.c       
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "crmultiq.h"
#include "mqdata.h"

#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    DR_Error("%s: Required condition (%s) fails!",routine,#cond);\
    goto RETURN;\
}} while(0)
 
/*----------------------------------------------------------------------------
 * return (exp(q*(x-x0))-1)/q
 *
 -----------------------------------------------------------------------------*/
static double qxq(double q, double x,double x0)
{
  if(fabs(q) < TINY) {
    return x-x0;
  }
  else {
    return (exp(q*(x-x0))-1.0)/q;
  }
}

/*----------------------------------------------------------------------------
 * MultiQ Black pricer     
 *
 * Uses Q parameters and the normal mean/stdev to compute call/put price
 * *mq must be initialized and updated beforehand
 * !!!C _MUST_ be computed beforehand to match the forward; otherwise, put-call
 * parity used to compute low puts and high calls will be incorrect !!!
 -----------------------------------------------------------------------------*/
int Q3MQPricerCR(
                 MQDATA            *mq,                /* (I) MQ data  */
                 long              optType,            /* (I) option type */
                 double            strike,             /* (I) option strike */
                 double            *price)             /* (O) option price & vega */
{
  static char      routine[]="Q3MQPricerCR";
  int              status = FAILURE;
 
  double fwdRate  = mq->fwdRate;
  double sigTotal   = mq->sigMQ * sqrt(mq->optExpy);
  
  long   nq;          /* number of intervals                    */
  double strkMult;    /* strike/fwdRate                         */
  double *q;          /* pointer to q parameters                */
  double *b;          /* pointer to b parameters                */
  double *k;          /* pointer to strike ratios               */
  double *x;
  double *d;
  double CK;

  double temp, xtemp, b2q;
  long   i, idx;
  int    sign;


  /* check vol */
  if(sigTotal < 0){
    DR_Error("%s: variance is negative", routine);
    goto RETURN;
  }

  /* check option type */
  if (optType != Q3_CALL && optType != Q3_PUT) 
  {
    DR_Error("%s: invalid option type", routine);
    goto RETURN;
  }

  /* write a special case for zero total vol! */
  if (sigTotal < Q3_MIN_VOL)
  {
    *price = MAX(Q3_COP_TO_COP(optType) * (fwdRate - strike), 0);
    status = SUCCESS;
    goto RETURN;
  }


  /* branch into low/high strikes (left/right); subtle      */
  /* choice needed in calibration: for ATM option, use left */
  /* intervals for put and high intervals for call          */
  strkMult = strike/fwdRate;

  CK = mq->C*mq->K;
  strkMult = strkMult - (1 - mq->K);
  /* now strkMult is the strike multiple for CK*H(X) */

  if (strkMult < CK) 
  {
    nq = mq->nbQL;
    q = mq->qL;
    b = mq->bL;
    d = mq->dL;
    k = mq->kL;
    x = mq->xL;
    sign = -1;
  } 
  else 
  {
    nq = mq->nbQR;
    q = mq->qR;
    b = mq->bR;
    d = mq->dR;
    k = mq->kR;
    x = mq->xR;
    sign = 1;
  }


  /* find strike location */
  idx = 0;
  while (idx < nq - 1 && sign * (strkMult - CK * k[idx+1]) >= 0.) idx++;

  /* find normal point corresponding to strike */

  if(fabs(q[idx]) < TINY) /* q[idx]=0, normal dist int */     
  {         
    temp = (strkMult/CK - k[idx])/b[idx];         
    xtemp = x[idx] + temp;        
  } 
  else 
  {
    b2q = b[idx]/q[idx];         
    temp = 1 + (strkMult/CK - k[idx])/b2q;         
    xtemp = (temp > 0 ? x[idx] + log(temp)/q[idx] :x[nq]);     
  }

  /* key step(!): price put for low strikes, call for high strikes */

  /* first interval */
  if((sign == 1 && xtemp >= x[idx+1])
     || (sign == -1 && xtemp <= x[idx+1]))
  {
    *price = 0.0;
  } 
  else 
  {
    *price = 
      BSQIntCR(xtemp,x[idx+1],k[idx], b[idx], q[idx], x[idx], 
               CK, sigTotal) - 
      strkMult * (NormCum(x[nq]/sigTotal) - NormCum(xtemp/sigTotal));
  }

  /* other intervals */
  for (i=idx+1; i<nq; i++) 
  {
    *price += 
      BSQIntCR(x[i],x[i+1],k[i], b[i], q[i], x[i], 
               CK, sigTotal);
  }

  /* we have priced call for low strikes and put for high strikes 
     (the latter includes ATM); if the goal is to price a high call or
     low put, need to use put-call parity. 
     !!! note that after the shift, strkMult is the strike multiple for CK*H(X), 
     not for the original model CK*H(X)-K+1; since C=1/E(H(X)), the
     forward of CK*H(X) is K, not 1 */
  if (strkMult < CK) 
  {
    /* put-call parity for low strikes */
    if (optType == Q3_CALL){
      *price += (mq->K - strkMult);
    }
  } 
  else 
  {
    /* put-call parity for high strikes */
    if (optType == Q3_PUT)  *price -= (mq->K - strkMult);
  }

  /* rescale by fwd */
  *price *= fwdRate;
  status = SUCCESS;

 RETURN:

  if (status == FAILURE) 
  {
    DR_Error("%s: Failed\n", routine);
  }

    
  return status;
} /* Q3MQPricer */



/*
	Risky discount factor given par spread
	For now, powerDel == volRatio
*/
double Zero(double parSpread, FADATA *fa)
{
	return pow(1+ fa->alphaDel * pow(parSpread, fa->powerDel)/fa->freqAnn, -fa->matDel*fa->freqAnn ); 
} /* Zero */

/* 
	Analytic approximation to the price of protection given a CDS spread and an
	equivalent interest rate, which makes the approximation correct for the 
	case where the input spread is the forward par spread
*/
double Protection(double parSpread, double IRate, FADATA *fa)
{
	double x;
	double l = fa->alphaAnn*pow(parSpread, fa->powerAnn); /* parSpread*fa->accrFactor/(1-fa->recovery);*/

	x = l/(l+IRate) * 
					(1 - pow(1+fa->alphaAnn*pow(parSpread, fa->powerAnn)/fa->freqAnn, -fa->matAnn*fa->freqAnn )
					* pow(1+IRate/fa->freqSwap, -fa->matAnn*fa->freqSwap));  
	return  x;
} /* Protection */


/*------------------------------------------------------------------------------
 * update mq coefficients
 *
 -----------------------------------------------------------------------------*/
int Q3MQUpdateCR(MQDATA     *mq)     /* (O) MQDATA structure             */
{
  static char      routine[] = "Q3MQUpdateCR";
  int              status    = FAILURE;
  int              i;
  double           sigTotal  = mq->sigMQ * sqrt(mq->optExpy);

  
  /* left side */
  for(i = 0;i < mq->nbQL;i++)
  {
    mq->xL[i] = sigTotal * NormCumInv(mq->dL[i]);
  }
  mq->xL[mq->nbQL] = - sigTotal * Q3_CUM_NUM_STDEV;
  
  mq->kL[0] = 1.0;
  mq->bL[0] = 1.0;
  for(i = 1;i < mq->nbQL;i++)
  {
    mq->bL[i] = mq->bL[i-1]*exp(mq->qL[i-1] * (mq->xL[i]-mq->xL[i-1]));
    mq->kL[i] = mq->kL[i-1] + mq->bL[i-1] * qxq(mq->qL[i-1], mq->xL[i], mq->xL[i-1]);
  }  
   
  /* right side */
  for(i = 0;i < mq->nbQR;i++)
  {
    mq->xR[i] = sigTotal * NormCumInv(mq->dR[i]);
  }
  mq->xR[mq->nbQR] = sigTotal * Q3_CUM_NUM_STDEV;

  mq->kR[0] = 1.0;
  mq->bR[0] = 1.0;
  for(i = 1;i < mq->nbQR;i++)
  {
    mq->bR[i] = mq->bR[i-1]*exp(mq->qR[i-1] * (mq->xR[i]-mq->xR[i-1]));
    mq->kR[i] = mq->kR[i-1] + mq->bR[i-1] * qxq(mq->qR[i-1], mq->xR[i], mq->xR[i-1]);
  }  

  if(Q3MQSetCandKCR(mq) == FAILURE) goto RETURN;
 
  status = SUCCESS;
  
 RETURN:

  if (status == FAILURE) 
  {
    DR_Error("%s: Failed\n", routine);
  }

  return status;
}

MQDATA* Q3MQMakeCR(
    double     fwd,
    double     vol,
    long       optType,
    double     time,
    CrxTQDist *qdist)
{
    static char routine[] = "Q3MQMakeCR";
    int         status    = FAILURE;

    int q3OptionType;
    int i;
    int nQs;
    double q[6];
    double d[6];
    
    MQDATA *mq = CrxMqdataMakeEmpty(3,3);
    if (mq == NULL)
       goto RETURN;

    switch (optType)
    {
    case CRX_OPTION_TYPE_CALL: q3OptionType = Q3_PUT; break;
    case CRX_OPTION_TYPE_PUT:  q3OptionType = Q3_CALL; break;
    default:
        DR_Error ("%s: Bad option type %ld", routine, optType);
        goto RETURN;
    }

    REQUIRE (qdist != NULL);
    REQUIRE (qdist->nQs == 6);
    REQUIRE (qdist->nDs == 5);

    nQs = qdist->nQs;
    for (i = 0; i < nQs-1; ++i)
    {
        q[i] = 1.0 - qdist->Qs[i];
        d[i] = qdist->Ds[i];
    }
    q[nQs-1] = 1.0 - qdist->Qs[nQs-1];
    d[nQs-1] = 1.0; /* this value should be ignored */

    if (Q3MQInitCR(mq,
                   fwd,
                   vol,
                   q3OptionType,
                   time,
                   nQs,
                   q,
                   nQs-1,
                   d) != SUCCESS)
    {
        DR_Error("%s: failed to initialise MQ structure", routine);
        goto RETURN;
    };

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        CrxMqdataFree (mq);
        mq = NULL;
        DR_Error("%s: Failed", routine);
    }

    return mq;
}
    
    


        

/*------------------------------------------------------------------------------
 * init mq for CR
 * 
 -----------------------------------------------------------------------------*/
int Q3MQInitCR(MQDATA     *mq,            /* (O) MQDATA structure             */
               double     fwdRate,        /* (I) fwd rate                     */
               double     sigATM,         /* (I) atm vol                      */
               long       optType,        /* (I) atm opt type                 */
               double     optExpy,        /* (I) atm opt expiry               */
               int        numQ,           /* (I) number of q's                */
               double     *q,             /* (I) a list of q's                */
               int        numD,           /* (I) number of d's                */
               double     *d)             /* (I) a list of d's                */
{
  int              i;
  static char      routine[] = "Q3MQInitCR";
  int              status    = FAILURE;

  if(numQ != 6 )
  {
    DR_Error("Number of Qs (%d) should be 6.\n",numQ);
    goto RETURN;
  }
  if(numD != 5)
  {
    DR_Error("Number of Ds (%d) should be 6.\n",numD);
    goto RETURN;
  }

  /* need to check that the order of deltas is correct 
     and no deltas are beyond the cutoff */
  for (i=0; i < numD-1; i++)
  {
    if (d[i+1]-d[i] < TINY)
    {
      DR_Error("Deltas are too close.\n");
      goto RETURN;
    }
  }
  /* note: add checks for delta cutoff and replace TINY by a minimal
     reasonable distance between deltas */

  /* misc */
  mq->fwdRate   = fwdRate;
  mq->sigATM    = sigATM;
  mq->optExpy   = optExpy;
  mq->optType   = Q3_CALL; 
  mq->calibType = Q3_MQ_TOL;

  mq->sigMQ     = sigATM;

  if (BSPricer(&mq->optATM,
               OPT_CALL,
               fwdRate,
               fwdRate,
               optExpy,
               sigATM,
               OPT_PRICE) == FAILURE)
    goto RETURN;
  


  /* left side */
  mq->nbQL = numQ/2;

  for(i = 0;i < mq->nbQL;i++)
  {
    mq->qL[i] = q[2-i];
    mq->dL[i] = d[2-i];
  }

  /* right side */
  mq->nbQR = numQ/2;
  for(i = 0;i < mq->nbQR;i++)
  {
    mq->qR[i] = q[3+i];
    mq->dR[i] = d[2+i];
  }

  // init, necessary
  mq->K = 1.0;

  // update misc coef
  if(Q3MQUpdateCR(mq) == FAILURE)
  {
    goto RETURN;
  }

  mq->optType = optType;

  status = SUCCESS;

 RETURN:

  if (status == FAILURE) 
  {
    DR_Error("%s: Failed\n", routine);
  }

  return status;
  
}

/*------------------------------------------------------------------------------
 * calculte E[H(x)] given q and vol, mq->C = 1/E[H(x)]
 *
 -----------------------------------------------------------------------------*/
int Q3MQSetCCR(MQDATA   *mq)       /* (I) MQ structure                  */
{
  int              i;
  static char      routine[] = "Q3MQSetCCR";
  int              status    = FAILURE;

  double           res;
  double           sigTotal  = mq->sigMQ * sqrt(mq->optExpy);

  if(!mq)
  {
    DR_Error("Credit MQ structure is null.");
    goto RETURN;
  }

  // special case
  if(fabs(sigTotal) < TINY)
  {
	mq->C = 1.0;
	status = SUCCESS;
	goto RETURN;
  }

  res = 0.0;  
   
  // left side
  for(i = mq->nbQL-1;i >= 0;i--)
  {
    res += BSQIntCR(
                    mq->xL[i+1],
                    mq->xL[i],
                    mq->kL[i],
                    mq->bL[i],
                    mq->qL[i],
                    mq->xL[i],
                    1.0, 
                    sigTotal);
  }
  

  // right side
  for(i = 0;i < mq->nbQR;i++)
  {
    res += BSQIntCR(
                    mq->xR[i],
                    mq->xR[i+1],
                    mq->kR[i],
                    mq->bR[i],
                    mq->qR[i],
                    mq->xR[i],
                    1.0, 
                    sigTotal);
  }

  if (fabs(res)<TINY) 
  {
    DR_Error("%s: failed in the normalization (residul = %lf)\n", routine, res);
    goto RETURN;
  }

  mq->C = 1.0/res;

  status = SUCCESS;

 RETURN:

  if (status == FAILURE) 
  {
    DR_Error("%s: Failed\n", routine);
  }

  return status;

} /* Q3MQSetCCR */
  
/*------------------------------------------------------------------------------
 * calculate C and then K to match forward and ATM vol
 *
 -----------------------------------------------------------------------------*/
int Q3MQSetCandKCR(MQDATA   *mq)    /* (I) MQ structure                  */
{
  static char      routine[] = "Q3MQSetCandKCR";
  int     status    = FAILURE;
  double  price     = 0.0;


  // set C first
  Q3MQSetCCR(mq);

  mq->K=1.0;

  // initialization
  if(Q3MQPricerCR(
                  mq,
                  mq->optType,
                  mq->fwdRate,
                  &price) == FAILURE)
    goto RETURN;

  if(price < 0) goto RETURN;

  if (fabs(price) < TINY) {
	
	if(fabs(mq->optATM) < TINY){
	  mq->K = 1;
	} else {
      DR_Error("%s: price = %f, optATM=%f\n", routine, price, mq->optATM);
	  goto RETURN;
	}
  } else {
	mq->K = mq->optATM/price;
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
 * Evaluates mapping function at normal point x given MQDATA data
 */
int Q3MQMapCR(
    MQDATA            *mq,      /* (I) MQ data       */
    double             xstd,    /* (I) normal point  */
    double            *yield)   /* (O) mapped yield  */
{
    static char      routine[]="Q3MQMapCR";
    int              status = FAILURE;

    /* distribution data */

    double fwdRate = mq->fwdRate;
    double sigTotal = mq->sigMQ * sqrt(mq->optExpy);

    long   nq; /* number of intervals      */
    double *q; /* pointer to q parameters  */
	double *b; /* pointer to b parameters  */
    double *k; /* pointer to strike ratios */
	double *x;
    
    long   idx;
    int    sign;

    double xval = xstd * sigTotal;

    /* check forward */
    if (fwdRate < Q3_MIN_FWD) 
    {
        DR_Error("%s: invalid fwdRate = %lf\n", routine, fwdRate);
        goto RETURN;
    }
	
    /* for zero volatility, return forward */
    if (mq->sigMQ < Q3_MIN_VOL) 
    {
        *yield = fwdRate;
        return (SUCCESS);
    }    

    /* branch into low/high strikes (left/right) */
    if (xval <= mq->xL[0])
	{
		nq = mq->nbQL;
		q = mq->qL;
		b = mq->bL;
		k = mq->kL;
		x = mq->xL;
		sign = -1;
	} 
	else 
	{
		nq = mq->nbQR;
		q = mq->qR;
		b = mq->bR;
		k = mq->kR;
		x = mq->xR;
		sign = 1;
	}
	
    /* find x location */
    idx = 0;
    while (idx < nq - 1 && sign * (xval - x[idx+1]) >= 0.) idx++;
	
    /* evaluate mapping function */

	*yield = k[idx] + b[idx] * qxq(q[idx], xval, x[idx]);

	*yield = mq->fwdRate*(mq->K*(mq->C*yield[0] - 1) + 1);
	
    status = SUCCESS;

 RETURN:
    
    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3MQMapCR */




/*f----------------------------------------------------------------------------
 * MultiQ: mapping derivative
 *
 * Evaluates mapping derivative at standard normal point x given MQDATA data.
 */
static int Q3MQMapDerivativeCR(
    MQDATA            *mq,      /* (I) MQ data       */
    double             xstd,    /* (I) normal point  */
    double            *dHdx)    /* (O) derivative    */
{
    static char      routine[]="Q3MQMapDerivativeCR";
    int              status = FAILURE;

    /* distribution data */

    double fwdRate = mq->fwdRate;
    double sigTotal = mq->sigMQ * sqrt(mq->optExpy);

    long   nq; /* number of intervals      */
    double *q; /* pointer to q parameters  */
	double *b; /* pointer to b parameters  */
    double *k; /* pointer to strike ratios */
	double *x;
    
    long   idx;
    int    sign;

    double xval = xstd * sigTotal;

    /* check forward */
    if (fwdRate < Q3_MIN_FWD) 
    {
        DR_Error("%s: invalid fwdRate = %lf\n", routine, fwdRate);
        goto RETURN;
    }

    /* branch into low/high strikes (left/right) */
    if (xval <= mq->xL[0])
	{
		nq = mq->nbQL;
		q = mq->qL;
		b = mq->bL;
		k = mq->kL;
		x = mq->xL;
		sign = -1;
	} 
	else 
	{
		nq = mq->nbQR;
		q = mq->qR;
		b = mq->bR;
		k = mq->kR;
		x = mq->xR;
		sign = 1;
	}
	
    /* find x location */
    idx = 0;
    while (idx < nq - 1 && sign * (xval - x[idx+1]) >= 0.) idx++;
	
    /* evaluate mapping function derivative */
    *dHdx = mq->fwdRate * mq->K * mq->C * b[idx] * sigTotal *
        exp(q[idx] * (xval-x[idx]));
	
    status = SUCCESS;

 RETURN:
    
    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3MQMapDerivativeCR */


/*f--------------------------------------------------------------------------
 * MultiQ: probability density function
 *
 * Computes PDF for mapped y-value.
 */
int Q3MQDensityCR(
    MQDATA *mq,
    double  yval,
    double *probDensity)
{
    static char routine[] = "Q3MQDensityCR";
    int         status    = FAILURE;

    /* x is standard normal ~ N(0,1) */
    double x;
    double dHdx;

    if (Q3MQMapInverseCR (mq, yval, &x) != SUCCESS)
        goto RETURN;

    if (Q3MQMapDerivativeCR (mq, x, &dHdx) != SUCCESS)
        goto RETURN;

    *probDensity = NormDens(x) / dHdx;
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: Failed", routine);

    return status;
}


/*f--------------------------------------------------------------------------
 * MultiQ: cumulative probability function
 *
 * Computes CPF for mapped y-value.
 */
int Q3MQCumCR(
    MQDATA *mq,
    double  yval,
    double *cumulativeProb)
{
    static char routine[] = "Q3MQCumCR";
    int         status    = FAILURE;

    /* x is standard normal ~ N(0,1) */
    double x;

    if (Q3MQMapInverseCR (mq, yval, &x) != SUCCESS)
        goto RETURN;

    *cumulativeProb = NormCum(x);
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: Failed", routine);

    return status;
}

/*f----------------------------------------------------------------------------
 * MultiQ: mapping inverse function
 *
 * Converts from strike value to standard normal point x given MQDATA.
 */
int Q3MQMapInverseCR(
    MQDATA            *mq,      /* (I) MQ data       */
    double             yval,    /* (I) Mapped point  */
    double            *xval)    /* (O) normal point  */
{
    static char      routine[]="Q3MQMapInverseCR";
    int              status = FAILURE;
 
    double fwdRate  = mq->fwdRate;
    double sigTotal = mq->sigMQ * sqrt(mq->optExpy);
    
    long   nq;          /* number of intervals                    */
    double strkMult;    /* strike/fwdRate                         */
    double *q;          /* pointer to q parameters                */
    double *b;          /* pointer to b parameters                */
    double *k;          /* pointer to strike ratios               */
    double *x;
    double *d;
    double CK;
    
    double temp, b2q;
    long   idx;
    int    sign;

    /* check vol */
    if(sigTotal < 0){
        DR_Error("%s: variance is negative", routine);
        goto RETURN;
    }
    
    /* write a special case for zero total vol! */
    if (sigTotal < Q3_MIN_VOL)
    {
        if (ARE_ALMOST_EQUAL(yval, fwdRate))
        {
            *xval = 0.0;
            status = SUCCESS;
        }
        else
        {
            DR_Error("%s: Cannot convert %f to normal point for zero "
                     "total volatility unless = forward %f\n",
                     routine, yval, fwdRate);
        }
        goto RETURN;
    }
    
    /* branch into low/high strikes (left/right); subtle      */
    /* choice needed in calibration: for ATM option, use left */
    /* intervals for put and high intervals for call          */
    strkMult = yval/fwdRate;
    
    CK = mq->C*mq->K;
    strkMult = strkMult - (1 - mq->K);
    /* now strkMult is the strike multiple for CK*H(X) */
    
    if (strkMult < CK) 
    {
        nq = mq->nbQL;
        q = mq->qL;
        b = mq->bL;
        d = mq->dL;
        k = mq->kL;
        x = mq->xL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        b = mq->bR;
        d = mq->dR;
        k = mq->kR;
        x = mq->xR;
        sign = 1;
    }
    
    /* find yval location */
    idx = 0;
    while (idx < nq - 1 && sign * (strkMult - CK * k[idx+1]) >= 0.) idx++;
    
    /* find normal point corresponding to yval */
    
    if(fabs(q[idx]) < TINY) /* q[idx]=0, normal dist int */     
    {         
        temp  = (strkMult/CK - k[idx])/b[idx];         
        *xval = x[idx] + temp;        
    } 
    else 
    {
        b2q   = b[idx]/q[idx];         
        temp  = 1 + (strkMult/CK - k[idx])/b2q;         
        *xval = (temp > 0 ? x[idx] + log(temp)/q[idx] : sign*Q3_CUM_NUM_STDEV);
    }
    
    *xval /= sigTotal;

    status = SUCCESS;
    
 RETURN:
    
    if (status == FAILURE) 
    {
        DR_Error("%s: Failed", routine);
    }
    return status;
} /* Q3MQMapInverseCR */
