/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File: stochvol.c     
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/

#include <math.h>
#include <ctype.h>                              /* toupper */
#include <stdio.h>

#include "q3.h"


/*f----------------------------------------------------------------------------
 * Q3SimpsIntegral    
 * 
 * Numerical integration : 1D Simpson (from data points).
 * <br><br>
 * Integrate function given by its values
 * using extended Simpson's rule
 * For less than 8 points use central Riemann sums
 * For a single point use value * step
 */

double  Q3SimpsIntegral
(
    double            step,             /* (I) */
    long              numPoints,        /* (I) */
    double           *fValues)          /* (I) */ 
{
    double c1 = 17./48.;
    double c2 = 59./48.;
    double c3 = 43./48.;
    double c4 = 49./48.;

    long i;
    long numEndPts; /*  for each end of the interval */
    double  intValue = 0.;
    double *fValue = NULL;

    if(numPoints >= 8) 
    {
        numEndPts = 4;
        intValue =
            (fValues[0] + fValues[numPoints-1]) * c1 +
            (fValues[1] + fValues[numPoints-2]) * c2 +
            (fValues[2] + fValues[numPoints-3]) * c3 +
            (fValues[3] + fValues[numPoints-4]) * c4;
    } 
    else if (numPoints >= 2) 
    {
        numEndPts = 1;
        intValue = 
            (fValues[0] + fValues[numPoints-1]) * 0.5;
    } 
    else 
    {
        numEndPts = 0;
    }
        
    for (i=2*numEndPts,fValue=fValues+numEndPts; i<numPoints;
         i++,fValue++) 
    {
        intValue += *fValue;
    }

    intValue *= step;

    return (intValue);
} /* Q3SimpsIntegral */

  
/*f----------------------------------------------------------------------------
 * Stoch Vol Black pricer     
 *
 * SV pricer for call/put. Assumes internal vol from SVDATA.
 * NOTE: price, vega on demand.
 */
int Q3SVPricer(
    SVDATA            *sv,           /* (I) SV data  */
    long              type,          /* (I) option type */
    double            strike,        /* (I) option strike */
    double            *price)        /* (O) option price & vega    */
{
    static char      routine[]="Q3SVPricer";
    int              status = FAILURE;

    double    step, normPt, vol, volStep;   /* vol integration vars    */
    double    expFact;                      /* geom brown fact         */
    double   *gridP = NULL, *gridV = NULL;  /* values of price & vega  */
    long      numGridPts = Q3_SV_INT_PTS;
    long      i;
    long      vOn;                          /* vega calc on or off     */
    
    double    expiry, fwdRate, volVolSV, sigSV, sigATM, q;   /* locals */


    /* set local variables from sm structure */
    expiry   = sv->expiry;
    fwdRate  = sv->fwdRate;
    volVolSV = sv->volVolSV; 
    sigSV    = sv->sigSV;     
    sigATM   = sv->sigATM;
    q        = 1-sv->q;

    /* input defense */
    if (Q3_COP_TYPE(type) != Q3_CALL && 
        Q3_COP_TYPE(type) != Q3_PUT  &&
        Q3_PAY_TYPE(type) != Q3_VNL  &&
        Q3_PAY_TYPE(type) != Q3_VEGA) goto RETURN;

    /* set integration for vega */
    vOn = (Q3_PAY_TYPE(type) == Q3_VEGA);

    /* check inputs: forward rate > 0, vol and vol of vol < max values */
    if (fwdRate                   < Q3_MIN_FWD     ||
        volVolSV * sqrt(expiry/3) > Q3_MAX_VOL_VOL ||
        sigATM   * sqrt(expiry)   > Q3_MAX_VOL       ) 
		
	{
		goto RETURN;
	}

    /* for sigATM = 0 compute intrinsic, no integration */
    if (fabs(sigATM * sqrt(expiry)) < Q3_MIN_VOL) 
    {
        if (BSQPricer(
            fwdRate,
            strike,
            expiry,
            sigATM,
            q,
            type,
            price) == FAILURE) goto RETURN;
        status = SUCCESS;
        goto RETURN;
    }

    /* allocate grid array */
    if ((gridP = DR_Array(DOUBLE,0,numGridPts-1)) == NULL) goto RETURN;
    if ((gridV = DR_Array(DOUBLE,0,numGridPts-1)) == NULL) goto RETURN;

    /* set grid  */ 
    volVolSV *= sqrt(expiry/3.);
    step      = 2.*Q3_SV_NUM_STDEV / numGridPts;
    normPt    = -Q3_SV_NUM_STDEV + step*0.5; 
    expFact   = exp((normPt-0.5*volVolSV)*volVolSV);
    vol       = sigSV * expFact;
    volStep   = exp(step*volVolSV);

    for (i=0; i<numGridPts; i++) 
    {
        double p[2],prob;

        if (BSQPricer(fwdRate,strike,expiry,vol,q,type,p) == FAILURE)
            goto RETURN;
        prob = NormDens(normPt);
        gridP[i] = prob * p[0];

        /* use chain rule to relate BSQPricer vega to SVPricer vega */
        if (vOn) gridV[i] = prob * p[1] * expFact;

        normPt  += step;
        vol     *= volStep;
        expFact *= volStep;
    }
    
    /* integrate */
    price[0] = Q3SimpsIntegral(step, numGridPts, gridP);
    if (vOn) price[1] = Q3SimpsIntegral(step, numGridPts, gridV);

    status = SUCCESS;

 RETURN:
    
    Free_DR_Array(gridP,DOUBLE,0,numGridPts-1);
    Free_DR_Array(gridV,DOUBLE,0,numGridPts-1);

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3SVPricer */


/*f----------------------------------------------------------------------------
 * Stoch Vol internal calibration    
 *
 * SV Newton-Raphson calibration for internal vol in SVDATA
 */
int Q3SVCalib(
    SVDATA             *sv)           /* (I) SV data  */
{
    static char      routine[]="Q3SVCalib";
    int              status = FAILURE;

    long      j;

    double    expiry, fwdRate, sigATM, q;  
    double    price[2], atmPrice;

    /* set local variables from sm structure */
    expiry  = sv->expiry;
    fwdRate = sv->fwdRate;
    sigATM  = sv->sigATM;
    q       = 1-sv->q;

    if (fwdRate < Q3_MIN_FWD) goto RETURN;

    /* calibration defense: for sigATM=0 no calibration */
    if (fabs(sigATM * sqrt(expiry)) < Q3_MIN_VOL) 
    {
        sv->sigSV = 0.;
        status = SUCCESS;
        goto RETURN;
    }

    /* price target option */
    if (BSQPricer(
        fwdRate,
        fwdRate,
        expiry,
        sigATM,
        1,
        Q3_CALL,
        price) == FAILURE) goto RETURN;
    atmPrice = price[0];

    /* solve for SV vol numerically */
    sv->sigSV = sigATM;

    for(j=0; j<Q3_NUM_SOLVER_ITER; j++) 
    {
        /* price and vega */
        if (Q3SVPricer(sv,Q3_CALL+Q3_VEGA,fwdRate,price) == FAILURE)
            goto RETURN;

        if(fabs(price[1]/(sv->sigSV)) < TINY) 
        {
            Q3ErrMsg("%s: Cannot calibrate average vol.\n", 
                     routine); goto RETURN;
        }
        sv->sigSV -= (price[0]-atmPrice) / price[1];
    }

    if(sv->sigSV < 0.) 
    {
        Q3ErrMsg("%s: Adjusted vol (%f) out of bounds.\n",
                 routine, sv->sigSV);
        goto RETURN;
    }

    status = SUCCESS;

 RETURN:
    
    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3SVCalib */


/*f----------------------------------------------------------------------------
 * Q3SVSmileInit
 *
 * Stoch Vol: initialize SV smile.
 */
int Q3SVSmileInit(
    double             fwdRate,         /* (I) fwd rate */
    double             sigATM,          /* (I) atm vol */
    double             expiry,          /* (I) in years */
    double            *smile,           /* (I) smile: SV smile */
    SVDATA            *sv)              /* (0) SV measure  */
{
    static char      routine[]="Q3SVSmileInit";
    int              status = FAILURE;

    double totVol, volVol, totVolVol;

    /* check inputs: expiry */
    if (expiry < 0.) 
    {
        Q3ErrMsg("%s: Expiry (%f) must be >= 0\n",
                 routine, expiry);
        goto RETURN;
    }

    /* fwd rate limit */
    if (fwdRate < Q3_MIN_FWD) 
    {
        Q3ErrMsg("%s: Fwd rate (%f) must be > %f.\n", 
                 routine, fwdRate, Q3_MIN_FWD);
        goto RETURN;
    }

    /* vol lower limit */
    if (sigATM < 0.) 
    {
        Q3ErrMsg("%s: Negative vol = %f\n", routine, sigATM);
        goto RETURN;
    }

    /* vol upper limit */
    totVol = sigATM * sqrt(expiry);
    if (totVol > Q3_MAX_VOL) 
    {
        Q3ErrMsg("%s: Atm vol out of bounds: vol = %f, T = %f\n", 
                 routine, sigATM, expiry);
        goto RETURN;
    }

    /* vol of vol lower limit */
    volVol = smile[1] * pow(sigATM, smile[2]) * pow(fwdRate, smile[3]);
    if (volVol < 0.) 
    {
        Q3ErrMsg("%s: Negative vol of vol = %f\n", routine, volVol);
        goto RETURN;
    }
   
    /* vol of vol upper limit */
    totVolVol = volVol * sqrt(expiry/3);
    if (totVolVol > Q3_MAX_VOL_VOL) 
    {
        Q3ErrMsg("%s Vol of vol out of bounds: vov = %f, T = %f\n", 
                 routine, volVol, expiry);
        goto RETURN;
    }

    /* skew limits */
    if (fabs (1. - smile[0]) * totVol > Q3_MAX_SKEW + TINY) 
    {
        Q3ErrMsg("%s: Stoch Vol skew (%6.4f) must be between %6.4f and "
                 "%6.4f.\n", routine, smile[0], 1. - Q3_MAX_SKEW / totVol, 
                 1. + Q3_MAX_SKEW / totVol);
        goto RETURN;
    }
    
    /* initialize Stoch Vol part of MQDATA structure */
    sv->fwdRate   = fwdRate;
    sv->expiry    = expiry;
    sv->sigATM    = sigATM;
    sv->q         = smile[0];
    sv->volVolMkt = smile[1];
    sv->bbV       = smile[2];
    sv->bbR       = smile[3];
    sv->volVolSV  = volVol;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: Failed\n", routine);
    }
  
    return status;
} /* Q3SVSmileInit */

