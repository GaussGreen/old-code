/******************************************************************************
 * Module:    CRXQ
 * File:      crxqfapricer.c
 * Function:    
 * Author:    Credit QRD, Charles Morcom, based on univariate work by JC Porras
 *****************************************************************************/
#include <math.h>
#include <stdio.h>

#include "crxq.h"
#include "q3.h"
#include "crxerror.h"

/* 
	Analytic approximation to the price of protection given a CDS spread and an
	equivalent interest rate, which makes the approximation correct for the 
	case where the input spread is the forward par spread
*/
double Protection(double parSpread, 
                  const CRXQFADATA* fa)
{
	double x;
    double R = fa->recovery;
    // approximate clean average/zero spread is VNFM approx fn of par spread
	double lambda = fa->alphaAnn*pow(parSpread, fa->powerAnn);
    // survival probability to maturity
    double zsp = pow(1 + lambda/fa->freqAnn, -fa->matAnn*fa->freqAnn);
    // zero IR
    double r = fa->IRate;
    // riskless discount to maturity
    double zdf = pow(1 + r/fa->freqSwap, -fa->matAnn*fa->freqSwap);

    /* This approximation assumes that the clean spread and
       the zero interest rate are constant from T to M, so
       that the protection leg may be integrated exactly.
       It also assumes that the rates are continuously
       compounded, which is not true, but the error is only
       in the lambda/(lambda+r) factor, and will be very small. */
	x = (1 - R) * (1 - zsp * zdf) * lambda/(lambda+r);  
	return  x;

} /* Protection */

/******************************************************************************
 * CRXQFAGridPricer
 * Computes the change of measure to FA, and then does 1D numerical
 * integration to find the option price. Bivariate options are priced
 * by the inner integral being explicitly evaluated inside the payoff
 * function. Assumes that the multi-q measure is already calibrated. 
 *****************************************************************************/
int CRXQFAGridPricer(
    CRXQPAYOFF  *optPayoff, /* (I) option payoff structure */
    CRXQFPAYOFF *payFunc,   /* (I) option payoff function  */
    CRXQFADATA   *fa,        /* (I) measure data            */
    double     *premium    /* (O) fwd premium/rate        */
    )
{
    static char routine[] = "CRXQFAGridPricer";
    int         status    = FAILURE;

    /* variables to hold numerical FA measure grid values and densities */
    double grid[CRXQ_1D_INT_PTS], yield[CRXQ_1D_INT_PTS], dens[CRXQ_1D_INT_PTS];
    double densNorm;

    /* calculate FA rate density */
    if (CRXQFADens(            
        fa,
        -CRXQ_1D_NUM_STDEV,   /* start */
        CRXQ_1D_NUM_STDEV,    /* end   */
        CRXQ_1D_INT_PTS,
        &(grid[0]),
        &(yield[0]),
        &(dens[0]),
        &densNorm) == FAILURE) goto RETURN;

    /* integrate density * payoff */
    if (CRXQSimpsonPricer1D(
        optPayoff,
        payFunc,
        &(grid[0]),
        &(yield[0]),
        &(dens[0]),
        densNorm,
        premium) == FAILURE) goto RETURN;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQFAGridPricer */


/* JCP: change to credit measure change (Eannuity/Annuity(L) * Z(L) */
/*f----------------------------------------------------------------------------
 * CRXQFADens
 *
 * Generates convexity and delay adjusted density for the adjusted-forward
 * (a.k.a. payment) measure.
 * The type of annuity adjustment depends on whether the rate is
 * an interest rate or a credit spread.
 *
 */
int CRXQFADens(
    const CRXQFADATA *fa,       /* (I) measure data              */
    double  start,      /* (I) start in std normal space */
    double  end,        /* (I) end in std normal space   */
    long    numGridPts, /* (I) num pts to generate       */
    double *grid,       /* (O) grid in normal space      */
    double *yields,     /* (O) grid in yield space       */
    double *faDens,     /* (O) probability density       */
    double *normC       /* (O) density normalization     */
    )
{
    static char routine[] = "CRXQFADens";
    int         status    = FAILURE;

	double delayFactor; // discount factor for payment delay
    double convexityFactor; // annuity adjustment for annuity convexity

    double p, z, x;
    const CRXQDATA* mq = fa->mq;

    /* low rate cutoff for Ann and Del functional forms */
    double yCutoff = CRXQ_FA_BDRY * fa->mq->fwdRate;

    /* Z-Grid step size */
    double step = (end - start)/(numGridPts-1);
    /* Current xval - incremented at head of loop */
    double xval = start-step;

    double sigMQ = mq->sigMQ;
    double muMQ = mq->muMQ;
    double y;  // place-holder for current rate point = yields[i]
    double yCut;
    long   i;

    /* Convexity and Delay change of measure params from VNFM power approx */
    /* Annuity Convexity */
    double freqAnn    = fa->freqAnn;
    double matAnn     = fa->matAnn;
    double powerAnn   = fa->powerAnn;
    double alphaAnn   = fa->alphaAnn;
    /* Delay */
    double freqDel    = fa->freqDel;
    double matDel     = fa->matDel;
    double powerDel   = fa->powerDel;
    double alphaDel   = fa->alphaDel;


    /* calculate mapping and FA density function */
    for (i = 0; i < numGridPts; i++)
    {
        xval += step;
        grid[i] = muMQ + sigMQ * xval;
        if (CRXQMap(
            mq,
            grid[i],
            &(yields[i])) == FAILURE) 
        {
            DR_Error("%s failed: invalid q-mapping for x-grid=%lf.\n", routine, grid[i]);
            goto RETURN;
        }
        
        y = yields[i];
        yCut = MAX (y, yCutoff);        // floor for very small y
        p    = MIN (y / yCutoff, 1.);

        /* Compute the convexity adjustment - annuity depends on
           whether y is credit or rate */
        if (CRXQ_CREDIT_SPREAD==fa->rateType)
        {
            /* Rate is a credit spread - use CDS annuity 
               The covexity adjustment is 1/annuity, and 
               the annuity is Protection/spread */
            convexityFactor = yCut / Protection(yCut, fa);

        } else {

            /* Rate is an interest rate: use swap annuity, which is
               just (1 - DF)/Y                                      */
            if (y > TINY) 
            {
                z    = pow(yCut, powerAnn);
                x    = 1. + alphaAnn * z / freqAnn;
                convexityFactor = yCut / (1.-pow(x, -freqAnn * matAnn));

                /* interpolate between 0 and yCutoff */
                convexityFactor = p * convexityFactor + (1 - p) / matAnn;
            } 
            else 
            {
                convexityFactor = 1. / matAnn;
            }
        }

        /* Compute the delay adjustment factor: 
           same computation for cr and rate y */
        if (y > TINY)
        {
            z = pow(yCut, powerDel);
            x = 1.0 + alphaDel * z / freqDel;
            delayFactor = p*pow(x, -freqDel*matDel) + (1-p);

        } else {
            delayFactor = 1.0;
        }
        
        faDens[i] = NormDens(xval) * delayFactor * convexityFactor;
            
    }
    
    /* normalize distribution (from q3.h) */
    *normC = Q3SimpsIntegral(
        step,
        numGridPts,
        faDens);
    if (*normC < TINY) 
    {
        DR_Error("%s: Non-positive total FA mass %lf.\n", routine, *normC);
        goto RETURN;
    }
    
    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQFADens */

/******************************************************************************
 * CRXQFACreditInit
 * If the rate is a credit spread rather than an interest-rate then, in order
 * to calculate the annuity for the measure change, you must know the IR
 * zero rate for the spread period, and the recovery rate.
 *****************************************************************************/
int CRXQFACreditInit(
    double     recoveryRate,
    double     accrualFactor,
    double     interestRate,
    double     interestRateFreq,
    double     riskyAnnuity,
    CRXQFADATA* fa
    ) 
{

    static char      routine[]="CRXQFACreditInit";
    int              status = FAILURE; 

    fa->rateType = CRXQ_CREDIT_SPREAD;
    fa->freqSwap = interestRateFreq;
    fa->recovery = recoveryRate;
    fa->IRate = interestRate;
    fa->accrFactor = accrualFactor;
    fa->ryAnnuity = riskyAnnuity;

    status = SUCCESS;
    return status;
}


/*f----------------------------------------------------------------------------
 * CRXQFASmileInit
 * 
 * Determines the functional form for the change of measure between the
 * FA and AA measures (convexity and delay), based on VNFM
 *
 */
int CRXQFASmileInit(
    double  expiry,        /* (I) observation date in yrs      */
    double  sigATM,        /* (I) ATM vol                      */
    double  start,         /* (I) swap/zero rate start in yrs  */
	long    freq,    
    double  swapMat,       /* (I) swap rate tenor in yrs       */
    double  swapRate,      /* (I) par swap rate                */
    double  fwdAnnuity,    /* (I) forward (swap) annuity       */
    double  zeroRateSwap,  /* (I) zero rate for same interval  */
    double  payDelay,      /* (I) payment delay                */   
    double  zeroRatePay,   /* (I) zero rate for delay interval */    
    long    numVnfmParams, /* (I) number of vnfm params        */
    double *vnfmParams,    /* (I) vnfm model parameters        */
    long    convexAdjSetl,  /* (I) settlement type              */
    long    delayAdjSetl,
    CRXQFADATA *fa             /* (I/O) FA data                    */
    )
{   
    static char      routine[]="CRXQFASmileInit";
    int              status = FAILURE; 

    double zeroSwapVolRatio, swapVol, zeroVol, zeroSwapCorr;

    /* This sets up the CRXQFADATA structure to be a rate */
    fa->rateType = CRXQ_INTEREST_RATE;
    fa->recovery = 0.0;
    fa->IRate = 0.0;
    fa->accrFactor = 0.0;

    /* rate frequencies and maturities */
    fa->freqAnn = fa->freqDel = fa->freqSwap = freq;
    fa->matAnn = swapMat;
    fa->matDel = payDelay;

    /* if cash settled, alpha = power = 1; otherwise, use VNFM */
    if (convexAdjSetl == CRXQ_CASH_SETL) 
    {
        /* no annuity adjustment needed */
        fa->alphaAnn = fa->powerAnn = 1.;
    } 
	else
	{
        /* From q3.h */
		if (Q3VNFMZero2Swap(
            expiry,
			0., // start from now
            start,
            freq,
            swapMat,
            swapRate,
            fwdAnnuity,
            swapMat,
            zeroRateSwap,
            numVnfmParams,
            vnfmParams,
            &swapVol,
            &zeroVol,
            &zeroSwapCorr) == FAILURE) 
        {
            DR_Error("%s failed: unable to calculate annuity VNFM power form for rate/spread.\n", routine);
            goto RETURN;
        }

        zeroSwapVolRatio = zeroVol / swapVol;
        fa->powerAnn = zeroSwapCorr * zeroSwapVolRatio;
        fa->alphaAnn = zeroRateSwap / pow(swapRate, fa->powerAnn) *
            exp(0.5 * sigATM * sigATM * expiry * fa->powerAnn *
                (1. - fa->powerAnn));      
    }

    /* VNFM calculation: functional form for delay adjustment */
    if (delayAdjSetl == CRXQ_CASH_SETL || payDelay < TINY) 
    {
        /* no delay adjustment needed */
        fa->powerDel = fa->alphaDel = 1.;
    } 
    else 
    {
	  if (Q3VNFMZero2Swap(
            expiry,
			0.,
            start,
            freq,
            swapMat,
            swapRate,
            fwdAnnuity,
            payDelay,
            zeroRatePay,
            numVnfmParams,
            vnfmParams,
            &swapVol,
            &zeroVol,
            &zeroSwapCorr) == FAILURE) 
      {
          DR_Error("%s failed: unable to calculate delay VNFM power form for rate/spread.\n", routine);
          goto RETURN;
      }


        zeroSwapVolRatio = zeroVol / swapVol;
        fa->powerDel = zeroSwapCorr * zeroSwapVolRatio;
        fa->alphaDel = zeroRatePay / pow(swapRate, fa->powerDel) *
            exp(0.5 * sigATM * sigATM * expiry * fa->powerDel *
                (1. - fa->powerDel));      
    }

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQFASmileInit */




