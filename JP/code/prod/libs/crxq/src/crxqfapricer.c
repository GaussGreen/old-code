/******************************************************************************
 * Module:      CRXQ
 * Submodule:
 * File:        crxqfapricer.c  
 * Function:    
 * Author:      Credit QRD
 * Revision:    $Header: $
 *****************************************************************************/

#include <math.h>
#include <stdio.h>

#include <crxmultiq/include/crmultiq.h>
#include <crxflow/include/crxerror.h>
#include <crxflow/include/optbsq.h>
#include "crxq.h"


/*-----------------------------------------------------------------------------
 * ProtectionCalc
 *
 * Analytic approximation to the price of protection given a CDS spread and an
 * equivalent interest rate, which makes the approximation correct for the
 * case where the input spread is the forward par spread
 *
 */ 
double ProtectionCalc(
       double parSpread,         /* (I) cds par spread */
       const CRXQFADATA* fa      /* (I) FA data        */
        )
{
    double x;                    /* (O) lambda/(lambda+r) *(1 - zsp * zdf)     */
    double r = fa->IRate;        /* zero IR                                    */
    double lambda;               /* approximate clean average/zero spread 
                                    is VNFM approx fn of par spread            */

    /* This approximation assumes that the clean spread and
       the zero interest rate are constant from T to M, so
       that the protection leg may be integrated exactly.
       It also assumes that the rates are continuously
       compounded, which is not true, but the error is only
       in the lambda/(lambda+r) factor, and will be very small.   */
    /* survival probability to maturity:  
       zsp = pow(1 + lambda/fa->freqAnn, -fa->matAnn*fa->freqAnn) */
    /* riskless discount to maturity:
       zdf = pow(1 + r/fa->freqSwap, -fa->matAnn*fa->freqSwap)    */

    lambda = fa->alphaAnn*pow(parSpread, fa->powerAnn);
    x = lambda/(lambda+r) *
            (1 - pow(1+fa->alphaAnn*pow(parSpread, fa->powerAnn)/fa->freqAnn,
			-fa->matAnn*fa->freqAnn )
            * pow(1+r/fa->freqSwap, -fa->matAnn*fa->freqSwap));  

    return  x;

} /* ProtectionCalc */


/*-----------------------------------------------------------------------------
 * CRXQFAGridPricer
 *
 * Computes the change of measure to FA, and then does 1D numerical
 * integration to find the option price. Bivariate options are priced
 * by the inner integral being explicitly evaluated inside the payoff
 * function. Assumes that the multi-q measure is already calibrated. 
 *
 */ 
int CRXQFAGridPricer(
    CRXQPAYOFF     *optPayoff, /* (I) option payoff structure */
    CRXQFPAYOFF    *payFunc,   /* (I) option payoff function  */
    CRXQFADATA     *fa,        /* (I) measure data            */
    double         *premium    /* (O) fwd premium/rate        */
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
    double  start,              /* (I) start in std normal space */
    double  end,                /* (I) end in std normal space   */
    long    numGridPts,         /* (I) num pts to generate       */
    double *grid,               /* (O) grid in normal space      */
    double *yields,             /* (O) grid in yield space       */
    double *faDens,             /* (O) probability density       */
    double *normC               /* (O) density normalization     */
    )
{
    static char routine[] = "CRXQFADens";
    int         status    = FAILURE;

    double delayFactor;         /* discount factor for payment delay */
    double convexityFactor;     /* annuity adjustment for annuity convexity */

    double p, z, x;
    const CRXQDATA* mq = fa->mq;

    /* low rate cutoff for Ann and Del functional forms */
    double yCutoff = CRXQ_FA_BDRY * fa->mq->fwdRate;

    double step = (end - start)/(numGridPts-1);
    double xval = start-step;

    double sigMQ = mq->sigMQ;
    double muMQ = mq->muMQ;
    double y;  /* place-holder for current rate point = yields[i] */
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
            DR_Error("%s failed: invalid q-mapping for x-grid=%lf.\n",
                       routine, grid[i]);
            goto RETURN;
        }
        
        y = yields[i];
        yCut = MAX (y, yCutoff);        /* floor for very small y */
        p    = MIN (y / yCutoff, 1.);

        /* Compute the convexity adjustment - annuity depends on
           whether y is credit or rate */
        if (CRXQ_CREDIT_SPREAD==fa->rateType)
        {
            /* Rate is a credit spread - use CDS annuity 
            The covexity adjustment is 1/annuity, and the 
            annuity is Protection/spread */
            convexityFactor = (yCut / (1.-fa->recovery)) 
					/ ProtectionCalc(yCut, fa);
        }
        else
        {
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
    
    /* normalize distribution (from crxq.h) */
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


/*f----------------------------------------------------------------------------
 * CRXQFACreditInit
 * 
 * If the rate is a credit spread rather than an interest-rate then, in order
 * to calculate the annuity for the measure change, you must know the IR
 * zero rate for the spread period, and the recovery rate.
 *
 */
int CRXQFACreditInit(
    double     recoveryRate,
    double     accrualFactor,
    double     interestRate,
    double     interestRateFreq,
    double     riskyAnnuity,
    CRXQFADATA* fa
    ) 
{

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
    long    rateType,       /* (I) credit spread or IR          */
    double  expiry,         /* (I) observation date in yrs      */
    double  sigATM,         /* (I) ATM vol                      */
    double  start,          /* (I) swap/zero rate start in yrs  */
    long    freq,           /* (I) index frequency              */
    double  swapMat,        /* (I) swap rate tenor in yrs       */
    double  cdsParSpread,   /* (I) cds par spread               */
    double  swapRate,       /* (I) par swap rate                */
    double  fwdAnnuity,     /* (I) forward (swap) annuity       */
    double  zeroRateSwap,   /* (I) zero rate for same interval  */
    double  payDelay,       /* (I) payment delay                */   
    double  zeroRatePay,    /* (I) zero rate for delay interval */    
    long    numVnfmParams,  /* (I) number of vnfm params        */
    double *vnfmParams,     /* (I) vnfm model parameters        */
    long    convexAdjSetl,  /* (I) settlement type              */
    long    delayAdjSetl,
    CRXQFADATA *fa          /* (I/O) FA data                    */
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
    fa->freqAnn = fa->freqDel = freq;        /* index frequency */
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
        if (CRXQVNFMZero2Swap(
            expiry,
            0., /* start from now  */
            start,
            freq,          /* frequency of index  */
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
            DR_Error("%s failed: unable to calculate annuity VNFM power form \
                       for rate/spread.\n", routine);
            goto RETURN;
        }

        zeroSwapVolRatio = zeroVol / swapVol;
        fa->powerAnn = zeroSwapCorr * zeroSwapVolRatio;

        if (rateType == CRXQ_INTEREST_RATE)
        {
            /* compute alphaAnn based on swap rate */
            fa->alphaAnn = zeroRateSwap / pow(swapRate, fa->powerAnn) *    
                exp(0.5 * sigATM * sigATM * expiry * fa->powerAnn *
                    (1. - fa->powerAnn));      
        }
        else
        {
            /* compute alphaAnn based on cdsParSpread */
            fa->alphaAnn = zeroRateSwap / pow(cdsParSpread, fa->powerAnn) *    
                exp(0.5 * sigATM * sigATM * expiry * fa->powerAnn *
                    (1. - fa->powerAnn));      
        }

    }

    /* VNFM calculation: functional form for delay adjustment */
    if (delayAdjSetl == CRXQ_CASH_SETL || payDelay < TINY) 
    {
        /* no delay adjustment needed */
        fa->powerDel = fa->alphaDel = 1.;
    } 
    else 
    {
        if (CRXQVNFMZero2Swap(
            expiry,
            0.,
            start,
            freq,              /* frequency of index */
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
          DR_Error("%s failed: unable to calculate delay VNFM power form \
                    for rate/spread.\n", routine);
          goto RETURN;
        }

        zeroSwapVolRatio = zeroVol / swapVol;
        fa->powerDel = zeroSwapCorr * zeroSwapVolRatio;
        if (rateType == CRXQ_INTEREST_RATE)
        {
           /* compute alphaDel based on swap rate */
	        fa->alphaDel = zeroRatePay / pow(swapRate, fa->powerDel) *
            exp(0.5 * sigATM * sigATM * expiry * fa->powerDel *
                (1. - fa->powerDel));      
        }
        else
        {
            /* compute alphaDel based on cdsParSpread */
            fa->alphaDel = zeroRatePay / pow(cdsParSpread, fa->powerDel) *
                exp(0.5 * sigATM * sigATM * expiry * fa->powerDel *
                    (1. - fa->powerDel));
		}
    }

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQFASmileInit */


#define  BOND_PRICE(c,f,m,y) \
         (fabs(y)<TINY) ? (1 + (m) * (c)) : ((c) / (y) + \
         (1- (c) / (y)) * pow(1 + (y)/(f), -(m) * (f))) 
/*-----------------------------------------------------------------------------
 * CRXQVNFMZero2Swap
 * 
 * Given VNFM parameters (weights, mean reversions, correlations), evaluates
 * annualized zero and swap vols and the correlation between coupon and zero 
 * rates.
 *
 */
int CRXQVNFMZero2Swap(
    double  expiry,        /* (I) observation date in yrs         */
    double  volStart,      /* (I) start of vol observation in yrs */
    double  start,         /* (I) swap/zero rate start in yrs     */
    long    freq,          /* (I) swap/zero rate frequency        */
    double  swapMat,       /* (I) swap rate tenor in yrs          */
    double  swapRate,      /* (I) par swap rate                   */
    double  fwdAnnuity,    /* (I) forward (swap) annuity          */
    double  zeroMat,       /* (I) zero rate tenor in yrs          */   
    double  zeroRate,      /* (I) zero coupon rate                */    
    long    numVnfmParams, /* (I) number of vnfm params           */
    double *vnfmParams,    /* (I) vnfm model parameters           */
    double *swapVol,       /* (O) swap vol                        */
    double *zeroVol,       /* (O) zero vol (given freq)           */
    double *zeroSwapCorr   /* (O) zero / swap correlation         */   
    )
{   
    static char      routine[]="CRXQVNFMZero2Swap";
    int              status = FAILURE; 

    double beta[3], wght[3], rho[3][3];
    double swapBCoeff[3], zeroBCoeff[3];
    double sVol, zVol, coVar;
    long   numFact;
    int    i,j;

    /* Adjust swapMat */
    if(freq == 12 || freq == 4 || freq == 2 || freq == 1)
    {
        swapMat= floor(freq*swapMat + 0.5) / (double) freq;
    }
    else if (freq == 365 || freq == 52 || freq == 26)
    {
        swapMat = 1.0/(double) freq;
    }
    else
    {
        DR_Error("%s: Invalid frequence: %ld \n", routine, freq);
        goto RETURN;
    }
    
    /* check swpMat, zeroMat */
    if (swapMat < TINY || zeroMat < TINY) goto RETURN;

    /* read VNFM parameters */
    switch(numVnfmParams) 
    {
    case 1:         /* 1 fact, alpha not given */
        numFact   = 1;
        beta[0]   = MAX(vnfmParams[0], CRXQ_BETA_SHIFT);
        wght[0]   = 1.;
        rho[0][0] = 1.;
        break;
    case 2:         /* 1 fact, alpha given */
        numFact   = 1;
        beta[0]   = MAX(vnfmParams[0], CRXQ_BETA_SHIFT);
        wght[0]   = vnfmParams[1];
        rho[0][0] = 1.;
        break;
    case 5:         /* 2 fact */
        numFact   = 2;
        beta[0]   = MAX(vnfmParams[0], CRXQ_BETA_SHIFT);
        beta[1]   = MAX(vnfmParams[1], CRXQ_BETA_SHIFT);
        wght[0]   = vnfmParams[2];
        wght[1]   = vnfmParams[3];
        rho[0][0] = rho[1][1] = 1.;
        rho[1][0] = rho[0][1] = vnfmParams[4];
        break;
    case 9:         /* 3 fact */
        numFact   = 3;
        beta[0]   = MAX(vnfmParams[0], CRXQ_BETA_SHIFT);
        beta[1]   = MAX(vnfmParams[1], CRXQ_BETA_SHIFT);
        beta[2]   = MAX(vnfmParams[2], CRXQ_BETA_SHIFT);
        wght[0]   = vnfmParams[3];
        wght[1]   = vnfmParams[4];
        wght[2]   = vnfmParams[5];
        rho[0][0] = rho[1][1] = rho[2][2] = 1.;
        rho[1][0] = rho[0][1] = vnfmParams[6];
        rho[2][0] = rho[0][2] = vnfmParams[7];
        rho[2][1] = rho[1][2] = vnfmParams[8];
        break;
    default:
        DR_Error("%s: Number of parameters (%d) must be 1, 2, 5, or 9.\n",
                 routine, numVnfmParams);
        goto RETURN;
    }

    /* setup B coefficients */
    for (i=0; i < numFact; i++) 
    {
        double swapPlusBeta, bondPrice1, bondPrice2;
    
        /* swap rate + mean reversion, on the swap rate compounding basis 
         * extra hassle to ensure cvx adj = delay adj when appropriate */
        swapPlusBeta = freq * ((1 + swapRate / freq) * 
                               exp (beta[i] / freq) - 1.);

        /* intermediate step: bond calculations */
        bondPrice1 = BOND_PRICE(swapRate, freq, swapMat, swapRate);
        bondPrice2 = BOND_PRICE(swapRate, freq, swapMat, swapPlusBeta);

        /* B coefficient for swap rate */
        swapBCoeff[i] = wght[i] * (bondPrice1 - bondPrice2) /
            (fwdAnnuity * beta[i]); 
            
        /* B coefficient for zero rate */
        zeroBCoeff[i] = wght[i] * (1 + zeroRate / freq) *
            (1 - exp( -beta[i] * zeroMat)) / (beta[i] * zeroMat);
    }

    /* total covariance */
    sVol = zVol = coVar = 0.;

    for (i=0; i< numFact; i++) 
    {
        for (j=0; j < numFact; j++) 
        {
            double covarGauss, betaSum, betaDur;
        
            betaSum = beta[i] + beta[j];
            betaDur = (betaSum*(expiry-volStart)>TINY)?
                (1-exp(-betaSum*(expiry-volStart))) / 
                (betaSum*(expiry-volStart)): 1.;

            covarGauss = rho[i][j] * exp(-betaSum*(start-expiry)) * betaDur;

            sVol  += swapBCoeff[i] * swapBCoeff[j] * covarGauss;
            zVol  += zeroBCoeff[i] * zeroBCoeff[j] * covarGauss;
            coVar += swapBCoeff[i] * zeroBCoeff[j] * covarGauss;
        }
    }

    *swapVol = sqrt(sVol);
    *zeroVol = sqrt(zVol);

    *zeroSwapCorr = coVar / sqrt(sVol*zVol);

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQVNFMZero2Swap */





