/******************************************************************************
 * Module:    Q3
 * Submodule:
 * File:      fapricer.c
 * Function:    
 * Author:    Interest Rates DR
 * Revision:  $Header$   
 *****************************************************************************/
#include <math.h>
#include <stdio.h>

#include "crmultiq.h"


/*f----------------------------------------------------------------------------
 * Q3FAGridPricer
 *
 * Change of measure: computes option prices after changing the measure
 * from AA to FA. Carries out both the CMS and delay adjustments. Assumes
 * Multiq measure has already been calibrated.
 */
int Q3FAGridPricer(
    PAYOFF    *optPayoff, /* (I) option payoff structure */
    FPAYOFF   *payFunc,   /* (I) option payoff function  */
    FADATA    *fa,        /* (I) measure data            */
    double    *premium    /* (O) fwd premium/rate        */
    )
{
    static char routine[] = "Q3FAGridPricer";
    int         status    = FAILURE;

    double grid[Q3_1D_INT_PTS], yield[Q3_1D_INT_PTS], dens[Q3_1D_INT_PTS];
    double densNorm;

    /* calculate FA yield density */
    if (Q3FADens(            
        fa,
        -Q3_1D_NUM_STDEV,   /* start */
        Q3_1D_NUM_STDEV,    /* end   */
        Q3_1D_INT_PTS,
        &(grid[0]),
        &(yield[0]),
        &(dens[0]),
        &densNorm) == FAILURE) goto RETURN;

    /* integrate density * payoff */
    if (Q3SimpsonPricer1D(
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

} /* Q3FAGridPricer */


/* JCP: change to credit measure change (Eannuity/Annuity(L) * Z(L) */
/*f----------------------------------------------------------------------------
 * Q3FADens
 *
 * Generates convexity and delay adjusted density
 *
 */
int Q3FADens(
    FADATA *fa,         /* (I) measure data              */
    double  start,      /* (I) start in std normal space */
    double  end,        /* (I) end in std normal space   */
    long    numGridPts, /* (I) num pts to generate       */
    double *grid,       /* (O) grid in normal space      */
    double *yields,     /* (O) grid in yield space       */
    double *faDens,     /* (O) density                   */
    double *normC       /* (O) density normalization     */
    )
{
    static char routine[] = "Q3FADens";
    int         status    = FAILURE;

	double zero, p;
    MQDATA *mq        = fa->mq;

    /* low rate cutoff for Ann and Del functional forms */
    double yCutoff = Q3_FA_BDRY * fa->mq->fwdRate;

    double step, xval, y,  yCut;
    long   i;

    step  = (end - start)/(numGridPts-1);
    xval  = start - step;

    /* calculate mapping and FA density function */
    for (i = 0; i < numGridPts; i++)
    {
        xval += step;
        grid[i] = mq->sigMQ * xval *sqrt(mq->optExpy); /* + mq->muMQ;   */
        if (Q3MQMapCR(
            mq,
            xval,
            &(yields[i])) == FAILURE) goto RETURN;
        
        y = yields[i];
        if (y > TINY)
        {
            yCut = MAX (y, yCutoff);
            p    = MIN (y / yCutoff, 1.);
            
            zero = p*Zero(yCut, fa) + (1-p);
        } 
        else {
            zero = 1;
            yCut = yCutoff;
        }
        
        faDens[i] = NormDens(xval) * zero 
            * (yCut/(1-fa->recovery)) 
            / Protection(yCut, fa->IRate, fa);
    }
    
    /* normalize distribution */
    *normC = Q3SimpsIntegral(
        step,
        numGridPts,
        faDens);
    if (*normC < TINY) 
    {
        DR_Error("%s: Non-positive total FA mass.\n", routine);
        goto RETURN;
    }
    
    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* Q3FADens */


/*f----------------------------------------------------------------------------
 * Q3FATailStrikes
 *
 * Calculate strikes corresponding to tail probabilities.
 * 
 */
int Q3FATailStrikes(
    FADATA *fa,        /* (I) measure data       */
    double *tailProb,  /* (I) tail probabilities */
    double *tailStrike /* (O) tail strikes       */
    )
{
    int status = FAILURE;

    double grid[Q3_1D_INT_PTS], yields[Q3_1D_INT_PTS], dens[Q3_1D_INT_PTS];
    double loProb, hiProb, cumProb, normDens;

    long   i, iLo, iHi;

    /* calculate FA yield density */
    if (Q3FADens(            
        fa,
        -Q3_1D_NUM_STDEV,   /* start */
        Q3_1D_NUM_STDEV,    /* end   */
        Q3_1D_INT_PTS,
        &(grid[0]),
        &(yields[0]),
        &(dens[0]),
        &normDens) == FAILURE) goto RETURN;

    /* find strikes for tail probabilities */
    cumProb = 0.;
    hiProb  = tailProb[0] * normDens / Q3_1D_STEP;
    loProb  = tailProb[1] * normDens / Q3_1D_STEP;        
    iLo     = 0;
    iHi     = Q3_1D_INT_PTS - 1;

    for (i = 0; i < Q3_1D_INT_PTS; i++) 
    {
        /* tail probs */
        cumProb += dens[i];
        if (cumProb < loProb)
            iLo = i;
        else if (cumProb < hiProb)
            iHi = i;
    }

    tailStrike[0] = yields[iHi];
    tailStrike[1] = yields[iLo];

    status = SUCCESS;

 RETURN:

    return status;

} /* Q3FATailStrikes */


/*f----------------------------------------------------------------------------
 * Q3FASmileInit
 * 
 * Determines the functional form for the change of measure between the
 * FA and AA measures (convexity and delay), based on VNFM
 *
 */
int Q3FASmileInit(
    double  expiry,        /* (I) observation date in yrs      */
    double  sigATM,         /* (I) ATM vol                      */
    double  start,         /* (I) swap/zero rate start in yrs  */
	long    freqSwap,
    long    freqCDS,       /* (I) swap/zero rate frequency     */
    double  swapMat,       /* (I) swap rate tenor in yrs       */
    double  swapRate,      /* (I) par swap rate                */
	double  cdsParSpread,
	double  recovery,
    double  fwdAnnuity,    /* (I) forward (swap) annuity       */
    double  zeroRateSwap,  /* (I) zero rate for same interval  */
    double  payDelay,      /* (I) payment delay                */   
    double  zeroRatePay,   /* (I) zero rate for delay interval */    
    long    numVnfmParams, /* (I) number of vnfm params        */
    double *vnfmParams,    /* (I) vnfm model parameters        */
    long    cashPhysSetl,  /* (I) settlement type              */
    FADATA *fa             /* (I/O) FA data                    */
    )
{   
    static char      routine[]="Q3FASmileInit";
    int              status = FAILURE; 

    double zeroSwapVolRatio, swapVol, zeroVol, zeroSwapCorr;

    /* rate frequencies and maturities */
    fa->freqAnn = freqCDS;
	fa->freqSwap = freqSwap;

    fa->matAnn = swapMat;
    fa->matDel = payDelay;

    /* if cash settled, alpha = power = 1; otherwise, use VNFM */
    if (cashPhysSetl == Q3_CASH_SETL) 
    {
        /* no annuity adjustment needed */
        fa->alphaAnn = fa->powerAnn = 1.;
    } 
	else
	{
		if (Q3VNFMZero2SwapCR(
            expiry,
			0.,
            start,
            freqCDS,
            swapMat,
            swapRate,
			recovery,
            fwdAnnuity,
            swapMat,
            zeroRateSwap,
            numVnfmParams,
            vnfmParams,
            &swapVol,
            &zeroVol,
            &zeroSwapCorr) == FAILURE) goto RETURN;

        zeroSwapVolRatio = zeroVol / swapVol;
        fa->powerAnn = zeroSwapCorr * zeroSwapVolRatio;
        fa->alphaAnn = zeroRateSwap / pow(cdsParSpread, fa->powerAnn) *
            exp(0.5 * sigATM * sigATM * expiry * fa->powerAnn *
                (1. - fa->powerAnn));      
    }

    /* VNFM calculation: functional form for delay adjustment */
    if (cashPhysSetl == Q3_CASH_SETL || payDelay < TINY) 
    {
        /* no delay adjustment needed */
        fa->powerDel = fa->alphaDel = 1.;
    } 
    else 
    {
	  if (Q3VNFMZero2SwapCR(
            expiry,
			0.,
            start,
            freqCDS,
            swapMat,
            swapRate,
			recovery,
            fwdAnnuity,
            payDelay,
            zeroRatePay,
            numVnfmParams,
            vnfmParams,
            &swapVol,
            &zeroVol,
            &zeroSwapCorr) == FAILURE) goto RETURN;

        zeroSwapVolRatio = zeroVol / swapVol;
        fa->powerDel = zeroSwapCorr * zeroSwapVolRatio;
        fa->alphaDel = zeroRatePay / pow(cdsParSpread, fa->powerDel) *
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

} /* Q3FASmileInit */

#define  BOND_PRICE(c,f,m,y) \
         (fabs(y)<TINY) ? (1 + (m) * (c)) : ((c) / (y) + \
         (1- (c) / (y)) * pow(1 + (y)/(f), -(m) * (f))) 
/*-----------------------------------------------------------------------------
 * Q3VNFMZero2Swap
 * 
 * Given VNFM parameters (weights, mean reversions, correlations), evaluates
 * annualized zero and swap vols and the correlation between coupon and zero 
 * rates.
 *
 */
int Q3VNFMZero2SwapCR(
    double  expiry,        /* (I) observation date in yrs         */
    double  volStart,      /* (I) start of vol observation in yrs */
    double  start,         /* (I) swap/zero rate start in yrs     */
    long    freq,          /* (I) swap/zero rate frequency        */
    double  swapMat,       /* (I) swap rate tenor in yrs          */
    double  swapRate,      /* (I) par swap rate                   */
	double  recovery,
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
    static char      routine[]="Q3VNFMZero2Swap";
    int              status = FAILURE; 

    double beta[3], wght[3], rho[3][3];
    double swapBCoeff[3], zeroBCoeff[3];
    double sVol, zVol, coVar;
    long   numFact;
    int    i,j;

    recovery = recovery; /* removed unused parameter warning */

    /* check swpMat, zeroMat */
    if (swapMat < TINY || zeroMat < TINY) goto RETURN;

    /* read VNFM parameters */
    switch(numVnfmParams) 
    {
    case 1:         /* 1 fact, alpha not given */
        numFact   = 1;
        beta[0]   = MAX(vnfmParams[0], Q3_BETA_SHIFT);
        wght[0]   = 1.;
        rho[0][0] = 1.;
        break;
    case 2:         /* 1 fact, alpha given */
        numFact   = 1;
        beta[0]   = MAX(vnfmParams[0], Q3_BETA_SHIFT);
        wght[0]   = vnfmParams[1];
        rho[0][0] = 1.;
        break;
    case 5:         /* 2 fact */
        numFact   = 2;
        beta[0]   = MAX(vnfmParams[0], Q3_BETA_SHIFT);
        beta[1]   = MAX(vnfmParams[1], Q3_BETA_SHIFT);
        wght[0]   = vnfmParams[2];
        wght[1]   = vnfmParams[3];
        rho[0][0] = rho[1][1] = 1.;
        rho[1][0] = rho[0][1] = vnfmParams[4];
        break;
    case 9:         /* 3 fact */
        numFact   = 3;
        beta[0]   = MAX(vnfmParams[0], Q3_BETA_SHIFT);
        beta[1]   = MAX(vnfmParams[1], Q3_BETA_SHIFT);
        beta[2]   = MAX(vnfmParams[2], Q3_BETA_SHIFT);
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

} /* Q3VNFMZero2Swap */



