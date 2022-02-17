/***************************************************************
 * Module:        crxvnfm
 * Submodule:     crxvnfm.c
 * File:          crxvnfm.c
 * Function:      Implementation of VNFM for credit
 * Author:        Charles Morcom, based on original CRX code by
 *                Dapeng Guan and others.
 ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include <alib/cashflow.h>
#include <alib/ldate.h>
#include <alib/tcurve.h>
#include <alib/gtobf.h>

#include <crxflow/include/crxerror.h>
#include <crxflow/include/optbsq.h>

#ifndef    MIN_SPOT_VOL_RATIO
#define    MIN_SPOT_VOL_RATIO    0.2
#endif

/**  Exponential decay function. */
double CrxExpDecay(double a, double t)
{ 
    double at = a*t;
    return (fabs (a*t) < 1E-10 ? 1e0 : (1e0 - exp(- at)) / (at)); 
}

/******************************************************************************
 * CrxBinarySearchLong
 * Performs binary search on an array of longs (a.k.a TDate) to
 * find the location of x in an ordered array aryX. If
 * x matches aryX exactly then x=aryX[pos]. If not,
 * r = -pos-1 is the index in aryX such that
 * aryX[r]>x, or r=n if x>aryX[n-1].
 * The only reason why this function ever returns FAILURE is if n<=0. Even
 * if the array is out of order, it will still return a value: it just
 * won't mean very much. You must ensure that aryX is sorted.
 *****************************************************************************/
int CrxBinarySearchLong(
    long        x,      /**<(I) Value to search for                          */
    int         n,      /**<(I) Array size                                   */
    const long* aryX,   /**<(I) Ordered array of values to search            */
    int*        pos     /**<(O) Location of x in xAry                        */    
    ) {

    static const char * function = "CrxBinarySearch";
    int  status = FAILURE;
    long xLo, xHi, xMid; // x-values for search
    int  iLo, iHi, iMid; // array indices for search

    if (n<=0)
    {
        DR_Error(function,
            "The array size (%d) should be strictly positive", n);
        goto RETURN;
    }

    // initialize and check within array bounds
    iLo = 0;
    xLo = aryX[iLo];
    iHi = n-1;
    xHi = aryX[iHi];
    if (x<xLo) 
    {
        /* x<aryX[0] */
        *pos = -1;
        status = SUCCESS;
        goto RETURN;

    } else if (x==xLo)
    {
        /* x=aryX[0] */
        *pos = 0;
        status = SUCCESS;
        goto RETURN;

    } else if (x>xHi)
    {
        /* x>aryX[n-1] */
        *pos = -n-1;
        status = SUCCESS;
        goto RETURN;

    } else if (x==xHi)
    {
        /* x=aryX[n-1] */
        *pos = n-1;
        status = SUCCESS;
        goto RETURN;

    }

    while (iHi-iLo>1) {
        /* You know xm<x<xp: halve interval size and try again, 
           checking midpoint */
        iMid = (iLo + iHi)/2; // integer division
        xMid = aryX[iMid];
        if (xMid==x)
        {
            *pos = -iMid-1;
            status = SUCCESS;
            goto RETURN;
        } else if (xMid>x)
        {
            /* Assume x in (xLo, xMid) */
            xHi = xMid;
            iHi = iMid;
        } else 
        {
            /* Assume x in (xMid, xHi) */
            xLo = xMid;
            iLo = iMid;
        }
    }

    /* Now you know that iHi is the lowest index s.t. x<aryX[iHi] */
    *pos = -iHi - 1;
    status = SUCCESS;

RETURN:
    return status;
}

/*****************************************************************************
 * CrxInterpLinearSquaredLong
 * Interpolates linearly on the square of the y-values - equivalent
 * to squaring all the y-values, linearly interpolating, and returning
 * the positive square-root of the result.
 *****************************************************************************/
int CrxInterpLinearSquaredLong(
    long            x,      /**<(I)Value to interpolate                      */
    int             n,      /**<(I) No. of x and y values                    */
    const long*     aryX,   /**<(I) Array of x-values                        */
    const TDate*    aryY,   /**<(I) Array of y-values                        */
    double*         result  /**<(O) Interpolated output y-value              */
                           ) {

    static const char * function = "CrxInterpLinearSquaredLong";
    int  status = FAILURE;
    int xLoc;
    long xLo, xHi;
    double yLo, yHi;
    
    if (FAILURE==CrxBinarySearchLong(x, n, aryX, &xLoc))
    {
        DR_Error(function, "Binary search in x-values failed.");
        goto RETURN;
    }

    if (xLoc==-1)
    {
        // flat extrapolation at low end
        *result = fabs(aryY[0]);
        status = SUCCESS;

    } else if (xLoc==-n-1) 
    {
        // flat extrapolation at high end
        *result = fabs(aryY[n-1]);
        status = SUCCESS;

    } else if (xLoc>=0)
    {
        // exact match
        *result = fabs(aryY[xLoc]);
        status = SUCCESS;
    } else 
    {
        xLo = aryX[-xLoc-2];
        xHi = aryX[-xLoc-1];
        yLo = aryY[-xLoc-2];
        yHi = aryY[-xLoc-1];
        if (xLo>=xHi) 
        {
            DR_Error(function,
                "The array of x-values is either not sorted, or not distinct at (%ld,%ld)",
                xLo,xHi);
        } else
        {    
            // linearly interpolate on the squares
            *result = sqrt(((xHi-x)*yLo*yLo + (x-xLo)*yLo*yLo)/(xHi-xLo));
            status = SUCCESS;
        }
    }

RETURN:
    return status;
}

/******************************************************************************
 * CrxVNFMQFactors
 * Calculates VNFM 'Q' factors for a specified CDS
 * Q[0] is the IR Q-component
 * Q[1] is the credit Q-component
 * Also computes the par spread and the annuity. 
 *****************************************************************************/
int CrxCreditVNFMQFactors(
    TDate         valueDate, /**<(I) CDS value Date                          */
    TDate         cdsStart,  /**<(I) Start date of CDS                       */
    TDate         cdsEnd,    /**<(I) Maturity of CDS                         */
    long          cdsDCC,    /**<(I) ALIB DCC for CDS coupon                 */
    TDateInterval cdsFreq,   /**<(I) CDS coupon Frequency                    */
    double        recovery,  /**<(I) Recovery Rate                           */
    double        betaIR,    /**<(I) Interest rate mean revers.              */
    double        backboneIR,    /**<(I) IR backbone: 0=normal, 1=lognormal  */
    const TCurve* irCurve,   /**<(I) IR Discount curve                       */
    double        betaCR,    /**<(I) Credit mean reversion                   */
    double        backboneCR,    /**<(I) IR backbone: 0=normal, 1=lognormal  */
    const TCurve* crCurve,   /**<(I) Credit clean-spread curve               */
    double*       Q,         /**<(O) credit VNFM q-factors                   */
    double*       parSpread, /**<(O) Forward par yield                       */
    double*       annuity    /**<(O) Annuity                                 */
    )
{    

    /*=========================================================================
     * VARIABLE DECLARATIONS - N.B. COMMENTS WHERE INITIALIZED
     *=======================================================================*/
    static const char* function = "CrxCreditVNFMCDSQ";
    int status = FAILURE;
    double irDiscToValue, irDisc, irDiscLast, irDiscStart;
    double crSProbToValue, crSProb, crSProbLast, crSProbStart, S, prevT;
    double ircrDisc, annuityL, protection, AIR, AIRLAST, BIR, CIR, ACR, BCR, CCR, DCR;
    int k;
    TDate prevCpnDate;

    /*=========================================================================
     * CREATE COUPON DATES AND PAYMENTS
     *=======================================================================*/
    TCashFlowList* cdsCpnCFL = GtoMakeCFL(
        1.0, /* Coupon rate */
        cdsStart,
        &cdsFreq,
        cdsEnd,
        cdsDCC,
        GTO_STUB_SIMPLE,
        FALSE, /* Stub at front */
        0, /* Don't add or subtract initial/final */
        GTO_BAD_DAY_NONE, /* Accruals */
        GTO_BAD_DAY_NONE, /* Payments */
        NULL /* No holidays */
        );
    if (cdsCpnCFL==NULL) 
    {
        DR_Error(function, "Unable to generate CDS cashflows from %ld to %ld", 
            cdsStart, cdsEnd);
    }

    /*=========================================================================
     * SET UP INITIAL DAYCOUNTS/DISCOUNTS/SURVIVAL PROBABILITIES ETC
     *=======================================================================*/
    /* DISCOUNT TO VALUE DATE */
    if (FAILURE==GtoInterpPV(valueDate, (TCurve*)irCurve, 
        GTO_FLAT_FORWARDS, &irDiscToValue)
        || irDiscToValue<=0)
    {
        DR_Error("%s failed: unable to compute discount to value date %ld", 
            function, valueDate);
        goto RETURN;
    }
    /* SURVIVAL PROBABILITY TO VALUE DATE */
    if (FAILURE==GtoInterpPV(valueDate, (TCurve*)crCurve, 
        GTO_FLAT_FORWARDS, &crSProbToValue)
        || crSProbToValue<=0)
    {
        DR_Error("%s failed: can't compute survival prob. to value date %ld", 
            function, valueDate);
        goto RETURN;
    }
    /* TIME TO CDS START - S */
    if (FAILURE==GtoDayCountFraction(valueDate, cdsStart, GTO_ACT_365, &S)) 
    {
        DR_Error("%s failed: can't compute CDS DCF to start %ld to %ld", 
            function, valueDate, cdsStart);
        goto RETURN;
    }
    prevT = S;
    /* DISCOUNT FACTOR TO S (START TIME OF CDS) */
    if (FAILURE==GtoInterpPV(cdsStart, (TCurve*)irCurve, 
        GTO_FLAT_FORWARDS, &irDiscStart)
        || irDiscStart<=0)
    {
        DR_Error("%s failed: Unable to compute +ve DF to CDS start date %ld",
            function, cdsStart);
        goto RETURN;
    }
    /* SURVIVAL PROBABILITY TO S */
    if (FAILURE==GtoInterpPV(cdsStart, (TCurve*)crCurve, 
        GTO_FLAT_FORWARDS, &crSProbStart)
        || crSProbStart<=0)
    {
        DR_Error("%s failed: Can't compute +ve s. prob. to CDS start date %ld",
            function, cdsStart);
        goto RETURN;
    }
    irDiscLast = irDiscStart;  // keep as previous value in loop
    crSProbLast = crSProbStart;// keep as previous value in loop

    /*=========================================================================
     * Coupon and protection contributions to A,B VNFM factors and annuity. 
     * Note that the annuity value does not include recovery of accruals
     * on default.
     * Note also that the integration interval for the protection leg is
     * implicitly taken to be the same as the coupon frequency. This is faster
     * and, given the other approximations, should be good enough since, in 
     * practice, it would probably be set to 3M or 6M anyway which are likely
     * to be the same as the coupon frequency.
     *=======================================================================*/
    annuityL = 0.0;
    protection = 0.0;
    AIR = 0.0; // various components of CR and IR Q factor FOR CDS
    AIRLAST = 0.0;
    BIR = 0.0;
    ACR = 0.0;
    BCR = 0.0;
    CIR = 0.0;
    CCR = 0.0;
    DCR = 0.0;
    /* prevCpnDate is used to calculate the mid-point for the accrual 
       recoveries. Assuming hasn't happened by the valueDate, it makes
       sense to go halfway between max(now,cpnDate) and cpnDate in terms
       of making an intelligent approximation to the default time */
    prevCpnDate = MAX(valueDate, cdsStart);
    for(k = 0; k < cdsCpnCFL->fNumItems; k++)
    {
        TDate cpnDate = cdsCpnCFL->fArray[k].fDate;
        double cpnAmt = cdsCpnCFL->fArray[k].fAmount;
        double cpnDCF = 0.0; // true day-count of coupon
        double timeToCpn = 0.0; // years from valueDate to coupon payment
        double cpnPV = 0.0, protPV = 0.0, acrTemp;
        double irMidDisc; // Zr to the middle of the cpn pd.
        double AIRMID = 0.0;
        TDate cpnMidDate = (cpnDate+prevCpnDate)/2;

        /* SKIP COUPONS WHICH HAVE PASSED */
        if (cpnDate<valueDate) continue;
        
        /* Calculate survival probability and discount factor at end of period */
        if (FAILURE==GtoInterpPV(cpnDate, (TCurve*)irCurve, 
            GTO_FLAT_FORWARDS, &irDisc)
            || irDisc<=0)
        {
            DR_Error("%s failed: Can't compute +ve DF to coupon %lf on %ld",
                function, cpnAmt, cpnDate);
            goto RETURN;
        }
        if (FAILURE==GtoInterpPV(cpnDate, (TCurve*)crCurve, 
            GTO_FLAT_FORWARDS, &crSProb)
            || crSProb<=0)
        {
            DR_Error("%s failed: Can't compute +ve s. prob. to cpn %lf on %ld",
                function, cpnAmt, cpnDate);
            goto RETURN;
        }
        irDisc /= irDiscToValue;
        crSProb /= crSProbToValue;

        /* Compute discount to middle of period */
        if (FAILURE==GtoInterpPV(cpnMidDate, (TCurve*)irCurve, 
            GTO_FLAT_FORWARDS, &irMidDisc)
            || irMidDisc<=0)
        {
            DR_Error("%s failed: Can't compute +ve DF to coupon period midpoint %lf on %ld",
                function, cpnMidDate);
            goto RETURN;
        }
        irMidDisc /= irDiscToValue;

        // calculate time to coupon and true (A/365) coupon day-count fraction
        if (FAILURE==GtoDayCountFraction(valueDate, cpnDate, 
            GTO_ACT_365, &timeToCpn)) 
        {
            DR_Error("%s failed: Can't compute cpn DC frac. from %ld to %ld",
                function,valueDate,cpnDate);
            goto RETURN;
        }
        cpnDCF = timeToCpn - prevT;
        

        ircrDisc = irDisc*crSProb;
        cpnPV = cpnAmt * ircrDisc;
        protPV = (log(crSProb/crSProbLast)/log(crSProb*irDisc/(crSProbLast*irDiscLast))) 
            * (crSProbLast*irDiscLast - crSProb*irDisc);
        protection += protPV;
        /* ANNUITY IS SCHEDULED AMOUNT + PROT ON HALF COUPON */
        annuityL += cpnPV + protPV * cpnAmt/2.0;
        
         /*=====================================================================
         * VNFM A,B,C,D factor computations for IR and CR
         *===================================================================*/
        /* MID COMPUTATION FIRST FOR B VALUES */
        /* BIR += AR.CPNPV + 0.5 CPN.ARMID.PROT */
        /** AIR(u) = \int_S^u r exp(-beta(t-S))dt is the VNFM 'A' 
            ACR is analogously defined for credit */
        AIRMID = AIR + exp (-betaIR * (prevT-S)) 
            * (backboneIR*log (irDiscLast/irMidDisc) + (1-backboneIR)/cpnDCF)
            * CrxExpDecay (betaIR, cpnDCF/2.0);
        AIR += exp (-betaIR * (prevT-S)) 
            * (backboneIR*log (irDiscLast/irDisc) + (1-backboneIR)/cpnDCF)
            * CrxExpDecay (betaIR, cpnDCF);
        
        acrTemp = exp (-betaCR * (prevT-S)) 
            * (backboneCR*log (crSProbLast/crSProb) + (1-backboneCR)/cpnDCF) 
            * CrxExpDecay (betaCR, cpnDCF) ;
        ACR  += acrTemp;
        /** payment contribution to vol. This is VNFM 'B' 
            with analogous definition for CR */
        /* FEE LEG */
        BIR  += AIR * cpnPV + AIRMID * protPV*cpnAmt/2.0;
        BCR  += ACR * cpnPV + (AIRLAST*crSProbLast - AIR*crSProb)*irMidDisc*cpnAmt/2.0;
        /* PROTECTION LEG */
        CIR += AIR * protPV;
        CCR += ACR * protPV;
        DCR -= exp (-betaCR * (prevT-S)) * CrxExpDecay (betaCR, cpnDCF) * protPV;


        AIRLAST = AIR;
        prevT = timeToCpn;
        irDiscLast = irDisc;
        crSProbLast = crSProb;
    }

    BIR /= annuityL;
    BCR /= annuityL;    
    CIR /= protection;
    CCR /= protection;
    DCR /= protection;

    Q[0] = BIR - CIR;
    Q[1] = BCR - CCR - DCR;
    *parSpread = (1.0 - recovery) * protection / annuityL;
    *annuity  = annuityL;

    status = SUCCESS;

 RETURN:
    GtoFreeCFL(cdsCpnCFL);
    return status;
}  /* QFactor */

/******************************************************************************
 * CrxCreditVNFMBlackScholesVol
 * Given credit and interest rate spot vols, mean reversions and
 * correlations, returns the Black-Scholes implied volatility of a given
 * swaption.
 *****************************************************************************/
int CrxCreditVNFMBlackScholesVol (
    TDate valueDate,
    TDateInterval cdsFreq,       /**<(I) Underlying BM CDS fee frequency     */
    long          cdsDCC,        /**<(I) Underlying BM CDS fee day-count     */
    double        recovery,      /**<(I) Recovery Rate                       */
    TDate         cdsOptionExpiry, /**<(I) Vol expiry date for swaption      */
    TDate         cdsStartDate,  /**<(I) Underlying BM CDS start dates       */
    TDate         cdsEndDate,    /**<(I) Underlying BM CDS maturity dates    */
    /* CREDIT VOL AND CURVE PARAMETERS */
    int           numVolCR,      /**<(I) Num. of input credit vol points     */
    const TDate*  volDateCR,     /**<(I) Credit spot vol dates               */
    double*       volCR,         /**<(I) Credit spot vols                    */
    double        betaCR,        /**<(I) Credit mean-reversion               */
    double        backboneCR,    /**<(I) CR backbone: 0=normal, 1=lognormal  */
    double        factorWeightCR,/**<(I) CR alpha factor weight              */
    const TCurve* crCurve,       /**<(I) Credit clean spread curve           */
    double        leftCR2Q,      /**<(I) Left Credit 2Q smile (1=lognormal)  */
    double        rightCR2Q,     /**<(I) Right Credit 2Q smile (1=lognormal) */
    double        CR2QF,         /**<(I) Forward shift 'F' parameter for 2Q  */
    /* IR VOL AND CURVE PARAMETERS */
    int           numVolIR,      /**<(I) Num. of IR spot vol points          */
    const TDate*  volDateIR,     /**<(I) IR vol dates                        */
    const double* volIR,         /**<(I) IR spot vols                        */
    double        betaIR,        /**<(I) IR mean-reversion                   */
    double        backboneIR,    /**<(I) CR backbone: 0=normal, 1=lognormal  */
    double        factorWeightIR,/**<(I) CR alpha factor weight              */
    const TCurve* irCurve,       /**<(I) IR zero curve                       */
    /* IR/CR CORRELATION */
    double        corrIRCR,      /**<(I) IR/CR spot-vol correlation          */
    double*       bsVol          /**<(O)Output Black-Scholes volatility      */
    ) {

    static const char* function = "CrxCreditVNFMBlackScholesVol";
    int status = FAILURE;
    double Q[2]; // array to store IR/CR Q factors
    double T; // A/365 time to expiry of option
    double parSpread;
    double annuity, atmPx;
    double irQ, crQ, irVar, crVar, ircrCovar, lastT, totalVol, A, fwd2Q;
    int iP, cP, i;

    /*=========================================================================
     * MERGE DATES AND CREATE IR VOLS TO MATCH CR VOL DATES
     * Vol interpolation is linear in the variance, so 0.5 power in the
     * the vol. This interpolation is only for IR spot vols. The credit
     * vol interpolation is flat, so that the benchmarks stand some chance
     * of being consistent!
     * create array big enough for dates even if all different
     *=======================================================================*/
    int numVolIRCR = 0;
    double* volIRMerged = 0;
    double* volCRMerged = 0;
    double alphaIR_2, alphaCR_2;
    TDate* volDateMerged = (TDate*)calloc((1+numVolIR+numVolCR), sizeof(TDate));
    TDate lastDate = valueDate;
    if (!volDateMerged){
        DR_Error("%s failed: can't allocate memory for %d combined dates",
            function, (1+numVolIR+numVolCR));
        goto RETURN;
    }
    volIRMerged = (double*) calloc((1+numVolIR+numVolCR), sizeof(double));
    if (!volIRMerged)
    {
        DR_Error("%s failed: can't allocate memory for %d combined IR vols",
            function, (1+numVolIR+numVolCR));
        goto RETURN;
    }
    volCRMerged = (double*) calloc((1+numVolIR+numVolCR), sizeof(double));
    if (!volCRMerged)
    {
        DR_Error("%s failed: can't allocate memory for %d combined CR vols",
            function, (1+numVolIR+numVolCR));
        goto RETURN;
    }
    /* A FEW MORE SIMPLE CHECKS */
    if (lastDate>=cdsOptionExpiry) {
        DR_Error("%s failed: option expiry date (%d) is after value date (%d).",
            function,cdsOptionExpiry,valueDate);
        goto RETURN;
    } else if (numVolIR<1 || numVolCR<1) {
        DR_Error("%s failed: You must provide at least one CR and at least one IR spot vol (no. IR=%d, no. CR=%d).",
            function, numVolIR, numVolCR);
        goto RETURN;
    }
    iP=0;
    cP=0;
    while (lastDate<cdsOptionExpiry && (iP<numVolIR || cP<numVolCR)) {
        if (iP==numVolIR) {
            // no more IR dates - consume CR dates
            // keep IR vols from last IR spot vol
            volDateMerged[numVolIRCR] = volDateCR[cP];
            volIRMerged[numVolIRCR] = volIR[numVolIR-1];
            volCRMerged[numVolIRCR] = volCR[cP];
            cP++;
        } else if (cP==numVolCR) {
            // no more CR dates - consume IR dates, carry forward CR vols
            volDateMerged[numVolIRCR] = volDateIR[iP];
            volIRMerged[numVolIRCR] = volIR[iP];
            volCRMerged[numVolIRCR] = volCR[numVolCR-1];
            iP++;
        } else {
            TDate irDate = volDateIR[iP];
            TDate crDate = volDateCR[cP];
            if (irDate==crDate) {
                // both same - take each & increment both
                volDateMerged[numVolIRCR] = irDate;
                volIRMerged[numVolIRCR] = volIR[iP];
                volCRMerged[numVolIRCR] = volCR[cP];
                iP++;
                cP++;
            } else if (irDate<crDate) {
                // consume IR date - before CR date
                // CR vol is now in the next interval: pick up from this, not last
                volDateMerged[numVolIRCR] = irDate;
                volIRMerged[numVolIRCR] = volIR[iP];
                volCRMerged[numVolIRCR] = volCR[cP];
                iP++;
            } else {
                // consume CR date - lower than IR date
                volDateMerged[numVolIRCR] = crDate;
                volCRMerged[numVolIRCR] = volCR[cP];
                /* Interpolate vol on variance, not vol */
                if (iP==0) {
                    volIRMerged[numVolIRCR] = volIR[0];
                } else {
                    TDate lastIRDate = volDateIR[iP-1];
                    TDate nextIRDate = volDateIR[iP];
                    double lastIRVar = volIR[iP-1]*volIR[iP-1];
                    double nextIRVar = volIR[iP]*volIR[iP];
                    volIRMerged[numVolIRCR] = sqrt(((nextIRDate-crDate)*lastIRVar + 
                        (crDate-lastIRDate)*nextIRVar)
                        /(double)(nextIRDate-lastIRDate));
                }
                cP++;
            }
        }
        lastDate = volDateMerged[numVolIRCR];
        numVolIRCR++;
    }
    /* MULTIPLY BY THE SQUARES OF THE FACTOR WEIGHTS */
    // slightly inefficient, but a little clearer, I think
    alphaIR_2 = factorWeightIR*factorWeightIR;
    alphaCR_2 = factorWeightCR*factorWeightCR;
    for (i=0; i<numVolIRCR; i++) {
        volIRMerged[i] *= alphaIR_2;
        volCRMerged[i] *= alphaCR_2;
    }
    /* NOW, MAKE SURE THAT THE LAST DATE IS THE EXPIRY OF THE OPTION */
    if (lastDate<cdsOptionExpiry) {
        int lastIdx = (numVolIRCR>0 ? numVolIRCR-1 : 0);
        TDate lastIRDate = volDateMerged[lastIdx];
        TDate nextIRDate = (iP<numVolIR ? volDateIR[iP] : cdsOptionExpiry);
        double lastIRVar = volIRMerged[lastIdx]*volIRMerged[lastIdx];
        double nextIRVar = (iP<numVolIR ? volIR[iP]*volIR[iP] : lastIRVar);
        volIRMerged[numVolIRCR] = sqrt(((nextIRDate-cdsOptionExpiry)*lastIRVar + 
                        (cdsOptionExpiry-lastIRDate)*nextIRVar)
                        /(double)(nextIRDate-lastIRDate));
        volDateMerged[numVolIRCR] = cdsOptionExpiry;
        volCRMerged[numVolIRCR] = volCRMerged[(cP<numVolCR ? cP : numVolCR-1)];
        numVolIRCR++;
    }


    /* Compute time from value to CDS start date */
    if (FAILURE==GtoDayCountFraction(
        valueDate, cdsOptionExpiry, GTO_ACT_365, &T)) {
        DR_Error("%s failed - can't compute d-c frac from %ld to %ld",
            function, valueDate, cdsOptionExpiry);
        goto RETURN;
    }
    
    /* Q factor */
    if(FAILURE==
        CrxCreditVNFMQFactors(
            valueDate,
            cdsStartDate,
            cdsEndDate,
            cdsDCC,
            cdsFreq,
            recovery,
            betaIR,
            backboneIR,
            irCurve,
            betaCR,
            backboneCR,
            crCurve,
            Q, 
            &parSpread,
            &annuity))
    {
        DR_Error("%s failed - could not calculate Q factors for CDS from %d to %d",
            function, cdsStartDate, cdsEndDate);
        goto RETURN;
    }
    irQ = Q[0];
    crQ = Q[1];
        
    /* NOW INTEGRATE THE TOTAL VARIANCE FROM VALUE DATE TO OPTION EXPIRY DATE */
    lastT = 0.0;
    irVar = 0;
    crVar = 0;
    ircrCovar = 0;
    for (i=0; i<numVolIRCR; i++) {
        double tI; // day-fraction for time-step
        if (FAILURE==
            GtoDayCountFraction(
                (i>0 ? volDateMerged[i-1] : valueDate),
                volDateMerged[i], GTO_ACT_365, &tI))
        {
            DR_Error("%s failed - can't compute DCF %d from %ld to %ld",
                function, i, 
                (i>0 ? volDateMerged[i-1] : valueDate),
                volDateMerged[i]);
            goto RETURN;
        }
        lastT += tI;
        irVar += volIRMerged[i]*volIRMerged[i]*exp(2*betaIR*lastT)*CrxExpDecay(2*betaIR, tI)*tI;
        crVar += volCRMerged[i]*volCRMerged[i]*exp(2*betaCR*lastT)*CrxExpDecay(2*betaCR, tI)*tI;
        ircrCovar += volCRMerged[i]*volIRMerged[i]*corrIRCR
            *exp((betaCR+betaIR)*lastT)*CrxExpDecay((betaCR+betaIR), tI)*tI;
    }
    irVar *= exp(-2*betaIR*T);
    crVar *= exp(-2*betaCR*T);
    ircrCovar *= exp(-(betaCR+betaIR)*T);
    totalVol = sqrt(irQ*irQ*irVar + 2*crQ*irQ*ircrCovar + crQ*crQ*crVar);

    /* NOW YOU NEED TO COVERT THE Q-VOL INTO A BS VOL */
    /* FIRST CALIBRATE THE Q DISTRIBUTION TO THE FORWARD */
    A=1.0;
    if (FAILURE==Crx2QForward(A,CR2QF,totalVol,leftCR2Q,rightCR2Q,&fwd2Q))
    {
        DR_Error("%s failed - cannot calibrate 2Q forward to "
            "total vol=%lf, qL=%lf, qR=%lf, F=%lf",
            totalVol, leftCR2Q, rightCR2Q, CR2QF);
        goto RETURN;
    }
    if (fwd2Q<DBL_EPSILON)
    {
        DR_Error("%s failed - cannot calibrate 2Q forward to "
            "total vol=%lf, qL=%lf, qR=%lf, F=%lf - negative fwd2Q %lf",
            totalVol, leftCR2Q, rightCR2Q, CR2QF, fwd2Q);
        goto RETURN;
    }
    A = parSpread/fwd2Q;

    // FIND THE 2Q ATM OPTION PRICE
    if (FAILURE==Crx2QOptionPrice(
        GtoOPTION_CALL,
        parSpread,
        A,    
        CR2QF,  
        totalVol,  
        leftCR2Q,  
        rightCR2Q,
        &atmPx)
        ) {
        DR_Error("%s failed - cannot find ATM 2Q option price to "
            "total vol=%lf, qL=%lf, qR=%lf, F, A=%lf, par spread %lf.",
            totalVol, leftCR2Q, rightCR2Q, CR2QF, A, parSpread);
        goto RETURN;
    }

    // invert to find BS vol
    if (FAILURE==BSImpVol(
           bsVol,
           parSpread,
           parSpread,
            T,
            atmPx,
            OPT_CALL,
            totalVol/sqrt(T))) {
        DR_Error("%s failed - compute ATM BS vol for "
            " strike=%lf, expiry=%lf, price=%lf.", parSpread, T, atmPx);
        goto RETURN;
    }
  
    status = SUCCESS;
  
 RETURN:
    if (volDateMerged) free(volDateMerged);
    if (volIRMerged) free(volIRMerged);
    if (volCRMerged) free(volCRMerged);
    return status;
}
            
/******************************************************************************
 * CrxCreditVNFMBootstrapSpotVols
 * Bootstraps a set of benchmark CDS and corresponding Black-Scholes vols
 * to an equivalent set of spot vols using the credit VNFM approximation.
 * The vol must be an implied BS vol - if your vols are multi-q, then you
 * should convert them before calling this function.
 * Note that the annuities and par spreads calculated are approximate - 
 * they do not include accrual recovery on default - so do not use them other
 * than as approximations.
 *****************************************************************************/
int CrxCreditVNFMBootstrapSpotVols (
    TDate valueDate,
    /* CREDIT VOL PARAMETERS */
    int           numVolCR,      /**<(I) Num. of input credit vol points     */
    const TDate*  volDateCR,     /**<(I) Credit vol dates                    */
    double*       volCR,         /**<(I/O) Input BS vols, output spot if cal */
    int*          volUsedCR,     /**<(I/O) Input/output 1 if should/did cal  */
    int           skipFlagCR,    /**<(I) If TRUE, skip & go on if fail BM cal*/
    TDateInterval cdsFreq,       /**<(I) Underlying BM CDS fee frequency     */
    long          cdsDCC,        /**<(I) Underlying BM CDS fee day-count     */
    double        recovery,      /**<(I) Recovery Rate                       */
    const TDate*  cdsStartDates, /**<(I) Underlying BM CDS start dates       */
    const TDate*  cdsEndDates,   /**<(I) Underlying BM CDS maturity dates    */
    double        betaCR,        /**<(I) Credit mean-reversion               */
    double        backboneCR,    /**<(I) CR backbone: 0=normal, 1=lognormal  */
    double        factorWeightCR,/**<(I) CR alpha factor weight              */
    const TCurve* crCurve,       /**<(I) Credit clean spread curve           */
    double        leftCR2Q,      /**<(I) Left Credit 2Q smile (1=lognormal)  */
    double        rightCR2Q,     /**<(I) Right Credit 2Q smile (1=lognormal) */
    double        CR2QF,         /**<(I) Forward shift 'F' parameter for 2Q  */
    double        corrIRCR,      /**<(I) IR/CR spot-vol correlation          */
    int           numVolIR,      /**<(I) Num. of IR spot vol points          */
    const TDate*  volDateIR,     /**<(I) IR vol dates                        */
    const double* volIR,         /**<(I) IR spot vols                        */
    double        betaIR,        /**<(I) IR mean-reversion                   */
    double        backboneIR,    /**<(I) IR backbone: 0=normal, 1=lognormal  */
    double        factorWeightIR,/**<(I) IR alpha factor weight              */
    const TCurve* irCurve,       /**<(I) IR zero curve                       */
    double*       parSpreadAry,  /**<(O)Output BM par spreads, where cal     */
    double*       annuityAry     /**<(O)Output BM annuity value, where cal   */
    )                  
     
{
    
    static const char* function = "CrxCreditVNFMBootstrapSpotVols";
    int status = FAILURE;
    double* volIRMerged = 0;
    int iP=0, cP=0, numVolIRCR=0;
    int idxIR, lastIdxIR, lastI, i;
    double lastExpY, lastSpotVol; //, bsATMPrice, A, fwd2Q, vol2Q;
    double varSumII, varSumCC, varSumIC;

    /*=========================================================================
     * MERGE DATES AND CREATE IR VOLS TO MATCH CR VOL DATES
     * Vol interpolation is linear in the variance, so 0.5 power in the
     * the vol.
     * create array big enough for dates even if all different
     *=======================================================================*/
    TDate* volDateIRMerged = (TDate*)calloc(numVolIR+numVolCR, sizeof(TDate));
    if (!volDateIRMerged)
    {
        DR_Error("%s failed: can't allocate memory for %d combined dates",
            function, (numVolIR+numVolCR));
        goto RETURN;
    }
    volIRMerged = (double*) calloc(numVolIR+numVolCR, sizeof(double));
    if (!volIRMerged)
    {
        DR_Error("%s failed: can't allocate memory for %d combined IR vols",
            function, (numVolIR+numVolCR));
        goto RETURN;
    }
    while (iP<numVolIR || cP<numVolCR) {
        if (iP==numVolIR) {
            // no more IR dates - consume CR dates
            // keep IR vols from last IR spot vol
            volDateIRMerged[numVolIRCR] = volDateCR[cP];
            volIRMerged[numVolIRCR] = volIR[numVolIR-1];
            cP++;
        } else if (cP==numVolCR) {
            // no more CR dates - consume IR dates
            volDateIRMerged[numVolIRCR] = volDateIR[iP];
            volIRMerged[numVolIRCR] = volIR[iP];
            iP++;
        } else {
            TDate irDate = volDateIR[iP];
            TDate crDate = volDateCR[cP];
            if (irDate==crDate) {
                // both same - doesn't matter & increment both
                volDateIRMerged[numVolIRCR] = irDate;
                volIRMerged[numVolIRCR] = volIR[iP];
                iP++;
                cP++;
            } else if (irDate<crDate) {
                // consume IR date - lower than CR date
                volDateIRMerged[numVolIRCR] = irDate;
                volIRMerged[numVolIRCR] = volIR[iP];
                iP++;
            } else {
                // consume CR date - lower than IR date
                volDateIRMerged[numVolIRCR] = crDate;
                /* Interpolate vol on variance, not vol */
                if (iP==0) {
                    volIRMerged[numVolIRCR] = volIR[0];
                } else {
                    TDate lastIRDate = volDateIR[iP-1];
                    TDate nextIRDate = volDateIR[iP];
                    double lastIRVar = volIR[iP-1]*volIR[iP-1];
                    double nextIRVar = volIR[iP]*volIR[iP];
                    volIRMerged[numVolIRCR] = sqrt(((nextIRDate-crDate)*lastIRVar + 
                        (crDate-lastIRDate)*nextIRVar)
                        /(double)(nextIRDate-lastIRDate));
                }
                cP++;
            }
        }
        numVolIRCR++;
    }
    /* Adjust for IR factor weight */
    for (i=0; i<numVolIRCR; i++) {
        volIRMerged[i] *= factorWeightIR*factorWeightIR;
    }

    /*=========================================================================
     * START BOOTSTRAPPING
     *=======================================================================*/
    idxIR = 0;
    lastIdxIR = 0;
    lastI = -1;
    lastExpY = 0; // time to expiry of last calibrated BM
    lastSpotVol = 0.0;
    varSumII = 0.0; // running total vol/covar integrals
    varSumCC = 0.0;
    varSumIC = 0.0;
    /* LOOP THROUGH CR BENCHMARKS AND BOOTSTRAP SPOT VOLS */
    for(i = 0; i < numVolCR; i++)
    {
        double Q[2]; // q-factor array
        double irQ;
        double crQ;
        double expY; // time from value to expiry
        double t; // time since last calibrated vol
        double y; // RHS 'constant' of quadratic eqn in s
        double a,b,c, discriminant, A, B, C;
        double newSpotVol, bsATMPrice, vol2Q, A2Q;

        if (FAILURE==GtoDayCountFraction(
            valueDate, volDateCR[i], GTO_ACT_365, &expY)) {
            DR_Error("%s failed - can't compute d-c frac from %ld to %ld",
                function, valueDate, volDateCR[i]);
            goto RETURN;
        }

        /* ONLY CALIBRATE INSTRUMENTS IF SPECIFIED */
        if (!volUsedCR[i]) continue;

        /* A/365 Time since last calibrated vol or base-date */
        t = (lastI>=0 ? expY - lastExpY : expY);

        /* Q factor */
        if(FAILURE==
            CrxCreditVNFMQFactors(valueDate,
                cdsStartDates[i],
                cdsEndDates[i],
                cdsDCC,
                cdsFreq,
                recovery,
                betaIR,
                backboneIR,
                irCurve,
                betaCR,
                backboneCR,
                crCurve,
                Q, 
                &parSpreadAry[i],
                &annuityAry[i]))
        {
            DR_Error("%s failed - could not calculate Q factors for BM idx %d.",
                function, i);
            goto RETURN;
        }
    
        irQ = Q[0];
        crQ = Q[1];
    
        /*=====================================================================
         * CONVERT QUOTED STANDARD BS VOLS INTO Q-SPACE VOLS FOR VNFM                 
         * ASSIGN TOTAL VOLATILITY (sigma*sigma*T) TO y                      
         *===================================================================*/
        y = volCR[i]*sqrt(expY); // use this for sT^0.5 just temporarily
        bsATMPrice = 2*parSpreadAry[i]*(GtoNormalCum(0.5*y)-0.5);
        if (bsATMPrice<=0) {
            DR_Error("%s failed - non positive BS ATM call price %lf "  
                "at time %lf and vol %lf", function, bsATMPrice, 
                expY, volCR[i]);
            goto RETURN;
        }

        /* CALIBRATE Q A AND sigmaQ=total vol PARAMETERS */
        if (FAILURE==Crx2QCalibrate(parSpreadAry[i],bsATMPrice,
            leftCR2Q,rightCR2Q,CR2QF,&A2Q,&vol2Q)) {
            DR_Error("%s failed - failed to calibrate 2Q A and sigmaQ "
                " forward=%lf, T=%lf, ATM px=%lf, lQ=%lf, rQ=%lf, muQ=%lf", 
                function, parSpreadAry[i], expY, bsATMPrice, 
                leftCR2Q, rightCR2Q, CR2QF);
            goto RETURN;
        }
        
        y = vol2Q*vol2Q; // this is the total forward variance = sigma^2 T
        /*=====================================================================
         * SOLVE QUADRATIC EQUATION TO FIND CURRENT SPOT VOL
         *===================================================================*/
        // y is the part of variance for the current time interval only
        // deduct variance due to earlier time intervals
        y -= irQ*irQ * varSumII * exp(-2.0 * betaIR * t) + 
              crQ*crQ * varSumCC * exp(-2.0 * betaCR * t) + 
              2 * irQ*crQ * /*corrIRCR **/ varSumIC * exp(-(betaCR + betaIR) * t); 
    
        /* The bootstrapped vol satisfies a quadratic equation
            y = as^2 + 2bs + c */
        a = CrxExpDecay(2.0 * betaCR, t) * t;// * exp(-2.0 * betaCR * );
        b = 0.0;
        c = 0.0;

        /* ADD IR VOL COMPONENTS - MAY BE SEVERAL IF DATES MISALIGNED */
        // step back in case failed calibration of last stage
        if (lastI<0) {
            idxIR = 0;
        } else {
            while (idxIR>0 && volDateIRMerged[idxIR]>volDateCR[lastI]) 
                idxIR--;
        }
        while (idxIR<numVolIRCR && volDateIRMerged[idxIR]<=volDateCR[i])
        {
            double tIR; // day-fraction for IR time-step
            if (FAILURE==
                GtoDayCountFraction(
                    (idxIR>0 ? volDateIRMerged[idxIR-1] : valueDate),
                    volDateIRMerged[idxIR], GTO_ACT_365, &tIR))
            {
                DR_Error("%s failed - can't compute IR DCF %d from %ld to %ld",
                    function, idxIR, 
                    (idxIR>0 ? volDateIRMerged[idxIR-1] : valueDate),
                    volDateIRMerged[idxIR]);
                goto RETURN;
            }
      
            b *= exp(-(betaCR + betaIR) * tIR);
            b += corrIRCR* volIRMerged[idxIR] 
                * CrxExpDecay(betaCR + betaIR, tIR) * tIR;
      
            c *= exp(-2.0 * betaIR * tIR);
            c += volIRMerged[idxIR] * volIRMerged[idxIR]
                * CrxExpDecay(2.0 * betaIR, tIR) * tIR;
      
            idxIR++;
        }

        A = a*crQ*crQ;
        B = b*crQ*irQ;
        C = c*irQ*irQ;
        discriminant = (B*B - A*(C-y));
        /* CHECK FOR IMAGINARY ROOTS - FAIL IF SO */
        if (discriminant<0)
        {
            DR_Error("%s failed - imaginary spot vol implied at i=%d",
                function, i);
            goto RETURN;
        }
        /* CHECK THAT NEW VOL IS POSITIVE - FAIL IF NOT */
        newSpotVol = (sqrt(discriminant) - B)/A;
        if (newSpotVol<TINY)
        {
            DR_Error("%s failed - non-positive spot vol %lf impled at i=%d", 
                function, newSpotVol, i);
            goto RETURN;
        }

        /* CHECK THAT VOL HAS NOT CHANGED TOO MUCH FROM LAST ONE 
         * IF SO, SKIP IT OR ABORT IF !skipFlagCR */
        if (lastI>=0 && 
            (newSpotVol/lastSpotVol<MIN_SPOT_VOL_RATIO 
            || newSpotVol/lastSpotVol>1./MIN_SPOT_VOL_RATIO)) {

            DR_Error("Skipping vol %d - vol too different from previous value", 
                i);
            volUsedCR[i] = 0; // mark that skipped by clearing this
            if (!skipFlagCR) {
                goto RETURN;
            } else {
                continue; // BM i loop
            }
            
        } else {
            /* CALIBRATION OK - SET SPOT VOL FOR ALL DATE INTERVALS UP TO i */
            lastI++;
            while (lastI<=i) 
            {
                volCR[lastI] = newSpotVol;
                lastI++;
            }
            lastI = i;
            /*  UPDATE VARIANCE INTEGRALS */
            varSumCC *= exp(-2.0 * betaCR * t);
            varSumCC += a * newSpotVol * newSpotVol; 
            varSumII *= exp(-2.0 * betaIR * t);
            varSumII += c;
            varSumIC *= exp(-(betaCR + betaIR) * t);
            varSumIC += b * newSpotVol;
            lastSpotVol = newSpotVol;
            lastExpY = expY;
        }
    }

    /* Adjust resulting credit vols for factor weight */
    for (i = 0; i < numVolCR; i++) {
        volCR[i] /= factorWeightCR*factorWeightCR;
    }
  
    status = SUCCESS;
  
 RETURN:
    if (volDateIRMerged) free(volDateIRMerged);
    if (volIRMerged) free(volIRMerged);
  
    return status;
} /* CrxCreditVNFMBootstrapSpotVols */

