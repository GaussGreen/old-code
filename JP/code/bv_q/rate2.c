/******************************************************************************
 * Module:	Q3
 * Submodule:
 * File: rate.c	
 * Function:	
 * Author:	Interest Rates DR
 * Revision:	$Header$
 *****************************************************************************/
#include <math.h>
#include <ctype.h>                      /* toupper */
#include <stdio.h>
#include <string.h>

#include "q3.h"

#include "cgeneral.h"
#include "macros.h"                     /* MAX, MIN */
#include "bastypes.h"
#include "cerror.h"
#include "ldate.h"                      /* GTO_ACT_365F */
#include "date_sup.h"                   /* GtoDateIntervalToFreq */
#include "convert.h"                    /* GtoFormatDate */
#include "tcurve.h"                     /* GtoConvertTCurveWrap */
#include "cashflow.h"                   /* GtoMakeCFL */
#include "presval.h"                    /* GtoCashFlowPVFromYield */
#include "optprop.h"                    /* GTO_OPTION_CALL, GTO_OPTION_PUT */ 
#include "zr2simp.h"                    /* GtoZerosToSimplePoint */
#include "zr2coup.h"                    /* GtoZerosToCouponsPoint */
#include "swaprate.h"                   /* GtoSwapRate2 */
#include "swappv.h"                     /* GtoSwapFloatPV */
#include "dtivlo.h"                     /* GtoDateIntervalNew*/

static int BusDayFrac(
           TDate    startDate,
           TDate    expiryDate,
           long     BusDaysPerYr,
           char     *HolidayFile,
           double   *BusDaysDisc);


/*f----------------------------------------------------------------------------
 * Q3SwapRateCalc2
 *
 * Swap rate, annuity, etc calculator from curve (for quasi-vanilla 1d)
 */
int Q3SwapRateCalc2 (
    long              VolBaseDate,       /* 1  (I) ATM vol base date       */
    long              SmileBaseDate,     /* 2  (I) Smile base date         */
    long              CurveBaseDate,     /* 3  (I) Index curve base date   */
    long              numIndxPts,        /* 4  (I) Nb of index points      */
    long              *indxDates,        /* 5  (I) Index dates             */
    double            *indxRates,        /* 6  (I) Act/365F index rates    */
    long              numDiscPts,        /* 7  (I) Nb of discount points   */ 
    long              *discDates,        /* 8  (I) Discount dates          */
    double            *discRates,        /* 9  (I) Act/365F discount rates */
    long              rateFreq,          /* 10 (I) rate frequency          */
    char              *rateDCC,          /* 11 (I) rate day count conv     */
    long              expiryDate,        /* 12 (I) rate reset date         */
    long              startDate,         /* 13 (I) rate effective date     */
    long              matDate,           /* 14 (I) rate maturity date      */
    long              payDate,           /* 15 (I) option payment date     */
    char              *HolidayVol,       /* 16 (I) holiday file name       */
    char              *BusVolDCC,        /* 17 (I) bus day disc convention */
    double            *output)           /* 18 (O) fwd, annuity etc        */
{
    static char routine[] = "Q3SwapRateCalc2";
    int         status    = FAILURE;

    double expiry, expiryVolvol, maturity, maturityVNFM, delay;
    double fwdRate;
    double zeroRateSwap, zeroRatePay, fwdRateVNFM, fwdAnnVNFM;
    int    matMonths, days, freq;

    TDateInterval interval, intervalVNFM;
    TDate         matDateVNFM;
    long          dayCntConv;
    long          dccSwapVNFM = GTO_B30_360_FIXED;
    long          dccZeroVNFM = GTO_B30_360;

    double        discStart, discPay;
    double        expToStart;

    TCurve        *indxZC = NULL;
    TCurve        *discZC = NULL;

    long           BusExpToStart    = 0;
    long           BusDaysDiff      = 0;
    long           BusDaysPerYr     = Q3_BUS_DAYS_YR;
    TMonthDayYear  MDYFirstDayLastYr;
    TMonthDayYear  MDYFirstDayNext2LastYr;
    TDate          FirstDayLastYr;
    TDate          FirstDayNext2LastYr;

    /* check date ordering */
    if (expiryDate < VolBaseDate){
            Q3ErrMsg("%s: option expiry date (%s) < ATM vol base date(%s)\n", 
                  routine, 
                  GtoFormatDate(expiryDate),
                  GtoFormatDate(VolBaseDate));
        goto RETURN;
    }
    
    if (expiryDate < SmileBaseDate){
            Q3ErrMsg("%s: option expiry date (%s) < Smile base date(%s)\n", 
                  routine, 
                  GtoFormatDate(expiryDate),
                  GtoFormatDate(SmileBaseDate));
        goto RETURN;
    }

    
    if (expiryDate > startDate || expiryDate > payDate) {
        Q3ErrMsg("%s: Effective date %s and payment date %s must be \
            >= expiration date %s\n", routine, GtoFormatDate(startDate),
            GtoFormatDate(payDate), GtoFormatDate(expiryDate));
        goto RETURN;
    }
    
    /* Note: as in smile library, shift payDate to 
     *       curve base date if between today and
     *       curve base date */
    payDate = MAX (payDate, CurveBaseDate);

    if (startDate < CurveBaseDate) {
        Q3ErrMsg("%s: Start date %s < Curve base date %s\n", routine,
            GtoFormatDate(startDate),GtoFormatDate(CurveBaseDate));
        goto RETURN;
    }

    /* compute rate maturity - 
     * if unadjusted matDate use 30/360 to get an integer 
     * if adjusted matDate, adjust the fraction number in 
     * fapricer.c 
     */
    if (GtoDayCountFraction(
        startDate,
        matDate,
        GTO_B30_360,
        &maturity) == FAILURE) goto RETURN;
    if (maturity < TINY) {
        Q3ErrMsg("%s: Rate maturity date %s must be after \
            rate start date %s\n", routine, GtoFormatDate(matDate),
            GtoFormatDate(startDate));
        goto RETURN;
    }

    
    /* exclude options on forward rates if (startDate-expiryDate)> 6 mos */
    if (GtoDayCountFraction(
        expiryDate,
        startDate,
        GTO_ACT_365F,
        &expToStart) == FAILURE) goto RETURN;
    if (expToStart > 0.502) {
        Q3ErrMsg("%s: Effective date %s  must be within 6 months \
            of expiration date %s\n", routine, GtoFormatDate(startDate), 
            GtoFormatDate(expiryDate));
        goto RETURN;
    }

    /* holiday file input checking*/
    if (BusVolDCC != NULL)
    {
        if (mystrcmp(BusVolDCC, FixBusDayPerYr) == 0) 
        {
            BusDaysPerYr = Q3_BUS_DAYS_YR;
        }
        else if (mystrcmp(BusVolDCC, ActBusDayPerYr) == 0)
        {
            if(GtoDateToMDY(
               expiryDate,
               &MDYFirstDayLastYr) == FAILURE) goto RETURN;
            
            MDYFirstDayLastYr.month         = 1;
            MDYFirstDayLastYr.day           = 1;
            MDYFirstDayNext2LastYr.month    = 1;
            MDYFirstDayNext2LastYr.day      = 1;
            MDYFirstDayNext2LastYr.year     = MDYFirstDayLastYr.year + 1;
            
            if(GtoMDYToDate(
                &MDYFirstDayLastYr, 
                &FirstDayLastYr) == FAILURE) goto RETURN;
            
            if(GtoMDYToDate(
                &MDYFirstDayNext2LastYr, 
                &FirstDayNext2LastYr) == FAILURE) goto RETURN;
            
            if(GtoBusinessDaysDiff(
                FirstDayLastYr, 
                FirstDayNext2LastYr, 
                HolidayVol, 
                &BusDaysPerYr) == FAILURE) goto RETURN;
        }
        else
        {
            Q3ErrMsg("%s: Invalid bus day vol DCC: %s.\n", 
                     routine, 
                     BusVolDCC);
            goto RETURN;
        }
    }

    /* expiry time */
    if (HolidayVol == NULL)
    {
        if (GtoDayCountFraction(
            VolBaseDate, 
            expiryDate,
            GTO_ACT_365F, 
            &expiry) == FAILURE) goto RETURN;
        if (GtoDayCountFraction(
            SmileBaseDate, 
            expiryDate,
            GTO_ACT_365F, 
            &expiryVolvol) == FAILURE) goto RETURN;
    }
    else
    {
        /* BUS/251F*/
        if((BusVolDCC == NULL) ||
           ((BusVolDCC != NULL)&&(mystrcmp(BusVolDCC, FixBusDayPerYr) == 0))) 
        {
            if(GtoBusinessDaysDiff(
               VolBaseDate,
               expiryDate,
               HolidayVol,
               &BusDaysDiff) == FAILURE) goto RETURN;
            if(GtoBusinessDaysDiff(
               expiryDate,
               startDate,
               HolidayVol,
               &BusExpToStart)==FAILURE) goto RETURN;
            expiry     = (double) (BusDaysDiff) / (double) (BusDaysPerYr);
            expToStart = (double) (BusExpToStart) / (double) (BusDaysPerYr);

            if(GtoBusinessDaysDiff(
               SmileBaseDate,
               expiryDate,
               HolidayVol,
               &BusDaysDiff) == FAILURE) goto RETURN;
            expiryVolvol = (double) (BusDaysDiff) / (double) (BusDaysPerYr);
        }
        else  /*BUS/BUS*/
        {
            if(BusDayFrac(
               VolBaseDate,
               expiryDate,
               BusDaysPerYr,
               HolidayVol,
               &expiry) == FAILURE) goto RETURN;
            if(BusDayFrac(
               expiryDate,
               startDate,
               BusDaysPerYr,
               HolidayVol,
               &expToStart) == FAILURE) goto RETURN;
            if(BusDayFrac(
               SmileBaseDate,
               expiryDate,
               BusDaysPerYr,
               HolidayVol,
               &expiryVolvol) == FAILURE) goto RETURN;
        }
        
    }

    /* compute expiry*/
    /* rate maturity in months */
    if (GtoMakeDateInterval(1, 'M', &interval) == FAILURE ||
        GtoCountDates(
            startDate,
            matDate, 
            &interval,
            &matMonths,
            &days) == FAILURE) goto RETURN;
    
    /* stub is allowed for Xm MM rates and swap rates */
    if (matMonths != 0)
    {
        /* MM rates: only 1,3,6,12 month maturities allowed */
        /* coupon rates: frequency must be 1,2,4,12         */

        if(rateFreq == 0)
        {
            if( (matMonths == 1 || matMonths == 3 || 
                 matMonths == 6 || matMonths == 12) && 
                 (days <= 10))
            {
                freq = 12 / (int) matMonths;
            }
            else if((matMonths == 2 || matMonths == 5 ||
                     matMonths == 11) &&
                     (days >=21)) 
            {
                matMonths++;
                freq = 12 / (int) matMonths;
            }
            else
            {                
                Q3ErrMsg("%s: Requested maturity %dm, %dd.\n"
                         "Only 1m, 3m, 6m, 12m money market rates allowed.\n",
                         routine,matMonths, days);
                goto RETURN;
            }

        }
        else 
        {
            if (rateFreq != 1 && rateFreq != 2 &&
                rateFreq != 4 && rateFreq != 12) {
                Q3ErrMsg("%s: Requested rate frequency %d.\n Only 1,2,4,12 "
                    "allowed.\n", routine, rateFreq);
                    goto RETURN;
            }
            if(matMonths < 11 || (matMonths == 11 && days < 21)){
                Q3ErrMsg("%s: The tenor %dm %dd is less than a year, \
                    rate frequence should be 0. \n", 
                    routine, matMonths, days);
                    goto RETURN;
            }
            if ( !((matMonths % (12 / (int)rateFreq) == 0 && days <= 10) ||
                   (matMonths % (12 /(int) rateFreq) == 
                   (12/(int) rateFreq - 1) && days >= 21)))
            {
                Q3ErrMsg("%s: The tenor %dm and %dd does not match to the \
                    rate frequency %d.\n ",
                    routine, matMonths, days, rateFreq);
                goto RETURN;
            }
            
            /* adjust the number of months to an interger number of 
             * coupon payment intervals*/
            if(days >=21) matMonths++;

            freq = (int) rateFreq;
        }
    }
    else /* 1D, 1W, 2W rates */
    {
        if (rateFreq != 0){
            Q3ErrMsg("%s: The tenor is less than a month.  \
                The rate frequence shold be 0 \n", routine);
            goto RETURN;
        }
        if (days >=1 && days <= 4) freq = 365;
        else if (days >=5 && days <= 11) freq = 52;
        else if (days >=12 && days <= 20) freq = 26;
        else if (days >=21) {
            freq = 12;
            matMonths = 1;
        }
        else{
            Q3ErrMsg("%s: The tenor cannot be zero.\n", routine);
            goto RETURN;
        }
        /* Be consistent to the freq definition of 1D, 1W, 2W*/
        dccSwapVNFM = GTO_ACT_365F; 
        dccZeroVNFM = GTO_ACT_365F;
    }
   
    /* convert zero curves */ 
    if ((discZC = GtoMakeTCurve(
        CurveBaseDate,
        discDates,
        discRates,
        numDiscPts,
        1,
        GTO_ACT_365F)) == NULL) goto RETURN;

    if ((indxZC = GtoMakeTCurve(
        CurveBaseDate,
        indxDates,
        indxRates,
        numIndxPts,
        1,
        GTO_ACT_365F)) == NULL) goto RETURN;

    if (GtoStringToDayCountConv(
        rateDCC,
        &dayCntConv) == FAILURE) goto RETURN;

    /* par rate */
    if (rateFreq != 0)
    {
        /* compute forward rate */
        if (GtoFreqAndTypeToInterval(
            freq,
            'N',
            &interval) == FAILURE) goto RETURN;

        if (GtoZerosToCouponsPoint(
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            &interval,
            matDate,
            dayCntConv,
            GTO_STUB_SIMPLE,
            FALSE,
            &fwdRate) == FAILURE) goto RETURN; 
    }
    else
    {
        if (GtoZerosToSimplePoint (
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            matDate,
            dayCntConv,
            &fwdRate) == FAILURE) goto RETURN;
    }

    
    /* 30/360 par swap rate/ MM rate for Vnfm calculation */
    if (rateFreq != 0)
    {
        /* par swap rate, tenor >= 1yr*/
        /*Addjust rate mature date for Vnfm calculation */
        if (GtoFreqAndTypeToInterval(
            freq,
            'N',
            &intervalVNFM) == FAILURE) goto RETURN;

        if(GtoDateAddMonths(
            startDate,        
            matMonths, 
            0,
            1,  
            &matDateVNFM) == FAILURE) goto RETURN;
        
        if (GtoZerosToCouponsPoint(
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            &intervalVNFM,
            matDateVNFM,
            dccSwapVNFM,
            GTO_STUB_SIMPLE,
            FALSE,
            &fwdRateVNFM) == FAILURE) goto RETURN;
    }
    else if (freq <= 12)
    {
        /* MM rate, tenor <= 1yr.*/
        
        /*Addjust rate mature date for Vnfm calculation */
        if(GtoDateAddMonths(
            startDate,        
            matMonths, 
            0,
            1,  
            &matDateVNFM) == FAILURE) goto RETURN;
        
        if (GtoZerosToSimplePoint (
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            matDateVNFM,
            dccSwapVNFM,
            &fwdRateVNFM) == FAILURE) goto RETURN;
    }
    else
    {
        /* MM rate, tenor = 1D, 1W, 2W */
        /* Addjust rate maturity date for Vnfm calculation, */
        switch(freq){
        case 365:
            if (GtoMakeDateInterval(
                1, 
                'D', 
                &intervalVNFM) == FAILURE) goto RETURN;
            break;
        case 52:
            if (GtoMakeDateInterval(
                7, 
                'D', 
                &intervalVNFM) == FAILURE) goto RETURN;
            break;
        case 26:
            if (GtoMakeDateInterval(
                14, 
                'D', 
                &intervalVNFM) == FAILURE) goto RETURN;
            break;
        default: 
            Q3ErrMsg("%s: Error: Invalid frequency%d.\n ",
                    routine);
        }

        if( GtoDtFwdAny(
            startDate,
            &intervalVNFM,
            &matDateVNFM) == FAILURE) goto RETURN;

        if (GtoZerosToSimplePoint (
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            matDateVNFM,
            dccSwapVNFM,
            &fwdRateVNFM) == FAILURE) goto RETURN;
    }

    /* fwd annuity and zero rate */
    if (GtoDiscountDate(
        matDateVNFM,
        indxZC,
        GTO_LINEAR_INTERP,
        &discPay) == FAILURE) goto RETURN;

    
    if (GtoDiscountDate(
        startDate,
        indxZC,
        GTO_LINEAR_INTERP,
        &discStart) == FAILURE) goto RETURN;
    
    if (GtoDayCountFraction(
        startDate,
        matDateVNFM,
        dccZeroVNFM,
        &maturityVNFM) == FAILURE) goto RETURN;
    zeroRateSwap = freq * (pow(discStart / discPay, 1./ freq /maturityVNFM) - 1.);
    fwdAnnVNFM = (1. - discPay / discStart) / fwdRateVNFM;

    /* compute delay interval */
    if (GtoDayCountFraction(
        startDate,
        payDate,
        dccZeroVNFM,
        &delay) == FAILURE) goto RETURN;

    /* delay zero rate (same compounding frequency as swap rate) */
    if (delay > TINY) {
        if (GtoDiscountDate(
            payDate,
            discZC,
            GTO_LINEAR_INTERP,
            &discPay) == FAILURE) goto RETURN;

        if (GtoDiscountDate(
            startDate,
            discZC,
            GTO_LINEAR_INTERP,
            &discStart) == FAILURE) goto RETURN;
         zeroRatePay = freq * (pow(discStart / discPay, 1./freq /delay) - 1.);
    } 
    else 
    {
        /* this rate not defined */
        zeroRatePay = -999;
    }
    
    /* write outputs */
    output[0]  = fwdRate;
    output[1]  = fwdRateVNFM;
    output[2]  = fwdAnnVNFM;
    output[3]  = zeroRateSwap;
    output[4]  = zeroRatePay;
    output[5]  = expiry;
    output[6]  = expiryVolvol;
    output[7]  = expiry + expToStart;
    output[8]  = maturity;
    output[9]  = delay;
    output[10] = freq;

    status = SUCCESS;

  RETURN:
    
    GtoFreeTCurve(discZC);
    GtoFreeTCurve(indxZC);

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }
  
    return status;
} /* Q3SwapRateCalc2 */

int Q3BasisLegCalc (
    long              CurveBaseDate,     /* 1  (I) Index curve base date   */
    long              numIndxPts,        /* 2  (I) Nb of index points      */
    long              *indxDates,        /* 3  (I) Basis index dates       */
    double            *indxRates,        /* 4  (I) Act/365F index rates    */
    long              numDiscPts,        /* 5  (I) Nb of discount points   */ 
    long              *discDates,        /* 6  (I) Discount dates          */
    double            *discRates,        /* 7  (I) Act/365F discount rates */
    long              fixLegFreq,        /* 8  (I) BMA and CMS fixed leg frequency      */
    char              *fixLegDCC,        /* 9  (I) BMA and CMS fixed Leg DCC            */
    long              rateFreq,          /* 10 (I) rate frequency          */
    char              *rateDCC,          /* 11 (I) rate day count conv     */
    long              expiryDate,        /* 12 (I) rate reset date         */
    long              startDate,         /* 13 (I) rate effective date     */
    long              matDate,           /* 14 (I) rate maturity date      */
    long              payDate,           /* 15 (I) rate maturity date      */
    char              *Holiday,          /* 16 (I) holiday file name       */
    double            *output            /* 17 (O) fwd, annuity etc        */
    )
{
    static char routine[] = "Q3BasisLegCalc2";
    int         status    = FAILURE;
 
    char          *fixedType    = "N";
    char          *floatType    = "N";
    char          *stubType     = "N";
    char          *stubAtEnd    = "B";
    char          *accBDC       = "N";
    char          *pmtBDC       = "N";
    char          *rstBDC       = "N";

    TBoolean       firstPmtFixed = FALSE; 
    double         firstPmtRate = 0.;

    TCurve        *indxZC = NULL;
    TCurve        *discZC = NULL;

    TDateInterval  fixedIvl;
    TDateInterval  floatIvl;
    long           fixedDCCA;
    long           floatDCCA;

    TBoolean       stubAtEndA;
    long           stubTypeA;

    long           accBDCA;
    long           pmtBDCA;
    long           rstBDCA;


    double         parRate, fltLegFv;
    double         delay;
    double         zeroDelay = 1.,
                   discPay = 1. , 
                   discStart = 1.;


    if (expiryDate > startDate) {
        Q3ErrMsg("%s: Effective date %s must be >= expiration date %s\n", 
            routine, 
            GtoFormatDate(startDate),
            GtoFormatDate(expiryDate));
        goto RETURN;
    }

    if (startDate < CurveBaseDate) {
        Q3ErrMsg("%s: Start date %s < Curve base date %s\n", routine,
            GtoFormatDate(startDate),GtoFormatDate(CurveBaseDate));
        goto RETURN;
    }

    /* convert zero curves */ 
    if ((discZC = GtoMakeTCurve(
        CurveBaseDate,
        discDates,
        discRates,
        numDiscPts,
        1,
        GTO_ACT_365F)) == NULL) goto RETURN;

    if ((indxZC = GtoMakeTCurve(
        CurveBaseDate,
        indxDates,
        indxRates,
        numIndxPts,
        1,
        GTO_ACT_365F)) == NULL) goto RETURN;

    /* Define BMA swap convention */
    if (GtoFreqAndTypeToInterval(
        fixLegFreq, 
        *fixedType, 
        &fixedIvl) == FAILURE) goto RETURN;

    if (GtoStringToDayCountConv(
        fixLegDCC, 
        &fixedDCCA) == FAILURE) goto RETURN;

    if (GtoFreqAndTypeToInterval(
        rateFreq, 
        *floatType, 
        &floatIvl) == FAILURE) goto RETURN;

    if (GtoStringToDayCountConv(
        rateDCC, 
        &floatDCCA) == FAILURE) goto RETURN;

    if (GtoStringToStubType(
            stubType, 
            &stubTypeA) == FAILURE ||
        GtoStubConvValid(
            routine, 
            stubTypeA) == FAILURE) {
        Q3ErrMsg("%s: Invalid stub type.\n", routine);
        goto RETURN;
    }

    if (GtoStringToStubAtEnd(
        stubAtEnd, 
        &stubAtEndA) == FAILURE) goto RETURN;

    accBDCA = toupper(accBDC[0]);
    pmtBDCA = toupper(pmtBDC[0]);
    rstBDCA = toupper(rstBDC[0]);

    /* compute par BMA swap rate */
    if (GtoSwapRate2(
        discZC,
        GTO_LINEAR_INTERP,
        startDate,
        matDate,
        &fixedIvl,
        fixedDCCA,
        1,   /* value floating */
        1.0, /*floatingFV,*/
        indxZC,
        GTO_LINEAR_INTERP,
        &floatIvl,
        floatDCCA,
        firstPmtFixed,
        firstPmtRate,
        FALSE,
        NULL,
        stubTypeA,
        stubAtEndA,
        accBDCA,
        pmtBDCA,
        rstBDCA,
        Holiday,
        &parRate) == FAILURE) goto RETURN;

    /* compute basis leg PV */
    if (GtoSwapFloatPV(
        discZC,
        GTO_LINEAR_INTERP,
        1.,
        0.,
        indxZC,
        GTO_LINEAR_INTERP,
        startDate,
        &floatIvl,
        matDate,
        floatDCCA,
        stubTypeA,
        stubAtEndA,
        FALSE,
        accBDCA,
        pmtBDCA,
        rstBDCA,
        Holiday,
        FALSE,
        0,
        startDate,
        FALSE,
        NULL,
        &fltLegFv) == FAILURE) goto RETURN;

    /* compute delay interval */
    if (GtoDayCountFraction(
        startDate,
        payDate,
        GTO_B30_360,
        &delay) == FAILURE) goto RETURN;

    /* delay zero rate (same compounding frequency as swap rate) */
    if (delay > TINY) {
        if (GtoDiscountDate(
            payDate,
            discZC,
            GTO_LINEAR_INTERP,
            &discPay) == FAILURE) goto RETURN;

        if (GtoDiscountDate(
            startDate,
            discZC,
            GTO_LINEAR_INTERP,
            &discStart) == FAILURE) goto RETURN;
         zeroDelay = discPay / discStart;
    } 

    output[0] = parRate;
    output[1] = fltLegFv / zeroDelay;   /* forward float leg */
    output[2] = fltLegFv / parRate;     /* Annuity           */
    output[3] = zeroDelay;              /* Z(start, pay)     */

    status = SUCCESS;
  RETURN:
    
    GtoFreeTCurve(discZC);
    GtoFreeTCurve(indxZC);

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }
  
    return status;
} /*Q3BasisLegCalc*/

/*f----------------------------------------------------------------------------
 * Q3SwapRateCalc3
 *
 * Swap rate, annuity, etc calculator from curve (for mid-curve products)
 */
int Q3SwapRateCalc3 (
    long              VolBaseDate,       /* 1  (I) ATM vol base date       */
    long              SmileBaseDate,     /* 2  (I) Smile base date         */
    long              CurveBaseDate,     /* 3  (I) Index curve base date   */
    long              numIndxPts,        /* 4  (I) Nb of index points      */
    long              *indxDates,        /* 5  (I) Index dates             */
    double            *indxRates,        /* 6  (I) Act/365F index rates    */
    long              numDiscPts,        /* 7  (I) Nb of discount points   */ 
    long              *discDates,        /* 8  (I) Discount dates          */
    double            *discRates,        /* 9  (I) Act/365F discount rates */
    long              rateFreq,          /* 10 (I) rate frequency          */
    char              *rateDCC,          /* 11 (I) rate day count conv     */
    char              *stubType,         /* 12 (I) stub type (S,B,D,N)     */
    char              *stubAtEnd,        /* 13 (I) (F) or (B)              */
    long              expiryDate,        /* 14 (I) rate reset date         */
    long              startDate,         /* 15 (I) rate effective date     */
    long              matDate,           /* 16 (I) rate maturity date      */
    long              payDate,           /* 17 (I) option payment date     */
    double            *output)           /* 18 (O) fwd, annuity etc        */
{
    static char routine[] = "Q3SwapRateCalc3";
    int         status    = FAILURE;

    double expiry, expiryVolvol, maturity, maturityVNFM, delay;
    double fwdRate, fwdRateVNFM, fwdAnnVNFM;
    double zeroRateSwap, zeroRatePay;
    double fwdAnn, fwdZero;

    int    matMonths, days;
    double shortTenor = -1; /* shortTenor = 1 if matMonths < 12/rateFreq */

    TDateInterval interval;
    long          dayCntConv;
    long          dccSwapVNFM = GTO_B30_360_FIXED;
    long          dccZeroVNFM = GTO_B30_360;

    double        discStart, discMat, discPay;
    double        expToStart;

    TCurve        *indxZC = NULL;
    TCurve        *discZC = NULL;

    TBoolean       stubAtEndA;
    long           stubTypeA;

    /* convert data as in ALIB function */
    if (GtoStringToStubType(
        stubType,
        &stubTypeA) == FAILURE ||
        GtoStubConvValid(
        routine,
        stubTypeA) == FAILURE)
    {
        Q3ErrMsg("%s: Invalid stub type.\n", routine);
        goto RETURN;
    }

    if (GTO_STUB_TYPE(stubTypeA) == GTO_STUB_DECOMPOUND)
    {
        Q3ErrMsg("%s: Decompoundend stubs not handled.\n", routine);
        goto RETURN;
    }

    if (GtoStringToStubAtEnd(
        stubAtEnd,
        &stubAtEndA) == FAILURE) goto RETURN;

    if (GtoStringToDayCountConv(
        rateDCC,
        &dayCntConv) == FAILURE) goto RETURN;

    /* check date ordering */
    if (expiryDate < VolBaseDate)
    {
        Q3ErrMsg("%s: option expiry date (%s) < ATM vol base date(%s)\n", 
            routine, 
            GtoFormatDate(expiryDate),
            GtoFormatDate(VolBaseDate));
        goto RETURN;
    }
    
    if (expiryDate < SmileBaseDate)
    {
        Q3ErrMsg("%s: option expiry date (%s) < Smile base date(%s)\n", 
            routine, 
            GtoFormatDate(expiryDate),
            GtoFormatDate(SmileBaseDate));
        goto RETURN;
    }
    
    if (expiryDate > startDate || expiryDate > payDate)
    {
        Q3ErrMsg("%s: Effective date %s and payment date %s must be \
            >= expiration date %s\n", routine, GtoFormatDate(startDate),
            GtoFormatDate(payDate), GtoFormatDate(expiryDate));
        goto RETURN;
    }
    
    /* Note: as in smile library, shift payDate to 
     *       curve base date if between today and
     *       curve base date */
    payDate = MAX (payDate, CurveBaseDate);

    if (startDate < CurveBaseDate)
    {
        Q3ErrMsg("%s: Start date %s < Curve base date %s\n", routine,
            GtoFormatDate(startDate),GtoFormatDate(CurveBaseDate));
        goto RETURN;
    }

    /* check rate frequency */
    if (rateFreq != 1 && rateFreq != 2 && rateFreq != 4 && rateFreq != 12)
    {
        Q3ErrMsg("%s: Requested rate frequency %d.\n Only 1,2,4,12 "
                 "allowed.\n", routine, rateFreq);
        goto RETURN;
    }

    /* compute rate maturity */ 
    if (GtoDayCountFraction(
        startDate,
        matDate,
        dayCntConv,
        &maturity) == FAILURE) goto RETURN;

    if (maturity < TINY)
    {
        Q3ErrMsg("%s: Rate maturity date %s must be after \
            rate start date %s\n", routine, GtoFormatDate(matDate),
            GtoFormatDate(startDate));
        goto RETURN;
    }

    if (GtoDayCountFraction(
        startDate,
        matDate,
        dccZeroVNFM,
        &maturityVNFM) == FAILURE) goto RETURN;

    /* compute maturity in months and days */
    if (GtoMakeDateInterval(1, 'M', &interval) == FAILURE ||
        GtoCountDates(
            startDate,
            matDate,
            &interval,
            &matMonths,
            &days) == FAILURE) goto RETURN;

    /* exclude options on forward rates if (startDate-expiryDate)> 6 mos */
    if (GtoDayCountFraction(
        expiryDate,
        startDate,
        GTO_ACT_365F,
        &expToStart) == FAILURE) goto RETURN;
    if (expToStart > 0.502)
    {
        Q3ErrMsg("%s: Effective date %s  must be within 6 months \
            of expiration date %s\n", routine, GtoFormatDate(startDate), 
            GtoFormatDate(expiryDate));
        goto RETURN;
    }

    /* expiry time */
    if (GtoDayCountFraction(
        VolBaseDate, 
        expiryDate,
        GTO_ACT_365F, 
        &expiry) == FAILURE) goto RETURN;
    if (GtoDayCountFraction(
        SmileBaseDate, 
        expiryDate,
        GTO_ACT_365F, 
        &expiryVolvol) == FAILURE) goto RETURN;

    /* convert zero curves */ 
    if ((discZC = GtoMakeTCurve(
        CurveBaseDate,
        discDates,
        discRates,
        numDiscPts,
        1,
        GTO_ACT_365F)) == NULL) goto RETURN;

    if ((indxZC = GtoMakeTCurve(
        CurveBaseDate,
        indxDates,
        indxRates,
        numIndxPts,
        1,
        GTO_ACT_365F)) == NULL) goto RETURN;

    /* par rate */
    if (GtoFreqAndTypeToInterval(
        rateFreq,
        'N',
        &interval) == FAILURE) goto RETURN;

    /* if matMonths < 12/rateFreq, we simply use GtoZerosToSimplePoint()
     * and return simple rate with ACT/365
     */

    if (matMonths < 12/rateFreq)
    {
        /* limiting cases */
        shortTenor = 1;

        if (GtoZerosToSimplePoint(
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            matDate,
            GTO_ACT_365F,
            &fwdRate) == FAILURE) goto RETURN;

        if (GtoZerosToSimplePoint(
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            matDate,
            GTO_ACT_365F,
            &fwdRateVNFM) == FAILURE) goto RETURN;
    }
    else
    {
        if (GtoZerosToCouponsPoint(
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            &interval,
            matDate,
            dayCntConv,
            stubTypeA,
            stubAtEndA,
            &fwdRate) == FAILURE) goto RETURN; 

        if (GtoZerosToCouponsPoint(
            indxZC,
            GTO_LINEAR_INTERP,
            startDate,
            &interval,
            matDate,
            dccSwapVNFM,
            GTO_STUB_SIMPLE,
            FALSE,
            &fwdRateVNFM) == FAILURE) goto RETURN;
    }

    /* discount factors */
    if (GtoDiscountDate(
        startDate,
        indxZC,
        GTO_LINEAR_INTERP,
        &discStart) == FAILURE) goto RETURN;

    if (GtoDiscountDate(
        matDate,
        indxZC,
        GTO_LINEAR_INTERP,
        &discMat) == FAILURE) goto RETURN;
   
    /* zeros, annuities */
    fwdZero = discMat/discStart;
    fwdAnn = (1. - fwdZero) / fwdRate;
    fwdAnnVNFM = (1. - fwdZero)/fwdRateVNFM;
    zeroRateSwap = rateFreq*(pow(1/fwdZero, 1./rateFreq/maturityVNFM) - 1.);

    /* compute delay interval */
    if (GtoDayCountFraction(
        startDate,
        payDate,
        dccZeroVNFM,
        &delay) == FAILURE) goto RETURN;

    /* delay zero rate (same compounding frequency as swap rate) */
    if (delay > TINY) {
        if (GtoDiscountDate(
            payDate,
            discZC,
            GTO_LINEAR_INTERP,
            &discPay) == FAILURE) goto RETURN;

        if (GtoDiscountDate(
            startDate,
            discZC,
            GTO_LINEAR_INTERP,
            &discStart) == FAILURE) goto RETURN;

         zeroRatePay = rateFreq * (pow(discStart / discPay, 1./rateFreq /delay) - 1.);
         fwdZero *= discStart/discPay; /* forward to paydate */
         fwdAnn  *= discStart/discPay; /* fowward to paydate */
    } 
    else 
    {
        /* this rate not defined */
        zeroRatePay = -999;
    }
    
    /* write outputs */

    output[0]  = fwdRate;
    output[1]  = fwdRateVNFM;
    output[2]  = fwdAnnVNFM;
    output[3]  = zeroRateSwap;
    output[4]  = zeroRatePay;
    output[5]  = expiry;
    output[6]  = expiryVolvol;
    output[7]  = expiry + expToStart;
    output[8]  = maturity;
    output[9]  = delay;
    output[10] = fwdAnn;
    output[11] = fwdZero;
    output[12] = shortTenor;

    status = SUCCESS;

  RETURN:
    
    GtoFreeTCurve(discZC);
    GtoFreeTCurve(indxZC);

    if (status == FAILURE) Q3ErrMsg("%s: Failed\n", routine);
  
    return status;
} /* Q3SwapRateCalc3 */

/*f------------------------------------------------------------------------
 * a copy from the vanilla product: used for midcurve only
 */
int Q3AdjustFreqDCCATMVol_BVQ(
    long              CurveBaseDate,    /* 1  (I) Base date of rate curve */
    long              numDiscPts,       /* 2  (I) Number of discount points */
    long              *discDates,       /* 3  (I) Discount dates */
    double            *discRates,       /* 4  (I) Act/365F discount rates */
    long              startDate,        /* 5  (I) Swap start date */
    long              matDate,          /* 6  (I) Swap maturity date */
    long              fixedFreq1,       /* 7  (I) non standard freq */
    long              fixedFreq2,       /* 8  (I) standard freq */
    char              *fixedType,       /* 9  (I) Type of payment(N,E,I,J) */
    char              *fixedDCC1,       /* 10 (I) non standard DCC */
    char              *fixedDCC2,       /* 11 (I) standard DCC */
    char              *stubType,        /* 12 (I) Stub type(S,B,D,N)*/
    char              *stubAtEnd,       /* 13 (I) (F)ront or (B)ack  */
    double            *output)          /* 14 (O) */
{
    static char routine[] = "Q3AdjustFreqDCCATMVol_BVP";
    int status = FAILURE;
    TCurve *discZC = NULL;
    long fixedDCCA1, fixedDCCA2;
    TBoolean stubAtEndA;
    long stubTypeA;
    TDateInterval *interval1 = NULL, 
        *interval2 = NULL;
    TDate oneyearlater;
    
    double d1, d2, f1, f2, y1, y2; /* refer to the doc */

    if (GtoStringToStubType(stubType, &stubTypeA) == FAILURE
            || GtoStubConvValid(routine, stubTypeA) == FAILURE) {
        Q3ErrMsg("%s: Invalid stub type.\n", routine);
        goto RETURN;
    }
    if (GTO_STUB_TYPE(stubTypeA) == GTO_STUB_DECOMPOUND) {
        Q3ErrMsg("%s: Decompoundend stubs not handled.\n", routine);
        goto RETURN;
    }

    if (startDate >= matDate) {
        Q3ErrMsg("%s: start date (%s) >= maturity date (%s)\n",
                routine, GtoFormatDate(startDate),
                GtoFormatDate(matDate));
        goto RETURN;
    }

    /* convert zero curves */
    if ((discZC = GtoMakeTCurve(
                    CurveBaseDate,
                    discDates,
                    discRates,
                    numDiscPts,
                    1,
                    GTO_ACT_365F)) == NULL) goto RETURN;
    
    /* consider the special case when freq == 13 */
    interval1 = GtoDateIntervalNew(28, "D");
    if (fixedFreq1 != 13 && GtoFreqAndTypeToInterval(
                fixedFreq1,
                *fixedType,
                interval1) == FAILURE) goto RETURN;

    /* consider the special case when freq == 13 */
    interval2 = GtoDateIntervalNew(28, "D");
    if (fixedFreq2 != 13 && GtoFreqAndTypeToInterval(
                fixedFreq2,
                *fixedType,
                interval2) == FAILURE) goto RETURN;

    if (GtoStringToDayCountConv(
                fixedDCC1,
                &fixedDCCA1) == FAILURE) goto RETURN;

    if (GtoStringToDayCountConv(
                fixedDCC2,
                &fixedDCCA2) == FAILURE) goto RETURN;

    if (GtoStringToStubAtEnd(
                stubAtEnd,
                &stubAtEndA) == FAILURE) goto RETURN;

    /* calculation y1, y2*/
    if (GtoZerosToCouponsPoint(
                discZC,
                GTO_LINEAR_INTERP,
                startDate,
                interval1,
                matDate,
                fixedDCCA1,
                stubTypeA,
                stubAtEndA,
                &y1) == FAILURE)
        goto RETURN;

    if (GtoZerosToCouponsPoint(
                discZC,
                GTO_LINEAR_INTERP,
                startDate,
                interval2,
                matDate,
                fixedDCCA2,
                stubTypeA,
                stubAtEndA,
                &y2) == FAILURE)
        goto RETURN;

    /* f1, f2 */
    f1 = (double) fixedFreq1;
    f2 = (double) fixedFreq2;

    /* calculate d1, d2 */
    if (GtoMakeDateInterval(
                1, 'A', interval1) == FAILURE)
        goto RETURN;
    if (GtoDtFwdAny(
                startDate,
                interval1,
                &oneyearlater) == FAILURE)
        goto RETURN;

    if (GtoDayCountFraction(
                startDate,
                oneyearlater,
                fixedDCCA1,
                &d1) == FAILURE)
        goto RETURN;
    d1 /= f1;

    if (GtoDayCountFraction(
                startDate,
                oneyearlater,
                fixedDCCA2,
                &d2) == FAILURE)
                goto RETURN;

    d2 /= f2;
    
    *output = d2*f2*(1.+d1*y1)*y2/(d1*f1*(1.+d2*y2)*y1);

    status = SUCCESS;
RETURN:
    GtoFreeTCurve(discZC);
    GtoDateIntervalDelete(interval1);
    GtoDateIntervalDelete(interval2);

    if (status == FAILURE)
        Q3ErrMsg("%s: Failed.\n", routine);

    return status;

} /* Q3AdjustFreqDCCATMVol_BVQ */
  
/* ------------------------------------------------------------
 * BusDayFrac
 */
static int BusDayFrac(
           TDate    startDate,
           TDate    expiryDate,
           long     BusDaysPerYr,
           char     *HolidayFile,
           double   *BusDaysDisc)
{
    int            status = FAILURE;
    long           BusDaysDiff = 0;
    TMonthDayYear  MDYstartDate;
    TMonthDayYear  MDYexpiryDate;
    TMonthDayYear  MDYtmpDate;
    TDate          tmpDate;
    
    if(GtoDateToMDY(
       expiryDate,
       &MDYexpiryDate) == FAILURE) goto RETURN;
    if(GtoDateToMDY(
       startDate,
       &MDYstartDate) == FAILURE) goto RETURN;

    MDYtmpDate.month  = MDYstartDate.month;
    MDYtmpDate.day    = MDYstartDate.day;
    MDYtmpDate.year   = MDYexpiryDate.year;
    
    if(GtoMDYToDate(
       &MDYtmpDate, 
       &tmpDate) == FAILURE) goto RETURN;

    if(GtoBusinessDaysDiff(
       tmpDate,
       expiryDate,
       HolidayFile,
       &BusDaysDiff)==FAILURE) goto RETURN;

    if (BusDaysPerYr == 0)
    {
        *BusDaysDisc = 0.;
    }
    else
    {
        if (tmpDate == expiryDate)
        {
            *BusDaysDisc = (double) (MDYexpiryDate.year - MDYstartDate.year);
        }
        else
        {
            *BusDaysDisc = (double) (MDYexpiryDate.year - MDYstartDate.year) + 
                          (double) (BusDaysDiff) / (double) (BusDaysPerYr);
        }
    }
    status = SUCCESS;

  RETURN:
    return(status);
}
