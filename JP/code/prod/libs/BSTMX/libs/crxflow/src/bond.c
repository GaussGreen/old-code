/***************************************************************
 * risky bond pricer
 ***************************************************************/

#include "bond.h"

/*********************************************************************************
 *    risky bond pricer
 *    
 ********************************************************************************/
int RiskyBondPV_O(
    double              *result,          /* (O) bond price or duration         */
    TDate               valDate,          /* (I) value date                     */
    TDate               issueDate,        /* (I) issue date                     */
    TDate               maturityDate,     /* (I) maturity date                  */
    TDateInterval       cpnInterval,      /* (I) coupon interval                */
    TBoolean            stubAtEnd,        /* (I) T=Stub at end; F=Stub at beg.  */
    long                stubConv,         /* (I) stub conv                      */
    TDateInterval       delay,            /* (I) default payment delay          */
    long                DCC,              /* (I) payment daycount convention    */
    double              notional,         /* (I) notional                       */
    double              coupon,           /* (I) coupon rate                    */
    double              recovery,         /* (I) recovery                       */
    KAccrualConv        accrualConv,      /* (I) accrual conv                   */
    KStubType           priceConv,        /* (I) price conv                     */
    TCurve              *discCurve,       /* (I) ir curve                       */
    TCurve              *spdCurve,        /* (I) cr curve                       */
    char                choice)           /* (I) 'P': Price, 'D': Duration      */
{
    int                 status   = FAILURE;
    static char         routine[] = "RiskyBondPV_O";

    TDateList           *cpnDateList = NULL;
    KFeeLeg_D           *feeLeg      = NULL;
    KProtLeg_D          *protLeg     = NULL;
    int                 i;
    double              feeLegPrice;
    double              protLegPrice;
    double              bondPriceL;
    double              bondYield;
    double              bondDuration;
    double              freq;
    double              dtmp;
    double              rlDf, ryDf;
    int                 NbCF;
    TDate               *feeAccStDates  = NULL, 
                        *feeAccEndDates = NULL, 
                        *feePayDates    = NULL;
    double              *feeNotionals   = NULL, 
                        *feeCoupons     = NULL;
    TDate               *protDates      = NULL;
    double              *protNotionals  = NULL, 
                        *protRecoveries = NULL;
    TDateInterval       DealFrequency;
    long                DrDates[200];
    
    GtoMakeDateInterval(1,'W',&DealFrequency);
    GtoDateIntervalToFreq(&cpnInterval,&freq);
    GtoDiscountDate(valDate, discCurve, INTERP_METHOD, &rlDf);
    GtoDiscountDate(valDate, spdCurve, INTERP_METHOD, &ryDf);

    if(valDate == maturityDate){
        *result = 0.0;
        status = SUCCESS;
        goto RETURN;
    }
    
    /* set up cpn date list */
    cpnDateList = GtoNewDateList(
        issueDate,
        maturityDate,
        &cpnInterval,        
        stubAtEnd);

    if(cpnDateList == NULL) {
        goto RETURN;
    }

    /* set up fee leg */
    NbCF           = cpnDateList->fNumItems - 1;
    
    feeAccStDates     = (TDate*) malloc(sizeof(TDate)*NbCF);
    feeAccEndDates    = (TDate*) malloc(sizeof(TDate)*NbCF);
    feePayDates       = (TDate*) malloc(sizeof(TDate)*NbCF);

    feeNotionals      = (double*) malloc(sizeof(double)*NbCF);
    feeCoupons        = (double*) malloc(sizeof(double)*NbCF);

    if((feeAccStDates  == NULL) ||
       (feeAccEndDates == NULL) ||
       (feePayDates    == NULL) ||
       (feeNotionals   == NULL) ||
       (feeCoupons     == NULL))
    {
        DR_Error("%s: faild to allocate memory", routine);
        goto RETURN;
    }
    

    for(i = 0;i < NbCF;i++)
    {
        feeNotionals[i]    = 1.0;
        feeCoupons[i]      = coupon;
        feePayDates[i]     = cpnDateList->fArray[i+1];
        feeAccStDates[i]   = cpnDateList->fArray[i];
        feeAccEndDates[i]  = cpnDateList->fArray[i+1];
        DrDates[i] = CrxTDate2DrDate(cpnDateList->fArray[i]);
    }
    

    if(!(feeLeg = RiskyFeeCreate(
           NbCF,
           feeAccStDates,
           feeAccEndDates,
           feePayDates,
           feeNotionals,
           feeCoupons,
           DCC,
           accrualConv,
           DealFrequency)))
    {
        DR_Error("Error creating fee leg.");
        goto RETURN;
    }
    
    /* pricing feeleg */
    if(RiskyFeePV_O(&feeLegPrice,
                    valDate,
                    valDate,
                    feeLeg,
                    priceConv,
                    discCurve,
                    spdCurve) == FAILURE)
    {
        DR_Error("Error calculating fee leg price.");
        goto RETURN;
    }

    /* add principal payment */
    feeLegPrice += RiskyDiscountFactor(valDate,maturityDate,discCurve,spdCurve);
    
    /* set up protection leg */
    protDates       =  (TDate*) malloc(sizeof(TDate)*1);
    protNotionals   =  (double *) malloc(sizeof(double)*1);
    protRecoveries  =  (double *) malloc(sizeof(double)*1);
    
    if((protDates      == NULL) ||
       (protNotionals  == NULL) ||
       (protRecoveries == NULL) )
    {
        DR_Error("%s: faild to allocate memory", routine);
        goto RETURN;
    }

    protDates[0]      = maturityDate;
    protNotionals[0]  = 1.0;
    protRecoveries[0] = 0.0;
    
    if(!(protLeg = ProtectionCreate(
           issueDate,
           1,
           protDates,
           protNotionals,
           protRecoveries,
           PAY_DEF,     /* Warning! hard-coded pay-def        */
           delay,
           DealFrequency)))
    {
        DR_Error("Error creating protection leg.\n");
        goto RETURN;
    }
           

    /* pricing protection leg */
    if(ProtectionPV_O(
           &protLegPrice,
           valDate,
           valDate,
           protLeg,
           discCurve,
           spdCurve) == FAILURE)
    {
        DR_Error("Error pricing protection leg.\n");
        goto RETURN;
    }

    
    bondPriceL = feeLegPrice + recovery * protLegPrice;

    if(choice == 'D'|| choice == 'Y'){
        GtoDayCountFraction(valDate, maturityDate,DCC, &dtmp);
        if(GtoBondYieldToMaturity(
               coupon,
               bondPriceL,
               (long)freq, /* Warning! */
               dtmp,
               0.1,
               stubConv,
               &bondYield) == FAILURE)
        {
            goto RETURN;
        }
        
        
        if(GtoBondModDuration(
               coupon,
               bondYield,
               (long)freq, /* Warning! */
               dtmp,
               stubConv,
               &bondDuration) == FAILURE)
        {
            goto RETURN;
        }
    }

    
    bondPriceL *= notional;

    switch(choice)
    {
    case 'P':
        *result = bondPriceL;
        break;
    case 'D':
        *result = bondDuration;
        break;
    case 'Y':
        *result = bondYield;
        break;
    default:
        DR_Error("Unknown output choice in risky bond pricer.\n");
        goto RETURN;
    }

    status = SUCCESS;
RETURN:

    if(cpnDateList    != NULL) GtoFreeDateList(cpnDateList);
    
    if(feeAccStDates  != NULL) free(feeAccStDates);
    if(feeAccEndDates != NULL) free(feeAccEndDates);
    if(feePayDates    != NULL) free(feePayDates);
    if(feeNotionals   != NULL) free(feeNotionals);
    if(feeCoupons     != NULL) free(feeCoupons);

    if(protDates      != NULL) free(protDates);
    if(protNotionals  != NULL) free(protNotionals);
    if(protRecoveries != NULL) free(protRecoveries);

    CrxFeeLegFree(feeLeg);
    CrxProtectionFree(protLeg);

    if (status != SUCCESS)
        DR_Error("%s: failed", routine);

    return status;
}
