/***************************************************************
 * Module:      cranalytics
 * File:        crcrv.c
 * Description: curve analytics
 ***************************************************************/

#include <math.h>

#include <common/include/drutils.h>
#include "crxutil.h"
#include "crcrv.h"

#include <alib/datelist.h>

#define MAX_COUPONS 100
#define MAX_CDS  25

/* Tolerance for Newton-Rapson */
#define TOL 1e-8

/* Minimum abs of first derivative */
#define SMALL 1e-14

#define MAX_ITER 20

/* DateFloor
 *
 * Computes result such that  baseDate + n*frequency = result <= date < baseDate + (n+1) * frequency
 *
 */
static int DateFloorSafe(TDate *result, long *nperiods, TDate baseDate, TDate date, TDateInterval frequency, int backup);

static int DateFloor(TDate *result, long *nperiods, TDate baseDate, TDate date, TDateInterval frequency)
{
    return DateFloorSafe(result, nperiods, baseDate, date, frequency, 0);
}

int DateFloorSafe(TDate *result, long *nperiods,  TDate baseDate, TDate date, TDateInterval frequency, int backup)
{

    char routine[] = "DateFloorSafe";
    long status = FAILURE;
    
    double dayFreq;
    
    TDate tmpDate, date1;

    if (GtoDateFromDateAndOffset( 
            date, 
            &frequency, 
            10, 
            &tmpDate) != SUCCESS)
    {
        DR_Error("%s: failed to advance date.", routine);
        goto RETURN;
    }

    dayFreq = (tmpDate-date) / 10.0;
    *nperiods =  (int) floor( (date  - baseDate)/dayFreq ) - backup;
    GtoDateFromDateAndOffset(baseDate, &frequency, *nperiods, &date1);

    if (date1 < date)
    {
        tmpDate = date1;
        do {
            nperiods[0]++;
            date1 = tmpDate;
            GtoDateFromDateAndOffset(baseDate, &frequency, *nperiods, &tmpDate);
        } while(tmpDate <= date);

        nperiods[0]--;
    }
    else if (date1 > date)
    {
        DateFloorSafe(&date1, nperiods, baseDate, date, frequency, backup + 2);
    }

    *result = date1;

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);
    return status;
}

/** Internal function */
int FowardToValueDate(
    double        *fwdFactor,      /**<(O) forward factor                  */
	TDate         valDate,         /**<(I) valuation date                  */
    const TCurve  *discZC,         /**<(I) Interest rate zero curve        */
    const TCurve  *sprdZC          /**<(I) Clean spread zero curve         */
    )
{

    double temp;
    long status = FAILURE;
    char routine[] = "FowardToValueDate";
      
    if (GtoDiscountDate(
        valDate, 
        (TCurve*)discZC, 
        INTERP_METHOD, 
        &temp) != SUCCESS)
    {
        DR_Error("%s: failed to compute discount factor.");
        goto RETURN;
    }

    *fwdFactor = 1./temp;

    if (sprdZC != NULL)
    {
        if (GtoDiscountDate(
                valDate, 
                (TCurve*)sprdZC, 
                INTERP_METHOD, 
                &temp) != SUCCESS)
        {
            DR_Error("%s: failed to compute discount factor.");
            goto RETURN;
        }

        *fwdFactor /= temp;
    }

    status = SUCCESS;
RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);
    return status;
}

KProtLeg_D* ProtectionCreate(
    TDate         stDate,        /* (I)   Start date of protection           */
    int           nbPrd,         /* (I)   Number of protection periods       */
    const TDate   *dates,        /* (I)   End dates for protection schedules */
    const double  *notionals,    /* (I)   Notional for each interval         */
    const double  *recoveries,   /* (I)   Recovery for each interval         */
    KProtPayConv  payType,       /* (I)   Timing of protection payment       */
    TDateInterval delay,         /* (I)   Payment delay (as offset from def) */
    TDateInterval frequency)     /* (I)   Integration interval               */
{

    char routine[] = "ProtectionCreate";
    int status = FAILURE;

    int         i; // loop counting var.
    
    KProtLeg_D  *protLeg = NULL; // the object that is being created

    /* ALLOCATE MEMORY FOR OVERALL STRUCTURE */
    protLeg = (KProtLeg_D *) calloc(1, sizeof(KProtLeg_D));
    if (protLeg == NULL)
    {
        DR_Error("%s: failed to allocate memory.", routine);
        goto RETURN;    
    }

    /* SET NON-POINTER MEMBERS */
    protLeg->mNbPrd     = nbPrd;        // No. periods
    protLeg->mStDate    = stDate;       // start date
    protLeg->mDelay     = delay;        // payment delay
    protLeg->mFrequency = frequency;    // integration freq.
    protLeg->mPayType   = payType;      // pay on default or at mty     

    /* ALLOCATE MEMORY FOR POINTER MEMBERS */
    if (!(protLeg->mDates = (TDate*) malloc(sizeof(TDate)*nbPrd))) 
    {
        DR_Error("%s failed: Unable to allocate memory for dates.\n");
        goto RETURN;
    }
    if (!(protLeg->mNotionals    = (double*) malloc(sizeof(double)*nbPrd))) 
    {
        DR_Error("%s failed: Unable to allocate memory for notionals.\n");
        goto RETURN;
    }
    if (!(protLeg->mRecoveries   = (double*) malloc(sizeof(double)*nbPrd)))
    {
        DR_Error("%s failed: Unable to allocate memory for recoveries.\n");
        goto RETURN;
    }

    /* FILL ARRAYS WITH CORRECT VALUES */
    for (i = 0; i < nbPrd; i++)
    {
       protLeg->mDates[i]       = dates[i];
       protLeg->mNotionals[i]   = notionals[i];
       protLeg->mRecoveries[i]  = recoveries[i];
    }
       
    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
    {
        DR_Error ("%s: failed", routine);
        return NULL;
    };

    return protLeg;

}

/**Free protection leg*/
int CrxProtectionFree(
    KProtLeg_D *protLeg)
{
    if(protLeg != NULL)
    {
        SafeFree(protLeg->mDates);
        SafeFree(protLeg->mNotionals);
        SafeFree(protLeg->mRecoveries);

        free(protLeg);
    }

    return SUCCESS;
}

/**Create a risky-fee leg */
KFeeLeg_D* RiskyFeeCreate(
    int           nbCF,
    const TDate   *accStDates,
    const TDate   *accEndDates,
    const TDate   *payDates,
    const double  *notionals,
    const double  *coupons,
    TDayCountConv dcc,
    KAccrualConv  payAccr,
	TDateInterval frequency)
{
    /* LOCAL VARIABLES */
    int status = FAILURE;
    char routine[] = "RiskyFeeCreate";
    int i;                          // loop variable
	KFeeLeg_D     *feeLeg = NULL;   // object that is being created

    /* ALLOCATE OBJECT */
    feeLeg = (KFeeLeg_D *) calloc(1, sizeof(KFeeLeg_D));
    if (feeLeg == NULL)
    {
        DR_Error("%s: failed to allocate memory for KFeeLeg_D\n", routine);
        goto RETURN;    
    }
    
    /* INITIALIZE NON-POINTER MEMBERS */
    feeLeg->mNbCF           = nbCF;         // No. cash-flows
    feeLeg->mDCC            = dcc;          // fee daycount convention
    feeLeg->mFrequency      = frequency;    // integration freq. for accruals
    feeLeg->mAccrualConv    = payAccr;      // accroed on default convention

    /* ALLOCATE POINTER MEMBERS */
    feeLeg->mAccStDates     = (TDate*) malloc(sizeof(TDate)*nbCF);
    if (!feeLeg->mAccStDates) {
        DR_Error("%s: failed to allocate memory for accrual start dates.\n", routine);
        goto RETURN;  
    }
    feeLeg->mAccEndDates    = (TDate*) malloc(sizeof(TDate)*nbCF);
    if (!feeLeg->mAccEndDates) {
        DR_Error("%s: failed to allocate memory for accrual end dates.\n", routine);
        goto RETURN;  
    }
    feeLeg->mPayDates       = (TDate*) malloc(sizeof(TDate)*nbCF);
    if (!feeLeg->mPayDates) {
        DR_Error("%s: failed to allocate memory for fee payment dates.\n", routine);
        goto RETURN;  
    }
    feeLeg->mNotionals      = (double*) malloc(sizeof(double)*nbCF);
    if (!feeLeg->mNotionals) {
        DR_Error("%s: failed to allocate memory for fee notionals.\n", routine);
        goto RETURN;
    }
    feeLeg->mCoupons        = (double*) malloc(sizeof(double)*nbCF);
    if (!feeLeg->mCoupons) {
        DR_Error("%s: failed to allocate memory for fee coupon rates.\n", routine);
        goto RETURN;
    }

    /* FILL ARRAY VALUES */
    for (i = 0; i < nbCF; i++)
    {
        feeLeg->mAccStDates[i]  = accStDates[i];
        feeLeg->mAccEndDates[i] = accEndDates[i];
        feeLeg->mPayDates[i]    = payDates[i];
        feeLeg->mNotionals[i]   = notionals[i];
        feeLeg->mCoupons[i]     = coupons[i];
    }

    status = SUCCESS;
RETURN:

    if (status != SUCCESS)
    {
        DR_Error ("%s: failed", routine);
        CrxFeeLegFree(feeLeg);

        return NULL;
    };
    return feeLeg;
}

/**Create a risky fee leg from a frequency, rather than from explicit arrays */
KFeeLeg_D* CrxFeeLegCreateFromFreq(
    TDate         sttDate,       /* (I) accrual start dates.                 */
    TDate         endDate,
    TDateInterval couponFrequency,
    KStubLocation stub,
    double        notional,
    double        coupon,
    TDayCountConv dcc,
    KAccrualConv  payAccr,
	TDateInterval frequency)
{
    long status    = FAILURE;
    char routine[] = "CrxFeeLegCreateFromFreq";

    /* LOCAL VARIABLES */
    TDateList     *dateList = NULL;
    TDate         *dates = NULL, 
                  *useDates = NULL;
    TDate         newDate;
    int           i, numDates;
    double        *notionals = NULL, 
                  *coupons   = NULL;
    TBoolean      stubFlag;
    KFeeLeg_D     *feeLeg = NULL;

    /* TRUE   =  stub at end */
    stubFlag  = (stub == SHORT_BACK) || (stub == LONG_BACK);
    /* Construct list of fee dates */
    dateList  = GtoNewDateList(sttDate, endDate, &couponFrequency, stubFlag);
    if (!dateList) {
        DR_Error("%s failed: Could not create date list.", routine);
        goto RETURN;
    }
    /* ALLOCATE SPACE FOR dates, notionals and coupons */
    dates     = (TDate *) malloc( (dateList->fNumItems+1) * sizeof(TDate) );
    if (!dates) {
        DR_Error("%s failed: Unable to allocate memory for dates.", routine);
        goto RETURN;
    }
    notionals = (double *) malloc( (dateList->fNumItems+1) * sizeof(double) );
    if (!notionals) {
        DR_Error("%s failed: Unable to allocate memory for notionals.", routine);
        goto RETURN;
    }
    coupons   = (double *) malloc( (dateList->fNumItems+1) * sizeof(double) );
    if (!coupons) {
        DR_Error("%s failed: Unable to allocate memory for coupon rates.", routine);
        goto RETURN;
    }

    /* FILL ARRAYS OF DATES, COUPONS & NOTIONALS */
    for (i= 0; i< dateList->fNumItems; i++)
    {
        dates[i]     = dateList->fArray[i];
        notionals[i] = notional;
        coupons[i]   = coupon;
    };

    numDates = dateList->fNumItems;
    useDates = dates;                   // copy so you can mess with the start and still delete OK
    
    /* ADJUST DATES FOR STUB */
    switch (stub)
    {
        /* Default date construction is OK for these stub types */
        case SHORT_BACK: 
        case SHORT_FRONT: break;

        /* Copy maturity back one index and decrement date size */
        case LONG_BACK:
            if (numDates < 2) break;
            GtoDateFromDateAndOffset(dates[numDates-1], &couponFrequency, -1, &newDate);
            if (newDate == dates[numDates-2]) break;
            dates[numDates-2] = dates[numDates-1];
            numDates--;
            break;

        /* Copy start forward one index and increment dates pointers */
        case LONG_FRONT:
            if (numDates < 3) break;
            GtoDateFromDateAndOffset(dates[0], &couponFrequency, 1, &newDate);
            if (newDate == dates[1]) break;
			dates[1] = dates[0];
            useDates++;
            numDates--;
            break;
    }

    if (!(feeLeg = RiskyFeeCreate(
            numDates-1, 
            useDates,       // implicitly assumes no payment delay
            useDates+1, 
            useDates+1, 
            notionals, 
            coupons, 
            dcc, 
            payAccr, 
            frequency)))
    {
        DR_Error("%s: failed to create risky fee leg\n", routine);
        goto RETURN;
    };
    
    status = SUCCESS;
RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    /* Free local pointer variables */
    SafeFree(notionals);
    SafeFree(coupons);
    SafeFree(dates);            // NB also frees memory for useDates
    GtoFreeDateList(dateList);

    if (status != SUCCESS)
    {
        DR_Error ("%s: failed", routine);
        CrxFeeLegFree(feeLeg);
        return NULL;
    };
	return feeLeg;
}


/* Free the fee leg
 */
int CrxFeeLegFree(
    KFeeLeg_D *feeLeg)
{
    if(feeLeg != NULL)
    {
        SafeFree(feeLeg->mAccStDates);
        SafeFree(feeLeg->mAccEndDates);
        SafeFree(feeLeg->mPayDates);
        SafeFree(feeLeg->mNotionals);
        SafeFree(feeLeg->mCoupons);
    }

    free(feeLeg);

    return SUCCESS;
}
      
/**Create a CDS object*/
KCDS_D* CDSCreate(
    int           nbFeePeriods,    /* (I) Number of fee payments           */
    const TDate   *accStDates,     /* (I) Accrual start dates for payments */
    const TDate   *accEndDates,    /* (I) Accrual end dates for payments   */
    const TDate   *payDates,       /* (I) Payment dates for fee payments   */
    const double  *notionals,      /* (I) Notionals for each period        */
    const double  *coupons,        /* (I) Coupon rates for each period     */   
	const double  *recoveries,     /* (I) recovery rates for coupon period */
    TDayCountConv dcc,             /* (I) */
    KAccrualConv  payAccr,         /* (I) */
    TDate         protStDate,      /* (I) */
    KProtPayConv  payType,         /* (I) */
	TDateInterval delay,           /* (I) */
	TDateInterval frequency)       /* (I) Contingent leg integration frequency, I think (CM) */
{
    
    char *routine = "CDSCreate";
    int  status   = FAILURE;

    KCDS_D *cds   = NULL; // object to create

    /** ALLOCATE MEMORY FOR MAIN OBJECT */
    cds = (KCDS_D*) calloc(1, sizeof(KCDS_D));
    if (cds == NULL)
    {
        DR_Error("%s: failed to allocate memory for KCDS_D.\n", routine);
        goto RETURN;
    }

    /* CREATE RISKY FEE LEG */
    if ( !(cds->feeLeg = RiskyFeeCreate(
        nbFeePeriods,
        accStDates,
        accEndDates, 
        payDates, 
        notionals, 
        coupons, 
        dcc, 
        payAccr, 
        frequency)))
    {
        DR_Error("%s: error creating fee leg", routine);
        goto RETURN;
    }

    /* CREATE PROTECTION LEG */
    if ( !(cds->protLeg = ProtectionCreate(
        protStDate, 
        nbFeePeriods, 
        accEndDates, 
        notionals, 
        recoveries, 
        payType, 
        delay, 
        frequency)))
    {
        DR_Error("%s: error creating protection leg", routine);
        goto RETURN;
    }

    status = SUCCESS;
RETURN:
    if (status != SUCCESS)
    {
        DR_Error ("%s: failed", routine);
        CDSFree(cds);
        return NULL;
    };

    return cds;
}


/** 
 *  Free CDS structure 
 */
int    CDSFree(
	KCDS_D        *cds)            /* (I) Pointer to CDS structure       */
{
    if(cds != NULL)
    {
        CrxFeeLegFree(cds->feeLeg);
        CrxProtectionFree(cds->protLeg);

        free(cds);
    }

    return SUCCESS;
}

/*********************************************************************************
 *    risky discount factor with flat forward interp
 ********************************************************************************/
double RiskyDiscountFactor(
    TDate               startDate,        /* (I) start date                     */
    TDate               endDate,          /* (I) end   date                     */
    const TCurve        *irCurve,         /* (I) ir curve                       */
    const TCurve        *crCurve)         /* (I) cr curve                       */
{
    /* Local vars */
    double irTmp1, irTmp2, crTmp1, crTmp2;

    /* IR Discount */
    GtoInterpPV(startDate,(TCurve*)irCurve,INTERP_METHOD,&irTmp1);
    GtoInterpPV(endDate,(TCurve*)irCurve,INTERP_METHOD,&irTmp2);

    /* Survival probability */
    GtoInterpPV(startDate,(TCurve*)crCurve,INTERP_METHOD,&crTmp1);
    GtoInterpPV(endDate,(TCurve*)crCurve,INTERP_METHOD,&crTmp2);

    return irTmp2 * crTmp2 / irTmp1 / crTmp1;
}

/* One-period protection leg value.  Computed by assuming flat forward
 * for credit AND interest rate within the period.
 */
#define IntgrlProt(ryDf1,ryDf,rlDf1,rlDf)  ((IS_EQUAL((ryDf),(ryDf1)))?(0e0): \
                                           ((ryDf1*rlDf1-ryDf*rlDf)           \
                                            *log(ryDf1/ryDf)/(log(ryDf1/ryDf) \
                                            +log(rlDf1/rlDf))))

/*-----------------------------------------------------------------------
 *
 * ProtectionSinglePV - LOCAL FUNCTION
 *
 * Computes price of protection for a period with constant recovery and notional
 *
 */
int    ProtectionSinglePV(
    double        *pv,
    TDate         today,
    TDate         stDate,
    TDate         endDate,
    double        notional,
    double        recovery,
    KProtPayConv  payType,
    TDate         payDate,
    TDateInterval delay,
    TDateInterval frequency,
    const TCurve  *discZC,
    const TCurve  *sprdZC)
{    
    long status = FAILURE;
    char routine[] = "ProtectionSinglePV";
    
    /* LOCAL VARIABLES */
    TDate date1, date, datePmt, datePmt1;
    double ryDf, rlDf, ryDf1, rlDf1;
    long np1;
    double result;
    

    result = 0.0;
    
    if (payType == PAY_MAT)
    {
        if (FAILURE==GtoDiscountDate(payDate, (TCurve*)discZC, GTO_FLAT_FORWARDS,  &rlDf))
        {
            DR_Error("%s failed: unable to calculate discount to maturity.\n");
            goto RETURN;
        }
        rlDf1 = rlDf;
    }
    
    date = MAX(today, stDate);
    
    /* date = start of protection as priced 
     * We need to subtract protection from grid point to protection start 
     * date >= today >= baseDate 
     * date1  = discZC->fBaseDate + frequency * ((date - discZC->fBaseDate) 
     *           / frequency); 
     */
    
    DateFloor(&date1, &np1, discZC->fBaseDate, date, frequency);
    
    /* date1 = highest i*freq + basedate
     * date >= date1 >= baseDate 
     */
    
    if (FAILURE==GtoDiscountDate(date1, (TCurve*)sprdZC, GTO_FLAT_FORWARDS, &ryDf1)) {
        goto RETURN;
    }
    if (FAILURE==GtoDiscountDate(date,  (TCurve*)sprdZC, GTO_FLAT_FORWARDS, &ryDf)) {
        goto RETURN;
    }
    
    if (payType == PAY_DEF)
    {
        /* Start date */
        if (FAILURE==GtoDtFwdAny(date1, &delay, &datePmt1))
            goto RETURN;
        if (FAILURE==GtoDiscountDate(datePmt1, (TCurve*)discZC, GTO_FLAT_FORWARDS, &rlDf1))
            goto RETURN;
        
        /* End date */
        if (FAILURE==GtoDtFwdAny(date, &delay, &datePmt))
            goto RETURN;
        if (FAILURE==GtoDiscountDate(datePmt, (TCurve*)discZC, GTO_FLAT_FORWARDS, &rlDf))
            goto RETURN;
    }
    
    /* Risky discount, one-period forward IR rate and credit spread
     * for the calculation of protection
     */
    
    /* One period protection value */
    result -= notional * (1. - recovery) * IntgrlProt(ryDf1, ryDf, 
                                                      rlDf1, rlDf);
    
    if (FAILURE==GtoDateFromDateAndOffset(discZC->fBaseDate, &frequency, ++np1, &date))
        goto RETURN;
    
    date = MIN(date , endDate);		    
    
    while (1)
    {
        if (payType == PAY_DEF)
        {
            /* Start date only */
            if (FAILURE==GtoDtFwdAny(date, &delay, &datePmt))
                goto RETURN;
            if (FAILURE==GtoDiscountDate(datePmt, (TCurve*)discZC, GTO_FLAT_FORWARDS, &rlDf))
                goto RETURN;
        };
        
        if (FAILURE==GtoDiscountDate(date, (TCurve*)sprdZC, GTO_FLAT_FORWARDS, &ryDf))
            goto RETURN;
        
        /* Add one period protection value */
        result += notional * (1. - recovery) * IntgrlProt(ryDf1, ryDf, 
                                                          rlDf1, rlDf);
        
        if (date == endDate) break;
        
        ryDf1 = ryDf;
        rlDf1 = rlDf;
        
        if (FAILURE==GtoDateFromDateAndOffset(discZC->fBaseDate, &frequency, ++np1, &date))
            goto RETURN;
        
        date = MIN(date, endDate);		    
    }
    
    *pv = result;
    
    status = SUCCESS;
RETURN:
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}

/**---------------------------------------------------------------------------
 *
 * ProtectionPV
 *
 * Computes price of protection for a schedule of recovery rates and notionals
 * (Non-object version)
 *
 */
int    ProtectionPV(
    double        *pv,              /* (O) pv of protection leg             */
    TDate         today,            /* (I) compute price as of date         */
    TDate         valDate,          /* (I) return price as of date          */
    TDate         stDate,           /* (I) start date of the protection     */
    int           nbPrd,            /* (I) number of periods                */
    const TDate   *dates,           /* (I) date schedule                    */
    const double  *notionals,       /* (I) notional schedule.               */
    const double  *recoveries,      /* (I) recovery rate schedule           */
    KProtPayConv  payType,          /* (I) PAY_DEF or PAY_MAT.              */
    TDateInterval delay,            /* (I) payment delay for protection leg */
    TDateInterval frequency,        /* (I) Integration interval             */
    const TCurve  *discZC,          /* (I) Interest rate zero curve         */
    const TCurve  *sprdZC)          /* (I) Clean spread zero curve          */
{
	long status = FAILURE;
	char routine[] = "ProtectionPV";

    /* LOCAL VARIABLES */
	int i, j;
	double pv1, factor;

	*pv = 0.0;
	i = 0;
	
	if (today > stDate)
	{
		while (today > dates[i])
		{
			i++;
			if (i == nbPrd)
				goto RETURN;
		};
	}

	for (j = i; j < nbPrd; j++)
	{
        if (ProtectionSinglePV(
            &pv1, 
            today, 
            j?dates[j-1]:stDate, 
            dates[j], 
            notionals[j],
            recoveries[j],
            payType,
            dates[nbPrd-1],      
            delay,               
            frequency,
            discZC,
            sprdZC) != SUCCESS)
        {
            DR_Error("%s: failed to compute protection single.", routine);
            goto RETURN;
        };

        *pv += pv1;
    }

    if ( FowardToValueDate(
            &factor, 
            valDate, 
            discZC, 
            sprdZC) != SUCCESS)
    {
        DR_Error("%s: failed to compute fwd factor", routine);
        goto RETURN;
    }

    *pv *= factor;

	status = SUCCESS;
RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

	return status;
}



/*------------------------------------------------------------------------
 *
 * ProtectionPV_O
 *
 * Computes the price of a protection leg (object version) 
 *
 * Protection from x to x is zero
 * Protection from 0 to x = sum [Prob[default happens on day i], i = 1 to x   (timing factors aside)
 *
 */
int ProtectionPV_O(
    double              *pv,              /* (O) protection price                */
    TDate               today,           /* (I) compute price as of date        */
	TDate               valDate,         /* (I) return price as of date         */
    const KProtLeg_D    *protLeg,         /* (I) Protection leg object           */
    const TCurve        *discZC,          /* (I) Interest rate zero curve        */
    const TCurve        *sprdZC)          /* (I) Clean spread zero curve         */
{
    int status = FAILURE;
	char routine[] = "ProtectionPV_O";

	int i, j;
	double pv1, factor;

	*pv = 0.0;
    
    i = 0;

    /* find first interval to be priced */
	if (today > protLeg->mStDate)
	{
		while (today > protLeg->mDates[i])
		{
			i++;
			if (i == protLeg->mNbPrd)
            {
                status = SUCCESS;
				goto RETURN;
            }
		};
	}

    for (j = i; j < protLeg->mNbPrd; j++)
	{
        if (ProtectionSinglePV(&pv1, 
            today, 
            (j > 0)?protLeg->mDates[j-1]:protLeg->mStDate, 
            protLeg->mDates[j], 
            protLeg->mNotionals[j],
            protLeg->mRecoveries[j],
            protLeg->mPayType,
            protLeg->mDates[protLeg->mNbPrd-1],      
            protLeg->mDelay,               
            protLeg->mFrequency,
            discZC,
            sprdZC
            ) == FAILURE) goto RETURN;

        *pv += pv1;
	}

    if (FowardToValueDate(
        &factor, 
        valDate, 
        discZC, 
        sprdZC) != SUCCESS)
    {
        DR_Error("%s: failed to compute fwd factor", routine);
        goto RETURN;
    }

    *pv *= factor;

	status = SUCCESS;
 RETURN:
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);
	return status;

} /* ProtectionPV_O */ 


/*------------------------------------------------------------------------
 *
 * RiskyFeePV_O
 *
 * Computes the price of a fee leg (object version) 
 *
 */

int    RiskyFeePV_O(
    double       *pv,
    TDate         today,           /* (I) compute price as of date        */
	TDate         valDate,         /* (I) return price as of date         */
    const KFeeLeg_D    *feeLeg,
    KStubType     priceConv,
    const TCurve       *discZC,          /* (I) Interest rate zero curve        */
    const TCurve       *sprdZC)          /* (I) Clean spread zero curve         */
{
	int i, j;
    long np1;
	double ryDf, rlDf, ryDf1;
	double dt, dt1, ntnlFactor, factor;
	TDate date, date1;

    int status = FAILURE;
    char routine[] = "RiskyFeePV_O";
    
    *pv = 0;

    i = 0;    
    /* if today < feeLeg->mAccStDates[0], we start with the first coupon */
    if (today >= feeLeg->mAccStDates[0])
    {
        while (i < feeLeg->mNbCF)
		    if ((today <feeLeg->mAccEndDates[i]) && (today >= feeLeg->mAccStDates[i]))
			    break;
		    else
			    i++;
    }

	if (i == feeLeg->mNbCF)
    {
        status = SUCCESS;
		goto RETURN;
    }

    ntnlFactor = 1.0;

    if (today >= feeLeg->mAccStDates[i])
    {
        if (BOND == priceConv)
        {

            GtoDayCountFraction(feeLeg->mAccStDates[i], today, feeLeg->mDCC, &dt);

            *pv = -dt * feeLeg->mCoupons[i] * feeLeg->mNotionals[i];
        }

        if (SIMPLE == priceConv)
        {
            GtoDayCountFraction(today, feeLeg->mAccEndDates[i], feeLeg->mDCC, &dt);

            GtoDayCountFraction(feeLeg->mAccStDates[i], feeLeg->mAccEndDates[i], feeLeg->mDCC, &dt1);

            ntnlFactor = dt/dt1;
        }
    }

	for (j= i; j < feeLeg->mNbCF; j++) 
    {
		if(feeLeg->mAccrualConv == ACCRUAL_PAY_ALL)
		{
		    date = MAX(today, feeLeg->mAccStDates[j]);

            GtoDiscountDate(date,  (TCurve*)sprdZC, INTERP_METHOD, &ryDf1);

            DateFloor(&date1, &np1, discZC->fBaseDate, date, feeLeg->mFrequency);

            GtoDayCountFraction(feeLeg->mAccStDates[j], date, feeLeg->mDCC, &dt1);

            if (date1 < date) { /* treat stub between date and date2 before we go to frequency intervals */

               GtoDateFromDateAndOffset(discZC->fBaseDate, (TDateInterval*)&(feeLeg->mFrequency), ++np1, &date);

               date = MIN(date, feeLeg->mAccEndDates[j]);

               GtoDiscountDate(date,  (TCurve*)sprdZC, INTERP_METHOD, &ryDf);

	           GtoDiscountDate(date,  (TCurve*)discZC, INTERP_METHOD, &rlDf);

               GtoDayCountFraction(feeLeg->mAccStDates[j], date, feeLeg->mDCC, &dt);

               *pv += 0.5 * (dt1+dt) * ntnlFactor * feeLeg->mNotionals[j] * feeLeg->mCoupons[j] * rlDf * (ryDf1-ryDf);

               ryDf1 = ryDf;

               dt1 = dt;

            }

            GtoDateFromDateAndOffset(discZC->fBaseDate, (TDateInterval*)&(feeLeg->mFrequency), ++np1, &date);
        
            date = MIN(date , feeLeg->mAccEndDates[j]);
 
		    while (date <= feeLeg->mAccEndDates[j])
		    {

                GtoDiscountDate(date, (TCurve*)discZC, INTERP_METHOD, &rlDf);
    
                GtoDiscountDate(date, (TCurve*)sprdZC, INTERP_METHOD, &ryDf);
                
                GtoDayCountFraction(feeLeg->mAccStDates[j], date, feeLeg->mDCC, &dt);
						        
			    *pv += 0.5 * (dt1+ dt) * ntnlFactor * feeLeg->mNotionals[j] * feeLeg->mCoupons[j] * rlDf * (ryDf1-ryDf);

                if (date == feeLeg->mAccEndDates[j]) break;

			    ryDf1 = ryDf;

                dt1 = dt;

                GtoDateFromDateAndOffset(discZC->fBaseDate, (TDateInterval*)&(feeLeg->mFrequency), ++np1, &date);
        
		        date = MIN(date , feeLeg->mAccEndDates[j]);		    
            }

        }  /* if accrualConv == PAY_ALL */
		
	    GtoDiscountDate(feeLeg->mAccEndDates[j], (TCurve*)discZC, INTERP_METHOD, &rlDf);
  
        GtoDiscountDate(feeLeg->mAccEndDates[j], (TCurve*)sprdZC, INTERP_METHOD, &ryDf);

        GtoDayCountFraction(feeLeg->mAccStDates[j], feeLeg->mAccEndDates[j], feeLeg->mDCC, &dt);

		*pv += feeLeg->mNotionals[j]*ntnlFactor*feeLeg->mCoupons[j]*rlDf*ryDf*dt;

        ntnlFactor = 1.0;

	}

    if ( FowardToValueDate(
            &factor, 
            valDate, 
            discZC, 
            sprdZC) != SUCCESS)
    {
        DR_Error("%s: failed to compute fwd factor", routine);
        goto RETURN;
    }

    *pv *= factor;

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}

/* Assign pointers to a fee leg, rather than deep copies for speed - local function only */
int FeeLegInit(
    KFeeLeg_D     *feeLeg,
    int           nbCF,              /* (I) nb of fee payment periods.           */
    const TDate         *accStDates,       /* (I) accrual start dates.                 */
    const TDate         *accEndDates,
    const TDate         *payDates,
    const double        *notionals,
    const double        *coupons,
    TDayCountConv dcc,
    KAccrualConv  payAccr,
	TDateInterval frequency)
{
    feeLeg->mAccStDates  = (TDate*)accStDates;
    feeLeg->mAccEndDates = (TDate*)accEndDates;
    feeLeg->mAccrualConv = payAccr;
    feeLeg->mCoupons     = (double*)coupons;
    feeLeg->mDCC         = dcc;
    feeLeg->mNbCF        = nbCF;
    feeLeg->mPayDates    = (TDate*)payDates;
    feeLeg->mNotionals   = (double*)notionals;
    feeLeg->mFrequency   = frequency;

    return SUCCESS;
}

/** 
 * Compute the value of a risky fee leg 
 * (Non-object version)
 *
 */
int    RiskyFeePV(
    double        *pv,              /* (O) price of fee leg                */
    TDate         today,            /* (I) compute price as of date        */
    TDate         valDate,          /* (I) return price as of date         */
    int           nbCF,             /* (I) nb of fee payment periods.      */
    const TDate         *accStDates,      /* (I) accrual start dates.            */
    const TDate         *accEndDates,
    const TDate         *payDates,
    const double        *notionals,
    const double        *coupons,
    TDayCountConv dcc,
    KAccrualConv  payAccr,
    KStubType     priceConv,
    TDateInterval frequency,        /* Integration interval? */
    const TCurve		  *discZC,          /* (I) Interest rate zero curve        */
    const TCurve		  *sprdZC)          /* (I) Clean spread zero curve         */
{
    int status = FAILURE;
    char routine[] = "RiskyFeePV";

    KFeeLeg_D feeLeg;
         
    if (FeeLegInit(
            &feeLeg,
            nbCF, 
            accStDates,                
            accEndDates,
            payDates,
            notionals,
            coupons,
            dcc,
            payAccr,
            frequency) != SUCCESS)
    {
        DR_Error("%s: failed to create fee leg.", routine);
        goto RETURN;
    }
    
    if (RiskyFeePV_O(
            pv, 
            today, 
            valDate,
            &feeLeg, 
            priceConv, 
            discZC, 
            sprdZC) != SUCCESS)
    {
        DR_Error("%s: failed to compute risky pv.", routine);
        goto RETURN;
    }
    
    status = SUCCESS;
    
RETURN:
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);
    return status;
    
} 




/**
 * Compute the adjustment factor to par spread for pay accrue upon default 
 * for a given accrual period.  
 * - adjFactor:   (O) adj. factor expressed as % of fee coupon
 * - valDate:     (I) valuation date
 * - accStDate:   (I) accrual start date
 * - accEndDate:  (I) accrual end date
 * - dcc:         (I) DCC (30/360, ACT/360, ACT/365)
 * - recovery:    (I) recovery rate 
 *
 * Return SUCCESS/FAILURE.
 */
int    PayAccrOnDefaultAdj(
    double        *adjFactor,
    TDate         today,           /* (I) compute price as of date        */
    TDate         valDate,         /* (I) return price as of date         */
    TDate         accStDate,
    TDate         accEndDate,
    TDayCountConv dcc,
    double        recovery,
    TDateInterval frequency,
    const TCurve       *discZC,          /* (I) Interest rate zero curve        */
    const TCurve       *sprdZC           /* (I) Clean spread zero curve         */
)
{
    long status = FAILURE;
    char routine[] = "PayAccrOnDefaultAdj";
    double one, pv1, pv2;

    one = 1.0;
    recovery = recovery; /* avoid warning on unused parameters */

    if (RiskyFeePV(
            &pv1, 
            today, 
            valDate, 
            1, 
            &accStDate, 
            &accEndDate, 
            &accEndDate, 
            &one, 
            &one, 
            dcc,
            ACCRUAL_PAY_NONE,
            SIMPLE, 
            frequency,         
            discZC, 
            sprdZC))
        goto RETURN;

    if (RiskyFeePV(
            &pv2, 
            today, 
            valDate, 
            1, 
            &accStDate, 
            &accEndDate, 
            &accEndDate, 
            &one, 
            &one, 
            dcc,
            ACCRUAL_PAY_ALL,
            SIMPLE, 
            frequency, 
            discZC, 
            sprdZC))
        goto RETURN;

    *adjFactor = pv1/pv2;

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
} 


/**
 * Compute the adjustment factor to par spread for pay accrue upon default 
 * for a given accrual period - INTERNAL FUNCTION
 * - adjFactor:   (O) adj. factor expressed as % of fee coupon
 * - valDate:     (I) valuation date
 * - accStDate:   (I) accrual start date
 * - accEndDate:  (I) accrual end date
 * - dcc:         (I) DCC (30/360, ACT/360, ACT/365)
 * - recovery:    (I) recovery rate 
 *
 * Return SUCCESS/FAILURE.
 */
int    RiskyDayCountFraction(
    double        *dcf,
    TDate         today,           /* (I) compute price as of date        */
    TDate         valDate,         /* (I) return price as of date         */
    TDate         accStDate,
    TDate         accEndDate,
    TDayCountConv dcc,
    double        recovery,
    TDateInterval frequency,
    const TCurve       *discZC,          /* (I) Interest rate zero curve        */
    const TCurve       *sprdZC           /* (I) Clean spread zero curve         */
    )
{
    long status = FAILURE;
    char routine[] = "RiskyDayCountFraction";
    double one, pv1, pv2;
    
    one = 1.0;
    recovery = recovery; /* avoid warning on unused parameters */

    if (RiskyFeePV(
            &pv1, 
            today, 
            valDate, 
            1, 
            &accStDate, 
            &accEndDate, 
            &accEndDate, 
            &one, 
            &one, 
            dcc,
            ACCRUAL_PAY_NONE,
            SIMPLE, 
            frequency, 
            discZC, 
            sprdZC))
        goto RETURN;

    if (RiskyFeePV(
            &pv2, 
            today, 
            valDate, 
            1, 
            &accStDate, 
            &accEndDate, 
            &accEndDate, 
            &one, 
            &one, 
            dcc,
            ACCRUAL_PAY_ALL,
            SIMPLE, 
            frequency, 
            discZC, 
            sprdZC))
        goto RETURN;

    if (GtoDayCountFraction(accStDate, accEndDate, dcc, dcf))
        goto RETURN;

    *dcf = pv1/pv2;

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;

} 

/**
 * Object version of CDS pricing
 */

int    CDSPV_O(
    double		*pv,
    TDate       today,
	TDate		valDate,
    const KCDS_D		*CDS,
	KStubType   priceConv,
    const TCurve		*discZC,
    const TCurve		*sprdZC
    )
{

    int status = FAILURE;
    char routine[] = "CDSPV_O";

    double price, price1;

    if (ProtectionPV_O(
        &price, 
        today, 
        valDate, 
        CDS->protLeg, 
        discZC, sprdZC) != SUCCESS)
    {
        DR_Error("%s: error pricing protection leg", routine);
        goto RETURN;
    }


    if (RiskyFeePV_O(
        &price1, 
        today, 
        valDate, 
        CDS->feeLeg, 
        priceConv, 
        discZC, 
        sprdZC) != SUCCESS)
    {
        DR_Error("%s: error pricing fee leg", routine);
        goto RETURN;
    }

    *pv = price - price1;

    status = SUCCESS;

RETURN:    

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
}


/*----------------------------------------------------------------------------------
 *
 * CDS Price - INTERNAL FUNCTION
 *
 * CDS price and first derivative wrt zero rate to maturity (last point on the input curve)
 *
 * Used by the solver inside the bootstrapping algorithm
 *
 */

int CDSPrice(
    double        *price,
    double        *dprice,
    double        *addedAnn,
    double        *addedProt,
    TDate         today,
    TDate         valDate,
    double        lastZero,
    double        parFee,
    double        recovery,
    double        overlapAnn,
    double        overlapProt,
    const KFeeLeg_D     *extraCoupons,
    const KProtLeg_D    *extraProtection,
    const TCurve        *irCurve, 
    TCurve        *crCurve
    )
{
	char routine[] = "CDSPrice";
    int status = FAILURE;

    double addedAnn1, addedProt1;
    
    
    crCurve->fArray[crCurve->fNumItems++].fDate = extraCoupons->mAccEndDates[extraCoupons->mNbCF - 1];
    crCurve->fArray[crCurve->fNumItems-1].fRate = lastZero;

    *price = -overlapAnn * parFee;
    
    RiskyFeePV_O(addedAnn, today, valDate, extraCoupons, NONE, irCurve, crCurve);

    *price -= *addedAnn * parFee;
    *price += overlapProt * (1 - recovery);

    ProtectionPV_O(addedProt, today, valDate, extraProtection, (TCurve*)irCurve, crCurve);

    *price += *addedProt * (1 - recovery);

    crCurve->fArray[crCurve->fNumItems-1].fRate += lastZero*0.01;

    if (RiskyFeePV_O(
        &addedAnn1,
        today,
        valDate, 
        extraCoupons,
        NONE,
        irCurve, 
        crCurve) != SUCCESS)
    {
        goto RETURN;
    }

    if (ProtectionPV_O(
        &addedProt1, 
        today,
        valDate, 
        extraProtection, 
        irCurve, 
        crCurve) != SUCCESS)
    {
        goto RETURN;
    }
    
    *dprice  = -*price;
    *dprice += (overlapProt * (1-recovery)  - overlapAnn * parFee);
    *dprice += (addedProt1 * (1-recovery)  - addedAnn1 * parFee);
    *dprice /= lastZero*0.01;

    status = SUCCESS;

RETURN:
    crCurve->fNumItems--;
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;

}

/*----------------------------------------------------------------------------------
 *
 * SolveParToClean - INTERNAL FUNCTION
 *
 * Takes in IR and credit zero curves, and the details of a par cds (coupon dates, last of which is the CDS maturity, 
 * and fee). Computes the zero rate to the CDS maturity which, when added to the input credit curve, makes the CDS par
 *
 * Coupon bank contains the first discount factors for the cds coupon dates,and the price of protection to the maturity 
 * of the previous cds. (the df's which depend on the part of the credit curve before its last date
 *
 */

int SolveParToClean(
    double          *zeroRate,        /* (O) output rate for couponDates  */
    double          *newAnnuity,
    double          *newProtection,
    TDate           today,
    TDate           valDate,
    int             doneCoupons,
    int             nbCoupons,        /* (O) total number of coupons      */
    const TDate           *accSttDates, 
    const TDate           *accEndDates, 
    const TDate           *payDates, 
    double          parFee, 
    double          recovery,
    TDayCountConv   dcc,
    KAccrualConv    payAccr,
    KProtPayConv    payType,
	TDateInterval   delay,
	TDateInterval   frequency, /* Integration frequency? */
    double          annuity,     
    double          protection,
    const TCurve          *irCurve, 
    TCurve          *crCurve
    )
{
    long status = FAILURE;
    char routine[] = "SolveParToClean";
    
    int i, newCoupons;
    KProtLeg_D *extraProt = NULL;
    KFeeLeg_D  *extraCoupons = NULL;

    double *notionals = NULL;
    double price;
    double dprice;
    double recoveryL;
    
    newCoupons = nbCoupons-doneCoupons;

	/* Allocate memory for notionals, and set them all to 1.0 */
    if (! (notionals  = (double*) malloc(newCoupons * sizeof(double)) ))
    {
        DR_Error("%s: failed to allocate memory.", routine);
        goto RETURN;
    }
    for (i = 0; i < newCoupons; i++) notionals[i] = 1;

    if (!(extraCoupons = RiskyFeeCreate(
        nbCoupons-doneCoupons, 
        accSttDates+doneCoupons, 
        accEndDates+doneCoupons, 
        payDates+doneCoupons, 
        notionals,                 /* notional is 1 */
        notionals,                 /* coupon rates are also 1 */
        dcc, 
        payAccr, 
        frequency)))
    {
        DR_Error("%s: failed to compute extra fee leg\n", routine);
        goto RETURN;
    };

    recoveryL = 0.0;

    if (!(extraProt = ProtectionCreate(
        accSttDates[doneCoupons], 
        1,                         /* only 1 period */
        accEndDates+nbCoupons-1,   /* last accrual date */
        notionals,                 /* notional is 1 */
        &recoveryL,                /* "pure" protection leg, recovery is 0 */
        payType,
        delay,
        frequency)))
    {
        DR_Error("%s: failed to compute extra protection leg\n", routine);
        goto RETURN;
    };

    *zeroRate = (crCurve->fNumItems == 0)? 0.01 : crCurve->fArray[crCurve->fNumItems-1].fRate;

    price = 1+TOL;

    i = 0;

    while (fabs(price) > SMALL && (i++ < MAX_ITER)) {
        
        if(CDSPrice(
            &price, 
            &dprice,
            newAnnuity,
            newProtection,
            today,
            valDate, 
            *zeroRate, 
            parFee,
            recovery,
            annuity,
            protection,
            extraCoupons,
            extraProt,
            irCurve,
            crCurve) != SUCCESS)
        {
			DR_Error("%s failed: could not price CDS for par fee %lf\n", routine, parFee);
            goto RETURN;
        };

        if (fabs(dprice) < TOL)
        {
			DR_Error("%s failed: derivative too small for Newton's method (%lf) on iteration %d with par fee %lf, recovery %lf, current price %lf, current zero-rate %lf\n", 
				routine, dprice, i, parFee, recovery, price, *zeroRate);
            goto RETURN;
        };

        *zeroRate -= price/dprice;
        
        if (-0.25 >= *zeroRate)
            *zeroRate = -0.25 + 1e-6;
    }

    if (fabs(price) > SMALL)
	{
		DR_Error("%s failed: price did not converge to zero (%lf)\n", routine, price);
        goto RETURN;
	}

    status = SUCCESS;
RETURN:  
	/* Free locally allocated objects */
    if (notionals) free(notionals);
    CrxProtectionFree(extraProt);
    CrxFeeLegFree(extraCoupons);
    
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;

}


/*------------------------------------------------------------------------------------
 *
 *  CreditZCBootstrap
 *
 *  This function creates a zero curve which prices to par a series of 
 *  CDS's with aligned cash-flows. It can also be used to bootstrap CDS's
 *  with non-aligned cash-flows, by calling it for one such CDS at a time, 
 *  passing the curve constructed by the previous call as input.
 *
 *  crCurveStub either NULL or points to a valid TCurve
 *
 *  Returns NULL (if fails), or a pointer to a valid TCurve (which the caller 
 *  is responsible for freeing)
 */

TCurve* CreditZCBootstrap(
    TDate         today,
	TDate         valDate,
	TDayCountConv dcc,
	KAccrualConv  accrual,
    KProtPayConv  payType,
	TDateInterval delay,
    int           nbCouponDates,
    const TDate         *accSttDates, 
    const TDate         *accEndDates, 
    const TDate         *payDates, 
    int           nbCDS,
    const int           *maturityIdx,
	const double        *parFees,
	double        recovery,
    TDateInterval frequency, /* Integration frequency */
	const TCurve        *irCurve,
	const TCurve        *crCurveStub
)
{

    
    char *routine = "CreditZCBootstrap";


    int cdsidx, npoints, i; 

    double zeroRate;
    double annuity, protection, newAnnuity, newProtection;

    TCurve *cveTmp = NULL;
    TCurve *crCurve = NULL;

	/* Do some basic consistency checks */
    for (i = 0; i < nbCouponDates; i++)
    {
        if (accSttDates[i] >= accEndDates[i])
        {
            DR_Error("%s: accrual start must be before accrual end.", routine);
            goto RETURN;
        }

        if (accEndDates[i] > accEndDates[i])
        {
            DR_Error("%s: payment date cannot be before accrual end.", routine);
            goto RETURN;
        }

        if (i)
        {
            if (accSttDates[i] <= accSttDates[i-1])
            {
                DR_Error("%s: accrual start dates must be increasing.", routine);
                goto RETURN;
            }

            if (accEndDates[i] <= accEndDates[i-1])
            {
                DR_Error("%s: accrual end dates must be increasing.", routine);
                goto RETURN;
            }

            if (accEndDates[i] <= accEndDates[i-1])
            {
                DR_Error("%s: payment dates must be increasing.", routine);
                goto RETURN;
            }
        }
    }

    for (i =0; i < nbCDS; i++)
    {
        if (maturityIdx[i] >= nbCouponDates)
        {
            DR_Error("%s : cds idx exceeds number of coupon dates", routine);
            goto RETURN;
        }

        if ( i> 0)
            if (maturityIdx[i] <= maturityIdx[i-1])
        {
            DR_Error("%s : cds indices must be increasing.", routine);
            goto RETURN;
        }
    }

    npoints = 0;
    
    if (crCurveStub != NULL)
        npoints = crCurveStub->fNumItems;
    
    npoints += nbCDS;


	/* Create a TCurve to hold the new credit-spread information */
	crCurve = GtoNewTCurve(irCurve->fBaseDate, npoints, 1, GTO_ACT_365F);
	if (crCurve==NULL) {
		DR_Error("%s failed: Could not create credit TCurve with %d points for base-date %d", 
			routine, npoints, irCurve->fBaseDate);
        goto RETURN;
	}
	/* Because the number of points is used during the bootstrapping, you have to set
		it to zero here */
	crCurve->fNumItems = 0;

    
    if (crCurveStub)
    {
        for (i=0; i < crCurveStub->fNumItems;i++)
        {
            crCurve->fArray[i] = crCurveStub->fArray[i];
        }

        crCurve->fNumItems = crCurveStub->fNumItems;
    }
    
    /* initialize crCurve with crCurveStub, but capacity for (crCurveStub->fNumItems + nbCDS) points  */
    annuity = protection = 0;

    for (cdsidx = 0; cdsidx < nbCDS; cdsidx++)
    {

        if(SolveParToClean(
            &zeroRate,
            &newAnnuity,
            &newProtection,
            today,
            valDate,
            (cdsidx==0) ? 0 : (maturityIdx[cdsidx-1]+1),
            maturityIdx[cdsidx]+1, 
            accSttDates, 
            accEndDates, 
            payDates,
            parFees[cdsidx],
            recovery, 
            dcc,
            accrual, 
            payType, 
            delay, 
            frequency, 
            annuity, 
            protection,
            irCurve, 
            crCurve) == FAILURE)
        {
			DR_Error("%s failed: failed to solve for clean spread %d of %d with par fee %lf.\n", routine, cdsidx, nbCDS, parFees[cdsidx]);
            goto RETURN;
        }

        crCurve->fArray[crCurve->fNumItems].fDate = accEndDates[maturityIdx[cdsidx]];
        crCurve->fArray[crCurve->fNumItems].fRate = zeroRate;
		crCurve->fNumItems++; // increment for bootstrapper

        annuity     += newAnnuity;
        protection  += newProtection;

    }

	/* Shift curve to CDS spot value-date */
    cveTmp      = GtoZCShift(crCurve, valDate);

RETURN:    
    GtoFreeTCurve(crCurve);
    return cveTmp;
}    
    

/** Creates payment dates for an array of CDS. The payment dates are put in dates[0][]. The number of payment
 *  dates for CDS i is index[0][i], so that maturities[j] = dates[0][index[0][j]-1]. The caller
 * is responsible for freeing dates and index when finished with them. */
int CDSCouponSchedule(
    TDate         **dates,
    int           **index,     
    TDate         sttDate,     // start Date
    TDateInterval frequency,   // coupon frequency
    int           nbCDS,       // number of CDS to generate coupon schedules for
    const TDate   *maturities) // array of CDS maturities
{
    char routine[] = "CDSCouponSchedule";
    long status = FAILURE;

    int i, j;
    long num;

    TDate tmp, couponDate = sttDate;

	/* This finds the final date and the number of dates starting from sttDate
	 * and going to maturities[nbCDS-1] in jumps of frequency */
    DateFloor(&tmp, &num, sttDate, maturities[nbCDS-1], frequency);
	/* If off by one (e.g. stub), increase number by one */
    if (tmp < maturities[nbCDS-1])
        num++;

    *dates = (TDate*) malloc( (num+1)* sizeof(TDate) ); // coupon dates
    *index = (int*) malloc( nbCDS* sizeof(int) );       // index[0][i] is the number of coupon dates for the ith CDS

    j = 0;

    for (i = 0; i < nbCDS; i++)
    {
        while (couponDate <= maturities[i])
        {
            dates[0][j] = couponDate;

			/* If you have reached the maturity, then stop */
            if (couponDate == maturities[i]) break;

			/* Add frequency*(j+1) to the base date */
            if (GtoDateFromDateAndOffset(dates[0][0], &frequency, ++j, &couponDate)) goto RETURN;

			/* If get stub at end, make sure that coupon dates can't overrun maturities */
            if (couponDate >= maturities[i]) couponDate = maturities[i];

        }

        index[0][i] = j-1; // only index 0 is ever filled. Why should index be [][], rather than just []?
    }

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);
    
    return status;
}

/**
* Builds a clean spread curve by offsetting the par-spreads by a fixed amount 
* relative to those implied by a given credit curve (crCurve). If crCurve is null,
* instead of offsetting by a fixed basis, a flat par-spread curve is built
* This uses the same contingent-leg integration frequency as the CDS coupon interval - i.e. usually 3M
*/
TCurve* CrxBuildCurveFromFlatParSpreadOffset (
	TDate today,
	TDate valueDate,
	TDayCountConv dcc,
	KAccrualConv accrual,
	KProtPayConv payType,
	TDateInterval delay,
	double parSpreadOffset,
	double recovery,
	TDateInterval frequency, // this is the coupon interval
    TDateInterval integFreq,  // Protection leg integration frequency
	const TCurve *irCurve,
	const TCurve *crCurve
	) {

	char          routine[] = "CrxBuildCurveFromFlatParSpreadOffset";               
	TCurve        *crCurveOut = NULL; // the resulting credit curve, which is returned.

    int           i;                  // loop counter
    TDate         *coupons = NULL;    // coupon payment dates
	int           *index = NULL;      // index into coupons showing end of ith CDS
	TDate         tenors[MAX_CDS];    // maturities of CDS
	double        parSpread[MAX_CDS];  // par-spreads of CDS
	/* Base new curve on spreads at these tenors */
	int ntenors = 14;
    char *tenorNames[] = {"3M", "6M", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y"};

    TDateInterval matInterval;
    KStubLocation stub = LONG_FRONT;
    double annuity;
    
	

	/* Compute par spreads for each tenor and add offset */
    for (i = 0; i < ntenors; i++) {


		if (FAILURE==GtoStringToDateInterval(tenorNames[i], tenorNames[i] /* for error msg */, &matInterval)) 
		{
			DR_Error("%s: failed to compute maturity interval for cds %d, tenor %s.\n", routine, i, tenorNames[i]);
			goto RETURN;
		}
		if (FAILURE==GtoDateFromDateAndOffset((crCurve ? crCurve->fBaseDate : valueDate), &matInterval, 1, tenors+i))
		{
			DR_Error("%s: failed to compute maturity date for cds %d, tenor %s.\n", routine, i, tenorNames[i]);
			goto RETURN;
		}

		
		if (crCurve) {

			
			/* If you have a valid curve, the new par spreads are the old ones plus the offset */
			if (CrxFwdParCDSSpread(
				&parSpread[i], // write par spread to parSpread[i]
				&annuity,      // do nothing with this, but you get it anyway
				today, 
				valueDate, 
				valueDate,
				tenors[i],
				frequency,
				dcc, 
				stub, 
				accrual, 
				payType, 
				delay, 
				recovery,
				integFreq, 
				irCurve, 
				crCurve
				) != SUCCESS)
			{
				DR_Error("%s: failed to compute par spread for tenor %s at date %d.\n", routine, tenorNames[i], tenors[i]);
				goto RETURN;
			}
		

			parSpread[i] += parSpreadOffset;

		} else {
			/* If no crCurve was passed, then you are creating a flat par-spread curve */
			parSpread[i] = parSpreadOffset;
		}
    }

	/* Next, create all the payment dates for the CDS (allocate and initialize coupons & index) */
    if (CDSCouponSchedule(
        &coupons, 
        &index, 
        valueDate, 
        frequency, 
        ntenors, 
        tenors) != SUCCESS)
    {
        DR_Error("%s: failed to generate coupon schedule.", routine);
        goto RETURN;
    };

    
	/* Finally, create the new credit curve */
	crCurveOut = CreditZCBootstrap(
        today,
        valueDate,
        dcc, 
        accrual,
        payType, 
        delay, 
        index[ntenors-1]+1, // this is the number of coupon dates in the longest CDS
        coupons,   // start dates
        coupons+1, // end dates
        coupons+1, // pay dates
        ntenors, 
        index, 
        parSpread, // new, offset, par spreads
        recovery, 
        integFreq, 
        irCurve, 
        NULL);
    if (crCurveOut==NULL)
    {
        DR_Error("%s: failed to bootstrap new offset credit curve.\n", routine);
        goto RETURN;
    }


RETURN:
    if (coupons!=NULL) free(coupons);
    if (index!=NULL) free(index);
    return crCurveOut;

}

int    CrxFwdParCDSSpread(
    double        *parSprd,
    double        *annuity,
    TDate         today,
	TDate         valDate,
    TDate         cdsStDate,
    TDate         cdsMatDate,
    TDateInterval cdsFreq,
    TDayCountConv cdsDCC,
	KStubLocation dateStub,
    KAccrualConv  idxPayAccFlag,
    KProtPayConv  protPayType,
	TDateInterval delay,
    double        recovery,        /* (I) Recovery rate                    */
	TDateInterval frequency,       /* (I) Frequency of input curve (I think this is the integration frequency CM) */ 
    const TCurve       *discZC,          /* (I) Interest rate zero curve         */
    const TCurve       *sprdZC)          /* (I) Clean spread zero curve          */
{
    char routine[] = "CrxFwdParCDSSpread";
    long status = FAILURE;

    double        one = 1;
    KFeeLeg_D     *feeLeg = NULL;


    if ((feeLeg = CrxFeeLegCreateFromFreq(
                        cdsStDate, 
                        cdsMatDate, 
                        cdsFreq, 
                        dateStub, 
                        1.0,
                        1.0, 
                        cdsDCC, 
                        idxPayAccFlag, 
                        frequency) ) == NULL)
    {
        DR_Error("%s: failed to create fee leg.", routine);
        goto RETURN;
    }

    if (ProtectionPV(
        parSprd,
        today,
        valDate, 
        cdsStDate, 
        1, 
        &cdsMatDate, 
        &one, 
        &recovery, 
        protPayType, 
        delay, 
        frequency, 
        discZC, 
        sprdZC) != SUCCESS)
    {
        DR_Error("%s: error computing protection leg", routine);
        goto RETURN;
    }


    if (RiskyFeePV_O(
        annuity,
        today,
        valDate, 
        feeLeg,
        NONE,
        discZC, 
        sprdZC) != SUCCESS)
    {
        DR_Error("%s: error computing fee leg", routine);
        goto RETURN;
    }

    *parSprd /=  *annuity;

    status  = SUCCESS;

RETURN:

    CrxFeeLegFree(feeLeg);
    
    if (status != SUCCESS)
        DR_Error ("%s: failed", routine);

    return status;
    
} 



/** Converts a curve calibrated at one contingent-leg integration frequency to another, new one. 
 Does so that calibration to various (hard-coded) par-fee tenors is consistent. 
 Assumes that recovery is always 50%, which does not seem like a good idea (CM 3/23/05) */
TCurve* CrxCleanCDSCurveConvert(
 	TDateInterval newFrequency,    /* (I) Frequency for output curve       */
 	TDateInterval frequency,       /* (I) Frequency of input curve         */
    double        recovery,        /* Assumed recovery rate */
    const TCurve       *discZC,          /* (I) Interest rate zero curve         */
    const TCurve       *sprdZC)          /* (I) Clean spread zero curve          */
{
    
	char routine[] = "CrxCleanCDSCurveConvert";
	TCurve          *adjSprdZC = NULL; // the resulting adjusted curve, if succeeds
	int             i; // loop counting variable
	int             ntenors = 14;
	TDate           *coupons = NULL;
	int             *index = NULL;

    TDate tenors[MAX_CDS];

    TDateInterval matInterval, noDelay, couponInterval;

    KStubLocation stub = LONG_FRONT;

    double parSpread[MAX_CDS], annuity;
    long dcc;
    
    char *tenorNames[] = {"3M", "6M", "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "15Y", "20Y"};
    dcc = GTO_ACT_360;
    //recovery = 0.5;

    GtoStringToDateInterval("3M", "", &couponInterval);
    GtoStringToDateInterval("0D", "", &noDelay);

    for (i = 0; i < ntenors; i++)
    {

        GtoStringToDateInterval(tenorNames[i], tenorNames[i], &matInterval);
        
        GtoDateFromDateAndOffset(sprdZC->fBaseDate, &matInterval, 1, tenors+i);

        if (CrxFwdParCDSSpread(
            parSpread+i, 
            &annuity, 
            sprdZC->fBaseDate, 
            sprdZC->fBaseDate, 
            sprdZC->fBaseDate,
            tenors[i],
            couponInterval, 
            dcc, 
            stub, 
            ACCRUAL_PAY_ALL, 
            PAY_DEF, 
            noDelay, 
            recovery,
            frequency, 
            discZC, 
            sprdZC
            ) != SUCCESS)
        {
            DR_Error("%s: failed to compute par spread curve.", routine);
            goto RETURN;
        }
    }

	/* Create coupon schedule - initialize coupons and index */
    if (CDSCouponSchedule(
        &coupons, 
        &index, 
        sprdZC->fBaseDate, 
        couponInterval, 
        ntenors, 
        tenors) != SUCCESS)
    {
        DR_Error("%s: failed to generate coupon schedule.", routine);
        goto RETURN;
    };
    
	adjSprdZC = CreditZCBootstrap(
        sprdZC->fBaseDate,
        sprdZC->fBaseDate,
        dcc, 
        ACCRUAL_PAY_ALL,
        PAY_DEF, 
        noDelay, 
        index[ntenors-1]+1, 
        coupons, 
        coupons+1, 
        coupons+1, 
        ntenors, 
        index, 
        parSpread, 
        recovery, 
        newFrequency, 
        discZC, 
        NULL);
    if (adjSprdZC==NULL)
    {
        DR_Error("%s: failed to bootstrap adjusted curve.", routine);
        goto RETURN;
    }

RETURN:

    if (coupons!=NULL) free(coupons);
    if (index!=NULL) free(index);

    return adjSprdZC;

}
