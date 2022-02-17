/*
***************************************************************************
** FILE NAME: creditcurve.c
**
** Credit curve functions.
**
** $Header$
***************************************************************************
*/

#include "creditcurve.h"

#include <cxutils/include/alibconv.h>

#include <math.h>

#include <alib/dtivlo.h>
#include <alib/gtomathp.h>
#include <alib/ldate.h>
#include <alib/tcurve.h>
#include <cxutils/include/cxmacros.h>
#include <cxutils/include/zerocurve.h>

/*f
***************************************************************************
** Validates (and potentially changes) a credit curve as part of its
** construction.
***************************************************************************
*/
int CxCreditCurveValidate
(CxTCreditCurve *creditCurve)
{
    static char routine[] = "CxCreditCurveValidate";
    int         status    = FAILURE;

    REQUIRE (creditCurve != NULL);

    if (creditCurve->timestep == NULL)
    {
        /* default of 1W chosen to match CDS option pricer default */
        creditCurve->timestep = GtoDateIntervalNew (1, "W");
        if (creditCurve->timestep == NULL)
            goto done; /* failure */
    }
    
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*f
***************************************************************************
** Constructs a credit curve. Easier to use this in a spreadsheet
** environment than first constructing a TCurve.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveNew
(CxTCreditCurveType type,
 TDate             baseDate,
 int               numDates,
 TDate            *dates,
 double           *spreads,
 long              rateType,
 CxTDayCountConv    dayCountConv,
 TDateInterval    *ivl)
{
    static char routine[] = "CxCreditCurveNew";

    TCurve        *tc = NULL;
    CxTCreditCurve *cc = NULL;

    tc = GtoMakeTCurve (baseDate,
                        dates,
                        spreads,
                        numDates,
                        rateType,
                        (long)dayCountConv);
    if (tc == NULL)
        goto done; /* failure */

    cc = CxCreditCurveMake (type, tc, ivl);

 done:

    GtoFreeTCurve (tc);

    if (cc == NULL)
        GtoErrMsgFailure (routine);

    return cc;
}

CxTCreditCurve* CxCreditCurveMakeSmoother
(CxTCreditCurve *cc,
 TDateList    *dl)
{
    static char routine[] = "CxCreditCurveMakeSmoother";

    TCurve *tc = NULL;
    CxTCreditCurve *out = NULL;

    REQUIRE (cc != NULL);
    
    tc  = CxZeroCurveMakeSmoother (cc->tc, dl);
    if (tc == NULL)
        goto done; /* failure */
    out = CxCreditCurveMake (cc->type, tc, cc->timestep);

 done:

    GtoFreeTCurve (tc);
    if (out == NULL)
        GtoErrMsgFailure (routine);

    return out;
}
 
/*f
***************************************************************************
** Coercion from TCurve to CxTCreditCurve.
**
** Some of our existing CX routines use TCurve internally (this will
** change), but externally we represent via CxTCreditCurve with annual
** rates instead of continuous rates.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveFromZeroCurve
(TCurve *zc)
{
    CxTCreditCurve *cc = CxCreditCurveMakeEmpty();

    cc->type = CX_CURVE_TYPE_FLOW;
    cc->tc   = GtoCopyCurve (zc);
    cc->timestep = NULL;
    if (CxCreditCurveValidate (cc) != SUCCESS)
    {
        CxCreditCurveFree (cc);
        cc = NULL;
    }

    return cc;
}

/*f
***************************************************************************
** Coercion from TCurve to CxTCreditCurve
**
** The main purpose of this routine is to enable existing test cases to
** work where previously the credit curve was represented as TCurve.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveFromTCurve
(TCurve *tc)
{
    return CxCreditCurveMake (CX_CURVE_TYPE_EXOTIC, tc, NULL);
}

#if 0
TCurve* CxCreditCurveToTCurve(CxTCreditCurve *cc)
{
    static char routine[] = "CxCreditCurveToZeroCurve";
    
    TCurve *zc = NULL;

    REQUIRE (cc != NULL);
    REQUIRE (cc->type == CX_CURVE_TYPE_FLOW);
    
    zc = CxZeroCurveFromTCurve (cc->tc, GTO_FLAT_FORWARDS);

 done:

    if (zc == NULL)
        GtoErrMsgFailure (routine);

    return zc;
}
#endif

TZeroCurve* CxCreditCurveToZeroCurve
(CxTCreditCurve *cc)
{
    return GtoZeroCurveFromTCurveNew (GtoCopyCurve (cc->tc),
                                      GTO_FLAT_FORWARDS);
}

/*f
***************************************************************************
** Routine for converting the compounding basis of a credit curve.
**
** The curve is amended in place.
***************************************************************************
*/
int CxCreditCurveConvertRateType
(CxTCreditCurve *cc, 
 long            newRateType)
{
    static char routine[] = "CxCreditCurveConvertRateType";
    int         status    = FAILURE;

    int i;

    REQUIRE (cc != NULL);

    if (IS_EQUAL(newRateType, cc->tc->fBasis))
    {
        /* do nothing */
    }
    else
    {
        for (i = 0; i < cc->tc->fNumItems; ++i)
        {
            if (CxConvertCompoundRate (cc->tc->fArray[i].fRate,
                                       cc->tc->fBasis,
                                       cc->tc->fDayCountConv,
                                       newRateType,
                                       cc->tc->fDayCountConv,
                                       &cc->tc->fArray[i].fRate) != SUCCESS)
                goto done; /* failure */
        }
        cc->tc->fBasis = newRateType;
    }
    
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

 
/**
***************************************************************************
** Returns the hazard rate for a particular date from a credit curve.
**
** Since survival probabilities in the credit curve are valued at the end
** of day, this calculation involves computing the survival probability from
** one day before and the current date, computing the one day hazard,
** and then converting this into an annualized exponential rate.
**
** For dates on or before the start date of the curve, the hazard rate is
** zero.
***************************************************************************
*/
int CxCreditCurveHazardRate
(CxTCreditCurve *creditCurve,
 TDate           date,
 double         *hazardRate)
{
    static char routine[] = "CxCreditCurveHazardRate";
    int         status    = FAILURE;

    double      oneDaySurvival;

    REQUIRE (creditCurve != NULL);
    REQUIRE (date >= 0);
    REQUIRE (hazardRate != NULL);

    if (CxCreditCurveConditionalSurvivalProb (creditCurve, 
                                              date-1, 
                                              date,
                                              &oneDaySurvival) != SUCCESS) 
        goto done; /* failure */

    if (ARE_ALMOST_EQUAL(oneDaySurvival,1.0))
    {
        *hazardRate = 0.0;
    }
    else
    {
        *hazardRate = -365.0 * log(oneDaySurvival);
    }
    
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}
    

/**
***************************************************************************
** Returns the conditional survival probability computed at the end of the
** day for endDate, conditional on survival to the end of the day for
** startDate.
***************************************************************************
*/
int CxCreditCurveConditionalSurvivalProb
(CxTCreditCurve *creditCurve,
 TDate           startDate,
 TDate           endDate,
 double         *survivalProb)
{
    static char routine[] = "CxCreditCurveConditionalSurvivalProb";
    int         status    = FAILURE;

    double      startSurvival;
    double      endSurvival;
 
    REQUIRE (creditCurve != NULL);
    REQUIRE (creditCurve->tc != NULL);
    REQUIRE (startDate >= 0);
    REQUIRE (endDate > startDate);
    REQUIRE (survivalProb != NULL);
    
    if (CxCreditCurveSurvivalProb (creditCurve, startDate, 
                                   &startSurvival) != SUCCESS)
        goto done; /* failure */

    if (CxCreditCurveSurvivalProb (creditCurve, endDate, 
                                   &endSurvival) != SUCCESS)
        goto done; /* failure */
 
    *survivalProb = endSurvival / startSurvival;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/**
***************************************************************************
** Returns the survival probability computed at the end of the day.
** This is computed relative to the base date of the credit curve.
***************************************************************************
*/
int CxCreditCurveSurvivalProb
(CxTCreditCurve *creditCurve,
 TDate           date,
 double         *survivalProb)
{
    static char routine[] = "CxCreditCurveSurvivalProb";
    int         status    = FAILURE;

    REQUIRE (creditCurve != NULL);
    REQUIRE (creditCurve->tc != NULL);
    REQUIRE (date >= 0);
    REQUIRE (survivalProb != NULL);
    
    if (date <= creditCurve->tc->fBaseDate)
    {
        *survivalProb = 1.0;
    }
    else
    {
        *survivalProb = CxZeroPrice (creditCurve->tc, date);
        if (GTO_ISNAN(*survivalProb))
            goto done; /* failure */
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}


/**
***************************************************************************
** Slides the credit curve to a new date in the future.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveSlide
(CxTCreditCurve *creditCurve,
 TDate           slideDate)
{
    static char routine[] = "CxCreditCurveSlide";
    int         status    = FAILURE;

    CxTCreditCurve *slideCurve = NULL;
    long            offset = 0;
    int             i;

    REQUIRE (creditCurve != NULL);
    REQUIRE (creditCurve->tc != NULL);

    slideCurve = CxCreditCurveCopy(creditCurve);

    offset     = slideDate - slideCurve->tc->fBaseDate;

    slideCurve->tc->fBaseDate = slideDate;

    for (i = 0; i < slideCurve->tc->fNumItems; ++i)
    {
        slideCurve->tc->fArray[i].fDate += offset;
    }
    
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        CxCreditCurveFree (slideCurve);
        slideCurve = NULL;
        GtoErrMsgFailure (routine);
    }

    return slideCurve;
}
    


