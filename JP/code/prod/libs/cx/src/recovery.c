/*
***************************************************************************
** SOURCE FILE: recovery.c
**
** Defines an extremely simple recovery curve structure.
**
** $Header$
***************************************************************************
*/

#include "recovery.h"

#include <assert.h>

#include "cxutils/include/cxmacros.h"
#include "cxutils/include/bsearch.h"
#include "cxutils/include/date.h"
#include "cxutils/include/datelist.h"

#include <alib/sintrp.h> 
#include <alib/today.h>


/*
***************************************************************************
** Constructs a recovery curve from a single recovery rate.
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveMakeFromRecoveryRate
(double recoveryRate)
{
    static char routine[] = "CxRecoveryCurveMakeFromRecoveryRate";

    CxTRecoveryCurve* recoveryCurve = NULL;
    TDate today;

    if (GtoToday(&today) != SUCCESS)
        goto done; /* failure */

    recoveryCurve = CxRecoveryCurveMake (
        1,
        &today,
        &recoveryRate,           
        CX_RECOVERY_LEFT_INTERP, /*interpType*/
        TRUE,                    /*extrapBefore*/
        TRUE);                   /*extrapAfter*/

 done:

    if (recoveryCurve == NULL)
        GtoErrMsgFailure (routine);

    return recoveryCurve;
}


/*
***************************************************************************
** Ensure that the dates in the recovery curve are in the correct order,
** plus various other validations.
***************************************************************************
*/
int CxRecoveryCurveValidate
(CxTRecoveryCurve* curve)
{
    static char routine[] = "CxRecoveryCurveValidate";
    int         status    = FAILURE;

    int i;

    REQUIRE(curve->numItems > 0);
    REQUIRE(curve->dates != NULL);
    REQUIRE(curve->recoveryRates != NULL);

    for (i = 0; i < curve->numItems; ++i)
    {
        if (i > 0)
        {
            if (curve->dates[i] <= curve->dates[i-1])
            {
                char d1[16];
                char d2[16];
                GtoErrMsg ("%s: Dates %s and %s are not in ascending order\n",
                           routine,
                           CxDateFormat(curve->dates[i-1], "DD-MMM-YYYY", d1),
                           CxDateFormat(curve->dates[i],   "DD-MMM-YYYY", d2));
                goto done; /* failure */
            }
        }
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure(routine);
    return status;
}

/*
***************************************************************************
** Interpolates a value from a recovery curve.
***************************************************************************
*/
int CxRecoveryCurveInterp
(const CxTRecoveryCurve* curve,        /* (I) Recovery curve */
 TDate                   date,         /* (I) Interpolation date */
 double                 *interpValue)  /* (O) Interpolated value */
{
    static char routine[] = "CxRecoveryCurveInterp";
    int         status    = FAILURE;
    
    long        exact;
    long        lo;
    long        hi;
    double      recoveryRate = 0.0;
    
    REQUIRE (curve != NULL);
    REQUIRE (curve->numItems > 0);
    REQUIRE (curve->dates != NULL);
    REQUIRE (curve->recoveryRates != NULL);

    if (CxBinarySearchLong (date,
                            curve->dates,
                            sizeof(TDate),
                            curve->numItems,
                            &exact,
                            &lo,
                            &hi) != SUCCESS) 
        goto done; /* failure */

    if (exact >= 0)
    {
        /* date found in dates */
        recoveryRate = curve->recoveryRates[exact];
    }
    else if (lo < 0)
    {
        /* date before start of dates */
        if (curve->extrapBefore)
        {
            recoveryRate = curve->recoveryRates[0];
        }
        else
        {
            char d1[16];
            char d2[16];
            GtoErrMsg ("%s: Requested date %s is before start date %s\n",
                       routine,
                       CxDateFormat(date,            "DD-MMM-YYYY", d1),
                       CxDateFormat(curve->dates[0], "DD-MMM-YYYY", d2));
            goto done; /* failure */
        }
    }
    else if (hi >= curve->numItems)
    {
        /* date after end of dates */
        if (curve->extrapAfter)
        {
            recoveryRate = curve->recoveryRates[curve->numItems-1];
        }
        else
        {
            char d1[16];
            char d2[16];
            GtoErrMsg ("%s: Requested date %s is after end date %s\n",
                       routine,
                       CxDateFormat(date, "DD-MMM-YYYY", d1),
                       CxDateFormat(curve->dates[curve->numItems-1],
                                    "DD-MMM-YYYY", d2));
            goto done; /* failure */
        }
    }
    else
    {
        /* interpolation required */
        /* everything else we can do directly */
        switch (curve->interpType)
        {
        case CX_RECOVERY_LINEAR_INTERP:
        {
            double loval = curve->recoveryRates[lo];
            double hival = curve->recoveryRates[hi];
            double lodate = curve->dates[lo];
            double hidate = curve->dates[hi];
            double fract;

            fract = (double)(date-lodate) / (double)(hidate-lodate);
            recoveryRate = loval + (hival - loval) * fract;
            break;
        }
        case CX_RECOVERY_RIGHT_INTERP:
            recoveryRate = curve->recoveryRates[hi];
            break;
        case CX_RECOVERY_LEFT_INTERP:
            recoveryRate = curve->recoveryRates[lo];
            break;
        /* no default - intelligent compilers can spot missing enum values */
        }
    }
        
    *interpValue = recoveryRate;
    status       = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*
***************************************************************************
** Computes an unweighted average loss between two dates.
***************************************************************************
*/
int CxRecoveryCurveAverage
(const CxTRecoveryCurve *recoveryCurve, /* (I) recovery curve */
 TDate                   startDate,     /* (I) Start date */
 TDate                   endDate,       /* (I) End date */
 double                 *averageValue)  /* (O) Average value */
{
    static char routine[] = "CxRecoveryCurveAverage";
    int         status    = FAILURE;

    double  total = 0.0;
    int     i;

    double loss0;
    double loss1;

    TDateList  *tl = NULL;

    REQUIRE (endDate > startDate);
    REQUIRE (averageValue != NULL);
    REQUIRE (recoveryCurve != NULL);

    tl = GtoNewDateListFromDates (recoveryCurve->dates, 
                                  recoveryCurve->numItems);
    if (tl == NULL) goto done; /* failure */
    tl = CxDateListAddDatesFreeOld (tl, 1, &startDate);
    if (tl == NULL) goto done; /* failure */
    tl = CxDateListAddDatesFreeOld (tl, 1, &endDate);
    if (tl == NULL) goto done; /* failure */
    tl = CxDateListTruncate (tl, startDate, TRUE, TRUE, TRUE);
    tl = CxDateListTruncate (tl, endDate, TRUE, FALSE, TRUE);

    assert (tl->fArray[0] == startDate);
    assert (tl->fArray[tl->fNumItems-1] == endDate);

    /* the integration - our interpolation functions are very simple, so
       we can average it out very easily between points.
    */

    if (CxRecoveryCurveInterp (recoveryCurve, startDate, &loss1) != SUCCESS)
        goto done; /* failure */

    for (i = 1; i < tl->fNumItems; ++i)
    {
        double thisAverage = 0.0;

        loss0 = loss1;
        if (CxRecoveryCurveInterp(recoveryCurve, tl->fArray[i], 
                                  &loss1) != SUCCESS)
            goto done; /* failure */
        
        switch (recoveryCurve->interpType)
        {
        case CX_RECOVERY_LINEAR_INTERP:
            thisAverage = (loss0 + loss1) * 0.5;
            break;
        case CX_RECOVERY_RIGHT_INTERP:
            thisAverage = loss1;
            break;
        case CX_RECOVERY_LEFT_INTERP:
            thisAverage = loss0;
            break;
        /* no default - intelligent compilers can spot missing enum values */
        }
        total   += thisAverage * (tl->fArray[i] - tl->fArray[i-1]);
    }

    status = SUCCESS;
    *averageValue = total / (endDate - startDate);
        
 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    GtoFreeDateList (tl);
    return status;
}


