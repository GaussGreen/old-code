/*
***************************************************************************
** FILENAME: crxconv.c
**
** Converts between CX and CRX data structures.
**
** $Header$
***************************************************************************
*/

#include "crxconv.h"

#include "recovery.h"
#include <cxutils/include/cxmacros.h>

KProtLeg_D* CxProtectionLegConvert
(CxTContingentLeg *contingentLeg,       /* (I) Contingent leg       */
 CxTRecoveryCurve *recoveryCurve,       /* (I) Recovery curve       */
 TDateInterval    *integrationInterval  /* (I) Integration interval */
)
{
    static char routine[] = "CxProtectionLegConvert";
    KProtLeg_D *out = NULL;

    int           i;
    double       *recoveries = NULL;
    TDateInterval delayInterval;

    /*
    ** There is a case to be made for extending the list of dates to match
    ** the integration interval if the loss curve is not flat. (TBD)
    */

    REQUIRE (contingentLeg != NULL);
    REQUIRE (integrationInterval != NULL);

    recoveries = NEW_ARRAY (double, contingentLeg->nbDates);
    
    if (recoveryCurve != NULL)
    {
        for (i = 0; i < contingentLeg->nbDates; ++i)
        {
            double recovery;
            
            if (CxRecoveryCurveInterp (recoveryCurve, contingentLeg->dates[i], 
                                       &recovery) != SUCCESS)
                goto done; /* failure */

            recoveries[i] = recovery;
        }
    }

    SET_TDATE_INTERVAL (delayInterval, contingentLeg->payDelay, 'D');

    out = ProtectionCreate (contingentLeg->startDate,
                            contingentLeg->nbDates,
                            contingentLeg->dates,
                            contingentLeg->notionals,
                            recoveries,
                            (KProtPayConv)contingentLeg->payType,
                            delayInterval,
                            *integrationInterval);
done:

    FREE (recoveries);
    if (out == NULL)
        GtoErrMsgFailure (routine);

    return out;
}
        

KFeeLeg_D* CxFeeLegConvert
(CxTFeeLeg     *feeLeg,             /* (I) Fee leg              */
 TDateInterval *integrationInterval /* (I) Integration interval */
)
{
    static char routine[] = "CxFeeLegConvert";
    KFeeLeg_D  *out = NULL;

    REQUIRE (feeLeg != NULL);
    REQUIRE (integrationInterval != NULL);

    out = RiskyFeeCreate (feeLeg->nbDates,
                          feeLeg->accStartDates,
                          feeLeg->accEndDates,
                          feeLeg->payDates,
                          feeLeg->notionals,
                          feeLeg->couponRates,
                          feeLeg->dcc,
                          (KAccrualConv)feeLeg->accrualPayConv,
                          *integrationInterval);
done:

    if (out == NULL)
        GtoErrMsgFailure (routine);

    return out;
}
    
