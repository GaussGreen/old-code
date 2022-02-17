/*
***************************************************************************
** SOURCE FILE: cx.c
***************************************************************************
*/

#include "cx.h"
#include "creditcurve.h"
#include "recovery.h"

#include <ctype.h>
#include <alib/convert.h>
#include <alib/dtivlo.h>
#include <alib/tcurve.h>
#include <cxutils/include/cxmacros.h>

/**
***************************************************************************
** Converts CreditCurveType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxCreditCurveTypeToString (CxTCreditCurveType value)
{
    static char routine[] = "CxCreditCurveTypeToString";

    switch (value)
    {
    case CX_CURVE_TYPE_EXOTIC:          return "X";
    case CX_CURVE_TYPE_FLOW:            return "F";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts CreditCurveType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxCreditCurveTypeFromString (const char* str, CxTCreditCurveType *val)
{
    static char routine[] = "CxCreditCurveTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CxTCreditCurveType)(CX_CURVE_TYPE_FLOW);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'F':
        if (strcmp (buf, "F") == 0) *val = CX_CURVE_TYPE_FLOW;
        else goto done;
        break;
    case 'X':
        if (strcmp (buf, "X") == 0) *val = CX_CURVE_TYPE_EXOTIC;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Converts RecoveryInterp to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxRecoveryInterpToString (CxTRecoveryInterp value)
{
    static char routine[] = "CxRecoveryInterpToString";

    switch (value)
    {
    case CX_RECOVERY_LINEAR_INTERP:     return "Linear";
    case CX_RECOVERY_RIGHT_INTERP:      return "Right";
    case CX_RECOVERY_LEFT_INTERP:       return "Left";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts RecoveryInterp from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxRecoveryInterpFromString (const char* str, CxTRecoveryInterp *val)
{
    static char routine[] = "CxRecoveryInterpFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[7];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CxTRecoveryInterp)(CX_RECOVERY_LINEAR_INTERP);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'L':
        if (strcmp (buf, "LEFT") == 0) *val = CX_RECOVERY_LEFT_INTERP;
        else if (strcmp (buf, "LINEAR") == 0) *val = CX_RECOVERY_LINEAR_INTERP;
        else goto done;
        break;
    case 'R':
        if (strcmp (buf, "RIGHT") == 0) *val = CX_RECOVERY_RIGHT_INTERP;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Converts ProtPayConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxProtPayConvToString (CxTProtPayConv value)
{
    static char routine[] = "CxProtPayConvToString";

    switch (value)
    {
    case CX_PROT_PAY_DEF:               return "Default";
    case CX_PROT_PAY_MAT:               return "Maturity";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts ProtPayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxProtPayConvFromString (const char* str, CxTProtPayConv *val)
{
    static char routine[] = "CxProtPayConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CxTProtPayConv)(CX_PROT_PAY_DEF);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'D':
        if (strcmp (buf, "D") == 0) *val = CX_PROT_PAY_DEF;
        else goto done;
        break;
    case 'M':
        if (strcmp (buf, "M") == 0) *val = CX_PROT_PAY_MAT;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Converts AccrualPayConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CxAccrualPayConvToString (CxTAccrualPayConv value)
{
    static char routine[] = "CxAccrualPayConvToString";

    switch (value)
    {
    case CX_ACCRUAL_PAY_NONE:           return "None";
    case CX_ACCRUAL_PAY_ALL:            return "All";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts AccrualPayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CxAccrualPayConvFromString (const char* str, CxTAccrualPayConv *val)
{
    static char routine[] = "CxAccrualPayConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CxTAccrualPayConv)(CX_ACCRUAL_PAY_ALL);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'A':
        if (strcmp (buf, "A") == 0) *val = CX_ACCRUAL_PAY_ALL;
        else goto done;
        break;
    case 'N':
        if (strcmp (buf, "N") == 0) *val = CX_ACCRUAL_PAY_NONE;
        else goto done;
        break;
    default:
        goto done;
    }
    status = SUCCESS;

  done:

    if (status != SUCCESS && str != NULL && val != NULL)
        GtoErrMsg("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Constructor for CxTCreditCurve
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveMake(
CxTCreditCurveType type,                /* (I) */
TCurve*         tc,                  /* (I) */
TDateInterval*  timestep             /* (I) */
)
{
    static char routine[] = "CxCreditCurveMake";
    int status = FAILURE;

    CxTCreditCurve* p = NULL;

    REQUIRE(tc != NULL);

    p = NEW(CxTCreditCurve);
    if (p==NULL) goto done;

    p->type            = type;
    p->tc              = GtoCopyCurve(tc);
    if (p->tc == NULL) goto done;

    if (timestep != NULL)
    {
        p->timestep = GtoDateIntervalMakeCopy(timestep);
        if (p->timestep == NULL) goto done;
    }
    else
    {
        p->timestep = NULL;
    }


    if (CxCreditCurveValidate(p) != SUCCESS) goto done; /* failure */

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCreditCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCreditCurve
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveMakeEmpty(void)
{
    static char routine[] = "CxCreditCurveMakeEmpty";
    int status = FAILURE;

    CxTCreditCurve* p = NULL;


    p = NEW(CxTCreditCurve);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCreditCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCreditCurve
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveCopy(CxTCreditCurve* src)
{
    CxTCreditCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCreditCurveMake(src->type,
                            src->tc,
                            src->timestep);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCreditCurve
***************************************************************************
*/
void CxCreditCurveFree(CxTCreditCurve *p)
{
    if (p != NULL)
    {
        GtoFreeTCurve(p->tc);
        GtoDateIntervalDelete(p->timestep);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTIrMarketData
***************************************************************************
*/
CxTIrMarketData* CxIrMarketDataMake(
char*           ccy,                 /* (I) */
TDate           valueDate,           /* (I) */
TCurve*         discCurve,           /* (I) */
TCurve*         liborCurve           /* (I) */
)
{
    static char routine[] = "CxIrMarketDataMake";
    int status = FAILURE;

    CxTIrMarketData* p = NULL;

    REQUIRE(ccy != NULL);
    REQUIRE(discCurve != NULL);

    p = NEW(CxTIrMarketData);
    if (p==NULL) goto done;

    p->ccy             = STRDUP(ccy);
    if (p->ccy == NULL) goto done;

    p->valueDate       = valueDate;
    p->discCurve       = GtoCopyCurve(discCurve);
    if (p->discCurve == NULL) goto done;

    if (liborCurve != NULL)
    {
        p->liborCurve = GtoCopyCurve(liborCurve);
        if (p->liborCurve == NULL) goto done;
    }
    else
    {
        p->liborCurve = NULL;
    }


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxIrMarketDataFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTIrMarketData
***************************************************************************
*/
CxTIrMarketData* CxIrMarketDataMakeEmpty(void)
{
    static char routine[] = "CxIrMarketDataMakeEmpty";
    int status = FAILURE;

    CxTIrMarketData* p = NULL;


    p = NEW(CxTIrMarketData);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxIrMarketDataFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTIrMarketData
***************************************************************************
*/
CxTIrMarketData* CxIrMarketDataCopy(CxTIrMarketData* src)
{
    CxTIrMarketData* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxIrMarketDataMake(src->ccy,
                             src->valueDate,
                             src->discCurve,
                             src->liborCurve);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTIrMarketData
***************************************************************************
*/
void CxIrMarketDataFree(CxTIrMarketData *p)
{
    if (p != NULL)
    {
        FREE(p->ccy);
        GtoFreeTCurve(p->discCurve);
        GtoFreeTCurve(p->liborCurve);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTRecoveryCurve
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveMake(
int             numItems,            /* (I) */
TDate*          dates,               /* (I) [numItems] */
double*         recoveryRates,       /* (I) [numItems] */
CxTRecoveryInterp interpType,          /* (I) */
TBoolean        extrapBefore,        /* (I) */
TBoolean        extrapAfter          /* (I) */
)
{
    static char routine[] = "CxRecoveryCurveMake";
    int status = FAILURE;

    CxTRecoveryCurve* p = NULL;

    REQUIRE(numItems > 0);
    REQUIRE(dates != NULL);
    REQUIRE(recoveryRates != NULL);

    p = NEW(CxTRecoveryCurve);
    if (p==NULL) goto done;

    p->numItems        = numItems;
    p->dates = NEW_ARRAY(TDate, p->numItems);
    if (p->dates == NULL) goto done;
    COPY_ARRAY (p->dates, dates, TDate, p->numItems);

    p->recoveryRates = NEW_ARRAY(double, p->numItems);
    if (p->recoveryRates == NULL) goto done;
    COPY_ARRAY (p->recoveryRates, recoveryRates, double, p->numItems);

    p->interpType      = interpType;
    p->extrapBefore    = extrapBefore;
    p->extrapAfter     = extrapAfter;

    if (CxRecoveryCurveValidate(p) != SUCCESS) goto done; /* failure */

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxRecoveryCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTRecoveryCurve
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveMakeEmpty(
int             numItems             /* (I) */
)
{
    static char routine[] = "CxRecoveryCurveMakeEmpty";
    int status = FAILURE;

    CxTRecoveryCurve* p = NULL;

    REQUIRE(numItems > 0);

    p = NEW(CxTRecoveryCurve);
    if (p==NULL) goto done;

    p->numItems        = numItems;

    p->dates = NEW_ARRAY(TDate, p->numItems);
    if (p->dates == NULL) goto done;

    p->recoveryRates = NEW_ARRAY(double, p->numItems);
    if (p->recoveryRates == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxRecoveryCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTRecoveryCurve
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveCopy(CxTRecoveryCurve* src)
{
    CxTRecoveryCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxRecoveryCurveMake(src->numItems,
                              src->dates,
                              src->recoveryRates,
                              src->interpType,
                              src->extrapBefore,
                              src->extrapAfter);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTRecoveryCurve
***************************************************************************
*/
void CxRecoveryCurveFree(CxTRecoveryCurve *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->recoveryRates);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTCreditMarketData
***************************************************************************
*/
CxTCreditMarketData* CxCreditMarketDataMake(
char*           id,                  /* (I) */
TDate           defaultDate,         /* (I) */
CxTCreditCurve* cleanSpreadCurve,    /* (I) */
CxTRecoveryCurve* recoveryCurve        /* (I) */
)
{
    static char routine[] = "CxCreditMarketDataMake";
    int status = FAILURE;

    CxTCreditMarketData* p = NULL;

    REQUIRE(id != NULL);
    REQUIRE(cleanSpreadCurve != NULL);
    REQUIRE(recoveryCurve != NULL);

    p = NEW(CxTCreditMarketData);
    if (p==NULL) goto done;

    p->id              = STRDUP(id);
    if (p->id == NULL) goto done;

    p->defaultDate     = defaultDate;
    p->cleanSpreadCurve = CxCreditCurveCopy(cleanSpreadCurve);
    if (p->cleanSpreadCurve == NULL) goto done;

    p->recoveryCurve   = CxRecoveryCurveCopy(recoveryCurve);
    if (p->recoveryCurve == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCreditMarketDataFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCreditMarketData
***************************************************************************
*/
CxTCreditMarketData* CxCreditMarketDataMakeEmpty(void)
{
    static char routine[] = "CxCreditMarketDataMakeEmpty";
    int status = FAILURE;

    CxTCreditMarketData* p = NULL;


    p = NEW(CxTCreditMarketData);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCreditMarketDataFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCreditMarketData
***************************************************************************
*/
CxTCreditMarketData* CxCreditMarketDataCopy(CxTCreditMarketData* src)
{
    CxTCreditMarketData* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCreditMarketDataMake(src->id,
                                 src->defaultDate,
                                 src->cleanSpreadCurve,
                                 src->recoveryCurve);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCreditMarketData
***************************************************************************
*/
void CxCreditMarketDataFree(CxTCreditMarketData *p)
{
    if (p != NULL)
    {
        FREE(p->id);
        CxCreditCurveFree(p->cleanSpreadCurve);
        CxRecoveryCurveFree(p->recoveryCurve);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTContingentLeg
***************************************************************************
*/
CxTContingentLeg* CxContingentLegMake(
TDate           startDate,           /* (I) */
int             nbDates,             /* (I) */
TDate*          dates,               /* (I) [nbDates] */
double*         notionals,           /* (I) [nbDates] */
CxTProtPayConv  payType,             /* (I) */
long            payDelay             /* (I) */
)
{
    static char routine[] = "CxContingentLegMake";
    int status = FAILURE;

    CxTContingentLeg* p = NULL;


    p = NEW(CxTContingentLeg);
    if (p==NULL) goto done;

    p->startDate       = startDate;
    p->nbDates         = nbDates;
    if (p->nbDates > 0 && dates != NULL)
    {
        p->dates = NEW_ARRAY(TDate, p->nbDates);
        if (p->dates == NULL) goto done;
        COPY_ARRAY (p->dates, dates, TDate, p->nbDates);
    }
    else
    {
        p->dates = NULL;
    }

    if (p->nbDates > 0 && notionals != NULL)
    {
        p->notionals = NEW_ARRAY(double, p->nbDates);
        if (p->notionals == NULL) goto done;
        COPY_ARRAY (p->notionals, notionals, double, p->nbDates);
    }
    else
    {
        p->notionals = NULL;
    }

    p->payType         = payType;
    p->payDelay        = payDelay;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxContingentLegFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTContingentLeg
***************************************************************************
*/
CxTContingentLeg* CxContingentLegMakeEmpty(
int             nbDates              /* (I) */
)
{
    static char routine[] = "CxContingentLegMakeEmpty";
    int status = FAILURE;

    CxTContingentLeg* p = NULL;


    p = NEW(CxTContingentLeg);
    if (p==NULL) goto done;

    p->nbDates         = nbDates;

    if (p->nbDates > 0)
    {
        p->dates = NEW_ARRAY(TDate, p->nbDates);
        if (p->dates == NULL) goto done;
    }
    else
    {
        p->dates = NULL;
    }

    if (p->nbDates > 0)
    {
        p->notionals = NEW_ARRAY(double, p->nbDates);
        if (p->notionals == NULL) goto done;
    }
    else
    {
        p->notionals = NULL;
    }


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxContingentLegFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTContingentLeg
***************************************************************************
*/
CxTContingentLeg* CxContingentLegCopy(CxTContingentLeg* src)
{
    CxTContingentLeg* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxContingentLegMake(src->startDate,
                              src->nbDates,
                              src->dates,
                              src->notionals,
                              src->payType,
                              src->payDelay);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTContingentLeg
***************************************************************************
*/
void CxContingentLegFree(CxTContingentLeg *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->notionals);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTFeeLeg
***************************************************************************
*/
CxTFeeLeg* CxFeeLegMake(
int             nbDates,             /* (I) */
TDate*          accStartDates,       /* (I) [nbDates] */
TDate*          accEndDates,         /* (I) [nbDates] */
TDate*          payDates,            /* (I) [nbDates] */
double*         notionals,           /* (I) [nbDates] */
double*         couponRates,         /* (I) [nbDates] */
CxTDayCountConv dcc,                 /* (I) */
CxTAccrualPayConv accrualPayConv,      /* (I) */
TBoolean        obsStartOfDay        /* (I) */
)
{
    static char routine[] = "CxFeeLegMake";
    int status = FAILURE;

    CxTFeeLeg* p = NULL;


    p = NEW(CxTFeeLeg);
    if (p==NULL) goto done;

    p->nbDates         = nbDates;
    if (p->nbDates > 0 && accStartDates != NULL)
    {
        p->accStartDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accStartDates == NULL) goto done;
        COPY_ARRAY (p->accStartDates, accStartDates, TDate, p->nbDates);
    }
    else
    {
        p->accStartDates = NULL;
    }

    if (p->nbDates > 0 && accEndDates != NULL)
    {
        p->accEndDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accEndDates == NULL) goto done;
        COPY_ARRAY (p->accEndDates, accEndDates, TDate, p->nbDates);
    }
    else
    {
        p->accEndDates = NULL;
    }

    if (p->nbDates > 0 && payDates != NULL)
    {
        p->payDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->payDates == NULL) goto done;
        COPY_ARRAY (p->payDates, payDates, TDate, p->nbDates);
    }
    else
    {
        p->payDates = NULL;
    }

    if (p->nbDates > 0 && notionals != NULL)
    {
        p->notionals = NEW_ARRAY(double, p->nbDates);
        if (p->notionals == NULL) goto done;
        COPY_ARRAY (p->notionals, notionals, double, p->nbDates);
    }
    else
    {
        p->notionals = NULL;
    }

    if (p->nbDates > 0 && couponRates != NULL)
    {
        p->couponRates = NEW_ARRAY(double, p->nbDates);
        if (p->couponRates == NULL) goto done;
        COPY_ARRAY (p->couponRates, couponRates, double, p->nbDates);
    }
    else
    {
        p->couponRates = NULL;
    }

    p->dcc             = dcc;
    p->accrualPayConv  = accrualPayConv;
    p->obsStartOfDay   = obsStartOfDay;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxFeeLegFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTFeeLeg
***************************************************************************
*/
CxTFeeLeg* CxFeeLegMakeEmpty(
int             nbDates              /* (I) */
)
{
    static char routine[] = "CxFeeLegMakeEmpty";
    int status = FAILURE;

    CxTFeeLeg* p = NULL;


    p = NEW(CxTFeeLeg);
    if (p==NULL) goto done;

    p->nbDates         = nbDates;

    if (p->nbDates > 0)
    {
        p->accStartDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accStartDates == NULL) goto done;
    }
    else
    {
        p->accStartDates = NULL;
    }

    if (p->nbDates > 0)
    {
        p->accEndDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->accEndDates == NULL) goto done;
    }
    else
    {
        p->accEndDates = NULL;
    }

    if (p->nbDates > 0)
    {
        p->payDates = NEW_ARRAY(TDate, p->nbDates);
        if (p->payDates == NULL) goto done;
    }
    else
    {
        p->payDates = NULL;
    }

    if (p->nbDates > 0)
    {
        p->notionals = NEW_ARRAY(double, p->nbDates);
        if (p->notionals == NULL) goto done;
    }
    else
    {
        p->notionals = NULL;
    }

    if (p->nbDates > 0)
    {
        p->couponRates = NEW_ARRAY(double, p->nbDates);
        if (p->couponRates == NULL) goto done;
    }
    else
    {
        p->couponRates = NULL;
    }


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxFeeLegFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTFeeLeg
***************************************************************************
*/
CxTFeeLeg* CxFeeLegCopy(CxTFeeLeg* src)
{
    CxTFeeLeg* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxFeeLegMake(src->nbDates,
                       src->accStartDates,
                       src->accEndDates,
                       src->payDates,
                       src->notionals,
                       src->couponRates,
                       src->dcc,
                       src->accrualPayConv,
                       src->obsStartOfDay);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTFeeLeg
***************************************************************************
*/
void CxFeeLegFree(CxTFeeLeg *p)
{
    if (p != NULL)
    {
        FREE(p->accStartDates);
        FREE(p->accEndDates);
        FREE(p->payDates);
        FREE(p->notionals);
        FREE(p->couponRates);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTCdsType
***************************************************************************
*/
CxTCdsType* CxCdsTypeMake(
TDateInterval*  couponInterval,      /* (I) */
CxTDayCountConv paymentDcc,          /* (I) */
CxTStubType     stubType,            /* (I) */
TBoolean        extraDay,            /* (I) */
TBoolean        payAccOnDefault,     /* (I) */
CxTBadDayConv   badDayConv           /* (I) */
)
{
    static char routine[] = "CxCdsTypeMake";
    int status = FAILURE;

    CxTCdsType* p = NULL;


    p = NEW(CxTCdsType);
    if (p==NULL) goto done;

    if (couponInterval != NULL)
    {
        p->couponInterval = GtoDateIntervalMakeCopy(couponInterval);
        if (p->couponInterval == NULL) goto done;
    }
    else
    {
        p->couponInterval = NULL;
    }

    p->paymentDcc      = paymentDcc;
    p->stubType        = stubType;
    p->extraDay        = extraDay;
    p->payAccOnDefault = payAccOnDefault;
    p->badDayConv      = badDayConv;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsTypeFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCdsType
***************************************************************************
*/
CxTCdsType* CxCdsTypeMakeEmpty(void)
{
    static char routine[] = "CxCdsTypeMakeEmpty";
    int status = FAILURE;

    CxTCdsType* p = NULL;


    p = NEW(CxTCdsType);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsTypeFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCdsType
***************************************************************************
*/
CxTCdsType* CxCdsTypeCopy(CxTCdsType* src)
{
    CxTCdsType* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCdsTypeMake(src->couponInterval,
                        src->paymentDcc,
                        src->stubType,
                        src->extraDay,
                        src->payAccOnDefault,
                        src->badDayConv);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCdsType
***************************************************************************
*/
void CxCdsTypeFree(CxTCdsType *p)
{
    if (p != NULL)
    {
        GtoDateIntervalDelete(p->couponInterval);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTCds
***************************************************************************
*/
CxTCds* CxCdsMake(
CxTCdsType*     type,                /* (I) */
TDate           startDate,           /* (I) */
TDate           endDate,             /* (I) */
double          notional             /* (I) */
)
{
    static char routine[] = "CxCdsMake";
    int status = FAILURE;

    CxTCds* p = NULL;

    REQUIRE(type != NULL);

    p = NEW(CxTCds);
    if (p==NULL) goto done;

    p->type            = CxCdsTypeCopy(type);
    if (p->type == NULL) goto done;

    p->startDate       = startDate;
    p->endDate         = endDate;
    p->notional        = notional;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCds
***************************************************************************
*/
CxTCds* CxCdsMakeEmpty(void)
{
    static char routine[] = "CxCdsMakeEmpty";
    int status = FAILURE;

    CxTCds* p = NULL;


    p = NEW(CxTCds);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCds
***************************************************************************
*/
CxTCds* CxCdsCopy(CxTCds* src)
{
    CxTCds* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCdsMake(src->type,
                    src->startDate,
                    src->endDate,
                    src->notional);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCds
***************************************************************************
*/
void CxCdsFree(CxTCds *p)
{
    if (p != NULL)
    {
        CxCdsTypeFree(p->type);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTCdsCurveDateAdj
***************************************************************************
*/
CxTCdsCurveDateAdj* CxCdsCurveDateAdjMake(
long            nbDate,              /* (I) */
TDate*          dates,               /* (I) [nbDate] */
double*         spreads,             /* (I) [nbDate] */
TMatrix2D*      tweakMatrix          /* (I) */
)
{
    static char routine[] = "CxCdsCurveDateAdjMake";
    int status = FAILURE;

    CxTCdsCurveDateAdj* p = NULL;

    REQUIRE(nbDate > 0);
    REQUIRE(dates != NULL);
    REQUIRE(spreads != NULL);

    p = NEW(CxTCdsCurveDateAdj);
    if (p==NULL) goto done;

    p->nbDate          = nbDate;
    p->dates = NEW_ARRAY(TDate, p->nbDate);
    if (p->dates == NULL) goto done;
    COPY_ARRAY (p->dates, dates, TDate, p->nbDate);

    p->spreads = NEW_ARRAY(double, p->nbDate);
    if (p->spreads == NULL) goto done;
    COPY_ARRAY (p->spreads, spreads, double, p->nbDate);

    if (tweakMatrix != NULL)
    {
        p->tweakMatrix = GtoMatrixCopy(tweakMatrix);
        if (p->tweakMatrix == NULL) goto done;
    }
    else
    {
        p->tweakMatrix = NULL;
    }


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsCurveDateAdjFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCdsCurveDateAdj
***************************************************************************
*/
CxTCdsCurveDateAdj* CxCdsCurveDateAdjMakeEmpty(
long            nbDate               /* (I) */
)
{
    static char routine[] = "CxCdsCurveDateAdjMakeEmpty";
    int status = FAILURE;

    CxTCdsCurveDateAdj* p = NULL;

    REQUIRE(nbDate > 0);

    p = NEW(CxTCdsCurveDateAdj);
    if (p==NULL) goto done;

    p->nbDate          = nbDate;

    p->dates = NEW_ARRAY(TDate, p->nbDate);
    if (p->dates == NULL) goto done;

    p->spreads = NEW_ARRAY(double, p->nbDate);
    if (p->spreads == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsCurveDateAdjFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCdsCurveDateAdj
***************************************************************************
*/
CxTCdsCurveDateAdj* CxCdsCurveDateAdjCopy(CxTCdsCurveDateAdj* src)
{
    CxTCdsCurveDateAdj* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCdsCurveDateAdjMake(src->nbDate,
                                src->dates,
                                src->spreads,
                                src->tweakMatrix);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCdsCurveDateAdj
***************************************************************************
*/
void CxCdsCurveDateAdjFree(CxTCdsCurveDateAdj *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->spreads);
        GtoMatrixFree(p->tweakMatrix);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTCdsCarryOutputs
***************************************************************************
*/
CxTCdsCarryOutputs* CxCdsCarryOutputsMake(
double          couponsPV,           /* (I) */
double          horizonAccrual,      /* (I) */
double          carry                /* (I) */
)
{
    static char routine[] = "CxCdsCarryOutputsMake";
    int status = FAILURE;

    CxTCdsCarryOutputs* p = NULL;


    p = NEW(CxTCdsCarryOutputs);
    if (p==NULL) goto done;

    p->couponsPV       = couponsPV;
    p->horizonAccrual  = horizonAccrual;
    p->carry           = carry;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsCarryOutputsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCdsCarryOutputs
***************************************************************************
*/
CxTCdsCarryOutputs* CxCdsCarryOutputsMakeEmpty(void)
{
    static char routine[] = "CxCdsCarryOutputsMakeEmpty";
    int status = FAILURE;

    CxTCdsCarryOutputs* p = NULL;


    p = NEW(CxTCdsCarryOutputs);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsCarryOutputsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCdsCarryOutputs
***************************************************************************
*/
CxTCdsCarryOutputs* CxCdsCarryOutputsCopy(CxTCdsCarryOutputs* src)
{
    CxTCdsCarryOutputs* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCdsCarryOutputsMake(src->couponsPV,
                                src->horizonAccrual,
                                src->carry);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCdsCarryOutputs
***************************************************************************
*/
void CxCdsCarryOutputsFree(CxTCdsCarryOutputs *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CxTCdsConventions
***************************************************************************
*/
CxTCdsConventions* CxCdsConventionsMake(
TBoolean        payAccOnDefault,     /* (I) */
TDateInterval*  couponInterval,      /* (I) */
CxTDayCountConv paymentDCC,          /* (I) */
CxTStubType     stubType,            /* (I) */
CxTCreditCurveType curveType,           /* (I) */
TDateInterval*  timestep,            /* (I) */
TDateInterval*  smoothInterval,      /* (I) */
TBoolean        protectStart,        /* (I) */
TBoolean        isPriceClean,        /* (I) */
long            delay,               /* (I) */
CxTBadDayConv   badDayConv           /* (I) */
)
{
    static char routine[] = "CxCdsConventionsMake";
    int status = FAILURE;

    CxTCdsConventions* p = NULL;


    p = NEW(CxTCdsConventions);
    if (p==NULL) goto done;

    p->payAccOnDefault = payAccOnDefault;
    if (couponInterval != NULL)
    {
        p->couponInterval = GtoDateIntervalMakeCopy(couponInterval);
        if (p->couponInterval == NULL) goto done;
    }
    else
    {
        p->couponInterval = NULL;
    }

    p->paymentDCC      = paymentDCC;
    p->stubType        = stubType;
    p->curveType       = curveType;
    if (timestep != NULL)
    {
        p->timestep = GtoDateIntervalMakeCopy(timestep);
        if (p->timestep == NULL) goto done;
    }
    else
    {
        p->timestep = NULL;
    }

    if (smoothInterval != NULL)
    {
        p->smoothInterval = GtoDateIntervalMakeCopy(smoothInterval);
        if (p->smoothInterval == NULL) goto done;
    }
    else
    {
        p->smoothInterval = NULL;
    }

    p->protectStart    = protectStart;
    p->isPriceClean    = isPriceClean;
    p->delay           = delay;
    p->badDayConv      = badDayConv;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsConventionsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CxTCdsConventions
***************************************************************************
*/
CxTCdsConventions* CxCdsConventionsMakeEmpty(void)
{
    static char routine[] = "CxCdsConventionsMakeEmpty";
    int status = FAILURE;

    CxTCdsConventions* p = NULL;


    p = NEW(CxTCdsConventions);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CxCdsConventionsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CxTCdsConventions
***************************************************************************
*/
CxTCdsConventions* CxCdsConventionsCopy(CxTCdsConventions* src)
{
    CxTCdsConventions* dst = NULL;
    if (src==NULL) return NULL;

    dst = CxCdsConventionsMake(src->payAccOnDefault,
                               src->couponInterval,
                               src->paymentDCC,
                               src->stubType,
                               src->curveType,
                               src->timestep,
                               src->smoothInterval,
                               src->protectStart,
                               src->isPriceClean,
                               src->delay,
                               src->badDayConv);

    return dst;
}

/**
***************************************************************************
** Destructor for CxTCdsConventions
***************************************************************************
*/
void CxCdsConventionsFree(CxTCdsConventions *p)
{
    if (p != NULL)
    {
        GtoDateIntervalDelete(p->couponInterval);
        GtoDateIntervalDelete(p->timestep);
        GtoDateIntervalDelete(p->smoothInterval);
        FREE(p);
    }
}

