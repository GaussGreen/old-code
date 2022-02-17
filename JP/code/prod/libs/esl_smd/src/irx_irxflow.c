/*
***************************************************************************
** SOURCE FILE: irxflow.c
**
** Defines data structures and enums used in the irxflow library.
***************************************************************************
*/

#include "irx/irxflow.h"
#include "irx/mktconv.h"
#include "irx/rate.h"
#include "irx/swap.h"

#include <ctype.h>
#include <irx/macros.h>
#include <irx/dateutils.h>

/**
***************************************************************************
** Converts BootstrapMethod to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxBootstrapMethodToString (IrxTBootstrapMethod value)
{
    static char routine[] = "irxBootstrapMethodToString";

    switch (value)
    {
    case IRX_BOOTSTRAP_FLAT:            return "Flat";
    case IRX_BOOTSTRAP_SMOOTH:          return "Smooth";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts BootstrapMethod from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxBootstrapMethodFromString (const char* str, IrxTBootstrapMethod *val)
{
    static char routine[] = "irxBootstrapMethodFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (IrxTBootstrapMethod)(IRX_BOOTSTRAP_FLAT);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'F':
        if (strcmp (buf, "F") == 0) *val = IRX_BOOTSTRAP_FLAT;
        else goto RETURN;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = IRX_BOOTSTRAP_SMOOTH;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Converts RateTypeType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxRateTypeTypeToString (IrxTRateTypeType value)
{
    static char routine[] = "irxRateTypeTypeToString";

    switch (value)
    {
    case IRX_RT_SIMPLE:                 return "Simple";
    case IRX_RT_DISCOUNT:               return "Discount";
    case IRX_RT_CONTINUOUS:             return "Continuous";
    case IRX_RT_FREQUENCY:              return "Frequency";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts RateTypeType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxRateTypeTypeFromString (const char* str, IrxTRateTypeType *val)
{
    static char routine[] = "irxRateTypeTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'C':
        if (strcmp (buf, "C") == 0) *val = IRX_RT_CONTINUOUS;
        else goto RETURN;
        break;
    case 'D':
        if (strcmp (buf, "D") == 0) *val = IRX_RT_DISCOUNT;
        else goto RETURN;
        break;
    case 'F':
        if (strcmp (buf, "F") == 0) *val = IRX_RT_FREQUENCY;
        else goto RETURN;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = IRX_RT_SIMPLE;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Converts StubPayment to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxStubPaymentToString (IrxTStubPayment value)
{
    static char routine[] = "irxStubPaymentToString";

    switch (value)
    {
    case IRX_STUB_SIMPLE:               return "Simple";
    case IRX_STUB_BOND:                 return "Bond";
    case IRX_STUB_NONE:                 return "None";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts StubPayment from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxStubPaymentFromString (const char* str, IrxTStubPayment *val)
{
    static char routine[] = "irxStubPaymentFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'B':
        if (strcmp (buf, "B") == 0) *val = IRX_STUB_BOND;
        else goto RETURN;
        break;
    case 'N':
        if (strcmp (buf, "N") == 0) *val = IRX_STUB_NONE;
        else goto RETURN;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = IRX_STUB_SIMPLE;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Constructor for IrxTCashFlowList
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMake(
int             numItems,            /* (I) */
IrxTDate const* dates,               /* (I) [numItems] */
double const*   amounts              /* (I) [numItems] */
)
{
    static char routine[] = "irxCashFlowListMake";
    int status = FAILURE;

    IrxTCashFlowList* p = NULL;

    REQUIRE(numItems > 0);
    REQUIRE(dates != NULL);
    REQUIRE(amounts != NULL);

    p = NEW(IrxTCashFlowList);
    if (p==NULL) goto RETURN;

    p->numItems        = numItems;
    p->dates = NEW_ARRAY(IrxTDate, p->numItems);
    if (p->dates == NULL) goto RETURN;
    COPY_ARRAY (p->dates, dates, IrxTDate, p->numItems);

    p->amounts = NEW_ARRAY(double, p->numItems);
    if (p->amounts == NULL) goto RETURN;
    COPY_ARRAY (p->amounts, amounts, double, p->numItems);


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxCashFlowListFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTCashFlowList
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMakeEmpty(
int             numItems             /* (I) */
)
{
    static char routine[] = "irxCashFlowListMakeEmpty";
    int status = FAILURE;

    IrxTCashFlowList* p = NULL;

    REQUIRE(numItems > 0);

    p = NEW(IrxTCashFlowList);
    if (p==NULL) goto RETURN;

    p->numItems        = numItems;

    p->dates = NEW_ARRAY(IrxTDate, p->numItems);
    if (p->dates == NULL) goto RETURN;

    p->amounts = NEW_ARRAY(double, p->numItems);
    if (p->amounts == NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxCashFlowListFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTCashFlowList
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListCopy(IrxTCashFlowList const* src)
{
    IrxTCashFlowList* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxCashFlowListMake(src->numItems,
                              src->dates,
                              src->amounts);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTCashFlowList
***************************************************************************
*/
void irxCashFlowListFree(IrxTCashFlowList *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->amounts);
        FREE(p);
    }
}



/**
***************************************************************************
** Set legacy default values for IrxTZeroCurve
***************************************************************************
*/
int IrxTZeroCurveSetDefault(IrxTZeroCurve *zc)
{
    static char routine[] = "irxZeroCurveSetDefault";

    if (zc != NULL)
    {
        zc->Today     = zc->baseDate; /**< Today's date                       */
        zc->SpotDays  = 0;            /**< Spot days                          */
        zc->ValueDate = zc->baseDate; /**< Value date, = baseDate             */

        /* Underlying yield curve conventions */
        zc->SwapFreq = 'S';           /**< Benchmark swap frequency           */
        strcpy(zc->SwapDCC, "360");   /**< Benchmark swap day count convention*/
        zc->MMB = 360;

        return SUCCESS;
    }
    else
    {
        irxError("%s: invalid input zero curve.\n", routine);
        return FAILURE;
    }
}



/**
***************************************************************************
** Set legacy values for IrxTZeroCurve from market
***************************************************************************
*/
int IrxTZeroCurveSetFromMarket(IrxTZeroCurve            *zc,
                               IrxTDate                 today,
                               const IrxTMarketConv     *marketConv)
{
    static char routine[] = "irxZeroCurveSetFromMarket";
    int status = FAILURE;

    if (zc != NULL)
    {
        zc->Today     = today;
        zc->SpotDays  = marketConv->daysToSpot;   
        zc->ValueDate = zc->baseDate;

        /* Underlying yield curve conventions */
        if (irxDateIntervalToFreq(marketConv->fixedIvl,
                                  &(zc->SwapFreq)) != SUCCESS) goto RETURN;

        switch (marketConv->fixedDcc)
        {
        case IRX_ACT_365F:
            strcpy(zc->SwapDCC, "365");
            break;
        case IRX_ACT_360:            
            strcpy(zc->SwapDCC, "360");
            break;
        case IRX_B30_360:
        case IRX_B30E_360:
            strcpy(zc->SwapDCC, "ACT"); /* conform to DR wrapper and Fix 3 */
            break;
        default:
            irxError("%s: invalid swap DCC %s.\n", 
                     routine,
                     irxDayCountConvToString(marketConv->fixedDcc));
            goto RETURN;
        }

        switch (marketConv->mmDcc)  /* could be floatDcc? */
        {
        case IRX_ACT_365F:
            zc->MMB = 365;
            break;
        case IRX_ACT_360:            
            zc->MMB = 360;
            break;
        default:
            irxError("%s: invalid MM DCC %s.\n", 
                     routine,
                     irxDayCountConvToString(marketConv->mmDcc));
            goto RETURN;
        }

    }
    else
    {
        irxError("%s: invalid input zero curve.\n", routine);
    }

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
    }

    return status;
}





/**
***************************************************************************
** Constructor for IrxTZeroCurve
***************************************************************************
*/
IrxTZeroCurve* irxZeroCurveMake(
IrxTDate        baseDate,            /* (I) */
IrxTDate        maxDate,             /* (I) */
int             numItems,            /* (I) */
IrxTDate const* startDates,          /* (I) [numItems] */
double const*   prices,              /* (I) [numItems] */
double const*   a0,                  /* (I) [numItems] */
double const*   a1,                  /* (I) [numItems] */
double const*   a2                   /* (I) [numItems] */
)
{
    static char routine[] = "irxZeroCurveMake";
    int status = FAILURE;

    IrxTZeroCurve* p = NULL;

    REQUIRE(numItems > 0);
    REQUIRE(startDates != NULL);
    REQUIRE(prices != NULL);
    REQUIRE(a0 != NULL);
    REQUIRE(a1 != NULL);
    REQUIRE(a2 != NULL);

    p = NEW(IrxTZeroCurve);
    if (p==NULL) goto RETURN;

    p->baseDate        = baseDate;
    p->maxDate         = maxDate;
    p->numItems        = numItems;
    p->startDates = NEW_ARRAY(IrxTDate, p->numItems);
    if (p->startDates == NULL) goto RETURN;
    COPY_ARRAY (p->startDates, startDates, IrxTDate, p->numItems);

    p->prices = NEW_ARRAY(double, p->numItems);
    if (p->prices == NULL) goto RETURN;
    COPY_ARRAY (p->prices, prices, double, p->numItems);

    p->a0 = NEW_ARRAY(double, p->numItems);
    if (p->a0 == NULL) goto RETURN;
    COPY_ARRAY (p->a0, a0, double, p->numItems);

    p->a1 = NEW_ARRAY(double, p->numItems);
    if (p->a1 == NULL) goto RETURN;
    COPY_ARRAY (p->a1, a1, double, p->numItems);

    p->a2 = NEW_ARRAY(double, p->numItems);
    if (p->a2 == NULL) goto RETURN;
    COPY_ARRAY (p->a2, a2, double, p->numItems);


    if (IrxTZeroCurveSetDefault(p) != SUCCESS) goto RETURN;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxZeroCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTZeroCurve
***************************************************************************
*/
int irxZeroCurveConstructEmpty(
IrxTZeroCurve*  p,
int             numItems             /* (I) */
)
{
    static char routine[] = "irxZeroCurveConstructEmpty";
    int status = FAILURE;

    REQUIRE(numItems > 0);

    if (p==NULL) goto RETURN;

    p->numItems        = numItems;

    p->startDates = NEW_ARRAY(IrxTDate, p->numItems);
    if (p->startDates == NULL) goto RETURN;

    p->prices = NEW_ARRAY(double, p->numItems);
    if (p->prices == NULL) goto RETURN;

    p->a0 = NEW_ARRAY(double, p->numItems);
    if (p->a0 == NULL) goto RETURN;

    p->a1 = NEW_ARRAY(double, p->numItems);
    if (p->a1 == NULL) goto RETURN;

    p->a2 = NEW_ARRAY(double, p->numItems);
    if (p->a2 == NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxZeroCurveDestroy(p);
    }
    return status;
}


IrxTZeroCurve* irxZeroCurveMakeEmpty(
int             numItems             /* (I) */
)
{
    static char routine[] = "irxZeroCurveMakeEmpty";

    IrxTZeroCurve* p = NULL;

    p = NEW(IrxTZeroCurve);
    if (p==NULL || irxZeroCurveConstructEmpty(p, numItems) != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxZeroCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTZeroCurve
***************************************************************************
*/
IrxTZeroCurve* irxZeroCurveCopy(IrxTZeroCurve const* src)
{
    IrxTZeroCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxZeroCurveMake(src->baseDate,
                           src->maxDate,
                           src->numItems,
                           src->startDates,
                           src->prices,
                           src->a0,
                           src->a1,
                           src->a2);

    /* For legacy wrapper backward compatibility purpose */
    dst->Today     = src->Today;
    dst->SpotDays  = src->SpotDays;
    dst->ValueDate = src->ValueDate;

    dst->SwapFreq = src->SwapFreq;
    strcpy(dst->SwapDCC, src->SwapDCC);
    dst->MMB = src->MMB;

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTZeroCurve
***************************************************************************
*/
void irxZeroCurveDestroy(IrxTZeroCurve *p)
{
    if (p != NULL)
    {
        FREE(p->startDates);
        FREE(p->prices);
        FREE(p->a0);
        FREE(p->a1);
        FREE(p->a2);
    }
}

/**
***************************************************************************
** Destructor + Delete for IrxTZeroCurve
***************************************************************************
*/
void irxZeroCurveFree(IrxTZeroCurve *p)
{
    if (p != NULL)
    {
        irxZeroCurveDestroy(p);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTMarketConv
***************************************************************************
*/
IrxTMarketConv* irxMarketConvMake(
IrxTDateInterval fixedIvl,            /* (I) */
IrxTDateInterval floatIvl,            /* (I) */
IrxTDateInterval cbsIvl,              /* (I) */
IrxTDayCountConv fixedDcc,            /* (I) */
IrxTDayCountConv floatDcc,            /* (I) */
IrxTDayCountConv cbsDcc,              /* (I) */
IrxTDayCountConv mmDcc,               /* (I) */
IrxTBadDayConv  paymentBdc,          /* (I) */
IrxTBadDayConv  accrualBdc,          /* (I) */
IrxTBadDayConv  resetBdc,            /* (I) */
IrxTBadDayConv  mmBdc,               /* (I) */
int             daysToSpot           /* (I) */
)
{
    static char routine[] = "irxMarketConvMake";
    int status = FAILURE;

    IrxTMarketConv* p = NULL;


    p = NEW(IrxTMarketConv);
    if (p==NULL) goto RETURN;

    p->fixedIvl        = fixedIvl;
    p->floatIvl        = floatIvl;
    p->cbsIvl          = cbsIvl;
    p->fixedDcc        = fixedDcc;
    p->floatDcc        = floatDcc;
    p->cbsDcc          = cbsDcc;
    p->mmDcc           = mmDcc;
    p->paymentBdc      = paymentBdc;
    p->accrualBdc      = accrualBdc;
    p->resetBdc        = resetBdc;
    p->mmBdc           = mmBdc;
    p->daysToSpot      = daysToSpot;

    if (irxMarketConvValidate(p) != SUCCESS) goto RETURN; /* failure */

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxMarketConvFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTMarketConv
***************************************************************************
*/
IrxTMarketConv* irxMarketConvMakeEmpty(void)
{
    static char routine[] = "irxMarketConvMakeEmpty";
    int status = FAILURE;

    IrxTMarketConv* p = NULL;


    p = NEW(IrxTMarketConv);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxMarketConvFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTMarketConv
***************************************************************************
*/
IrxTMarketConv* irxMarketConvCopy(IrxTMarketConv const* src)
{
    IrxTMarketConv* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxMarketConvMake(src->fixedIvl,
                            src->floatIvl,
                            src->cbsIvl,
                            src->fixedDcc,
                            src->floatDcc,
                            src->cbsDcc,
                            src->mmDcc,
                            src->paymentBdc,
                            src->accrualBdc,
                            src->resetBdc,
                            src->mmBdc,
                            src->daysToSpot);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTMarketConv
***************************************************************************
*/
void irxMarketConvFree(IrxTMarketConv *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTRateType
***************************************************************************
*/
IrxTRateType* irxRateTypeMake(
IrxTRateTypeType type,                /* (I) */
int             frequency            /* (I) */
)
{
    static char routine[] = "irxRateTypeMake";
    int status = FAILURE;

    IrxTRateType* p = NULL;


    p = NEW(IrxTRateType);
    if (p==NULL) goto RETURN;

    p->type            = type;
    p->frequency       = frequency;

    if (irxRateTypeValidate(p) != SUCCESS) goto RETURN; /* failure */

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxRateTypeFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTRateType
***************************************************************************
*/
IrxTRateType* irxRateTypeMakeEmpty(void)
{
    static char routine[] = "irxRateTypeMakeEmpty";
    int status = FAILURE;

    IrxTRateType* p = NULL;


    p = NEW(IrxTRateType);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxRateTypeFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTRateType
***************************************************************************
*/
IrxTRateType* irxRateTypeCopy(IrxTRateType const* src)
{
    IrxTRateType* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxRateTypeMake(src->type,
                          src->frequency);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTRateType
***************************************************************************
*/
void irxRateTypeFree(IrxTRateType *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTSwapZeroCurve
***************************************************************************
*/
IrxTSwapZeroCurve* irxSwapZeroCurveMake(
IrxTDate        today,               /* (I) */
IrxTMarketConv const* marketConv,          /* (I) */
IrxTCalendar const* calendar,            /* (I) */
IrxTZeroCurve const* discountCurve,       /* (I) */
IrxTZeroCurve const* indexCurve           /* (I) */
)
{
    static char routine[] = "irxSwapZeroCurveMake";
    int status = FAILURE;

    IrxTSwapZeroCurve* p = NULL;

    REQUIRE(marketConv != NULL);
    REQUIRE(calendar != NULL);
    REQUIRE(discountCurve != NULL);

    p = NEW(IrxTSwapZeroCurve);
    if (p==NULL) goto RETURN;

    p->today           = today;
    p->marketConv      = irxMarketConvCopy(marketConv);
    if (p->marketConv == NULL) goto RETURN;

    p->calendar        = irxCalendarCopy(calendar);
    if (p->calendar == NULL) goto RETURN;

    p->discountCurve   = irxZeroCurveCopy(discountCurve);
    if (p->discountCurve == NULL) goto RETURN;

    if (indexCurve != NULL)
    {
        p->indexCurve = irxZeroCurveCopy(indexCurve);
        if (p->indexCurve == NULL) goto RETURN;
    }
    else
    {
        p->indexCurve = NULL;
    }


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxSwapZeroCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTSwapZeroCurve
***************************************************************************
*/
IrxTSwapZeroCurve* irxSwapZeroCurveMakeEmpty(void)
{
    static char routine[] = "irxSwapZeroCurveMakeEmpty";
    int status = FAILURE;

    IrxTSwapZeroCurve* p = NULL;


    p = NEW(IrxTSwapZeroCurve);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxSwapZeroCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTSwapZeroCurve
***************************************************************************
*/
IrxTSwapZeroCurve* irxSwapZeroCurveCopy(IrxTSwapZeroCurve const* src)
{
    IrxTSwapZeroCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxSwapZeroCurveMake(src->today,
                               src->marketConv,
                               src->calendar,
                               src->discountCurve,
                               src->indexCurve);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTSwapZeroCurve
***************************************************************************
*/
void irxSwapZeroCurveFree(IrxTSwapZeroCurve *p)
{
    if (p != NULL)
    {
        irxMarketConvFree(p->marketConv);
        irxCalendarFree(p->calendar);
        irxZeroCurveFree(p->discountCurve);
        irxZeroCurveFree(p->indexCurve);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTSwap
***************************************************************************
*/
IrxTSwap* irxSwapMake(
IrxTMarketConv const* marketConv,          /* (I) */
double          couponRate,          /* (I) */
double          spread,              /* (I) */
IrxTDate        startDate,           /* (I) */
IrxTDate        maturityDate,        /* (I) */
IrxTDate        rollDate,            /* (I) */
IrxTStubLocation stubLocation,        /* (I) */
IrxTStubPayment fixedStubPayment,    /* (I) */
IrxTStubPayment floatStubPayment     /* (I) */
)
{
    static char routine[] = "irxSwapMake";
    int status = FAILURE;

    IrxTSwap* p = NULL;

    REQUIRE(marketConv != NULL);

    p = NEW(IrxTSwap);
    if (p==NULL) goto RETURN;

    p->marketConv      = irxMarketConvCopy(marketConv);
    if (p->marketConv == NULL) goto RETURN;

    p->couponRate      = couponRate;
    p->spread          = spread;
    p->startDate       = startDate;
    p->maturityDate    = maturityDate;
    p->rollDate        = rollDate;
    p->stubLocation    = stubLocation;
    p->fixedStubPayment = fixedStubPayment;
    p->floatStubPayment = floatStubPayment;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxSwapFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTSwap
***************************************************************************
*/
IrxTSwap* irxSwapMakeEmpty(void)
{
    static char routine[] = "irxSwapMakeEmpty";
    int status = FAILURE;

    IrxTSwap* p = NULL;


    p = NEW(IrxTSwap);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxSwapFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTSwap
***************************************************************************
*/
IrxTSwap* irxSwapCopy(IrxTSwap const* src)
{
    IrxTSwap* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxSwapMake(src->marketConv,
                      src->couponRate,
                      src->spread,
                      src->startDate,
                      src->maturityDate,
                      src->rollDate,
                      src->stubLocation,
                      src->fixedStubPayment,
                      src->floatStubPayment);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTSwap
***************************************************************************
*/
void irxSwapFree(IrxTSwap *p)
{
    if (p != NULL)
    {
        irxMarketConvFree(p->marketConv);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTCouponPayment
***************************************************************************
*/
IrxTCouponPayment* irxCouponPaymentMake(
double          notional,            /* (I) */
double          couponRate,          /* (I) */
double          spread,              /* (I) */
int             indexCurveNumber,    /* (I) */
IrxTDayCountConv dcc,                 /* (I) */
IrxTDate        accrueStartDate,     /* (I) */
IrxTDate        accrueEndDate,       /* (I) */
IrxTDate        rateStartDate,       /* (I) */
IrxTDate        rateEndDate          /* (I) */
)
{
    static char routine[] = "irxCouponPaymentMake";
    int status = FAILURE;

    IrxTCouponPayment* p = NULL;


    p = NEW(IrxTCouponPayment);
    if (p==NULL) goto RETURN;

    p->notional        = notional;
    p->couponRate      = couponRate;
    p->spread          = spread;
    p->indexCurveNumber = indexCurveNumber;
    p->dcc             = dcc;
    p->accrueStartDate = accrueStartDate;
    p->accrueEndDate   = accrueEndDate;
    p->rateStartDate   = rateStartDate;
    p->rateEndDate     = rateEndDate;

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxCouponPaymentFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTCouponPayment
***************************************************************************
*/
IrxTCouponPayment* irxCouponPaymentMakeEmpty(void)
{
    static char routine[] = "irxCouponPaymentMakeEmpty";
    int status = FAILURE;

    IrxTCouponPayment* p = NULL;


    p = NEW(IrxTCouponPayment);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxCouponPaymentFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTCouponPayment
***************************************************************************
*/
IrxTCouponPayment* irxCouponPaymentCopy(IrxTCouponPayment const* src)
{
    IrxTCouponPayment* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxCouponPaymentMake(src->notional,
                               src->couponRate,
                               src->spread,
                               src->indexCurveNumber,
                               src->dcc,
                               src->accrueStartDate,
                               src->accrueEndDate,
                               src->rateStartDate,
                               src->rateEndDate);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTCouponPayment
***************************************************************************
*/
void irxCouponPaymentFree(IrxTCouponPayment *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTSwapPayments
***************************************************************************
*/
IrxTSwapPayments* irxSwapPaymentsMake(
int             numFlows,            /* (I) */
IrxTDate const* dates,               /* (I) [numFlows] */
double const*   amounts,             /* (I) [numFlows] */
IrxTCouponPayment const* couponPayments       /* (I) [numFlows] */
)
{
    static char routine[] = "irxSwapPaymentsMake";
    int status = FAILURE;

    IrxTSwapPayments* p = NULL;

    REQUIRE(numFlows > 0);
    REQUIRE(dates != NULL);
    REQUIRE(amounts != NULL);
    REQUIRE(couponPayments != NULL);

    p = NEW(IrxTSwapPayments);
    if (p==NULL) goto RETURN;

    p->numFlows        = numFlows;
    p->dates = NEW_ARRAY(IrxTDate, p->numFlows);
    if (p->dates == NULL) goto RETURN;
    COPY_ARRAY (p->dates, dates, IrxTDate, p->numFlows);

    p->amounts = NEW_ARRAY(double, p->numFlows);
    if (p->amounts == NULL) goto RETURN;
    COPY_ARRAY (p->amounts, amounts, double, p->numFlows);

    p->couponPayments = NEW_ARRAY(IrxTCouponPayment, p->numFlows);
    if (p->couponPayments == NULL) goto RETURN;
    COPY_ARRAY (p->couponPayments, couponPayments, IrxTCouponPayment, p->numFlows);


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxSwapPaymentsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTSwapPayments
***************************************************************************
*/
IrxTSwapPayments* irxSwapPaymentsMakeEmpty(
int             numFlows             /* (I) */
)
{
    static char routine[] = "irxSwapPaymentsMakeEmpty";
    int status = FAILURE;

    IrxTSwapPayments* p = NULL;

    REQUIRE(numFlows > 0);

    p = NEW(IrxTSwapPayments);
    if (p==NULL) goto RETURN;

    p->numFlows        = numFlows;

    p->dates = NEW_ARRAY(IrxTDate, p->numFlows);
    if (p->dates == NULL) goto RETURN;

    p->amounts = NEW_ARRAY(double, p->numFlows);
    if (p->amounts == NULL) goto RETURN;

    p->couponPayments = NEW_ARRAY(IrxTCouponPayment, p->numFlows);
    if (p->couponPayments == NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxSwapPaymentsFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTSwapPayments
***************************************************************************
*/
IrxTSwapPayments* irxSwapPaymentsCopy(IrxTSwapPayments const* src)
{
    IrxTSwapPayments* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxSwapPaymentsMake(src->numFlows,
                              src->dates,
                              src->amounts,
                              src->couponPayments);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTSwapPayments
***************************************************************************
*/
void irxSwapPaymentsFree(IrxTSwapPayments *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->amounts);
        FREE(p->couponPayments);
        FREE(p);
    }
}

