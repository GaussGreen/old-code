/*
***************************************************************************
** SOURCE FILE: crxdata.c
**
** Defines data structures for crxflow library.
***************************************************************************
*/

#include "crxdata.h"
#include "crxmacros.h"

#include <ctype.h>
#include <alib/bondcnst.h>
#include <alib/convert.h>
#include <alib/dtivlo.h>

/**
***************************************************************************
** Converts CdsOptionPayoff to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxCdsOptionPayoffToString (CrxTCdsOptionPayoff value)
{
    static char routine[] = "CrxCdsOptionPayoffToString";

    switch (value)
    {
    case CRX_CDS_OPTION_PAYOFF_MARKET:  return "Market";
    case CRX_CDS_OPTION_PAYOFF_STRIKE:  return "Strike";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts CdsOptionPayoff from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxCdsOptionPayoffFromString (const char* str, CrxTCdsOptionPayoff *val)
{
    static char routine[] = "CrxCdsOptionPayoffFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CrxTCdsOptionPayoff)(CRX_CDS_OPTION_PAYOFF_MARKET);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'M':
        if (strcmp (buf, "M") == 0) *val = CRX_CDS_OPTION_PAYOFF_MARKET;
        else goto done;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = CRX_CDS_OPTION_PAYOFF_STRIKE;
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
** Converts CdsOptionStrikeType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxCdsOptionStrikeTypeToString (CrxTCdsOptionStrikeType value)
{
    static char routine[] = "CrxCdsOptionStrikeTypeToString";

    switch (value)
    {
    case CRX_CDS_OPTION_STRIKE_TYPE_SPREAD: return "Spread";
    case CRX_CDS_OPTION_STRIKE_TYPE_PRICE: return "Price";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts CdsOptionStrikeType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxCdsOptionStrikeTypeFromString (const char* str, CrxTCdsOptionStrikeType *val)
{
    static char routine[] = "CrxCdsOptionStrikeTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        *val = (CrxTCdsOptionStrikeType)(CRX_CDS_OPTION_STRIKE_TYPE_SPREAD);
        return SUCCESS;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'P':
        if (strcmp (buf, "P") == 0) *val = CRX_CDS_OPTION_STRIKE_TYPE_PRICE;
        else goto done;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = CRX_CDS_OPTION_STRIKE_TYPE_SPREAD;
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
** Converts DistType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxDistTypeToString (long value)
{
    static char routine[] = "CrxDistTypeToString";

    switch (value)
    {
    case CRX_DIST_TYPE_LOGNORMAL:       return "Lognormal";
    case CRX_DIST_TYPE_NORMAL:          return "Normal";
    case CRX_DIST_TYPE_Q:               return "Q";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts DistType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxDistTypeFromString (const char* str, long *val)
{
    static char routine[] = "CrxDistTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'L':
        if (strcmp (buf, "L") == 0) *val = CRX_DIST_TYPE_LOGNORMAL;
        else goto done;
        break;
    case 'N':
        if (strcmp (buf, "N") == 0) *val = CRX_DIST_TYPE_NORMAL;
        else goto done;
        break;
    case 'Q':
        if (strcmp (buf, "Q") == 0) *val = CRX_DIST_TYPE_Q;
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
** Converts OptionType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxOptionTypeToString (long value)
{
    static char routine[] = "CrxOptionTypeToString";

    switch (value)
    {
    case CRX_OPTION_TYPE_CALL:          return "Call";
    case CRX_OPTION_TYPE_PUT:           return "Put";
    case CRX_OPTION_TYPE_STRADDLE:      return "Straddle";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts OptionType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxOptionTypeFromString (const char* str, long *val)
{
    static char routine[] = "CrxOptionTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'C':
        if (strcmp (buf, "C") == 0) *val = CRX_OPTION_TYPE_CALL;
        else goto done;
        break;
    case 'P':
        if (strcmp (buf, "P") == 0) *val = CRX_OPTION_TYPE_PUT;
        else goto done;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = CRX_OPTION_TYPE_STRADDLE;
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
** Converts DayCountConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* CrxDayCountConvToString (long value)
{
    static char routine[] = "CrxDayCountConvToString";

    switch (value)
    {
    case CRX_ACT_ACT:                   return "ACT/ACT";
    case CRX_ACT_365F:                  return "ACT/365F";
    case CRX_ACT_360:                   return "ACT/360";
    case CRX_B30_360:                   return "30/360";
    case CRX_B30E_360:                  return "30E/360";
    }

    GtoErrMsg("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts DayCountConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int CrxDayCountConvFromString (const char* str, long *val)
{
    static char routine[] = "CrxDayCountConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[9];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        GtoErrMsg ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case '3':
        if (strcmp (buf, "30/360") == 0) *val = CRX_B30_360;
        else if (strcmp (buf, "30E/360") == 0) *val = CRX_B30E_360;
        else goto done;
        break;
    case 'A':
        if (strcmp (buf, "ACT/360") == 0) *val = CRX_ACT_360;
        else if (strcmp (buf, "ACT/365F") == 0) *val = CRX_ACT_365F;
        else if (strcmp (buf, "ACT/ACT") == 0) *val = CRX_ACT_ACT;
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
** Constructor for CrxTQDist
***************************************************************************
*/
CrxTQDist* CrxQDistMake(
int             nQs,                 /* (I) */
double*         Qs,                  /* (I) [nQs] */
int             nDs,                 /* (I) */
double*         Ds                   /* (I) [nDs] */
)
{
    static char routine[] = "CrxQDistMake";
    int status = FAILURE;

    CrxTQDist* p = NULL;

    REQUIRE(nQs > 0);
    REQUIRE(Qs != NULL);
    REQUIRE(nDs > 0);
    REQUIRE(Ds != NULL);

    p = NEW(CrxTQDist);
    if (p==NULL) goto done;

    p->nQs             = nQs;
    p->Qs = NEW_ARRAY(double, p->nQs);
    if (p->Qs == NULL) goto done;
    COPY_ARRAY (p->Qs, Qs, double, p->nQs);

    p->nDs             = nDs;
    p->Ds = NEW_ARRAY(double, p->nDs);
    if (p->Ds == NULL) goto done;
    COPY_ARRAY (p->Ds, Ds, double, p->nDs);


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxQDistFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTQDist
***************************************************************************
*/
CrxTQDist* CrxQDistMakeEmpty(
int             nQs,                 /* (I) */
int             nDs                  /* (I) */
)
{
    static char routine[] = "CrxQDistMakeEmpty";
    int status = FAILURE;

    CrxTQDist* p = NULL;

    REQUIRE(nQs > 0);
    REQUIRE(nDs > 0);

    p = NEW(CrxTQDist);
    if (p==NULL) goto done;

    p->nQs             = nQs;

    p->Qs = NEW_ARRAY(double, p->nQs);
    if (p->Qs == NULL) goto done;

    p->nDs             = nDs;

    p->Ds = NEW_ARRAY(double, p->nDs);
    if (p->Ds == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxQDistFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTQDist
***************************************************************************
*/
CrxTQDist* CrxQDistCopy(CrxTQDist* src)
{
    CrxTQDist* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxQDistMake(src->nQs,
                       src->Qs,
                       src->nDs,
                       src->Ds);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTQDist
***************************************************************************
*/
void CrxQDistFree(CrxTQDist *p)
{
    if (p != NULL)
    {
        FREE(p->Qs);
        FREE(p->Ds);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTCdsOption
***************************************************************************
*/
CrxTCdsOption* CrxCdsOptionMake(
TDate           issueDate,           /* (I) */
TDate           maturityDate,        /* (I) */
TDateInterval*  feeInterval,         /* (I) */
long            dcc,                 /* (I) */
TBoolean        payAccOnDefault,     /* (I) */
TBoolean        isIndex,             /* (I) */
double          lossSoFar,           /* (I) */
long            optionType,          /* (I) */
TDate           exerciseDate,        /* (I) */
TDate           paymentDate,         /* (I) */
TBoolean        koOnDefault,         /* (I) */
double          strike,              /* (I) */
CrxTCdsOptionPayoff optionPayoff,        /* (I) */
double          coupon,              /* (I) */
CrxTCdsOptionStrikeType strikeType           /* (I) */
)
{
    static char routine[] = "CrxCdsOptionMake";
    int status = FAILURE;

    CrxTCdsOption* p = NULL;

    REQUIRE(feeInterval != NULL);

    p = NEW(CrxTCdsOption);
    if (p==NULL) goto done;

    p->issueDate       = issueDate;
    p->maturityDate    = maturityDate;
    p->feeInterval     = GtoDateIntervalMakeCopy(feeInterval);
    if (p->feeInterval == NULL) goto done;

    p->dcc             = dcc;
    p->payAccOnDefault = payAccOnDefault;
    p->isIndex         = isIndex;
    p->lossSoFar       = lossSoFar;
    p->optionType      = optionType;
    p->exerciseDate    = exerciseDate;
    p->paymentDate     = paymentDate;
    p->koOnDefault     = koOnDefault;
    p->strike          = strike;
    p->optionPayoff    = optionPayoff;
    p->coupon          = coupon;
    p->strikeType      = strikeType;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTCdsOption
***************************************************************************
*/
CrxTCdsOption* CrxCdsOptionMakeEmpty(void)
{
    static char routine[] = "CrxCdsOptionMakeEmpty";
    int status = FAILURE;

    CrxTCdsOption* p = NULL;


    p = NEW(CrxTCdsOption);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTCdsOption
***************************************************************************
*/
CrxTCdsOption* CrxCdsOptionCopy(CrxTCdsOption* src)
{
    CrxTCdsOption* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxCdsOptionMake(src->issueDate,
                           src->maturityDate,
                           src->feeInterval,
                           src->dcc,
                           src->payAccOnDefault,
                           src->isIndex,
                           src->lossSoFar,
                           src->optionType,
                           src->exerciseDate,
                           src->paymentDate,
                           src->koOnDefault,
                           src->strike,
                           src->optionPayoff,
                           src->coupon,
                           src->strikeType);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTCdsOption
***************************************************************************
*/
void CrxCdsOptionFree(CrxTCdsOption *p)
{
    if (p != NULL)
    {
        GtoDateIntervalDelete(p->feeInterval);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTCdsOptionCalc
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalcMake(
double          price,               /* (I) */
double          vol,                 /* (I) */
double          fwdSpread,           /* (I) */
double          fwdSpreadUnadj,      /* (I) */
double          annuity,             /* (I) */
double          timeToExpiry,        /* (I) */
double          fwdPrice             /* (I) */
)
{
    static char routine[] = "CrxCdsOptionCalcMake";
    int status = FAILURE;

    CrxTCdsOptionCalc* p = NULL;


    p = NEW(CrxTCdsOptionCalc);
    if (p==NULL) goto done;

    p->price           = price;
    p->vol             = vol;
    p->fwdSpread       = fwdSpread;
    p->fwdSpreadUnadj  = fwdSpreadUnadj;
    p->annuity         = annuity;
    p->timeToExpiry    = timeToExpiry;
    p->fwdPrice        = fwdPrice;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionCalcFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTCdsOptionCalc
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalcMakeEmpty(void)
{
    static char routine[] = "CrxCdsOptionCalcMakeEmpty";
    int status = FAILURE;

    CrxTCdsOptionCalc* p = NULL;


    p = NEW(CrxTCdsOptionCalc);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionCalcFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTCdsOptionCalc
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalcCopy(CrxTCdsOptionCalc* src)
{
    CrxTCdsOptionCalc* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxCdsOptionCalcMake(src->price,
                               src->vol,
                               src->fwdSpread,
                               src->fwdSpreadUnadj,
                               src->annuity,
                               src->timeToExpiry,
                               src->fwdPrice);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTCdsOptionCalc
***************************************************************************
*/
void CrxCdsOptionCalcFree(CrxTCdsOptionCalc *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTCdsOptionQOptimization
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationMake(
CrxTQDist*      qdist,               /* (I) */
double          vol,                 /* (I) */
TMatrix2D*      direction,           /* (I) */
int             iter,                /* (I) */
double          value,               /* (I) */
double          vdiff                /* (I) */
)
{
    static char routine[] = "CrxCdsOptionQOptimizationMake";
    int status = FAILURE;

    CrxTCdsOptionQOptimization* p = NULL;

    REQUIRE(qdist != NULL);

    p = NEW(CrxTCdsOptionQOptimization);
    if (p==NULL) goto done;

    p->qdist           = CrxQDistCopy(qdist);
    if (p->qdist == NULL) goto done;

    p->vol             = vol;
    if (direction != NULL)
    {
        p->direction = GtoMatrixCopy(direction);
        if (p->direction == NULL) goto done;
    }
    else
    {
        p->direction = NULL;
    }

    p->iter            = iter;
    p->value           = value;
    p->vdiff           = vdiff;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionQOptimizationFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTCdsOptionQOptimization
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationMakeEmpty(void)
{
    static char routine[] = "CrxCdsOptionQOptimizationMakeEmpty";
    int status = FAILURE;

    CrxTCdsOptionQOptimization* p = NULL;


    p = NEW(CrxTCdsOptionQOptimization);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionQOptimizationFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTCdsOptionQOptimization
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationCopy(CrxTCdsOptionQOptimization* src)
{
    CrxTCdsOptionQOptimization* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxCdsOptionQOptimizationMake(src->qdist,
                                        src->vol,
                                        src->direction,
                                        src->iter,
                                        src->value,
                                        src->vdiff);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTCdsOptionQOptimization
***************************************************************************
*/
void CrxCdsOptionQOptimizationFree(CrxTCdsOptionQOptimization *p)
{
    if (p != NULL)
    {
        CrxQDistFree(p->qdist);
        GtoMatrixFree(p->direction);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTCdsOptionQOptModel
***************************************************************************
*/
CrxTCdsOptionQOptModel* CrxCdsOptionQOptModelMake(
TBoolean        optimizeVol,         /* (I) */
int             numQsNoOpt,          /* (I) */
char**          qsNoOpt,             /* (I) [numQsNoOpt] */
int             maxIter,             /* (I) */
double          tolerance,           /* (I) */
double          maxTime              /* (I) */
)
{
    static char routine[] = "CrxCdsOptionQOptModelMake";
    int status = FAILURE;

    CrxTCdsOptionQOptModel* p = NULL;


    p = NEW(CrxTCdsOptionQOptModel);
    if (p==NULL) goto done;

    p->optimizeVol     = optimizeVol;
    p->numQsNoOpt      = numQsNoOpt;
    if (p->numQsNoOpt > 0 && qsNoOpt != NULL)
    {
        int i;
        p->qsNoOpt = NEW_ARRAY(char*, p->numQsNoOpt);
        if (p->qsNoOpt == NULL) goto done;
        for (i = 0; i < p->numQsNoOpt; ++i)
        {
            if (qsNoOpt[i] != NULL)
            {
                p->qsNoOpt[i] = STRDUP(qsNoOpt[i]);
                if (p->qsNoOpt[i] == NULL) goto done;
            }
        }
    }
    else
    {
        p->qsNoOpt = NULL;
    }

    p->maxIter         = maxIter;
    p->tolerance       = tolerance;
    p->maxTime         = maxTime;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionQOptModelFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTCdsOptionQOptModel
***************************************************************************
*/
CrxTCdsOptionQOptModel* CrxCdsOptionQOptModelMakeEmpty(
int             numQsNoOpt           /* (I) */
)
{
    static char routine[] = "CrxCdsOptionQOptModelMakeEmpty";
    int status = FAILURE;

    CrxTCdsOptionQOptModel* p = NULL;


    p = NEW(CrxTCdsOptionQOptModel);
    if (p==NULL) goto done;

    p->numQsNoOpt      = numQsNoOpt;

    if (p->numQsNoOpt > 0)
    {
        p->qsNoOpt = NEW_ARRAY(char*, p->numQsNoOpt);
        if (p->qsNoOpt == NULL) goto done;
    }
    else
    {
        p->qsNoOpt = NULL;
    }


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxCdsOptionQOptModelFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTCdsOptionQOptModel
***************************************************************************
*/
CrxTCdsOptionQOptModel* CrxCdsOptionQOptModelCopy(CrxTCdsOptionQOptModel* src)
{
    CrxTCdsOptionQOptModel* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxCdsOptionQOptModelMake(src->optimizeVol,
                                    src->numQsNoOpt,
                                    src->qsNoOpt,
                                    src->maxIter,
                                    src->tolerance,
                                    src->maxTime);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTCdsOptionQOptModel
***************************************************************************
*/
void CrxCdsOptionQOptModelFree(CrxTCdsOptionQOptModel *p)
{
    if (p != NULL)
    {
        if(p->qsNoOpt != NULL)
        {
            int i;
            for (i = 0; i < p->numQsNoOpt; ++i)
                FREE(p->qsNoOpt[i]);
            FREE(p->qsNoOpt);
        }
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceMake(
TDate           settleDate,          /* (I) */
double          cleanPrice           /* (I) */
)
{
    static char routine[] = "CrxBondPriceMake";
    int status = FAILURE;

    CrxTBondPrice* p = NULL;


    p = NEW(CrxTBondPrice);
    if (p==NULL) goto done;

    p->settleDate      = settleDate;
    p->cleanPrice      = cleanPrice;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceMakeEmpty(void)
{
    static char routine[] = "CrxBondPriceMakeEmpty";
    int status = FAILURE;

    CrxTBondPrice* p = NULL;


    p = NEW(CrxTBondPrice);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceCopy(CrxTBondPrice* src)
{
    CrxTBondPrice* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondPriceMake(src->settleDate,
                           src->cleanPrice);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondPrice
***************************************************************************
*/
void CrxBondPriceFree(CrxTBondPrice *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveMake(
TDate           spotSettleDate,      /* (I) */
long            dcc,                 /* (I) */
int             numDates,            /* (I) */
TDate*          dates,               /* (I) [numDates] */
double*         rates                /* (I) [numDates] */
)
{
    static char routine[] = "CrxBondRepoCurveMake";
    int status = FAILURE;

    CrxTBondRepoCurve* p = NULL;

    REQUIRE(numDates > 0);
    REQUIRE(dates != NULL);
    REQUIRE(rates != NULL);

    p = NEW(CrxTBondRepoCurve);
    if (p==NULL) goto done;

    p->spotSettleDate  = spotSettleDate;
    p->dcc             = dcc;
    p->numDates        = numDates;
    p->dates = NEW_ARRAY(TDate, p->numDates);
    if (p->dates == NULL) goto done;
    COPY_ARRAY (p->dates, dates, TDate, p->numDates);

    p->rates = NEW_ARRAY(double, p->numDates);
    if (p->rates == NULL) goto done;
    COPY_ARRAY (p->rates, rates, double, p->numDates);


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondRepoCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveMakeEmpty(
int             numDates             /* (I) */
)
{
    static char routine[] = "CrxBondRepoCurveMakeEmpty";
    int status = FAILURE;

    CrxTBondRepoCurve* p = NULL;

    REQUIRE(numDates > 0);

    p = NEW(CrxTBondRepoCurve);
    if (p==NULL) goto done;

    p->numDates        = numDates;

    p->dates = NEW_ARRAY(TDate, p->numDates);
    if (p->dates == NULL) goto done;

    p->rates = NEW_ARRAY(double, p->numDates);
    if (p->rates == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondRepoCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveCopy(CrxTBondRepoCurve* src)
{
    CrxTBondRepoCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondRepoCurveMake(src->spotSettleDate,
                               src->dcc,
                               src->numDates,
                               src->dates,
                               src->rates);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondRepoCurve
***************************************************************************
*/
void CrxBondRepoCurveFree(CrxTBondRepoCurve *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->rates);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveMake(
TDate           today,               /* (I) */
int             numDates,            /* (I) */
TDate*          dates,               /* (I) [numDates] */
double*         vols                 /* (I) [numDates] */
)
{
    static char routine[] = "CrxBondPriceVolCurveMake";
    int status = FAILURE;

    CrxTBondPriceVolCurve* p = NULL;

    REQUIRE(numDates > 0);
    REQUIRE(dates != NULL);
    REQUIRE(vols != NULL);

    p = NEW(CrxTBondPriceVolCurve);
    if (p==NULL) goto done;

    p->today           = today;
    p->numDates        = numDates;
    p->dates = NEW_ARRAY(TDate, p->numDates);
    if (p->dates == NULL) goto done;
    COPY_ARRAY (p->dates, dates, TDate, p->numDates);

    p->vols = NEW_ARRAY(double, p->numDates);
    if (p->vols == NULL) goto done;
    COPY_ARRAY (p->vols, vols, double, p->numDates);


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceVolCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveMakeEmpty(
int             numDates             /* (I) */
)
{
    static char routine[] = "CrxBondPriceVolCurveMakeEmpty";
    int status = FAILURE;

    CrxTBondPriceVolCurve* p = NULL;

    REQUIRE(numDates > 0);

    p = NEW(CrxTBondPriceVolCurve);
    if (p==NULL) goto done;

    p->numDates        = numDates;

    p->dates = NEW_ARRAY(TDate, p->numDates);
    if (p->dates == NULL) goto done;

    p->vols = NEW_ARRAY(double, p->numDates);
    if (p->vols == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceVolCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveCopy(CrxTBondPriceVolCurve* src)
{
    CrxTBondPriceVolCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondPriceVolCurveMake(src->today,
                                   src->numDates,
                                   src->dates,
                                   src->vols);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondPriceVolCurve
***************************************************************************
*/
void CrxBondPriceVolCurveFree(CrxTBondPriceVolCurve *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->vols);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveMake(
TDate           today,               /* (I) */
int             numDates,            /* (I) */
TDate*          dates,               /* (I) [numDates] */
double*         vols                 /* (I) [numDates] */
)
{
    static char routine[] = "CrxBondSpreadVolCurveMake";
    int status = FAILURE;

    CrxTBondSpreadVolCurve* p = NULL;

    REQUIRE(numDates > 0);
    REQUIRE(dates != NULL);
    REQUIRE(vols != NULL);

    p = NEW(CrxTBondSpreadVolCurve);
    if (p==NULL) goto done;

    p->today           = today;
    p->numDates        = numDates;
    p->dates = NEW_ARRAY(TDate, p->numDates);
    if (p->dates == NULL) goto done;
    COPY_ARRAY (p->dates, dates, TDate, p->numDates);

    p->vols = NEW_ARRAY(double, p->numDates);
    if (p->vols == NULL) goto done;
    COPY_ARRAY (p->vols, vols, double, p->numDates);


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondSpreadVolCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveMakeEmpty(
int             numDates             /* (I) */
)
{
    static char routine[] = "CrxBondSpreadVolCurveMakeEmpty";
    int status = FAILURE;

    CrxTBondSpreadVolCurve* p = NULL;

    REQUIRE(numDates > 0);

    p = NEW(CrxTBondSpreadVolCurve);
    if (p==NULL) goto done;

    p->numDates        = numDates;

    p->dates = NEW_ARRAY(TDate, p->numDates);
    if (p->dates == NULL) goto done;

    p->vols = NEW_ARRAY(double, p->numDates);
    if (p->vols == NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondSpreadVolCurveFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveCopy(CrxTBondSpreadVolCurve* src)
{
    CrxTBondSpreadVolCurve* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondSpreadVolCurveMake(src->today,
                                    src->numDates,
                                    src->dates,
                                    src->vols);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondSpreadVolCurve
***************************************************************************
*/
void CrxBondSpreadVolCurveFree(CrxTBondSpreadVolCurve *p)
{
    if (p != NULL)
    {
        FREE(p->dates);
        FREE(p->vols);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondPriceOption
***************************************************************************
*/
CrxTBondPriceOption* CrxBondPriceOptionMake(
TBond*          bond,                /* (I) */
long            optionType,          /* (I) */
TDate           exerciseDate,        /* (I) */
TDate           paymentDate,         /* (I) */
double          strikePrice          /* (I) */
)
{
    static char routine[] = "CrxBondPriceOptionMake";
    int status = FAILURE;

    CrxTBondPriceOption* p = NULL;

    REQUIRE(bond != NULL);

    p = NEW(CrxTBondPriceOption);
    if (p==NULL) goto done;

    p->bond            = GtoCopyTBond(bond);
    if (p->bond == NULL) goto done;

    p->optionType      = optionType;
    p->exerciseDate    = exerciseDate;
    p->paymentDate     = paymentDate;
    p->strikePrice     = strikePrice;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceOptionFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondPriceOption
***************************************************************************
*/
CrxTBondPriceOption* CrxBondPriceOptionMakeEmpty(void)
{
    static char routine[] = "CrxBondPriceOptionMakeEmpty";
    int status = FAILURE;

    CrxTBondPriceOption* p = NULL;


    p = NEW(CrxTBondPriceOption);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceOptionFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondPriceOption
***************************************************************************
*/
CrxTBondPriceOption* CrxBondPriceOptionCopy(CrxTBondPriceOption* src)
{
    CrxTBondPriceOption* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondPriceOptionMake(src->bond,
                                 src->optionType,
                                 src->exerciseDate,
                                 src->paymentDate,
                                 src->strikePrice);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondPriceOption
***************************************************************************
*/
void CrxBondPriceOptionFree(CrxTBondPriceOption *p)
{
    if (p != NULL)
    {
        GtoFreeTBond(p->bond);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondPriceOptionCalc
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalcMake(
double          optionPrice,         /* (I) */
double          repoRate,            /* (I) */
double          fwdPrice,            /* (I) */
double          fwdPremium           /* (I) */
)
{
    static char routine[] = "CrxBondPriceOptionCalcMake";
    int status = FAILURE;

    CrxTBondPriceOptionCalc* p = NULL;


    p = NEW(CrxTBondPriceOptionCalc);
    if (p==NULL) goto done;

    p->optionPrice     = optionPrice;
    p->repoRate        = repoRate;
    p->fwdPrice        = fwdPrice;
    p->fwdPremium      = fwdPremium;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceOptionCalcFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondPriceOptionCalc
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalcMakeEmpty(void)
{
    static char routine[] = "CrxBondPriceOptionCalcMakeEmpty";
    int status = FAILURE;

    CrxTBondPriceOptionCalc* p = NULL;


    p = NEW(CrxTBondPriceOptionCalc);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondPriceOptionCalcFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondPriceOptionCalc
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalcCopy(CrxTBondPriceOptionCalc* src)
{
    CrxTBondPriceOptionCalc* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondPriceOptionCalcMake(src->optionPrice,
                                     src->repoRate,
                                     src->fwdPrice,
                                     src->fwdPremium);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondPriceOptionCalc
***************************************************************************
*/
void CrxBondPriceOptionCalcFree(CrxTBondPriceOptionCalc *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondSpreadOption
***************************************************************************
*/
CrxTBondSpreadOption* CrxBondSpreadOptionMake(
TBond*          bond,                /* (I) */
TBond*          refBond,             /* (I) */
long            optionType,          /* (I) */
TDate           exerciseDate,        /* (I) */
TDate           paymentDate,         /* (I) */
double          strikeSpread         /* (I) */
)
{
    static char routine[] = "CrxBondSpreadOptionMake";
    int status = FAILURE;

    CrxTBondSpreadOption* p = NULL;

    REQUIRE(bond != NULL);
    REQUIRE(refBond != NULL);

    p = NEW(CrxTBondSpreadOption);
    if (p==NULL) goto done;

    p->bond            = GtoCopyTBond(bond);
    if (p->bond == NULL) goto done;

    p->refBond         = GtoCopyTBond(refBond);
    if (p->refBond == NULL) goto done;

    p->optionType      = optionType;
    p->exerciseDate    = exerciseDate;
    p->paymentDate     = paymentDate;
    p->strikeSpread    = strikeSpread;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondSpreadOptionFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondSpreadOption
***************************************************************************
*/
CrxTBondSpreadOption* CrxBondSpreadOptionMakeEmpty(void)
{
    static char routine[] = "CrxBondSpreadOptionMakeEmpty";
    int status = FAILURE;

    CrxTBondSpreadOption* p = NULL;


    p = NEW(CrxTBondSpreadOption);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondSpreadOptionFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondSpreadOption
***************************************************************************
*/
CrxTBondSpreadOption* CrxBondSpreadOptionCopy(CrxTBondSpreadOption* src)
{
    CrxTBondSpreadOption* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondSpreadOptionMake(src->bond,
                                  src->refBond,
                                  src->optionType,
                                  src->exerciseDate,
                                  src->paymentDate,
                                  src->strikeSpread);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondSpreadOption
***************************************************************************
*/
void CrxBondSpreadOptionFree(CrxTBondSpreadOption *p)
{
    if (p != NULL)
    {
        GtoFreeTBond(p->bond);
        GtoFreeTBond(p->refBond);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for CrxTBondSpreadOptionCalc
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalcMake(
double          optionPrice,         /* (I) */
double          repoRate,            /* (I) */
double          refRepoRate,         /* (I) */
double          fwdPrice,            /* (I) */
double          refFwdPrice,         /* (I) */
double          fwdSpread,           /* (I) */
double          dpdy,                /* (I) */
double          spreadPremium        /* (I) */
)
{
    static char routine[] = "CrxBondSpreadOptionCalcMake";
    int status = FAILURE;

    CrxTBondSpreadOptionCalc* p = NULL;


    p = NEW(CrxTBondSpreadOptionCalc);
    if (p==NULL) goto done;

    p->optionPrice     = optionPrice;
    p->repoRate        = repoRate;
    p->refRepoRate     = refRepoRate;
    p->fwdPrice        = fwdPrice;
    p->refFwdPrice     = refFwdPrice;
    p->fwdSpread       = fwdSpread;
    p->dpdy            = dpdy;
    p->spreadPremium   = spreadPremium;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondSpreadOptionCalcFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for CrxTBondSpreadOptionCalc
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalcMakeEmpty(void)
{
    static char routine[] = "CrxBondSpreadOptionCalcMakeEmpty";
    int status = FAILURE;

    CrxTBondSpreadOptionCalc* p = NULL;


    p = NEW(CrxTBondSpreadOptionCalc);
    if (p==NULL) goto done;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxBondSpreadOptionCalcFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for CrxTBondSpreadOptionCalc
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalcCopy(CrxTBondSpreadOptionCalc* src)
{
    CrxTBondSpreadOptionCalc* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxBondSpreadOptionCalcMake(src->optionPrice,
                                      src->repoRate,
                                      src->refRepoRate,
                                      src->fwdPrice,
                                      src->refFwdPrice,
                                      src->fwdSpread,
                                      src->dpdy,
                                      src->spreadPremium);

    return dst;
}

/**
***************************************************************************
** Destructor for CrxTBondSpreadOptionCalc
***************************************************************************
*/
void CrxBondSpreadOptionCalcFree(CrxTBondSpreadOptionCalc *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

