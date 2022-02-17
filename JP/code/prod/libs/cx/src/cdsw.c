/*
***************************************************************************
** FILE NAME: cdsw.c
**
** CDSW functions. These are functions which emulate the CDSW page on
** Bloomberg.
**
** $Header$
***************************************************************************
*/

#include "cdsw.h"

#include "cds.h"
#include "cdsbootstrap.h"
#include "recovery.h"

#include <alib/cerror.h>
#include <alib/rtbrent.h>

typedef struct
{
    TDate           today;
    TDate           valueDate;
    TDate           startDate;
    TDate           endDate;
    double          couponRate;
    TBoolean        payAccruedOnDefault;
    TDateInterval  *dateInterval;
    CxTStubType     stubType;
    CxTDayCountConv accrueDCC;
    CxTBadDayConv   badDayConv;
    CxTCalendar    *calendar;
    TCurve         *discCurve;
    double          upfrontCharge;
    double          recoveryRate;
    TBoolean        inclusiveProtection;
    TBoolean        payAccruedAtStart;
} CDSW_SPREAD_CONTEXT;


/* static function declarations */
static int cdswSpreadSolverFunction
(double               bbgSpread,
 CDSW_SPREAD_CONTEXT *context,
 double              *diff);

/*
 * Computes the upfront charge for a flat spread par curve.
 */
int CxCdswUpfrontCharge
(TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 double          couponRate,
 TBoolean        payAccruedOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 CxTDayCountConv accrueDCC,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TCurve         *discCurve,
 double          bbgSpread,
 double          recoveryRate,
 TBoolean        inclusiveProtection,
 TBoolean        payAccruedAtStart,
 double         *upfrontCharge)
{
    static char routine[] = "CxCdswUpfrontCharge";
    int         status    = FAILURE;

    CxTCreditCurve *flatSpreadCurve = NULL;
    CxTRecoveryCurve *recoveryCurve = NULL;

    recoveryCurve = CxRecoveryCurveMakeFromRecoveryRate(recoveryRate);
    if (recoveryCurve == NULL)
        goto done; /* failure */

    flatSpreadCurve = CxCdsBootstrap (
        today,
        discCurve,
        startDate,
        valueDate,
        1,
        &endDate,
        &bbgSpread,
        NULL,
        NULL,
        recoveryCurve,
        payAccruedOnDefault,
        dateInterval,
        accrueDCC,
        stubType,
        CX_CURVE_TYPE_FLOW,
        NULL,
        NULL,
        inclusiveProtection,
        FALSE, /* isPriceClean - irrelevant in this context */
        0,     /* delay */
        badDayConv,
        calendar);

    if (flatSpreadCurve == NULL)
        goto done; /* failure */

    if (CxCdsPrice (today,
                    valueDate,
                    startDate,
                    endDate,
                    0, /* delay */
                    couponRate,
                    payAccruedOnDefault,
                    dateInterval,
                    stubType,
                    accrueDCC,
                    badDayConv,
                    calendar,
                    discCurve,
                    flatSpreadCurve,
                    recoveryCurve,
                    inclusiveProtection,
                    payAccruedAtStart,
                    upfrontCharge) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    CxCreditCurveFree (flatSpreadCurve);
    CxRecoveryCurveFree (recoveryCurve);

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*
 * Computes the flat spread required to match the upfront charge.
 */
int CxCdswSpread
(TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 double          couponRate,
 TBoolean        payAccruedOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 CxTDayCountConv accrueDCC,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TCurve         *discCurve,
 double          upfrontCharge,
 double          recoveryRate,
 TBoolean        inclusiveProtection,
 TBoolean        payAccruedAtStart,
 double         *bbgSpread)
{
    static char routine[] = "CxCdswSpread";
    int         status    = FAILURE;

    CDSW_SPREAD_CONTEXT context;

    context.today               = today;
    context.valueDate           = valueDate;
    context.startDate           = startDate;
    context.endDate             = endDate;
    context.couponRate          = couponRate;
    context.payAccruedOnDefault = payAccruedOnDefault;
    context.dateInterval        = dateInterval;
    context.stubType            = stubType;
    context.accrueDCC           = accrueDCC;
    context.badDayConv          = badDayConv;
    context.calendar            = calendar;
    context.discCurve           = discCurve;
    context.upfrontCharge       = upfrontCharge;
    context.recoveryRate        = recoveryRate;
    context.inclusiveProtection = inclusiveProtection;
    context.payAccruedAtStart   = payAccruedAtStart;

    if (GtoRootFindBrent ((TObjectFunc)cdswSpreadSolverFunction,
                          &context,
                          0.0,    /* boundLo */
                          100.0,  /* boundHi */
                          100,    /* numIterations */
                          0.01,   /* guess */
                          0.0001, /* initialXStep */
                          0.0,    /* initialFDeriv */
                          1e-8,   /* xacc */
                          1e-8,   /* facc */
                          bbgSpread) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);
    return status;
}


/* static function declarations */
static int cdswSpreadSolverFunction
(double               bbgSpread,
 CDSW_SPREAD_CONTEXT *context,
 double              *diff)
{
    double upfrontCharge;

    if (CxCdswUpfrontCharge (context->today,
                             context->valueDate,
                             context->startDate,
                             context->endDate,
                             context->couponRate,
                             context->payAccruedOnDefault,
                             context->dateInterval,
                             context->stubType,
                             context->accrueDCC,
                             context->badDayConv,
                             context->calendar,
                             context->discCurve,
                             bbgSpread,
                             context->recoveryRate,
                             context->inclusiveProtection,
                             context->payAccruedAtStart,
                             &upfrontCharge) != SUCCESS)
        return FAILURE;

    *diff = upfrontCharge - context->upfrontCharge;
    return SUCCESS;
}
