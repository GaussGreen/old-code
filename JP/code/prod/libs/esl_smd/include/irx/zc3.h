#ifndef IRX_ZC3_H
#define IRX_ZC3_H

#include "irxflow.h"

/**
 * Zero curve bootstrapping without currency basis.
 */
IrxTSwapZeroCurve* irxZeroCurve3
(IrxTBootstrapMethod method,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             prices,
 long*               benchmarkFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTBool            valueFloating,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate);
 

/**
 * Zero curve bootstrapping with currency basis.
 *
 * Same as ZeroCurve3 but with three extra parameters (cbsRates, cbsCalendar,
 * shortCbsSpreads) and one less parameter (valueFloating).
 *
 * Also uses more information from the marketConv.
 */
IrxTSwapZeroCurve* irxZeroCurve3CB
(IrxTBootstrapMethod method,
 IrxTBool            reuseDiscountCurve,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             cbsRates,
 double*             prices,
 long*               benchmarkFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTCalendar*       cbsCalendar,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate);

#endif
