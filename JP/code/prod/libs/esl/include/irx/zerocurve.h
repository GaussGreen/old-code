#ifndef IRX_ZERO_CURVE_H
#define IRX_ZERO_CURVE_H

#include "irxflow.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Constructs a zero curve from a set of input dates and rates. The
 * resulting curve will use flat forward interpolation between dates.
 *
 * Note that we do not actually store the zero curve as dates and rates.
 */
IrxTZeroCurve* irxZeroCurveMakeFromRates(
    IrxTDate          baseDate,
    int               numDates,
    const IrxTDate   *zeroDates,
    const double     *zeroRates,
    IrxTRateType      rateType,
    IrxTDayCountConv  dcc);

/**
 * Like one above but requires empty curve object
 */
int irxZeroCurveConstructFromRates(
    IrxTZeroCurve* crv,
    IrxTDate          baseDate,
    int               numDates,
    const IrxTDate   *zeroDates,
    const double     *zeroRates,
    IrxTRateType      rateType,
    IrxTDayCountConv  dcc);

/**
 * Like the ones above, except set the interpMethod explicitly
 */
IrxTZeroCurve* irxZeroCurveMakeFromRates2(
    IrxTDate             baseDate,
    int                  numDates,
    const IrxTDate*      zeroDates,
    const double*        zeroRates,
    IrxTRateType         rateType,
    IrxTBootstrapMethod  interpMethod,
    IrxTDayCountConv     dcc);

int irxZeroCurveConstructFromRates2(
    IrxTZeroCurve*       crv,
    IrxTDate             baseDate,
    int                  numDates,
    const IrxTDate*      zeroDates,
    const double*        zeroRates,
    IrxTRateType         rateType,
    IrxTBootstrapMethod  interpMethod,
    IrxTDayCountConv     dcc);

/**
 * Constructs a zero curve from a set of input dates and prices. The
 * resulting curve will use flat forward interpolation between dates.
 *
 * You can provide the baseDate as one of the zeroDates, but in this
 * case the price for that date must be 1.
 *
 * Note that we do not actually store the zero curve solely as dates
 * and prices.
 */
IrxTZeroCurve* irxZeroCurveMakeFromPrices(
    IrxTDate          baseDate,
    int               numDates,
    const IrxTDate   *zeroDates,
    const double     *zeroPrices);


/**---------------------------------------------------------
 * I/O: read IrxTZeroCurve from wrapper file (skip comments).
     * (1) London format:
     * # Start date
     * 19970411
     * # Money Market basis (360 or 365)
     * 360
     * # Annual or semi-annual curve ("A" or "S")
     * S
     * # Year basis for benchmark swaps ("ACT", "365" or "360")
     * ACT
     * # No of entries
     * 84
     * #zero maturity yyyymmdd rates (ACT/365F annual)
     * 19970412 5.734349908886
     * 19970511 5.913551526763
     * etc
 */
int  irxZeroCurveConstructFromDRWFile(IrxTZeroCurve *zc, FILE *fp);

/**
 * Accessor to get base date from the supplied zero curve.
 */
IrxTDate irxZeroCurveBaseDate(const IrxTZeroCurve *zeroCurve);

/**
 * Accessor to get last date from the supplied zero curve.
 */
IrxTDate irxZeroCurveLastDate(const IrxTZeroCurve *zeroCurve);

/**
 * Calculates the zero price for a given date. This is the price as seen
 * from the base date of the zero curve.
 */
int irxZeroPrice(
    const IrxTZeroCurve *zc,
    IrxTDate             date,
    double              *price);

/**
 * Calculates the zero price for a given start date and end date.
 */
int irxFwdZeroPrice(
    const IrxTZeroCurve *zc,
    IrxTDate             startDate,
    IrxTDate             endDate,
    double              *fwdPrice);

/**
 * Calculates the spot zero rate from the base date to the given date
 */
int irxZeroRate(
    const IrxTZeroCurve *zc,
    IrxTDate             endDate,
    IrxTDayCountConv     dcc,
    IrxTRateType         rateType,
    double              *rate);

/**
 * Calculates a forward zero rate for a given start date and end date.
 */
int irxFwdZeroRate(
    const IrxTZeroCurve *zc,
    IrxTDate             startDate,
    IrxTDate             endDate,
    IrxTDayCountConv     dcc,
    IrxTRateType         rateType,
    double              *fwdRate);

/**
 * Extend the zero curve to Today (from ValueDate)
 * look at  QLib/rates/libs/srm3/src/paryield.c
 * Returns the discount factor between today and baseDate
 */
double irxZeroCurveExtendToToday(IrxTZeroCurve* zc);

/**
 * Set the zero curve's notion of Today and extend the zero curve to Today
 * (from ValueDate); look at  QLib/rates/libs/srm3/src/paryield.c
 * Returns the discount factor between today and baseDate
 */
double irxZeroCurveSetAndExtendToToday(IrxTZeroCurve* zc, IrxTDate newToday);

#ifdef __cplusplus
}
#endif

#endif /* IRX_ZERO_CURVE_H */

