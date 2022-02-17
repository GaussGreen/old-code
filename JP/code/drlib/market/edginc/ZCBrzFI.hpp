//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZCBrzFI.hpp
//
//   Description : Curve factory implementation ported from ALIB BRZ_ZC_FI.
//                 Creates a zero curve from Brazilian benchmarks and COPOM
//                 date information.  The input data consists of:
//                 1) CDI overnight rate (quoted as 30D overnight rate)
//                 2) next CDI overnoght rate (quoted as 30D overnight rate)
//                 3) 1 day CD futures prices for a sequence of benchmark
//                    instruments in increasing order of maturities (quoted as 
//                    the present value of on PU [price unit] at maturity)
//                 4) CD versus prefixed swap rates for a sequence of benchmark 
//                    swaps in increasing order of maturities where the first 
//                    swap maturity comes after the last futures maturity
//                 5) COPOM meeting dates and possibly an initial forward rate
//                    forecast.
//
//   Author      : Richard Appleton
//
//   Date        : 24th March 2006
//
//----------------------------------------------------------------------------

#ifndef ZC_BRZ_FI_HPP
#define ZC_BRZ_FI_HPP

#include "edginc/DateTime.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include <string>



DRLIB_BEGIN_NAMESPACE
using std::string;


class MARKET_DLL ZCBrzFI : public IZeroCurveFactory
{
public:
    static CClassConstSP const TYPE;

    ZCBrzFI();

    virtual ~ZCBrzFI();

    int          hashCode() const;

    virtual bool equals(const IZeroCurveFactory* zcm) const;

    virtual bool isFloatLegValued() const;

    virtual void validatePop2Object();

    void validate
    (
        const DateTime&                today,
        const DateTime&                spotDate,
        const ZeroCurveBenchmarkArray& benchmarks,
        const DayCountConvention*      moneyMarketDcc,
        const DayCountConvention*      swapFixedDcc,
        const DayCountConvention*      swapFloatDcc,
        const DayCountConvention*      swapBasisDcc,
        const MaturityPeriod*          swapFixedIvl,
        const MaturityPeriod*          swapFloatIvl,
        const MaturityPeriod*          swapBasisIvl,
        const BadDayConvention*        swapBadDayConvention,
        const bool                     floatRateFixed,
        const ExpiryArray*             fixDates,
        const DoubleArray*             fixRates,
        IRVolBaseConstSP               volModelIR,
        HolidayConstSP                 holidays,
        const CurrencyBasisWrapper&    ccyBasis,
        const Expiry*                  futuresMaturity
    ) const;

    ZeroPairSP bootstrap
    (
        const DateTime&                today,
        const DateTime&                spotDate,
        const ZeroCurveBenchmarkArray& benchmarks,
        const DayCountConvention*      moneyMarketDcc,
        const DayCountConvention*      swapFixedDcc,
        const DayCountConvention*      swapFloatDcc,
        const DayCountConvention*      swapBasisDcc,
        const MaturityPeriod*          swapFixedIvl,
        const MaturityPeriod*          swapFloatIvl,
        const MaturityPeriod*          swapBasisIvl,
        const BadDayConvention*        swapBadDayConvention,
        const bool                     floatRateFixed,
        const ExpiryArray*             fixDates,
        const DoubleArray*             fixRates,
        IRVolBaseConstSP               volModelIR,
        HolidayConstSP                 holidays,
        const CurrencyBasisWrapper&    ccyBasis,
        const Expiry*                  futuresMaturity
    ) const;

private:
    static void load(CClassSP& clazz);

    static double BRZ_FUTURES_PRICE_UNIT;

    //fields
    DateTimeArraySP copomDates;
};


DRLIB_END_NAMESPACE

#endif // ZC_BRZ_FI_HPP
