//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZeroCurve3.hpp
//
//   Description : Curve factory implementation for ALIB zero curve 3 method
//                 Implementation of the ALIB zero curve "best of breed" 
//                 methodology. The code has been extracted from ALIB and 
//                 converted to make use of native EDR classes and data types
//
//   Author      : Richard Appleton
//
//   Date        : 25th April 2005
//
//----------------------------------------------------------------------------

#ifndef ZERO_CURVE3_HPP
#define ZERO_CURVE3_HPP

#include "edginc/DateTime.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include <string>



DRLIB_BEGIN_NAMESPACE
using std::string;

class RateData;
class SwapData;
class Turn;
class IRVolBase;



class MARKET_DLL ZeroCurve3 : public IZeroCurveFactory
{
public:
    static CClassConstSP const TYPE;

    ZeroCurve3();

    virtual ~ZeroCurve3();

    int          hashCode() const;

    virtual bool equals(const IZeroCurveFactory* zcm) const;

    virtual bool isFloatLegValued() const;

    virtual void validatePop2Object();

	void validate(
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

	ZeroPairSP bootstrap(
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

    ZeroPairSP bootstrap(
        const DateTime&                    valueDate,
        const ZC3ZeroCurve&                stub,
        const ZeroCurveBenchmarkArray&     benchmarks,
        const DayCountConvention&          moneyMarketDcc,
        const DateTime*                    futuresMMDt,
        const Expiry*                      futuresMaturity,
        const DayCountConvention&          swapFixedDcc,
        const DayCountConvention*          swapFloatDcc,
        const DayCountConvention*          swapBasisDcc,
        const MaturityPeriod&              swapFixedIvl,
        const MaturityPeriod*              swapFloatIvl,
        const MaturityPeriod*              swapBasisIvl,
        const BadDayConvention&            swapBadDayConvention,
        const bool                         floatRateFixed,
        const ZeroCurve*                   fixingCurve,
        const IRVolBase*                   volModelIR,
        const Holiday&                     holidays,
        const CurrencyBasisWrapper&        ccyBasis,
        const ZC3ZeroInterpolationArraySP& interpolators,
        const DateTime*                    extrapDt
        ) const;

    ZC3ZeroCurveSP getStubCurve(
        const DateTime& baseDate, 
        const string& stubInterpType) const;

    /**
     * Parse input instrument data.
     */
    void splitData(
        ZC3RateDataArray&              ratesData,
        ZC3SwapDataArray&              swapsData,
        ZC3TurnArray&                  turnsData,
        const ZC3ZeroCurve&            szcStub,
        const ZeroCurveBenchmarkArray& benchmarks,
        const DayCountConvention&      mmDcc,
        const Expiry*                  futuresMaturity,
        const bool                     valueFloating,
        const DayCountConvention&      fixedDcc,
        const DayCountConvention*      floatDcc,       // only used if value floating
        const DayCountConvention*      basisDcc,       // only used if basis rates not NULL
        const MaturityPeriod&          fixedIvl,
        const MaturityPeriod*          floatIvl,       // only used if value floating
        const MaturityPeriod*          basisIvl,       // only used if basis rates not NULL
        const bool                     floatRateFixed, // only used if value floating
        const ZeroCurve*               fixCurve,       // only used if float rate fixed
        const BadDayConvention&        badDayConv,
        const Holiday&                 holidays,
        const CurrencyBasisWrapper&    ccyBasis
        ) const;

    /*
     * Select parameters for segment bootstrap where curve has multiple 
     * interpolation types.
     */
    void benchmarksExtractInstruments(
        const ZC3RateDataArray& ratesData,
        const ZC3SwapDataArray& swapsData,
        const ZC3TurnArray&     turnsData,
        const DateTime&         from,
        const DateTime&         to,
        bool                    recursive,
        ZC3RateDataArray&       ratesDataOut,
        ZC3SwapDataArray&       swapsDataOut,
        ZC3TurnArray&           turnsDataOut,
        DateTime&               lastDateUsed,
        bool&                   doSomething,
        bool&                   lastSection,
        bool&                   futuresUsed
        ) const;

    //fields
    bool     valueFloating;       // is the floating leg of the swap considered when building
    bool     mtmAdjustmentFlag;   // futures flag #1
    string   futuresStubType;     // futures flag #2
    string   futuresGapMethod;    // futures flag #3
    int      futuresGapDays;      // futures flag #4
    ExpirySP futuresMMDate;
    string   futuresMMDateBadDayConvention;
    ExpirySP extrapDate;
    string   stubType;
    string   couponInterpolationMethod; // controls how the zc is built
    string   zeroInterpolationMethod;   // controls how the zc is interpolated
};


DRLIB_END_NAMESPACE

#endif // ZERO_CURVE3_HPP
