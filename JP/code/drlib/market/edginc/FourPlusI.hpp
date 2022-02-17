//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourPlusI.hpp
//
//   Description : Defines how swap rates are added to a zero curve
//
//   Author      : Andrew J Swain
//
//   Date        : 13 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef FOURPLUSI_HPP
#define FOURPLUSI_HPP

#include <string>
#include "edginc/DateTime.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/FourPlusIZeroCurve.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/ZeroCurve.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE


/** 
 * This curve factory builds a 4+i style zero curve.
 *
 */
class MARKET_DLL FourPlusI : public IZeroCurveFactory {
public:
    static CClassConstSP const TYPE;

	FourPlusI();

    virtual ~FourPlusI();

    int          hashCode() const;

    virtual bool equals(const IZeroCurveFactory* zcm) const;

    virtual bool isFloatLegValued() const;

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

    /**
     * Use default rates to convert curve to a risky curve.
     */
    virtual IYieldCurve* makeRiskyCurve(
        const BootstrappedYieldCurve& riskless, 
        const CashFlowArray&          defaultRates, 
        double                        recovery) const;

private:

		//build standard curve, using pre-built discount curve
		void buildZeroCurve(
                            FourPlusIZeroCurve*       ref,
                            const ZeroCurve*          discount,
                            DateTime                  today,
                            DateTime                  spotDate,
                            HolidayConstSP            rateHolidays,
                            const DayCountConvention& moneyMarketDayCount,
                            const ExpiryArray&        moneyMarketExpiries,
                            const DoubleArray&        moneyMarketRates,
                            const BadDayConvention*   swapBadDayConvention,
                            const DayCountConvention* swapFixedDayCount,
                            const DayCountConvention* swapFloatDayCount,
                            int                       swapFixedLegFrequency,
                            int                       swapFloatLegFrequency,
                            ExpiryArray&              swapExpiries,
                            DoubleArray&              swapRates,
                            const ExpiryArray*        futureExpiries,    //expect NULL if not being used
                            const DoubleArray*        futureRates,       //expect NULL if not being used
                            const Expiry*             futureMaturity,
                            HolidayConstSP            basisHolidays,
                            const DoubleArray*        basisRates,        //expect NULL if not being used
                            bool                      buildProjectionCurve
                            ) const;

    // parse input data
    void parseData(
        const ZeroCurveBenchmarkArray& benchmarks,
        StringArray&                   instruments,
        ExpiryArray&                   expiries,
        DoubleArray&                   rates,
        ExpiryArray&                   moneyMarketExpiries,
        ExpiryArraySP&                 futExpiries,
        ExpiryArray&                   swapExpiries,
        DoubleArray&                   moneyMarketRates,
        DoubleArraySP&                 futRates,
        DoubleArray&                   swapRates,
        bool&                          useFutures) const;

    // check inputs
    void validate(const DateTime&      today,
                  const DateTime&      spotDate,
                  const int            swapFixedLegFrequency,
                  const DoubleArray&   moneyMarketRates,
                  const DoubleArray&   swapRates,
                  const DateTimeArray& bMarks,
                  const bool           badDayAdjust,
                  const ExpiryArray*   futureExpiries,
                  const DoubleArray*   futureRates) const;

    double getBasisRate(
        const ZeroCurveBenchmark&   benchmark,
        const DayCountConvention&   moneyMarketDcc,
        const DayCountConvention&   swapFixedDcc,
        const MaturityPeriod&       swapFixedIvl,
        const MaturityPeriod&       swapFloatIvl,
        const CurrencyBasisWrapper& ccyBasis) const;

    /** add swaps to a zero curve */
    void addSwaps(
        const DateTime&           valueDate,
        const DateTime&           minDate,   // none before here
        const DateTimeArray&      swapDates,
        const DoubleArray&        swapRates,
        int                       swapFrequency,
        const DayCountConvention* dcc,
        FourPlusIZeroCurve*       zc) const;

    // advanced version that values both fixed and floating leg of swap
    void addSwaps(
        const ZeroCurve*          discount,         // (I) 
        const DateTime&           valueDate,
        const DateTime&           minDate,          // none before here
        const DateTimeArray&      swapDates,
        const DoubleArray&        swapRates,
        int                       fixedSwapFreq,    // (I) 
        int                       floatSwapFreq,    // (I) 
        const DayCountConvention* fixedDCC,         // (I) 
        const DayCountConvention* floatDCC,         // (I) 
        FourPlusIZeroCurve*       zc) const;

    friend class FourPlusIHelper;

    // fields
    string zeroInterpolationMethod;   // how to interp on zero curve
    bool   badDayAdjust; // $unregistered
};

DRLIB_END_NAMESPACE

#endif
