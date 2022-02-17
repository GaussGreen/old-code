//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids Derivatives Research
//
//   Filename    : IZeroCurveFactory.hpp
//
//   Description : Defines interface for zero curve factories.
//
//   Author      : Gordon C Stephens
//
//   Date        : 18 January 2005
//
//----------------------------------------------------------------------------

#ifndef ICURVEFACTORY_HPP
#define ICURVEFACTORY_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/ZeroPair.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/ZeroCurveBenchmark.hpp"


DRLIB_BEGIN_NAMESPACE

class BootstrappedYieldCurve;



/**
 * Interface for zero curve factory classes.
 */
class MARKET_DLL IZeroCurveFactory : public CObject
{
public:
    static CClassConstSP const TYPE;

    virtual bool equals(const IZeroCurveFactory* zcm) const = 0;

    /**
     * Test if floating leg is to be valued.
     */
    virtual bool isFloatLegValued() const = 0;

    /**
     * Validate input data, and throw exception if not valid.
     */
	virtual void validate(
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
        ) const = 0;

    /**
     * Create zero curves.
     * 
     * It is part of the contract for this factory that this method creates
     * both the discounting and estimating curves in a single call.  These
     * could be the same.
     *
     * validate() must be called sometime before calling this method : for 
     * performance reasons it is not called from within this method.
     */
	virtual ZeroPairSP bootstrap(
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
        ) const = 0;
    
    virtual void setAdjustDates(bool adjust);

    /**
     * Use default rates to convert curve to a risky curve.
     */
    virtual IYieldCurve* makeRiskyCurve(
        const BootstrappedYieldCurve& riskless, 
        const CashFlowArray&          defaultRates, 
        double                        recovery) const;

protected:
    IZeroCurveFactory(CClassConstSP clazz);

    bool adjustDates; // bad day adjust cash end / futures start dates?

private:
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<IZeroCurveFactory> IZeroCurveFactoryConstSP;
typedef smartPtr<IZeroCurveFactory>      IZeroCurveFactorySP;

DRLIB_END_NAMESPACE


#endif  // ICURVEFACTORY_HPP
