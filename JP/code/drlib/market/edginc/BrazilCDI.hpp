//----------------------------------------------------------------------------
//
//   Group       : QR Equities NY
//
//   Filename    : BrazilCDI.hpp
//
//   Description : Barilian CDI Curve
//
//   Author      : Dmytro Zhuravytsky
//
//   Date        : 15 March 2006
//
//----------------------------------------------------------------------------

#ifndef BRAZIL_CDI_HPP
#define BRAZIL_CDI_HPP

#include "edginc/ZeroCurve.hpp"
#include "edginc/ZeroCurveFactory.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL BrazilCDICurve : public ZeroCurve
{
public:
    static CClassConstSP const TYPE;

    BrazilCDICurve();
    BrazilCDICurve(
        const ZeroCurveBenchmarkArray & benchmarks,
        const DateTime & refDate,
        const HolidayConstSP & holidays,
        int spotOffset );

    /** Calculates discount factor for a date */
    virtual double discountFactor(const DateTime& date) const;

    /** Calculates an zero coupon rate from a ZCurve at some date */
    virtual double zeroCouponRate(const DateTime& date) const;

    /** how long is the curve ? */
    virtual int length() const { return dates.size(); }
    
    /**
     * The date returned is the first date for which we have genuine
     * information regarding rates and discount factors. It is possible that
     * using the curve for dates before this date will give answers, but they
     * might not be based on real information.
     */
    virtual const DateTime& firstDate() const { return dates.front(); }

    /** when does it end ? */
    virtual const DateTime & endDate() const { return dates.back(); }
    
    /** strip out the dates */
    virtual DateTimeArray getDates() const { return DateTimeArray( dates ); }

    /** strip out the rates and dates */
    virtual CashFlowArraySP getRatesAndDates() const { return CashFlow::createCashFlowArray( dates, rates ); }

    virtual const DateTime& getBaseDate() const { return baseDate; }

    /**
     * Returns a key used to optimize repeated calculations of discount
     * factors/forward rate. The calc method for this key returns the 
     * natural logarithm of the discount factor.
     */
    virtual YieldCurve::IKey* logOfDiscFactorKey() const;

    // Interpolates log of discount for a date
    double logDiscFactor(const DateTime& date) const;

private:
    friend class BrazilCDICurveHelper;

    BrazilCDICurve(const BrazilCDICurve & BrazilCDICurve);
    BrazilCDICurve & operator =(const BrazilCDICurve & BrazilCDICurve);

    DateTime baseDate;
    DateTimeArray dates;
    DoubleArray rates;
    DayCountConventionConstSP dcc;

    DoubleArray yearFrac;
    DoubleArray logDisc;
};

class MARKET_DLL BrazilCDIFactory : public IZeroCurveFactory
{
public:
    static CClassConstSP const TYPE;

    BrazilCDIFactory();

    virtual bool equals(const IZeroCurveFactory* zcm) const;

    virtual bool isFloatLegValued() const;

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
        const Expiry*                  futuresMaturity ) const;

    /**
     * Create zero curves.
     *
     * It is part of the contract for this factory that this method creates
     * both the discounting and estimating curves in a single call.  These
     * could be the same.
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
        const Expiry*                  futuresMaturity ) const;

    friend class BrazilCDIFactoryHelper;

protected:
    ZeroCurveBenchmarkArraySP adjustBenchmarks(
        const ZeroCurveBenchmarkArray & benchmarks,
        const DateTime                & refDate,
        const Holiday                 & holidays ) const;
};

DRLIB_END_NAMESPACE

#endif
