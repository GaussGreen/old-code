//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : YieldCurveAdapter.hpp
//
//   Description : Adapter to convert yield curve to zero curve.
//                 Enables YieldCurves to be used in SwapTool methods.
//
//   Author      : Richard Appleton
//
//   Date        : 27th February 2006
//
//----------------------------------------------------------------------------

#ifndef YIELD_CURVE_ADAPTER_HPP
#define YIELD_CURVE_ADAPTER_HPP

#include "edginc/YieldCurve.hpp"
#include "edginc/ZeroCurve.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL YieldCurveAdapter : public ZeroCurve
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    YieldCurveAdapter(const IYieldCurve& yc, bool projection = false);
    virtual ~YieldCurveAdapter();
    
    /** Calculates discount factor for a date */
    virtual double discountFactor(const DateTime& date) const;

    /** Calculates an zero coupon rate from a ZCurve at some date */
    virtual double zeroCouponRate(const DateTime& date) const;

    /** how long is the curve ? */
    virtual int length() const;
    
    /**
     * The date returned is the first date for which we have genuine
     * information regarding rates and discount factors. It is possible that
     * using the curve for dates before this date will give answers, but they
     * might not be based on real information.
     */
    virtual const DateTime& firstDate() const;

    /** when does it end ? */
    virtual const DateTime& endDate() const;
    
    /** strip out the dates */
    virtual DateTimeArray getDates() const;

    /** strip out the rates and dates */
    virtual CashFlowArraySP getRatesAndDates() const;

    virtual const DateTime& getBaseDate() const;

    /**
     * Returns a key used to optimize repeated calculations of discount
     * factors/forward rate. The calc method for this key returns the 
     * natural logarithm of the discount factor.  The default implementation
     * on this class just takes the log of the result from calling pv().
     */
    virtual YieldCurve::IKey* logOfDiscFactorKey() const;

private:
    class LogOfDiscFactorKey;

    YieldCurveAdapter();

    const IYieldCurve& get() const;

    IYieldCurveSP      yc;  // not const beacuse of reflection macros - but treated as const
    bool               projection;
    mutable DateTime   baseDate; // $unregistered
    mutable DateTime   firstDt; // $unregistered
};



DRLIB_END_NAMESPACE
#endif // YIELD_CURVE_ADAPTER_HPP
