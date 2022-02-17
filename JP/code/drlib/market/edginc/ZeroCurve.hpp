//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZeroCurve.hpp
//
//   Description : Abstract base class for zero curves
//
//   Author      : Richard Appleton
//
//   Date        : 4th July 2005
//
//----------------------------------------------------------------------------

#ifndef ZERO_CURVE_HPP
#define ZERO_CURVE_HPP

#include <string>
#include "edginc/DayCountConvention.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

class Stub;

/**
 * Abstract base class for all zero curve implementations.
 *
 * The methods implemented within this class are independent of how the zero
 * curves is constructed or its internal data structures.
 */
class MARKET_DLL ZeroCurve : public CObject
{
public:
    static CClassConstSP const TYPE;

    virtual ~ZeroCurve();
    
    /** Calculates discount factor for a date */
    virtual double discountFactor(const DateTime& date) const = 0;

    /** Calculates an zero coupon rate from a ZCurve at some date */
    virtual double zeroCouponRate(const DateTime& date) const = 0;

    /** how long is the curve ? */
    virtual int length() const = 0;
    
    /**
     * The date returned is the first date for which we have genuine
     * information regarding rates and discount factors. It is possible that
     * using the curve for dates before this date will give answers, but they
     * might not be based on real information.
     */
    virtual const DateTime& firstDate() const = 0;

    /** when does it end ? */
    virtual const DateTime& endDate() const = 0;
    
    /** strip out the dates */
    virtual DateTimeArray getDates() const = 0;

    /** strip out the rates and dates */
    virtual CashFlowArraySP getRatesAndDates() const = 0;

    virtual const DateTime& getBaseDate() const = 0;

    /**
     * Returns a key used to optimize repeated calculations of discount
     * factors/forward rate. The calc method for this key returns the 
     * natural logarithm of the discount factor.  The default implementation
     * on this class just takes the log of the result from calling pv().
     */
    virtual YieldCurve::IKey* logOfDiscFactorKey() const;

    // useful methods applicable to any type of zero curve

    /** 
     * Compute discount factor between two dates
     * @param lodate Lower date
     * @param hidate Upper date
     * @return Discount factor between lodate & hidate
     */
    double pv(const DateTime& lodate, const DateTime& hidate) const;

    double fwd(
        const DateTime&           lodate, 
        const DateTime&           hidate,
        const DayCountConvention* dcc,
        int                       basis) const;


    double parSwapRate(
        const DateTime&           startDate, 
        const DateTime&           endDate,
        const MaturityPeriod&     period, 
        const DayCountConvention& dcc,
        const Stub&               stubType,
        bool                      stubAtEnd,
        const BadDayConvention&   accBadDayConv,
        const BadDayConvention&   payBadDayConv,
        const Holiday&            holidays) const ;

    /**
     * Compute PV for cashflow array on curve base date.
     */
    double pv(const CashFlowArray& cfl) const;

    /**
     * Compute value for cashflow array on specified date.
     */
    double fv(const CashFlowArray& cfl, const DateTime& date) const;

protected:
    ZeroCurve(CClassConstSP clazz);

private:
    class LogOfDiscFactorKey;
};


typedef smartConstPtr<ZeroCurve> ZeroCurveConstSP;
typedef smartPtr<ZeroCurve> ZeroCurveSP;
typedef array<ZeroCurveSP, ZeroCurve> ZeroCurveArray;


DRLIB_END_NAMESPACE
#endif // ZERO_CURVE_HPP
