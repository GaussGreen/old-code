//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : SimpleZeroCurve.hpp
//
//   Description : Zero curve with dates/rate explicitly defined.
//
//   Author      : Richard Appleton
//
//   Date        : 27 June 2005
//
//----------------------------------------------------------------------------

#ifndef SIMPLE_ZERO_CURVE_HPP
#define SIMPLE_ZERO_CURVE_HPP

#include <string>
#include "edginc/DayCountConvention.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/YieldCurve.hpp"


using namespace std;    // string

DRLIB_BEGIN_NAMESPACE


class SimpleZeroCurveHelper;


 /**
  * Zero curve where the dates and rates are explicitly specified, rather than
  * calculated via some bootstrapping process.
  */
class MARKET_DLL SimpleZeroCurve : public ZeroCurve
{
public:
    static CClassConstSP const TYPE;

    /** Construct zero curve from raw zero dates/rates */
    SimpleZeroCurve(
        const DateTime&           baseDate, 
        int                       basis,
        const DayCountConvention* dcc,
        const string&             interpolation,
        const ExpiryArray&        dates,
        const DoubleArray&        values);

    ~SimpleZeroCurve();

    /** Calculates discount factor for a date */
    double discountFactor(const DateTime& date) const;

    /** Calculates an zero coupon rate from a ZCurve at some date */
    double zeroCouponRate(const DateTime& date) const;

    /** how long is the curve ? */
    int length() const;
    
    /** when is first date for which we have genuine information? */
    const DateTime& firstDate() const;

    /** when does it end ? */
    const DateTime& endDate() const;
    
    /** strip out the dates */
    DateTimeArray getDates() const;

    /** strip out the rates and dates */
    CashFlowArraySP getRatesAndDates() const;

    const DateTime& getBaseDate() const;

    /** Returns a key used to optimize repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative 
        year fraction (Act/365F) between the two dates. */
    YieldCurve::IKey* logOfDiscFactorKey() const;

    /** Convert a this risk free zero curve to a risky zero curve */
    // TBD!! remove need for this!! (ALIB curves throw exception if called)
    void CDSriskyZeroCurve(double recovery, CashFlowArray& defaultRates);

private:
    class LogOfDiscFactorKey;
    friend class LogOfDiscFactorKey;
    friend class SimpleZeroCurveHelper;

    SimpleZeroCurve();
    
    void addRate(const DateTime&, double rate);
    void insertValue(const DateTime&, double rate);

    double fastInterpRate(int date, int* idx) const;
    double interpolate(
        const DateTime& date,
        const DateTime& loDate,
        double          loRate,
        const DateTime& hiDate,
        double          hiRate) const;

    // fields
    DateTime                  baseDate;
    int                       basis;
    string                    interpolation;
    DayCountConventionConstSP dayCountConv;
    DateTimeArray             dates;
    DoubleArray               rates;
    DoubleArray               discount; // $unregistered

    // used internally
    mutable int loBound; // $unregistered
    mutable int hiBound; // $unregistered
};


DRLIB_END_NAMESPACE

#endif // SIMPLE_ZERO_CURVE_HPP
