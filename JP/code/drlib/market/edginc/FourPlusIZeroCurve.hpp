//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ZeroCurve.hpp
//
//   Description : Zero curve
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef FOURPLUSI_ZEROCURVE_HPP
#define FOURPLUSI_ZEROCURVE_HPP

#include "edginc/ZeroCurve.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include <string>

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** 
 * 4+i style zero curve.
 *
 * This class is not meant to be used directly by clients - it is intended to
 * be accessed via a yield curve which acts as a proxy onto the methods 
 * detailed below.
 */
class MARKET_DLL FourPlusIZeroCurve : public ZeroCurve {
public:
    static CClassConstSP const TYPE;
    friend class FourPlusIZeroCurveHelper;

    /** set up empty shell */
    FourPlusIZeroCurve(const DateTime& valueDate, int length = 0 /* unused */ );
	FourPlusIZeroCurve(const DateTime& valueDate, string zeroInterpMethod, int length = 0 /* unused */ );

    /** Constructs a zero curve from raw zero dates/rates */
    FourPlusIZeroCurve(const DateTime&           valueDate, 
              DayCountConventionConstSP dayCountConv,
              int                       basis,
              const DateTimeArray&      dates,
              const DoubleArray&        rates);

    virtual ~FourPlusIZeroCurve();

    virtual void validatePop2Object();

    /** Convert a this risk free zero curve to a risky zero curve */
    void CDSriskyZeroCurve(double recovery, const CashFlowArray& defaultRates);

    /** Adds simple-interest money market bond to ZCurve */
    void addMoneyMarket(const DateTimeArray&      dates,
                        const DoubleArray&        rates,
                        int                       numMoneyMarket,
                        const DayCountConvention* dcc);


    /** Adds a discount factor (at a specified date) to a ZCurve */
    void addDiscountFactor(const DateTime& date,
                           double          disc);

    /** Adds a zero rate to a ZCurve */
    void addZeroRate(const DateTime& date,
                     double          rate);
    
    /** Adds a general zero rate to a ZCurve */
    void addGeneralRate(const DateTime&           date,
                        double                    rate,
                        int                       basis,
                        const DayCountConvention* dcc);

    /** Adds a zero rate and discount to a ZCurve */
    void addRateAndDiscount(const DateTime& date,
                            double          rate,
                            double          disc);

    /** Adds information represented by a list-of-cash-flows to a zero
        curve.  Any cash flows which are already covered by the zero curve are
        discounted at rates derived from the zero curve.  Cash flows beyond the
        zero curve imply discount factors, which are added to the zero curve.  If
        there is more than one such cash flow, several points are added to the curve,
        which are calculated by using an interative root-finding secant method,
        where the discount factor for the last cash flow is guessed (and the other
        discount factors are implied by interpolation) where the current price = 
        net present value of all the cash flows.
    */
    void addCashFlows(
        const CashFlowArray& cfl,     // cash flows to add to ZCurve
        double               price,   // current price of cash-flows
        const DateTime&      date);   // date to add last point at


    /** Adds a forward rate to a ZCurve - requires stub rate in place already */
    void addForward(
        const DateTime& dateStart,    // starting date of forward
        const DateTime& dateEnd,      // ending date of forward
        double          discount);    // forward discount factor

    /** Adds a strip of futures to a ZCurve using a 3M stub */
    void addFutures(
        const DateTimeArray&      dateStart,  // start of future
        const DateTimeArray&      dateEnd,    // end date of future
        const DoubleArray&        rates,      // e.g. .06 for 9400
        double                    vol,        // of future
        const DayCountConvention* dcc);       // of money mkt

    /** add stub for futures based on 3M rate */
    void add3mFutureStub(const DateTimeArray&      futureStart,
                         const DateTimeArray&      futureEnd,
                         const DoubleArray&        futRates,
                         double                    rate3M,
                         double                    vol3M,
                         const DayCountConvention* dcc);

    /** Adds points to a zero curve from another zero curve, but only
        those dates before first date of the other curve */
    void addPrefixCurve(const FourPlusIZeroCurve* zc);

    /** Adds points to a zero curve from another zero curve, but only
        those dates after last date of the other curve */
    void addSuffixCurve(const FourPlusIZeroCurve* zc);
    
    /** Converts zc-style rate into a discount factor */
    double computeDiscount(const DateTime& date,
                                   double          rate) const;

    /** Calculates discount factor for a date */
    double discountFactor(const DateTime& date) const;

    /** Calculates an interpolated rate from a ZCurve at some date */
    virtual double zeroCouponRate(const DateTime& date) const;

    /** Calculates an interpolated rate from a ZCurve at some date 
        with option on how to handle dates past end of curve - either
        go flat or linearly extrapolate */
    double interpolate(const DateTime& date, bool extendFlat) const;

    /** Calculates coupon rate for a swap from a zero curve.
        algorithm:
        We know the par price of the swap should be equal to the sum of the
        present value of the principal plus present value all of the coupons.
        Note: coupon[i] = couponFrequency * yearFraction(cDate[i] - cDate[i-1])

        par price = PV(principal) + sum of PV(coupon[i])
        ""     =      ""        + sum of PV(couponRate * yearFrac[i])
        ""     =      ""        + sum of discount[i] * couponRate * yearFrac[i]
        ""     =      ""        + couponRate * (sum of discount[i]*yearFrac[i])
        
        solving for couponRate:
        couponRate=(par price-PV(principal))/( sum of discount[i]*yearFrac[i])
    */
    double swapRate(
        const DateTime&           startDate,
        const DateTime&           matDate,
        bool                      stubAtEnd,
        int                       count,           // interval = count periods
        const string&             period,          // e.g. Y, M, W, D
        const DayCountConvention* dcc) const ;

    /** how long is the curve ? */
    int length() const;
    
    const DateTime& getBaseDate() const;

    /** when is first date for which we have genuine information? */
    const DateTime& firstDate() const;

    /** when does it end ? */
    const DateTime& endDate() const;
    
    /**
     * Returns a key used to optimize repeated calculations of discount
     * factors/forward rate. The calc method for this key returns the 
     * natural logarithm of the discount factor.
     */
    virtual YieldCurve::IKey* logOfDiscFactorKey() const;

    /** strip out the rates and dates */
    CashFlowArraySP getRatesAndDates()const;

    /** strip out the dates */
    DateTimeArray getDates()const;

protected:
    FourPlusIZeroCurve(CClassConstSP clazz);
    FourPlusIZeroCurve(CClassConstSP clazz, const DateTime& valueDate);

    DateTime      valueDate;

    // if ZeroCurve is going to be derived from, this had better be virtual too 
	/* encapsulate the fast linear interpolation in its own function 
	   so that other implementations can override it */
	virtual double doFastInterp(int date,
						        const int idx) const;

private:
    FourPlusIZeroCurve();
    FourPlusIZeroCurve(const FourPlusIZeroCurve& rhs);
//    FourPlusIZeroCurve& operator=(const FourPlusIZeroCurve& rhs); // undefined to prevent usage

    double fastInterpRate(int        date,       /* (I) interp date */
                          int*       idx) const; /* (M) where to search from */

	double calcRateInterp( const DateTime& date,
							const DateTime& loDate, 
							double          loRate, 
							const DateTime& hiDate, 
							double          hiRate) const;

    void zapDate(const DateTime& date);

    // fields
    DateTimeArray dates;
    DoubleArray   rates;
    DoubleArray   discount;
    int           basis;
    DayCountConventionConstSP dayCountConv;
	string zeroInterpolationMethod;

    // used internally
    mutable int loBound;
    mutable int hiBound;

    class FourPlusILogOfDiscFactorKey;
    friend class FourPlusILogOfDiscFactorKey;
};


typedef smartPtr<FourPlusIZeroCurve> FourPlusIZeroCurveSP;


DRLIB_END_NAMESPACE
#endif // FOURPLUSI__ZEROCURVE_HPP
