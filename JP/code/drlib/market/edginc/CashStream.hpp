//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : CashStream.hpp
//
//   Description : Stream of fixed or floating cashflows.
//
//   Author      : Richard Appleton
//
//   Date        : 25th April 2005
//
//----------------------------------------------------------------------------

#ifndef CASH_STREAM_HPP
#define CASH_STREAM_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/RiskFreeCashFlow.hpp"


DRLIB_BEGIN_NAMESPACE

class IRVolBase;
class ZeroCurve;
class BadDayConvention;
class DayCountConvention;
class Holiday;
class Stub;


class MARKET_DLL CashStream : public CObject
{
public:
    static CClassConstSP const TYPE;

    CashStream(int size = 0);
    ~CashStream();

    void add(const CashFlowArray& cfl);
    void add(RiskFreeCashFlow* point);
    void addFixedFlow(const DateTime& date, double amount);
    void addFloatFlow(
        const DateTime&           date, 
        double                    notional, 
        double                    spread,           // additive
        const DayCountConvention* dcc, 
        const DateTime&           accrueStartDate, 
        const DateTime&           accrueEndDate, 
        const DateTime&           rateStartDate, 
        const DateTime&           rateEndDate);

    bool empty() const;
    int  size() const;
    const RiskFreeCashFlow& operator[](int i) const;

    /**
     * Sort cash stream by date.
     */
    void sort();

    /**
     * Computes last effective date for a cash stream. This is the last date
     * from a zero curve that would affect the value of the cash stream.
     */
    DateTime maturityDate() const;

    /**
     * Enumerates critical dates; the critical points are the set of points 
     * that are required to correctly reprice the benchmark instruments.
     */
    DateTimeArraySP getCriticalDates() const;

    /**
     * Convert cash stream to a cashflow list.  Note that cashflows on the
     * same date are combined, and the cashflow list is returned in date order.
     */
    CashFlowArraySP flows(const ZeroCurve& curve, const IRVolBase* volModelIR);

    /**
     * Adds a vanilla floating leg to a cash stream.
     */
    void addFloatLegVanilla(
        double                    notional,
        double                    spread,
        bool                      isAdditive,
        const DateTime&           startDate,
        const DateTime&           valueDate,
        const DateTime&           rollDate,
        const MaturityPeriod&     interval,
        const DateTime&           maturityDate,
        const Stub&               stubType,
        bool                      stubAtEnd,
        bool                      subtractInitial,
        bool                      addFinal,
        const BadDayConvention&   accBadDayConv,
        const BadDayConvention&   payBadDayConv,
        const BadDayConvention&   resetBadDayConv,
        bool                      firstPaymentFixed,
        double                    firstPaymentRate,
        const DayCountConvention& dayCountConv,
        const Holiday&            holidays);

private:
    static void load(CClassSP& clazz);

    double stubPaymentAmount(
        const DateTime&           prevCouponDate,
        const DateTime&           nextCouponDate,
        const DateTime&           stubStart,
        const DateTime&           stubEnd,
        double                    rate,
        const DayCountConvention& couponDayCountConv
        // stub type hard-coded to GTO_STUB_SIMPLE
        ) const;

    RiskFreeCashFlowArray cashflows; // $unregistered
};


typedef smartPtr<CashStream>      CashStreamSP;
typedef array<CashStream>         CashStreamArray;


/**
 * Specializations of arrayObjectCast (needed as arrays are not array of pointers)
 */
template <> class MARKET_DLL arrayObjectCast<CashStream>
{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const CashStream& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(CashStream& value);

    /** Turns the IObjectSP into a CashStream */
    static CashStream fromIObject(IObjectSP& value);
};


DRLIB_END_NAMESPACE

#endif // CASH_STREAM_HPP
