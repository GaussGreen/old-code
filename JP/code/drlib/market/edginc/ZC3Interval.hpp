//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3Interval.hpp
//
//   Description : Helper class for bootstrapping ALIB style zero curve
//
//   Author      : Richard Appleton
//
//   Date        : 9th May 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_INTERVAL_HPP
#define ZC3_INTERVAL_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"


DRLIB_BEGIN_NAMESPACE


class ZC3Interval;
typedef smartPtr<ZC3Interval> ZC3IntervalSP;
typedef array<ZC3Interval>    ZC3IntervalArray;


class MARKET_DLL ZC3Interval : public CObject
{
public:
    /**
     * Perform interval merging.
     */
    static void mergeAndSort(ZC3IntervalArray& intervals);

    static bool lessThan(const ZC3Interval& first, const ZC3Interval& second);

    static CClassConstSP const TYPE;
    
    ZC3Interval();   // for STL use only
    ZC3Interval(const CClassConstSP& clazz);
    ZC3Interval(const DateTime& start, const DateTime& end, double rate, const CClassConstSP& clazz = TYPE);

    string toString() const;

    bool            isForSamePeriodAs(const ZC3Interval& other) const;
    const DateTime& getStartDate() const;
    const DateTime& getEndDate() const;
    double          getContRate() const;

    void            setEndDate(const DateTime& date);
    void            setContRate(double rate);

    // inclusive of start date, exclusive of end date
    bool startsWithin(const DateTime& start, const DateTime& end) const;

protected:
    DateTime startDate;
    DateTime endDate;
    double   contRate;   // actually an integral over [startDate, endDate]

private:
    static void load(CClassSP& clazz);
};


class MARKET_DLL ZC3Adjustment : public ZC3Interval
{
public:
    static CClassConstSP const TYPE;

    ZC3Adjustment();     // for STL use only
    ZC3Adjustment(const DateTime& start, const DateTime& end, double rate, double discount);

    bool   isAbsolute() const;
    double getDiscount() const;
    void   modify(double  contRate);

private:
    static void load(CClassSP& clazz);

    double discount;
};


typedef smartPtr<ZC3Adjustment> ZC3AdjustmentSP;
typedef array<ZC3Adjustment>    ZC3AdjustmentArray;


/**
 * Specializations of arrayObjectCast (needed as arrays are not array of pointers)
 */
template <> class MARKET_DLL arrayObjectCast<ZC3Interval>
{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ZC3Interval& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ZC3Interval& value);

    /** Turns the IObjectSP into a ZC3Interval */
    static ZC3Interval fromIObject(IObjectSP& value);
};


template <> class MARKET_DLL arrayObjectCast<ZC3Adjustment>
{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ZC3Adjustment& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ZC3Adjustment& value);

    /** Turns the IObjectSP into a ZC3Adjustment */
    static ZC3Adjustment fromIObject(IObjectSP& value);
};


DRLIB_END_NAMESPACE
#endif // ZC3_INTERVAL_HPP
