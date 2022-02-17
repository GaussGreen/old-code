//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Interval.hpp
//
//   Description : A Date interval
//
//   Author      : André Segger
//
//   Date        : 05 October 2001
//
//
//----------------------------------------------------------------------------

#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include "edginc/DateTime.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL Interval: public CObject {
public:
    friend class IntervalHelper;
    static CClassConstSP const TYPE;

    //// Constructor - inline for performance
    Interval(const DateTime& startDate, const DateTime& endDate): CObject(TYPE), 
                              startDate(startDate), endDate(endDate) {}
    
    virtual ~Interval(){} // inline for performance;

    /** returns the start date of the interval */
    virtual void validatePop2Object();

    /** returns the start date of the interval */
    const DateTime& getStartDate() const;

    /** returns the end date of the interval */
    const DateTime& getEndDate() const;
    
    /** returns the number of days between the two dates */
    int daysDiff() const;

    /** public as DateTimeArray is an array of structures */
    Interval();

    //// Added to allow array template to be instantiated
    //// Implementation just compares two dates.
    bool operator==(const Interval& rhs) const;

private:
    DateTime startDate;
    DateTime endDate;
};



typedef smartPtr<Interval> IntervalSP;

/** specialisations of arrayObjectCast */
template <> class TOOLKIT_DLL arrayObjectCast<Interval>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const Interval& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(Interval& value);

    /** Turns the IObjectSP into a DateTime */
    static const Interval& fromIObject(IObjectSP& value);
};


typedef array<Interval, Interval>    IntervalArray;
typedef smartPtr<IntervalArray>      IntervalArraySP;
typedef smartConstPtr<IntervalArray> IntervalArrayConstSP;

#ifndef QLIB_INTERVAL_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<Interval _COMMA_ Interval>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<Interval _COMMA_ Interval>);
#endif
typedef Interval              CInterval;
typedef IntervalArray         CIntervalArray;
typedef IntervalArraySP       CIntervalArraySP;
typedef IntervalArrayConstSP  CIntervalArrayConstSP;

DRLIB_END_NAMESPACE

#endif


