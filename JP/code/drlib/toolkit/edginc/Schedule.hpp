//---------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Schedule.hpp
//
//   Description : A schedule of dates & levels you can interpolate on
//
//   Author      : Andrew J Swain
//
//   Date        : 9 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef SCHEDULE_HPP
#define SCHEDULE_HPP

#include <string>
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** A schedule of dates & levels you can interpolate on */

class TOOLKIT_DLL Schedule : public CObject {
public:
    static CClassConstSP const TYPE;

    static const int ADJUST_NONE;
    static const int ADJUST_CONSTANT;
    static const int ADJUST_INTERPOLATE;

    /** how you can interpolate */
    static const string INTERP_NONE;
    static const string INTERP_LINEAR;
    static const string INTERP_STAIRS;

    
    Schedule(const DateTimeArray& dates,
             const DoubleArray&   values,
             const string&        interp);

    Schedule(const Schedule& rhs,
             double          value);

    virtual ~Schedule();

    /** overridden for performance */
    virtual IObject* clone() const;

    /** convenience access to last date */
    DateTime lastDate() const;
    /** convenience access to last value */
    double   lastValue() const;

    /** convenient access to first date */
    DateTime firstDate() const;
    /** convenience access to first value */
    double   firstValue() const;
  
    /** how long is the schedule ? */
    int length() const;

    /** v. basic interpolator */
    double interpolate(const DateTime& date) const;

    double interpolate(const DateTime::Date& date) const;

    bool   interpolateCVBSchedule(int      adjustmentType,
                                  double   redemptionValue,
                                  DateTime bondMaturity,
                                  DateTime interpDate,
                                  double   &level) const;

    /** scale rates by a double */
    /* (wouldn't need such an inane function if we could specify schedules
        as percentage even after start date (gasp!) ) */
    void scale(double factor);

    /** validation */
    virtual void validatePop2Object();

    /** due to a weak validation, empty schedules are allowed to be created
        adding a validation in validatePop2Object() is not the solution as
        many instruments are creating empty schedule without using them 
        This (poor) solution prevents to use empty schedule to prevent edr
        to core */
    void isEmpty(string method) const;

    /** return interp type */
    string getInterp() const;

    /** return date list (deep copy) */
    DateTimeArray getDates() const;

    /** return value list (deep copy) */
    DoubleArray getValues() const;

    /** returns reference to date array */
    const DateTimeArray& getDateArray() const;

    /** returns reference to value array */
    const DoubleArray& getValueArray() const;

    /** Determine if a date is contained inside a schedule */
    bool coversDateRange(const DateTime& lowDate, const DateTime& highDate, bool datesDefineLimit)const;

    /* This tests the rates of a schedule to see if they are all
       greater than or equal to zero. */
    bool isNonNegative()const;
    
    /* This tests the rates of a schedule to see if they are all strictly greater than zero. */
    bool isPositive() const;

    /** This tests the rates of a schedule to see if they are all
        almost equal to zero. */
    bool isAlmostZero()const;

    /** Determines if all rates in the schedule are equal */
    bool isFlat()const;

    bool timesAreAll(int timeToCompare) const;

    /** This function returns the next schedule date for a given (value) date and a maturity date 
        up to which the schedules will automatically extend after the last date in the schedule's
        date list */
    bool getNextDate(const DateTime& fromDate,
                     const DateTime& bondMaturity, 
                     DateTime& nextDate) const;
    
    /** returns a subset of the schedule that covers the supplied date range.
        If INTERP_NONE, returns any values that lie inside the given dates
        otherwise retuns a value per day for each day between supplied values.
        Advances from lower bound in 1 day increments  using best guess at
    time of day to give most complete results */
    CashFlowArray* subset(const DateTime& lowDate, const DateTime& highDate) const;

    /** check the dates in schedule are same to refDates*/
    bool isSameSchedule(const DateTimeArray refDates, bool compTime = true) const;

private:
    friend class ScheduleHelper;

    Schedule();
    Schedule(const Schedule &rhs);
    Schedule& operator=(const Schedule& rhs);

    DateTimeArrayConstSP dates;
    DoubleArrayConstSP   values;
    string               interp;
};

typedef smartConstPtr<Schedule> ScheduleConstSP;
typedef smartPtr<Schedule> ScheduleSP;

#ifndef QLIB_SCHEDULE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<Schedule>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Schedule>);
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<ScheduleSP>(ScheduleSP* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<ScheduleSP>(ScheduleSP* t,
                                                  IObjectSP o));
#endif

typedef array<ScheduleSP, Schedule> ScheduleArray;
typedef smartPtr<ScheduleArray> ScheduleArraySP;

DRLIB_END_NAMESPACE
#endif


