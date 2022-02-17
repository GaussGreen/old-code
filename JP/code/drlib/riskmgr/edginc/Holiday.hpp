//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Holiday.hpp
//
//   Description : Holiday representation
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef HOLIDAY_HPP
#define HOLIDAY_HPP

#include <map>
#include "edginc/DateTime.hpp"
#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE
class WrapperNameCollector;

typedef map<int, int, less<int> > intMap;

/** A class whch defines holidays and weekends and provides a range of functions 
    related to holidays and business days */
class RISKMGR_DLL Holiday: public MarketObject {
public:
    static CClassConstSP const TYPE;
    friend class HolidayHelper;

    ~Holiday();

    /** Holidays are essentially immutable objects so clone just returns this
        - more or less */
    virtual IObject* clone() const;
    
    /** Constructs a Holiday with full flexibility
          @param name        identifier for printing
          @param holidays    non-business days (may be null)
          @param useWeekends weekends (Saturday & Sunday) are 
        not business days */
    Holiday(const string& name, 
            const DateTimeArray& holidays,
            bool useWeekends);

    /** Returns true if the suppplied hols are same as this hol */
    bool equals(const Holiday* hols) const;

    /** Returns a simple hash code value: the hash code of the holiday name */
    virtual int hashCode() const;

    /** Create a Holiday where only weekends are non-business days */
    static Holiday* weekendsOnly();
 
    /** Create a Holiday where all days are business days */
    static Holiday* noHolidays();

    /** is a date a business day? */
    bool isBusinessDay(const DateTime& date) const;
    
    /** is a date a holiday (includes weekends) */
    bool isHoliday(const DateTime& date) const;

    /** add a number of business days to a date */
    DateTime addBusinessDays(const DateTime& from, int busDays) const;

    /** return copy of holiday name */
    string getName() const;

    /** is every day a business day? */
    bool isAllWorkingDays() const; 

    /** Calculates the number of business days between two dates (FROM & TO).
        
       Algorithm:
       1. if FROM = TO, the result is 0 
       2. if FROM < TO, the result is the number of business days in the CLOSED
       interval of [FROM+1,TO]
       3. if FROM > TO, the result is negated number of business days in the 
          CLOSED interval of [TO,FROM-1] */
    int businessDaysDiff(const DateTime&  fromDate,
                         const DateTime&  toDate) const;

 /* returns two dates defining the next holiday on or after/before
    start date.  When fwdSearch is TRUE the two dates returned are
    holStart which is MAX(startDate, start of holiday) and holEnd which
    is the date of the last day in the holiday. When
    fwdSearch is FALSE then the algorithm works as if you travelling backwards
    in time so that holStart will be MIN(startDate, last day in holiday) etc.
    NB The time of day in the dates returned is undefined */
    void findNextHol(const DateTime&   startDate, 
                     bool              fwdSearch, 
                     DateTime&         holStart,
                     DateTime&         holEnd,  
                     bool&             foundHol) const;

    /** Calculates the number of business days for the given year and whether 
        or not it is a leap year  */
    void numBusDaysInYear(int   year, 
                          int&  numBusDays, 
                          bool& isLeapYear) const;

    /** used for translating between EDR and ALIB holiday objects only.
        Simply returns the array of holiday dates and whether
        weekends are included as bad days or not */
    DateTimeArrayConstSP toALIB(bool& areWeekendsHolidays) const;

    /**Checks dates are in order. Removes dates at w/ends if counting w/ends */
    void validatePop2Object();

private:
    Holiday();
    Holiday(const Holiday &rhs);
    Holiday& operator=(const Holiday& rhs);

    // fields - essentially Holiday object is immutable so the same object
    // can be reused multiple times. The cache obviously is not immutable
    // but a value once written in there is never changed. If we were doing
    // multithreaded then we would need to put a lock around access to the 
    // cache
    const string                name;
    const bool                  useWeekends;
    const DateTimeArrayConstSP  holidays;
    mutable intMap              cache;  // $unregistered

    static void acceptWrapperNameCollector(Holiday* holiday,
                                           WrapperNameCollector* collector);

    typedef vector<DateTime>::const_iterator HolIter;
    typedef vector<DateTime>::const_reverse_iterator HolRevIter;

    void findNextWeekEndHol(const DateTime&   startDate, 
                            bool              fwdSearch, 
                            DateTime&         holStart,
                            DateTime&         holEnd,  
                            bool&             foundHol) const;

   /** Compute the number of week-day holidays starting at the input holiday 
       index and continuing until the holiday list is at or beyond toDate. 
       Week-day holidays are holidays which occur on a non-weekend day. */
    int  calcNumWeekdayHolidays (
        const DateTime&   toDate,        /* End date.                */
        HolIter&          startPt,       /* Starting holiday index.  */ 
        const HolIter&    endPt) const;  /* ending holiday index.  */ 

    int  calcNumWeekdayHolidays (
        const DateTime&   toDate,        /* End date.                */
        HolRevIter&       startPt,       /* Starting holiday index.  */ 
        const HolRevIter& endPt) const;  /* ending holiday index.  */ 

    /** This function forwards a date by given number of business days under 
        the following conditions: 

        (1) Saturdays and Sundays are both holidays.
        (2) The direction is forwards in time.                            */ 
    DateTime forwardStandardWeekends ( 
        const DateTime& fromDate,           /* (I) start date             */ 
        int             numBusDays) const;  /* (I) how many business days */

    /** This function handles the special case where the number of business 
        days to forward by is very large. */ 
    DateTime forwardNonStandardWeekends (
        const DateTime& fromDate,              /* (I) start date */ 
        int             numBusDaysLeft,        /* (I) abs. num. bus. days */ 
        int             direction,             // (I) +1=forwards, -1=backwards
        int             busDaysPerWeek) const; // (I) num. bus. days per week 
};

typedef smartConstPtr<Holiday> HolidayConstSP;
typedef smartPtr<Holiday> HolidaySP;
#ifndef QLIB_HOLIDAY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<Holiday>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<Holiday>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<Holiday>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<Holiday>);
#endif

// support for wrapper class
typedef MarketWrapper<Holiday> HolidayWrapper;
#ifndef QLIB_HOLIDAY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL MarketWrapper<Holiday>);
EXTERN_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetInLine<HolidayWrapper>(HolidayWrapper* t));
EXTERN_TEMPLATE(void RISKMGR_DLL FieldSetInLine<HolidayWrapper>(HolidayWrapper* t,
                                                    IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL MarketWrapper<Holiday>);
INSTANTIATE_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetInLine<HolidayWrapper>(
                         HolidayWrapper* t));
INSTANTIATE_TEMPLATE(void RISKMGR_DLL FieldSetInLine<HolidayWrapper>(HolidayWrapper* t,
                                                         IObjectSP o));
#endif

DRLIB_END_NAMESPACE
#endif
