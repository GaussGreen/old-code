//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DateTime.hpp
//
//   Description : Date-time representation
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef DATETIME_HPP
#define DATETIME_HPP

#include "edginc/Object.hpp"
#include "edginc/Array.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/CombinableResult.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/VirtualDestructorBase.hpp"
using namespace std;

DRLIB_BEGIN_NAMESPACE

#define EDG_WEEKEND_SATURDAY    0x0020
#define EDG_WEEKEND_SUNDAY      0x0040
#define EDG_WEEKEND_STANDARD    (EDG_WEEKEND_SATURDAY | EDG_WEEKEND_SUNDAY)

class DateTime; // predeclaration in order to declare array

typedef array<DateTime> DateTimeArray;
typedef smartPtr<DateTimeArray> DateTimeArraySP;
typedef smartConstPtr<DateTimeArray> DateTimeArrayConstSP;


#ifndef QLIB_DATETIME_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<DateTime _COMMA_ DateTime>);
#else
// normally this would go in the source file but MSVC warns that it has already
// been instantiated if you do that - unclear why
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<DateTime _COMMA_ DateTime>);
#endif

/** DateTimeArray is an array of DateTime structures - not an array of
    pointers. Note that this approach forces the no argument DateTime
    constructor to be public */
class TOOLKIT_DLL DateTime: public CObject {
public:
    friend class DateTimeHelper;
    static CClassConstSP const TYPE;
    static const int TIME_IN_DAY;
    static const int TIME_MAX;
    static const int TIME_MIN;
    static const int START_OF_DAY_TIME;
    static const int BEFORE_EX_DIV_TIME;
    static const int EX_DIV_BEFORE_START_TIME;
    static const int END_OF_DAY_TIME;
    static const int EFFECTIVE_TIME_IN_DAY;
    static const int DAYS_PER_YEAR;
    static const int TIME_IN_YEAR; 
    static const int DAYS_PER_WEEK;
    static const int VOL_TIME_IN_DAY;
    static const int VOL_TIME_IN_YEAR;
    
    static const int SATURDAY;
    static const int SUNDAY;
    static const int MONDAY;
    static const int TUESDAY;
    static const int WEDNESDAY;
    static const int THURSDAY;
    static const int FRIDAY;

    enum DayOfWeek { Saturday, Sunday, Monday, Tuesday, Wednesday, Thursday, Friday };
    
    static const string START_OF_DAY;
    static const string END_OF_DAY;
    static const string BEFORE_EX_DIV;
    static const string EX_DIV_BEFORE_START;

    static bool  isLeapYear(long year);

    //// Constructor - inline for performance
    DateTime(int date, int time);

    
    //// Copy Constructor - inline for performance
    DateTime(const DateTime& dt);

    /** takes date in dd-mmm-yyyy format, time is either
        START_OF_DAY, END_OF_DAY or hh-mm-ss format
    */
    DateTime(const string& date, const string& time);

    DateTime& operator=(const DateTime& rhs);

    virtual ~DateTime();
    
    /** Overridden for performance */
    virtual IObject* clone() const;

    /** returns a hashcode for the object based upon the date and time.
        Overridden for performance */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one. Overridden
        for performance */
    virtual bool equalTo(const IObject* obj) const;
    
    /** Returns true if the datetime has not been initialsed ie it has been
        built through the default constructor */
    bool empty() const;

    bool equals(const DateTime& date, bool compTime=true) const;

    int compareTo(const DateTime& date) const;

    bool isGreater(const DateTime& date) const;
    bool isGreaterOrEqual(const DateTime& date) const;
    bool isLess(const DateTime& date) const;
    bool isLessOrEqual(const DateTime& date) const;

    // for weak (ignoring time) comparison
    static bool isLessWeak(const DateTime& date1, const DateTime& date2);

    /** Returns MAX(this, date) */
    const DateTime& max(const DateTime& date) const;

    /** Returns MAX(this, date1, date2) ie maximum of all three dates */
    const DateTime& max(const DateTime& date1, const DateTime& date2) const;

    /** Returns MIN(this, date) */
    const DateTime& min(const DateTime& date) const;

    bool operator < (const DateTime& date) const;
    bool operator > (const DateTime& date) const;

    bool operator <= (const DateTime& date) const;
    bool operator >= (const DateTime& date) const;
    
    bool operator == (const DateTime& date) const;
    bool operator != (const DateTime& date) const;

    /** convertToYearFracsUsingRefDate applies yearFrac to every element of
    gDateArray using 'this' as the reference date and returns DoubleArraySP */
    template<class TContainer4DateTime>
        DoubleArraySP convertToYearFracsUsingRefDate(
        const TContainer4DateTime & gDateTimeContainter) const
    {
        const int _size = gDateTimeContainter.size();
        DoubleArraySP _result(new DoubleArray(_size));
        const DateTime _this = (*this); // to save time
        typename TContainer4DateTime::const_iterator 
            _cit = gDateTimeContainter.begin();
        int s=0; for (; s<_size; s++, _cit++)
        {
            (*_result)[s] = _this.yearFrac(*_cit);
        }
        return _result;
    }
    /* Returns a DateTime with the number of days increased by offset */
    DateTime rollDate(const int daysOffset) const;

    /* Returns a DateTime with the number of months increased by offset */
    DateTime rollDateInMonths(const int monthsOffset) const;

    /* Attempts to answer the nth dayOfWeek (e.g. 3rd Thursday) of the month
       of this DateTime */
    DateTime nthWeekdayOfMonth(int n, int dayOfWeek) const;

    /* Answers the last dayOfWeek of the month of this DateTime */
    DateTime lastWeekdayOfMonth(int dayOfWeek) const;

    /** Returns the number of days between this and loDate. Equivalent to
        this.getDate() - loDate.getDate() */
    int daysDiff(const DateTime& loDate) const;

    int getDate() const;
        
    int getTime() const;

    /** returns the day of the week with 0 = Saturday ... 6 = Friday */
    int getWeekday() const;
    
    /** is day on a weekend */
    bool isWeekend() const;

    /** Defined such that if dt is a DateTime, 
        dt.isWeekend() == DateTime::isWeekend(dt.getDate()) */
    static bool isWeekend(int date);

    /** is day in a leap year ? */
    bool isLeap() const;

    /** is day last day in month ? */
    bool isEndOfMonth() const;

    /** returns the last day of the month that the date falls in.
        Time of day remains the same. ignoreLeapYears = true sets
        to February 28th in leap years */
    DateTime returnEndOfMonth(bool ignoreLeapYears) const;

    /** returns the first day of the month that this date falls in. */
    DateTime firstOfMonth() const;

    /** Given a dateTime, if the time is < EDG_START_OF_DAY_TIME
        makes the time = EDG_START_OF_DAY_TIME.
        If the time is > EDG_END_OF_DAY_TIME 
        makes the time = EDG_END_OF_DAY_TIME. */
    void moveIntoDay();


    /** Calculates simple year fraction between 2 dates. Assumes 365F */
    double yearFrac(const DateTime &dateTo) const;

    /** returns string representation of date - good for error messages */
    string toString() const;

    /** format a date in DD-MMM-YYYY format */
    static string dateFormat(int date);

    /** format a time either as SOD, EOD or in hh:mm:ss format */
    static string timeFormat(int time);

    /** convert hh:mm:ss into a time */
    static int timeConvert(string time);

    /** to do: add debug field showing day of the week */

    /** public as DateTimeArray is an array of structures */
    DateTime();

    /** Interval class allows DateTimes to be subtracted. The resulting
        interval can then be added to a DateTime. Note no type support */
    class TOOLKIT_DLL Interval{
    public:
        friend class DateTime;
        //// comparison methods
        bool isGreaterOrEqual(const Interval& interval) const;
        bool isLess(const Interval& interval) const;
        bool isGreater(const Interval& interval) const;
        
    private:
        Interval(int dateDiff, int timeDiff);
        int dateDiff;
        int timeDiff;
    };

    /** Returns a date-time midway between this date and the supplied date */
    DateTime midPoint(const DateTime& date) const;

    /** Returns the interval between this date and date ie "this - date" */
    Interval subtract(const DateTime& date) const;

    /** Adds given interval to date */
    DateTime add(const Interval& interval) const;


    /** Calculates the expected end date of a list of future dates which 
        are weighted */
    static const DateTime expectedFutureDate(
        const DateTime& today,
        const array<DateTime>& futDates,
        const DoubleArray& futWeights);

    class TOOLKIT_DLL MonthDayYear {
        friend class DateTime;
    public:
        int day;
        int month;
        int year;

        DateTime toDateTime(int time = START_OF_DAY_TIME);
        int wholeMonthsBetween(MonthDayYear mdyTo);
        MonthDayYear(int day, int month, int year);

        /** returns the number of days in this date's month */
        int daysInMonth(bool ignoreLeapYears = false) const;

    private:
        void normalize();
    };

    MonthDayYear toMDY() const;

    /** Creates DateTime from IR (Alib) style date ie yyyymmdd. Hard coded to
        start of day */
    static DateTime fromIrDate(long irDate);
    
    /** converts to an IR (Alib) style date ie yyyymmdd */
    long toIrDate(void) const;

    /** Returns the year (eg 2003) in which this date falls */
    int getYear() const;

    /** Returns DateTime for Dec 31st EOD in supplied year */
    static DateTime endOfYear(int year);

    /** Counts the number of leap year days between a start date
        (inclusive) and an end date (exclusive).  Needed for
        Japanese Government Bonds and some floating rate notes
    */
    static int countLeapYearDays(const DateTime& d1, const DateTime& d2);

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object in from reader  */
    virtual void import(Reader::Node* elem, Reader* reader); 
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Validates that supplied dates are increasing order. If failIfEmpty is
        true, then throws exception if dates is empty. The description string
        is used in the exception messages should describe what the dates 
        represent */
    static void ensureIncreasing(const DateTimeArray& dates,
                                 const string&        description,
                                 bool                 failIfEmpty);
    /** As above but takes a const char* */
    static void ensureIncreasing(const DateTimeArray& dates,
                                 const char*          description,
                                 bool                 failIfEmpty);

    /** Validates that supplied dates are strictly increasing order. 
        If failIfEmpty is true, then throws exception if dates is empty. 
        The description string is used in the exception messages should
        describe what the dates represent */
    static void ensureStrictlyIncreasing(const DateTimeArray& dates,
                                         const string&        description,
                                         bool                 failIfEmpty);
    
    /** As above but takes a const char* */
    static void ensureStrictlyIncreasing(const DateTimeArray& dates,
                                         const char*          description,
                                         bool                 failIfEmpty);

    /** Returns the number of dates strictly greater than this date. Assumes
        dates are in increasing order */
    int numFutureDates(const DateTimeArray& allDates) const;

    /** Returns an array of dates strictly greater than this date. Assumes
        dates are in increasing order */
    DateTimeArray getFutureDates(const DateTimeArray& allDates) const;

    /** Returns an array of dates greater or equal to this date. Assumes
        dates are in increasing order */
    DateTimeArray getNotPastDates(const DateTimeArray& allDates) const;

    /** Returns the number of dates <= this date. Assumes
        dates are in increasing order */
    int numPastDates(const DateTimeArray& allDates) const;

    /** Returns an array of dates <= this date. Assumes
        dates are in increasing order */
    DateTimeArray getPastDates(const DateTimeArray& allDates) const;

    /** Returns an array of dates from allDates that are bounded 
        inclusively by start and end dates.  start has to be strictly
        smaller than end. */
    static DateTimeArray getInclusiveDates(const DateTime & start,
                                           const DateTime & end,
                                           const DateTimeArray & allDates);

    /** Creates an array of integers of length dates.size()+1. Each
        element in the array is the offset to the index of the array
        which corresponds to the next date in datesForMap. Eg
        Dates:              A B C D E F G H I J
        datesForMap:          B     E       I
        map:                1 0 2 1 0 3 2 1 0 1 0
        The output parameter isTrivial is true if datesForMap coincides with
        the dates. An exception is thrown if any of the dates
        in datesForMap are not in the supplied list of dates. The
        supplied dates must be in order. Typical use of such a map, dateMap:
        for (iStep += dateMap[iStep]; iStep < endIdx; 
        iStep++, iStep += dateMap[iStep]){
        ...
        }
        Note the fact that the dateMap is of length dates.size()+1 means that
        the increment statement "iStep++, iStep += dateMap[iStep]" is safe */
    static IntArray  createMapping(const DateTimeArray&  dates,
                                   const DateTimeArray&  datesForMap,
                                   bool&                 isTrivial);
    
    /** Creates an array of integers of length datesForOffset.size(). Each
        element in the array is the offset to the index of the dates array
        which corresponds to the next date in datesForOffset. Eg
        Index:              0 1 2 3 4 5 6 7 8 9
        Dates:              A B C D E F G H I J
        datesForOffset:       B     E       I
        result returned:    1 4 8
        Also, lastPastDateIdx is set to the largest idx such that 
        dates[idx] < valueDate.
    */
    static vector<int> createMapping2( 
        const DateTimeArray& dates,          // assumed sorted
        const DateTimeArray& datesForOffset, // assumed sorted
        const DateTime& valueDate,
        int& lastPastDateIdx );

    /** Creates an array of bools of length dates.size(). Each
        element in the array is set to true if the date exists in the 
        datesForMap. Eg
        Dates:              A B C D E F G H I J
        datesForMap:          B     E       I
        map:                F T F F T F F F T F 
        Typical use of such a map is if we want to do averaging over certain days etc
        ...
        } Arrays need not be sorted */

    static vector<bool> createmapping3(
        const DateTimeArray& dates,                          
        const DateTimeArray& datesForMapping );

    /** returns an array of length part.size() giving the index of
        each element of 'part' in 'full'. Both arrays must be ordered
        [strictly increasing], and 'part' must be a subset of
        'full' */
    static vector<int> getIndexes(const DateTimeArray& full,
                                  const DateTimeArray& part);
    /** returns an array of size full.size() giving the index of full[i]
        inside part[] or -1 if this element is out of part().
        Preconditions are the same as for getIndexes.
    */                                  
    static vector<int> getProjection(const DateTimeArray& full,
                                     const DateTimeArray& part);

    /** returns an array of size full.size(), giving an index j for full[i]
        such that: j = min(k: part[k] >= full[i]) if full[i] <= part[last],
        otherwise j = last.  In other words, the index j for full[i] is 
        the index of the smallest element in part that is bigger than or 
        equal to full[i], given full[i] <= part[last].  If full[i] > 
        part[last], then j = last. */
    static vector<int> getCeilingProjection(const DateTimeArray & full,
                                            const DateTimeArray & part);

    /** getFloorProjection is similar to getCeilingProjection except it 
        returns indices j's for full[i] such that: j = max(k: part[k] <= 
        full[i]), if full[i] >= part[first], otherwise j = first. */
    static vector<int> getFloorProjection(const DateTimeArray & full,
                                          const DateTimeArray & part);

    /** Rather nasty function which forces the time of day in each of the dates
        supplied to the specified value */
    static void setTimeOfDay(DateTimeArray& dates, // (M)
                             int            timeOfDay);

    /** Dates must be increasing (but does not have to be strictly). 
        Returns the index of the [last] date which is greater
        or equal than the supplied target date. Returns dates.size()
        if target date > all dates. */
    int findUpper(const DateTimeArray& dates) const;

    /** Dates must be increasing (but does not have to be strictly).
        Returns the index of the [last] date which is less or equal
        than the supplied target date. Returns -1 if target date < all
        dates.  */
    int findLower(const DateTimeArray& dates) const;

	/** Dates must be strictly increasing.
		Returns the index of the date array that is nearest to the date.
		Returns the left index if it is exactly midway between two indices.
		Returns -1 if the array is empty
	*/
	int findNearest(const DateTimeArray& dates) const;

    /** Return iterator to the date.  The search is done on the interval
        [dateStart, dateEnd), ie. excluding dateEnd as in stl */
    DateTimeArray::const_iterator findIterator(
        DateTimeArray::const_iterator dateStart,
        DateTimeArray::const_iterator dateEnd) const;

    /** Returns the index of the first instance of this date in the array
        of dates. Fails if the date does not exist. The supplied dates
        must be weakly increasing (a binary search is used) */
    int find(const DateTimeArray& dates) const;

    /** Same as above find if used dates.begin(), dates.end() */
    int find(DateTimeArray::const_iterator dateStart,
             DateTimeArray::const_iterator dateEnd) const;

    /** Returns the index of the first instance of this date in the array
        of dates. Fails if the date does not exist. The supplied dates
        must be weakly increasing (a binary search is used) 
        NOTE THIS METHOD IGNORES THE TIME*/
    int findIgnoringTime(const DateTimeArray& dates) const;

    /** Returns true if this date is between upper and lower
        - the bounds are inclusive */
    bool within(const DateTime& lower,
                const DateTime& upper) const;

    /** Returns true if the two lists of dates are identical */
    static bool equals(const DateTimeArray& dates1,
                       const DateTimeArray& dates2);

    /** merge n date arrays into one removing duplicate dates. All date arrays
        must be be sorted in increasing order and contain no duplicates
        themselves. May want to review parameter type */
    static DateTimeArray merge(const vector<DateTimeArray>& dtArrayVector);

    /** merge n date arrays into one removing duplicate dates. All date arrays
        must be be sorted in increasing order and contain no duplicates
        themselves. May want to review parameter type */
    static DateTimeArray merge(
        const vector<const DateTimeArray*>& dtArrayVector);

    /** merge 2 date arrays into one removing duplicate dates. Both date arrays
        must be be sorted in increasing order and contain no duplicates
        themselves. May want to review parameter type */
    static DateTimeArray merge(const DateTimeArray& dtArray1,
                               const DateTimeArray& dtArray2);

    //// another variant on the merge theme, same rules apply
    static void merge(
        const vector<DateTimeArray::const_iterator>& dtStart,     // (I)
        const vector<DateTimeArray::const_iterator>& dtEnd,       // (I)
        DateTimeArray&                               mergedDates);// (O)

    //// another variant on the merge theme, same rules apply
    static void merge(const vector<const DateTimeArray*>& dtArrayVector, // (I)
                      DateTimeArray&                      mergedDates);  // (O)

    /** Another merge except those dates <= lowerLimit or > upperLimit are
        omitted */
    static void merge(const DateTime&                     lowerLimit,
                      const DateTime&                     upperLimit,
                      const vector<const DateTimeArray*>& dtArrayVector, // (I)
                      DateTimeArray&                      mergedDates);  // (O)

    //// another variant, mainly when you have many shared DateTimeArrays
    static DateTimeArray merge(const vector<DateTimeArrayConstSP>& v);

    /** What do I do? Why aren't the parameters const? RA to document */
    static DateTimeArraySP subtract(DateTimeArray& first, 
                                    DateTimeArray& second);


    /** Removes duplicate dates in a DateTimeArray provided that the duplicate
        dates are next to each other */
    static void removeDuplicates(DateTimeArray& dates,
                                 bool           ignoreTimeOfDay);

    /** Removes dates which lie outside lowerDate or upperDate */
    static void removeOutliers(DateTimeArray&  dates,
                               const DateTime& lowerDate,
                               const DateTime& upperDate);

    /** Report if first is a subset of second */
    static bool isSubset(const DateTimeArray& full,
                         const DateTimeArray& part);
                         
    /** Takes an array which is possibly not in order and with duplicates. In-place sorts and removes duplicates and returns itself. Sorting will use "DateTime::operator < ()". Worst case complexity is the same as of sort (i.e. O(n log n)) */
    static DateTimeArray& doSortUniq(DateTimeArray& dates);

    /** Date portion of DateTime - may want to make a top level class. 
        Currently exists for addin purposes */
    class TOOLKIT_DLL Date: public CObject{
    public:
        static CClassConstSP const TYPE;
        friend class DateTimeHelper;
        /** returns string representation of date */
        string toString() const;
        /** write object out to writer */
        virtual void write(const string& tag, Writer* writer) const;
        /** populate an empty object in from reader  */
        virtual void import(Reader::Node* elem, Reader* reader); 
        /** write object out in 'output' format - ie suitable for comparing
            regression files with */
        virtual void outputWrite(const string& linePrefix,
                                 const string& prefix, ostream& stream) const;
        //// construct object using date as an integer
        Date(int date);
        //// same as empty constructor to DateTime
        Date();

        //// Copy Constructor - inline for performance
        Date(const Date& dt);

        //// returns date as an integer
        int getDate() const;
        /** returns true if dates are the same */
        bool equals(const Date& date) const;
        /** returns true if this date is greater than date supplied */
        bool isGreater(const Date& date) const;
        /** returns true if this date is greater than or equal to date
            supplied */
        bool isGreaterOrEqual(const Date& dateToCompare) const;
    private:
        //// fields
        int date;
    };
    typedef smartPtr<Date> DateSP;
    typedef array<DateSP, Date> DateArray;
    typedef smartPtr<DateArray> DateArraySP;

    /** Time portion of DateTime - may want to make a top level class. 
        Currently exists for addin purposes */
    class TOOLKIT_DLL Time: public CObject{
    public:
        static CClassConstSP const TYPE;
        friend class DateTimeHelper;
        /** returns string representation of date */
        string toString() const;
        /** write object out to writer */
        virtual void write(const string& tag, Writer* writer) const;
        /** populate an empty object in from reader  */
        virtual void import(Reader::Node* elem, Reader* reader); 
        /** write object out in 'output' format - ie suitable for comparing
            regression files with */
        virtual void outputWrite(const string& linePrefix,
                                 const string& prefix, ostream& stream) const;
        //// construct object using time as an integer
        Time(int time);
        //// returns time as an integer
        int getTime() const;
    private:
        //// fields
        int time;
    };
    typedef smartPtr<Time> TimeSP;
    typedef array<TimeSP, Time> TimeArray;
    typedef smartPtr<TimeArray> TimeArraySP;


private:
    struct GTOper;
#ifdef DEBUG
public: // ideally not but want to call from function in DateTime.cpp
    //// handy for debugging. DO NOT USE ANYWHERE ELSE
    const char* p() const;
private:
#endif
    /// fields
    int date;
    int time;
};

typedef smartPtr<DateTime> DateTimeSP;
typedef smartConstPtr<DateTime> DateTimeConstSP;
#ifndef QLIB_DATETIME_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<DateTime>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DateTime>);
#endif

#if !defined(DEBUG) || defined(QLIB_DATETIME_CPP)
//// inlined functions for performance
OPT_INLINE DateTime::DateTime(int date, int time): IObject(), CObject(TYPE), date(date), time(time) {}
OPT_INLINE DateTime::DateTime(const DateTime& dt): IObject(), CObject(TYPE), date(dt.date), time(dt.time){}
OPT_INLINE DateTime::~DateTime(){}
OPT_INLINE int DateTime::getDate() const{return date;}
OPT_INLINE DateTime::Date::Date(const Date& dt): CObject(TYPE), date(dt.date){}
OPT_INLINE DateTime& DateTime::operator=(const DateTime& rhs){date=rhs.date;time=rhs.time; return *this;}
#endif

/** specialisations of arrayObjectCast */
template <> class TOOLKIT_DLL arrayObjectCast<DateTime>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const DateTime& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(DateTime& value);

    /** Turns the IObjectSP into a DateTime */
    static const DateTime& fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class TOOLKIT_DLL arrayClone<DateTime>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

/** DateTimeArray is an array of DateTime structures - not an array of
    pointers. Note that this approach forces the no argument DateTime
    constructor to be public */

//typedef smartConstPtr<DateTimeArray> DateTimeArrayConstSP;

#ifndef QLIB_DATETIME_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<DateTimeArray _COMMA_ DateTimeArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<DateTime::DateSP _COMMA_ DateTime::Date>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<DateTime::TimeSP _COMMA_ DateTime::Time>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<DateTimeArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DateTimeArray>);
// again, avoid code bloat by declaring common implementation of function
// templates as extern
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<DateTime>(DateTime* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<DateTimeArray>(DateTimeArray* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<DateTime>(DateTime* t, IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<DateTimeArray>(DateTimeArray* t,
                                                               IObjectSP o));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<DateTimeSP>(DateTimeSP* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<DateTimeArraySP>(
                    DateTimeArraySP* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<DateTimeSP>(DateTimeSP* t,
                                                              IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<DateTimeArraySP>(DateTimeArraySP* t,
                                                                   IObjectSP o));
#else
//// seem to need this before the rest of this file otherwise (MSVC) the 
//// template is instantiated
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<DateTimeArray _COMMA_ DateTimeArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<DateTime::DateSP _COMMA_ DateTime::Date>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<DateTime::TimeSP _COMMA_ DateTime::Time>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<DateTime>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DateTime>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<DateTimeArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DateTimeArray>);
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<DateTime>(DateTime* t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<DateTimeArray>(DateTimeArray* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<DateTime>(DateTime* t, IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<DateTimeArray>(DateTimeArray* t,
                                                                    IObjectSP o));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<DateTimeSP>(DateTimeSP* t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<DateTimeArraySP>(
                         DateTimeArraySP* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<DateTimeSP>(DateTimeSP* t,
                                                                   IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<DateTimeArraySP>(DateTimeArraySP* t,
                                                                        IObjectSP o));
#endif

typedef DateTime CDateTime;
typedef DateTimeArray CDateTimeArray;
typedef DateTimeArraySP CDateTimeArraySP;
typedef DateTimeArrayConstSP  CDateTimeArrayConstSP;

/** Defining a DateTimeCluster to be a DateTimeArrayArray. Note an array of
    structure not an array of pointers */
typedef array<DateTimeArray> DateTimeArrayArray;
// DateTimeArrayArray is the 'correct' name, cluster is just a convenience
typedef DateTimeArrayArray  DateTimeCluster;

/** specialisations of arrayObjectCast (needed as the array is not an array
    of pointers) */
template <> class TOOLKIT_DLL arrayObjectCast<DateTimeArray>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const DateTimeArray& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(DateTimeArray& value);

    /** Turns the IObjectSP into a DateTime */
    static const DateTimeArray& fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class TOOLKIT_DLL arrayClone<DateTimeArray>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

DECLARE(DateTimeCluster);

/** Class for Addin functions that require a single optional datetime */
class TOOLKIT_DLL DateTimeAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /** the single [optional] parameter that the addin takes */
    DateTimeSP  dateTime;

private:
    /** for reflection */
    DateTimeAddin();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    static IObject* defaultDateTimeAddin();
    
};

/** Class for Addin functions that require an optional DateTimeArray */
class TOOLKIT_DLL DateTimeArrayAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /** the single [optional] parameter that the addin takes */
    DateTimeArraySP  dateTimeArray;

private:
    /** for reflection */
    DateTimeArrayAddin();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    static IObject* defaultDateTimeArrayAddin();
};

/** wrapper around DateTimeArray that supports smart adding & scaling */
class TOOLKIT_DLL DateTimeList: public CObject,
                                virtual public CombinableResult {
public:
    static CClassConstSP const TYPE;

    /** stores unique, sorted version of dl */
    DateTimeList(const DateTimeArray* dl);

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add DateTimeList to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

private:
    friend class DateTimeListHelper;
    DateTimeList();

    DateTimeArraySP dl;
};

DECLARE(DateTimeList);

#ifndef QLIB_DATETIME_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<DateTimeCluster>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DateTimeCluster>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DateTimeList>);
#endif

DRLIB_END_NAMESPACE

#endif


