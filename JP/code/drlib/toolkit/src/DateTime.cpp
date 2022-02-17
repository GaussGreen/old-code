//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DateTime.cpp
//
//   Description : Prototype date-time representation
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_DATETIME_CPP
#include "edginc/DateTime.hpp"
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Null.hpp"
#include "edginc/Writer.hpp"
#include "edginc/BoxedEnum.hpp"

#include <float.h>   /* DBL_EPSILON cannot use Maths:: methods because
                        util is before toolkit in the include path */
#include <algorithm>
#include <queue>

DRLIB_BEGIN_NAMESPACE

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<DateTimeCluster>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DateTimeCluster>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DateTimeList>);


typedef vector<DateTime::Interval> DateTimeIntervalArray;

const int DateTime::TIME_IN_DAY = 86400;
const int DateTime::TIME_MAX    = TIME_IN_DAY - 1;
const int DateTime::TIME_MIN    = 0;
const int DateTime::START_OF_DAY_TIME = 5;
const int DateTime::BEFORE_EX_DIV_TIME = 3;
const int DateTime::EX_DIV_BEFORE_START_TIME = 4;
const int DateTime::END_OF_DAY_TIME = TIME_MAX;
const int DateTime::EFFECTIVE_TIME_IN_DAY = TIME_MAX - START_OF_DAY_TIME;
const int DateTime::DAYS_PER_YEAR = 365;
const int DateTime::TIME_IN_YEAR  = DAYS_PER_YEAR * EFFECTIVE_TIME_IN_DAY;
const int DateTime::DAYS_PER_WEEK = 7;
const int DateTime::VOL_TIME_IN_DAY = END_OF_DAY_TIME - START_OF_DAY_TIME;
const int DateTime::VOL_TIME_IN_YEAR = DAYS_PER_YEAR * VOL_TIME_IN_DAY;

// IMPORTANT:  Be sure to keep enum DateTime::DayOfWeek (in header) synchronized with these values!
const int DateTime::SATURDAY      = 0;
const int DateTime::SUNDAY        = 1;
const int DateTime::MONDAY        = 2;
const int DateTime::TUESDAY       = 3;
const int DateTime::WEDNESDAY     = 4;
const int DateTime::THURSDAY      = 5;
const int DateTime::FRIDAY        = 6;

START_PUBLIC_ENUM_DEFINITION(DateTime::DayOfWeek, "Names the days of the week");
ENUM_VALUE_AND_NAME(DateTime::Saturday, "Saturday", "Saturday");
ENUM_VALUE_AND_NAME(DateTime::Sunday, "Sunday", "Sunday");
ENUM_VALUE_AND_NAME(DateTime::Monday, "Monday", "Monday");
ENUM_VALUE_AND_NAME(DateTime::Tuesday, "Tuesday", "Tuesday");
ENUM_VALUE_AND_NAME(DateTime::Wednesday, "Wednesday", "Wednesday");
ENUM_VALUE_AND_NAME(DateTime::Thursday, "Thursday", "Thursday");
ENUM_VALUE_AND_NAME(DateTime::Friday, "Friday", "Friday");
END_ENUM_DEFINITION(DateTime::DayOfWeek);

const string DateTime::START_OF_DAY        = "SOD";
const string DateTime::END_OF_DAY          = "EOD";
const string DateTime::BEFORE_EX_DIV       = "BEX";
const string DateTime::EX_DIV_BEFORE_START = "XBS";

// convert dd-mmm-yyyy into a date
static int convertDate(string date);

/** takes date in dd-mmm-yyyy format, time is either
    START_OF_DAY, END_OF_DAY or hh-mm-ss format */
DateTime::DateTime(const string& date, const string& time): CObject(TYPE){
    try{
        this->date = convertDate(date);
        this->time = DateTime::timeConvert(time);
    } catch (exception){ // fix for gcc
        throw;
    }
}

bool DateTime::isLeapYear(long year) {
    return ( ((year)%4   == 0) && ((year)%100 != 0)) ||
        ((year)%400 == 0);
}

/** Calculates the expected end date of a list of future dates which
    are weighted */
const DateTime DateTime::expectedFutureDate(
    const DateTime& today,
    const array<DateTime>& futDates,
    const DoubleArray& futWeights)
{
    static const string method = "DateTime::expectedFutureDate";
    double weighting  = 0.0;
    int i=0;

    // check that the dates and weights arrays are not empty
    if (futDates.size() == 0 ||
        futWeights.size() == 0)
    {
        throw ModelException(method,
                             "future dates and/or future weights arrays "
                             "cannot be empty");
    }

    // check the date and weights arrays are the same length
    if (futWeights.size() != futDates.size())
    {
        throw ModelException(method,
                             "future dates and future weights arrays must be "
                             "the same length");
    }

    int numDates = futDates.size();

    // check that future dates are strictly increasing and are after today
    for (i=1; i < numDates; i++)
    {
        if (futDates[i-1].isGreaterOrEqual(futDates[i]))
        {
            throw ModelException(method,
                                 "future dates not in ascending order.\n" +
                                 futDates[i-1].toString() +
                                 " occurs before " +
                                 futDates[i].toString() +
                                 " in the dates list.");
        }
    }

    // check that future dates are greater than today
    if (today.isGreaterOrEqual(futDates[0]))
    {
        throw ModelException(method,
                             "today " +
                             today.toString() +
                             " is after first future date " +
                             futDates[0].toString());
    }

    // check that weights are all positive
    for (i = 0; i < numDates; i++)
    {
        if (futWeights[i] < - DBL_EPSILON)
        {
            throw ModelException(method,
                                 "future weight " +
                                 Format::toString(futWeights[i]) +
                                 " is negative");
        }
    }

    DateTimeIntervalArray intervals;
    for (i = 0; i < numDates; i++)
    {
        intervals.push_back(futDates[i].subtract(today));
        weighting += futWeights[i];
    }

    // if weighting = zero
    if (weighting < DBL_EPSILON && weighting > -DBL_EPSILON)
    {
        throw ModelException(method,
                             "sum of future futWeights must be > 0\n");
    }

    double days = 0.0;
    double time = 0.0;

    for (i = 0; i < numDates; i++)
    {
        days += intervals[i].dateDiff * futWeights[i];
        time += intervals[i].timeDiff * futWeights[i];
    }

    days /= weighting;
    time /= weighting;

    // add any fractional days onto time
    time += DateTime::TIME_IN_DAY * modf(days, &days);

    // split time into any whole days & leftover time
    int extraDays = (int)(time / DateTime::TIME_IN_DAY);
    int extraTime = (int)(time - extraDays * DateTime::TIME_IN_DAY);

    DateTime::Interval expectedLife(((int)(days + extraDays)), extraTime);

    // add expected life onto today to get expected end
    DateTime expectedEndDate = today.add(expectedLife);

    return expectedEndDate;
}

// DateTime::DateTime(int date, int time) constructor in header file for
// performance reasons

/** Returns true if the datetime has not been initialsed ie it has been
    built through the default constructor */
bool DateTime::empty() const{
    return (date == 0 && time == 0);
}


// for reflection
DateTime::DateTime(): CObject(TYPE), date(0), time(0){
    // empty
}

// ~DateTime() is in the .hpp for performance

/** Overridden for performance */
IObject* DateTime::clone() const{
    return new DateTime(*this);
}

bool DateTime::equals(const DateTime &date, bool compTime) const {
    return (this->date == date.date && (!compTime || this->time == date.time));
}

int DateTime::compareTo(const DateTime &date) const {
    if (isGreater(date)) {
        return 1;
    }
    else if (equals(date)) {
        return 0;
    }
    else {
        return -1;
    }
}

bool DateTime::operator < (const DateTime& date) const
{
    return isLess(date);
}
bool DateTime::operator <= (const DateTime& date) const
{
    return !isGreater(date);
}
bool DateTime::operator > (const DateTime& date) const
{
    return isGreater(date);
}
bool DateTime::operator >= (const DateTime& date) const
{
    return !isLess(date);
}
bool DateTime::operator == (const DateTime& date) const
{
    return (this->date == date.date && this->time == date.time);
}
bool DateTime::operator != (const DateTime& date) const
{
    return (this->date != date.date || this->time != date.time);
}


bool DateTime::isGreater(const DateTime &date) const {
    return (this->date > date.date ||
            ((this->date == date.date) && (this->time > date.time)));
}

bool DateTime::isGreaterOrEqual(const DateTime &date) const {
    return (this->date > date.date ||
            ((this->date == date.date) && (this->time >= date.time)));
}

bool DateTime::isLess(const DateTime &date) const {
    return (this->date < date.date ||
            ((this->date == date.date) && (this->time < date.time)));
}

bool DateTime::isLessOrEqual(const DateTime &date) const {
    return (this->date < date.date ||
            ((this->date == date.date) && (this->time <= date.time)));
}

// for weak (ignoring time) comparison
bool DateTime::isLessWeak(const DateTime& date1, const DateTime& date2) {
    return (date1.date < date2.date);
}

/** Returns MAX(this, date) */
const DateTime& DateTime::max(const DateTime& date) const{
    return isGreater(date)? *this: date;
}

/** Returns MAX(this, date1, date2) ie maximum of all three dates */
const DateTime& DateTime::max(const DateTime& date1,
                              const DateTime& date2) const{
    if (*this > date1) {
        return max(date2);
    }
    return date1.max(date2);
}

/** Returns MIN(this, date) */
const DateTime& DateTime::min(const DateTime& date) const{
    return isLess(date)? *this: date;
}

DateTime DateTime::rollDate(const int offset) const {
    return DateTime(date + offset, time);
}

/* Returns a DateTime with the number of months increased by offset */
DateTime DateTime::rollDateInMonths(const int monthsOffset) const {
    MonthDayYear mdy = toMDY();
    mdy.month += monthsOffset;

    return mdy.toDateTime(time);
}

// Attempts to answer the nth specified weekday of the month of this DateTime.
// For example, to get the 3rd Thursday of the month, set n to 3, and
// dayOfWeek to DateTime::THURSDAY (or DateTime::Thursday).
DateTime DateTime::nthWeekdayOfMonth(int n, int dayOfWeek) const {
    try {
        MonthDayYear mdy(toMDY());
        int firstWeekdayOfMonth = firstOfMonth().getWeekday();

        mdy.day = 1 + (n - 1) * 7 + dayOfWeek - firstWeekdayOfMonth;
        if (dayOfWeek < firstWeekdayOfMonth)
            mdy.day += 7;

        // Attempt to convert to DateTime (but this may fail if we've overshot the days in the month)
        return mdy.toDateTime(time);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "Failed");
    } 
}

// Answers the last specified weekday of the month of this DateTime.
DateTime DateTime::lastWeekdayOfMonth(int dayOfWeek) const {
    DateTime endOfMonth(returnEndOfMonth(false));

    int daysFromLastWeekday = endOfMonth.getWeekday() - dayOfWeek;
    if (daysFromLastWeekday < 0)
        daysFromLastWeekday += 7;

    return endOfMonth.rollDate(-daysFromLastWeekday);
}

/** Validates that supplied dates are increasing order. If failIfEmpty is
    true, then throws exception if dates is empty. The description string
    is used in the exception messages should describe what the dates
    represent */
void DateTime::ensureIncreasing(const DateTimeArray& dates,
                                const char*          description,
                                bool                 failIfEmpty){
    static const string routine("DateTime::ensureIncreasing");
    if (failIfEmpty && dates.empty()){
        throw ModelException(routine, string(description) +
                             " (array of dates) is empty");
    }
    for (int i = 1; i < dates.size(); i++){
        if (dates[i-1].isGreater(dates[i])){
            throw ModelException(routine, "Dates in " + string(description) +
                                 " (array of dates) are not increasing : [" +
                                 Format::toString(i) + "] = " + dates[i-1].toString() +
                                 " is after [" + Format::toString(i+1) + "] = " +
                                 dates[i].toString() + "]");
        }
    }
}

/** Validates that supplied dates are increasing order. If failIfEmpty is
    true, then throws exception if dates is empty. The description string
    is used in the exception messages should describe what the dates
    represent */
void DateTime::ensureStrictlyIncreasing(const DateTimeArray& dates,
                                        const char*          description,
                                        bool                 failIfEmpty){
    static const string routine("DateTime::ensureStrictlyIncreasing");
    if (failIfEmpty && dates.empty()){
        throw ModelException(routine, string(description) +
                             " (array of dates) is empty");
    }
    for (int i = 1; i < dates.size(); i++){
        if (!(dates[i].isGreater(dates[i-1]))) {
            throw ModelException(routine, "Dates in " + string(description) +
                                 " (array of dates) are not striclt increasing : [" +
                                 Format::toString(i) + "] = " + dates[i-1].toString() +
                                 " is after or equal to [" + Format::toString(i+1) + "] = " +
                                 dates[i].toString() + "]");
        }
    }
}
//// same as above but takes a string
void DateTime::ensureStrictlyIncreasing(const DateTimeArray& dates,
                                        const string&        description,
                                        bool                 failIfEmpty){
    ensureStrictlyIncreasing(dates, description.c_str(), failIfEmpty);
}

//// same as above but takes a string
void DateTime::ensureIncreasing(const DateTimeArray& dates,
                                const string&        description,
                                bool                 failIfEmpty){
    ensureIncreasing(dates, description.c_str(), failIfEmpty);
}

/** Returns the number of dates strictly greater than this date. Assumes
    dates are in increasing order */
int DateTime::numFutureDates(const DateTimeArray& allDates) const{
    return distance(upper_bound(allDates.begin(), allDates.end(), *this), allDates.end());
//    return (allDates.size()-(findLower(allDates)+1));

}

/** Returns an array of dates strictly greater than baseDate. Assumes
    dates are in strictly increasing order */
DateTimeArray DateTime::getFutureDates(const DateTimeArray& allDates) const
{
//     int i;
//     for (i = 0; i < allDates.size(); ++i)
//         if (allDates[i]>*this)
//             break;
//
//     return DateTimeArray(allDates.begin()+i, allDates.end());
    return DateTimeArray(upper_bound(allDates.begin(), allDates.end(), *this), allDates.end());
}

DateTimeArray DateTime::getNotPastDates(const DateTimeArray& allDates) const{
    return DateTimeArray(lower_bound(allDates.begin(), allDates.end(), *this), allDates.end());
//     DateTimeArray futureDates;
//     futureDates.reserve(allDates.size()/2);
//     int i;
//     for (i = 0; i < allDates.size(); ++i)
//         if (allDates[i]>=*this)
//             break;
//
//     for (; i < allDates.size(); ++i)
//         futureDates.push_back(allDates[i]);
//
//     return futureDates;
}

/** Returns the number of dates <= this date. Assumes
    dates are in increasing order */
int DateTime::numPastDates(const DateTimeArray& allDates) const{
    return (findLower(allDates)+1);
#if 0
    // to do: replace with a binary search
    int numDates = allDates.size();
    int i;
    for (i = 0; i < numDates && isGreaterOrEqual(allDates[i]); i++); // empty
    return i;
#endif
}

/** Returns an array of dates <= baseDate. Assumes dates are in order */
DateTimeArray DateTime::getPastDates(const DateTimeArray& allDates) const{
    int i = numPastDates(allDates);
    DateTimeArray pastDates(i);
    for (int j = 0; j < i; j++){
        pastDates[j] = allDates[j];
    }
    return pastDates;
}

/** Returns an array of dates from allDates that are bounded
    inclusively by start and end dates */
DateTimeArray DateTime::getInclusiveDates(const DateTime & start,
                                          const DateTime & end,
                                          const DateTimeArray & allDates)
{
    ASSERT(end.isGreater(start));
    DateTimeArray result;

    if (end.isGreaterOrEqual(allDates.front()) || (allDates.back()).isGreaterOrEqual(start))
        // only do it if [start, end] overlaps with allDates
    {
        int startIdx = start.findUpper(allDates);
        int endIdx = end.findLower(allDates);

        if (startIdx <= endIdx) {
            // only do it if [start, end] actually bracket a date in allDates
            result.resize(endIdx - startIdx + 1);
            for (int i = startIdx, j = 0; i <= endIdx; ++i, ++j)
                result[j] = allDates[i];
        }
    }
    return result;
}


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
IntArray  DateTime::createMapping(const DateTimeArray&  dates,
                                  const DateTimeArray&  datesForMap,
                                  bool&                 isTrivial){
    IntArray myMap(dates.size()+1);
    // start by working backwards through the 2 arrays
    int iMapDate = datesForMap.size()-1;
    int offset = 0;
    isTrivial = true; // default
    for (int iDate = dates.size()-1; iDate >= 0; iDate--){
        if (iMapDate >= 0 && datesForMap[iMapDate].equals(dates[iDate])){
            offset = 0;
            iMapDate--;
        } else{
            offset++;
            isTrivial = false;
        }
        myMap[iDate] = offset;
    }
    if (iMapDate >= 0){
        throw ModelException("DateTime::createMap",
                             "Supplied date ("+
                             datesForMap[iMapDate].toString()+")"
                             " not found in list of all dates");
    }
    return myMap;
}

vector<int> DateTime::createMapping2(
    const DateTimeArray& dates,          // assumed sorted
    const DateTimeArray& datesForOffset, // assumed sorted
    const DateTime& valueDate,
    int& lastPastDateIdx )
{
    // FIX to take into account dates are sorted.
    const DateTimeArray& D = dates;
    vector<int> keyDate( datesForOffset.size() );
    array<DateTime>::const_iterator I;
    for ( size_t i = 0; i < keyDate.size(); ++i ) {
        const DateTime& V = datesForOffset[i];
        I = std::find( D.begin(), D.end(), V );
        if ( I == D.end() ) {
            throw ModelException("Could not find date '"
                                 + V.toString() + "' in the input date list");
        }
        keyDate[i] = I - D.begin();
    }
    I = std::lower_bound( datesForOffset.begin(), datesForOffset.end(), valueDate );
    lastPastDateIdx = I - datesForOffset.begin();
    return keyDate;

}
/** Creates an array of bools of length dates.size(). Each
    element in the array is set to true if the date exists in the
    datesForMap. Eg
    Dates:              A B C D E F G H I J
    datesForMap:          B     E       I
    map:                F T F F T F F F T F
    Typical use of such a map is if we want to do averaging over certain days etc
    ...
    } Arrays need not be sorted */

vector<bool> DateTime::createmapping3(
    const DateTimeArray& dates,
    const DateTimeArray& datesForMapping )
{
    const DateTimeArray& D = datesForMapping;
    vector<bool> keyDate(dates.size() );
    array<DateTime>::const_iterator I;
    for (size_t i = 0; i < keyDate.size(); ++i) {
        const DateTime& V = dates[i];
        I = std::find( D.begin(), D.end(), V);
        if ( I == D.end() ) {
            keyDate[i] = false;
        } else {
            keyDate[i] = true;
        }
    }
    return keyDate;
}

/** Returns true if the two lists of dates are identical */
bool DateTime::equals(const DateTimeArray& dates1,
                      const DateTimeArray& dates2){
    if (&dates1 == &dates2){
        return true; // could be from the same smart pointer
    }
    int len = dates1.size();
    if (len != dates2.size()){
        return false;
    }
    for (int i = 0; i < len; i++){
        if (!dates1[i].equals(dates2[i])){
            return false;
        }
    }
    return true;
}

/** merge n date arrays into one removing duplicate dates. All date arrays
    must be be sorted in increasing order and contain no duplicates
    themselves */
DateTimeArray DateTime::merge(const vector<DateTimeArray>& dtArrayVector){
    int numArrays = dtArrayVector.size();
    vector<const DateTimeArray*>   arrayByPtr(numArrays);
    for (int i = 0; i < numArrays; i++){
        arrayByPtr[i] = &dtArrayVector[i];
    }
    return merge(arrayByPtr);
}


/** merge n date arrays into one removing duplicate dates. All date arrays
    must be be sorted in increasing order and contain no duplicates
    themselves */
DateTimeArray DateTime::merge(
    const vector<const DateTimeArray*>& dtArrayVector){
    DateTimeArray mergedDates;
    merge(dtArrayVector, mergedDates);
    return mergedDates;
}

/** Same as above, put reference to DateTimeArray supplied to write results
    to */
void DateTime::merge(const vector<const DateTimeArray*>& dtArrayVector, // (I)
                     DateTimeArray&                      mergedDates){  // (O)
    vector<DateTimeArray::const_iterator> dtStart;
    vector<DateTimeArray::const_iterator> dtEnd;

    dtStart.reserve(dtArrayVector.size());
    dtEnd.reserve(dtArrayVector.size());

    for (size_t i = 0; i < dtArrayVector.size(); i++){
        dtStart.push_back(dtArrayVector[i]->begin());
        dtEnd.push_back(dtArrayVector[i]->end());
    }
    merge(dtStart, dtEnd, mergedDates);
}

/** Another merge except those dates <= lowerLimit or > upperLimit are
    omitted */
void DateTime::merge(const DateTime&                     lowerLimit,
                     const DateTime&                     upperLimit,
                     const vector<const DateTimeArray*>& dtArrayVector, // (I)
                     DateTimeArray&                      mergedDates){  // (O)
    //ensure lowerLimit is indeed earlier than upperLimit
    if (lowerLimit.isGreater(upperLimit))
    {
        throw ModelException("DateTime::merge","lowerLimit is greater than upperLimit");
    }
    vector<DateTimeArray::const_iterator> dtStart;
	dtStart.reserve(dtArrayVector.size());
    vector<DateTimeArray::const_iterator> dtEnd;
	dtEnd.reserve(dtArrayVector.size());
    for (size_t i = 0; i != dtArrayVector.size(); i++){
        dtStart.push_back( dtArrayVector[i]->begin()+
						   lowerLimit.numPastDates(*dtArrayVector[i]));
        dtEnd.push_back( dtArrayVector[i]->end()-
						 upperLimit.numFutureDates(*dtArrayVector[i]));
    }
    merge(dtStart, dtEnd, mergedDates);
}


/** Aux class that says date1 < date2 iff date2 is actually before the date1. Used in priority queue below.
 */
class DateComparator : public binary_function<size_t, size_t, bool>
{
    const vector<DateTimeArray::const_iterator> & curr_iterators;
public:
    DateComparator(const vector<DateTimeArray::const_iterator> & iters) : curr_iterators(iters)
        {
        }

    bool operator()(const size_t & l, const size_t & r)
        {
            return *curr_iterators[r] < *curr_iterators[l]; // we have to have "top" to be the least in the order
        }
};

/*
   The previously used method had two performance drawbacks:
   1. it searches for the smallest active element in a loop => very bad for a big number of arrays
   2. it inserts elements one by one => prevents fast copies of vector's ranges

   We overcome the first issue by using a heap (priority_queue over a
   vector) -- just note that priority_queue returns the largest element,
   not the smallest, so we need an aux class to reverse comparison.

   The second issue is solved by postponing push_back-s until linear
   chain of inserts is interrupted. This should work fine for the common
   situation when one element is added to a big array.

   Overall, we spent O(log(dtStart.size())) time on each element.
*/

/** merge n date arrays into one removing duplicate dates. All date arrays
    must be be sorted in increasing order and contain no duplicates
    themselves */
void DateTime::merge(
    const vector<DateTimeArray::const_iterator>& dtStart,     // (I)
    const vector<DateTimeArray::const_iterator>& dtEnd,       // (I)
    DateTimeArray&                               mergedDates) // (O)
{
    if (dtStart.size() != dtEnd.size())
        throw ModelException("dtStart.size() != dtEnd.size()");

    const size_t n = dtStart.size();
    vector<DateTimeArray::const_iterator> pos(dtStart); // current positions in dtArrayVector[]

    // the "largest" element should be the smallest in our ordering
    DateComparator dateComparator(pos); // compares dates pointed to by pos[i] and pos[j]

    /* keep an active element from each of the dtArrayVector[] */
    priority_queue<size_t, vector<size_t>, DateComparator> heap(dateComparator);

    // populate heap
    for(size_t i=0; i < n; ++i) {
        // initialize heap: remember that dtArrayVector[i] is assumed to be sorted already
        if (pos[i] != dtEnd[i])
            heap.push(i); // insert positions into the heap (queue) from all non-empty ranges
    }

    // this function is frequently called with the second array being
    // a 1-elem vector; do some more optimizations to insert
    // consequitive ranges at a time.

    size_t last_top = n+1; // (n+1) is the "out of range" marker
    DateTimeArray::const_iterator last_start; // we also keep the start of the range to be copied

    // Main loop; at each iteration at least one elem is popped out and at most one is pushed in
    while(!heap.empty()) {
        size_t top = heap.top(); // found the current least element
        heap.pop(); // pop the top element out of the queue

        if (top != last_top) { // active element is from array different than the prev. one
            if (last_top != n+1) {
                // copy consecutive range dtArrayVector[last_top] onto result
                mergedDates.insert(mergedDates.end(), last_start, pos[last_top]);
                last_top = n+1;  // mark that we copied the range
            }
            // check if what we found should be added to the answer
            if (mergedDates.empty() || (mergedDates.back() != *pos[top])) {
                last_top = top; // notice the source of the new range
                last_start = pos[last_top]; // notice location in the new range
            }
        }
        ++pos[top]; //advance to the next element
        // if there is next element, submit it to the heap
        if (pos[top] != dtEnd[top])
            heap.push(top);
    } // end of while(!heap.empty())

    if (last_top != n+1) // insert the last range
        mergedDates.insert(mergedDates.end(), last_start, pos[last_top]);

    return;
}

/** merge 2 date arrays into one removing duplicate dates. Both date arrays
    must be be sorted in increasing order and contain no duplicates
    themselves */
DateTimeArray DateTime::merge(const DateTimeArray& dtArray1,
                              const DateTimeArray& dtArray2){
    if (dtArray1.empty()){
        return dtArray2;
    }
    if (dtArray2.empty()){
        return dtArray1;
    }
    vector<const DateTimeArray*>   arrayByPtr(2);
    arrayByPtr[0] = &dtArray1;
    arrayByPtr[1] = &dtArray2;
    return merge(arrayByPtr);
}

/** Simple extension to merge via SPs */
DateTimeArray DateTime::merge(const vector<DateTimeArrayConstSP>& v)
{
    vector<const DateTimeArray *>  l;
    for(size_t i=0; i != v.size(); ++i)
        l.push_back(v[i].get());
    return DateTime::merge(l);
}

/** Removes duplicate dates in a DateTimeArray provided that the duplicate
    dates are next to each other */
/* Attn: possibly calls erase() in the middle of a vector in a loop: O(n^2) worst case complexity */
void DateTime::removeDuplicates(DateTimeArray& dates,
                                bool           ignoreTimeOfDay){
    for (vector<DateTime>::iterator iter(dates.begin()); iter != dates.end();
         /* inc in loop */){
        if (iter != dates.begin() &&
            iter->equals(*(iter-1), !ignoreTimeOfDay)){
            iter = dates.erase(iter);
        } else {
            ++iter;
        }
    }
}

/** Removes dates which lie outside lowerDate or upperDate */
/* Attn: possibly calls erase() in the middle of a vector in a loop: O(n^2) worst case complexity */
void DateTime::removeOutliers(DateTimeArray&  dates,
                              const DateTime& lowerDate,
                              const DateTime& upperDate){
    for (vector<DateTime>::iterator iter(dates.begin()); iter != dates.end();
         /* inc in loop */){
        if (lowerDate.isGreater(*iter) || upperDate.isLess(*iter)){
            iter = dates.erase(iter);
        } else {
            ++iter;
        }
    }
}


/** What do I do? Why aren't the parameters const? RA to document */
DateTimeArraySP DateTime::subtract(DateTimeArray& first, DateTimeArray& second)
{
    // ensure date lists are sorted
    sort(first.begin(), first.end());
    sort(second.begin(), second.end());

    DateTimeArraySP result(new DateTimeArray());
    set_difference(first.begin(), first.end(),
                   second.begin(), second.end(),
                   result->back_inserter());
    return result;
}


/** returns an array of length part.size() giving the index of each element of
    'part' in 'full'. Both arrays must be ordered [strictly increasing], and
    'part' must be a subset of 'full' */
vector<int> DateTime::getIndexes(const DateTimeArray& full,
                                 const DateTimeArray& part){
    vector<int> indexes(part.size());
    int j = 0;
    for (int i = 0; i < part.size(); i++){
        while (!part[i].equals(full[j]) && j < full.size()){
            j++;
        }
        if (j == full.size()){
            throw ModelException("DateTime::getIndexes", "invalid inputs");
        }
        indexes[i] = j;
        j++;
    }
    return indexes;
}

/* returns indexes of projection
   t0 t1 t2 t3 t4
   t2    t4
   yields:
   -1 -1  0 -1  1
*/
vector<int> DateTime::getProjection(const DateTimeArray& full,
                                    const DateTimeArray& part)
{
    vector<int> locations = DateTime::getIndexes(full, part);
    vector<int> result(full.size(), -1);
    for(int i=0; i < part.size(); ++i)
        result[locations[i]] = i;
    return result;
}

/** returns an array of size full.size(), giving an index j for full[i]
    such that: j = min(k: part[k] >= full[i]) if full[i] <= part[last],
    otherwise, j = last.  In other words, the index j for full[i] is the
    index of the smallest element in part that is bigger than or equal to
    full[i], given full[i] <= part[last].  If full[i] >     part[last], then
    j = last.

    Example 1:
    t0 t1 t2 t3 t4 t5               full
    t1    t3       t6    part
    yields:
    0  0  1  1  2  2

    Example 2:
    t0 t1 t2 t3 t4 t5 t6    full
    t1    t3     t4                      part
    yields:
    0  0  1  1  2  2  2

*/
vector<int> DateTime::getCeilingProjection(const DateTimeArray & full,
                                           const DateTimeArray & part)
{
    ASSERT(part.size() > 0);
    size_t fullSize = full.size();
    size_t partSize = part.size();
    vector<int> result(fullSize);

    size_t partIdx = 0;
    for (size_t i = 0; i < fullSize; ++i)
    {
        while (full[i].isGreater(part[partIdx])) {
            if (partIdx < partSize - 1)
                ++partIdx;
            else
                break;
        }
        result[i] = partIdx;
    }
    return result;
}

/** getFloorProjection is similar to getCeilingProjection except it returns
    indices j's for full[i] such that: j = max(k: part[k] <= full[i]), if
    full[i] >= part[first], otherwise j = first.

    Example 1:
    t0 t1 t2 t3 t4 t5       full
    t0    t2 t3             part
    yields:
    0  0  1  2  2  2

    Example 2:
    t0 t1 t2 t3 t4 t5 t6    full
    t2    t4    t6    part
    yields:
    0  0  0  0  1  1  2
*/
vector<int> DateTime::getFloorProjection(const DateTimeArray & full,
                                         const DateTimeArray & part)
{
    ASSERT(part.size() > 0);
    size_t fullSize = full.size();
    size_t partSize = part.size();
    vector<int> result(fullSize);

    size_t partIdx = partSize - 1;
    for (int i = fullSize - 1; i >= 0; --i) // i has to be int!
    {
        while (part[partIdx].isGreater(full[i])) {
            if (partIdx > 0)
                --partIdx;
            else
                break;
        }
        result[i] = partIdx;
    }
    return result;
}


/** Rather nasty function which forces the time of day in each of the dates
    supplied to the specified value */
void DateTime::setTimeOfDay(DateTimeArray& dates, // (M)
                            int            timeOfDay){
    for (int i = 0; i < dates.size(); i++){
        dates[i].time = timeOfDay;
    }
}

//// silly class required by STL
struct DateTime::GTOper{
    bool operator()(const DateTime& d1, const DateTime& d2){
        return d1.isGreater(d2);
    }
};

/** Dates must be increasing (but does not have to be strictly).
    Returns the index of the [last] date which is greater
    or equal than the supplied target date. Returns dates.size()
    if target date > all dates. */
int DateTime::findUpper(const DateTimeArray& dates) const{
    // The impenetrable STL again ....
    // Hopefully it is obvious to you that this is the correct algorithm
    // can't say it was obvious to me :-)
    return dates.rend() -
        upper_bound(dates.rbegin(), dates.rend(), *this, GTOper());
}

/** Dates must be increasing (but does not have to be strictly).
    Returns the index of the [last] date which is less or equal
    than the supplied target date. Returns -1 if target date < all
    dates.  */
int DateTime::findLower(const DateTimeArray& dates) const{
    // The impenetrable STL again ....
    // Hopefully it is obvious to you that this is the correct algorithm
    // can't say it was obvious to me :-)
    return upper_bound(dates.begin(), dates.end(), *this) - dates.begin()-1;
}

/** Dates must be strictly increasing.
	Returns the index of the array date that is nearest to the date.
	Returns -1 if the array is empty
*/
int DateTime::findNearest(const DateTimeArray& dates) const
{
	if (dates.size() == 0) return -1;
	if (*this <= dates[0]) return 0;
	int size = dates.size();
	if (*this >= dates.back()) return (size-1);

	for (int i=1; i< size; ++i)
	{
		if ((dates[i-1] < *this) && (*this <= dates[i]))
		{
			Interval left = this->subtract(dates[i-1]);
			Interval right = dates[i].subtract(*this);
			if (right.isGreaterOrEqual(left))
				return i-1;
			else
				return i;
		}
	}
	return -1;// it will never reach here
}

/** Return iterator to the date.  The search is done on the interval
    [dateStart, dateEnd), ie. excluding dateEnd as in stl */
DateTimeArray::const_iterator DateTime::findIterator(
    DateTimeArray::const_iterator dateStart,
    DateTimeArray::const_iterator dateEnd) const
{
    pair<DateTimeArray::const_iterator, DateTimeArray::const_iterator> loc =
        equal_range(dateStart, dateEnd, *this);
    if (loc.first == loc.second){
        return dateEnd;
    }
    return loc.first;
}

/** Returns the index of the first instance of this date in the array
    of dates. Fails if the date does not exist. The supplied dates
    must be weakly increasing (a binary search is used) */
int DateTime::find(const DateTimeArray& dates) const{
    return find(dates.begin(), dates.end());
}

/** Same as above find if used dates.begin(), dates.end() */
int DateTime::find(DateTimeArray::const_iterator dateStart,
                   DateTimeArray::const_iterator dateEnd) const{
    pair<DateTimeArray::const_iterator, DateTimeArray::const_iterator> loc =
        equal_range(dateStart, dateEnd, *this);
    if (loc.first == loc.second){
        throw ModelException("DateTime::find", "Date "+toString()+" not in "
                             "supplied dates");
    }
    return (loc.first-dateStart);
}

/** Returns the index of the first instance of this date in the array
    of dates. Fails if the date does not exist. The supplied dates
    must be weakly increasing (a binary search is used)
    NOTE THIS METHOD IGNORES THE TIME*/
int DateTime::findIgnoringTime(const DateTimeArray& dates) const{
    DateTimeArray::const_iterator dateStart = dates.begin();
    DateTimeArray::const_iterator dateEnd = dates.end();
    pair<DateTimeArray::const_iterator, DateTimeArray::const_iterator> loc =
        equal_range(dateStart, dateEnd, *this, isLessWeak);
    if (loc.first == loc.second){
        throw ModelException("DateTime::find", "Date "+toString()+" not in "
                             "supplied dates");
    }
    return (loc.first-dateStart);
}

/** Returns true if this date is between upper and lower
    - the bounds are inclusive */
bool DateTime::within(const DateTime& lower,
                      const DateTime& upper) const
{
    return isGreaterOrEqual(lower) && isLessOrEqual(upper);
}

bool DateTime::isSubset(const DateTimeArray& full,
                        const DateTimeArray& part) {
    int p=0;
    for(int f=0; f<full.size() && p<part.size(); f++) {
        if (full[f]==part[p]) {
            p++;
        } else if (full[f]>part[p]) {
            return false;
        }
    }
    return (p==part.size());
}

/** Takes an array which is possibly not in order and with duplicates. In-place sorts and removes duplicates and returns itself. Sorting will use "DateTime::operator < ()",  unique() will use "operator ==", so both date and time are compared. Worst case complexity is the same as of std::sort (i.e. O(n log n)). */
DateTimeArray& DateTime::doSortUniq(DateTimeArray& dates)
{
    if (dates.size()==0) // do nothing
        return dates;
    std::sort(dates.begin(), dates.end()); // will use DateTime::operator < ()
    dates.erase(std::unique(dates.begin(), dates.end()), dates.end()); // move non-unique elements at the end and erase them
    return dates;
}

/** Returns the number of days between this and loDate. Equivalent to
    this.getDate() - loDate.getDate() */
int DateTime::daysDiff(const DateTime& loDate) const{
    return (date - loDate.date);
}

// DateTime::getDate() const is in the .hpp for performance reasons

int DateTime::getTime() const {
    return time;
}

void DateTime::moveIntoDay(){
    if (time  < START_OF_DAY_TIME){
        time = START_OF_DAY_TIME;
    } else if (time > END_OF_DAY_TIME){
        time = END_OF_DAY_TIME;
    }
}

/** returns a hashcode for the object based upon the date and time */
int DateTime::hashCode() const{
    return (date ^ time); // this is an "exclusive or"
}

/** Indicates whether some other object is "equal to" this one. Overridden
    for performance */
bool DateTime::equalTo(const IObject* obj) const{
    if (this == obj){
        return true;
    }
    if (!obj || obj->getClass() != TYPE){
        return false;
    }
    const DateTime* dt = STATIC_CAST(DateTime, obj);
    return (*dt == *this);
}

double DateTime::yearFrac(const DateTime &dateTo) const {
    return ((double)(dateTo.date - date))/ ((double) DAYS_PER_YEAR) +
        ((double)(dateTo.time - time))/ ((double) TIME_IN_YEAR);
}

typedef struct                         /* One entry for each month in cache */
{
    int             date;              /* Date of first day in cached month  */
    short int       month;             /* month= 1..12 */
    short int       year;              /* E.g. 1996  (NOT 96) */
} TDateCacheEntry;

static TDateCacheEntry gDateCacheArray[] = {
    {143905,1,1995},{143936,2,1995},{143964,3,1995},{143995,4,1995},
    {144025,5,1995},{144056,6,1995},{144086,7,1995},{144117,8,1995},
    {144148,9,1995},{144178,10,1995},{144209,11,1995},{144239,12,1995},
    {144270,1,1996},{144301,2,1996},{144330,3,1996},{144361,4,1996},
    {144391,5,1996},{144422,6,1996},{144452,7,1996},{144483,8,1996},
    {144514,9,1996},{144544,10,1996},{144575,11,1996},{144605,12,1996},
    {144636,1,1997},{144667,2,1997},{144695,3,1997},{144726,4,1997},
    {144756,5,1997},{144787,6,1997},{144817,7,1997},{144848,8,1997},
    {144879,9,1997},{144909,10,1997},{144940,11,1997},{144970,12,1997},
    {145001,1,1998},{145032,2,1998},{145060,3,1998},{145091,4,1998},
    {145121,5,1998},{145152,6,1998},{145182,7,1998},{145213,8,1998},
    {145244,9,1998},{145274,10,1998},{145305,11,1998},{145335,12,1998},
    {145366,1,1999},{145397,2,1999},{145425,3,1999},{145456,4,1999},
    {145486,5,1999},{145517,6,1999},{145547,7,1999},{145578,8,1999},
    {145609,9,1999},{145639,10,1999},{145670,11,1999},{145700,12,1999},
    {145731,1,2000},{145762,2,2000},{145791,3,2000},{145822,4,2000},
    {145852,5,2000},{145883,6,2000},{145913,7,2000},{145944,8,2000},
    {145975,9,2000},{146005,10,2000},{146036,11,2000},{146066,12,2000},
    {146097,1,2001},{146128,2,2001},{146156,3,2001},{146187,4,2001},
    {146217,5,2001},{146248,6,2001},{146278,7,2001},{146309,8,2001},
    {146340,9,2001},{146370,10,2001},{146401,11,2001},{146431,12,2001},
    {146462,1,2002},{146493,2,2002},{146521,3,2002},{146552,4,2002},
    {146582,5,2002},{146613,6,2002},{146643,7,2002},{146674,8,2002},
    {146705,9,2002},{146735,10,2002},{146766,11,2002},{146796,12,2002},
    {146827,1,2003},{146858,2,2003},{146886,3,2003},{146917,4,2003},
    {146947,5,2003},{146978,6,2003},{147008,7,2003},{147039,8,2003},
    {147070,9,2003},{147100,10,2003},{147131,11,2003},{147161,12,2003},
    {147192,1,2004},{147223,2,2004},{147252,3,2004},{147283,4,2004},
    {147313,5,2004},{147344,6,2004},{147374,7,2004},{147405,8,2004},
    {147436,9,2004},{147466,10,2004},{147497,11,2004},{147527,12,2004},
    {147558,1,2005},{147589,2,2005},{147617,3,2005},{147648,4,2005},
    {147678,5,2005},{147709,6,2005},{147739,7,2005},{147770,8,2005},
    {147801,9,2005},{147831,10,2005},{147862,11,2005},{147892,12,2005},
    {147923,1,2006},{147954,2,2006},{147982,3,2006},{148013,4,2006},
    {148043,5,2006},{148074,6,2006},{148104,7,2006},{148135,8,2006},
    {148166,9,2006},{148196,10,2006},{148227,11,2006},{148257,12,2006},
    {148288,1,2007},{148319,2,2007},{148347,3,2007},{148378,4,2007},
    {148408,5,2007},{148439,6,2007},{148469,7,2007},{148500,8,2007},
    {148531,9,2007},{148561,10,2007},{148592,11,2007},{148622,12,2007},
    {148653,1,2008},{148684,2,2008},{148713,3,2008},{148744,4,2008},
    {148774,5,2008},{148805,6,2008},{148835,7,2008},{148866,8,2008},
    {148897,9,2008},{148927,10,2008},{148958,11,2008},{148988,12,2008},
    {149019,1,2009},{149050,2,2009},{149078,3,2009},{149109,4,2009},
    {149139,5,2009},{149170,6,2009},{149200,7,2009},{149231,8,2009},
    {149262,9,2009},{149292,10,2009},{149323,11,2009},{149353,12,2009},
    {149384,1,2010},{149415,2,2010},{149443,3,2010},{149474,4,2010},
    {149504,5,2010},{149535,6,2010},{149565,7,2010},{149596,8,2010},
    {149627,9,2010},{149657,10,2010},{149688,11,2010},{149718,12,2010},
    {149749,1,2011},{149780,2,2011},{149808,3,2011},{149839,4,2011},
    {149869,5,2011},{149900,6,2011},{149930,7,2011},{149961,8,2011},
    {149992,9,2011},{150022,10,2011},{150053,11,2011},{150083,12,2011},
    {150114,1,2012},{150145,2,2012},{150174,3,2012},{150205,4,2012},
    {150235,5,2012},{150266,6,2012},{150296,7,2012},{150327,8,2012},
    {150358,9,2012},{150388,10,2012},{150419,11,2012},{150449,12,2012},
    {150480,1,2013},{150511,2,2013},{150539,3,2013},{150570,4,2013},
    {150600,5,2013},{150631,6,2013},{150661,7,2013},{150692,8,2013},
    {150723,9,2013},{150753,10,2013},{150784,11,2013},{150814,12,2013},
    {150845,1,2014},{150876,2,2014},{150904,3,2014},{150935,4,2014},
    {150965,5,2014},{150996,6,2014},{151026,7,2014},{151057,8,2014},
    {151088,9,2014},{151118,10,2014},{151149,11,2014},{151179,12,2014},
    {151210,1,2015},{151241,2,2015},{151269,3,2015},{151300,4,2015},
    {151330,5,2015},{151361,6,2015},{151391,7,2015},{151422,8,2015},
    {151453,9,2015},{151483,10,2015},{151514,11,2015},{151544,12,2015},
    {151575,1,2016},{151606,2,2016},{151635,3,2016},{151666,4,2016},
    {151696,5,2016},{151727,6,2016},{151757,7,2016},{151788,8,2016},
    {151819,9,2016},{151849,10,2016},{151880,11,2016},{151910,12,2016},
    {151941,1,2017},{151972,2,2017},{152000,3,2017},{152031,4,2017},
    {152061,5,2017},{152092,6,2017},{152122,7,2017},{152153,8,2017},
    {152184,9,2017},{152214,10,2017},{152245,11,2017},{152275,12,2017},
    {152306,1,2018},{152337,2,2018},{152365,3,2018},{152396,4,2018},
    {152426,5,2018},{152457,6,2018},{152487,7,2018},{152518,8,2018},
    {152549,9,2018},{152579,10,2018},{152610,11,2018},{152640,12,2018},
    {152671,1,2019},{152702,2,2019},{152730,3,2019},{152761,4,2019},
    {152791,5,2019},{152822,6,2019},{152852,7,2019},{152883,8,2019},
    {152914,9,2019},{152944,10,2019},{152975,11,2019},{153005,12,2019},
    {153036,1,2020},{153067,2,2020},{153096,3,2020},{153127,4,2020},
    {153157,5,2020},{153188,6,2020},{153218,7,2020},{153249,8,2020},
    {153280,9,2020},{153310,10,2020},{153341,11,2020},{153371,12,2020},
    {153402,1,2021},{153433,2,2021},{153461,3,2021},{153492,4,2021},
    {153522,5,2021},{153553,6,2021},{153583,7,2021},{153614,8,2021},
    {153645,9,2021},{153675,10,2021},{153706,11,2021},{153736,12,2021},
    {153767,1,2022},{153798,2,2022},{153826,3,2022},{153857,4,2022},
    {153887,5,2022},{153918,6,2022},{153948,7,2022},{153979,8,2022},
    {154010,9,2022},{154040,10,2022},{154071,11,2022},{154101,12,2022},
    {154132,1,2023},{154163,2,2023},{154191,3,2023},{154222,4,2023},
    {154252,5,2023},{154283,6,2023},{154313,7,2023},{154344,8,2023},
    {154375,9,2023},{154405,10,2023},{154436,11,2023},{154466,12,2023},
    {154497,1,2024},{154528,2,2024},{154557,3,2024},{154588,4,2024},
    {154618,5,2024},{154649,6,2024},{154679,7,2024},{154710,8,2024},
    {154741,9,2024},{154771,10,2024},{154802,11,2024},{154832,12,2024},
    {154863,1,2025},{154894,2,2025},{154922,3,2025},{154953,4,2025},
    {154983,5,2025},{155014,6,2025},{155044,7,2025},{155075,8,2025},
    {155106,9,2025},{155136,10,2025},{155167,11,2025},{155197,12,2025},
    {155228,1,2026},{155259,2,2026},{155287,3,2026},{155318,4,2026},
    {155348,5,2026},{155379,6,2026},{155409,7,2026},{155440,8,2026},
    {155471,9,2026},{155501,10,2026},{155532,11,2026},{155562,12,2026},
    {155593,1,2027},{155624,2,2027},{155652,3,2027},{155683,4,2027},
    {155713,5,2027},{155744,6,2027},{155774,7,2027},{155805,8,2027},
    {155836,9,2027},{155866,10,2027},{155897,11,2027},{155927,12,2027},
    {155958,1,2028},{155989,2,2028},{156018,3,2028},{156049,4,2028},
    {156079,5,2028},{156110,6,2028},{156140,7,2028},{156171,8,2028},
    {156202,9,2028},{156232,10,2028},{156263,11,2028},{156293,12,2028},
    {156324,1,2029},{156355,2,2029},{156383,3,2029},{156414,4,2029},
    {156444,5,2029},{156475,6,2029},{156505,7,2029},{156536,8,2029},
    {156567,9,2029},{156597,10,2029},{156628,11,2029},{156658,12,2029}
};

static int  leapCumDays[] = {
    -1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};

static int  cumDays[] = {
    -1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364};


static int  leapDays[] = {
    0, 31,  29,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31};
/* JAN  FEB  MAR  APR  MAY  JUN  JUL  AUG  SEP  OCT  NOV  DEC */

static int  days[] = {
    0, 31,  28,  31,  30,  31,  30,  31,  31,  30,  31,  30,  31};
/* JAN  FEB  MAR  APR  MAY  JUN  JUL  AUG  SEP  OCT  NOV  DEC */


static char* months[] = {"Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec"};

#define DATE_CACHE_NUM_ITEMS  ((int) (sizeof(gDateCacheArray)/ \
                               sizeof(TDateCacheEntry)))
/** we are using ALIB's TDate convention */
#define TDATE_BASE_YEAR     1601
#define DAYS_IN_1_YEAR      365
#define DAYS_IN_4_YEARS     1461
#define DAYS_IN_100_YEARS   (DAYS_IN_4_YEARS * 25 - 1)
#define DAYS_IN_400_YEARS   (DAYS_IN_100_YEARS * 4 + 1)
#define MONTHS_PER_YEAR     12
#define MAX_DAYS_PER_MONTH  31

/* A generally useful macro. */
#define IS_LEAP(year) (                 \
    (((year)%4 == 0) && ((year)%100 != 0)) || \
    ((year)%400 == 0) \
    )

/* Convert a date into its corresponding month - day - year
 * Note the month, day and year arguments are references and
 * WILL BE MODIFIED to contain the appropriate values */
static void dateToMDY(int date,
                      int &month,
                      int &day,
                      int &year)
{
    int         fourYearBlocks;
    int         count;
    int        *cumDaysp;

    /* Check if date is covered by cache  */
    if (gDateCacheArray[0].date <= date &&
        date <= gDateCacheArray[DATE_CACHE_NUM_ITEMS-1].date)
    {
        int i = (int) (date - gDateCacheArray[0].date);
        /* i = 12*(i/365) + (i%365)/29; */ /* make guess at proper index */
        i /= 29;                /* make guess at proper index */
        if (i>DATE_CACHE_NUM_ITEMS-1)
        {
            i = DATE_CACHE_NUM_ITEMS-1;
        }
        while (gDateCacheArray[i].date > date) /* can guess too high */
        {
            i--;
        }

        year  = gDateCacheArray[i].year;
        month = gDateCacheArray[i].month;
        day   = date - gDateCacheArray[i].date + 1;
    } else {
        if (date < 0){
            throw ModelException("DateToMDY", "Date ("+Format::toString(date)+
                                 ") is negative");
        }

        year = TDATE_BASE_YEAR;

        /* Get year */
        year += 400 * (date / DAYS_IN_400_YEARS); // Note: Integer division
        date %= DAYS_IN_400_YEARS;  // Note: date is in [0, DAYS_IN_400_YEARS-1]


        /* Go through this loop at most 3 times so that Dec 31 in the
         * year 2000, 2400, etc doesn't get moved to the year 2001, 2401. */
        for (count = 3; date >= DAYS_IN_100_YEARS && count > 0; count--){
            date -= DAYS_IN_100_YEARS;
            year += 100;
        }

        /* Dont have to make sure we go through at most 24 times since
         * the last 4 years (of 100) has *less* (or the same number of)
         * days than the other groups of 4 years.
         */
        if (date >= DAYS_IN_4_YEARS){
            fourYearBlocks = date/DAYS_IN_4_YEARS;
            date -= DAYS_IN_4_YEARS * fourYearBlocks;
            year += (int)fourYearBlocks << 2;   /* Multiply by 4 */
        }

        /* Go through this loop at most 3 times so that Dec 31 in a leap
         * year does not get moved to the next year. */
        for (count = 3; date >= DAYS_IN_1_YEAR && count > 0; count--)
        {
            date -= DAYS_IN_1_YEAR;
            year += 1;
        }

        /* Get month and date */

        /* date/32 is a good lower bound for month. */
        month = (date >> 5) + 1;

        if (IS_LEAP(year)){
            cumDaysp = leapCumDays + month;
        } else {
            cumDaysp = cumDays + month;
        }
        /* There is an extra increment and decrement of cumDaysp here, but
           it's necessary in order to set month correctly. */
        for ( ; date > *cumDaysp; month++){
            cumDaysp++;
        }
        day = date - *(--cumDaysp);
    }
}

/* Convert a month - day - year combination into the corresponding
 * date. Note the date argument is a reference and
 * WILL BE MODIFIED to contain the appropriate value */
static void mdyToDate(const int &monthToConvert,
                      const int &dayToConvert,
                      const int &yearToConvert,
                      int &date)
{
    static const string method ="mdyToDate";

    int      dt = 0;
    int      fourYearBlocks;
    int      year  = yearToConvert;
    int      month = monthToConvert;
    int      day   = dayToConvert;
    bool     isLeap;
    char     buffer[128];

    // Check if date is covered by cache
    if (gDateCacheArray[0].year <= year &&
        year <= gDateCacheArray[DATE_CACHE_NUM_ITEMS-1].year)
    {
        int i = 12*(year - gDateCacheArray[0].year) + month-1;  /* index */

        if (day < 1  ||  31 < day || month < 1  ||  12 < month)
        {
            sprintf(buffer, "invalid date %d-%d-%d", day, month, year);
            throw ModelException(method, buffer);
        }

        date = gDateCacheArray[i].date + day - 1;
        if (i < DATE_CACHE_NUM_ITEMS-1 && date > gDateCacheArray[i+1].date)
        {   // note: don't have to check last, as it's december
            sprintf(buffer, "invalid date %d-%d-%d", day, month, year);
            throw ModelException(method, buffer);
        }
    }
    else
    {
        year   = year - TDATE_BASE_YEAR;
        isLeap = IS_LEAP(yearToConvert);

        // Make sure day is in range
        if (day >= 1 && day <= 28) {
            /*EMPTY*/;                      /* Guaranteed to be OK */
            /* Avoid doing check below */
        }
        else if (day < 1 ||
                 (isLeap ? day > leapDays[month] : day > (days[month])))
        {
            sprintf(buffer, "invalid date %d-%d-%d",
                    dayToConvert, monthToConvert, yearToConvert);
            throw ModelException(method, buffer);
        }

        // Make sure month and year are in range
        if (month < 1 || month > MONTHS_PER_YEAR || yearToConvert < TDATE_BASE_YEAR)
        {
            sprintf(buffer, "invalid date %d-%d-%d",
                    dayToConvert, monthToConvert, yearToConvert);
            throw ModelException(method, buffer);
        }

        // Take years into account
        dt += DAYS_IN_400_YEARS * (year / 400); // Note: Integer division
        year %= 400; // Note: year is in [0, 399]


        if (year >= 100) {
            dt += DAYS_IN_100_YEARS * (year / 100); // Note: Integer division
            year %= 100; // Note: year is in [0, 99]
        }

        if (year >= 4) {
            fourYearBlocks = (int)(year>>2);       /* Divide by 4 */
            year -= (int)(fourYearBlocks<<2);       /* Multiply by 4 */
            dt   += fourYearBlocks * DAYS_IN_4_YEARS;
        }

        dt += year * DAYS_IN_1_YEAR;

        if (isLeap) {
            dt += leapCumDays[month-1] + day;
        }
        else {
            dt += cumDays[month-1] + day;
        }

        date = dt;
    }
}


// convert dd-mmm-yyyy into a date
static int convertDate(string date) {
    static const string method = "convertDate";

    try {
        int    day;
        int    month;
        int    year;
        int    dt;
        string d;
        string m;
        string y;

        /* see if string is remotely in date format */
        if (date.find("-") == std::string::npos || date.length() != 11)
        {
            string msg="date string (" + date + ") not in dd-mmm-yyyy format";
            throw ModelException(method, msg);
        }

        // chop it into dd, mmm & yyyy
        int  dash1 = date.find("-");
        int  dash2 = date.rfind("-");
        bool found;

        d = date.substr(0, 2);
        m = date.substr(dash1+1, 3);
        y = date.substr(dash2+1, 4);

        // day & year are easy
        sscanf(d.c_str(), "%d", &day);
        sscanf(y.c_str(), "%d", &year);

        month = 0;
        found = false;

        while (month < MONTHS_PER_YEAR && !found) {
            if (CString::equalsIgnoreCase(m, months[month])) {
                found = true;
            }
            month++; // if found, this converts from 0-1 to 1-12
        }

        if (!found) {
            throw ModelException(method, "couldn't find month ("+m+")");
        }

        mdyToDate(month, day, year, dt);

        // all OK
        return dt;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}


typedef struct  _THourMinSec
{
    int hour;                         /* In range [0,23] */
    int min;                          /* In range [0-59] */
    int sec;                          /* In range [0-59-] */
} THourMinSec;

static void timeToHMS(int          aTime,
                      THourMinSec *hms)   // (O)
{
    const int   halfsec   = DateTime::TIME_IN_DAY / (2 * 24 * 60 * 60);
    if (aTime < 0 || aTime >= DateTime::TIME_IN_DAY){
        throw ModelException("timeToHMS", "Invalid time: "+
                             Format::toString(aTime));
    }
    /* add on half a second to round to the nearest second */
    int time      = aTime + halfsec;
    hms->hour     = (time * 24L)/DateTime::TIME_IN_DAY;
    int unitHours = (hms->hour * DateTime::TIME_IN_DAY)/24L;
    hms->min      =  ((time - unitHours) * 24L * 60L)/DateTime::TIME_IN_DAY;
    int unitMins  = (hms->min * DateTime::TIME_IN_DAY) / (24L * 60L);
    hms->sec      = ((time - unitHours - unitMins) *
                     (24L * 60L * 60L)) / DateTime::TIME_IN_DAY;
    return;
}

#define MAX_LINE_LENGTH 25
string DateTime::timeFormat(int time) {
    string tm;

    if (time == DateTime::START_OF_DAY_TIME) {
        tm = DateTime::START_OF_DAY;
    }
    else if (time == DateTime::END_OF_DAY_TIME) {
        tm = DateTime::END_OF_DAY;
    }
    else if (time == DateTime::BEFORE_EX_DIV_TIME) {
        tm = DateTime::BEFORE_EX_DIV;
    }
    else if (time == DateTime::EX_DIV_BEFORE_START_TIME) {
        tm = DateTime::EX_DIV_BEFORE_START;
    }
    else {
        THourMinSec hms = {0};
        timeToHMS(time, &hms);
        char   buffer[MAX_LINE_LENGTH];

        sprintf(buffer, "%02d:%02d:%02d", hms.hour, hms.min, hms.sec);
        tm = buffer;
    }
    return tm;
}

string DateTime::dateFormat(int date) {
    char   buffer[MAX_LINE_LENGTH];
    int    month;
    int    day;
    int    year;

    dateToMDY(date, month, day, year);
    sprintf(buffer, "%02d-%s-%d", day, months[month-1], year);
    return string(buffer);
}

/** returns string representation of date - good for error messages */
string DateTime::toString() const{
    char   buffer[MAX_LINE_LENGTH];
    sprintf(buffer, "%s %s",
            DateTime::dateFormat(date).c_str(),  // date portion
            DateTime::timeFormat(time).c_str()); // time portion
    return string(buffer);
}

#ifdef DEBUG
//// handy for debugging. DO NOT USE ANYWHERE ELSE
const char* DateTime::p() const{
    // certainly not threadsafe! - I tried using new but it seemed to cause
    // problems with the debugger
    static char buffer[25];
    try{
        const string& str = toString();
        strcpy(buffer, str.c_str());
    } catch (exception&){
        strcpy(buffer, "invalid datetime");
    }
    return buffer;
}
#endif

/* Interval constructor */
DateTime::Interval::Interval(int dateDiff, int timeDiff):
    dateDiff(dateDiff), timeDiff(timeDiff){}

//// comparison methods
bool DateTime::Interval::isGreaterOrEqual(const Interval& interval) const{
    return (this->dateDiff > interval.dateDiff ||
            ((this->dateDiff == interval.dateDiff) &&
             (this->timeDiff >= interval.timeDiff)));
}
bool DateTime::Interval::isLess(const Interval& interval) const{
    return (this->dateDiff < interval.dateDiff ||
            ((this->dateDiff == interval.dateDiff) &&
             (this->timeDiff < interval.timeDiff)));
}
bool DateTime::Interval::isGreater(const Interval& interval) const{
    return (this->dateDiff > interval.dateDiff ||
            ((this->dateDiff == interval.dateDiff) &&
             (this->timeDiff > interval.timeDiff)));
}

/** Returns a date-time midway between this date and the supplied date */
DateTime DateTime::midPoint(const DateTime& date) const{
    int daysDiff = date.date - this->date;
    int timeDiff = date.time - this->time;
    int dayMidPoint = daysDiff/2;
    if (dayMidPoint * 2 != daysDiff){
        /* this is open to interpretation (eg midpoint between 2004/1/1 SOD and
           2004/1/2 EOD is what?) */
        if (this->time+(timeDiff+TIME_IN_DAY)/2 >= TIME_IN_DAY){
            dayMidPoint++;
        } else {
            timeDiff += TIME_IN_DAY;
        }
    }
    return DateTime(this->date+dayMidPoint, this->time+timeDiff/2);
}

/** Returns the interval between this date and date
    ie "this - dateToSubtract" */
DateTime::Interval DateTime::subtract(const DateTime& dateToSubtract) const{
    int dateDiff = date - dateToSubtract.date;
    int timeDiff = time -  dateToSubtract.time;
    if (time < dateToSubtract.time){
        dateDiff--;
        timeDiff = TIME_IN_DAY + timeDiff; // as timeDiff is < 0 here
    }
    return Interval(dateDiff, timeDiff);
}

/** Adds given interval to date */
DateTime DateTime::add(const Interval& interval) const{
    int newDate = date + interval.dateDiff;
    int newTime = time + interval.timeDiff;
    if (newTime >= TIME_IN_DAY){
        newDate++;
        newTime -= TIME_IN_DAY;
    }
    return DateTime(newDate, newTime);
}

/* returns the day of the week with 0 = Saturday ... 6 = Friday */
int DateTime::getWeekday() const{
    return ((date + 2 ) % 7);
}

/* is day on a weekend? */
bool DateTime::isWeekend() const{
    return isWeekend(date);
}

/** Defined such that if dt is a DateTime,
    dt.isWeekend() == DateTime::isWeekend(dt.getDate()) */
bool DateTime::isWeekend(int date){
    return ( ( (1 << ((date) % 7)) & EDG_WEEKEND_STANDARD ) != 0 );
}

/** is day in a leap year ? */
bool DateTime::isLeap() const {
    try {
        DateTime::MonthDayYear mdy = toMDY();

        return DateTime::isLeapYear(mdy.year);
    }
    catch (exception &e) {
        throw ModelException(&e, "DateTime::isLeap");
    }
}

/** is day last day in month ? */
bool DateTime::isEndOfMonth() const
{
    try {
        bool isEOM;
        DateTime::MonthDayYear mdy = toMDY();

        if (isLeap()) {
            isEOM = (mdy.day == leapDays[mdy.month]);
        } else {
            isEOM = (mdy.day == days[mdy.month]);
        }
        return isEOM;
    }
    catch (exception &e) {
        throw ModelException(&e, "DateTime::isLeap");
    }
}

/** returns the last day of the month that the date falls in.
    Time of day remains the same. ignoreLeapYears = true sets
    to February 28th in leap years */
DateTime DateTime::returnEndOfMonth(bool ignoreLeapYears) const
{
    DateTime::MonthDayYear mdy(toMDY());
    mdy.day = mdy.daysInMonth(ignoreLeapYears);
    return mdy.toDateTime(getTime());
}

/** returns the first day of the month that this date falls in. */
DateTime DateTime::firstOfMonth() const {
    MonthDayYear mdy = toMDY();
    mdy.day = 1;
    return mdy.toDateTime(getTime());
}

// public interface onto month-day-year
DateTime::MonthDayYear::MonthDayYear(int day, int month, int year) :
    day(day), month(month), year(year)
{
    if ( day < 1 || day > 31 || month < 1 || month > 12 ||
         year < TDATE_BASE_YEAR ) {
        throw ModelException("MonthDayYear::MonthDayYear",
                             "Date (" + Format::toString(day) +
                             "/"      + Format::toString(month) +
                             "/"      + Format::toString(year) + ") is invalid");
    }
}

DateTime DateTime::MonthDayYear::toDateTime(int time) {
    static const string method = "DateTime::MonthDayYear::toDateTime";
    try {
        int date;

        normalize();

        mdyToDate(month, day, year, date);

        return DateTime(date, time);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    }
}

/* Returns the number of whole months between this dateTime and
 * the mdyTo parameter. The result is > 0 if mdyTo is greater
 * than this dateTime, and < 0 otherwise
 *
 * CAUTION:
 * This function is NOT simmetrical: In general,
 * dateA.wholeMonthsBetween(dateB) is the same as dateB.wholeMonthsBetween(dateA)
 * BUT if we are going forwards, going from the last day of one month
 * (e.g. 31-Dec-2003) to the last day of another month (e.g. 30-Jun-2006) does
 * count the last month in (ie, result = 30 month) - However, going backwards it
 * does not, as if one day was missing (ie, result = -29).
 * This behaviour is consistent with existing code in other areas and has
 * been maintained until proven wrong */
int DateTime::MonthDayYear::wholeMonthsBetween(MonthDayYear mdyTo) {
    static const string method = "DateTime::MonthDayYear::wholeMonthsTill";
    try {
        int numMonths;

        numMonths  = (mdyTo.year - year) * 12;
        numMonths += (mdyTo.month - month);
        if ((numMonths > 0) && // Going forwards
            (mdyTo.day < day) &&  // This month is not complete
            // unless we are on the last day of the month
            (DateTime::isLeapYear(mdyTo.year) ? !(leapDays[mdyTo.month] == mdyTo.day) :
             !(days[mdyTo.month] == mdyTo.day)))
        {
            // Going forwards this month is not complete and we are
            // not going from 29 Feb to the 28 Feb of a non-leap year
            numMonths--;
        }
        else if ((numMonths < 0) && (mdyTo.day > day)) {
            // Going backwards this month is not complete
            numMonths++;
        }
        // else nothing to do
        return numMonths;
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    }
}

/**  Answers the number of days in the month of the date given by this MonthDayYear */
int DateTime::MonthDayYear::daysInMonth(bool ignoreLeapYears) const {
    if (ignoreLeapYears || !DateTime::isLeapYear(year))
        return days[month];
    else
        return leapDays[month];      
}

void DateTime::MonthDayYear::normalize() {
    static const string method = "DateTime::MonthDayYear::normalize";
    try {
        year += month / MONTHS_PER_YEAR; // Note: Integer division
        month %= MONTHS_PER_YEAR;  // Note: if MONTHS_PER_YEAR=12, month is in [-11,11]

        if (month < 1) {
            month += MONTHS_PER_YEAR;
            year--;
        }

        if (day < 1 || day > MAX_DAYS_PER_MONTH) {
            throw ModelException(method, "invalid day");
        }

        if (IS_LEAP(year)) {
            if (day > leapDays[month]) {
                day = leapDays[month];
            }
        }
        else {
            if (day > days[month]) {
                day = days[month];
            }
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}


DateTime::MonthDayYear DateTime::toMDY() const {
    static const string method = "DateTime::toMDY";
    try {
        int month;
        int day;
        int year;

        dateToMDY(date, month, day, year);

        return MonthDayYear(day, month, year);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    }
}

/* Creates DateTime from IR (Alib) style date ie yyyymmdd. Hard coded to
   start of day */
DateTime DateTime::fromIrDate(long irDate){
    int date;
    int month;
    int day;
    int year;

    year  =  irDate/10000; /* Truncation forced */
    irDate -= 10000 * year;
    month = irDate/100;
    day   = (irDate - 100 * month);

    mdyToDate(month, day, year, date);
    return DateTime(date, START_OF_DAY_TIME);
}

/* converts to an IR (Alib) style date ie yyyymmdd */
long DateTime::toIrDate(void) const {
    MonthDayYear mdy = toMDY();
    return (mdy.year * 100 + mdy.month )* 100 + mdy.day;
}

/** Returns the year (eg 2003) in which this date falls.
    To do: reuse in conjunction with toMDY */
int DateTime::getYear() const{
    /** Quantify suggests that the cache is actually slower than this
        simple integer arithmetic */
    int date = this->date;
    int year = TDATE_BASE_YEAR;
    int count;

    while (date >= DAYS_IN_400_YEARS){
        date -= DAYS_IN_400_YEARS;
        year += 400;
    }
    /* Go through this loop at most 3 times so that Dec 31 in the
     * year 2000, 2400, etc doesn't get moved to the year 2001, 2401. */
    for (count = 3; date >= DAYS_IN_100_YEARS && count > 0; count--){
        date -= DAYS_IN_100_YEARS;
        year += 100;
    }
    /* Dont have to make sure we go through at most 24 times since
     * the last 4 years (of 100) has *less* (or the same number of)
     * days than the other groups of 4 years.
     */
    if (date >= DAYS_IN_4_YEARS){
        int fourYearBlocks = date/DAYS_IN_4_YEARS;
        date -= DAYS_IN_4_YEARS * fourYearBlocks;
        year += (int)fourYearBlocks << 2;   /* Multiply by 4 */
    }
    /* Go through this loop at most 3 times so that Dec 31 in a leap
     * year does not get moved to the next year. */
    for (count = 3; date >= DAYS_IN_1_YEAR && count > 0; count--) {
        date -= DAYS_IN_1_YEAR;
        year += 1;
    }
    return year;
}

/** Returns DateTime for Dec 31st EOD in supplied year */
DateTime DateTime::endOfYear(int year) {
    /** Quantify suggests that the cache is actually slower than this
        simple integer arithmetic */
    const int month = 12;
    const int day =   31;
    bool isLeap = IS_LEAP(year);
    year   = year - TDATE_BASE_YEAR;
    int dt = 0;

    // Take years into account
    dt += DAYS_IN_400_YEARS * (year / 400); // Note: Integer division
    year %= 400; // Note: year is in [0, 399]


    if (year >= 100) {
        dt += DAYS_IN_100_YEARS * (year / 100); // Note: Integer division
        year %= 100; // Note: year is in [0, 99]
    }

    if (year >= 4) {
        int fourYearBlocks = (int)(year>>2);       /* Divide by 4 */
        year -= (int)(fourYearBlocks<<2);       /* Multiply by 4 */
        dt   += fourYearBlocks * DAYS_IN_4_YEARS;
    }

    dt += year * DAYS_IN_1_YEAR;

    if (isLeap) {
        dt += leapCumDays[month-1] + day;
    }
    else {
        dt += cumDays[month-1] + day;
    }
    return DateTime(dt,  END_OF_DAY_TIME);
}

/** Counts the number of leap year days between a start date
    (inclusive) and an end date (exclusive).  Needed for
    Japanese Government Bonds and some floating rate notes
*/
int DateTime::countLeapYearDays(const DateTime& d1, const DateTime& d2) {
    static const string method = "DateTime::countLeapYearDays";
    try {
        int      days;
        DateTime date1;
        DateTime date2;

        if (d1.isGreater(d2)) {
            date1 = d2;
            date2 = d1;
        }
        else {
            date1 = d1;
            date2 = d2;
        }

        DateTime::MonthDayYear mdy1 = date1.toMDY();
        DateTime::MonthDayYear mdy2 = date2.toMDY();

        if (date1.equals(date2)) {
            if (mdy1.month == 2 && mdy1.day == 29) {
                days = 1;
            }
            else {
                days = 0;
            }
        }
        else {
            days = 0;

            // If the first year is a leap year, count the extra day if it
            // falls on or after the start date
            if (DateTime::isLeapYear(mdy1.year)) {
                mdy1.month = 2;
                mdy1.day = 29;

                DateTime temp = mdy1.toDateTime(d1.getTime());

                if (!date1.isGreater(temp)) {
                    days++;
                }
            }

            // If the last year is a leap year, count the extra day if it
            // falls before the end date

            // However in addition, if the years are the same and a leap year,
            // then subtract one from the result.

            if (DateTime::isLeapYear(mdy2.year)) {
                mdy2.month = 2;
                mdy2.day = 29;

                DateTime temp = mdy2.toDateTime(d2.getTime());

                if (temp.isLess(date2)) {
                    days++;
                }

                if (mdy1.year == mdy2.year)
                {
                    --days;
                }
            }

            // Otherwise, loop through the intervening years and check each one
            // to see if it's a leap year

            for (int year = mdy1.year + 1; year < mdy2.year; year++) {
                if (DateTime::isLeapYear(year)) {
                    days++;
                }
            }
        }

        return days;
    }
    catch (exception& e) {
        throw ModelException(&e, method, "Failed");
    }
}

/** write object out to writer */
void DateTime::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get() && !empty()){
            // don't write out empty DateTimes
            writer->write(toString());
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "DateTime::write");
    }
}


// convert hh:mm:ss into a time
int DateTime::timeConvert(string time) {
    static const string method = "DateTime::timeConvert";

    int  hr;
    int  min;
    int  sec;
    int  tm;
    char buffer[128];

    // try special strings for start & end of day first
    if (time == DateTime::START_OF_DAY) {
        tm = DateTime::START_OF_DAY_TIME;
    }
    else if (time == DateTime::END_OF_DAY) {
        tm = DateTime::END_OF_DAY_TIME;
    }
    else if (time == DateTime::BEFORE_EX_DIV) {
        tm = DateTime::BEFORE_EX_DIV_TIME;
    }
    else if (time == DateTime::EX_DIV_BEFORE_START) {
        tm = DateTime::EX_DIV_BEFORE_START_TIME;
    }
    else {
        /* see if string is remotely in time format */
        if (time.find(":") == std::string::npos)
        {
            string msg = "time string (" + time + ") not in hh:mm:ss format";
            throw ModelException(method, msg);
        }

        sscanf(time.c_str(), "%d:%d:%d", &hr, &min, &sec);

        /* because you can't trust anybody */
        if (hr < 0 || hr > 23)
        {
            sprintf(buffer, "hour (%d) out of bounds (0 <= hr < 24) ", hr);
            throw ModelException(method, buffer);
        }

        if (min < 0 || min > 59)
        {
            sprintf(buffer, "minutes (%d) out of bounds (0 <= min < 60)", min);
            throw ModelException(method, buffer);
        }

        if (sec < 0 || sec > 59)
        {
            sprintf(buffer, "seconds (%d) out of bounds (0 <= sec < 60)", sec);
            throw ModelException(method, buffer);
        }

        tm = hr*60*60 + min*60 + sec;

        if (tm < 0)
        {
            sprintf(buffer, "time (%d) is negative !", tm);
            throw ModelException(method, buffer);
        }

        if (tm > DateTime::TIME_MAX)
        {
            sprintf(buffer, "time (%d) too big !", tm);
            throw ModelException(method, buffer);
        }
    }
    // all OK
    return tm;
}

/** populate an empty object from reader */
void DateTime::import(Reader::Node* elem, Reader* reader) {
    try {
        string    formatted = elem->value();
        if (!formatted.empty()){
            // empty datetimes => no value specified ie empty() was true
            string    dt = formatted.substr(0, 11);
            string    tm = formatted.substr(12);
            date = convertDate(dt);
            time = DateTime::timeConvert(tm);
        }
    }
    catch (exception& e) {
        throw ModelException(&e, "DateTime::import");
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void DateTime::outputWrite(const string& linePrefix,
                           const string& prefix, ostream& stream) const{
    stream << linePrefix << prefix << ": " << toString() << endl;
}

DateTime::Date::Date(int date): CObject(TYPE), date(date){}

//// same as empty constructor to DateTime
DateTime::Date::Date(): CObject(TYPE), date(0){} // to do: move to source file

int DateTime::Date::getDate() const{
    return date;
}

/** returns true if dates are the same */
bool DateTime::Date::equals(const Date& dateToCompare) const{
    return date == dateToCompare.date;
}

/** returns true if this date is greater than date supplied */
bool DateTime::Date::isGreater(const Date& dateToCompare) const{
    return date > dateToCompare.date;
}

/** returns true if this date is greater than or equal to date supplied */
bool DateTime::Date::isGreaterOrEqual(const Date& dateToCompare) const{
    return date >= dateToCompare.date;
}

/** returns string representation of date */
string DateTime::Date::toString() const{
    int    month;
    int    day;
    int    year;
    char   buffer[MAX_LINE_LENGTH];

    dateToMDY(date, month, day, year);
    sprintf(buffer, "%02d-%s-%d", day, months[month-1], year);
    return string(buffer);
}

/** write object out to writer */
void DateTime::Date::write(const string& tag, Writer* writer) const{
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(toString());
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "DateTime::Date::write");
    }
}

/** populate an empty object from reader */
void DateTime::Date::import(Reader::Node* elem, Reader* reader){
    static const string routine("DateTime::Date::import");
    try {
        string formatted = elem->value();
        date = convertDate(formatted);
    }
    catch (exception& e) {
        throw ModelException(e, routine, "Failed");
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void DateTime::Date::outputWrite(const string& linePrefix,
                                 const string& prefix,
                                 ostream& stream) const{
    stream << linePrefix << prefix << ": " << toString() << endl;
}

DateTime::Time::Time(int time): CObject(TYPE), time(time){}

int DateTime::Time::getTime() const{
    return time;
}

/** returns string representation of time */
string DateTime::Time::toString() const{
    return DateTime::timeFormat(time);
}

/** write object out to writer */
void DateTime::Time::write(const string& tag, Writer* writer) const{
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(toString());
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(e, "DateTime::Time::write");
    }
}

/** populate an empty object from reader */
void DateTime::Time::import(Reader::Node* elem, Reader* reader){
    static const string routine("DateTime::Time::import");
    try {
        string formatted = elem->value();
        time = DateTime::timeConvert(formatted);
    }
    catch (exception& e) {
        throw ModelException(e, routine, "Failed");
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void DateTime::Time::outputWrite(const string& linePrefix,
                                 const string& prefix,
                                 ostream& stream) const{
    stream << linePrefix << prefix << ": " << toString() << endl;
}

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<DateTime>::toIObject(
    const DateTime& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<DateTime>::toIObject(DateTime& value){
    return DateTimeSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DateTime */
const DateTime& arrayObjectCast<DateTime>::fromIObject(IObjectSP& value){
    DateTime *dtPtr = DYNAMIC_CAST(DateTime, value.get());
    return *dtPtr;
}

// explicit clone for arrays of datetimes - for performance
IObject* arrayClone<DateTime>::clone(const CArray* arrayToClone){
    const DateTimeArray& theArray =
        static_cast<const DateTimeArray&>(*arrayToClone);
    return new DateTimeArray(theArray);
}

DEFINE_TEMPLATE_TYPE(DateTimeArray);

/** for reflection */
DateTimeAddin::DateTimeAddin():  CObject(TYPE){}

/** Invoked when this class is 'loaded' */
void DateTimeAddin::load(CClassSP& clazz){
    REGISTER(DateTimeAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultDateTimeAddin);
    FIELD(dateTime, "A date-time");
    FIELD_MAKE_OPTIONAL(dateTime);
}

IObject* DateTimeAddin::defaultDateTimeAddin(){
    return new DateTimeAddin();
}


CClassConstSP const DateTimeAddin::TYPE = CClass::registerClassLoadMethod(
    "DateTimeAddin", typeid(DateTimeAddin), load);

/** specialisations of arrayObjectCast for DateTimeCluster */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<DateTimeArray>::toIObject(
    const DateTimeArray& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<DateTimeArray>::toIObject(DateTimeArray& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DateTimeArray */
const DateTimeArray& arrayObjectCast<DateTimeArray>::fromIObject(
    IObjectSP& value){
    DateTimeArray *dtPtr = DYNAMIC_CAST(DateTimeArray, value.get());
    return *dtPtr;
}

// explicit clone for arrays of datetime arrays - for performance
IObject* arrayClone<DateTimeArray>::clone(const CArray* arrayToClone){
    const DateTimeCluster& theCluster =
        static_cast<const DateTimeCluster&>(*arrayToClone);
    return new DateTimeCluster(theCluster);
}

DEFINE_TEMPLATE_TYPE_WITH_NAME("DateTimeArrayArray", DateTimeCluster);

/** for reflection */
DateTimeArrayAddin::DateTimeArrayAddin():  CObject(TYPE){}

/** Invoked when Class is 'loaded' */
void DateTimeArrayAddin::load(CClassSP& clazz){
    REGISTER(DateTimeArrayAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultDateTimeArrayAddin);
    FIELD(dateTimeArray, "An array of date-times");
    FIELD_MAKE_OPTIONAL(dateTimeArray);
}

IObject* DateTimeArrayAddin::defaultDateTimeArrayAddin(){
    return new DateTimeArrayAddin();
}

CClassConstSP const DateTimeArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "DateTimeArrayAddin", typeid(DateTimeArrayAddin), load);

///// class for addin function that takes in a separate date and time
class DateTimeAddin2: public CObject{
public:
    static CClassConstSP const TYPE;
    // parameters for addin
    DateTime::DateSP   date;
    DateTime::TimeSP   time;

    DateTimeAddin2(): CObject(TYPE), date(0), time(0){}

    static IObjectSP makeDateTime(DateTimeAddin2* params){
        if (!params->date && !params->time){
            // this is useful for optional datetimes
            return CNull::create();
        }
        if (!params->date || !params->time){
            throw ModelException("EDateTimeAddin2::makeDateTime",
                                 "Both or neither parameters must be supplied");
        }
        return IObjectSP(new DateTime(params->date->getDate(),
                                      params->time->getTime()));
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DateTimeAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDateTimeAddin2);
        FIELD(date, "Date");
        FIELD_MAKE_OPTIONAL(date);
        FIELD(time, "Time");
        FIELD_MAKE_OPTIONAL(time);
        Addin::registerClassObjectMethod("DATE_TIME",
                                         Addin::UTILITIES,
                                         "Builds date-time from date and time",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)makeDateTime);
    }

    static IObject* defaultDateTimeAddin2(){
        return new DateTimeAddin2();
    }
};

CClassConstSP const DateTimeAddin2::TYPE = CClass::registerClassLoadMethod(
    "DateTimeAddin2", typeid(DateTimeAddin2), load);

typedef array<smartPtr<DateTime::Date>, DateTime::Date> DateArray;
typedef array<smartPtr<DateTime::Time>, DateTime::Time> TimeArray;

DEFINE_TEMPLATE_TYPE_WITH_NAME("DateTime::DateArray", DateArray);

DEFINE_TEMPLATE_TYPE_WITH_NAME("DateTime::TimeArray", TimeArray);

///// class for addin function that takes in a separate dates and times
class DateTimeArrayAddin2: public CObject{
public:
    static CClassConstSP const TYPE;
    // parameters for addin
    DateTime::DateArraySP   dates;
    DateTime::TimeArraySP   times;

    DateTimeArrayAddin2(): CObject(TYPE) {}

    static IObjectSP makeDateTime(DateTimeArrayAddin2* params){
        if (params->times->size() != 1 &&
            params->times->size() != params->dates->size()){
            throw ModelException("DateTimeArrayAddin2::makeDateTime",
                                 "Number of cells indicating the time must be"
                                 " 1 or equal to the number of dates");
        }
        DateTimeArraySP dtArray(new DateTimeArray(params->dates->size()));
        for (int i = 0; i < params->dates->size(); i++){
            DateTime::TimeSP& time = (params->times->size() == 1)?
                (*params->times)[0]: (*params->times)[i];
            (*dtArray)[i] =  DateTime((*params->dates)[i]->getDate(),
                                      time->getTime());
        }
        return dtArray;
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DateTimeArrayAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDateTimeArrayAddin2);
        FIELD(dates, "array of dates");
        FIELD(times, "a single time or an array");
        Addin::registerClassObjectMethod("DATE_TIME_ARRAY",
                                         Addin::UTILITIES,
                                         "Builds date-time array from "
                                         "separate dates and times",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)makeDateTime);
    }

    static IObject* defaultDateTimeArrayAddin2(){
        return new DateTimeArrayAddin2();
    }
};

CClassConstSP const DateTimeArrayAddin2::TYPE =CClass::registerClassLoadMethod(
    "DateTimeArrayAddin2", typeid(DateTimeArrayAddin2), load);

//// class for addin function that merges array of dates
class DTArrayMergeAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    // parameters for addin
    DateTimeArraySP   dates1;
    DateTimeArraySP   dates2;
    // be generous to the spreadsheet user
    DateTimeArraySP   dates3;
    DateTimeArraySP   dates4;

    DTArrayMergeAddin(): CObject(TYPE) {}

    static IObjectSP merge(DTArrayMergeAddin* params){
        vector<DateTimeArray> arrays(1, *params->dates1.get());
        if (params->dates2.get()){
            arrays.push_back(*params->dates2.get());
        }
        if (params->dates3.get()){
            arrays.push_back(*params->dates3.get());
        }
        if (params->dates4.get()){
            arrays.push_back(*params->dates4.get());
        }
        for (unsigned int i = 0; i < arrays.size(); i++){
            DateTime::ensureIncreasing(arrays[i], "list ", false);
            // move non-unique elements at the end and erase them
            // note this operates in-place, hence the local copies
            arrays[i].erase(std::unique(arrays[i].begin(), arrays[i].end()), arrays[i].end());
        }
        IObjectSP mergedDates(new DateTimeArray(DateTime::merge(arrays)));
        return mergedDates;
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DTArrayMergeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDTArrayMergeAddin);
        FIELD(dates1, "array of dates");
        FIELD(dates2, "optional array of dates");
        FIELD_MAKE_OPTIONAL(dates2);
        FIELD(dates3, "optional array of dates");
        FIELD_MAKE_OPTIONAL(dates3);
        FIELD(dates4, "optional array of dates");
        FIELD_MAKE_OPTIONAL(dates4);
        Addin::registerClassObjectMethod("DATETIME_MERGE",
                                         Addin::UTILITIES,
                                         "Merges up to four ordered "
                                         "date-time arrays",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)merge);
    }

    static IObject* defaultDTArrayMergeAddin(){
        return new DTArrayMergeAddin();
    }
};

CClassConstSP const DTArrayMergeAddin::TYPE =CClass::registerClassLoadMethod(
    "DTArrayMergeAddin", typeid(DTArrayMergeAddin), load);

/** stores unique, sorted version of dl */
DateTimeList::DateTimeList(const DateTimeArray* dates):
    CObject(TYPE), dl(copy(dates)) {

    DateTimeArraySP unique(new DateTimeArray(0));
    sort(dl->begin(), this->dl->end());

    // strip duplicates
    if (!dl->empty()) {
        unique->push_back((*dl)[0]);
    }

    for (int i = 1; i < dl->size(); i++) {
        if ((*dl)[i] > (*dl)[i-1]) {
            unique->push_back((*dl)[i]);
        }
    }

    dl = unique;
}

/** scale by factor x */
void DateTimeList::scale(double x) {
    // empty - just what would it do ?
}

/** add DateTimeList to this result (Implementation of
    CombinableResult) */
void DateTimeList::add(const CombinableResult& x, double scaleFactor) {
    const DateTimeList& toAdd = *DYNAMIC_CAST(DateTimeList, &x);
    DateTimeArray unique = DateTime::merge(*dl, *toAdd.dl);
    dl = DateTimeArraySP(copy(&unique));
}

DateTimeList::DateTimeList() : CObject(TYPE) {}

class DateTimeListHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DateTimeList, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultDateTimeList);
        FIELD(dl, "dates");
    }
    static IObject* defaultDateTimeList(){
        return new DateTimeList();
    }
};

CClassConstSP const DateTimeList::TYPE = CClass::registerClassLoadMethod(
    "DateTimeList", typeid(DateTimeList), DateTimeListHelper::load);

/** Class for Addin functions that require two datetime's */
class DateTimeCompareAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    DateTime dateTime1;
    DateTime dateTime2;

private:
    /** for reflection */
    DateTimeCompareAddin():
        CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DateTimeCompareAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(dateTime1, "'lhs' date-time");
        FIELD(dateTime2, "'rhs' date-time");
        Addin::registerClassBoolMethod("DATETIME_EQUAL",
                                       Addin::UTILITIES,
                                       "Returns true if 'dateTime1 == dateTime2'",
                                       TYPE,
                                       (Addin::BoolMethod*)equal);
        Addin::registerClassBoolMethod("DATETIME_LESS",
                                       Addin::UTILITIES,
                                       "Returns true if 'dateTime1 < dateTime2'",
                                       TYPE,
                                       (Addin::BoolMethod*)less);
        Addin::registerClassBoolMethod("DATETIME_GREATER",
                                       Addin::UTILITIES,
                                       "Returns true if 'dateTime1 > dateTime2'",
                                       TYPE,
                                       (Addin::BoolMethod*)greater);
    }

    static bool equal(DateTimeCompareAddin* params){
        return params->dateTime1.equals(params->dateTime2);
    }

    static bool greater(DateTimeCompareAddin* params){
        return params->dateTime1.isGreater(params->dateTime2);
    }

    static bool less(DateTimeCompareAddin* params){
        return params->dateTime1.isLess(params->dateTime2);
    }

    static IObject* defaultCtor(){
        return new DateTimeCompareAddin();
    }
};

CClassConstSP const DateTimeCompareAddin::TYPE = CClass::registerClassLoadMethod(
    "DateTimeCompareAddin", typeid(DateTimeCompareAddin), load);

// we never put a public/private interface on DateTime, yet we don't want
// to expose the two ints that make up its data. Too late to change, so for
// DR Interface purposes provide a proxy which acts as the interface we
// want to expose
class DateTimeProxy: public CObject {
public:
    static CClassConstSP const TYPE;

private:
    string date;
    string time;

    static IObjectSP toDateTime(const IObjectConstSP& obj) {
        const DateTimeProxy* proxy = dynamic_cast<const DateTimeProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("DateTimeProxy::toDateTime",
                                 "object is not a DateTimeProxy");
        }

        return IObjectSP(new DateTime(proxy->date, proxy->time));
    }

    static IObjectSP fromDateTime(const IObjectConstSP& obj) {
        const DateTime* dt = dynamic_cast<const DateTime*>(obj.get());
        if (!dt) {
            throw ModelException("DateTimeProxy::fromDateTime",
                                 "object is not a DateTime");
        }

        return IObjectSP(new DateTimeProxy(DateTime::dateFormat(dt->getDate()),
                                           DateTime::timeFormat(dt->getTime())));
    }

    DateTimeProxy(const string& date, const string& time) : CObject(TYPE),
                                                            date(date), time(time) {}

    /** for reflection */
    DateTimeProxy():  CObject(TYPE){}

    static IObject* defaultDateTimeProxy(){
        return new DateTimeProxy();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DateTimeProxy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDateTimeProxy);
        registerObjectProxy(DateTime::TYPE,
                            DateTimeProxy::TYPE,
                            fromDateTime,
                            toDateTime);
        FIELD(date, "date (in DD-MMM-YYYY format)");
        FIELD(time, "time (either HH-MM-SS format or BEX, XBS, SOD, EOD)");
    }
};

CClassConstSP const DateTimeProxy::TYPE = CClass::registerClassLoadMethod(
    "DateTimeProxy", typeid(DateTimeProxy), load);


class DateTimeHelper{
public:
    /** Invoked when Date class is 'loaded' */
    static void loadDate(CClassSP& clazz){
        REGISTER(DateTime::Date, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDate);
        FIELD(date, "date");
    }

    static IObject* defaultDate(){
        return new DateTime::Date(0);
    }

    /** Invoked when Time class is 'loaded' */
    static void loadTime(CClassSP& clazz){
        REGISTER(DateTime::Time, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTime);
        FIELD(time, "time");
    }

    static IObject* defaultTime(){
        return new DateTime::Time(0);
    }

    /** Invoked when DateTime class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDRIProxyType(DateTimeProxy::TYPE); // use proxy for dri
        REGISTER(DateTime, clazz);
        SUPERCLASS(CObject);
        clazz->enableCloneOptimisations(); /* avoid cloning DateTimes when
                                              using copy constructor */
        EMPTY_SHELL_METHOD(defaultDateTime);
        FIELD(date, "date");
        FIELD(time, "time");

        // register our addin functions
        Addin::registerClassObjectMethod(
            "DATETIME",
            Addin::UTILITIES,
            "Constructs a handle to a date-time",
            DateTimeAddin::TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createDateTime);

        Addin::registerClassObjectMethod(
            "DATETIME_ARRAY",
            Addin::UTILITIES,
            "Constructs a handle to an array of dateTimes",
            DateTimeArrayAddin::TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createDateTimeArray);
    }

    static IObject* defaultDateTime(){
        return new DateTime();
    }

    /** addin function - either creates a DateTime or a CNull if no datetime is
        supplied */
    static IObjectSP createDateTime(DateTimeAddin* params){
        if (!params->dateTime){
            return CNull::create();
        }
        return params->dateTime;
    }

    /** addin function - returns a DateTimeArray (creates empty one if
        no dates supplied) */
    static IObjectSP createDateTimeArray(DateTimeArrayAddin* params){
        if (!params->dateTimeArray){
            return IObjectSP(new DateTimeArray());
        }
        return params->dateTimeArray;
    }
};

CClassConstSP const DateTime::Date::TYPE = CClass::registerClassLoadMethod(
    "DateTime::Date", typeid(DateTime::Date), DateTimeHelper::loadDate);
CClassConstSP const DateTime::Time::TYPE = CClass::registerClassLoadMethod(
    "DateTime::Time", typeid(DateTime::Time), DateTimeHelper::loadTime);
CClassConstSP const DateTime::TYPE = CClass::registerClassLoadMethod(
    "DateTime", typeid(DateTime), DateTimeHelper::load);

DRLIB_END_NAMESPACE
