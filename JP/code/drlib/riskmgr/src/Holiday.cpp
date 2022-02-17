//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Holiday.cpp
//
//   Description : Holiday representation
//
//   Author      : Andrew J Swain 
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------
 
#include "edginc/config.hpp"
#define QLIB_HOLIDAY_CPP
#include "edginc/Holiday.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include <algorithm>
#include <set>
#include <fstream>

DRLIB_BEGIN_NAMESPACE

Holiday::~Holiday(){}

/** Holidays are essentially immutable objects so clone just returns this
    - more or less */
IObject* Holiday::clone() const{
    Holiday* hol;
    if (getRefCount() == 0){
        // not being accessed via smart pointer - object probably on stack
        hol = new Holiday(*this);
    } else {
        hol = const_cast<Holiday*>(this);
    }
    return hol;
}

/** Returns true if the suppplied hols are same as this hol */
bool Holiday::equals(const Holiday* hols) const{
    if (hols == this ||
        (hols->name == name && hols->useWeekends == useWeekends &&
         DateTime::equals(*hols->holidays, *holidays))){
        return true;
    }
    return false;
}

/** Returns a simple hash code value: the hash code of the holiday name */
int Holiday::hashCode() const {
    return hash_string(name);
}


//// constructor that just copies references over
Holiday::Holiday(const Holiday &rhs): 
    MarketObject(TYPE), name(rhs.name), useWeekends(rhs.useWeekends),
    holidays(rhs.holidays), cache(rhs.cache){}

/** Constructs a Holiday with full flexibility
    @param name        identifier for printing
    @param holidays    non-business days (may be null)
    @param useWeekends weekends (Saturday & Sunday) are not business days
*/
Holiday::Holiday(const string&        name, 
                 const DateTimeArray& holidays,
                 bool                 useWeekends) :
    MarketObject(TYPE),
    name(name), useWeekends(useWeekends), 
    holidays(new DateTimeArray(holidays)) {
    validatePop2Object();
}

void Holiday::validatePop2Object() {
    // check holidays are in ascending order
    // and unique by date. The time of day is irrelevant as we always
    // ignore it. Also remove any hols on weekends if useWeekends is true
    // validatePop2Object has privileged access
    DateTimeArray& holidays = const_cast<DateTimeArray&>(*this->holidays);
    for (vector<DateTime>::iterator iter = holidays.begin(); 
         iter != holidays.end(); /* inc in loop body */){
        if (iter != holidays.begin() && 
            iter->getDate() <= (iter-1)->getDate()){
            throw ModelException("Holiday::validatePop2Object",
                                 "holidays must be strictly increasing "
                                 "and on distinct dates: holiday (" +
                                 iter->toString() + 
                                 ") is on or before (" + 
                                 (iter-1)->toString() + ")");
        }
        if (useWeekends && iter->isWeekend()){
            iter = holidays.erase(iter);
        } else {
            ++iter;
        }
    }
}

/** Create a Holiday where only weekends are non-business days */
Holiday* Holiday::weekendsOnly() {
    DateTimeArray empty;
    return new Holiday("weekends only", empty, true);
}
    
/** Create a Holiday where all days are business days */
Holiday* Holiday::noHolidays() {
    DateTimeArray empty;
    return new Holiday("no holidays", empty, false);
}

/** is a date a business day? */
bool Holiday::isBusinessDay(const DateTime& date) const {
    return !isHoliday(date);
}

struct ComparatorLT{
    bool operator()(const DateTime& d1, const DateTime& d2){
        return (d1.getDate() < d2.getDate());
    }
};
struct ComparatorGT{
    bool operator()(const DateTime& d1, const DateTime& d2){
        return (d1.getDate() > d2.getDate());
    }
};

/** used for translating between EDR and ALIB holiday objects only.
    Simply returns the array of holiday dates and whether
    weekends are included as bad days or not */
DateTimeArrayConstSP Holiday::toALIB(bool& areWeekendsHolidays) const {
    areWeekendsHolidays = useWeekends;
    return holidays;
}

/** is a date a holiday (includes weekends) */
bool Holiday::isHoliday(const DateTime& date) const {
    if (useWeekends && date.isWeekend()){
        return true;
    }
    return binary_search(holidays->begin(), holidays->end(), 
                         date, ComparatorLT());
}

/** add a number of business days to a date */
DateTime Holiday::addBusinessDays(const DateTime& from, int busDays) const 
{
    if (busDays == 0){
        return from; // optimise return
    }
    if (holidays->empty() && !useWeekends) {
        // No adjustments at all.
        return from.rollDate(busDays);
    } 
    int      intervalSign   = Maths::sign(busDays);
    int      numBusDaysLeft = abs(busDays);
    if (intervalSign == 1 && useWeekends) {
        return forwardStandardWeekends(from, numBusDaysLeft);
    } 
    // Get the number of business days per week. In-line for speed.
    int  busDaysPerWeek = useWeekends? 5:  DateTime::DAYS_PER_WEEK;
    return forwardNonStandardWeekends(from,
                                      numBusDaysLeft, 
                                      intervalSign,
                                      busDaysPerWeek);
}  

/** return copy of holiday name */
string Holiday::getName() const{
    return name;
}

/** is holiday equivalent to ALIB "NO_WEEKENDS" ? */
bool Holiday::isAllWorkingDays() const{
    return (!useWeekends && holidays->size() == 0);
}  

/** Calculates the number of business days between two dates (FROM & TO).
    
   Algorithm:
   1. if FROM = TO, the result is 0 
   2. if FROM < TO, the result is the number of business days in the CLOSED
   interval of [FROM+1,TO]
   3. if FROM > TO, the result is negated number of business days in the 
      CLOSED interval of [TO,FROM-1] */
int Holiday::businessDaysDiff(const DateTime&  fromDate,
                              const DateTime&  toDate) const{
    int fromDt = fromDate.getDate();
    int toDt = toDate.getDate();
    // check whether both dates are the same
    if (fromDt == toDt){
        return 0;
    }
    /* Set the direction. */ 
    int         signum     = toDt < fromDt? -1: 1;

    /* First, compute the date difference adjusting for weekends only. */ 

    /* Get the number of business days per week. */
    int busDaysPerWeek = useWeekends? 5: DateTime::DAYS_PER_WEEK;

    int numWeeks = abs((toDt - fromDt) / DateTime::DAYS_PER_WEEK);
    int  curDate = fromDt + DateTime::DAYS_PER_WEEK * numWeeks * signum;

    // count number of extra business days
    int nrExtraDays = 0;
    if (!useWeekends){
        nrExtraDays = abs(curDate - toDt);
    } else {
        /* need to determine how many bus days there are between 
           (curDate, toDate] */
        int numDaysToCheck = abs(curDate - toDt);
        for (int i = 0; i < numDaysToCheck; i++){
            curDate += signum;
            if (!DateTime::isWeekend(curDate)){
                nrExtraDays++;
            } 
        }
    }
  
    int result = (numWeeks * busDaysPerWeek) + nrExtraDays; 

    /* Now count the number of weekday holidays and adjust the date
    ** difference by the result. The call to calcNumWeekdayHolidays
    ** searches the holiday list at most once. */
    if (holidays->size() > 0) {
        DateTime  nextDay = DateTime(fromDt + signum,
                                     fromDate.getTime());

        if (signum == 1){
            HolIter pos(lower_bound(holidays->begin(), holidays->end(), 
                                    nextDay, ComparatorLT()));
            if (pos != holidays->end()){
                result  -= calcNumWeekdayHolidays(toDate, pos, 
                                                  holidays->end());
            }
        } else {
            HolRevIter pos(lower_bound(holidays->rbegin(), holidays->rend(), 
                                       nextDay, ComparatorGT()));
            if (pos != holidays->rend()){
                result  -= calcNumWeekdayHolidays(toDate, pos, 
                                                  holidays->rend());
            }
        }
    }
    result *= signum;      /* change sign if we are going backwards. */
    return result;
}

/** Compute the number of week-day holidays starting at the input holiday 
    index and continuing until the holiday list is at or beyond toDate. 
    Week-day holidays are holidays which occur on a non-weekend day. */
int Holiday::calcNumWeekdayHolidays (
    const DateTime&   toDate,        /* End date.                */
    HolIter&          startPt,       /* Starting holiday index.  */ 
    const HolIter&    endPt) const { /* Starting holiday index.  */ 
    // make assumption that if weekends are hols all dates which are
    // weekends have been removed
    HolIter pos = upper_bound(startPt, endPt, toDate, ComparatorLT());
    int numHols = pos-startPt;
    startPt = pos;
    return numHols;
}

int  Holiday::calcNumWeekdayHolidays (
    const DateTime&   toDate,        /* End date.                */
    HolRevIter&       startPt,       /* Starting holiday index.  */ 
    const HolRevIter& endPt) const { /* Starting holiday index.  */ 
    // make assumption that if weekends are hols all dates which are
    // weekends have been removed
    HolRevIter pos = upper_bound(startPt, endPt, toDate, ComparatorGT());
    int numHols = pos-startPt;
    startPt = pos;
    return numHols;
}

/** Same as findNextHol but with only weekends as possible holidays */
void Holiday::findNextWeekEndHol(const DateTime&   startDate, 
                                 bool              fwdSearch, 
                                 DateTime&         holStart,
                                 DateTime&         holEnd,  
                                 bool&             foundHol) const{
    if (!useWeekends){
        foundHol = false;
    } else {
        foundHol = true;
        int weekDay = startDate.getWeekday();
        if (weekDay < 2){
            holStart = startDate; // on a weekend already
            holEnd = startDate.rollDate(fwdSearch?
                                        1-weekDay: -weekDay); // Sa/Su
        } else {
            // so move to next Sat (if fwds)
            holStart = startDate.rollDate(fwdSearch? 7-weekDay: 1-weekDay);
            holEnd = holStart.rollDate(fwdSearch? 1: -1);
        }
    }
}

/* returns two dates defining the next holiday on or after/before
   start date.  When fwdSearch is TRUE the two dates returned are
   holStart which is MAX(startDate, start of holiday) and holEnd which
   is the date of the last day in the holiday. When
   fwdSearch is FALSE then the algorithm works as if you travelling backwards
   in time so that holStart will be MIN(startDate, last day in holiday) etc.
   NB The time of day in the dates returned is undefined */
void Holiday::findNextHol(const DateTime&   startDate, 
                          bool              fwdSearch, 
                          DateTime&         holStart,
                          DateTime&         holEnd,  
                          bool&             foundHol) const
{
    /* Strategy: find next hol index. If useWeekends, roll from start
       date until hit the hol or a weekend. Roll across hol/weekend , then
       roll until next hol stopping when date is not a weekend. Repeat
       rolling until next hol etc */
    if (fwdSearch){
        // first find next hol index
        HolIter pos = lower_bound(holidays->begin(), holidays->end(), 
                                  startDate, ComparatorLT());
        if (pos == holidays->end()){
            // gone beyond end of hols
            findNextWeekEndHol(startDate, fwdSearch, holStart,
                               holEnd, foundHol);
            return;
        }
        foundHol = true;
        const HolIter& stopPt = holidays->end();
        if (useWeekends){
            // roll from start date until hit the hol or a weekend. 
            int date = startDate.getDate();
            while (date < pos->getDate() && !DateTime::isWeekend(date)){
                date++;
            }
            holStart = DateTime(date, DateTime::START_OF_DAY_TIME); // save
            //Roll across hol/weekend
            do{
                if (pos != stopPt && date == pos->getDate()){
                    // if date points to current holiday, move to next holiday
                    ++pos;
                }
                date++;
            } while (DateTime::isWeekend(date) ||
                     (pos != stopPt && date == pos->getDate()));
            holEnd = DateTime(date-1, DateTime::END_OF_DAY_TIME);
        } else {
            holStart = *pos; // save
            do{
                ++pos;
            } while (pos != stopPt && 
                     (pos-1)->getDate()+1 == pos->getDate());
            holEnd = *(pos-1);
        }
    } else {
        // first find next hol index
        HolRevIter pos = lower_bound(holidays->rbegin(), holidays->rend(), 
                                     startDate, ComparatorGT());
        if (pos == holidays->rend()){
            // gone beyond end of hols
            findNextWeekEndHol(startDate, fwdSearch, holStart,
                               holEnd, foundHol);
            return;
        }
        foundHol = true;
        const HolRevIter& stopPt = holidays->rend();
        if (useWeekends){
            // roll from start date until hit the hol or a weekend. 
            int date = startDate.getDate();
            while (date > pos->getDate() && !DateTime::isWeekend(date)){
                date--;
            }
            holStart = DateTime(date, DateTime::END_OF_DAY_TIME); // save
            //Roll across hol/weekend
            do{
                if (pos != stopPt && date == pos->getDate()){
                    // if date points to current holiday, move to next holiday
                    ++pos;
                }
                date--;
            } while (DateTime::isWeekend(date) ||
                     (pos != stopPt && date == pos->getDate()));
            holEnd = DateTime(date+1, DateTime::START_OF_DAY_TIME);
        } else {
            holStart = *pos; // save
            do{
                ++pos;
            } while (pos != stopPt && 
                     (pos-1)->getDate()-1 == pos->getDate());
            holEnd = *(pos-1);
        }
    }
}


/** Calculates the number of business days for the given year and whether 
    or not it is a leap year  */
void Holiday::numBusDaysInYear(int          year, 
                               int&         numBusDays, 
                               bool&        isLeapYear) const
{
    // Check if information is already in the cache 
    intMap::iterator iter(cache.find(year));
    if (iter == cache.end()) {
        // convert start and end dates of given year to TDate format
        DateTime::MonthDayYear        startOfYearMDY(31,12,year-1);
        DateTime::MonthDayYear        endOfYearMDY(31,12,year);
        DateTime startOfYear(startOfYearMDY.toDateTime());
        DateTime endOfYear(endOfYearMDY.toDateTime());

        numBusDays = businessDaysDiff(startOfYear, endOfYear);

        // insert the new value into the hash table
        cache.insert(intMap::value_type(year, numBusDays));
    }
    else  /* data is already in the cache */
    {
        // second is what normal people would call value in (key,value) pairs
        numBusDays = iter->second;
    }
    
    isLeapYear = DateTime::isLeapYear(year);
}


/** This function forwards a date by given number of business days under 
    the following conditions: 
    (1) Saturdays and Sundays are both holidays.
    (2) The direction is forwards in time.                            */ 
DateTime Holiday::forwardStandardWeekends ( 
    const DateTime& fromDate,       /* (I) start date */ 
    int             numBusDays) const
{
    if (numBusDays == 0) {
        return fromDate;
    }

    /*
     * Table containing the offset to add to
     * a date to skip a given number of business
     * days.  Only valid for less than 5 business
     * days, and when saturday and sunday are both
     * holidays.
     *
     * So:
     * offsetTable[dayOfWeek(date)][nDays]
     * is the number of days to add to date
     * to get a new business day.
     *
     * Note: The indexes are set up as DateTime mod 7
     */
    static int offsetTable[7][5] = 
    {/*n: 0  1  2  3  4 */
        {     0, 1, 2, 3, 4 },  /* Monday */
        {     0, 1, 2, 3, 6 },  /* Tuesday */
        {     0, 1, 2, 5, 6 },  /* Wednesday */
        {     0, 1, 4, 5, 6 },  /* Thursday */
        {     0, 3, 4, 5, 6 },  /* Friday */
        {    -1, 2, 3, 4, 5 },  /* Saturday */
        {    -2, 1, 2, 3, 4 }   /* Sunday */
    };


    /*
     * Calculate result if no holidays are involved.
     * 
     * Move forward over a week for every 5 business days.  
     * Use a table for moving 0..4 business days from each day of the week.
     */

    int fromDt = fromDate.getDate();
    int fromTm = fromDate.getTime();
    int dt  = fromDt +  7 * (numBusDays / 5) + 
                       offsetTable[ fromDt % 7 ][ numBusDays % 5 ];

    int numHols  = holidays->size();
    if (numHols > 0) {
        /* Handle any holidays there are. Use binary search() to find date in 
         * the datelist. */
        HolIter pos(upper_bound(holidays->begin(), holidays->end(), 
                                fromDate, ComparatorLT()));

        while (pos != holidays->end() &&  dt >= pos->getDate()) {
            dt += offsetTable[dt % 7 ][ 1 ];
            ++pos;        
        } 
    }
    return DateTime(dt, fromTm);
}

/** This function handles the special case where the number of business 
    days to go forward by is very large. */ 
DateTime Holiday::forwardNonStandardWeekends (
    const DateTime& fromDate,             /* (I) start date */ 
    int             numBusDaysLeft,       /* (I) abs. num. bus. days */ 
    int             direction,            /* (I) +1=forwards, -1=backwards */
    int             busDaysPerWeek) const /* (I) num. bus. days per week */ 
{ 
    /* First, adjust for weekends only, pretending there are no holidays. */ 
    int numWeeks = Maths::max(numBusDaysLeft/busDaysPerWeek, 0);
    DateTime curDate(fromDate.rollDate(
                         DateTime::DAYS_PER_WEEK*numWeeks*direction));
    /* Search the holiday list for the first holiday strictly after
    ** (if going forward in time) or strictly before (if going
    ** backward in time) the start date. Note that the holiday list is
    ** assumed to be sorted in increasing order. */  
    if (direction == 1){
        const HolIter& stopPt = holidays->end();
        HolIter pos = upper_bound(holidays->begin(), stopPt, fromDate,
                                  ComparatorLT());

        /* Count the number of weekday holidays starting from the current
        ** holiday index and going up to the curDate.  Count weekday
        ** holidays only because holidays occurring on week-end days have
        ** been skipped by the previous calculation. */ 
        numBusDaysLeft -= (busDaysPerWeek * numWeeks);
        if (pos != stopPt){
            numBusDaysLeft += calcNumWeekdayHolidays(curDate, pos, stopPt);
        }
        /* Now search one day at a time starting one day beyond
           the current date. */
        int date = curDate.getDate();
        while (numBusDaysLeft > 0) { 
            date++;
            /* Check if the day is a holiday first. If it is, don't decrement 
               numBusDaysLeft. */ 
            if (pos != stopPt  && 
                date == pos->getDate()){
                ++pos;
            } else {
                /* Not a holiday. */
                /* If the day is a weekday, then decrement numBusDaysLeft, 
                ** otherwise continue looping. */ 
                if (!(useWeekends && DateTime::isWeekend(date))) {
                    numBusDaysLeft--; 
                }
            } 
        }
        curDate = DateTime(date, curDate.getTime());
    } else { 
        HolRevIter stopPt = holidays->rend();
        HolRevIter pos = upper_bound(holidays->rbegin(), stopPt,
                                     fromDate, ComparatorGT());

        /* Count the number of weekday holidays starting from the current
        ** holiday index and going up to the curDate.  Count weekday
        ** holidays only because holidays occurring on week-end days have
        ** been skipped by the previous calculation. */ 
        numBusDaysLeft -= busDaysPerWeek * numWeeks;
        if (pos != stopPt){
            numBusDaysLeft += calcNumWeekdayHolidays(curDate, pos, stopPt);
        }
    
        /* Now search one day at a time starting one day beyond
        ** the current date. */    
        int date = curDate.getDate();
        while (numBusDaysLeft > 0) { 
            date--;
            /* Check if the day is a holiday first. If it is, don't decrement 
               numBusDaysLeft. */ 
            if (pos != stopPt  && 
                date == pos->getDate()){
                ++pos;
            } else {
                /* Not a holiday. */
                /* If the day is a weekday, then decrement numBusDaysLeft, 
                ** otherwise continue looping. */ 
                if (!(useWeekends && DateTime::isWeekend(date))) {
                    numBusDaysLeft--; 
                }
            } 
        } 
        curDate = DateTime(date, curDate.getTime());
    }
    return curDate; 
}

void Holiday::acceptWrapperNameCollector(Holiday* holiday, 
                                         WrapperNameCollector* collector)
{
    collector->addName(holiday->getName());
}


/** for reflection */
Holiday::Holiday():MarketObject(TYPE), useWeekends(true){}

class HolidayHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Holiday, clazz);
        SUPERCLASS(MarketObject);
        EMPTY_SHELL_METHOD(defaultHoliday);
        FIELD(name, "Name identifying holidays");
        FIELD(useWeekends, "Whether weekends are holidays");
        FIELD(holidays, "holidays");
        ClassSetAcceptMethod(Holiday::acceptWrapperNameCollector);

        Addin::registerConstructor("HOLIDAY",
                                   Addin::UTILITIES,
                                   "Creates a Holiday object",
                                   Holiday::TYPE);
    }
    
    static IObject* defaultHoliday(){
        return new Holiday();
    }
};

CClassConstSP const Holiday::TYPE = CClass::registerClassLoadMethod(
    "Holiday", typeid(Holiday), HolidayHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(HolidayWrapper);

class IsHolidayAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    // addin parameters
    HolidaySP         hols;
    DateTimeArraySP   dates;

    static IObjectSP isHoliday(IsHolidayAddin* params){
        CBoolArraySP results(new BoolArray(params->dates->size()));
        for (int i = 0; i < results->size(); i++){
            (*results)[i] = params->hols->isHoliday((*params->dates)[i]);
        }
        return results;
    }

    IsHolidayAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IsHolidayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIsHolidayAddin);
        FIELD(hols, "Holiday object");
        FIELD(dates, "dates to query");
        Addin::registerInstanceObjectMethod("IS_HOLIDAY",
                                            Addin::MARKET,
                                            "Determines if supplied dates are"
                                            " holidays",
                                            TYPE,
                                            false,
                                            Addin::expandSimple,
                                            (Addin::ObjMethod*)isHoliday);
    }
    
    static IObject* defaultIsHolidayAddin(){
        return new IsHolidayAddin();
    }
};
CClassConstSP const IsHolidayAddin::TYPE = CClass::registerClassLoadMethod(
    "IsHolidayAddin", typeid(IsHolidayAddin), load);

class BusDaysDiffAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    // addin parameters
    HolidaySP         hols;
    DateTimeArraySP   dates1;
    DateTimeArraySP   dates2;

    static IObjectSP busDaysDiff(BusDaysDiffAddin* params){
        if (params->dates2->size() != params->dates1->size()){
            throw ModelException("BusDaysDiffAddin", "Both arrays of dates"
                                 " must be the same length");
        }
        CIntArraySP results(new IntArray(params->dates1->size()));
        for (int i = 0; i < results->size(); i++){
            (*results)[i] = 
                params->hols->businessDaysDiff((*params->dates1)[i],
                                               (*params->dates2)[i]);
        }
        return results;
    }

    BusDaysDiffAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BusDaysDiffAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBusDaysDiffAddin);
        FIELD(hols, "Holiday object");
        FIELD(dates1, "set of first dates");
        FIELD(dates2, "set of second dates");
        Addin::registerInstanceObjectMethod("BUS_DAYS_DIFF",
                                            Addin::MARKET,
                                            "Determines number of business "
                                            " days between pairs of 2 dates",
                                            TYPE,
                                            false,
                                            Addin::expandSimple,
                                            (Addin::ObjMethod*)busDaysDiff);
    }
    
    static IObject* defaultBusDaysDiffAddin(){
        return new BusDaysDiffAddin();
    }
};
CClassConstSP const BusDaysDiffAddin::TYPE = CClass::registerClassLoadMethod(
    "BusDaysDiffAddin", typeid(BusDaysDiffAddin), load);

class BusDaysFwdAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    // addin parameters
    HolidaySP         hols;
    DateTimeArraySP   dates;
    CIntArraySP       numBusDays;

    static IObjectSP busDaysFwd(BusDaysFwdAddin* params){
        if (params->dates->size() != params->numBusDays->size()){
            throw ModelException("BusDaysFwdAddin", "Both arrays"
                                 " must be the same length");
        }
        DateTimeArraySP results(new DateTimeArray(params->dates->size()));
        for (int i = 0; i < results->size(); i++){
            (*results)[i] = 
                params->hols->addBusinessDays((*params->dates)[i], 
                                              (*params->numBusDays)[i]);
        }
        return results;
    }

    BusDaysFwdAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BusDaysFwdAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBusDaysFwdAddin);
        FIELD(hols, "Holiday object");
        FIELD(dates, "dates");
        FIELD(numBusDays, "number of business days to add");
        Addin::registerInstanceObjectMethod("BUS_DAYS_FWD",
                                            Addin::MARKET,
                                            "Calculates new dates using offset"
                                            " of given number of business "
                                            "dates",
                                            TYPE,
                                            false,
                                            Addin::expandSimple,
                                            (Addin::ObjMethod*)busDaysFwd);
    }
    
    static IObject* defaultBusDaysFwdAddin(){
        return new BusDaysFwdAddin();
    }
};
CClassConstSP const BusDaysFwdAddin::TYPE = CClass::registerClassLoadMethod(
    "BusDaysFwdAddin", typeid(BusDaysFwdAddin), load);

class FindNextHolAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    // addin parameters
    HolidaySP         hols;
    DateTimeArraySP   dates;
    CBoolArraySP      fwdSearch;

    static IObjectSP findNextHol(FindNextHolAddin* params){
        if (params->dates->size() != params->fwdSearch->size()){
            throw ModelException("FindNextHolAddin", "Both arrays"
                                 " must be the same length");
        }
        ObjectArraySP results(new ObjectArray(params->dates->size()));
        for (int i = 0; i < results->size(); i++){
            bool     foundHol;
            DateTimeArraySP boundary(new DateTimeArray(2));
            params->hols->findNextHol((*params->dates)[i],
                                      (*params->fwdSearch)[i],
                                      (*boundary)[0], (*boundary)[1],
                                      foundHol);
            if (foundHol){
                (*results)[i] = boundary;
            }
        }
        return results;
    }

    FindNextHolAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FindNextHolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFindNextHolAddin);
        FIELD(hols, "Holiday object");
        FIELD(dates, "date");
        FIELD(fwdSearch, "Forwards or backwards search");
        Addin::registerInstanceObjectMethod("FIND_NEXT_HOL",
                                            Addin::MARKET,
                                            "Finds next holiday period "
                                            "before/after given date",
                                            TYPE,
                                            false,
                                            Addin::expandMulti,
                                            (Addin::ObjMethod*)findNextHol);
    }
    
    static IObject* defaultFindNextHolAddin(){
        return new FindNextHolAddin();
    }
};
CClassConstSP const FindNextHolAddin::TYPE = CClass::registerClassLoadMethod(
    "FindNextHolAddin", typeid(FindNextHolAddin), load);


class LoadHolidaysAddin: public CObject
{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    string calendar;
    string holidays;
    string additionalHolidays;
    string workingDays;
    string format;

private:
    /*
     * Read holiday data file in ALIB format,
     * ie. YYYYMMDD comment text
     *
     * The file may have comment lines deliminated by a leading '#', and also 
     * weekend day directives [the '#' must be in the first column].
     *
     * cf. ALIB: buscache.c#1068 (GtoHolidayListRead)
     */
    static bool readAlibFile(const string& calendar, ifstream& is, set<DateTime>& dates)
    {
        static const string method = "LoadHolidaysAddin::readAlibFile";

        // read holiday data into supplied set
        char buffer[256];
        bool saturday = true;
        bool sunday = true;
        string NOT_SATURDAY = "# SATURDAY_NOT_ALWAYS_HOLIDAY";
        string NOT_SUNDAY = "# SUNDAY_NOT_ALWAYS_HOLIDAY";

        while (is.getline(buffer, 255))
        {
            if (buffer[0] == '#')
            {
                if (CString::equalsIgnoreCase(buffer, NOT_SATURDAY, NOT_SATURDAY.length()))
                {
                    saturday = false;
                }

                if (CString::equalsIgnoreCase(buffer, NOT_SUNDAY, NOT_SUNDAY.length()))
                {
                    sunday = false;
                }
            }
            else
            {
                long date = atol(buffer);
                if (date > 16010101)
                {
                    dates.insert(DateTime::fromIrDate(date));
                }
            }
        }

        return saturday && sunday;
    }


    /*
     * Convert string in MM/DD/YYYY format to a (long) date in YYYYMMDD format.
     */
    static int stringToDate(const char* str)
    {
        /* pull numeric values out of string */
        int month = atoi(&str[0]);
        int day = atoi(&str[3]);
        int year = atoi(&str[6]);

        /* convert days, months and years to date */
        int date = (long)year*10000 + (long)month*100 + (long)day;
        return date;
    }


    /*
     * Read holiday data file in ChaseTools / Catalyst format,
     * ie. MM/DD/YYYY,NOM (NOM = holiday calendar name).
     *
     * The file may have comments after the calendar name, and are deliminated
     * by whitespace or '#'.
     *
     * cf. ALIB: busday.c#2218 (readHolidaysAll)
     */
    static bool readChaseToolsFile(const string& calendar, const string& filename, ifstream& is, set<DateTime>& dates)
    {
        static const string method = "LoadHolidaysAddin::readChaseToolsFile";

        // read holiday data into supplied set
        int i;
        char ch;
        const int DATE_BUFFER_LEN = 20;
        const int NAME_BUFFER_LEN = 20;
        char dateBuffer[DATE_BUFFER_LEN];
        char nameBuffer[NAME_BUFFER_LEN];

        do
        {
            // read date
            dateBuffer[0] = '\0';
            for (i = 0 ; i < DATE_BUFFER_LEN-1 ; i++)
            {
                ch = is.get();
                if (ch == EOF || ch == ' ' || ch == ',' || ch == '\n' || ch == '\t')
                    break;
                dateBuffer[i] = ch;
            }

            if (i != 0)
            {
                dateBuffer[i] = '\0';
            }

            if (i == DATE_BUFFER_LEN)
            {
                string msg = Format::toString(
                    "Overflowed buffer while reading date from file %s",
                    filename.c_str());
                throw ModelException(method, msg);
            }

            // read calendar name
            nameBuffer[0] = '\0';
            i = 0;
            bool endOfCentreData = false;
            do
            {
                ch = toupper(is.get());
                if (ch == ' ' || ch == '#')
                {
                    endOfCentreData = true;
                }

                if (ch == EOF || ch == '\n')
                {
                    break;
                }

                if (!endOfCentreData)
                {
                    nameBuffer[i] = ch;
                    i++;
                }
            } while(ch != '\n' && i < NAME_BUFFER_LEN - 1);

            if (i != 0)
            {
                nameBuffer[i] = '\0';
            }

            // check for errors
            if (ch != EOF && i == 0 && dateBuffer[0] != '\0')
            {
                string msg = Format::toString(
                    "No holiday centre specified for %s in file %s", 
                    dateBuffer, filename.c_str());
                throw ModelException(method, msg);
            }

            if (i == NAME_BUFFER_LEN)
            {
                string msg = Format::toString(
                    "Overflowed buffer while reading holiday from file %s", 
                    filename.c_str());
                throw ModelException(method, msg);
            }

            // if data is for the calendar wanted add date to output set
            // (ignores case when comparing names)
            if (CString::equalsIgnoreCase(nameBuffer, calendar))
            {
                int dt = stringToDate(dateBuffer);
                dates.insert(DateTime::fromIrDate(dt));
            }
        } while(ch != EOF);

        return true; // data format assumes weekends are non-working days
    }

    /*
     * Read holiday data file.
     */
    static bool read(const string& format, const string& calendar, const string& filename, set<DateTime>& dates)
    {
        static const string method = "LoadHolidaysAddin::read";

        bool useWeekends = true; 

        // open data file
        ifstream is(filename.c_str());
        if (!is)
        {
            string msg = Format::toString(
                "Could not open %s holiday data file", filename.c_str());
            throw ModelException(method, msg);
        }

        try
        {
            if (CString::equalsIgnoreCase(format,"ALIB"))
            {
                useWeekends = readAlibFile(calendar, is, dates);
            }
            else if (CString::equalsIgnoreCase(format,"ChaseTools"))
            {
                useWeekends = readChaseToolsFile(calendar, filename, is, dates);
            }
            else
            {
                string msg = Format::toString("Unrecognized data format (%s)", format.c_str());
                throw ModelException(method, msg);
            }

            is.close();
        }
        catch(const ModelException& x)
        {
            is.close();
            throw x;
        }
        catch(...)
        {
            is.close();

            string msg = Format::toString(
                "Could not read holiday data from %s for %s", 
                filename.c_str(), calendar.c_str());
            throw ModelException(method, msg);
        }

        return useWeekends;
    }

public:
    static IObjectSP loadHolidays(LoadHolidaysAddin* params)
    {
        static const string method = "LoadHolidaysAddin::loadHolidays";

        if (params->calendar.empty())
        {
            throw ModelException(method, "holiday calendar name not provided");
        }

        if (params->holidays.empty())
        {
            throw ModelException(method, "holiday data file name not provided");
        }

        DateTimeArray holidays;
        set<DateTime> dates;
        bool useWeekends = read(params->format, params->calendar, params->holidays, dates);

        if (!params->additionalHolidays.empty())
        {
            read(params->format, params->calendar, params->additionalHolidays, dates);
        }

        if (!params->workingDays.empty())
        {
            set<DateTime> workingDays;
            read(params->format, params->calendar, params->workingDays, workingDays);
            set_difference(dates.begin(), dates.end(), workingDays.begin(), workingDays.end(), holidays.back_inserter());
        }
        else
        {
            std::copy(dates.begin(), dates.end(), holidays.back_inserter());
        }

        std::sort(holidays.begin(), holidays.end());
        return HolidaySP(new Holiday(params->calendar, holidays, useWeekends));
    }

    LoadHolidaysAddin(): CObject(TYPE), format("ChaseTools") {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LoadHolidaysAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLoadHolidaysAddin);
        FIELD(calendar,           "Name of holiday calendar to load");
        FIELD(holidays,           "Name of file containing holiday data");
        FIELD(additionalHolidays, "Name of file containing additional holidays");
        FIELD(workingDays,        "Name of file contining days that are not holidays");
        FIELD(format,             "Data file format - ALIB or ChaseTools (default)");
        FIELD_MAKE_OPTIONAL(additionalHolidays);
        FIELD_MAKE_OPTIONAL(workingDays);
        FIELD_MAKE_OPTIONAL(format);

        Addin::registerInstanceObjectMethod("HOLIDAY_LOAD",
                                            Addin::UTILITIES,
                                            "Creates Holiday object from file data",
                                            TYPE,
                                            false,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*) loadHolidays);
    }
    
    static IObject* defaultLoadHolidaysAddin(){
        return new LoadHolidaysAddin();
    }
};

CClassConstSP const LoadHolidaysAddin::TYPE = CClass::registerClassLoadMethod(
    "LoadHolidaysAddin", typeid(LoadHolidaysAddin), load);

bool HolidayLoad() {
    return Holiday::TYPE != NULL &&
           BusDaysDiffAddin::TYPE != NULL &&
           BusDaysFwdAddin::TYPE != NULL &&
           FindNextHolAddin::TYPE != NULL &&
           Holiday::TYPE != NULL &&
           HolidayWrapper::TYPE != NULL &&
           IsHolidayAddin::TYPE != NULL &&
           LoadHolidaysAddin::TYPE != NULL;
}

DRLIB_END_NAMESPACE
