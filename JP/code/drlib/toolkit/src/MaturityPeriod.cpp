//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MaturityPeriod.cpp
//
//   Description : Defines floating expiries used to define yield curve & 
//                 vol surface points e.g. 1M, 5Y
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MATURITYPERIOD_CPP
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Hashtable.hpp"

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<MaturityPeriod>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<MaturityPeriod>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<MaturityPeriodSP _COMMA_ MaturityPeriod>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<MaturityPeriodArray>);
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<MaturityPeriodSP>(
                         MaturityPeriodSP* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<MaturityPeriodSP>(
                         MaturityPeriodSP* t, IObjectSP o));

#define MAX_BUFFER 32
// could be a static field in the class but then I have to check the 
// header out ....
#define VALID_MATURITIES "DMAYSQWIJKCT"
/** Override clone method to copy our extra data over */
IObject* MaturityPeriod::clone() const{
    int& count = getRefCount();
    if (count == 0){
        return new MaturityPeriod(*this);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

// implemented for fast clone
MaturityPeriod::MaturityPeriod(const MaturityPeriod& matPeriod):
    Expiry(TYPE), period(matPeriod.period), 
    interval(matPeriod.interval), count(matPeriod.count){}

/** Returns true if given expiry matches this */
bool MaturityPeriod::equalTo(const IObject* expiry) const{
    if (this == expiry){
        return true;
    }
    if (!expiry || expiry->getClass() != TYPE){
        return false;
    }
    const MaturityPeriod* matPeriod = STATIC_CAST(MaturityPeriod, expiry);
    return (count == matPeriod->count && interval == matPeriod->interval);
}

/** Returns true if given expiry matches this */
bool MaturityPeriod::equals(const Expiry* expiry) const{
    return equalTo(expiry);
}

/** returns a hashcode for the object based upon the period.
    Overridden for performance */
int MaturityPeriod::hashCode() const{
    return (((size_t) TYPE) ^ hash_string(period));
}

/** this gets called after an object is constructed from a data dictionary.
    Not after an object has been copied (see override of clone method below) */
void MaturityPeriod::validatePop2Object(){
    static const string routine = "MaturityPeriod::validatePop2Object";
    // special case for ON
    if (period == "ON") {
        // validatePop2Object has privileged access
        const_cast<int&>(count)    = 1;
        const_cast<string&>(interval) = "D";
    } 
    else {           
        // chop string (e.g. 6M) into count (i.e. 6) and interval (i.e. M)

        const char *inp = period.c_str();
        char  numberBuff[MAX_BUFFER];
        char *nump = numberBuff;            /* Pointer to number */

        /* Copy sign,if any, to numberBuff */
        if (*inp == '-' || *inp == '+')
        {
            *nump++ = *inp++;
        }

        /* Copy digits, if any, to numberBuff */
        while (isdigit(*inp)) { 
            *nump++ = *inp++;
        }

        *nump = '\0';                       /* Null terminate */

        if (inp != period.c_str())   {
            // Found some digits 
            const_cast<int&>(count) = atoi(numberBuff);
        }
        else {
            const_cast<int&>(count) = 1;
        }
        // dangerous to build a string from a char which is 0
        char intervalAsChar = toupper(*inp);
        if (intervalAsChar == '\0' || *(inp+1) != '\0' ||
            !strchr(VALID_MATURITIES, intervalAsChar)){
            if (intervalAsChar != '\0'){
                const_cast<string&>(interval) = intervalAsChar;
            }
            throw ModelException(routine, "Invalid period: '"+interval+
                                 "' for MaturityPeriod "+period+". "+
                                 "Period must be D, M, A, Y, S, Q, W or C");
        }
        const_cast<string&>(interval) = intervalAsChar;
    }
}

MaturityPeriod::MaturityPeriod(const string& period) : 
    Expiry(TYPE), period(period), count(0){
    // populate count and interval
    validatePop2Object();
}

MaturityPeriod::MaturityPeriod(int count, const string& interval) :
    Expiry(TYPE), period(Format::toString(count)+interval), 
    count(count) {
    static const string routine = "MaturityPeriod::MaturityPeriod";
    if (interval.size() != 1 ||
        !strchr(VALID_MATURITIES, toupper(interval[0]))){
        throw ModelException(routine, "Invalid period: '"+interval+
                             "' for MaturityPeriod "+period+". "+
                             "Period must be D, M, A, Y, S, Q, W or C");
    }
    const_cast<string&>(this->interval) = toupper(interval[0]);
}


MaturityPeriod::MaturityPeriod(int frequency)  : 
	Expiry(TYPE) {
	switch(frequency) {
	case 1:
		period = "1Y";
		interval = "A";
		count = 1;
		break;

	case 2:
		period = "6M";
		interval = "S";
		count = 1;
		break;

    case 3:
        period = "4M";
        interval = "M";
        count = 4;
        break;

	case 4:
		period = "3M";
		interval = "M";
		count = 3;
		break;

    case 6:
        period = "2M";
        interval = "M";
        count = 2;
        break;

    case 12:
        period = "1M";
		interval = "M";
		count = 1;
		break;

	default:
        string msg = Format::toString("Invalid frequency: %d", frequency);
		throw ModelException("MaturityPeriod::MaturityPeriod", msg);
	}
}

MaturityPeriod::~MaturityPeriod() {
    // empty
}

string MaturityPeriod::toString() const {
    return period;
}

DateTime MaturityPeriod::toDate(const DateTime& aDate) const {
    static const string method = "MaturityPeriod::toDate";
    try {
        return toDate(count, interval, aDate);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    } 
}

// Using the year & month fields of a MonthDayYear,
// sets the day field to be the Nth (whichOne) occurence
// of a particular day of the week. Thus to get the 3rd
// Thursday in a month, set whichOne to 3, and dowDesired
// to DateTime::THURSDAY.
static void nthWeekDayOfMonth(
    int                     dowDesired,
    int                     whichOne,
    DateTime::MonthDayYear& mdy) {
    static const string method = "MaturityPeriod::nthWeekDayOfMonth";
    try {
        mdy.day = (whichOne-1)*7 + 1;
        // see if valid
        DateTime date = mdy.toDateTime();
        int      dow  = date.getWeekday();

        if (dow <= dowDesired) {
            mdy.day += (dowDesired - dow);
        }
        else {
            mdy.day += (dowDesired + 7 - dow);
        }
        // might go off end of month, but private and would get
        // trapped when we flip back to a date
    }
    catch (exception& e) {
        throw ModelException(e, method);
    } 
}

static DateTime immGeneralForward(
    DateTime startDate,
    int      numPeriods,
    int      monthsInPeriod,
    int      dayOfWeek,
    int      whichDow) {
    static const string method = "MaturityPeriod::immGeneralForward";
    try {
        if (numPeriods == 0) {
            return startDate;  // easy
        }

        DateTime::MonthDayYear mdy = startDate.toMDY();

        // Determine whether we are before or after the desired day
        // of the week in the next possible month. For example, if
        // monthsInPeriod is 3, then possible months are 3,6,9,12.
        // If monthsInPeriod is 1, then we use the current month.
        mdy.month = ((mdy.month - 1)/monthsInPeriod + 1)*monthsInPeriod;

        // Update the day field of the mdy to be the required dayOfWeek
        nthWeekDayOfMonth(dayOfWeek, whichDow, mdy);

        DateTime baseDate = mdy.toDateTime();
        int      jumpPeriods;

        // If baseDate is before(after) startDate when we're going
        // forward (backwards), still need to skip
        // numPeriods periods. If baseDate is after (before) startDate, 
        // we've already done one skip.
        if ((numPeriods > 0 && startDate.isGreaterOrEqual(baseDate)) ||
            (numPeriods < 0 && baseDate.isGreaterOrEqual(startDate))) {
            jumpPeriods = numPeriods ;
        }
        else {
            if (numPeriods > 0) {
                jumpPeriods = numPeriods - 1 ;
            }
            else {
                jumpPeriods = numPeriods + 1 ;
            }
        }
        
        // Jump forward the required # months
        mdy.month += monthsInPeriod*jumpPeriods;

        // Update the day field of the mdy to be the required dayOfWeek
        nthWeekDayOfMonth(dayOfWeek, whichDow, mdy);
        
        return mdy.toDateTime();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    } 
}

/**
 * Returns Credit IMM dates
 * Eg1: if monthsInPeriod=3 and dayOfMonth=20, resul will
 *      be nth (n=numPeriods) 20-Mar, 20-Jun, 20-Sep or 20-Dec following
 *      startDate.
 *      Eg1.1: if startDate=20-Jan-2006 and numPeriod=3, then function
 *             will return 20-Sep-20006
 *      Eg1.2: if startDate=20-Mar-2006 and numPeriod=1, then function
 *             will return 20-Jun-2006
 *      Eg1.3: if startDate=20-Mar-2006 and numPeriod=-1, then function
 *             will return 20-Dec-2005
 *      Eg1.4: if startDate="any date" and numPeriod=0, then function
 *             will return "any date"
 * */
static DateTime creditIMMDate(
    DateTime startDate,
    int      numPeriods,
    int      monthsInPeriod,
    int      dayOfMonth) {
    static const string method = "MaturityPeriod::creditIMMDate";
    try {
        if (numPeriods == 0) {
            return startDate;
        }

        DateTime::MonthDayYear mdy = startDate.toMDY();

        // Determine base date
        if (mdy.month % monthsInPeriod != 0) {
            mdy.month = ((mdy.month - 1)/monthsInPeriod + 1)*monthsInPeriod;
        } else {
            if (mdy.day >= dayOfMonth) {
                if (mdy.day != dayOfMonth || numPeriods > 0) {
                    mdy.month += monthsInPeriod;
                }
            }
        }
        mdy.day = dayOfMonth;
        
        // Jump the required number of periods
        if (numPeriods > 0) {
            mdy.month += (numPeriods-1) * monthsInPeriod;
        } else {
            mdy.month += numPeriods * monthsInPeriod;
        }
        
        return mdy.toDateTime();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    } 
}

/** Answers whether the specified interval produces dates purely by rolling some number of months */
/** static */
bool MaturityPeriod::isMonthBasedInterval(const string& interval) {
    switch (interval[0]) {
        case 'A':
        case 'M':
        case 'Q':
        case 'S':
        case 'Y':
            return true;
            break;
        default:
            return false;
            break;
    }
}

DateTime MaturityPeriod::toDate(
    int           count, 
    const char*   interval,
    const DateTime& aDate) 
{
    static const string method = "MaturityPeriod::toDate";
    try {
        switch (interval[0])
        {
        case 'D':                      
            return aDate.rollDate(count);
        case 'W':
            return aDate.rollDate(7 * count);
        case 'I':  // IMM dates (quarterly)
            return immGeneralForward(aDate, count, 3, DateTime::WEDNESDAY, 3);
        case 'J':  // IMM dates (monthly)
            return immGeneralForward(aDate, count, 1, DateTime::WEDNESDAY, 3);
        case 'K':  // Aussie IMM dates (quarterly)
            return immGeneralForward(aDate, count, 3, DateTime::FRIDAY, 2);
        case 'M':  // MONTHly increments
            return aDate.rollDateInMonths(count);
        case 'A':  // ANNUAL increments
        case 'Y':  // YEARly increments
            return aDate.rollDateInMonths(12 * count);
        case 'S':  // SEMIANNUAL increments
            return aDate.rollDateInMonths(6 * count);
        case 'Q':  // QUARTERly increments
            return aDate.rollDateInMonths(3 * count);
        case 'C':  // CREDIT IMM dates (quarterly)
            return creditIMMDate(aDate, count, 3, 20);
        case 'T':  // Equity derivatives expiry dates (monthly)
            return immGeneralForward(aDate, count, 1, DateTime::FRIDAY, 3);
        default:                
            throw ModelException(method, "invalid period " +
                                 Format::toString(count) + interval);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed for period " + 
                             Format::toString(count) + interval);
    } 
}

/** as above but takes a string */
DateTime MaturityPeriod::toDate(
    int             count, 
    const string&   interval,
    const DateTime& aDate){
    return toDate(count, interval.c_str(), aDate);
}

/** write object out to writer */
void MaturityPeriod::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(toString());
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(e, "MaturityPeriod::write");
    }
}

/** populate an empty object from reader */
void MaturityPeriod::import(Reader::Node* elem, Reader* reader) {
    static const string method = "MaturityPeriod::import";
    try {
        // import has privileged access
        const_cast<string&>(period) = elem->value();
        validatePop2Object();   
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
} 

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void MaturityPeriod::outputWrite(const string& linePrefix,
                                 const string& prefix,
                                 ostream&      stream) const{
    stream << linePrefix << prefix << ": " << period << endl;
}

/** for use with SwapTool methods */
void MaturityPeriod::decompose(int& count, string& interval) const {
    count    = this->count;
    interval = this->interval;
}

/** subtract two dates (d1 - d2) and return difference as a period */
MaturityPeriod* MaturityPeriod::dateSubtract(
    const DateTime& d1, const DateTime& d2) {
    static const string method = "MaturityPeriod::dateSubtract";
    try {
        DateTime::MonthDayYear mdy1 = d1.toMDY();
        DateTime::MonthDayYear mdy2 = d2.toMDY();

        // check for months 
        if (mdy1.day == mdy2.day) {
            if (mdy1.month == mdy2.month) {
                int years = mdy1.year - mdy2.year;
                return new MaturityPeriod(years, "Y");
            }
            else {
                int months = 12*(mdy1.year - mdy2.year) + mdy1.month - mdy2.month;
                return new MaturityPeriod(months, "M");
            }
        }


        // check for last day of month
        if (mdy2.day > mdy1.day) {
            if (d1.isEndOfMonth() && d2.isEndOfMonth()) {
                int months = 12*(mdy1.year-mdy2.year) + mdy1.month-mdy2.month;
                return new MaturityPeriod(months, "M");
            }
        }

        // failing that, just use days
        int days = d1.daysDiff(d2);
        if (days%7 == 0) {
            return new MaturityPeriod(days/7, "W");
        }
        else {
            return new MaturityPeriod(days, "D");
        }
    }    
    catch (exception& e) {
        throw ModelException(e, method);
    } 
}

/** convert to equivalent number of years */
double MaturityPeriod::toYears() const {
    static const string method = "MaturityPeriod::toYears";
    try {
        double years;
        switch (interval[0])
        {
        case 'M':                     // MONTHly increments 
        case 'J':                     // IMM dates (monthly)
            years = count/12.0;
            break;
        case 'A':                     /* ANNUAL increments */
        case 'Y':                     /* YEARly increments */
            years = count;
            break;
        case 'S':                     /* SEMIANNUALL increments */
            years = count/2.0;
            break;
        case 'Q':                     /* QUARTERly increments */
        case 'I':                     // IMM dates (quarterly)
        case 'K':                     // Aussie IMM dates (quarterly)
        case 'C':                     // CREDIT IMM dates (quarterly)
            years = count/4.0;
            break;
        case 'D':                      
            years = count/365.0;
            break;
        case 'W':
            years = count*7.0/365.0;
            break;
        default:
            throw ModelException(method, "invalid period "+ 
                                 Format::toString(count) + interval);
        }
        return years;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    } 
}

int MaturityPeriod::toMonths() const {
    try {
        switch (interval[0])
        {
        case 'M':                     // MONTHly increments 
        case 'J':                     // IMM dates (monthly)
            return count;
        case 'A':                     /* ANNUAL increments */
        case 'Y':                     /* YEARly increments */
            return 12*count;
        case 'S':                     /* SEMIANNUALL increments */
            return 6*count;
        case 'Q':                     /* QUARTERly increments */
        case 'I':                     // IMM dates (quarterly)
        case 'K':                     // Aussie IMM dates (quarterly)
        case 'C':                     // CREDIT IMM dates (quarterly)
            return 3*count;
        default:
            throw ModelException(toString()+" cannot be converted to a number of months");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "MaturityPeriod::toMonth");
    } 
}

int MaturityPeriod::annualFrequency() const {
	const char * method = "MaturityPeriod::annualFrequency";
	int months = toMonths();
	if (months>12 || (12%months!=0)) throw ModelException(method, Format::toString(months)+" cannot be converted to an annual frequency");
	return 12/months;
}

int MaturityPeriod::approxAnnualFrequency() const {
	return static_cast<int>(0.5 + 1.0 / toYears());
}

bool MaturityPeriod::isZeroLength() const {
    return count == 0;
}

/* for reflection */
MaturityPeriod::MaturityPeriod():Expiry(TYPE), count(0){}

class MaturityPeriodHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MaturityPeriod, clazz);
        SUPERCLASS(Expiry);
        EMPTY_SHELL_METHOD(defaultMaturityPeriod);
        clazz->enableCloneOptimisations(); // assignment operator clones
        FIELD(period, "period eg 1W");
        // interval and count not registered
    }
    
    static IObject* defaultMaturityPeriod(){
        return new MaturityPeriod();
    }
};

CClassConstSP const MaturityPeriod::TYPE = CClass::registerClassLoadMethod(
    "MaturityPeriod", typeid(MaturityPeriod), MaturityPeriodHelper::load);

//template<> CClassConstSP const array<MaturityPeriodSP, MaturityPeriod>::TYPE = 
//CClass::registerClassLoadMethod(
//    "MaturityPeriodArray", typeid(MaturityPeriodArray), load);
DEFINE_TEMPLATE_TYPE(MaturityPeriodArray);

class DateFwdAddin: public CObject{
    static CClassConstSP const TYPE;

    DateTime date;     
    string   period;

    /** the 'addin function' */
    static IObjectSP fwdDate(DateFwdAddin* params){
        static const string routine = "DateFwdAddin::fwdDate";
        try {
            MaturityPeriod period(params->period);
            return (IObjectSP(period.toDate(params->date).clone()));
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    DateFwdAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DateFwdAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDateFwdAddin);
        FIELD(date, "date");
        FIELD(period, "period");

        Addin::registerClassObjectMethod("DATE_FWD",
                                         Addin::UTILITIES,
                                         "Add an interval to a date",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)fwdDate);
    }

    static IObject* defaultDateFwdAddin(){
        return new DateFwdAddin();
    }
    
};

CClassConstSP const DateFwdAddin::TYPE = CClass::registerClassLoadMethod(
    "DateFwdAddin", typeid(DateFwdAddin), load);

class DateOffsetAddin: public CObject, virtual public ClientRunnable {
    static CClassConstSP const TYPE;

    DateTime date;     
    string   period;

    /** the 'addin function' */
    static IObjectSP offset(DateOffsetAddin* params){
        static const string routine = "DateOffsetAddin::offset";
        try {
            MaturityPeriod period(params->period);
            DateTime       date = period.toDate(params->date);
            return (IObjectSP(CInt::create(date.daysDiff(params->date))));
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    // EdrAction version of addin
    IObjectSP run() {
        return offset(this);
    }

    /** for reflection */
    DateOffsetAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setDescription("returns days between date & date + interval");
        REGISTER(DateOffsetAddin, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultDateOffsetAddin);
        FIELD(date, "date");
        FIELD(period, "period");

        Addin::registerClassObjectMethod("DATE_OFFSET",
                                         Addin::UTILITIES,
                                         "Offset of interval from a date",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)offset);
    }

    static IObject* defaultDateOffsetAddin(){
        return new DateOffsetAddin();
    }
    
};

CClassConstSP const DateOffsetAddin::TYPE = CClass::registerClassLoadMethod(
    "DateOffsetAddin", typeid(DateOffsetAddin), load);

class DateIntervalAddin: public CObject{
    static CClassConstSP const TYPE;

    DateTime dateFrom;     
    DateTime dateTo;
    bool     ignoreTimeOfDay;

    /** the 'addin function' */
    static IObjectSP dateInterval(DateIntervalAddin* params){
        static const string routine = "DateIntervalAddin::dateInterval";
        try {
            DateTime dateFrom(params->dateFrom);
            DateTime dateTo(params->dateTo);
            if (params->ignoreTimeOfDay){
                dateFrom = DateTime(dateFrom.getDate(), DateTime::START_OF_DAY_TIME);
                dateTo = DateTime(dateTo.getDate(), DateTime::START_OF_DAY_TIME);
            }
            MaturityPeriodSP period(MaturityPeriod::dateSubtract(dateTo, dateFrom));
            return (IObjectSP(CString::create(period->toString())));
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    DateIntervalAddin():
    CObject(TYPE),
    ignoreTimeOfDay(false){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DateIntervalAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDateIntervalAddin);
        FIELD(dateFrom, "dateFrom");
        FIELD(dateTo, "dateTo");
        FIELD(ignoreTimeOfDay, "ignoreTimeOfDay");
        FIELD_MAKE_OPTIONAL(ignoreTimeOfDay);

        Addin::registerClassObjectMethod("DATE_INTERVAL",
                                         Addin::UTILITIES,
                                         "Returns a date interval",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)dateInterval);
    }

    static IObject* defaultDateIntervalAddin(){
        return new DateIntervalAddin();
    }
    
};

CClassConstSP const DateIntervalAddin::TYPE = CClass::registerClassLoadMethod(
    "DateIntervalAddin", typeid(DateIntervalAddin), load);

DRLIB_END_NAMESPACE
