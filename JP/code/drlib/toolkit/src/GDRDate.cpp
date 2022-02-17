//----------------------------------------------------------------------------
//
//   Group       : Global Derivatives Research
//
//   Filename    : GDRDate.cpp
//
//   Description : Wrapper for the date in DR Interface
//
//   Author      : Mark A Robson
//
//   Date        : 5 Dec 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GDRDate.hpp"
#include "edginc/Format.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Null.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE
#define TDATE_BASE_YEAR     1601
#define DAYS_IN_1_YEAR      365
#define DAYS_IN_4_YEARS     1461
#define DAYS_IN_100_YEARS   (DAYS_IN_4_YEARS * 25 - 1)
#define DAYS_IN_400_YEARS   (DAYS_IN_100_YEARS * 4 + 1)
#define MONTHS_PER_YEAR     12
/* A generally useful macro. */
#define IS_LEAP(year) (                 \
    (((year)%4 == 0) && ((year)%100 != 0)) || \
    ((year)%400 == 0) \
    )

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

/** Converts ALIB date to a DRDate */
static DRDate toDRDate(int tDate){
    int         count;
    int        *cumDaysp;
    int         date = tDate;
    int         year = TDATE_BASE_YEAR;
    int         month;
    int         day;
    
    if (date < 0){
        throw ModelException("toDRDate", "Date ("+Format::toString(date)+ 
                             ") is negative");
    }
    
    /* Get year */
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
     * days than the other groups of 4 years. */
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
    DRDate theDate;
    theDate.year = year;
    theDate.month = month;
    theDate.day = day;
    return theDate;
}

/** Converts to ALIB date from DRDate */
static int toALIBDate(DRDate drDate){
    static const string method("toALIBDate");

    int      dt = 0;  
    int      fourYearBlocks;
    int      year =  drDate.year;
    int      month = drDate.month;
    int      day =   drDate.day;
    bool     isLeap;
    char     buffer[128];
    
    isLeap = IS_LEAP(year);
    year   = year - TDATE_BASE_YEAR;                
    
    // Make sure day is in range
    if (day >= 1 && day <= 28) {
        /*EMPTY*/;                      /* Guaranteed to be OK */
        /* Avoid doing check below */
    } else if (day < 1 ||
               (isLeap ? day > leapDays[month] : day > (days[month]))){
        sprintf(buffer, "invalid date %d-%d-%d", day, month, year);
        throw ModelException(method, buffer);
    }
    
    // Make sure month and year are in range
    if (month < 1 || month > MONTHS_PER_YEAR || year < 0)
    {
        sprintf(buffer, "invalid date %d-%d-%d", day, month, year);
        throw ModelException(method, buffer);
    }
    
    // Take years into account
    while (year >= 400) {
        year -= 400;
        dt   += DAYS_IN_400_YEARS;
    }
    
    while (year >= 100) {
        year -= 100;
        dt   += DAYS_IN_100_YEARS;
    }
    
    if (year >= 4) {
        fourYearBlocks = (int)(year>>2);       /* Divide by 4 */
        year -= (int)(fourYearBlocks<<2);       /* Multiply by 4 */
        dt   += fourYearBlocks * DAYS_IN_4_YEARS;
    }
    
    while (year >= 1) {
        year -= 1;
        dt   += DAYS_IN_1_YEAR;
    }
    
    if (isLeap) {
        dt += leapCumDays[month-1] + day;
    } else {
        dt += cumDays[month-1] + day;
    }
    return dt;
}

GDRDate::~GDRDate(){}

/** Are these objects equal (ie contain the same data) */
bool GDRDate::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const GDRDate* theDate = STATIC_CAST(GDRDate, drObj.get());
    return (theDate->date.year == date.year &&
            theDate->date.month == date.month &&
            theDate->date.day == date.day);
}

const static char* months[] = {"Jan","Feb","Mar","Apr","May","Jun",
                               "Jul","Aug","Sep","Oct","Nov","Dec"};
#define MONTHS_PER_YEAR     12

/** Returns a string representation of this date */
string GDRDate::toString() const{
    return Format::toString("%02d-%s-%d", (int)date.day, months[date.month-1], 
                            (int)date.year);
}

/** write object out to writer */
void GDRDate::write(const string& tag, Writer* writer) const{
    IObjectConstSP obj(writer->objectStart(tag, "", this, false));
    if (obj.get()){
        writer->write(toString());
    }
    writer->objectEnd(tag, this);
}
//// convert dd-mmm-yyyy into a date
DRDate GDRDate::convertDate(const string& date) {
    static const string method = "GDRDate::convertDate";
    try {
        /* see if string is remotely in date format */
        if (date.find("-") == std::string::npos || date.length() != 11) {
            string msg="date string (" + date + ") not in dd-mmm-yyyy format";
            throw ModelException(method, msg);
        }

        // chop it into dd, mmm & yyyy
        int  dash1 = date.find("-");
        int  dash2 = date.rfind("-");
        string d(date.substr(0, 2));
        string m(date.substr(dash1+1, 3));
        string y(date.substr(dash2+1, 4));

        DRDate drDate;
        // day & year are easy
        int day, year;
        sscanf(d.c_str(), "%d", &day);
        drDate.day = day; // we're going to validate the day/year
        sscanf(y.c_str(), "%d", &year);
        drDate.year = year;
        drDate.month = 0;

        bool found = false;
        while (drDate.month < MONTHS_PER_YEAR && !found) {
            if (CString::equalsIgnoreCase(m, months[(int)drDate.month])) {
                found = true;
            }
            drDate.month++; // if found, this converts from 0-1 to 1-12
        }
        if (!found) {
            throw ModelException(method, "couldn't find month ("+m+")");
        }
    
        return drDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns ALIB style date */
int GDRDate::tDate() const{
    return toALIBDate(date);
}

//// constructor
GDRDate* GDRDate::create(int tDate){ // ALIB style date
    return new GDRDate(toDRDate(tDate));
}

/** Create a new date offset from this date */
GDRDate* GDRDate::rollDate(int offset){
    // because we're working with the world's worst date convention ....
    int tDate = toALIBDate(date);
    return create(tDate+offset);
}

/** populate an empty object from description */
void GDRDate::import(Reader::Node* elem, Reader* reader){
    const_cast<DRDate&>(date) = convertDate(elem->value());
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void GDRDate::outputWrite(const string& linePrefix,
                          const string& prefix, ostream& stream) const{
    stream << linePrefix << prefix << ": " << toString() << endl;
}

// just returns this
IObject* GDRDate::clone() const{
    if (getRefCount() == 0){
        return new GDRDate(date);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}    

/** Returns the date being wrapped */
DRDate GDRDate::dateValue() const{
    return date;
}

//// constructors
GDRDate* GDRDate::create(DRDate date){
    return new GDRDate(date);
}
GDRDate* GDRDate::create(const char* date){
    return new GDRDate(convertDate(date));
}


GDRDate::GDRDate(DRDate date): CObject(TYPE), date(date){}

GDRDate::GDRDate(): CObject(TYPE), date() {
    DRDate& date = const_cast<DRDate&>(this->date);
    date.year = 0;
    date.month = 0;
    date.day = 0;
}

IObject* GDRDate::defaultConstructor(){
    return new GDRDate();
}

void GDRDate::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setNative();
    REGISTER(GDRDate, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDRObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
}
CClassConstSP const GDRDate::TYPE = CClass::registerClassLoadMethod(
    "GDRDate", typeid(GDRDate), load);

/** Class for addins that take a single optional date */
class GDRDate::AddinFnc: public CObject{
public:
    static CClassConstSP const TYPE;
    // default constructor
    AddinFnc(): CObject(TYPE){}
    // single optional field. hack it for now by taking an integer
    smartPtr<GDRDate> dateValue;

    ~AddinFnc(){}
private:
    /** addin function - either creates a GDRDate or a Null if no date is
        supplied */
    static IObjectSP createGDRDate(AddinFnc* params){
        if (!params->dateValue){
            return CNull::create();
        }
        return params->dateValue;
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AddinFnc, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAddinFnc);
        FIELD(dateValue, "date");
        FIELD_MAKE_OPTIONAL(dateValue);
        Addin::registerClassObjectMethod("GDRDATE",
                                         Addin::UTILITIES,
                                         "Constructs a handle to a GDR Date",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createGDRDate);
    }

    static IObject* defaultAddinFnc(){
        return new AddinFnc();
    }
};

CClassConstSP const GDRDate::AddinFnc::TYPE = CClass::registerClassLoadMethod(
    "GDRDate::AddinFnc", typeid(AddinFnc), load);


DRLIB_END_NAMESPACE
