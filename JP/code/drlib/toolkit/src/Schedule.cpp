//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Schedule.cpp
//
//   Description : A schedule of dates & levels you can interpolate on
//
//   Author      : Andrew J Swain
//
//   Date        : 9 March 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_SCHEDULE_CPP
#include "edginc/Schedule.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<Schedule>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Schedule>);
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<ScheduleSP>(ScheduleSP* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<ScheduleSP>(ScheduleSP* t,
                                                       IObjectSP o));

const string Schedule::INTERP_NONE   = "N";
const string Schedule::INTERP_LINEAR = "L";
const string Schedule::INTERP_STAIRS = "S";

const int Schedule::ADJUST_NONE            = 0;
const int Schedule::ADJUST_CONSTANT        = 1;
const int Schedule::ADJUST_INTERPOLATE     = 2;


Schedule::Schedule(
        const DateTimeArray& dates,
        const DoubleArray&   values,
        const string&        interp):
CObject(TYPE),
dates(new DateTimeArray(dates)),
values(new DoubleArray(values)),
interp(interp) {
    validatePop2Object();
}

Schedule::Schedule(const Schedule& rhs,
                   double          value):
CObject(TYPE),
dates(rhs.dates),
interp(rhs.interp){
    values = DoubleArrayConstSP(new DoubleArray(rhs.values->size(), value));
}


Schedule::~Schedule() {}

/** overridden for performance */
IObject* Schedule::clone() const{
    Schedule* copy = new Schedule();
    copy->dates = this->dates;   // shallow copy the SP's over
    copy->values = this->values; // shallow copy the SP's over
    copy->interp = this->interp;
    return copy;
}

DateTime Schedule::lastDate() const {
    this->isEmpty("Schedule::lastDate");
    return ((*dates)[length()-1]);
}

double Schedule::lastValue() const {
    this->isEmpty("Schedule::lastValue");
    return ((*values)[length()-1]);
}

DateTime Schedule::firstDate() const {
    this->isEmpty("Schedule::firstDate");
    return ((*dates)[0]);
}

double Schedule::firstValue() const {
    this->isEmpty("Schedule::firstValue");
    return ((*values)[0]);
}

int Schedule::length() const {
    return dates->getLength();
}

/** Determines if all rates in the schedule are equal */
bool Schedule::isFlat()const
{
    this->isEmpty("Schedule::isFlat");
    
    bool isFlat = true;

    for (int idx = 1; idx < dates->size(); idx++)
    {
        if (!Maths::equals((*values)[idx], (*values)[0]))
        {
            isFlat = false;
            break;
        }
    }
    return isFlat;
}

/** Determines if all times in schedule are timeToCompare */
bool Schedule::timesAreAll(int timeToCompare) const
{
    this->isEmpty("Schedule::timesAreAll");

    for (int idx = 0; idx < dates->size(); idx++)
    {
        if ((*dates)[idx].getTime() != timeToCompare)
        {
            return false;
        }
    }
    return true;
}

/** Determine if a date is contained inside a schedule */
bool Schedule::coversDateRange(const DateTime& lowDate, const DateTime& highDate, bool datesDefineLimit)const
{
    static const string method = "Schedule::coversDateRange";
    this->isEmpty(method);

    bool coversDateRange;
    
    if (lowDate.isGreater(highDate)) {
        throw ModelException(method,
                             "lowDate passed " + lowDate.toString() +
                             " is after high date " + highDate.toString());
    }

    if (!datesDefineLimit) {
        if (interp != INTERP_NONE) {
            if ((*dates)[0].isGreater(lowDate) ||
                highDate.isGreater((*dates)[dates->size()-1])) {
                coversDateRange = false;
            }
            else {
                coversDateRange = true;
            }
        }
        else {
            /* interpolation not allowed so find exact dates using 2 binary searches 
               loop twice looking for the highDate and lowDate */
            bool lowFound = false, highFound = false;
            
            if (dates->size() == 0) {
                coversDateRange = false;
            }
            else if (dates->size() == 1) {
                coversDateRange = (*dates)[0].equals(lowDate) &&  (*dates)[0].equals(highDate);
            }
            else if (dates->size() == 2) {
                coversDateRange = ((*dates)[0]==lowDate ||  (*dates)[1]==lowDate) 
                    && ((*dates)[0]==highDate ||  (*dates)[1]==highDate) ;
            }
            else {
                for (int loopIter = 0; loopIter < 2; loopIter++) {
                    int lo = 0, hi = dates->size()-1, mid;    
                    if (loopIter == 1 && highFound == false) {
                        // no need to do the second search
                        break;
                    }
                    while((hi - lo) > 1) {
                        mid = (lo + hi) >> 1;  // compute a mid point;
                        if ((*dates)[mid].equals(loopIter?lowDate : highDate)) {
                            if (loopIter) {
                                lowFound = true;
                            }
                            else {
                                highFound = true;
                            }
                            break;
                        }
                        if ((*dates)[mid].isGreater(lowDate)) {
                            hi = mid;
                        }
                        else {
                            lo = mid;
                        }
                    }
                }
            
                coversDateRange = lowFound && highFound;
            }
        }
    }
    else {
        coversDateRange = lowDate.equals((*dates)[0]) && highDate.equals((*dates)[dates->size()-1]);
    }
    
    return coversDateRange;
}



/** v. basic interpolator */
double Schedule::interpolate(const DateTime& date) const {
    static const string method = "Schedule::interpolate";
    try {
        this->isEmpty(method);

        double value = 0.0; // stop a compiler warning
        bool   hasDate;
        // first see if date is in schedule 
        if (interp == Schedule::INTERP_STAIRS || 
            interp == Schedule::INTERP_LINEAR)
        {
            // takes care of case of length() = 1
            if (date.equals((*dates)[0]))
            {
                value = (*values)[0];
            }
            else
            {
                // sufficient to be within first and last date 
                hasDate = date.isGreaterOrEqual((*dates)[0]) &&
                    !date.isGreater(lastDate());
                if (hasDate) {
                    // determine schedule date >= given date 
                    int i;
                    for (i = 1; /* sched->numItems > 1 if we are here */
                         i < length() && date.isGreater((*dates)[i]);
                         i++); /* empty loop */
                    
                    if (interp == Schedule::INTERP_STAIRS) {
                        if (date.equals((*dates)[i])) {
                            i++;
                        }
                        value = (*values)[i-1];
                    }
                    else  {
                        // simple linear interpolation 
                        double hi_lo = (*dates)[i].getDate() - (*dates)[i-1].getDate();
                        double dt_lo = date.getDate() - (*dates)[i-1].getDate();
                        
                    if (hi_lo <= DBL_EPSILON) {
                        value =  (*values)[i-1];
                    }
                    else {
                        value = (*values)[i-1] + ((*values)[i] - (*values)[i-1]) * dt_lo/hi_lo;
                    }
                    }
                }
                else {
                    throw ModelException(method, 
                                         "date " + date.toString() + 
                                         " out of schedule bounds [" + 
                                         (*dates)[0].toString() + ", " + 
                                         (*dates)[length()-1].toString() + "]");
                }
            }
        }
        else if (interp == Schedule::INTERP_NONE) {
            int i;
            for (i = 0, hasDate = false; i < length() && !hasDate; i++) {
                if (date.equals((*dates)[i], false)) {
                    hasDate = true;
                    value   = (*values)[i];
                }
            }

            if (!hasDate) {
                throw ModelException(method, 
                                     "date " + date.toString() + 
                                     " not found in schedule (NB interp type = None)");
            }
        }
        else {
            throw ModelException(method, "invalid interp type " + interp);
        }
        return value;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

// As above but allows just date (without time)
double Schedule::interpolate(const DateTime::Date& date) const {
    this->isEmpty("Schedule::interpolate");

    DateTime dt(date.getDate(), DateTime::START_OF_DAY_TIME);
    return interpolate(dt);
}

bool Schedule::interpolateCVBSchedule(int      adjustmentType,
                                      double   redemptionValue,
                                      DateTime bondMaturity,
                                      DateTime interpDate,
                                      double   &level) const
{
    static const string method = "Schedule::interpolateCVBSchedule";
    this->isEmpty(method);

    bool dateIsInSchedule = false;
    int i;
    level = 0;

    // cut out early if we can
    if (dates->size() <= 0) {
        dateIsInSchedule = false;
        level = 0.;
        return dateIsInSchedule;
    }

    if ( interp == Schedule::INTERP_NONE ) {
        // None - Same as regular schedule
        dateIsInSchedule = false;
        level = 0;
        for (i=0; i<dates->size(); i++) {
            if (interpDate == (*dates)[i]) {
                dateIsInSchedule = true;
                level = (*values)[i];
                break;
            }
        }
    } else if (interp == Schedule::INTERP_STAIRS) {
        // Stairs - Different from regular schedule in that after last date it returns
        //     last amount
        if (interpDate >= (*dates)[0]) {
            dateIsInSchedule = true;
            for (i=dates->size(); i>0; i--){
                if (interpDate >= (*dates)[i-1]) {
                    level = (*values)[i-1];
                    break;
                }
            }
        } else {
            dateIsInSchedule = false;
            level = 0.;
        }
    } else if (interp == Schedule::INTERP_LINEAR) {
        // Linear - Different from regular schedule in that it automatically extends to
        //     maturity. Maturity date level is the bonds redemption if redemption is
        //     greater than last level in the schedule, otherwise it is set to the 
        //     last level in schedule - matAdjust flag controls if redemption is used as well
        if (interpDate >= (*dates)[0]) {
            dateIsInSchedule = true;
            if (interpDate >= (*dates)[dates->size()-1]) { // handle extention to maturity
                if ( adjustmentType == ADJUST_NONE ) {
                    if (interpDate == (*dates)[dates->size()-1]) {
                        level = (*values)[dates->size()-1];
                    } else {
                        dateIsInSchedule = false;
                        level = 0.0;
                    }
                } else {
                    if ( adjustmentType == ADJUST_CONSTANT || (*values)[dates->size()-1] >= redemptionValue ||
                        bondMaturity.getDate() == (*dates)[dates->size()-1].getDate()) {
                        level = (*values)[dates->size()-1];
                    } else if ( adjustmentType == ADJUST_INTERPOLATE ) {
                        level = (*values)[dates->size()-1] + 
                            (redemptionValue - (*values)[dates->size()-1])/
                            (bondMaturity.getDate() - (*dates)[dates->size()-1].getDate())*
                            (interpDate.getDate() - (*dates)[dates->size()-1].getDate());
                    } else {
                        throw ModelException(method, "Unknown CVB schedule adjustment type " + Format::toString(adjustmentType));
                    }
                }
            } else { // we know we're in the middle
                for (i=dates->size()-1; i>0; i--){
                    if (interpDate >= (*dates)[i-1]) {
                        level = (*values)[i-1] + ((*values)[i] - (*values)[i-1])/
                            ((*dates)[i].getDate() - (*dates)[i-1].getDate())*
                            (interpDate.getDate() - (*dates)[i-1].getDate());
                        break;
                    }
                }
            }
        } else {
            dateIsInSchedule = false;
            level = 0.;
        }

    } else {
        throw ModelException(method, "Unknown interpolation type " + interp);
    }

    return dateIsInSchedule;
}


/** scale rates by a double */
void Schedule::scale(double factor)
{
    /* "copy on write" approach (see clone()) */
    int size = length();
    DoubleArraySP newValues(new DoubleArray(size));
    for (int i = 0; i < size; i++)
    {
        (*newValues)[i] = (*values)[i] * factor;
    }
    values = newValues;
}

/* This tests the rates of a schedule to see if they are all
   greater than or equal to zero. Slightly different from
   IsNegative() in that IsPositive() tests for greater than zero. */
bool Schedule::isNonNegative()const
{
    bool result = true;
    
    for (int idx = 0; idx < dates->size(); idx++)
    {
        if (Maths::isNegative((*values)[idx]))
        {
            result = false;
            break;
        }
    }
    return result;
}

/* This tests the rates of a schedule to see if they are all strictly 
   greater than zero. */
bool Schedule::isPositive() const
{
    this->isEmpty("Schedule::isPositive");

    bool positive = true;
    
    for (int idx = 0; idx < dates->size(); idx++) {
        if (!Maths::isPositive((*values)[idx])) {
            positive = false;
            break;
        }
    }
    return positive;
}


/** This tests the rates of a schedule to see if they are all
   almost equal to zero. */
bool Schedule::isAlmostZero()const
{
    this->isEmpty("Schedule::isAlmostZero");

    bool result = true;
    
    for (int idx = 0; idx < dates->size(); idx++)
    {
        if (!Maths::isZero((*values)[idx]))
        {
            result = false;
            break;
        }
    }
    return result;
}

bool Schedule::getNextDate(const DateTime& fromDate, 
                           const DateTime& bondMaturity,
                           DateTime&       nextDate) const
{
    static const string method = "Schedule::getNextDate";
    this->isEmpty(method);
    
    bool foundDate = false;
    try {
        if (dates->size() == 0 ) {
           return false;
        }

        if (interp == Schedule::INTERP_STAIRS ||
            interp == Schedule::INTERP_LINEAR)
        {
            if ( fromDate <= (*dates)[0] ) {
               nextDate = (*dates)[0];
               foundDate = true;
            } else if ( fromDate > (*dates)[dates->size()-1] && fromDate > bondMaturity) {
               foundDate = false;
            } else {
               nextDate = fromDate;
               foundDate = true;              
            }
        }
        else if (interp == Schedule::INTERP_NONE) {
            int i=0;
            if ( fromDate < (*dates)[0] ) {
               nextDate = (*dates)[0];
               foundDate = true;
            } else if ( fromDate > (*dates)[dates->size()-1] && fromDate > bondMaturity ) {
               foundDate = false;
            } else {
               while (i<dates->size() && (*dates)[i] < fromDate ) {
                 ++i;
               }
               if ( i == dates->size()) {
                  foundDate = false; 
               } else {
                  nextDate = (*dates)[i];
                  foundDate = true;
               }
            }
        }
        else {
            throw ModelException(method, "invalid interp type " + interp);
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }

    if ( foundDate && nextDate <= bondMaturity ) {
        return true;
    } else {
        return false;
    }
}

void Schedule::validatePop2Object() {
    static const string method = "Schedule::validatePop2Object";
    try {
        if ( dates->size() != values->size() ) { 
           throw ModelException(method,
                                "number of dates (" + 
                                Format::toString(dates->size()) + 
                                ") must be the same as number of values (" +
                                Format::toString(values->size()) + ").");
        }

        for (int i = 1; i < length(); i++) {
            if (!(*dates)[i].isGreater((*dates)[i-1])) {
                if ((*dates)[i] == (*dates)[i-1]) {
                    // if dates and values are identical for [i] and [i-1],
                    // then this is not an error
                    if ((*values)[i] != (*values)[i-1]) {
                        throw ModelException(method,
                            "dates["+Format::toString(i-1)
                            +"] and dates["+Format::toString(i)
                            +"] are identical ("
                            +(*dates)[i].toString()+") but values "
                            "are different ("+Format::toString((*values)[i-1])
                            +" and " + Format::toString((*values)[i])+")");
                    }
                }
                else {
                    throw ModelException(method,
                                        "dates not in increasing order (" +
                                        (*dates)[i].toString() + " <= " +
                                        (*dates)[i-1].toString() + ")");
                }
            }
        }

        if (interp != Schedule::INTERP_NONE &&
            interp != Schedule::INTERP_LINEAR &&
            interp != Schedule::INTERP_STAIRS) {
            throw ModelException(method, "invalid interp type " + interp+
                                 ". Must be (N)one, (L)inear or (S)tairs");
        }
    }            
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void Schedule::isEmpty(string method) const {
    if ( dates->size() == 0) { 
        throw ModelException(method+"::isEmpty",
                             "Empty schedule can not be used");
    }
}

// initialized dates and values to prevent core dump when using an empty schedule
Schedule::Schedule(): CObject(TYPE),
                      dates(new DateTimeArray(0)),
                      values(new DoubleArray(0)) {
    //nothing
}

    /** return interp type */
string Schedule::getInterp() const{
    return interp;
}

/** return date list (deep copy) */
DateTimeArray Schedule::getDates() const{
    return *dates;
}

/** return value list (deep copy) */
DoubleArray Schedule::getValues() const{
    return *values;
} 

/** returns reference to date array */
const DateTimeArray& Schedule::getDateArray() const
{
    return *dates;
}

const DoubleArray& Schedule::getValueArray() const
{
    return *values;
}

/** returns a subset of the schedule that covers the supplied date range.
    If INTERP_NONE, returns any values that lie inside the given dates
    otherwise retuns a value per day for each day between supplied values.
    Advances from lower bound in 1 day increments using best guess at
    time of day to give most complete results */
CashFlowArray* Schedule::subset(const DateTime& lowDate, 
                                const DateTime& highDate) const {
    static const string method("Schedule::subset");
    try {
        this->isEmpty(method);

        CashFlowArraySP cfl(new CashFlowArray(0));
        // do simple case first
        if (interp == INTERP_NONE) {
            for (int i = 0; i < dates->size(); i++) {
                if ((*dates)[i] >= lowDate && (*dates)[i] <= highDate) {
                    CashFlow cf((*dates)[i], (*values)[i]);
                    cfl->push_back(cf);                        
                }
            }
        }
        else {
            DateTime start = firstDate();
            DateTime end   = lastDate();
            if (lowDate <= end && highDate >= start) {
                // some portion of the schedule overlaps with the range
                DateTime interpDate(lowDate.getDate(), start.getTime());
                while (interpDate < highDate) {
                    if (coversDateRange(interpDate, interpDate, false)) {
                        CashFlow cf(interpDate, interpolate(interpDate));
                        cfl->push_back(cf); 
                    }
                    interpDate = interpDate.rollDate(1);
                }
                interpDate = DateTime(highDate.getDate(), end.getTime());
                if (coversDateRange(interpDate, interpDate, false)) { 
                    // pump in upper bound
                    CashFlow cf(interpDate, interpolate(interpDate));
                    cfl->push_back(cf);                 
                }
            }
        }
        return (cfl.release());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** check the dates in schedule are same to refDates*/
bool Schedule::isSameSchedule(const DateTimeArray refDates, bool compTime) const{
    double iSize = this->length();
    if (iSize != refDates.size())
        return false;
    else {
        for (int i=0; i<iSize; i++){
            if ((*dates)[i].equals(refDates[i], compTime) == false)
                return false;
        }
    }
    return true;
}

class ScheduleHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Schedule, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSchedule);
        FIELD(dates, "dates");
        FIELD(values, "values");
        FIELD(interp, "interpolation style");
        Addin::registerConstructor("SCHEDULE",
                                   Addin::UTILITIES,
                                   "Creates a Schedule",
                                   Schedule::TYPE);
    }

    static IObject* defaultSchedule(){
        return new Schedule();
    }
};

CClassConstSP const Schedule::TYPE = CClass::registerClassLoadMethod(
    "Schedule", typeid(Schedule), ScheduleHelper::load);


// this class is not compatible with IMS so do not use it in a public interface
// but only for the private interface (i.e. to have a better code)
DEFINE_TEMPLATE_TYPE(ScheduleArray);

/** Addin to calculate the probability of NO DEFAULT **/ 
class ScheduleAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    ScheduleSP         schedule;
    DateTime           interpDate;

    /** the 'addin function' */
    static double interpolate(ScheduleAddin* params) {
        static const string routine = "ScheduleAddin::interpolate";

        params->schedule->isEmpty(routine);

        if (!params->schedule) {
            throw ModelException(routine, "Schedule must not be null");
        }

        double level;
        double dateFound;
        dateFound = params->schedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                             100.0,
                                                             params->interpDate,
                                                             params->interpDate, 
                                                             level);

        return level;
    }

    /** for reflection */
    ScheduleAddin():  CObject(TYPE){}

    
 /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ScheduleAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultScheduleAddin);
        FIELD(schedule,             "Schedule");
        FIELD(interpDate,    "interpolation date");
        Addin::registerClassDoubleMethod("SCHEDULE_INTERPOLATE",
                                         Addin::CONV_BOND,
                                         "Interpolate a schedule at a given date",
                                         TYPE,
                                         (Addin::DoubleMethod*)interpolate);
    }

    static IObject* defaultScheduleAddin() {
        return new ScheduleAddin();
    }
};

CClassConstSP const ScheduleAddin::TYPE = CClass::registerClassLoadMethod(
    "ScheduleAddin", typeid(ScheduleAddin), load);


DRLIB_END_NAMESPACE
