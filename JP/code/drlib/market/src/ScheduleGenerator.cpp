//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ScheduleGenerator.cpp
//
//   Description : Various classes to support schedules and schedule
//                 generation.
//
//   Author      : Jerry Cohen
//
//   Date        : November 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SCHEDULEGENERATOR_CPP
#include "edginc/ScheduleGenerator.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Addin.hpp"
#include "edginc/BoxedEnum.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

#define DEFINE_arrayObjectCast(X) \
IObjectConstSP arrayObjectCast<X>::toIObject( \
    const X& value){ \
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value)); \
    return objValue; \
} \
\
IObjectSP arrayObjectCast<X>::toIObject(X& value){ \
    return OptionDatesSP::attachToRef(&value); \
} \
\
const OptionDates& arrayObjectCast<OptionDates>::fromIObject(IObjectSP& value){ \
    OptionDates *dtPtr = DYNAMIC_CAST(X, value.get()); \
    return *dtPtr; \
}

///////////////////////////////////////////////////////////////////////////////
// Some helper functions
///////////////////////////////////////////////////////////////////////////////
static void ensureIncreasing(const DateTime& date1,
                             const DateTime& date2,
                             const string& description) {
    if (date1 > date2) {
        string msg(description);
        msg += " dates must be in increasing order: date1 = "
            + date1.toString() + ", date2 = " + date2.toString();
        throw ModelException(__FUNCTION__, msg);
    }
}

static void ensureStrictlyIncreasing(const DateTime& date1,
                                     const DateTime& date2,
                                     const string& description) {
    if (date1 >= date2) {
        string msg(description);
        msg += " dates must be in strictly increasing order: date1 = "
            + date1.toString() + ", date2 = " + date2.toString();
        throw ModelException(__FUNCTION__, msg);
    }
}
                                    
static void validatePeriodDates(const DateTime& start, const DateTime& end) {
    if (start > end) {
        string msg("Period start date (");
        msg += start.toString() + ") must be <= end date ("
            + end.toString() + ")";
        throw ModelException(__FUNCTION__, msg);
    }
}


///////////////////////////////////////////////////////////////////////////////
// ISched - Top-level interface for all schedules.
///////////////////////////////////////////////////////////////////////////////
ISchedConstSP ISched::flatSchedule() const { return schedule(); }

ISched::ISched(const CClassConstSP& type) : CObject(type) {}

ISched::~ISched() {}

void ISched::load(CClassSP& clazz) {
    clazz->setPublic();     // So that we can create arrays of these.
    REGISTER(ISched, clazz);
    clazz->setDescription("Interface for various dated schedules");
    SUPERCLASS(CObject);
}

CClassConstSP const ISched::TYPE = CClass::registerClassLoadMethod(
    "ISched", typeid(ISched), ISched::load);


///////////////////////////////////////////////////////////////////////////////
// ScheduleGenerator - An executable object that allows clients to generate
//                     a specified schedule with a specified market (for
//                     Holidays).
///////////////////////////////////////////////////////////////////////////////
ScheduleGenerator::ScheduleGenerator(const CClassConstSP& type) :
CObject(type), sched(ISchedSP(0)), market(CMarketDataSP(0)), flatten(false) {}

IObjectSP ScheduleGenerator::runAddin(ScheduleGenerator* scheduleGenerator) {
    return scheduleGenerator->run();
}

IObject* ScheduleGenerator::defaultConstructor() { return new ScheduleGenerator(TYPE); }

void ScheduleGenerator::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(ScheduleGenerator, clazz);
    clazz->setDescription("Executable object that generates a schedule");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    IMPLEMENTS(ClientRunnable);
    FIELD(sched, "Object that defines the schedule to generate");
    FIELD(market, "A market environment that will be used to resolve references to holidays");
    FIELD(flatten, "Whether to answer a flattened schedule (default = false)");
    FIELD_MAKE_OPTIONAL(flatten);

    Addin::registerClassObjectMethod("ScheduleGenerator",
                                     Addin::UTILITIES,
                                     "Generates a schedule as specified via the \"sched\" parameter",
                                     ScheduleGenerator::TYPE,
                                     true,
                                     Addin::returnHandle,
                                     (Addin::ObjMethod*)runAddin);
}

CClassConstSP const ScheduleGenerator::TYPE = CClass::registerClassLoadMethod(
    "ScheduleGenerator", typeid(ScheduleGenerator), ScheduleGenerator::load);

class MarketWrappersAction: public ObjectIteration::IAction {
public:
    const IModel *model;
    const MarketData* market;

    MarketWrappersAction(const IModel *model, const MarketData* market)
        : model(model), market(market) {}

        virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
            MarketObjectWrapper* mow = dynamic_cast<MarketObjectWrapper*>(state.getObject().get());
            if (!mow->getName().empty()) {
                mow->getData(model, market, mow->getMOType());
            }
            return false; // do not go into the market object we just fetched
        }
};

IObjectSP ScheduleGenerator::run() {
    // Populate market object wrappers with actual objects (e.g. Holidays) from environment
    NonPricingModelSP modelForMarketFetching(new NonPricingModel());
    MarketWrappersAction action(modelForMarketFetching.get(), market.get());
    ObjectIteration iter(MarketObjectWrapper::TYPE);
    iter.recurse(action, ISchedSP::constCast(sched));

    return ISchedSP::constCast(flatten ? sched->flatSchedule() : sched->schedule());
}


///////////////////////////////////////////////////////////////////////////////
// IDateSched - Interface for classes that produce a schedule of dates (i.e. a
//              DateTimeArray).
///////////////////////////////////////////////////////////////////////////////
DateTimeArrayConstSP IDateSched::elements() const {
    if (!cachedElemsSP)
        setup();
    return cachedElemsSP;
}

ISchedConstSP IDateSched::schedule() const { return ISchedConstSP(new DateSched(elements())); }

IDateSched::IDateSched(const CClassConstSP& type) :
    ISched(type), cachedElemsSP(DateTimeArrayConstSP(0)) {}

void IDateSched::load(CClassSP& clazz) {
    clazz->setPublic();     // In case people want to create arrays of these.
    REGISTER(IDateSched, clazz);
    clazz->setDescription("Interface for a hierarchy of classes that represent a "
                          "schedule of dates");
    SUPERCLASS(ISched);
    FIELD(cachedElemsSP, "");
    FIELD_MAKE_TRANSIENT(cachedElemsSP);
}

CClassConstSP const IDateSched::TYPE = CClass::registerClassLoadMethod(
    "IDateSched", typeid(IDateSched), IDateSched::load);

DEFINE_TEMPLATE_TYPE(IDateSchedArray);


///////////////////////////////////////////////////////////////////////////////
// DateSched - An explicitly specified date schedule.
///////////////////////////////////////////////////////////////////////////////
void DateSched::setup() const { cachedElemsSP = explicitDates; }

IObject* DateSched::clone() const {
    if (getRefCount() == 0)
        return new DateSched(*this);
    else
        return const_cast<DateSched*>(this);
}

DateSched::DateSched(const DateTimeArrayConstSP& explicitDates) :
    IDateSched(TYPE), explicitDates(explicitDates) {
        if (!explicitDates)
            throw ModelException(__FUNCTION__, "Explicit DateTimeArray cannot be NULL");
}

DateSched::DateSched(const CClassConstSP& type) :
    IDateSched(type), explicitDates(DateTimeArrayConstSP(0)) {}

void DateSched::validatePop2Object() {
    DateTime::ensureIncreasing(*explicitDates, "explicit date schedule", false);
}

IObject* DateSched::defaultConstructor() { return new DateSched(TYPE); }

void DateSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(DateSched, clazz);
    clazz->setDescription("An explicitly specified array of dates");
    SUPERCLASS(IDateSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitDates, dates, "The array of dates");
}

CClassConstSP const DateSched::TYPE = CClass::registerClassLoadMethod(
        "DateSched", typeid(DateSched), DateSched::load);


///////////////////////////////////////////////////////////////////////////////
// DateSchedGen - A flexible, periodic date generator.  Dates are generated
//                "backwards" by applying multiples of a MaturityPeriod from
//                either endDate or backStubStartDate (if specified) to either
//                startDate or frontStubEndDate (if specified).  By default, a
//                short front stub will be produced if startDate is not on
//                cycle w.r.t. endDate, and an explicit frontStubEndDate or
//                backStubStartDate is not specified.  However, these may be
//                specified to produce any combination of short/long front
//                and/or back stubs.  For month-based MaturityPeriods (e.g.
//                'M', 'Q', 'S', 'A', 'Y', where dates cycle on a particular
//                day of the month, a.k.a. the roll day), an explicit rollDay
//                or rollOnLast may be specified.  Otherwise, the rollDay
//                defaults to the day of the endDate or an explicitly specified
//                backStubStartDate.  Checks are made to ensure consistency of
//                inputs and to guard against producing more than one front
//                and/or back stub.
///////////////////////////////////////////////////////////////////////////////
void DateSchedGen::setup() const {
    DateTimeArraySP dateListSP(new DateTimeArray);
    DateTimeArray& dateList = *dateListSP;

    // Dates are generated backwards (i.e. future to past/present)...
    if (includeEnd)
        dateList.push_back(endDate);

    const DateTime& periodicBase = periodicBaseDate();
    const DateTime& periodicCutoff = periodicCutoffDate();
    int i = 0;
    DateTime date;

    do {
        date = MaturityPeriod::toDate(i-- * periodCount, periodInterval, periodicBase);
        if (useExplicitRollDay)
            date = rollDayAdjustedDate(date);

        if (date > startDate && date < endDate)
            dateList.push_back(date);

    } while (date > periodicCutoff);

    if (!frontStubEndDate.empty() && date != frontStubEndDate) {
        // Guard against generating multiple front stubs.
        string msg("Invalid explicit frontStubEndDate: ");
        msg += frontStubEndDate.toString() + " (not on cycle)";
        throw ModelException(__FUNCTION__, msg);
    }

    if (includeStart)
        dateList.push_back(startDate);

    // Reverse the list to answer dates in ascending order...
    for (int low = 0, high = dateList.size() - 1; low < high; ++low, --high) {
        DateTime tmp(dateList[low]);
        dateList[low] = dateList[high];
        dateList[high] = tmp;
    }

    cachedElemsSP = dateListSP;
}

const DateTime& DateSchedGen::periodicBaseDate() const {
    return backStubStartDate.empty() ? endDate : backStubStartDate;
}

const DateTime& DateSchedGen::periodicCutoffDate() const {
    return frontStubEndDate.empty() ? startDate : frontStubEndDate;
}

DateTime DateSchedGen::rollDayAdjustedDate(const DateTime& srcDate) const {
    DateTime::MonthDayYear mdy(srcDate.toMDY());
    mdy.day = min(effectiveRollDay, mdy.daysInMonth(ignoreLeapYears));
    return mdy.toDateTime(srcDate.getTime());
}

// Answers true if the two dates are not separated by one whole period (i.e. MaturityPeriod).
bool DateSchedGen::isStubPeriod(const DateTime& start, const DateTime& end) const {
    // We assume that the period is a stub period unless proven otherwise.
    bool isStub = true;

    DateTime periodicDate(MaturityPeriod::toDate(-periodCount, periodInterval, end));

    if (effectiveRollDay > 0) {
        // When dates are produced via an interval that implies that a roll day is used
        // (i.e. an interval that is equivalent to a multiple of months), we must consider
        // that the days of the period start and end dates may be different when the roll
        // day is near the end of the month (e.g. the 29th, 30th, or 31st).
        if (start == rollDayAdjustedDate(periodicDate) &&
            end == rollDayAdjustedDate(MaturityPeriod::toDate(periodCount, periodInterval, start)))
                isStub = false;
    }
    else if (start == periodicDate && end == MaturityPeriod::toDate(periodCount, periodInterval, start)) {
        // For other types of intervals (e.g. daily, weekly, IMM, etc), we take the
        // conservative approach of ensuring that applying the interval to the start
        // date or applying the reverse interval to the end date yields the same result.
        isStub = false;
    }

    return isStub;
}

// If the maturityPeriod is equivalent to some number of months, then dates will be generated
// on a particular day of the month; which we call the roll day.  Normally, the roll day is
// given by the day of the month of the endDate (or explicit backStubStartDate).  However, an
// explicit rollDay or rollOnLast may be specified, but some checks are made to ensure
// consistency with an explicitly specified backStubStartDate.
//
// This method sets:
//      * "effectiveRollDay" to either a valid day of the month or zero (if a roll day is not
//        applicable for the specified maturityPeriod).
//      * "useExplicitRollDay" to true iff we need to explicitly adjust the day of dates
//        generated from endDate (or an explicit backStubStartDate) to be on
//        min(generatedDate.daysInMonth(ignoreLeapYears), effectiveRollDay).
void DateSchedGen::setEffectiveRollDay() {
    // Roll day only applies if we are using a month-based interval.
    if (MaturityPeriod::isMonthBasedInterval(periodInterval)) {

        DateTime::MonthDayYear baseMDY(periodicBaseDate().toMDY());

        if (rollOnLast || rollDay != 0) {

            if (rollOnLast && rollDay != 0 && rollDay != 31)
                throw ModelException(__FUNCTION__, "Cannot specify both an explicit rollDay and rollOnLast");

            int baseDaysInMonth = baseMDY.daysInMonth(ignoreLeapYears);

            // Always set effectiveRollDay, but only set useExplicitRollDay flag if necessary...
            if (rollOnLast) {
                if (!backStubStartDate.empty() && baseMDY.day != baseDaysInMonth) {
                    string msg("When using an explicit backStubStartDate and rollOnLast is true, "
                               "backStubStartDate (");
                    msg += backStubStartDate.toString() + ") must be the last day of the month.";
                    throw ModelException(__FUNCTION__, msg);
                }
                effectiveRollDay = 31;
            }
            else {
                if (!backStubStartDate.empty() && rollDay != baseMDY.day && rollDay <= baseDaysInMonth) {
                    string msg("When using an explicit backStubStartDate and rollDay, "
                               "rollDay must be either greater than number of days "
                               "in the month of backStubStartDate or the same day of "
                               "the month as backStubStartDate");
                    throw ModelException(__FUNCTION__, msg);
                }
                if (rollDay < 1 || rollDay > 31)
                    throw ModelException(__FUNCTION__, "Explicit rollDay must be between 1 and 31");
                
                effectiveRollDay = rollDay;
            }

            if (effectiveRollDay != baseMDY.day)
                useExplicitRollDay = true;
        }
        else
            effectiveRollDay = baseMDY.day;
    }
    else if (rollOnLast || rollDay != 0)
        throw ModelException(__FUNCTION__, "When using an explicit rollDay or rollOnLast is true, "
                                           "maturityPeriod must be equivalent to some multiple of months");
}

void DateSchedGen::validatePop2Object() {
    if (endDate < startDate)
    {
        string msg("Start date (");
        msg += startDate.toString() + ") must be on or before end date (" + endDate.toString() + ").";
        throw ModelException(__FUNCTION__, msg);
    }

    if (!frontStubEndDate.empty() && (frontStubEndDate < startDate || frontStubEndDate > endDate))
        throw ModelException(__FUNCTION__, "frontStubEndDate must be between the start and end dates (inclusive)");

    if (!backStubStartDate.empty() && (backStubStartDate < startDate || backStubStartDate > endDate))
        throw ModelException(__FUNCTION__, "backStubStartDate must be between the start and end dates (inclusive)");

    if (!frontStubEndDate.empty() && !backStubStartDate.empty() && backStubStartDate < frontStubEndDate)
    {
        string msg("Explicit front stub end date (");
        msg += frontStubEndDate.toString() + ") must be on or before explicit back stub start date ("
            + backStubStartDate.toString() + ").";
        throw ModelException(__FUNCTION__, msg);
    }

    // Decompose MaturityPeriod and canonicalize count so that we can apply negative multiples
    // of the count to generate successive periodic dates.
    maturityPeriod->decompose(periodCount, periodInterval);
    periodCount = abs(periodCount);
    if (periodCount == 0)
        throw ModelException(__FUNCTION__, "The maturity period cannot be zero length");

    // Set effective roll day (if applicable) and do some more input validation to ensure
    // consistency of inputs.
    setEffectiveRollDay();

    if (effectiveRollDay == 0 && !backStubStartDate.empty()) {
        // Say we're generating IMM dates, and the user specifies a backStubStartDate that is not
        // an IMM date.  We disallow this because it could produce two back stubs!
        DateTime date(MaturityPeriod::toDate(-periodCount, periodInterval, backStubStartDate));
        if (backStubStartDate != MaturityPeriod::toDate(periodCount, periodInterval, date)) {
            string msg("Invalid explicit backStubStartDate: ");
            msg += backStubStartDate.toString() + " (not on cycle)";
            throw ModelException(__FUNCTION__, msg);
        }
    }
}

IObject* DateSchedGen::clone() const {
    if (getRefCount() == 0)
        return new DateSchedGen(*this);
    else
        return const_cast<DateSchedGen*>(this);
}

DateSchedGen::DateSchedGen(
                    const DateTime& startDate,
                    const DateTime& endDate,
                    const MaturityPeriodConstSP& maturityPeriod,
                    const DateTime& frontStubEndDate,
                    const DateTime& backStubStartDate,
                    bool includeStart,
                    bool includeEnd,
                    bool rollOnLast,
                    bool ignoreLeapYears,
                    int rollDay) :
    IDateSched(TYPE), startDate(startDate), endDate(endDate), maturityPeriod(maturityPeriod),
    frontStubEndDate(frontStubEndDate), backStubStartDate(backStubStartDate),
    includeStart(includeStart), includeEnd(includeEnd), rollOnLast(rollOnLast),
    ignoreLeapYears(ignoreLeapYears), rollDay(rollDay), periodCount(0),
    effectiveRollDay(0), useExplicitRollDay(false) {
        if (!maturityPeriod)
            throw ModelException(__FUNCTION__, "maturityPeriod cannot be NULL");
        validatePop2Object();
}

DateSchedGen::DateSchedGen(const CClassConstSP& type) :
    IDateSched(type), maturityPeriod(MaturityPeriodConstSP(0)), includeStart(true),
    includeEnd(true), rollOnLast(false), ignoreLeapYears(false), rollDay(0), periodCount(0),
    effectiveRollDay(0), useExplicitRollDay(false) {}

IObject* DateSchedGen::defaultConstructor() { return new DateSchedGen(TYPE); }

void DateSchedGen::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(DateSchedGen, clazz);
    clazz->setDescription("A flexible date schedule generator");
    SUPERCLASS(IDateSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(startDate, "Earliest date to potentially include in output");
    FIELD(endDate, "Latest date to potentially include in output");
    FIELD(maturityPeriod, "Interval between periodic dates");
    FIELD(frontStubEndDate, "Explicit front stub end (periodic start) date");
    FIELD_MAKE_OPTIONAL(frontStubEndDate);
    FIELD(backStubStartDate, "Explicit back stub start (periodic end) date");
    FIELD_MAKE_OPTIONAL(backStubStartDate);
    FIELD(includeStart, "Specifies whether startDate is included in the output (default = true)");
    FIELD_MAKE_OPTIONAL(includeStart);
    FIELD(includeEnd, "Specifies whether endDate is included in the output (default = true)");
    FIELD_MAKE_OPTIONAL(includeEnd);
    FIELD(rollOnLast, "Specifies whether to generic periodic dates that are always the last day of the month (default = false)");
    FIELD_MAKE_OPTIONAL(rollOnLast);
    FIELD(ignoreLeapYears, "Specifies whether February should always be considered to have just 28 days");
    FIELD_MAKE_OPTIONAL(ignoreLeapYears);
    FIELD(rollDay, "Specifies an explicit day of the month for dates that are generated with a month-based MaturityPeriod");
    FIELD_MAKE_OPTIONAL(rollDay);

    FIELD(periodCount, "");
    FIELD_MAKE_TRANSIENT(periodCount);
    FIELD(periodInterval, "");
    FIELD_MAKE_TRANSIENT(periodInterval);
    FIELD(effectiveRollDay, "");
    FIELD_MAKE_TRANSIENT(effectiveRollDay);
    FIELD(useExplicitRollDay, "");
    FIELD_MAKE_TRANSIENT(useExplicitRollDay);
}

CClassConstSP const DateSchedGen::TYPE = CClass::registerClassLoadMethod(
    "DateSchedGen", typeid(DateSchedGen), DateSchedGen::load);


///////////////////////////////////////////////////////////////////////////////
// AdjustedDateSched - Answers a new date schedule that is the result of
//                     applying a specified date adjustment to each of the
//                     dates of a source date schedule.
///////////////////////////////////////////////////////////////////////////////
void AdjustedDateSched::setup() const {
    const DateTimeArrayConstSP& sourceDatesSP = dateSched->elements();
    const DateTimeArray& sourceDates = *sourceDatesSP;
    DateTimeArraySP adjustedSchedSP(new DateTimeArray);
    DateTimeArray& adjustedSched = *adjustedSchedSP;

    adjustedSched.reserve(sourceDates.size());

    const int maxIdx = adjustLast ? sourceDates.size() : sourceDates.size() - 1;
    DateTime adjustedDate;
    for (int i = adjustFirst ? 0 : 1; i < maxIdx; ++i) {
        adjustedDate = dateAdjustment->dateFor(sourceDates[i]);
        if (includeDuplicates || adjustedSched.empty() || adjustedDate != adjustedSched.back())
            adjustedSched.push_back(adjustedDate);
    }
    
    cachedElemsSP = adjustedSchedSP;
}

IObject* AdjustedDateSched::clone() const {
    if (getRefCount() == 0)
        return new AdjustedDateSched(*this);
    else
        return const_cast<AdjustedDateSched*>(this);
}

AdjustedDateSched::AdjustedDateSched(const IDateSchedConstSP& dateSched,
                                     const IDateAdjustmentConstSP& dateAdjustment,
                                     bool adjustFirst,
                                     bool adjustLast,
                                     bool includeDuplicates) :
    IDateSched(TYPE), dateSched(dateSched), dateAdjustment(dateAdjustment),
    adjustFirst(adjustFirst), adjustLast(adjustLast),
    includeDuplicates(includeDuplicates) {
        if (!dateSched)
            throw ModelException(__FUNCTION__, "Date schedule cannot be NULL");
        if (!dateAdjustment)
            throw ModelException(__FUNCTION__, "Date adjustment cannot be NULL");
}

AdjustedDateSched::AdjustedDateSched(const CClassConstSP& type) :
    IDateSched(type), dateSched(IDateSchedConstSP(0)),
    dateAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)),
    adjustFirst(true), adjustLast(true), includeDuplicates(false) {}

IObject* AdjustedDateSched::defaultConstructor() { return new AdjustedDateSched(TYPE); }

void AdjustedDateSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(AdjustedDateSched, clazz);
    clazz->setDescription("Produces an array of dates by applying one or more date "
                          "adjustments to a given array of dates");
    SUPERCLASS(IDateSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(dateSched, "Source date schedule for producing adjusted schedule");
    FIELD(dateAdjustment, "DateAdjustment(s) to apply to source date schedule");
    FIELD(adjustFirst, "Whether to apply date adjustment to the first date");
    FIELD_MAKE_OPTIONAL(adjustFirst);
    FIELD(adjustLast, "Whether to apply date adjustment to the last date");
    FIELD_MAKE_OPTIONAL(adjustLast);
    FIELD(includeDuplicates, "Whether to include consecutive duplicate dates in the output");
    FIELD_MAKE_OPTIONAL(includeDuplicates);
}

CClassConstSP const AdjustedDateSched::TYPE = CClass::registerClassLoadMethod(
    "AdjustedDateSched", typeid(AdjustedDateSched), AdjustedDateSched::load);

        
///////////////////////////////////////////////////////////////////////////////
// DateSchedUnion - Answers a new date schedule that is the union of zero or
//                  more date schedules where the order of dates is preserved.
///////////////////////////////////////////////////////////////////////////////
void DateSchedUnion::setup() const {
    DateTimeArraySP mergedDateArray(new DateTimeArray);

    const IDateSchedArray& sourceScheds = *dateScheds;

    for (int i = 0; i < sourceScheds.size(); ++i) {
        DateTimeArrayConstSP datesToMerge(sourceScheds[i]->elements());
        DateTimeArraySP result(new DateTimeArray);
        insert_iterator<DateTimeArray> resultInserter(*result, result->begin());

        set_union(mergedDateArray->begin(), mergedDateArray->end(),
                  datesToMerge->begin(), datesToMerge->end(),
                  resultInserter);
        mergedDateArray = result;
    }

    cachedElemsSP = mergedDateArray;
}

IObject* DateSchedUnion::clone() const {
    if (getRefCount() == 0)
        return new DateSchedUnion(*this);
    else
        return const_cast<DateSchedUnion*>(this);
}

DateSchedUnion::DateSchedUnion(const IDateSchedArrayConstSP& dateScheds) :
    IDateSched(TYPE), dateScheds(dateScheds) {
        if (!dateScheds)
            throw ModelException(__FUNCTION__, "Date schedules array cannot be NULL");
}

DateSchedUnion::DateSchedUnion(const CClassConstSP& type) :
    IDateSched(type), dateScheds(IDateSchedArrayConstSP(0)) {}

IObject* DateSchedUnion::defaultConstructor() { return new DateSchedUnion(TYPE); }

void DateSchedUnion::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(DateSchedUnion, clazz);
    clazz->setDescription("Produces an array of dates that is the union of "
                          "zero or more date schedules with the order of "
                          "dates preserved");
    SUPERCLASS(IDateSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(dateScheds, "Source date schedules for producing merged date schedule");
}

CClassConstSP const DateSchedUnion::TYPE = CClass::registerClassLoadMethod(
    "DateSchedUnion", typeid(DateSchedUnion), DateSchedUnion::load);


///////////////////////////////////////////////////////////////////////////////
// Period - Represents an explicitly specified date/time interval.  The isStub
//          attribute can be used to mark periods that are of irregular length.
///////////////////////////////////////////////////////////////////////////////
bool Period::operator ==(const Period& period) const {
    return &period == this || (period.start == start && period.end == end && period.isStub == isStub);
}

bool Period::operator !=(const Period& period) const {
    return !(*this == period);
}

string Period::toString() const {
    string s("(start=");
    s += start.toString() + ", end=" + end.toString() + ", isStub=" + Format::toString(isStub) + ")";
    return s;
}

Period::Period(const CClassConstSP& type) : CObject(type), isStub(false) {}

void Period::validatePop2Object() {
    validatePeriodDates(start, end);
}

IObject* Period::defaultConstructor() { return new Period(TYPE); }

void Period::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(Period, clazz);
    clazz->setDescription("A time period specified by explicit dates");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(start, "Period start date");
    FIELD(end, "Period end date");
    FIELD(isStub, "Indicates whether this is a stub (i.e. irregular length) period");
    FIELD_MAKE_OPTIONAL(isStub);
}

CClassConstSP const Period::TYPE = CClass::registerClassLoadMethod(
    "Period", typeid(Period), Period::load);

DEFINE_TEMPLATE_TYPE(PeriodArray);

START_PUBLIC_ENUM_DEFINITION(Period::Boundary, "Specifies a period boundary date");
ENUM_VALUE_AND_NAME(Period::None, "None", "Indicates that no boundary date is specified");
ENUM_VALUE_AND_NAME(Period::Start, "PeriodStart", "Indicates the period start");
ENUM_VALUE_AND_NAME(Period::End, "PeriodEnd", "Indicates the period end");
END_ENUM_DEFINITION(Period::Boundary);


///////////////////////////////////////////////////////////////////////////////
// IPeriodSched - Interface for classes that answer a schedule of periods (i.e.
//                a PeriodArray.
///////////////////////////////////////////////////////////////////////////////
PeriodArrayConstSP IPeriodSched::elements() const {
    if (!cachedElemsSP)
        setup();
    return cachedElemsSP;
}

ISchedConstSP IPeriodSched::schedule() const { return ISchedConstSP(new PeriodSched(elements())); }

IPeriodSched::IPeriodSched(const CClassConstSP& type) :
    ISched(type), cachedElemsSP(PeriodArrayConstSP(0)) {}

void IPeriodSched::load(CClassSP& clazz) {
    clazz->setPublic();     // In case people want to create arrays of these.
    REGISTER(IPeriodSched, clazz);
    clazz->setDescription("Interface for a hierarchy of classes that represent a "
                          "schedule of Periods");
    SUPERCLASS(ISched);
    FIELD(cachedElemsSP, "");
    FIELD_MAKE_TRANSIENT(cachedElemsSP);
}

CClassConstSP const IPeriodSched::TYPE = CClass::registerClassLoadMethod(
    "IPeriodSched", typeid(IPeriodSched), IPeriodSched::load);


///////////////////////////////////////////////////////////////////////////////
// DateSchedFromPeriodSched - Answers a new date schedule that is the result of
//                            selecting either the period start or end dates
//                            from a source period schedule.  Optionally, one
//                            may specify a date adjustment to apply to the
//                            selected dates.
///////////////////////////////////////////////////////////////////////////////
void DateSchedFromPeriodSched::setup() const {
    const PeriodArrayConstSP& sourcePeriodsSP = periodSched->elements();
    const PeriodArray& sourcePeriods = *sourcePeriodsSP;
    DateTimeArraySP datesSP(new DateTimeArray(sourcePeriods.size()));
    DateTimeArray& dates = *datesSP;

    for (int i = 0; i < sourcePeriods.size(); ++i)
        if (startOrEnd == Period::Start)
            dates[i] = dateAdjustment->dateFor(sourcePeriods[i]->start);
        else
            dates[i] = dateAdjustment->dateFor(sourcePeriods[i]->end);

    cachedElemsSP = datesSP;
}

void DateSchedFromPeriodSched::validatePop2Object() {
    if (startOrEnd == Period::None)
            throw ModelException(__FUNCTION__, "startOrEnd cannot be None");
}

IObject* DateSchedFromPeriodSched::clone() const {
    if (getRefCount() == 0)
        return new DateSchedFromPeriodSched(*this);
    else
        return const_cast<DateSchedFromPeriodSched*>(this);
}

DateSchedFromPeriodSched::DateSchedFromPeriodSched(const IPeriodSchedConstSP& periodSched,
                                                   Period::Boundary startOrEnd,
                                                   const IDateAdjustmentConstSP& dateAdjustment) :
    IDateSched(TYPE), periodSched(periodSched), startOrEnd(startOrEnd), dateAdjustment(dateAdjustment) {
        if (!periodSched)
            throw ModelException(__FUNCTION__, "Period schedule cannot be NULL");
        if (!dateAdjustment)
            throw ModelException(__FUNCTION__, "Date adjustment cannot be NULL");
        validatePop2Object();
}

DateSchedFromPeriodSched::DateSchedFromPeriodSched(const CClassConstSP& type) :
    IDateSched(type), periodSched(IPeriodSchedConstSP(0)), startOrEnd(Period::Start),
    dateAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)) {}

IObject* DateSchedFromPeriodSched::defaultConstructor() { return new DateSchedFromPeriodSched(TYPE); }

void DateSchedFromPeriodSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(DateSchedFromPeriodSched, clazz);
    clazz->setDescription("Creates a date schedule from either the ("
                          "optionally adjusted) start or end dates of a "
                          "period schedule");
    SUPERCLASS(IDateSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(periodSched, "The source period schedule");
    FIELD(startOrEnd, "Specifies whether to select period start or end dates");
    FIELD(dateAdjustment, "Date adjustment(s) to apply to the selected dates");
    FIELD_MAKE_OPTIONAL(dateAdjustment);
}

CClassConstSP const DateSchedFromPeriodSched::TYPE = CClass::registerClassLoadMethod(
        "DateSchedFromPeriodSched", typeid(DateSchedFromPeriodSched), DateSchedFromPeriodSched::load);


///////////////////////////////////////////////////////////////////////////////
// PeriodSched - An explicitly specified period schedule.
///////////////////////////////////////////////////////////////////////////////
void PeriodSched::setup() const { cachedElemsSP = explicitPeriods; }

IObject* PeriodSched::clone() const {
    if (getRefCount() == 0)
        return new PeriodSched(*this);
    else
        return const_cast<PeriodSched*>(this);
}

PeriodSched::PeriodSched(const PeriodArrayConstSP& explicitPeriods) :
    IPeriodSched(TYPE), explicitPeriods(explicitPeriods) {
        if (!explicitPeriods)
            throw ModelException(__FUNCTION__, "Explicit PeriodArray cannot be NULL");
}

PeriodSched::PeriodSched(const CClassConstSP& type) :
    IPeriodSched(type), explicitPeriods(PeriodArrayConstSP(0)) {}

void PeriodSched::validatePop2Object() {
    const PeriodArray& periodArray = *explicitPeriods;
    for (int i = 1; i < periodArray.size(); ++i) {
        const Period& cur = *(periodArray[i - 1]);
        const Period& nxt = *(periodArray[i]);
        ensureStrictlyIncreasing(cur.start, nxt.start, "Period start");
        ensureStrictlyIncreasing(cur.end, nxt.end, "Period end");
    }
}

IObject* PeriodSched::defaultConstructor() { return new PeriodSched(TYPE); }

void PeriodSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(PeriodSched, clazz);
    clazz->setDescription("An explicitly specified array of Periods");
    SUPERCLASS(IPeriodSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitPeriods, periods, "The array of Periods");
}

CClassConstSP const PeriodSched::TYPE = CClass::registerClassLoadMethod(
    "PeriodSched", typeid(PeriodSched), PeriodSched::load);


///////////////////////////////////////////////////////////////////////////////
// PeriodSchedGen - Generates a schedule of periods from the consective dates
//                  that are produced by the specified DateSchedGen.  The first
//                  and last periods may be marked as stubs if the first period
//                  start date and last period end date are not periodic dates,
//                  respectively.
///////////////////////////////////////////////////////////////////////////////
void PeriodSchedGen::setup() const {
    PeriodArraySP periodListSP(new PeriodArray);
    PeriodArray& periodList = *periodListSP;

    const DateTimeArrayConstSP& dateListSP = dateSchedGen->elements();
    const DateTimeArray& dateList = *dateListSP;

    const int lastPeriod = dateList.size() - 1;
    for (int i = 1; i <= lastPeriod; ++i)
    {
        const DateTime& start = dateList[i - 1];
        const DateTime& end = dateList[i];
        bool isStub;
        if (i == 1 || i == lastPeriod)
            // We only need to check the first and last periods since the DateSchedGen
            // ensures that at most one front and/or back stub will exist, and we don't
            // want to waste time checking intermediate periods for no reason.
            isStub = dateSchedGen->isStubPeriod(start, end);
        else
            isStub = false;
        periodList.push_back(PeriodSP(new Period(start, end, isStub)));
    }

    cachedElemsSP = periodListSP;
}

IObject* PeriodSchedGen::clone() const {
    if (getRefCount() == 0)
        return new PeriodSchedGen(*this);
    else
        return const_cast<PeriodSchedGen*>(this);
}

PeriodSchedGen::PeriodSchedGen(const DateSchedGenConstSP& dateSchedGen) :
    IPeriodSched(TYPE), dateSchedGen(dateSchedGen) {
        if (!dateSchedGen)
            throw ModelException(__FUNCTION__, "Date schedule cannot be NULL");
}

PeriodSchedGen::PeriodSchedGen(const CClassConstSP& type) :
    IPeriodSched(type), dateSchedGen(DateSchedGenConstSP(0)) {}

IObject* PeriodSchedGen::defaultConstructor() { return new PeriodSchedGen(TYPE); }

void PeriodSchedGen::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(PeriodSchedGen, clazz);
    clazz->setDescription("A flexible Period schedule generator");
    SUPERCLASS(IPeriodSched);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(dateSchedGen, "A date schedule generator that will be used to produce an array of contiguous Periods");
}

CClassConstSP const PeriodSchedGen::TYPE = CClass::registerClassLoadMethod(
    "PeriodSchedGen", typeid(PeriodSchedGen), PeriodSchedGen::load);


///////////////////////////////////////////////////////////////////////////////
// AdjustedPeriodSched - Generates a schedule of periods that is the result of
//                       applying a date adjustment to each period's start and
//                       end date.  One can independently specify whether to
//                       adjust the first period start date and the last period
//                       end date.
///////////////////////////////////////////////////////////////////////////////
void AdjustedPeriodSched::setup() const {
    const PeriodArrayConstSP& sourcePeriods = periodSched->elements();
    PeriodArraySP adjustedSchedSP(new PeriodArray(*sourcePeriods));
    PeriodArray& adjustedSched = *adjustedSchedSP;

    bool adjustStart, adjustEnd;
    const int maxIdx = adjustedSched.size() - 1;
    for (int i = 0; i <= maxIdx; ++i)
    { 
        adjustStart = adjustEnd = true;
        if (i == 0)
            adjustStart = adjustFirst;
        else if (i == maxIdx)
            adjustEnd = adjustLast;
            
        Period& period = *(adjustedSched[i]);
        if (adjustStart)
            period.start = dateAdjustment->dateFor(period.start);
        if (adjustEnd)
            period.end = dateAdjustment->dateFor(period.end);
    }

    cachedElemsSP = adjustedSchedSP;
}

IObject* AdjustedPeriodSched::clone() const {
    if (getRefCount() == 0)
        return new AdjustedPeriodSched(*this);
    else
        return const_cast<AdjustedPeriodSched*>(this);
}

AdjustedPeriodSched::AdjustedPeriodSched(const IPeriodSchedConstSP& periodSched,
                                         const IDateAdjustmentConstSP& dateAdjustment,
                                         bool adjustFirst,
                                         bool adjustLast) :
    IPeriodSched(TYPE), periodSched(periodSched), dateAdjustment(dateAdjustment),
    adjustFirst(adjustFirst), adjustLast(adjustLast) {
        if (!periodSched)
            throw ModelException(__FUNCTION__, "Period schedule cannot be NULL");
        if (!dateAdjustment)
            throw ModelException(__FUNCTION__, "Date adjustment cannot be NULL");
}

AdjustedPeriodSched::AdjustedPeriodSched(const CClassConstSP& type) :
    IPeriodSched(type), periodSched(IPeriodSchedConstSP(0)),
    dateAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)),
    adjustFirst(true), adjustLast(true) {}

IObject* AdjustedPeriodSched::defaultConstructor() { return new AdjustedPeriodSched(TYPE); }

void AdjustedPeriodSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(AdjustedPeriodSched, clazz);
    clazz->setDescription("Produces an array of Periods whose dates will be adjusted by "
                          "applying one or more date adjustments to each Period's start "
                          "and end dates");
    SUPERCLASS(IPeriodSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(periodSched, "Source Period schedule for producing adjusted schedule");
    FIELD(dateAdjustment, "DateAdjustment(s) to apply to source Period schedule");
    FIELD(adjustFirst, "Whether to apply date adjustment to the first Period start date");
    FIELD_MAKE_OPTIONAL(adjustFirst);
    FIELD(adjustLast, "Whether to apply date adjustment to the last Period end date");
    FIELD_MAKE_OPTIONAL(adjustLast);
}

CClassConstSP const AdjustedPeriodSched::TYPE = CClass::registerClassLoadMethod(
    "AdjustedPeriodSched", typeid(AdjustedPeriodSched), AdjustedPeriodSched::load);


///////////////////////////////////////////////////////////////////////////////
// PeriodObsDates - Associates an array of (observation) dates with a Period.
///////////////////////////////////////////////////////////////////////////////
PeriodObsDates::PeriodObsDates(const PeriodConstSP& period, const DateTimeArrayConstSP& obsDates) :
    CObject(TYPE), period(period), obsDates(obsDates) {
    if (!period)
        throw ModelException(__FUNCTION__, "Period cannot be NULL");
    if (!obsDates)
        throw ModelException(__FUNCTION__, "Observation dates array cannot be NULL");
}

PeriodObsDates::PeriodObsDates(const PeriodConstSP& period, const DateTime& obsDate) :
    CObject(TYPE), period(period) {
    if (!period)
        throw ModelException(__FUNCTION__, "Period cannot be NULL");
    DateTimeArraySP anObsDatesSP(new DateTimeArray(1));
    (*anObsDatesSP)[0] = obsDate;
    obsDates = anObsDatesSP;
}

PeriodObsDates::PeriodObsDates(const CClassConstSP& type) : CObject(type), period(PeriodConstSP(0)),
    obsDates(DateTimeArrayConstSP(0)) {}

void PeriodObsDates::validatePop2Object() {
    validatePeriodDates(period->start, period->end);
    DateTime::ensureIncreasing(*obsDates, "period observation dates", false);
}

IObject* PeriodObsDates::defaultConstructor() { return new PeriodObsDates(TYPE); }

void PeriodObsDates::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(PeriodObsDates, clazz);
    clazz->setDescription("An object that gives observation dates associated with a particular period");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(period, "The period for the associated observation dates");
    FIELD(obsDates, "The observation dates for the associated period");
}

CClassConstSP const PeriodObsDates::TYPE = CClass::registerClassLoadMethod(
    "PeriodObsDates", typeid(PeriodObsDates), PeriodObsDates::load);

DEFINE_TEMPLATE_TYPE(PeriodObsDatesArray);


///////////////////////////////////////////////////////////////////////////////
// FlatPeriodObsDates - A "flattened" version of PeriodObsDates.
///////////////////////////////////////////////////////////////////////////////
FlatPeriodObsDate::FlatPeriodObsDate(const DateTime& periodStart, const DateTime& periodEnd,
                                     bool isStub, const DateTime& obsDate) :
    CObject(TYPE), periodStart(periodStart), periodEnd(periodEnd), isStub(isStub), obsDate(obsDate) {}

FlatPeriodObsDate::FlatPeriodObsDate(const CClassConstSP& type) : CObject(type), isStub(false) {}

void FlatPeriodObsDate::validatePop2Object() {
    validatePeriodDates(periodStart, periodEnd);
}

IObject* FlatPeriodObsDate::defaultConstructor() { return new FlatPeriodObsDate(TYPE); }

void FlatPeriodObsDate::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FlatPeriodObsDate, clazz);
    clazz->setDescription("An object that gives an observation date associated with a particular period, "
                          "which is represented in a flattened form");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(periodStart, "The period start date for the associated observation date");
    FIELD(periodEnd, "The period end date for the associated observation date");
    FIELD(isStub, "Whether the period is an irregular length period");
    FIELD(obsDate, "An observation date for the associated period");
}

CClassConstSP const FlatPeriodObsDate::TYPE = CClass::registerClassLoadMethod(
    "FlatPeriodObsDate", typeid(FlatPeriodObsDate), FlatPeriodObsDate::load);

DEFINE_TEMPLATE_TYPE(FlatPeriodObsDateArray);


///////////////////////////////////////////////////////////////////////////////
// IObsSched - Interface for classes that represent a schedule of
//             PeriodObsDates.
///////////////////////////////////////////////////////////////////////////////
PeriodObsDatesArrayConstSP IObsSched::elements() const {
    if (!cachedElemsSP)
        setup();
    return cachedElemsSP;
}

ISchedConstSP IObsSched::schedule() const { return ISchedConstSP(new ObsSched(elements())); }

ISchedConstSP IObsSched::flatSchedule() const {
    const PeriodObsDatesArrayConstSP& periodObsDatesArraySP = elements();
    const PeriodObsDatesArray& periodObsDatesArray = *periodObsDatesArraySP;

    FlatPeriodObsDateArraySP flatPeriodObsDatesSP(new FlatPeriodObsDateArray);
    FlatPeriodObsDateArray& flatPeriodObsDates = *flatPeriodObsDatesSP;

    for (int i = 0; i < periodObsDatesArray.size(); ++i) {
        const Period& period = *(periodObsDatesArray[i]->period);
        const DateTimeArray& obsDates = *(periodObsDatesArray[i]->obsDates);

        for (int j = 0; j < obsDates.size(); ++j)
            flatPeriodObsDates.push_back(
                FlatPeriodObsDateSP(new FlatPeriodObsDate(period.start, period.end, period.isStub, obsDates[j])));
    }

    return ISchedConstSP(new FlatObsSched(flatPeriodObsDatesSP));
}

IObsSched::IObsSched(const CClassConstSP& type) :
    ISched(type), cachedElemsSP(PeriodObsDatesArrayConstSP(0)) {}

void IObsSched::load(CClassSP& clazz) {
    clazz->setPublic();     // In case people want to create arrays of these.
    REGISTER(IObsSched, clazz);
    clazz->setDescription("Interface for a hierarchy of classes that represent a "
                          "schedule of observation dates for an associated period");
    SUPERCLASS(ISched);

    FIELD(cachedElemsSP, "");
    FIELD_MAKE_TRANSIENT(cachedElemsSP);
}

CClassConstSP const IObsSched::TYPE = CClass::registerClassLoadMethod(
    "IObsSched", typeid(IObsSched), IObsSched::load);


///////////////////////////////////////////////////////////////////////////////
// ObsSched - An explicitly specified schedule of PeriodObsDates.
///////////////////////////////////////////////////////////////////////////////
void ObsSched::setup() const { cachedElemsSP = explicitPeriodObsDates; }

IObject* ObsSched::clone() const {
    if (getRefCount() == 0)
        return new ObsSched(*this);
    else
        return const_cast<ObsSched*>(this);
}

ObsSched::ObsSched(const PeriodObsDatesArrayConstSP& explicitPeriodObsDates) :
    IObsSched(TYPE), explicitPeriodObsDates(explicitPeriodObsDates) {
        if (!explicitPeriodObsDates)
            throw ModelException(__FUNCTION__, "Explicit period observation "
                                               "dates array cannot be NULL");
}

ObsSched::ObsSched(const CClassConstSP& type) :
    IObsSched(type), explicitPeriodObsDates(PeriodObsDatesArrayConstSP(0)) {}

void ObsSched::validatePop2Object() {
    const PeriodObsDatesArray& periodObsDatesArray = *explicitPeriodObsDates;
    for (int i = 1; i < periodObsDatesArray.size(); ++i) {
        const PeriodObsDates& cur = *(periodObsDatesArray[i - 1]);
        const PeriodObsDates& nxt = *(periodObsDatesArray[i]);
        ensureStrictlyIncreasing(cur.period->start, nxt.period->start, "Period start");
        ensureStrictlyIncreasing(cur.period->end, nxt.period->end, "Period end");
    }
}

IObject* ObsSched::defaultConstructor() { return new ObsSched(TYPE); }

void ObsSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(ObsSched, clazz);
    clazz->setDescription("An explicitly specified array of PeriodObsDates");
    SUPERCLASS(IObsSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitPeriodObsDates, periodObsDates, "The array of PeriodObsDates");
}

CClassConstSP const ObsSched::TYPE = CClass::registerClassLoadMethod(
    "ObsSched", typeid(ObsSched), ObsSched::load);


///////////////////////////////////////////////////////////////////////////////
// FlatObsSched - An explicitly specified schedule of FlatPeriodObsDates.
///////////////////////////////////////////////////////////////////////////////
ISchedConstSP FlatObsSched::flatSchedule() const { return ISchedConstSP(this); }

void FlatObsSched::setup() const {
    PeriodObsDatesArraySP periodObsDatesArraySP(new PeriodObsDatesArray);
    PeriodObsDatesArray& periodObsDatesArray = *periodObsDatesArraySP;

    const FlatPeriodObsDateArray& flatPeriodObsDates = *explicitFlatPeriodObsDates;

    PeriodSP periodSP(0);
    DateTimeArraySP obsDatesSP(0);
    DateTime periodStart, periodEnd;
    bool isStub;

    for (int i = 0; i < flatPeriodObsDates.size(); ++i)
    {
        const FlatPeriodObsDate& flatPeriodObsDate = *(flatPeriodObsDates[i]);
        if (i == 0 ||
            flatPeriodObsDate.periodStart != periodStart ||
            flatPeriodObsDate.periodEnd != periodEnd)
        {
            if (!!periodSP)
                periodObsDatesArray.push_back(PeriodObsDatesSP(new PeriodObsDates(periodSP, obsDatesSP)));
            periodStart = flatPeriodObsDate.periodStart;
            periodEnd = flatPeriodObsDate.periodEnd;
            isStub = flatPeriodObsDate.isStub;
            periodSP.reset(new Period(periodStart, periodEnd, isStub));
            obsDatesSP.reset(new DateTimeArray);
        }
        else if (flatPeriodObsDate.isStub != isStub) {
            string msg("Encountered inconsistent stub specifications for period (start=");
            msg += periodStart.toString() + ", end=" + periodEnd.toString() + ")";
            throw ModelException(__FUNCTION__, msg);
        }

        obsDatesSP->push_back(flatPeriodObsDate.obsDate);
    }

    if (!!periodSP)
        periodObsDatesArray.push_back(PeriodObsDatesSP(new PeriodObsDates(periodSP, obsDatesSP)));

    cachedElemsSP = periodObsDatesArraySP;
}

IObject* FlatObsSched::clone() const {
    if (getRefCount() == 0)
        return new FlatObsSched(*this);
    else
        return const_cast<FlatObsSched*>(this);
}

FlatObsSched::FlatObsSched(const FlatPeriodObsDateArrayConstSP& explicitFlatPeriodObsDates) :
    IObsSched(TYPE), explicitFlatPeriodObsDates(explicitFlatPeriodObsDates) {
        if (!explicitFlatPeriodObsDates)
            throw ModelException(__FUNCTION__, "Explicit flat period observation "
                                               "dates array cannot be NULL");
}

FlatObsSched::FlatObsSched(const CClassConstSP& type) :
    IObsSched(type), explicitFlatPeriodObsDates(FlatPeriodObsDateArrayConstSP(0)) {}

void FlatObsSched::validatePop2Object() {
    const FlatPeriodObsDateArray& flatPeriodObsDateArray = *explicitFlatPeriodObsDates;
    for (int i = 1; i < flatPeriodObsDateArray.size(); ++i) {
        const FlatPeriodObsDate& cur = *(flatPeriodObsDateArray[i - 1]);
        const FlatPeriodObsDate& nxt = *(flatPeriodObsDateArray[i]);

        ensureIncreasing(cur.periodStart, nxt.periodStart, "Period start");
        ensureIncreasing(cur.periodEnd, nxt.periodEnd, "Period end");

        if (cur.periodStart == nxt.periodStart && cur.periodEnd == nxt.periodEnd)
            ensureIncreasing(cur.obsDate, nxt.obsDate, "Period observation");
    }
}

IObject* FlatObsSched::defaultConstructor() { return new FlatObsSched(TYPE); }

void FlatObsSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FlatObsSched, clazz);
    clazz->setDescription("An explicitly specified array of FlatPeriodObsDate");
    SUPERCLASS(IObsSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitFlatPeriodObsDates, flatPeriodObsDates, "The array of FlatPeriodObsDate");
}

CClassConstSP const FlatObsSched::TYPE = CClass::registerClassLoadMethod(
    "FlatObsSched", typeid(FlatObsSched), FlatObsSched::load);


///////////////////////////////////////////////////////////////////////////////
// RelativeObsSched - Generates a schedule of PeriodObsDates where there is one
//                    observation date for each period, and each observation
//                    date is either the period start or end date plus an
//                    optional date adjustment.
///////////////////////////////////////////////////////////////////////////////
void RelativeObsSched::setup() const {
    const PeriodArrayConstSP& sourcePeriodsSP = periodSched->elements();
    const PeriodArray& sourcePeriods = *sourcePeriodsSP;
    PeriodObsDatesArraySP periodObsDatesArraySP(new PeriodObsDatesArray(sourcePeriods.size()));
    PeriodObsDatesArray& periodObsDatesArray = *periodObsDatesArraySP;
    
    DateSchedFromPeriodSched obsSched(periodSched, relativeToDate, obsAdjustment);
    const DateTimeArrayConstSP& obsDatesSP = obsSched.elements();
    const DateTimeArray& obsDates = *obsDatesSP;

    for (int i = 0; i < sourcePeriods.size(); ++i)
        periodObsDatesArray[i] = PeriodObsDatesSP(new PeriodObsDates(sourcePeriods[i], obsDates[i]));

    cachedElemsSP = periodObsDatesArraySP;
}

IObject* RelativeObsSched::clone() const {
    if (getRefCount() == 0)
        return new RelativeObsSched(*this);
    else
        return const_cast<RelativeObsSched*>(this);
}

void RelativeObsSched::validatePop2Object() {
    if (relativeToDate == Period::None)
            throw ModelException(__FUNCTION__, "relativeToDate cannot be None");
}

RelativeObsSched::RelativeObsSched(const IPeriodSchedConstSP& periodSched,
                                   Period::Boundary relativeToDate,
                                   const IDateAdjustmentConstSP& obsAdjustment) :
    IObsSched(TYPE), periodSched(periodSched), relativeToDate(relativeToDate), obsAdjustment(obsAdjustment) {
        if (!periodSched)
            throw ModelException(__FUNCTION__, "Period schedule cannot be NULL");
        if (!obsAdjustment)
            throw ModelException(__FUNCTION__, "Observation date adjustment cannot be NULL");
        validatePop2Object();
}

RelativeObsSched::RelativeObsSched(const CClassConstSP& type) :
    IObsSched(type), periodSched(IPeriodSchedConstSP(0)), relativeToDate(Period::Start),
    obsAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)) {}

IObject* RelativeObsSched::defaultConstructor() { return new RelativeObsSched(TYPE); }

void RelativeObsSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(RelativeObsSched, clazz);
    clazz->setDescription("Produces an array of PeriodObsDates where the observation date for "
                          "each PeriodObsDates is given by some date adjustment (e.g. offset) "
                          "that is applied to each period start or end date of a given period schedule");
    SUPERCLASS(IObsSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(periodSched, "Source Period schedule for producing PeriodObsDatesArray");
    FIELD(relativeToDate, "Selects whether observation is relative to period start or end date");
    FIELD(obsAdjustment, "DateAdjustment(s) to apply to source period start or end date to produce an observation date");
    FIELD_MAKE_OPTIONAL(obsAdjustment);
}

CClassConstSP const RelativeObsSched::TYPE = CClass::registerClassLoadMethod(
    "RelativeObsSched", typeid(RelativeObsSched), RelativeObsSched::load);


///////////////////////////////////////////////////////////////////////////////
// WindowedObsSched - Generates a schedule of PeriodObsDates where there are
//                    zero or more observation dates associated with each
//                    period.  Observation dates are given by an independent
//                    date schedule and associated with each period via a
//                    rolling selection window that is specified relative to
//                    either or both the period start and end dates plus
//                    optional offsets.  A cutoff offset may be specified
//                    relative to either the period start or end date to
//                    indicate the date after which to stop selecting
//                    observations that otherwise fall within the selection
//                    window.  The repeatLastObs flag specifies whether to
//                    repeat the observation date that falls on the cutoff date
//                    for each of the observation dates that were skipped as a
//                    result of the cutoff.
///////////////////////////////////////////////////////////////////////////////
void WindowedObsSched::setup() const {
    const PeriodArrayConstSP& sourcePeriodsSP = periodSched->elements();
    const PeriodArray& sourcePeriods = *sourcePeriodsSP;
    const DateTimeArrayConstSP& obsDatesSP = obsSched->elements();
    const DateTimeArray& obsDates = *obsDatesSP;

    PeriodObsDatesArraySP periodObsDatesArraySP(new PeriodObsDatesArray(sourcePeriods.size()));
    PeriodObsDatesArray& periodObsDatesArray = *periodObsDatesArraySP;

    DateTime cutoffDate, windowStartDate, windowEndDate;
    int obsIdx = 0;

    for (int periodIdx = 0; periodIdx < sourcePeriods.size(); ++periodIdx) {
        DateTimeArraySP obsDateArraySP(new DateTimeArray);

        if (selectionStartRelativeToDate != Period::None) {
            const DateTime& selectionStartDate = selectionStartRelativeToDate == Period::End ?
                                                 sourcePeriods[periodIdx]->end : sourcePeriods[periodIdx]->start;
            windowStartDate = selectionStartOffset->dateFor(selectionStartDate);
        }
 
        const DateTime& selectionEndDate = selectionEndRelativeToDate == Period::Start ?
                                           sourcePeriods[periodIdx]->start : sourcePeriods[periodIdx]->end;
        windowEndDate = selectionEndOffset->dateFor(selectionEndDate);
         
        if (obsCutoffRelativeToDate != Period::None) {
            const DateTime& obsCutoffDate = obsCutoffRelativeToDate == Period::Start ?
                                            sourcePeriods[periodIdx]->start : sourcePeriods[periodIdx]->end;
            cutoffDate = min(obsCutoffOffset->dateFor(obsCutoffDate), windowEndDate);
        }
        else
            cutoffDate = windowEndDate;

        while (obsIdx < obsDates.size() && obsDates[obsIdx] < windowStartDate)
            ++obsIdx;

        while (obsIdx < obsDates.size() && obsDates[obsIdx] <= cutoffDate)
            obsDateArraySP->push_back(obsAdjustment->dateFor(obsDates[obsIdx++]));

        if (repeatLastObs && obsIdx > 0 && obsIdx < obsDates.size()) {
            DateTime repeatedObsDate(obsAdjustment->dateFor(obsDates[obsIdx - 1]));
            for (; obsIdx < obsDates.size() && obsDates[obsIdx] <= windowEndDate; ++obsIdx)
                obsDateArraySP->push_back(repeatedObsDate);
        }
        else
            for (; obsIdx < obsDates.size() && obsDates[obsIdx] <= windowEndDate; ++obsIdx);

        periodObsDatesArray[periodIdx] = PeriodObsDatesSP(new PeriodObsDates(sourcePeriods[periodIdx], obsDateArraySP));
    }

    cachedElemsSP = periodObsDatesArraySP;
}

IObject* WindowedObsSched::clone() const {
    if (getRefCount() == 0)
        return new WindowedObsSched(*this);
    else
        return const_cast<WindowedObsSched*>(this);
}

void WindowedObsSched::validatePop2Object() {
    if (selectionEndRelativeToDate == Period::None)
            throw ModelException(__FUNCTION__, "selectionEndRelativeToDate cannot be None");
}

WindowedObsSched::WindowedObsSched(const IPeriodSchedConstSP& periodSched,
                                   const IDateSchedConstSP& obsSched,
                                   const IDateAdjustmentConstSP& obsAdjustment,
                                   Period::Boundary selectionStartRelativeToDate,
                                   const IDateAdjustmentConstSP& selectionStartOffset,
                                   Period::Boundary selectionEndRelativeToDate,
                                   const IDateAdjustmentConstSP&selectionEndOffset,
                                   Period::Boundary obsCutoffRelativeToDate,
                                   const IDateAdjustmentConstSP& obsCutoffOffset,
                                   bool repeatLastObs) :
    IObsSched(TYPE), periodSched(periodSched), obsSched(obsSched), obsAdjustment(obsAdjustment),
    selectionStartRelativeToDate(selectionStartRelativeToDate), selectionStartOffset(selectionStartOffset),
    selectionEndRelativeToDate(selectionEndRelativeToDate), selectionEndOffset(selectionEndOffset),
    obsCutoffRelativeToDate(obsCutoffRelativeToDate), obsCutoffOffset(obsCutoffOffset), repeatLastObs(repeatLastObs) {
        if (!periodSched)
            throw ModelException(__FUNCTION__, "Period schedule cannot be NULL");
        if (!obsSched)
            throw ModelException(__FUNCTION__, "Observation date schedule cannot be NULL");
        if (!obsAdjustment)
            throw ModelException(__FUNCTION__, "Observation date adjustment cannot be NULL");
        if (!selectionStartOffset)
            throw ModelException(__FUNCTION__, "Selection window start offset cannot be NULL");
        if (!selectionEndOffset)
            throw ModelException(__FUNCTION__, "Selection window end offset cannot be NULL");
        if (!obsCutoffOffset)
            throw ModelException(__FUNCTION__, "Observation cutoff offset cannot be NULL");
        validatePop2Object();
}

WindowedObsSched::WindowedObsSched(const CClassConstSP& type) :
    IObsSched(type), periodSched(IPeriodSchedConstSP(0)), obsSched(IDateSchedConstSP(0)),
    obsAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)),
    selectionStartRelativeToDate(Period::None),
    selectionStartOffset(IDateAdjustmentConstSP(new NoDateAdjustment)),
    selectionEndRelativeToDate(Period::End),
    selectionEndOffset(IDateAdjustmentConstSP(new NoDateAdjustment)),
    obsCutoffRelativeToDate(Period::None),
    obsCutoffOffset(IDateAdjustmentConstSP(new NoDateAdjustment)), repeatLastObs(false) {}

IObject* WindowedObsSched::defaultConstructor() { return new WindowedObsSched(TYPE); }

void WindowedObsSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(WindowedObsSched, clazz);
    clazz->setDescription("Produces an array of PeriodObsDates where observation dates given by a "
                          "date schedule are assigned to each PeriodObsDates by a rolling selection "
                          "window that is specified relative to each period start and end date");
    SUPERCLASS(IObsSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(periodSched, "Source Period schedule for producing PeriodObsDatesArray");
    FIELD(obsSched, "Observation date schedule");
    FIELD(obsAdjustment, "Date adjustment to apply to selected observation dates");
    FIELD_MAKE_OPTIONAL(obsAdjustment);
    FIELD(selectionStartRelativeToDate, "Specifies the period date that is the base for the start of each selection window");
    FIELD_MAKE_OPTIONAL(selectionStartRelativeToDate);
    FIELD(selectionStartOffset, "Applied to the selectionStartRelativeToDate to produce the current selection window start date (inclusive)");
    FIELD_MAKE_OPTIONAL(selectionStartOffset);
    FIELD(selectionEndRelativeToDate, "Specifies the period date that is the base for the end of each selection window");
    FIELD_MAKE_OPTIONAL(selectionEndRelativeToDate);
    FIELD(selectionEndOffset, "Applied to the selectionEndRelativeToDate to produce the current selection window end date (inclusive)");
    FIELD_MAKE_OPTIONAL(selectionEndOffset);
    FIELD(obsCutoffRelativeToDate, "Specifies the period date that is the base for the cutoff date within each selection window");
    FIELD_MAKE_OPTIONAL(obsCutoffRelativeToDate);
    FIELD(obsCutoffOffset, "Applied to the obsCutoffRelativeToDate to determine the observation cutoff date (inclusive) "
                           "for the current selection window");
    FIELD_MAKE_OPTIONAL(obsCutoffOffset);
    FIELD(repeatLastObs, "Whether to repeat the observation on the cutoff date when the cutoff is before the "
                         "selection window end date");
    FIELD_MAKE_OPTIONAL(repeatLastObs);
}

CClassConstSP const WindowedObsSched::TYPE = CClass::registerClassLoadMethod(
    "WindowedObsSched", typeid(WindowedObsSched), WindowedObsSched::load);


///////////////////////////////////////////////////////////////////////////////
// OnOrBeforeObsSched - Generates a schedule of PeriodObsDates where there is
//                      zero or one observation date for each period.
//                      Observation dates are generated via an independent
//                      schedule and are associated with a period by selecting
//                      the last observation date that is on or before a
//                      rolling cutoff date.  The cutoff date may be either the
//                      period start or end date subject to an optional date
//                      adjustment.
///////////////////////////////////////////////////////////////////////////////
void OnOrBeforeObsSched::setup() const {
    const PeriodArrayConstSP& sourcePeriodsSP = periodSched->elements();
    const PeriodArray& sourcePeriods = *sourcePeriodsSP;
    const DateTimeArrayConstSP& obsDatesSP = obsSched->elements();
    const DateTimeArray& obsDates = *obsDatesSP;

    PeriodObsDatesArraySP periodObsDatesArraySP(new PeriodObsDatesArray(sourcePeriods.size()));
    PeriodObsDatesArray& periodObsDatesArray = *periodObsDatesArraySP;

    DateSchedFromPeriodSched obsCutoffSched(periodSched, cutoffRelativeToDate, cutoffOffset);
    const DateTimeArrayConstSP& obsCutoffDatesSP = obsCutoffSched.elements();
    const DateTimeArray& obsCutoffDates = *obsCutoffDatesSP;

    int obsIdx = 0;
    DateTime obsDate;
    for (int periodIdx = 0; periodIdx < sourcePeriods.size(); ++periodIdx) {
        DateTimeArraySP obsDateArraySP(new DateTimeArray);
        DateTime cutoffDate(obsCutoffDates[periodIdx]);

        // obsDates are expected to be monotonically increasing
        for (; obsIdx < obsDates.size() && obsDates[obsIdx] <= cutoffDate; ++obsIdx);

        if (obsIdx - 1 >= 0 && (obsDate = obsDates[obsIdx - 1]) <= cutoffDate)
            obsDateArraySP->push_back(obsDate);

        periodObsDatesArray[periodIdx] = PeriodObsDatesSP(
                                                new PeriodObsDates(sourcePeriods[periodIdx], obsDateArraySP));
    }

    cachedElemsSP = periodObsDatesArraySP;
}

IObject* OnOrBeforeObsSched::clone() const {
    if (getRefCount() == 0)
        return new OnOrBeforeObsSched(*this);
    else
        return const_cast<OnOrBeforeObsSched*>(this);
}

void OnOrBeforeObsSched::validatePop2Object() {
    if (cutoffRelativeToDate == Period::None)
            throw ModelException(__FUNCTION__, "cutoffRelativeToDate cannot be None");
}

OnOrBeforeObsSched::OnOrBeforeObsSched(const IPeriodSchedConstSP& periodSched,
                                       const IDateSchedConstSP& obsSched,
                                       Period::Boundary cutoffRelativeToDate,
                                       const IDateAdjustmentConstSP& cutoffOffset) :
    IObsSched(TYPE), periodSched(periodSched), obsSched(obsSched), cutoffRelativeToDate(cutoffRelativeToDate),
    cutoffOffset(cutoffOffset) {
        if (!periodSched)
            throw ModelException(__FUNCTION__, "Period schedule cannot be NULL");
        if (!obsSched)
            throw ModelException(__FUNCTION__, "Observation date schedule cannot be NULL");
        if (!cutoffOffset)
            throw ModelException(__FUNCTION__, "Observation cutoff offset cannot be NULL");
        validatePop2Object();
}

OnOrBeforeObsSched::OnOrBeforeObsSched(const CClassConstSP& type) :
    IObsSched(type), periodSched(IPeriodSchedConstSP(0)), obsSched(IDateSchedConstSP(0)),
    cutoffRelativeToDate(Period::Start), cutoffOffset(IDateAdjustmentConstSP(new NoDateAdjustment)) {}

IObject* OnOrBeforeObsSched::defaultConstructor() { return new OnOrBeforeObsSched(TYPE); }

void OnOrBeforeObsSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OnOrBeforeObsSched, clazz);
    clazz->setDescription("Produces an array of PeriodObsDates where the observation date for "
                          "each PeriodObsDates is the last observation date that occurs on or "
                          "before the cutoff date (i.e. period start or end), which may be "
                          "optionally adjusted by an offset");
    SUPERCLASS(IObsSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(periodSched, "Source Period schedule for producing PeriodObsDatesArray");
    FIELD(obsSched, "Source date schedule for observation dates");
    FIELD(cutoffRelativeToDate, "Specifies whether observation dates are selected relative to period start or end date");
    FIELD(cutoffOffset, "DateAdjustment(s) to apply to selected cutoff date");
    FIELD_MAKE_OPTIONAL(cutoffOffset);
}

CClassConstSP const OnOrBeforeObsSched::TYPE = CClass::registerClassLoadMethod(
    "OnOrBeforeObsSched", typeid(OnOrBeforeObsSched), OnOrBeforeObsSched::load);


///////////////////////////////////////////////////////////////////////////////
// OptionDates - Object to hold an option exercise date and associated
//               notification date.
///////////////////////////////////////////////////////////////////////////////
OptionDates::OptionDates(const DateTime& exerciseDate, const DateTime& notificationDate) :
    CObject(TYPE), exerciseDate(exerciseDate), notificationDate(notificationDate) {}

OptionDates::OptionDates(const CClassConstSP& type) : CObject(type) {}

IObject* OptionDates::defaultConstructor() { return new OptionDates(TYPE); }

void OptionDates::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OptionDates, clazz);
    clazz->setDescription("An option exercise date and associated notification date");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(exerciseDate, "Option exercise date");
    FIELD(notificationDate, "Option notification date");
}

CClassConstSP const OptionDates::TYPE = CClass::registerClassLoadMethod(
    "OptionDates", typeid(OptionDates), OptionDates::load);

DEFINE_TEMPLATE_TYPE(OptionDatesArray);
DEFINE_arrayObjectCast(OptionDates);

///////////////////////////////////////////////////////////////////////////////
// IOptionSched - Interface for classes that generate a schedule of
//                OptionDates.
///////////////////////////////////////////////////////////////////////////////
void IOptionSched::check(bool enforceStrictOrdering) const {
    for (int i = 0; i < size(); ++i) {
        const OptionDates& cur = optionDatesArray[i];
        ensureIncreasing(cur.notificationDate, cur.exerciseDate, "Option notification vs exercise");
        
        if (enforceStrictOrdering && i > 0)
            ensureStrictlyIncreasing(optionDatesArray[i-1].exerciseDate, cur.notificationDate,
                                     "Previous exercise vs current notification");
    }
}

void IOptionSched::fillColumns() const {
    colExerciseDateSP.reset(new DateTimeArray(size()));
    colNotificationDateSP.reset(new DateTimeArray(size()));
    for (int i=0; i<size(); ++i) {
        (*colExerciseDateSP)[i] = optionDatesArray[i].exerciseDate;
        (*colNotificationDateSP)[i] = optionDatesArray[i].notificationDate;
    }
}

IOptionSched::IOptionSched(const CClassConstSP& type) :
    ISched(type) {}

void IOptionSched::load(CClassSP& clazz) {
    clazz->setPublic();     // In case people want to create arrays of these.
    REGISTER(IOptionSched, clazz);
    clazz->setDescription("Interface for a hierarchy of classes that represent a "
                          "schedule of OptionDates");
    SUPERCLASS(ISched);

    FIELD(optionDatesArray, "");
    FIELD_MAKE_TRANSIENT(optionDatesArray);
    FIELD(colExerciseDateSP, "");
    FIELD_MAKE_TRANSIENT(colExerciseDateSP);
    FIELD(colNotificationDateSP, "");
    FIELD_MAKE_TRANSIENT(colNotificationDateSP);
}

CClassConstSP const IOptionSched::TYPE = CClass::registerClassLoadMethod(
    "IOptionSched", typeid(IOptionSched), IOptionSched::load);


///////////////////////////////////////////////////////////////////////////////
// OptionSched - An explicitly specified schedule of OptionDates.
///////////////////////////////////////////////////////////////////////////////
void OptionSched::setup() const { 
    optionDatesArray = *explicitOptionDates; 
    fillColumns();
}

OptionSched::OptionSched(const OptionDatesArrayConstSP& explicitOptionDates) :
    IOptionSched(TYPE), explicitOptionDates(explicitOptionDates) {}

OptionSched::OptionSched(const CClassConstSP& type) :
    IOptionSched(type) {}

IObject* OptionSched::clone() const {
    if (getRefCount() == 0)
        return new OptionSched(*this);

    return const_cast<OptionSched*>(this);
}

IObject* OptionSched::defaultConstructor() { return new OptionSched(TYPE); }

ISchedConstSP OptionSched::schedule() const { return ISchedConstSP(this); }

void OptionSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OptionSched, clazz);
    clazz->setDescription("An explicitly specified array of OptionDates");
    SUPERCLASS(IOptionSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitOptionDates, optionDates, "The array of OptionDates");
}

CClassConstSP const OptionSched::TYPE = CClass::registerClassLoadMethod(
    "OptionSched", typeid(OptionSched), OptionSched::load);


///////////////////////////////////////////////////////////////////////////////
// OptionSchedGen - Generates a schedule of OptionDates from a specified date
//                  schedule.  Separate date adjustments are applied to produce
//                  each exercise and notification date, and the notification
//                  date adjustment may be applied to either the unadjusted
//                  exercise date (i.e. input date schedule) or the adjusted
//                  exercise date.
///////////////////////////////////////////////////////////////////////////////
ISchedConstSP OptionSchedGen::schedule() const { 
    return ISchedConstSP(new OptionSched(
        OptionDatesArrayConstSP(new OptionDatesArray(optionDatesArray)))); 
}

void OptionSchedGen::setup() const {
    const DateTimeArrayConstSP& dateListSP = dateSched->elements();
    const DateTimeArray& dateList = *dateListSP;

    optionDatesArray.reserve(dateList.size());

    DateTime exerciseDate, notificationDate;
    for (int i = 0; i < dateList.size(); ++i)
    {
        exerciseDate = exerciseAdjustment->dateFor(dateList[i]);
        if (notificationRelativeToDate == Exercise)
            notificationDate = notificationAdjustment->dateFor(exerciseDate);
        else
            notificationDate = notificationAdjustment->dateFor(dateList[i]);
        optionDatesArray.push_back(OptionDates(exerciseDate, notificationDate));
    }
    fillColumns();
}

IObject* OptionSchedGen::clone() const {
    if (getRefCount() == 0)
        return new OptionSchedGen(*this);
    else
        return const_cast<OptionSchedGen*>(this);
}

OptionSchedGen::OptionSchedGen(const IDateSchedConstSP& dateSched,
                               const IDateAdjustmentConstSP& exerciseAdjustment,
                               OptionDate notificationRelativeToDate,
                               const IDateAdjustmentConstSP& notificationAdjustment) :
    IOptionSched(TYPE), dateSched(dateSched), exerciseAdjustment(exerciseAdjustment),
    notificationRelativeToDate(notificationRelativeToDate), notificationAdjustment(notificationAdjustment) {
        if (!dateSched)
            throw ModelException(__FUNCTION__, "Date schedule cannot be NULL");
        if (!exerciseAdjustment)
            throw ModelException(__FUNCTION__, "Excercise adjustment cannot be NULL");
        if (!notificationAdjustment)
            throw ModelException(__FUNCTION__, "Notification adjustment cannot be NULL");
}

OptionSchedGen::OptionSchedGen(const CClassConstSP& type) : IOptionSched(type), dateSched(IDateSchedConstSP(0)),
    exerciseAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)), notificationRelativeToDate(Exercise),
    notificationAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)) {}

IObject* OptionSchedGen::defaultConstructor() { return new OptionSchedGen(TYPE); }

void OptionSchedGen::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OptionSchedGen, clazz);
    clazz->setDescription("An option schedule generator");
    SUPERCLASS(IOptionSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(dateSched, "Date schedule that provides unadjusted option exercise dates");
    FIELD(exerciseAdjustment, "Date adjustment(s) to apply to the date schedule to produce option exercise dates");
    FIELD_MAKE_OPTIONAL(exerciseAdjustment);
    FIELD(notificationRelativeToDate, "Indicates whether notification dates are relative to their associated "
                                      "exercise or unadjusted exercise date");
    FIELD_MAKE_OPTIONAL(notificationRelativeToDate);
    FIELD(notificationAdjustment, "Date adjustment(s) to apply to the selected date to produce option notification dates");
    FIELD_MAKE_OPTIONAL(notificationAdjustment);
}

CClassConstSP const OptionSchedGen::TYPE = CClass::registerClassLoadMethod(
    "OptionSchedGen", typeid(OptionSchedGen), OptionSchedGen::load);

START_PUBLIC_ENUM_DEFINITION(OptionSchedGen::OptionDate, "Specifies an option date (e.g. Exercise or UnadjustedExercise)");
ENUM_VALUE_AND_NAME(OptionSchedGen::Exercise, "Exercise", "Indicates the adjusted option exercise date");
ENUM_VALUE_AND_NAME(OptionSchedGen::UnadjustedExercise, "UnadjustedExercise", "Indicates the unadjusted option exercise date");
END_ENUM_DEFINITION(OptionSchedGen::OptionDate);

///////////////////////////////////////////////////////////////////////////////
// FixedCouponDates - Object to hold dates associated with a fixed-rate
//                    coupon payment.
///////////////////////////////////////////////////////////////////////////////
FixedCouponDates::FixedCouponDates(const PeriodConstSP& accrualPeriod, const DateTime& paymentDate) :
    CObject(TYPE), accrualPeriod(accrualPeriod), paymentDate(paymentDate) {
        if (!accrualPeriod)
            throw ModelException(__FUNCTION__, "Accrual period cannot be NULL");
}

FixedCouponDates::FixedCouponDates(const CClassConstSP& type) : CObject(type), accrualPeriod(PeriodConstSP(0)) {}

FixedCouponDates::FixedCouponDates(const CClassConstSP& type, const PeriodConstSP& accrualPeriod, const DateTime& paymentDate) :
    CObject(type), accrualPeriod(accrualPeriod), paymentDate(paymentDate) {
        if (!accrualPeriod)
            throw ModelException(__FUNCTION__, "Accrual period cannot be NULL");
}

void FixedCouponDates::validatePop2Object() {
    validatePeriodDates(accrualPeriod->start, accrualPeriod->end);
}

IObject* FixedCouponDates::defaultConstructor() { return new FixedCouponDates(TYPE); }

void FixedCouponDates::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FixedCouponDates, clazz);
    clazz->setDescription("An object that gives dates associated with a particular fixed coupon");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(accrualPeriod, "Coupon accrual period");
    FIELD(paymentDate, "Coupon payment date");
}

CClassConstSP const FixedCouponDates::TYPE = CClass::registerClassLoadMethod(
    "FixedCouponDates", typeid(FixedCouponDates), FixedCouponDates::load);

DEFINE_TEMPLATE_TYPE(FixedCouponDatesArray);


///////////////////////////////////////////////////////////////////////////////
// FlatFixedCouponDates - A "flattened" version of FixedCouponDates.
///////////////////////////////////////////////////////////////////////////////
FlatFixedCouponDates::FlatFixedCouponDates(const DateTime& accrualStart, const DateTime& accrualEnd,
                                           bool isStub, const DateTime& paymentDate) :
    CObject(TYPE), accrualStart(accrualStart), accrualEnd(accrualEnd), isStub(isStub), paymentDate(paymentDate) {}

FlatFixedCouponDates::FlatFixedCouponDates(const CClassConstSP& type) : CObject(type), isStub(false) {}

FlatFixedCouponDates::FlatFixedCouponDates(const CClassConstSP& type, const DateTime& accrualStart,
                                           const DateTime& accrualEnd, bool isStub, const DateTime& paymentDate) :
    CObject(type), accrualStart(accrualStart), accrualEnd(accrualEnd), isStub(isStub), paymentDate(paymentDate) {}

void FlatFixedCouponDates::validatePop2Object() {
    validatePeriodDates(accrualStart, accrualEnd);
}

IObject* FlatFixedCouponDates::defaultConstructor() { return new FlatFixedCouponDates(TYPE); }

void FlatFixedCouponDates::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FlatFixedCouponDates, clazz);
    clazz->setDescription("A flattened form of FixedCouponDates");
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(accrualStart, "Coupon accrual period start date");
    FIELD(accrualEnd, "Coupon accrual period end date");
    FIELD(isStub, "Whether the accrual period is an irregular length period");
    FIELD(paymentDate, "Coupon payment date");
}

CClassConstSP const FlatFixedCouponDates::TYPE = CClass::registerClassLoadMethod(
    "FlatFixedCouponDates", typeid(FlatFixedCouponDates), FlatFixedCouponDates::load);

DEFINE_TEMPLATE_TYPE(FlatFixedCouponDatesArray);


///////////////////////////////////////////////////////////////////////////////
// IFixedCouponSched - Interface for classes that generate a schedule of
//                     FixedCouponDates.
///////////////////////////////////////////////////////////////////////////////
FixedCouponDatesArrayConstSP IFixedCouponSched::elements() const {
    if (!cachedElemsSP)
        setup();
    return cachedElemsSP;
}

ISchedConstSP IFixedCouponSched::schedule() const { return ISchedConstSP(new FixedCouponSched(elements())); }

ISchedConstSP IFixedCouponSched::flatSchedule() const {
    const FixedCouponDatesArrayConstSP& fixedCouponDatesArraySP = elements();
    const FixedCouponDatesArray& fixedCouponDatesArray = *fixedCouponDatesArraySP;

    FlatFixedCouponDatesArraySP flatFixedCouponDatesSP(new FlatFixedCouponDatesArray);
    FlatFixedCouponDatesArray& flatFixedCouponDates = *flatFixedCouponDatesSP;

    for (int i = 0; i < fixedCouponDatesArray.size(); ++i) {
        const Period& accrualPeriod = *(fixedCouponDatesArray[i]->accrualPeriod);
        const DateTime& paymentDate = fixedCouponDatesArray[i]->paymentDate;

        flatFixedCouponDates.push_back(
                FlatFixedCouponDatesSP(new FlatFixedCouponDates(accrualPeriod.start, accrualPeriod.end,
                                                                accrualPeriod.isStub, paymentDate)));
    }

    return ISchedConstSP(new FlatFixedCouponSched(flatFixedCouponDatesSP));
}

IFixedCouponSched::IFixedCouponSched(const CClassConstSP& type) :
    ISched(type), cachedElemsSP(FixedCouponDatesArrayConstSP(0)) {}

void IFixedCouponSched::load(CClassSP& clazz) {
    clazz->setPublic();     // In case people want to create arrays of these.
    REGISTER(IFixedCouponSched, clazz);
    clazz->setDescription("Interface for a hierarchy of classes that represent a "
                          "schedule of fixed coupons");
    SUPERCLASS(ISched);

    FIELD(cachedElemsSP, "");
    FIELD_MAKE_TRANSIENT(cachedElemsSP);
}

CClassConstSP const IFixedCouponSched::TYPE = CClass::registerClassLoadMethod(
    "IFixedCouponSched", typeid(IFixedCouponSched), IFixedCouponSched::load);


///////////////////////////////////////////////////////////////////////////////
// FixedCouponSched - An explicitly specified schedule of FixedCouponDates.
///////////////////////////////////////////////////////////////////////////////
void FixedCouponSched::setup() const { cachedElemsSP = explicitFixedCouponDates; }

IObject* FixedCouponSched::clone() const {
    if (getRefCount() == 0)
        return new FixedCouponSched(*this);
    else
        return const_cast<FixedCouponSched*>(this);
}

FixedCouponSched::FixedCouponSched(const FixedCouponDatesArrayConstSP& explicitFixedCouponDates) :
    IFixedCouponSched(TYPE), explicitFixedCouponDates(explicitFixedCouponDates) {
        if (!explicitFixedCouponDates)
            throw ModelException(__FUNCTION__, "Explicit FixedCouponDatesArray cannot be NULL");
}

FixedCouponSched::FixedCouponSched(const CClassConstSP& type) :
    IFixedCouponSched(type), explicitFixedCouponDates(FixedCouponDatesArrayConstSP(0)) {}

void FixedCouponSched::validatePop2Object() {
    const FixedCouponDatesArray& fixedCouponDatesArray = *explicitFixedCouponDates;
    for (int i = 1; i < fixedCouponDatesArray.size(); ++i) {
        const FixedCouponDates& cur = *(fixedCouponDatesArray[i - 1]);
        const FixedCouponDates& nxt = *(fixedCouponDatesArray[i]);

        ensureStrictlyIncreasing(cur.accrualPeriod->start, nxt.accrualPeriod->start,
                                 "Accrual period start");
        ensureStrictlyIncreasing(cur.accrualPeriod->end, nxt.accrualPeriod->end,
                                 "Accrual period end");
    }
}

IObject* FixedCouponSched::defaultConstructor() { return new FixedCouponSched(TYPE); }

void FixedCouponSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FixedCouponSched, clazz);
    clazz->setDescription("An explicitly specified array of FixedCouponDates");
    SUPERCLASS(IFixedCouponSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitFixedCouponDates, fixedCouponDates, "The array of FixedCouponDates");
}

CClassConstSP const FixedCouponSched::TYPE = CClass::registerClassLoadMethod(
    "FixedCouponSched", typeid(FixedCouponSched), FixedCouponSched::load);


///////////////////////////////////////////////////////////////////////////////
// FlatFixedCouponSched - An explicitly specified schedule of
//                        FlatFixedCouponDates.
///////////////////////////////////////////////////////////////////////////////
ISchedConstSP FlatFixedCouponSched::flatSchedule() const { return ISchedConstSP(this); }

void FlatFixedCouponSched::setup() const {
    FixedCouponDatesArraySP fixedCouponDatesArraySP(new FixedCouponDatesArray);
    FixedCouponDatesArray& fixedCouponDatesArray = *fixedCouponDatesArraySP;

    const FlatFixedCouponDatesArray& flatFixedCouponDates = *explicitFlatFixedCouponDates;

    PeriodSP accrualPeriodSP(0);
    DateTime accrualStart, accrualEnd, paymentDate;
    bool isStub;

    for (int i = 0; i < flatFixedCouponDates.size(); ++i)
    {
        const FlatFixedCouponDates& aFlatFixedCouponDates = *(flatFixedCouponDates[i]);
        if (i == 0 ||
            aFlatFixedCouponDates.accrualStart != accrualStart ||
            aFlatFixedCouponDates.accrualEnd != accrualEnd)
        {
            if (!!accrualPeriodSP)
                fixedCouponDatesArray.push_back(
                    FixedCouponDatesSP(new FixedCouponDates(accrualPeriodSP, paymentDate)));
            accrualStart = aFlatFixedCouponDates.accrualStart;
            accrualEnd = aFlatFixedCouponDates.accrualEnd;
            paymentDate = aFlatFixedCouponDates.paymentDate;
            isStub = aFlatFixedCouponDates.isStub;
            accrualPeriodSP.reset(new Period(accrualStart, accrualEnd, isStub));
        }
        else if (aFlatFixedCouponDates.isStub != isStub) {
            string msg("Encountered inconsistent stub specifications for accrual period (start=");
            msg += accrualStart.toString() + ", end=" + accrualEnd.toString() + ")";
            throw ModelException(__FUNCTION__, msg);
        }
    }

    if (!!accrualPeriodSP)
        fixedCouponDatesArray.push_back(
            FixedCouponDatesSP(new FixedCouponDates(accrualPeriodSP, paymentDate)));

    cachedElemsSP = fixedCouponDatesArraySP;
}

IObject* FlatFixedCouponSched::clone() const {
    if (getRefCount() == 0)
        return new FlatFixedCouponSched(*this);
    else
        return const_cast<FlatFixedCouponSched*>(this);
}

FlatFixedCouponSched::FlatFixedCouponSched(const FlatFixedCouponDatesArrayConstSP& explicitFlatFixedCouponDates) :
    IFixedCouponSched(TYPE), explicitFlatFixedCouponDates(explicitFlatFixedCouponDates) {
        if (!explicitFlatFixedCouponDates)
            throw ModelException(__FUNCTION__, "Explicit FlatFixedCouponDatesArray cannot be NULL");
}

FlatFixedCouponSched::FlatFixedCouponSched(const CClassConstSP& type) :
    IFixedCouponSched(type), explicitFlatFixedCouponDates(FlatFixedCouponDatesArrayConstSP(0)) {}

void FlatFixedCouponSched::validatePop2Object() {
    const FlatFixedCouponDatesArray& flatFixedCouponDatesArray = *explicitFlatFixedCouponDates;
    for (int i = 1; i < flatFixedCouponDatesArray.size(); ++i) {
        const FlatFixedCouponDates& cur = *(flatFixedCouponDatesArray[i - 1]);
        const FlatFixedCouponDates& nxt = *(flatFixedCouponDatesArray[i]);

        ensureIncreasing(cur.accrualStart, nxt.accrualStart, "Accrual period start");
        ensureIncreasing(cur.accrualEnd, nxt.accrualEnd, "Accrual period end");
    }
}

IObject* FlatFixedCouponSched::defaultConstructor() { return new FlatFixedCouponSched(TYPE); }

void FlatFixedCouponSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FlatFixedCouponSched, clazz);
    clazz->setDescription("An explicitly specified array of FlatFixedCouponDates");
    SUPERCLASS(IFixedCouponSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitFlatFixedCouponDates, flatFixedCouponDates, "The array of FlatFixedCouponDates");
}

CClassConstSP const FlatFixedCouponSched::TYPE = CClass::registerClassLoadMethod(
    "FlatFixedCouponSched", typeid(FlatFixedCouponSched), FlatFixedCouponSched::load);


///////////////////////////////////////////////////////////////////////////////
// FixedCouponSchedGen - Generates a schedule of FixedCouponDates where accrual
//                       periods are given by an IPeriodSched and payment dates
//                       are specified relative to each period start or end
//                       date plus an optional date adjustment.
///////////////////////////////////////////////////////////////////////////////
void FixedCouponSchedGen::setup() const {
    const PeriodArrayConstSP& accrualPeriodsSP = accrualSched->elements();
    const PeriodArray& accrualPeriods = *accrualPeriodsSP;

    DateSchedFromPeriodSched paymentSched(accrualSched, paymentRelativeToDate, paymentDateAdjustment);
    const DateTimeArrayConstSP& paymentDatesSP = paymentSched.elements();
    const DateTimeArray& paymentDates = *paymentDatesSP;

    FixedCouponDatesArraySP fixedCouponDatesArraySP(new FixedCouponDatesArray(accrualPeriods.size()));
    FixedCouponDatesArray& fixedCouponDatesArray = *fixedCouponDatesArraySP;

    for (int i = 0; i < accrualPeriods.size(); ++i)
        fixedCouponDatesArray[i] = FixedCouponDatesSP(new FixedCouponDates(accrualPeriods[i], paymentDates[i]));

    cachedElemsSP = fixedCouponDatesArraySP;
}

IObject* FixedCouponSchedGen::clone() const {
    if (getRefCount() == 0)
        return new FixedCouponSchedGen(*this);
    else
        return const_cast<FixedCouponSchedGen*>(this);
}

void FixedCouponSchedGen::validatePop2Object() {
    if (paymentRelativeToDate == Period::None)
            throw ModelException(__FUNCTION__, "paymentRelativeToDate cannot be None");
}

FixedCouponSchedGen::FixedCouponSchedGen(const IPeriodSchedConstSP& accrualSched,
                                         Period::Boundary paymentRelativeToDate,
                                         const IDateAdjustmentConstSP& paymentDateAdjustment) :
    IFixedCouponSched(TYPE), accrualSched(accrualSched), paymentRelativeToDate(paymentRelativeToDate),
    paymentDateAdjustment(paymentDateAdjustment) {
        if (!accrualSched)
            throw ModelException(__FUNCTION__, "Accrual period schedule cannot be NULL");
        if (!paymentDateAdjustment)
            throw ModelException(__FUNCTION__, "Payment date adjustment cannot be NULL");
        validatePop2Object();
}

FixedCouponSchedGen::FixedCouponSchedGen(const CClassConstSP& type) :
    IFixedCouponSched(type), accrualSched(IPeriodSchedConstSP(0)), paymentRelativeToDate(Period::End),
    paymentDateAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)) {}

IObject* FixedCouponSchedGen::defaultConstructor() { return new FixedCouponSchedGen(TYPE); }

void FixedCouponSchedGen::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FixedCouponSchedGen, clazz);
    clazz->setDescription("Generates an array of FixedCouponDates");
    SUPERCLASS(IFixedCouponSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(accrualSched, "Source of accrual periods");
    FIELD(paymentRelativeToDate, "Selects whether payment date is relative to period start or end date");
    FIELD_MAKE_OPTIONAL(paymentRelativeToDate);
    FIELD(paymentDateAdjustment, "DateAdjustment(s) to apply to source period start or end date to produce a payment date");
    FIELD_MAKE_OPTIONAL(paymentDateAdjustment);
}

CClassConstSP const FixedCouponSchedGen::TYPE = CClass::registerClassLoadMethod(
    "FixedCouponSchedGen", typeid(FixedCouponSchedGen), FixedCouponSchedGen::load);


///////////////////////////////////////////////////////////////////////////////
// FloatCouponDates - Object to hold dates associated with a floating-rate
//                    coupon payment.
///////////////////////////////////////////////////////////////////////////////
FloatCouponDates::FloatCouponDates(const PeriodConstSP& accrualPeriod, const DateTime& paymentDate,
                                   const DateTimeArrayConstSP& resetDates) :
    FixedCouponDates(TYPE, accrualPeriod, paymentDate), resetDates(resetDates) {
        if (!resetDates)
            throw ModelException(__FUNCTION__, "Reset dates array cannot be NULL");
}

FloatCouponDates::FloatCouponDates(const CClassConstSP& type) : FixedCouponDates(type), resetDates(DateTimeArrayConstSP(0)) {}

void FloatCouponDates::validatePop2Object() {
    FixedCouponDates::validatePop2Object();
    DateTime::ensureIncreasing(*resetDates, "reset observation dates", false);
}

IObject* FloatCouponDates::defaultConstructor() { return new FloatCouponDates(TYPE); }

void FloatCouponDates::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FloatCouponDates, clazz);
    clazz->setDescription("An object that gives dates associated with a particular floating coupon");
    SUPERCLASS(FixedCouponDates);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(resetDates, "Floating coupon reset date(s)");
}

CClassConstSP const FloatCouponDates::TYPE = CClass::registerClassLoadMethod(
    "FloatCouponDates", typeid(FloatCouponDates), FloatCouponDates::load);

DEFINE_TEMPLATE_TYPE(FloatCouponDatesArray);


///////////////////////////////////////////////////////////////////////////////
// FlatFloatCouponDates - A "flattened" version of FloatCouponDates.
///////////////////////////////////////////////////////////////////////////////
FlatFloatCouponDates::FlatFloatCouponDates(const DateTime& accrualStart, const DateTime& accrualEnd,
                                           bool isStub, const DateTime& paymentDate, const DateTime& resetDate) :
    FlatFixedCouponDates(TYPE, accrualStart, accrualEnd, isStub, paymentDate), resetDate(resetDate) {}

FlatFloatCouponDates::FlatFloatCouponDates(const CClassConstSP& type) : FlatFixedCouponDates(type) {}

IObject* FlatFloatCouponDates::defaultConstructor() { return new FlatFloatCouponDates(TYPE); }

void FlatFloatCouponDates::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FlatFloatCouponDates, clazz);
    clazz->setDescription("A flattened form of FloatCouponDates");
    SUPERCLASS(FlatFixedCouponDates);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(resetDate, "A floating coupon reset date");
}

CClassConstSP const FlatFloatCouponDates::TYPE = CClass::registerClassLoadMethod(
    "FlatFloatCouponDates", typeid(FlatFloatCouponDates), FlatFloatCouponDates::load);

DEFINE_TEMPLATE_TYPE(FlatFloatCouponDatesArray);


///////////////////////////////////////////////////////////////////////////////
// IFloatCouponSched - Interface for classes that generate a schedule of
//                     FloatCouponDates.
///////////////////////////////////////////////////////////////////////////////
ISchedConstSP IFloatCouponSched::schedule() const { return ISchedConstSP(new FloatCouponSched(elements())); }

FloatCouponDatesArrayConstSP IFloatCouponSched::elements() const {
    if (!cachedElemsSP)
        setup();
    return cachedElemsSP;
}

ISchedConstSP IFloatCouponSched::flatSchedule() const {
    const FloatCouponDatesArrayConstSP& floatCouponDatesArraySP = elements();
    const FloatCouponDatesArray& floatCouponDatesArray = *floatCouponDatesArraySP;

    FlatFloatCouponDatesArraySP flatFloatCouponDatesSP(new FlatFloatCouponDatesArray);
    FlatFloatCouponDatesArray& flatFloatCouponDates = *flatFloatCouponDatesSP;

    for (int i = 0; i < floatCouponDatesArray.size(); ++i) {
        const Period& accrualPeriod = *(floatCouponDatesArray[i]->accrualPeriod);
        const DateTime& paymentDate = floatCouponDatesArray[i]->paymentDate;
        const DateTimeArray& resetDates = *(floatCouponDatesArray[i]->resetDates);

        for (int j = 0; j < resetDates.size(); ++j)
            flatFloatCouponDates.push_back(
                FlatFloatCouponDatesSP(new FlatFloatCouponDates(accrualPeriod.start, accrualPeriod.end,
                                                                accrualPeriod.isStub, paymentDate,
                                                                resetDates[j])));
    }

    return ISchedConstSP(new FlatFloatCouponSched(flatFloatCouponDatesSP));
}

IFloatCouponSched::IFloatCouponSched(const CClassConstSP& type) :
    ISched(type), cachedElemsSP(FloatCouponDatesArrayConstSP(0)) {}

void IFloatCouponSched::load(CClassSP& clazz) {
    clazz->setPublic();     // In case people want to create arrays of these.
    REGISTER(IFloatCouponSched, clazz);
    clazz->setDescription("Interface for a hierarchy of classes that represent a "
                          "schedule of floating coupons");
    SUPERCLASS(ISched);

    FIELD(cachedElemsSP, "");
    FIELD_MAKE_TRANSIENT(cachedElemsSP);
}

CClassConstSP const IFloatCouponSched::TYPE = CClass::registerClassLoadMethod(
    "IFloatCouponSched", typeid(IFloatCouponSched), IFloatCouponSched::load);


///////////////////////////////////////////////////////////////////////////////
// FloatCouponSched - An explicitly specified schedule of FloatCouponDates.
///////////////////////////////////////////////////////////////////////////////
void FloatCouponSched::setup() const { cachedElemsSP = explicitFloatCouponDates; }

IObject* FloatCouponSched::clone() const {
    if (getRefCount() == 0)
        return new FloatCouponSched(*this);
    else
        return const_cast<FloatCouponSched*>(this);
}

FloatCouponSched::FloatCouponSched(const FloatCouponDatesArrayConstSP& explicitFloatCouponDates) :
    IFloatCouponSched(TYPE), explicitFloatCouponDates(explicitFloatCouponDates) {
        if (!explicitFloatCouponDates)
            throw ModelException(__FUNCTION__, "Explicit FloatCouponDatesArray cannot be NULL");
}

FloatCouponSched::FloatCouponSched(const CClassConstSP& type) :
    IFloatCouponSched(type), explicitFloatCouponDates(FloatCouponDatesArrayConstSP(0)) {}

void FloatCouponSched::validatePop2Object() {
    const FloatCouponDatesArray& floatCouponDatesArray = *explicitFloatCouponDates;
    for (int i = 1; i < floatCouponDatesArray.size(); ++i) {
        const FloatCouponDates& cur = *(floatCouponDatesArray[i - 1]);
        const FloatCouponDates& nxt = *(floatCouponDatesArray[i]);
        ensureStrictlyIncreasing(cur.accrualPeriod->start, nxt.accrualPeriod->start,
                                 "Accrual period start");
        ensureStrictlyIncreasing(cur.accrualPeriod->end, nxt.accrualPeriod->end,
                                 "Accrual period end");
    }
}

IObject* FloatCouponSched::defaultConstructor() { return new FloatCouponSched(TYPE); }

void FloatCouponSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FloatCouponSched, clazz);
    clazz->setDescription("An explicitly specified array of FloatCouponDates");
    SUPERCLASS(IFloatCouponSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitFloatCouponDates, floatCouponDates, "The array of FloatCouponDates");
}

CClassConstSP const FloatCouponSched::TYPE = CClass::registerClassLoadMethod(
    "FloatCouponSched", typeid(FloatCouponSched), FloatCouponSched::load);


///////////////////////////////////////////////////////////////////////////////
// FlatFloatCouponSched - An explicitly specified schedule of
//                        FlatFloatCouponDates.
///////////////////////////////////////////////////////////////////////////////
ISchedConstSP FlatFloatCouponSched::flatSchedule() const { return ISchedConstSP(this); }

void FlatFloatCouponSched::setup() const {
    FloatCouponDatesArraySP floatCouponDatesArraySP(new FloatCouponDatesArray);
    FloatCouponDatesArray& floatCouponDatesArray = *floatCouponDatesArraySP;

    const FlatFloatCouponDatesArray& flatFloatCouponDates = *explicitFlatFloatCouponDates;

    PeriodSP accrualPeriodSP(0);
    DateTimeArraySP resetDatesSP(0);
    DateTime accrualStart, accrualEnd, paymentDate;
    bool isStub;

    for (int i = 0; i < flatFloatCouponDates.size(); ++i)
    {
        const FlatFloatCouponDates& aFlatFloatCouponDates = *(flatFloatCouponDates[i]);
        if (i == 0 ||
            aFlatFloatCouponDates.accrualStart != accrualStart ||
            aFlatFloatCouponDates.accrualEnd != accrualEnd)
        {
            if (!!accrualPeriodSP)
                floatCouponDatesArray.push_back(
                    FloatCouponDatesSP(new FloatCouponDates(accrualPeriodSP, paymentDate, resetDatesSP)));
            accrualStart = aFlatFloatCouponDates.accrualStart;
            accrualEnd = aFlatFloatCouponDates.accrualEnd;
            paymentDate = aFlatFloatCouponDates.paymentDate;
            isStub = aFlatFloatCouponDates.isStub;
            accrualPeriodSP.reset(new Period(accrualStart, accrualEnd, isStub));
            resetDatesSP.reset(new DateTimeArray);
        }
        else if (aFlatFloatCouponDates.isStub != isStub) {
            string msg("Encountered inconsistent stub specifications for accrual period (start=");
            msg += accrualStart.toString() + ", end=" + accrualEnd.toString() + ")";
            throw ModelException(__FUNCTION__, msg);
        }

        resetDatesSP->push_back(aFlatFloatCouponDates.resetDate);
    }

    if (!!accrualPeriodSP)
        floatCouponDatesArray.push_back(
            FloatCouponDatesSP(new FloatCouponDates(accrualPeriodSP, paymentDate, resetDatesSP)));

    cachedElemsSP = floatCouponDatesArraySP;
}

IObject* FlatFloatCouponSched::clone() const {
    if (getRefCount() == 0)
        return new FlatFloatCouponSched(*this);
    else
        return const_cast<FlatFloatCouponSched*>(this);
}

FlatFloatCouponSched::FlatFloatCouponSched(const FlatFloatCouponDatesArrayConstSP& explicitFlatFloatCouponDates) :
    IFloatCouponSched(TYPE), explicitFlatFloatCouponDates(explicitFlatFloatCouponDates) {
        if (!explicitFlatFloatCouponDates)
            throw ModelException(__FUNCTION__, "Explicit FlatFloatCouponDatesArray cannot be NULL");
}

FlatFloatCouponSched::FlatFloatCouponSched(const CClassConstSP& type) :
    IFloatCouponSched(type), explicitFlatFloatCouponDates(FlatFloatCouponDatesArrayConstSP(0)) {}

void FlatFloatCouponSched::validatePop2Object() {
    const FlatFloatCouponDatesArray& flatFloatCouponDatesArray = *explicitFlatFloatCouponDates;
    for (int i = 1; i < flatFloatCouponDatesArray.size(); ++i) {
        const FlatFloatCouponDates& cur = *(flatFloatCouponDatesArray[i - 1]);
        const FlatFloatCouponDates& nxt = *(flatFloatCouponDatesArray[i]);

        ensureIncreasing(cur.accrualStart, nxt.accrualStart, "Accrual period start");
        ensureIncreasing(cur.accrualEnd, nxt.accrualEnd, "Accural period end");

        if (cur.accrualStart == nxt.accrualStart && cur.accrualEnd == nxt.accrualEnd)
            ensureIncreasing(cur.resetDate, nxt.resetDate, "Reset observation");
    }
}

IObject* FlatFloatCouponSched::defaultConstructor() { return new FlatFloatCouponSched(TYPE); }

void FlatFloatCouponSched::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FlatFloatCouponSched, clazz);
    clazz->setDescription("An explicitly specified array of FlatFloatCouponDates");
    SUPERCLASS(IFloatCouponSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD_USING_ALIAS(explicitFlatFloatCouponDates, flatFloatCouponDates, "The array of FlatFloatCouponDates");
}

CClassConstSP const FlatFloatCouponSched::TYPE = CClass::registerClassLoadMethod(
    "FlatFloatCouponSched", typeid(FlatFloatCouponSched), FlatFloatCouponSched::load);


///////////////////////////////////////////////////////////////////////////////
// FloatCouponSchedGen - Generates a schedule of FloatCouponDates where accrual
//                       periods are given by an IPeriodSched and payment dates
//                       are specified relative to each period start or end
//                       date plus an optional date adjustment.  Reset
//                       observation dates are specified by an IObsSched that
//                       must produce PeriodObsDates that correspond to the
//                       accrual periods (both in number and order).
///////////////////////////////////////////////////////////////////////////////
void FloatCouponSchedGen::setup() const {
    const PeriodArrayConstSP& accrualPeriodsSP = accrualSched->elements();
    const PeriodArray& accrualPeriods = *accrualPeriodsSP;
    const PeriodObsDatesArrayConstSP& resetPeriodObsSP = resetSched->elements();
    const PeriodObsDatesArray& resetPeriodObs = *resetPeriodObsSP;

    if (resetPeriodObs.size() != accrualPeriods.size())
    {
        throw ModelException(__FUNCTION__, "Mismatched number of reset observation "
                                           "periods with respect to accrual periods");
    }

    DateSchedFromPeriodSched paymentSched(accrualSched, paymentRelativeToDate, paymentDateAdjustment);
    const DateTimeArrayConstSP& paymentDatesSP = paymentSched.elements();
    const DateTimeArray& paymentDates = *paymentDatesSP;

    FloatCouponDatesArraySP floatCouponDatesArraySP(new FloatCouponDatesArray(accrualPeriods.size()));
    FloatCouponDatesArray& floatCouponDatesArray = *floatCouponDatesArraySP;

    for (int i = 0; i < accrualPeriods.size(); ++i) {
        const Period& accrualPeriod = *(accrualPeriods[i]);

        if (*(resetPeriodObs[i]->period) != accrualPeriod)
        {
            string msg("Mismatched reset observation period ");
            msg += resetPeriodObs[i]->period->toString() + " for accrual period "
                + accrualPeriod.toString() + ".  i = " + Format::toString(i);
            throw ModelException(__FUNCTION__, msg);
        }

        floatCouponDatesArray[i] = FloatCouponDatesSP(
                                        new FloatCouponDates(
                                                accrualPeriods[i],
                                                paymentDates[i],
                                                resetPeriodObs[i]->obsDates));
    }

    cachedElemsSP = floatCouponDatesArraySP;
}

IObject* FloatCouponSchedGen::clone() const {
    if (getRefCount() == 0)
        return new FloatCouponSchedGen(*this);
    else
        return const_cast<FloatCouponSchedGen*>(this);
}

void FloatCouponSchedGen::validatePop2Object() {
    if (paymentRelativeToDate == Period::None)
            throw ModelException(__FUNCTION__, "paymentRelativeToDate cannot be None");
}

FloatCouponSchedGen::FloatCouponSchedGen(const IPeriodSchedConstSP& accrualSched,
                                         Period::Boundary paymentRelativeToDate,
                                         const IDateAdjustmentConstSP& paymentDateAdjustment,
                                         const IObsSchedConstSP& resetSched) :
    IFloatCouponSched(TYPE), accrualSched(accrualSched), paymentRelativeToDate(paymentRelativeToDate),
    paymentDateAdjustment(paymentDateAdjustment), resetSched(resetSched) {
        if (!accrualSched)
            throw ModelException(__FUNCTION__, "Accrual period schedule cannot be NULL");
        if (!paymentDateAdjustment)
            throw ModelException(__FUNCTION__, "Payment date adjustment cannot be NULL");
        if (!resetSched)
            throw ModelException(__FUNCTION__, "Reset observation schedule cannot be NULL");
        validatePop2Object();
}

FloatCouponSchedGen::FloatCouponSchedGen(const CClassConstSP& type) :
    IFloatCouponSched(type), accrualSched(IPeriodSchedConstSP(0)), paymentRelativeToDate(Period::End),
    paymentDateAdjustment(IDateAdjustmentConstSP(new NoDateAdjustment)), resetSched(IObsSchedConstSP(0)) {}

IObject* FloatCouponSchedGen::defaultConstructor() { return new FloatCouponSchedGen(TYPE); }

void FloatCouponSchedGen::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FloatCouponSchedGen, clazz);
    clazz->setDescription("Generates an array of FloatCouponDates");
    SUPERCLASS(IFloatCouponSched);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(accrualSched, "Source of accrual periods");
    FIELD(paymentRelativeToDate, "Selects whether payment date is relative to period start or end date");
    FIELD_MAKE_OPTIONAL(paymentRelativeToDate);
    FIELD(paymentDateAdjustment, "DateAdjustment(s) to apply to source period start or end date to produce a payment date");
    FIELD_MAKE_OPTIONAL(paymentDateAdjustment);
    FIELD(resetSched, "Schedule of reset observation dates per accrual period");
}

CClassConstSP const FloatCouponSchedGen::TYPE = CClass::registerClassLoadMethod(
    "FloatCouponSchedGen", typeid(FloatCouponSchedGen), FloatCouponSchedGen::load);


///////////////////////////////////////////////////////////////////////////////
// Function to force linker to include this object file
///////////////////////////////////////////////////////////////////////////////
bool ScheduleGeneratorLoad() {
    return ISched::TYPE != 0;
}

DRLIB_END_NAMESPACE
