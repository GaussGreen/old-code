//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : DateAdjustment.cpp
//
//   Description : Various transformations that produce a date from a given date.
//
//   Author      : Jerry Cohen
//
//   Date        : November 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_DATEADJUSTMENT_CPP
#include "edginc/DateAdjustment.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/BoxedEnum.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// IDateAdjustment - Top-level interface for date adjustments (i.e.
//                   transformations on a date).
///////////////////////////////////////////////////////////////////////////////
void IDateAdjustment::load(CClassSP& clazz) {
    clazz->setPublic(); // Needed in order to be able to create IDateAdjustmentArray.
    REGISTER_INTERFACE(IDateAdjustment, clazz);
    clazz->setDescription("Interface for various date adjustments that produce a (new) date for a given date");
    EXTENDS(IObject);
}

CClassConstSP const IDateAdjustment::TYPE = CClass::registerInterfaceLoadMethod(
    "IDateAdjustment", typeid(IDateAdjustment), IDateAdjustment::load);

DEFINE_TEMPLATE_TYPE(IDateAdjustmentArray);


///////////////////////////////////////////////////////////////////////////////
// DateAdjustmentSequence - A collection of zero or more DateAdjustments that
//                          are applied in sequence to the base date.
///////////////////////////////////////////////////////////////////////////////
DateTime DateAdjustmentSequence::dateFor(const DateTime& baseDate) const {

    DateTime date(baseDate);
    for (int i = 0; i < adjustments->size(); ++i)
        date = (*adjustments)[i]->dateFor(date);

    return date;
}

DateAdjustmentSequence::DateAdjustmentSequence(const IDateAdjustmentArrayConstSP& adjustments) :
    CObject(TYPE), adjustments(adjustments) {
        if (!adjustments)
            throw ModelException(__FUNCTION__, "DateAdjustmentArray cannot be NULL");
}

DateAdjustmentSequence::DateAdjustmentSequence(const CClassConstSP& type) :
    CObject(type), adjustments(IDateAdjustmentArrayConstSP(0)) {}

IObject* DateAdjustmentSequence::defaultConstructor() { return new DateAdjustmentSequence(TYPE); }

void DateAdjustmentSequence::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(DateAdjustmentSequence, clazz);
    clazz->setDescription("A collection of zero or more date adjustments that are to be applied in sequence");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(adjustments, "A collection of date adjustments to apply to a base date produce a new date");
}

CClassConstSP const DateAdjustmentSequence::TYPE = CClass::registerClassLoadMethod(
    "DateAdjustmentSequence", typeid(DateAdjustmentSequence), DateAdjustmentSequence::load);


///////////////////////////////////////////////////////////////////////////////
// NoDateAdjustment - Answers the base date unchanged.  This can be a
//                    convenient default adjustment.
///////////////////////////////////////////////////////////////////////////////
IObject* NoDateAdjustment::defaultConstructor() { return new NoDateAdjustment; }

void NoDateAdjustment::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(NoDateAdjustment, clazz);
    clazz->setDescription("Does not perform any adjustment and simply answers the base date");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const NoDateAdjustment::TYPE = CClass::registerClassLoadMethod(
    "NoDateAdjustment", typeid(NoDateAdjustment), NoDateAdjustment::load);


///////////////////////////////////////////////////////////////////////////////
// FixedDate - Always answers a specified fixed date regardless of the base
//             date.
///////////////////////////////////////////////////////////////////////////////
DateTime FixedDate::dateFor(const DateTime& baseDate) const { return date; }

FixedDate::FixedDate(const DateTime& date) : CObject(TYPE), date(date) {}

FixedDate::FixedDate(const CClassConstSP& type) : CObject(type) {}

IObject* FixedDate::defaultConstructor() { return new FixedDate(TYPE); }

void FixedDate::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FixedDate, clazz);
    clazz->setDescription("Simply answers a fixed date regardless of the base date");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(date, "The fixed date to answer for any base date");
}

CClassConstSP const FixedDate::TYPE = CClass::registerClassLoadMethod(
    "FixedDate", typeid(FixedDate), FixedDate::load);


///////////////////////////////////////////////////////////////////////////////
// IntervalOffset - Answers a date that is a specified interval (e.g. "3M")
//                  from the base date.
///////////////////////////////////////////////////////////////////////////////
DateTime IntervalOffset::dateFor(const DateTime& baseDate) const {
    return maturityPeriod->toDate(baseDate);
}

IntervalOffset::IntervalOffset(const MaturityPeriodConstSP& maturityPeriod) :
    CObject(TYPE), maturityPeriod(maturityPeriod) {
        if (!maturityPeriod)
            throw ModelException(__FUNCTION__, "maturityPeriod cannot be NULL");
}

IntervalOffset::IntervalOffset(const CClassConstSP& type) :
    CObject(type), maturityPeriod(MaturityPeriodConstSP(0)) {}

IObject* IntervalOffset::defaultConstructor() { return new IntervalOffset(TYPE); }

void IntervalOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(IntervalOffset, clazz);
    clazz->setDescription("Produces a date that is a specified interval from a given base date");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(maturityPeriod, "The interval (e.g. 3 months) to offset a base date to produce a new date");
}

CClassConstSP const IntervalOffset::TYPE = CClass::registerClassLoadMethod(
    "IntervalOffset", typeid(IntervalOffset), IntervalOffset::load);


///////////////////////////////////////////////////////////////////////////////
// NthDayOffset - Answers a date that is the nth day of the month of the base
//                date (or throws an exception if n is not valid for that
//                month).
///////////////////////////////////////////////////////////////////////////////
DateTime NthDayOffset::dateFor(const DateTime& baseDate) const {
    DateTime::MonthDayYear mdy(baseDate.toMDY());
    mdy.day = n;
    return mdy.toDateTime();
}

NthDayOffset::NthDayOffset(int n) : CObject(TYPE), n(n) {}

void NthDayOffset::validatePop2Object() {
    if (n < 1 || n > 31) {
        throw ModelException(__FUNCTION__, "The day of the month must be between 1 and 31");
    }
}

IObject* NthDayOffset::defaultConstructor() { return new NthDayOffset(0); }

void NthDayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(NthDayOffset, clazz);
    clazz->setDescription("Answers the nth day of the month of the base date");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(n, "The desired day of the month of the base date")
}

CClassConstSP const NthDayOffset::TYPE = CClass::registerClassLoadMethod(
    "NthDayOffset", typeid(NthDayOffset), NthDayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// LastDayOffset - Answers the date of the last day of the month of the base
//                 date.
///////////////////////////////////////////////////////////////////////////////
DateTime LastDayOffset::dateFor(const DateTime& baseDate) const {
    return baseDate.returnEndOfMonth(ignoreLeapYears);
}

LastDayOffset::LastDayOffset(bool ignoreLeapYears) :
    CObject(TYPE), ignoreLeapYears(ignoreLeapYears) {}

IObject* LastDayOffset::defaultConstructor() { return new LastDayOffset; }

void LastDayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LastDayOffset, clazz);
    clazz->setDescription("Answers the last day of the month of the base date");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ignoreLeapYears, "Specifies whether February should always be considered to have just 28 days");
    FIELD_MAKE_OPTIONAL(ignoreLeapYears);
}

CClassConstSP const LastDayOffset::TYPE = CClass::registerClassLoadMethod(
    "LastDayOffset", typeid(LastDayOffset), LastDayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// WeekdayOffset - Abstract base class for DateAdjustments that attempt to
//                 answer the date of a specified day of the week (e.g. Monday)
//                 relative to the base date.
///////////////////////////////////////////////////////////////////////////////
WeekdayOffset::WeekdayOffset(const CClassConstSP& type, DateTime::DayOfWeek weekday) : CObject(type), weekday(weekday) {}

void WeekdayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(WeekdayOffset, clazz);
    clazz->setDescription("Common base class for adjustments that answer a particular weekday relative to the base date");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);

    FIELD(weekday, "A day of the week");
}

CClassConstSP const WeekdayOffset::TYPE = CClass::registerClassLoadMethod(
    "WeekdayOffset", typeid(WeekdayOffset), WeekdayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// NthWeekdayOffset - Attempts to answer the date of the nth specified day of
//                    the week (e.g. 3rd Wednesday) of the month of the base
//                    date.
///////////////////////////////////////////////////////////////////////////////
DateTime NthWeekdayOffset::dateFor(const DateTime& baseDate) const {
    return baseDate.nthWeekdayOfMonth(n, weekday);
}

NthWeekdayOffset::NthWeekdayOffset(int n, DateTime::DayOfWeek weekday) : WeekdayOffset(TYPE, weekday), n(n) {}

void NthWeekdayOffset::validatePop2Object() {
    if (n < 1 || n > 5) {
        throw ModelException(__FUNCTION__, "The ordinal number for the nth weekday must be between 1 and 5");
    }
}

IObject* NthWeekdayOffset::defaultConstructor() { return new NthWeekdayOffset(0, DateTime::Sunday); }

void NthWeekdayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(NthWeekdayOffset, clazz);
    clazz->setDescription("Answers the nth specified weekday of the month of the base date");
    SUPERCLASS(WeekdayOffset);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(n, "The ordinal number of the specified weekday to answer (e.g. 1 = first)")
}

CClassConstSP const NthWeekdayOffset::TYPE = CClass::registerClassLoadMethod(
    "NthWeekdayOffset", typeid(NthWeekdayOffset), NthWeekdayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// LastWeekdayOffset - Answers the date of the last specified day of the week
//                     (e.g. last Wednesday) of the month of the base date.
///////////////////////////////////////////////////////////////////////////////
DateTime LastWeekdayOffset::dateFor(const DateTime& baseDate) const {
    return baseDate.lastWeekdayOfMonth(weekday);
}

LastWeekdayOffset::LastWeekdayOffset(DateTime::DayOfWeek weekday) : WeekdayOffset(TYPE, weekday) {}

IObject* LastWeekdayOffset::defaultConstructor() { return new LastWeekdayOffset(DateTime::Sunday); }

void LastWeekdayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LastWeekdayOffset, clazz);
    clazz->setDescription("Answers the last specified weekday of the month of the base date");
    SUPERCLASS(WeekdayOffset);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const LastWeekdayOffset::TYPE = CClass::registerClassLoadMethod(
    "LastWeekdayOffset", typeid(LastWeekdayOffset), LastWeekdayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// OnOrBeforeWeekdayOffset - Answers the date of the specified day of the week
//                           that is on or before the base date.
///////////////////////////////////////////////////////////////////////////////
DateTime OnOrBeforeWeekdayOffset::dateFor(const DateTime& baseDate) const {
    int daysPastTarget = baseDate.getWeekday() - weekday;

    if (daysPastTarget < 0)
        daysPastTarget += 7;

    return baseDate.rollDate(-daysPastTarget);
}

OnOrBeforeWeekdayOffset::OnOrBeforeWeekdayOffset(DateTime::DayOfWeek weekday) : WeekdayOffset(TYPE, weekday) {}

IObject* OnOrBeforeWeekdayOffset::defaultConstructor() { return new OnOrBeforeWeekdayOffset(DateTime::Sunday); }

void OnOrBeforeWeekdayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OnOrBeforeWeekdayOffset, clazz);
    clazz->setDescription("Answers the nearest specified weekday that is on or before the base date");
    SUPERCLASS(WeekdayOffset);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const OnOrBeforeWeekdayOffset::TYPE = CClass::registerClassLoadMethod(
    "OnOrBeforeWeekdayOffset", typeid(OnOrBeforeWeekdayOffset), OnOrBeforeWeekdayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// OnOrAfterWeekdayOffset - Answers the date of the specified day of the week
//                          that is on or after the base date.
///////////////////////////////////////////////////////////////////////////////
DateTime OnOrAfterWeekdayOffset::dateFor(const DateTime& baseDate) const {
    int daysBeforeTarget = weekday - baseDate.getWeekday();

    if (daysBeforeTarget < 0)
        daysBeforeTarget += 7;

    return baseDate.rollDate(daysBeforeTarget);
}

OnOrAfterWeekdayOffset::OnOrAfterWeekdayOffset(DateTime::DayOfWeek weekday) : WeekdayOffset(TYPE, weekday) {}

IObject* OnOrAfterWeekdayOffset::defaultConstructor() { return new OnOrAfterWeekdayOffset(DateTime::Sunday); }

void OnOrAfterWeekdayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OnOrAfterWeekdayOffset, clazz);
    clazz->setDescription("Answers the nearest specified weekday that is on or after the base date");
    SUPERCLASS(WeekdayOffset);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const OnOrAfterWeekdayOffset::TYPE = CClass::registerClassLoadMethod(
    "OnOrAfterWeekdayOffset", typeid(OnOrAfterWeekdayOffset), OnOrAfterWeekdayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// SameWeekWeekdayOffset - Answers the date of the specified day of the week
//                         that is in the same week as the base date.
//                         NOTE:  Weeks are considered to begin on Monday and
//                                end the following Sunday.
///////////////////////////////////////////////////////////////////////////////
DateTime SameWeekWeekdayOffset::dateFor(const DateTime& baseDate) const {
    // TO REVIEW:  QLib starts weeks on Saturday.  To be consistent with Kapital (at least),
    // we need to adjust the weekday indices to make weeks start on Monday.
    const int weekdayIndexAdjustment = DateTime::Monday - DateTime::Saturday;
    int targetWeekdayIndex = weekday - weekdayIndexAdjustment;
    if (targetWeekdayIndex < 0)
        targetWeekdayIndex += 7;
    int weekdayIndex = baseDate.getWeekday() - weekdayIndexAdjustment;
    if (weekdayIndex < 0)
        weekdayIndex += 7;

    return baseDate.rollDate(targetWeekdayIndex - weekdayIndex);
}

SameWeekWeekdayOffset::SameWeekWeekdayOffset(DateTime::DayOfWeek weekday) : WeekdayOffset(TYPE, weekday) {}

IObject* SameWeekWeekdayOffset::defaultConstructor() { return new SameWeekWeekdayOffset(DateTime::Sunday); }

void SameWeekWeekdayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(SameWeekWeekdayOffset, clazz);
    clazz->setDescription("Answers the specified weekday that is in the same week as the base date");
    SUPERCLASS(WeekdayOffset);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const SameWeekWeekdayOffset::TYPE = CClass::registerClassLoadMethod(
    "SameWeekWeekdayOffset", typeid(SameWeekWeekdayOffset), SameWeekWeekdayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// BusDateAdjustment - Abstract base class for DateAdjustments that operate
//                     in the context of a business calendar (i.e. Holiday).
///////////////////////////////////////////////////////////////////////////////
BusDateAdjustment::BusDateAdjustment(const CClassConstSP& type,
                                     const HolidayWrapper& holidays) :
    CObject(type), holidays(holidays) {}

BusDateAdjustment::BusDateAdjustment(const CClassConstSP& type) : CObject(type) {}

void BusDateAdjustment::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(BusDateAdjustment, clazz);
    clazz->setDescription("Abstract base class for a hierarchy of date adjustments that operate "
                          "in the context of zero or more business calendars");
    SUPERCLASS(CObject);
    IMPLEMENTS(IDateAdjustment);

    FIELD(holidays, "The holidays that determine whether a date is a good business day");
}

CClassConstSP const BusDateAdjustment::TYPE = CClass::registerClassLoadMethod(
    "BusDateAdjustment", typeid(BusDateAdjustment), BusDateAdjustment::load);


///////////////////////////////////////////////////////////////////////////////
// BusDaysOffset - Answers a date that is a specified number of business days
//                 from the base date.
///////////////////////////////////////////////////////////////////////////////
DateTime BusDaysOffset::dateFor(const DateTime& baseDate) const
{
    return holidays->addBusinessDays(baseDate, n);
}

BusDaysOffset::BusDaysOffset(int n, const HolidayWrapper& holidays) :
    BusDateAdjustment(TYPE, holidays), n(n) {}

BusDaysOffset::BusDaysOffset(const CClassConstSP& type) :
    BusDateAdjustment(type), n(0) {}

IObject* BusDaysOffset::defaultConstructor() { return new BusDaysOffset(TYPE); }

void BusDaysOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(BusDaysOffset, clazz);
    clazz->setDescription("Produces a date that is a specified number of business days from a given base date");
    SUPERCLASS(BusDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(n, "The number of business days to add to the base date to produce a new date");
}

CClassConstSP const BusDaysOffset::TYPE = CClass::registerClassLoadMethod(
    "BusDaysOffset", typeid(BusDaysOffset), BusDaysOffset::load);


///////////////////////////////////////////////////////////////////////////////
// NthBusDayOffset - Answers the nth business day of the month of the base date
//                   (or throws an exception if no such day exists).
///////////////////////////////////////////////////////////////////////////////
DateTime NthBusDayOffset::dateFor(const DateTime& baseDate) const
{
    BadDayAdjustment badDayFollowing(BadDayAdjustment::Following, holidays);
    DateTime firstBusDay(badDayFollowing.dateFor(baseDate.firstOfMonth()));
    DateTime nthBusDay(holidays->addBusinessDays(firstBusDay, n - 1));

    // Ensure that we are still in the same month as baseDate.
    DateTime::MonthDayYear baseMDY(baseDate.toMDY());
    DateTime::MonthDayYear newMDY(nthBusDay.toMDY());
    if (newMDY.month != baseMDY.month)
    {
        string msg("Invalid business day number (");
        msg += Format::toString(n) + ") for base date: " + baseDate.toString();
        throw ModelException(__FUNCTION__, msg);
    }

    return nthBusDay;
}

NthBusDayOffset::NthBusDayOffset(int n, const HolidayWrapper& holidays) :
    BusDateAdjustment(TYPE, holidays), n(n) {}

NthBusDayOffset::NthBusDayOffset(const CClassConstSP& type) :
    BusDateAdjustment(type), n(0) {}

void NthBusDayOffset::validatePop2Object() {
    if (n < 1 || n > 31) {
        throw ModelException(__FUNCTION__, "The specified business day must be between 1 and 31");
    }
}

IObject* NthBusDayOffset::defaultConstructor() { return new NthBusDayOffset(TYPE); }

void NthBusDayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(NthBusDayOffset, clazz);
    clazz->setDescription("Produces a date that is the nth business day of the month of a given base date");
    SUPERCLASS(BusDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(n, "The desired business day of the month of the base date");
}

CClassConstSP const NthBusDayOffset::TYPE = CClass::registerClassLoadMethod(
    "NthBusDayOffset", typeid(NthBusDayOffset), NthBusDayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// LastBusDayOffset - Answers the last business day of the month of the base
//                    date (or throws an exception if no such day exists).
///////////////////////////////////////////////////////////////////////////////
DateTime LastBusDayOffset::dateFor(const DateTime& baseDate) const
{
    DateTime endOfMonth(baseDate.returnEndOfMonth(ignoreLeapYears));

    BadDayAdjustment badDayPrevious(BadDayAdjustment::Previous, holidays);
    DateTime lastBusDay(badDayPrevious.dateFor(endOfMonth));

    // Ensure that we are still in the same month as baseDate (e.g. if all days in the
    // month are marked as holidays).
    DateTime::MonthDayYear baseMDY(baseDate.toMDY());
    DateTime::MonthDayYear newMDY(lastBusDay.toMDY());
    if (newMDY.month != baseMDY.month)
    {
        string msg("There are no valid business days in the month of the base date: ");
        msg += baseDate.toString();
        throw ModelException(__FUNCTION__, msg);
    }

    return lastBusDay;
}

LastBusDayOffset::LastBusDayOffset(const HolidayWrapper& holidays, bool ignoreLeapYears) :
    BusDateAdjustment(TYPE, holidays), ignoreLeapYears(ignoreLeapYears) {}

LastBusDayOffset::LastBusDayOffset(const CClassConstSP& type) :
    BusDateAdjustment(type), ignoreLeapYears(false) {}

IObject* LastBusDayOffset::defaultConstructor() { return new LastBusDayOffset(TYPE); }

void LastBusDayOffset::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LastBusDayOffset, clazz);
    clazz->setDescription("Produces a date that is the last business day of the month of a given base date");
    SUPERCLASS(BusDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ignoreLeapYears, "Specifies whether February should always be considered to have just 28 days");
    FIELD_MAKE_OPTIONAL(ignoreLeapYears);
}

CClassConstSP const LastBusDayOffset::TYPE = CClass::registerClassLoadMethod(
    "LastBusDayOffset", typeid(LastBusDayOffset), LastBusDayOffset::load);


///////////////////////////////////////////////////////////////////////////////
// BadDayAdjustment - Answers the base date if it is a good business day.
//                    Otherwise, the specified bad day convention is applied
//                    to produce an appropriate good business day.
///////////////////////////////////////////////////////////////////////////////
DateTime BadDayAdjustment::dateFor(const DateTime& baseDate) const {
    static const char* conventionStrings[] = {"N", "P", "F", "M" };

    BadDayConventionSP _convention(BadDayConventionFactory::make(conventionStrings[convention]));
    return _convention->adjust(baseDate, holidays.get());
}

BadDayAdjustment::BadDayAdjustment(ConventionName convention,
                                   const HolidayWrapper& holidays) :
    BusDateAdjustment(TYPE, holidays), convention(convention) {}

BadDayAdjustment::BadDayAdjustment(const CClassConstSP& type) :
    BusDateAdjustment(TYPE), convention(None) {}

IObject* BadDayAdjustment::defaultConstructor() { return new BadDayAdjustment(TYPE); }

void BadDayAdjustment::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(BadDayAdjustment, clazz);
    clazz->setDescription("Adjusts a date that is not a good business day using a specified bad day convention");
    SUPERCLASS(BusDateAdjustment);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(convention, "The name of the convention to apply");
}

CClassConstSP const BadDayAdjustment::TYPE = CClass::registerClassLoadMethod(
    "BadDayAdjustment", typeid(BadDayAdjustment), BadDayAdjustment::load);

START_PUBLIC_ENUM_DEFINITION(BadDayAdjustment::ConventionName, "Specifies the name of a bad day convention");
ENUM_VALUE_AND_NAME(BadDayAdjustment::None, "None",
                    "Answers the base date regardless of whether it is a bad business day");
ENUM_VALUE_AND_NAME(BadDayAdjustment::Previous, "Previous",
                    "Answers the first good business day that is on or before the base date");
ENUM_VALUE_AND_NAME(BadDayAdjustment::Following, "Following",
                    "Answers the first good business day that is on or after the base date");
ENUM_VALUE_AND_NAME(BadDayAdjustment::ModifiedFollowing, "ModifiedFollowing",
                    "Answers the first good business day that is on or after the base date but within the same month.  "
                    "Otherwise, answers the previous good business day");
END_ENUM_DEFINITION(BadDayAdjustment::ConventionName);


///////////////////////////////////////////////////////////////////////////////
// Force linker to include this object file.
///////////////////////////////////////////////////////////////////////////////
bool DateAdjustmentLoad() {
    return IDateAdjustment::TYPE != 0;
}

DRLIB_END_NAMESPACE
