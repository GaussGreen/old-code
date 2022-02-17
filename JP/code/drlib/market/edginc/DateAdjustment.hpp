//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : DateAdjustment.hpp
//
//   Description : Various transformations that produce a date from a given date.
//
//   Author      : Jerry Cohen
//
//   Date        : November 2006
//
//----------------------------------------------------------------------------

#ifndef DateAdjustment_HPP
#define DateAdjustment_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/MarketData_forward.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// IDateAdjustment - Top-level interface for date adjustments (i.e.
//                   transformations on a date).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL IDateAdjustment : public virtual IObject {
public:
    static CClassConstSP const TYPE;

    // The main DateAdjustment interface.
    virtual DateTime dateFor(const DateTime& baseDate) const = 0;

protected:
    IDateAdjustment() {}

private:
    static void load(CClassSP& clazz);
};

DECLARE(IDateAdjustment);

///////////////////////////////////////////////////////////////////////////////
// DateAdjustmentSequence - A collection of zero or more DateAdjustments that
//                          are applied in sequence to the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL DateAdjustmentSequence: public CObject,
                                         public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

    DateAdjustmentSequence(const IDateAdjustmentArrayConstSP& adjustments);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    DateAdjustmentSequence(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IDateAdjustmentArrayConstSP adjustments;
};

///////////////////////////////////////////////////////////////////////////////
// NoDateAdjustment - Answers the base date unchanged.  This can be a
//                    convenient default adjustment.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL NoDateAdjustment: public CObject,
                                   public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

    NoDateAdjustment() : CObject(TYPE) {}

    virtual DateTime dateFor(const DateTime& baseDate) const { return baseDate; }

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// FixedDate - Always answers a specified fixed date regardless of the base
//             date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FixedDate : public CObject,
                             public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

    FixedDate(const DateTime& date);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    FixedDate(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    DateTime date;
};

///////////////////////////////////////////////////////////////////////////////
// IntervalOffset - Answers a date that is a specified interval (e.g. "3M")
//                  from the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL IntervalOffset: public CObject, 
                                 public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

    IntervalOffset(const MaturityPeriodConstSP& maturityPeriod);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    IntervalOffset(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    MaturityPeriodConstSP maturityPeriod;
};

///////////////////////////////////////////////////////////////////////////////
// NthDayOffset - Answers a date that is the nth day of the month of the base
//                date (or throws an exception if n is not valid for that
//                month).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL NthDayOffset: public CObject, 
                               public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

    NthDayOffset(int n);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    virtual void validatePop2Object();

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    int n;
};

///////////////////////////////////////////////////////////////////////////////
// LastDayOffset - Answers the date of the last day of the month of the base
//                 date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL LastDayOffset: public CObject, 
                                public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

    LastDayOffset(bool ignoreLeapYears = false);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    bool ignoreLeapYears;
};

///////////////////////////////////////////////////////////////////////////////
// WeekdayOffset - Abstract base class for DateAdjustments that attempt to
//                 answer the date of a specified day of the week (e.g. Monday)
//                 relative to the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL WeekdayOffset: public CObject, 
                                public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

protected:
    WeekdayOffset(const CClassConstSP& type, DateTime::DayOfWeek weekday);

    DateTime::DayOfWeek weekday;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// NthWeekdayOffset - Attempts to answer the date of the nth specified day of
//                    the week (e.g. 3rd Wednesday) of the month of the base
//                    date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL NthWeekdayOffset: public WeekdayOffset {
public:
    static CClassConstSP const TYPE;

    NthWeekdayOffset(int n, DateTime::DayOfWeek weekday);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    virtual void validatePop2Object();

    int n;
};

///////////////////////////////////////////////////////////////////////////////
// LastWeekdayOffset - Answers the date of the last specified day of the week
//                     (e.g. last Wednesday) of the month of the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL LastWeekdayOffset: public WeekdayOffset {
public:
    static CClassConstSP const TYPE;

    LastWeekdayOffset(DateTime::DayOfWeek weekday);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// OnOrBeforeWeekdayOffset - Answers the date of the specified day of the week
//                           that is on or before the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL OnOrBeforeWeekdayOffset: public WeekdayOffset {
public:
    static CClassConstSP const TYPE;

    OnOrBeforeWeekdayOffset(DateTime::DayOfWeek weekday);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// OnOrAfterWeekdayOffset - Answers the date of the specified day of the week
//                          that is on or after the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL OnOrAfterWeekdayOffset: public WeekdayOffset {
public:
    static CClassConstSP const TYPE;

    OnOrAfterWeekdayOffset(DateTime::DayOfWeek weekday);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// SameWeekWeekdayOffset - Answers the date of the specified day of the week
//                         that is in the same week as the base date.
//                         NOTE:  Weeks are considered to begin on Monday and
//                                end the following Sunday.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL SameWeekWeekdayOffset: public WeekdayOffset {
public:
    static CClassConstSP const TYPE;

    SameWeekWeekdayOffset(DateTime::DayOfWeek weekday);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// BusDateAdjustment - Abstract base class for DateAdjustments that operate
//                     in the context of a business calendar (i.e. Holiday).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL BusDateAdjustment : public CObject,
                                     public virtual IDateAdjustment {
public:
    static CClassConstSP const TYPE;

protected:
    BusDateAdjustment(const CClassConstSP& type,
                      const HolidayWrapper& holidays);

    BusDateAdjustment(const CClassConstSP& type);

    mutable HolidayWrapper holidays;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// BusDaysOffset - Answers a date that is a specified number of business days
//                 from the base date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL BusDaysOffset: public BusDateAdjustment {
public:
    static CClassConstSP const TYPE;

    BusDaysOffset(int n, const HolidayWrapper& holidays);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    BusDaysOffset(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    int n;
};

///////////////////////////////////////////////////////////////////////////////
// NthBusDayOffset - Answers the nth business day of the month of the base date
//                   (or throws an exception if no such day exists).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL NthBusDayOffset: public BusDateAdjustment {
public:
    static CClassConstSP const TYPE;

    NthBusDayOffset(int n, const HolidayWrapper& holidays);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    NthBusDayOffset(const CClassConstSP& type);

    virtual void validatePop2Object();

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    int n;
};

///////////////////////////////////////////////////////////////////////////////
// LastBusDayOffset - Answers the last business day of the month of the base
//                    date (or throws an exception if no such day exists).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL LastBusDayOffset: public BusDateAdjustment {
public:
    static CClassConstSP const TYPE;

    LastBusDayOffset(const HolidayWrapper& holidays, bool ignoreLeapYears = false);

    virtual DateTime dateFor(const DateTime& baseDate) const;

private:
    LastBusDayOffset(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    bool ignoreLeapYears;
};

///////////////////////////////////////////////////////////////////////////////
// BadDayAdjustment - Answers the base date if it is a good business day.
//                    Otherwise, the specified bad day convention is applied
//                    to produce an appropriate good business day.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL BadDayAdjustment : public BusDateAdjustment {
public:
    static CClassConstSP const TYPE;

    enum ConventionName { None, Previous, Following, ModifiedFollowing };

    BadDayAdjustment(ConventionName convention, const HolidayWrapper& holidays);

    virtual DateTime dateFor(const DateTime& baseDate) const;

protected:
    BadDayAdjustment(const CClassConstSP& type);

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    ConventionName convention;
};

DRLIB_END_NAMESPACE

#endif
