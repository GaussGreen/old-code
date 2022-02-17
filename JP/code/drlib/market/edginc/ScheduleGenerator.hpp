//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ScheduleGenerator.hpp
//
//   Description : Various classes to support schedules and schedule
//                 generation.
//
//   Author      : Jerry Cohen
//
//   Date        : November 2006
//
//----------------------------------------------------------------------------

#ifndef ScheduleGenerator_HPP
#define ScheduleGenerator_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/DateAdjustment.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// Experimental macros, to be moved somewhere else
///////////////////////////////////////////////////////////////////////////////
#define DECLARE_WITH_INLINE_ARRAY(T)                         \
    typedef smartConstPtr<T> T##ConstSP;                \
    typedef smartPtr<T> T##SP;                          \
    typedef array<T, T> T##Array;                       \
    typedef smartPtr<T##Array> T##ArraySP;              \
    typedef smartConstPtr<T##Array> T##ArrayConstSP;

#define DECLARE_arrayObjectCast(X) \
template <> class MARKET_DLL arrayObjectCast<X>{ \
public: \
    /** Casts array element to an IObject */ \
    static IObjectConstSP toIObject(const X& value); \
\
    /** Casts array element to an IObject */ \
    static IObjectSP toIObject(X& value); \
\
    /** Turns the IObjectSP into a X */ \
    static const X& fromIObject(IObjectSP& value); \
};

///////////////////////////////////////////////////////////////////////////////
// ISched - Top-level interface for all schedules.
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(ISched);

class MARKET_DLL ISched : public CObject {
public:
    static CClassConstSP const TYPE;

    virtual ~ISched();

    // The main schedule generation interface.
    virtual ISchedConstSP schedule() const = 0;

    // Subclasses can override to provide a flattened schedule that is more
    // easily consumed by certain clients (e.g. spreadsheets).
    virtual ISchedConstSP flatSchedule() const;

    virtual void setup() const = 0;

protected:
    ISched(const CClassConstSP& type);

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// ScheduleGenerator - An executable object that allows clients to generate
//                     a specified schedule with a specified market (for
//                     Holidays).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL ScheduleGenerator : public CObject,
                                     public virtual ClientRunnable {
public:
    static CClassConstSP const TYPE;

    virtual IObjectSP run();

protected:
    ScheduleGenerator(const CClassConstSP& type);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    static IObjectSP runAddin(ScheduleGenerator*);

    ISchedConstSP sched;
    MarketDataConstSP market;
    bool flatten;
};

///////////////////////////////////////////////////////////////////////////////
// IDateSched - Interface for classes that produce a schedule of dates (i.e. a
//              DateTimeArray).
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(IDateSched);

class MARKET_DLL IDateSched : public ISched {
public:
    static CClassConstSP const TYPE;
    
    virtual ISchedConstSP schedule() const;

    // Populates cache if needed.
    DateTimeArrayConstSP elements() const;

    // ASSUMES that one has previously called DateSched().
    const DateTimeArray& cachedElems() const { return *cachedElemsSP; }
    
protected:
    IDateSched(const CClassConstSP& type);

    mutable DateTimeArrayConstSP cachedElemsSP;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// DateSched - An explicitly specified date schedule where dates are expected
//             to be in increasing order, but duplicates are allowed.
//             For performance, this condition is only enforced when DateSched
//             objects are constructed via pop2Object.
///////////////////////////////////////////////////////////////////////////////
class DateSched;

class MARKET_DLL DateSched: public IDateSched {
public:
    static CClassConstSP const TYPE;

    DateSched(const DateTimeArrayConstSP& explicitDates);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    DateSched(const CClassConstSP& type);

    virtual void validatePop2Object();
    
    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    DateTimeArrayConstSP explicitDates;
};

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
class DateSchedGen;

class MARKET_DLL DateSchedGen: public IDateSched {
public:
    static CClassConstSP const TYPE;

    DateSchedGen(
        const DateTime& startDate,
        const DateTime& endDate,
        const MaturityPeriodConstSP& maturityPeriod,
        const DateTime& frontStubEndDate,
        const DateTime& backStubStartDate,
        bool includeStart = true,
        bool includeEnd = true,
        bool rollOnLast = false,
        bool ignoreLeapYears = false,
        int rollDay = 0);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

    // Answers true if the two dates are not separated by one whole period (i.e. MaturityPeriod).
    bool isStubPeriod(const DateTime& start, const DateTime& end) const;

private:
    DateSchedGen(const CClassConstSP& type);

    const DateTime& periodicBaseDate() const;
    const DateTime& periodicCutoffDate() const;

    void setEffectiveRollDay();
    DateTime rollDayAdjustedDate(const DateTime& srcDate) const;

    virtual void setup() const;

    virtual void validatePop2Object();

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    DateTime                startDate;
    DateTime                endDate;
    MaturityPeriodConstSP   maturityPeriod;
    DateTime                frontStubEndDate;
    DateTime                backStubStartDate;
    bool                    includeStart;
    bool                    includeEnd;
    bool                    rollOnLast;
    bool                    ignoreLeapYears;
    int                     rollDay;

    // Transient
    int                     periodCount;
    string                  periodInterval;
    int                     effectiveRollDay;
    bool                    useExplicitRollDay;
};

DECLARE(DateSchedGen);

///////////////////////////////////////////////////////////////////////////////
// AdjustedDateSched - Answers a new date schedule that is the result of
//                     applying a specified date adjustment to each of the
//                     dates of a source date schedule.
///////////////////////////////////////////////////////////////////////////////
class AdjustedDateSched;

class MARKET_DLL AdjustedDateSched: public IDateSched {
public:
    static CClassConstSP const TYPE;

    AdjustedDateSched(const IDateSchedConstSP& dateSched,
                      const IDateAdjustmentConstSP& dateAdjustment,
                      bool adjustFirst = true,
                      bool adjustLast = true,
                      bool includeDuplicates = false);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    AdjustedDateSched(const CClassConstSP& type);

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IDateSchedConstSP dateSched;
    IDateAdjustmentConstSP dateAdjustment;
    bool adjustFirst;
    bool adjustLast;
    bool includeDuplicates;
};

///////////////////////////////////////////////////////////////////////////////
// DateSchedUnion - Answers the union of zero or more date schedules with the
//                  order of dates preserved.
///////////////////////////////////////////////////////////////////////////////
class DateSchedUnion;

class MARKET_DLL DateSchedUnion: public IDateSched {
public:
    static CClassConstSP const TYPE;

    DateSchedUnion(const IDateSchedArrayConstSP& dateScheds);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    DateSchedUnion(const CClassConstSP& type);

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IDateSchedArrayConstSP  dateScheds;
};

///////////////////////////////////////////////////////////////////////////////
// Period - Represents an explicitly specified date/time interval.  The isStub
//          attribute can be used to mark periods that are of irregular length.
//          The period start date is expected to be less than or equal to the
//          period end date.  For performance, this condition is only enforced
//          when DateSched objects are constructed via pop2Object.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL Period: public CObject {
public:
    static CClassConstSP const TYPE;

    enum Boundary { None, Start, End };

    Period(const DateTime& start, const DateTime& end, bool isStub = false) :
        CObject(TYPE), start(start), end(end), isStub(isStub) {}

    bool operator ==(const Period& period) const;
    bool operator !=(const Period& period) const;

    string toString() const;

    DateTime start;
    DateTime end;
    bool isStub;

protected:
    Period(const CClassConstSP& type);

    virtual void validatePop2Object();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(Period);

///////////////////////////////////////////////////////////////////////////////
// IPeriodSched - Interface for classes that answer a schedule of periods (i.e.
//                a PeriodArray.
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(IPeriodSched);

class MARKET_DLL IPeriodSched : public ISched {
public:
    static CClassConstSP const TYPE;

    virtual ISchedConstSP schedule() const;

    // Populates cache if needed.
    PeriodArrayConstSP elements() const;

    // ASSUMES one has previously called DateSched().
    const PeriodArray& cachedElems() const { return *cachedElemsSP; }

protected:
    IPeriodSched(const CClassConstSP& type);

    mutable PeriodArrayConstSP cachedElemsSP;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// DateSchedFromPeriodSched - Answers a new date schedule that is the result of
//                            selecting either the period start or end dates
//                            from a source period schedule.  Optionally, one
//                            may specify a date adjustment to apply to the
//                            selected dates.
///////////////////////////////////////////////////////////////////////////////
class DateSchedFromPeriodSched;

class MARKET_DLL DateSchedFromPeriodSched : public IDateSched {
public:
    static CClassConstSP const TYPE;

    DateSchedFromPeriodSched(const IPeriodSchedConstSP& periodSched, Period::Boundary startOrEnd,
                             const IDateAdjustmentConstSP& dateAdjustment);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    virtual void validatePop2Object();

    virtual void setup() const;

    DateSchedFromPeriodSched(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP     periodSched;
    Period::Boundary        startOrEnd;
    IDateAdjustmentConstSP  dateAdjustment;
};

///////////////////////////////////////////////////////////////////////////////
// PeriodSched - An explicitly specified period schedule where periods are
//               expected to be in strictly increasing order such that:
//                   For all i: period[i].start < period[i+1].start and
//                              period[i].end < period[i+1].end)
//               However, periods MAY overlap such that:
//                   period[i+1].start < period[i].end
//               For performance, we only check for increasing periods when
//               PeriodSched objects are constructed via pop2Object.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL PeriodSched: public IPeriodSched {
public:
    static CClassConstSP const TYPE;

    PeriodSched(const PeriodArrayConstSP& explicitPeriods);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    PeriodSched(const CClassConstSP& type);

    virtual void validatePop2Object();
    
    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    PeriodArrayConstSP explicitPeriods;
};

///////////////////////////////////////////////////////////////////////////////
// PeriodSchedGen - Generates a schedule of periods from the consective dates
//                  that are produced by the specified DateSchedGen.  The first
//                  and last periods may be marked as stubs if the first period
//                  start date and last period end date are not periodic dates,
//                  respectively.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL PeriodSchedGen: public IPeriodSched {
public:
    static CClassConstSP const TYPE;

    PeriodSchedGen(const DateSchedGenConstSP& dateSchedGen);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    PeriodSchedGen(const CClassConstSP& type);
    
    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    DateSchedGenConstSP    dateSchedGen;
};

///////////////////////////////////////////////////////////////////////////////
// AdjustedPeriodSched - Generates a schedule of periods that is the result of
//                       applying a date adjustment to each period's start and
//                       end date.  One can independently specify whether to
//                       adjust the first period start date and the last period
//                       end date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL AdjustedPeriodSched: public IPeriodSched {
public:
    static CClassConstSP const TYPE;

    AdjustedPeriodSched(const IPeriodSchedConstSP& periodSched,
                        const IDateAdjustmentConstSP& dateAdjustment,
                        bool adjustFirst = true,
                        bool adjustLast = true);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    AdjustedPeriodSched(const CClassConstSP& type);
    
    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP periodSched;
    IDateAdjustmentConstSP dateAdjustment;
    bool adjustFirst;
    bool adjustLast;
};

///////////////////////////////////////////////////////////////////////////////
// PeriodObsDates - Associates an array of (observation) dates with a Period.
//                  The observation dates are expected to be in increasing
//                  order, but this condition is only enforced when
//                  PeriodObsDates are constructed via pop2Object.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL PeriodObsDates: public CObject {
public:
    static CClassConstSP const TYPE;

    PeriodObsDates(const PeriodConstSP& period, const DateTimeArrayConstSP& obsDates);
    PeriodObsDates(const PeriodConstSP& period, const DateTime& obsDate);

    PeriodConstSP           period;
    DateTimeArrayConstSP    obsDates;

protected:
    PeriodObsDates(const CClassConstSP& type);

    virtual void validatePop2Object();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(PeriodObsDates);

///////////////////////////////////////////////////////////////////////////////
// FlatPeriodObsDates - A "flattened" version of PeriodObsDates.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FlatPeriodObsDate : public CObject {
public:
    static CClassConstSP const TYPE;

    FlatPeriodObsDate(const DateTime& periodStart,
                      const DateTime& periodEnd,
                      bool isStub,
                      const DateTime& obsDate);

    DateTime        periodStart;
    DateTime        periodEnd;
    bool            isStub;
    DateTime        obsDate;

protected:
    FlatPeriodObsDate(const CClassConstSP& type);

    virtual void validatePop2Object();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(FlatPeriodObsDate);

///////////////////////////////////////////////////////////////////////////////
// IObsSched - Interface for classes that represent a schedule of
//             PeriodObsDates.
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(IObsSched);

class MARKET_DLL IObsSched : public ISched {
public:
    static CClassConstSP const TYPE;

    virtual ISchedConstSP schedule() const;
    virtual ISchedConstSP flatSchedule() const;

    // Populates cache if needed.
    PeriodObsDatesArrayConstSP elements() const;

    // ASSUMES one has previously called DateSched().
    const PeriodObsDatesArray& cachedElems() const { return *cachedElemsSP; }

protected:
    IObsSched(const CClassConstSP& type);

    mutable PeriodObsDatesArrayConstSP cachedElemsSP;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// ObsSched - An explicitly specified schedule of PeriodObsDates.  The periods
//            in the PeriodObsDatesArray are expected to appear in increasing
//            order like with an explicit PeriodSched.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL ObsSched: public IObsSched {
public:
    static CClassConstSP const TYPE;

    ObsSched(const PeriodObsDatesArrayConstSP& explicitPeriodObsDates);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    ObsSched(const CClassConstSP& type);

    virtual void validatePop2Object();
    
    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    PeriodObsDatesArrayConstSP explicitPeriodObsDates;
};

///////////////////////////////////////////////////////////////////////////////
// FlatObsSched - An explicitly specified schedule of FlatPeriodObsDates.
//                Periods are expected to be specified in increasing order.
//                Reset observation dates for the same period are also expected
//                to be in increasing order.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FlatObsSched: public IObsSched {
public:
    static CClassConstSP const TYPE;

    FlatObsSched(const FlatPeriodObsDateArrayConstSP& explicitFlatPeriodObsDates);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

    virtual ISchedConstSP flatSchedule() const;

private:
    FlatObsSched(const CClassConstSP& type);

    virtual void validatePop2Object();

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    FlatPeriodObsDateArrayConstSP explicitFlatPeriodObsDates;
};

///////////////////////////////////////////////////////////////////////////////
// RelativeObsSched - Generates a schedule of PeriodObsDates where there is one
//                    observation date for each period, and each observation
//                    date is either the period start or end date plus an
//                    optional date adjustment.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL RelativeObsSched: public IObsSched {
public:
    static CClassConstSP const TYPE;

    RelativeObsSched(const IPeriodSchedConstSP& periodSched,
                     Period::Boundary relativeToDate,
                     const IDateAdjustmentConstSP& obsAdjustment);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    virtual void validatePop2Object();

    virtual void setup() const;

    RelativeObsSched(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP     periodSched;
    Period::Boundary        relativeToDate;
    IDateAdjustmentConstSP  obsAdjustment;
};

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
class MARKET_DLL WindowedObsSched: public IObsSched {
public:
    static CClassConstSP const TYPE;

    WindowedObsSched(const IPeriodSchedConstSP& periodSched,
                     const IDateSchedConstSP& obsSched,
                     const IDateAdjustmentConstSP& obsAdjustment,
                     Period::Boundary selectionStartRelativeToDate,
                     const IDateAdjustmentConstSP& selectionStartOffset,
                     Period::Boundary selectionEndRelativeToDate,
                     const IDateAdjustmentConstSP& selectionEndOffset,
                     Period::Boundary obsCutoffRelativeToDate,
                     const IDateAdjustmentConstSP& obsCutoffOffset,
                     bool  repeatLastObs);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    virtual void validatePop2Object();

    virtual void setup() const;

    WindowedObsSched(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP     periodSched;
    IDateSchedConstSP       obsSched;
    IDateAdjustmentConstSP  obsAdjustment;
    Period::Boundary        selectionStartRelativeToDate;
    IDateAdjustmentConstSP  selectionStartOffset;
    Period::Boundary        selectionEndRelativeToDate;
    IDateAdjustmentConstSP  selectionEndOffset;
    Period::Boundary        obsCutoffRelativeToDate;
    IDateAdjustmentConstSP  obsCutoffOffset;
    bool                    repeatLastObs;
};

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
class MARKET_DLL OnOrBeforeObsSched: public IObsSched {
public:
    static CClassConstSP const TYPE;

    OnOrBeforeObsSched(const IPeriodSchedConstSP& periodSched,
                       const IDateSchedConstSP& obsSched,
                       Period::Boundary cutoffRelativeToDate,
                       const IDateAdjustmentConstSP& cutoffOffset);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    virtual void validatePop2Object();

    virtual void setup() const;

    OnOrBeforeObsSched(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP     periodSched;
    IDateSchedConstSP       obsSched;
    Period::Boundary        cutoffRelativeToDate;
    IDateAdjustmentConstSP  cutoffOffset;
};

///////////////////////////////////////////////////////////////////////////////
// OptionDates - Object to hold an option exercise date and associated
//               notification date.
///////////////////////////////////////////////////////////////////////////////

class MARKET_DLL OptionDates: public CObject {
public:
    static CClassConstSP const TYPE;

    OptionDates(const CClassConstSP& type=TYPE);
    OptionDates(const DateTime& exerciseDate, const DateTime& notificationDate);

    DateTime exerciseDate;
    DateTime notificationDate;

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE_WITH_INLINE_ARRAY(OptionDates);
DECLARE_arrayObjectCast(OptionDates);

///////////////////////////////////////////////////////////////////////////////
// IOptionSched - Interface for classes that generate a schedule of
//                OptionDates.
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(IOptionSched);

class MARKET_DLL IOptionSched : public ISched {
public:
    static CClassConstSP const TYPE;

    void check(bool ensureStrictOrdering) const;

    // access the schedule as an array of OptionDates
    OptionDates const operator[](int i) const { return optionDatesArray[i]; }
    int size()       const { return optionDatesArray.size(); }

    // to access the schedule with iterators
    typedef OptionDatesArray::const_iterator Iterator;
    Iterator begin() const { return optionDatesArray.begin(); }
    Iterator end()   const { return optionDatesArray.end(); }

    // to access the schedule colomn-wise
    DateTimeArrayConstSP exerciseDates()     const { return colExerciseDateSP; }
    DateTimeArrayConstSP notificationDates() const { return colNotificationDateSP; }

protected:
    IOptionSched(const CClassConstSP& type);
    void fillColumns() const; // to be called by the child's setup()

    // the computed schedule, to be fill by the child's setup()
    mutable OptionDatesArray optionDatesArray;

private:
    // the columns, computed from optionDatesArray by fillColumns()
    mutable DateTimeArraySP colExerciseDateSP;
    mutable DateTimeArraySP colNotificationDateSP;

    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// OptionSched - An explicitly specified schedule of OptionDates that are
//               expected to be specified in increasing order such that:
//                   For all i: exercise[i] < exercise[i+1] and
//                              notification[i] < notification[i+1]
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL OptionSched: public IOptionSched {
public:
    static CClassConstSP const TYPE;

    virtual void setup() const;
    virtual ISchedConstSP schedule() const;
    OptionSched(const OptionDatesArrayConstSP& explicitOptionDates);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    OptionSched(const CClassConstSP& type);
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    OptionDatesArrayConstSP explicitOptionDates;
};

///////////////////////////////////////////////////////////////////////////////
// OptionSchedGen - Generates a schedule of OptionDates from a specified date
//                  schedule.  Separate date adjustments are applied to produce
//                  each exercise and notification date, and the notification
//                  date adjustment may be applied to either the unadjusted
//                  exercise date (i.e. input date schedule) or the adjusted
//                  exercise date.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL OptionSchedGen: public IOptionSched {
public:
    static CClassConstSP const TYPE;

    enum OptionDate { Exercise, UnadjustedExercise };

    virtual void setup() const;

    virtual ISchedConstSP schedule() const;

    OptionSchedGen(const IDateSchedConstSP& dateSched,
                   const IDateAdjustmentConstSP& exerciseAdjustment,
                   OptionDate notificationRelativeToDate,
                   const IDateAdjustmentConstSP& notificationAdjustment);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    OptionSchedGen(const CClassConstSP& type);
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IDateSchedConstSP       dateSched;
    IDateAdjustmentConstSP  exerciseAdjustment;
    OptionDate              notificationRelativeToDate;
    IDateAdjustmentConstSP  notificationAdjustment;
};

///////////////////////////////////////////////////////////////////////////////
// FixedCouponDates - Object to hold dates associated with a fixed-rate
//                    coupon payment.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FixedCouponDates: public CObject {
public:
    static CClassConstSP const TYPE;

    FixedCouponDates(const PeriodConstSP& accrualPeriod,
                     const DateTime& paymentDate);

    PeriodConstSP               accrualPeriod;
    DateTime                    paymentDate;

protected:
    FixedCouponDates(const CClassConstSP& type);
    FixedCouponDates(const CClassConstSP& type, const PeriodConstSP& accrualPeriod, const DateTime& paymentDate);

    virtual void validatePop2Object();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(FixedCouponDates);

///////////////////////////////////////////////////////////////////////////////
// FlatFixedCouponDates - A "flattened" version of FixedCouponDates.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FlatFixedCouponDates : public CObject {
public:
    static CClassConstSP const TYPE;

    FlatFixedCouponDates(const DateTime& accrualStart,
                         const DateTime& accrualEnd,
                         bool isStub,
                         const DateTime& paymentDate);

    DateTime    accrualStart;
    DateTime    accrualEnd;
    bool        isStub;
    DateTime    paymentDate;

protected:
    FlatFixedCouponDates(const CClassConstSP& type);
    FlatFixedCouponDates(const CClassConstSP& type,
                         const DateTime& accrualStart,
                         const DateTime& accrualEnd,
                         bool isStub,
                         const DateTime& paymentDate);

    virtual void validatePop2Object();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(FlatFixedCouponDates);

///////////////////////////////////////////////////////////////////////////////
// IFixedCouponSched - Interface for classes that generate a schedule of
//                     FixedCouponDates.
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(IFixedCouponSched);

class MARKET_DLL IFixedCouponSched : public ISched {
public:
    static CClassConstSP const TYPE;

    virtual ISchedConstSP schedule() const;
    virtual ISchedConstSP flatSchedule() const;

    // Populates cache if needed.
    FixedCouponDatesArrayConstSP elements() const;

    // ASSUMES one has previously called setup().
    const FixedCouponDatesArray& cachedElems() const { return *cachedElemsSP; }

protected:
    IFixedCouponSched(const CClassConstSP& type);

    mutable FixedCouponDatesArrayConstSP cachedElemsSP;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// FixedCouponSched - An explicitly specified schedule of FixedCouponDates.
//                    Accrual periods are expected to be specified in
//                    increasing order, but overlapping periods are okay.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FixedCouponSched: public IFixedCouponSched {
public:
    static CClassConstSP const TYPE;

    FixedCouponSched(const FixedCouponDatesArrayConstSP& explicitFixedCouponDates);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    FixedCouponSched(const CClassConstSP& type);

    virtual void validatePop2Object();

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    FixedCouponDatesArrayConstSP explicitFixedCouponDates;
};

///////////////////////////////////////////////////////////////////////////////
// FlatFixedCouponSched - An explicitly specified schedule of
//                        FlatFixedCouponDates.  Accrual periods are expected
//                        to be specified in increasing order, but overlapping
//                        periods are okay.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FlatFixedCouponSched: public IFixedCouponSched {
public:
    static CClassConstSP const TYPE;

    FlatFixedCouponSched(const FlatFixedCouponDatesArrayConstSP& explicitFlatFixedCouponDates);

    virtual ISchedConstSP flatSchedule() const;

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    FlatFixedCouponSched(const CClassConstSP& type);

    virtual void validatePop2Object();

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    FlatFixedCouponDatesArrayConstSP explicitFlatFixedCouponDates;
};

///////////////////////////////////////////////////////////////////////////////
// FixedCouponSchedGen - Generates a schedule of FixedCouponDates where accrual
//                       periods are given by an IPeriodSched and payment dates
//                       are specified relative to each period start or end
//                       date plus an optional date adjustment.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FixedCouponSchedGen: public IFixedCouponSched {
public:
    static CClassConstSP const TYPE;

    FixedCouponSchedGen(const IPeriodSchedConstSP& accrualSched,
                        Period::Boundary paymentRelativeToDate,
                        const IDateAdjustmentConstSP& paymentDateAdjustment);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    virtual void validatePop2Object();

    virtual void setup() const;

    FixedCouponSchedGen(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP     accrualSched;
    Period::Boundary        paymentRelativeToDate;
    IDateAdjustmentConstSP  paymentDateAdjustment;
};

///////////////////////////////////////////////////////////////////////////////
// FloatCouponDates - Object to hold dates associated with a floating-rate
//                    coupon payment.  Reset dates are expected to be in
//                    increasing order.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FloatCouponDates: public FixedCouponDates {
public:
    static CClassConstSP const TYPE;

    FloatCouponDates(const PeriodConstSP& accrualPeriod,
                     const DateTime& paymentDate,
                     const DateTimeArrayConstSP& resetDates);

    DateTimeArrayConstSP        resetDates;

protected:
    FloatCouponDates(const CClassConstSP& type);

    virtual void validatePop2Object();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(FloatCouponDates);

///////////////////////////////////////////////////////////////////////////////
// FlatFloatCouponDates - A "flattened" version of FloatCouponDates.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FlatFloatCouponDates : public FlatFixedCouponDates {
public:
    static CClassConstSP const TYPE;

    FlatFloatCouponDates(const DateTime& accrualStart,
                         const DateTime& accrualEnd,
                         bool isStub,
                         const DateTime& paymentDate,
                         const DateTime& resetDate);

    DateTime    resetDate;

protected:
    FlatFloatCouponDates(const CClassConstSP& type);

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
};

DECLARE(FlatFloatCouponDates);

///////////////////////////////////////////////////////////////////////////////
// IFloatCouponSched - Interface for classes that generate a schedule of
//                     FloatCouponDates.
///////////////////////////////////////////////////////////////////////////////
FORWARD_DECLARE(IFloatCouponSched);

class MARKET_DLL IFloatCouponSched : public ISched {
public:
    static CClassConstSP const TYPE;

    virtual ISchedConstSP schedule() const;
    virtual ISchedConstSP flatSchedule() const;

    // Populates cache if needed.
    FloatCouponDatesArrayConstSP elements() const;

    // ASSUMES one has previously called DateSched().
    const FloatCouponDatesArray& cachedElems() const { return *cachedElemsSP; }

protected:
    IFloatCouponSched(const CClassConstSP& type);

    mutable FloatCouponDatesArrayConstSP cachedElemsSP;

private:
    static void load(CClassSP& clazz);
};

///////////////////////////////////////////////////////////////////////////////
// FloatCouponSched - An explicitly specified schedule of FloatCouponDates.
//                    Accrual periods are expected to be in increasing order.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FloatCouponSched: public IFloatCouponSched {
public:
    static CClassConstSP const TYPE;

    FloatCouponSched(const FloatCouponDatesArrayConstSP& explicitFloatCouponDates);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    FloatCouponSched(const CClassConstSP& type);

    virtual void validatePop2Object();

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    FloatCouponDatesArrayConstSP explicitFloatCouponDates;
};

///////////////////////////////////////////////////////////////////////////////
// FlatFloatCouponSched - An explicitly specified schedule of
//                        FlatFloatCouponDates.  Accrual periods are expected
//                        to be in increasing order, and reset dates for a
//                        particular accrual period are also expected in
//                        increasing order.
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FlatFloatCouponSched: public IFloatCouponSched {
public:
    static CClassConstSP const TYPE;

    FlatFloatCouponSched(const FlatFloatCouponDatesArrayConstSP& explicitFlatFloatCouponDates);

    virtual ISchedConstSP flatSchedule() const;

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    FlatFloatCouponSched(const CClassConstSP& type);

    virtual void validatePop2Object();

    virtual void setup() const;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    FlatFloatCouponDatesArrayConstSP explicitFlatFloatCouponDates;
};

///////////////////////////////////////////////////////////////////////////////
// FloatCouponSchedGen - Generates a schedule of FloatCouponDates where accrual
//                       periods are given by an IPeriodSched and payment dates
//                       are specified relative to each period start or end
//                       date plus an optional date adjustment.  Reset
//                       observation dates are specified by an IObsSched that
//                       must produce PeriodObsDates that correspond to the
//                       accrual periods (both in number and order).
///////////////////////////////////////////////////////////////////////////////
class MARKET_DLL FloatCouponSchedGen: public IFloatCouponSched {
public:
    static CClassConstSP const TYPE;

    FloatCouponSchedGen(const IPeriodSchedConstSP& accrualSched,
                        Period::Boundary paymentRelativeToDate,
                        const IDateAdjustmentConstSP& paymentDateAdjustment,
                        const IObsSchedConstSP& resetSched);

    // These objects are intended to be immutable.  Override clone to support
    // caching.
    virtual IObject* clone() const;

private:
    virtual void validatePop2Object();

    virtual void setup() const;

    FloatCouponSchedGen(const CClassConstSP& type);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    IPeriodSchedConstSP     accrualSched;
    Period::Boundary        paymentRelativeToDate;
    IDateAdjustmentConstSP  paymentDateAdjustment;
    IObsSchedConstSP        resetSched;
};

DRLIB_END_NAMESPACE

#endif
