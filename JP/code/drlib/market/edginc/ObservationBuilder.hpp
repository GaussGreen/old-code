//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ObservationBuilder.hpp
//
//   Description : Parameterized Observation Builder
//
//   Author      : Manos Venardos
//
//   Date        : 23 February 2006
//
//
//----------------------------------------------------------------------------

#ifndef OBSERVATION_BUILDER_HPP
#define OBSERVATION_BUILDER_HPP

#include "edginc/DateBuilder.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/RevertTypeConvert.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for parameterized timeline */
class MARKET_DLL IObservationBuilder: virtual public IObject {
public:
   
    static CClassConstSP const TYPE;

    /** retrieve dates */
    virtual DateTimeArraySP dateList() const = 0;
    /** retrieve obs types */
    virtual ObservationTypeArraySP obsTypes() const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IObservationBuilder> IObservationBuilderSP;

// ObservationBuilderDaily builds a sample schedule to include all days except
// (optionally) weekends and a set of explicit excludeDates.
// Observation times are defaulted from ObservationTypes. 
class MARKET_DLL ObservationBuilderDaily:
    public CObject, 
    virtual public IObservationBuilder {
public:
    static CClassConstSP const TYPE;
    friend class ObservationBuilderDailyHelper;

    /** retrieve dates */
    virtual DateTimeArraySP dateList() const;
    /** retrieve obs types */
    virtual ObservationTypeArraySP obsTypes() const;
    virtual void validatePop2Object();
        
    // explicit constructor
    ObservationBuilderDaily(const DateTime& startDate,
                            const DateTime& endDate,
                            bool excludeWeekends,
                            const DateTimeArray& excludeDates,
                            const string& startObsType,
                            const string& obsType,
                            const string& endObsType);

private:
    DateTime startDate;          //!< Start date for timeline. Time ignored
    DateTime endDate;            //!< End date for timeline. Time ignored
    bool     excludeWeekends;    //!< True: exclude weekends, False: include weekends
    DateTimeArray excludeDates; //!< Dates to exclude, optional - defaults to empty date list
    HolidaySP       hols;        // Transient hols built from excludeWeekends and excludeDates
    
    string         startObsType;  //!< Observation type for start date
    string         obsType;       //!< Observation type for intermediate dates
    string         endObsType;    //!< Observation type for end date

    // transient
    DateTimeArray dates;
    ObservationTypeArray types;

    ObservationBuilderDaily();

protected:
    ObservationBuilderDaily(CClassConstSP clazz);
};

// SimpleObservationSchedule builds a sample schedule from an array of dates
// and an ObservationType.
// It's intended to be used as an intermediate for use when upgrading products
// with primitive sampling schedules to support more sophisticated 
// IObservationBuilder sampling schedules.
//
// It registers a method for converting DateTimeArray to IObservationBuilder
// (as a SimpleObservationBuilder wrapping the array and obsType 'NotUsed') with
// the reflection mechanism, to enable drop-in replacement of DateTimeArrays
// when upgrading products.
class SimpleObservationSchedule: public CObject,
                                 virtual public IObservationBuilder,
                                 virtual public IRevertTypeConvert {
public:
    static CClassConstSP const TYPE;
    friend class ObservationBuilderDailyHelper;

    /** retrieve dates */
    virtual DateTimeArraySP dateList() const;
    /** retrieve obs types */
    virtual ObservationTypeArraySP obsTypes() const;

    // Conversion function for DateTimeArray.
    static IObjectSP fromDateTimeArray(const IObjectSP& object, 
                                       CClassConstSP    requiredType);

    // <<IRevertTypeConvert>> interface. Revert object to type for interface.
    virtual IObjectSP revert(const string& interfaceType) const;

    // Explicit constructor
    SimpleObservationSchedule(IDateBuilderSP theDateBuilder,
                              ObservationTypeSP theObsType,
                              const string& theExchange);

private:
    IDateBuilderSP    dates;    //!< Observation dates
    ObservationTypeSP obsType;  //!< Type of observation for this schedule
    string            source;   //!< Name of the exchange

    static void load(CClassSP& clazz);
    static IObject* defaultSimpleObservationSchedule();

    SimpleObservationSchedule();

protected:
    SimpleObservationSchedule(CClassConstSP clazz);
};

DRLIB_END_NAMESPACE

#endif
