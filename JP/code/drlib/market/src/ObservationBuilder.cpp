//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ObservationBuilder.cpp
//
//   Description : Observation Builder
//
//   Author      : Manos Venardos
//
//   Date        : 23 February 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE


void IObservationBuilder::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IObservationBuilder, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}


CClassConstSP const IObservationBuilder::TYPE = CClass::registerInterfaceLoadMethod(
    "IObservationBuilder", typeid(IObservationBuilder), IObservationBuilder::load);


/////////////////////////////////////////////////////////////////////////////////////


/** Addin class wrapping date builder */
class ObservationBuilderAddin: public CObject {
public:
    static CClassConstSP const TYPE;
    
    /** Creates a DateTimeArray */
    DateTimeArraySP buildObsDates() {
        return observationBuilder->dateList();
    }

    /** Creates a DateTimeArray */
    ObservationTypeArraySP buildObsTypes() {
        return observationBuilder->obsTypes();
    }

private:
    /** Registration */
    static void load(CClassSP& clazz) {
        REGISTER(ObservationBuilderAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObservationBuilderAddin);
        FIELD(observationBuilder, "Observation builder");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        
        Addin::registerObjectMethod("BUILD_OBSERVATION_DATES",
                                    Addin::RISK,
                                    "Produces date list from date builder",
                                    false,
                                    Addin::expandSimple,
                                    &ObservationBuilderAddin::buildObsDates);

        Addin::registerObjectMethod("BUILD_OBSERVATION_TYPES",
                                    Addin::RISK,
                                    "Produces list of obs types from date builder",
                                    false,
                                    Addin::expandSimple,
                                    &ObservationBuilderAddin::buildObsTypes);
    }

    /** Empty shell method */
    static IObject* defaultObservationBuilderAddin() {
        return new ObservationBuilderAddin();
    }

    /** Empty constructor */
    ObservationBuilderAddin(): CObject(TYPE) {}
    
    IObservationBuilderSP observationBuilder;     //!< Observation builder object
};

CClassConstSP const ObservationBuilderAddin::TYPE = CClass::registerClassLoadMethod(
    "ObservationBuilderAddin", typeid(ObservationBuilderAddin), load);


/////////////////////////////////////////////////////////////////////////////////////


/** Parameterized date sequence with observation types */
class ObservationBuilder: public CObject, virtual public IObservationBuilder {
public:
    static CClassConstSP const TYPE;

    /** Get date list */
    virtual DateTimeArraySP dateList() const {
        static const string routine = "ObservationBuilder::dates";
        try {
            // Build dates
            return dateBuilder->dates();
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Get list of obs types */
    virtual ObservationTypeArraySP obsTypes() const {
        static const string routine = "ObservationBuilder::obsTypes";
        try {
            ObservationTypeArraySP copyTypes(copy(&types));
            return copyTypes;
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    virtual void validatePop2Object() {
        static const string routine = "ObservationBuilder::validatePop2Object";
        try {
            ObservationTypeSP interm = ObservationType::make(obsType);
            types.push_back(ObservationType::make(startObsType));
            for (int i = 0; i < dateBuilder->size()-2; ++i) {
                types.push_back(interm);
            }
            types.push_back(ObservationType::make(endObsType));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    IDateBuilderSP dateBuilder;   //!< Date builder
    
    string         startObsType;  //!< Observation type for start date
    string         endObsType;    //!< Observation type for end date
    string         obsType;       //!< Observation type for intermediate dates

    //transient
    ObservationTypeArray types;
    
    static void load(CClassSP& clazz) {
        REGISTER(ObservationBuilder, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IObservationBuilder);
        EMPTY_SHELL_METHOD(defaultObservationBuilder);
        FIELD(dateBuilder, "Date builder");
        FIELD(startObsType, "Observation type for start date");
        FIELD(endObsType, "Observation type for end date");
        FIELD(obsType, "Observation type for intermediate dates");
        clazz->setPublic(); // make visible to EAS/spreadsheet        
    }

    static IObject* defaultObservationBuilder() {
        return new ObservationBuilder();
    }

    ObservationBuilder(): CObject(TYPE), types(0) {}

protected:
    ObservationBuilder(CClassConstSP clazz): 
         CObject(clazz) {}
};

CClassConstSP const ObservationBuilder::TYPE = CClass::registerClassLoadMethod(
    "ObservationBuilder", typeid(ObservationBuilder), ObservationBuilder::load);

// For IMS purposes only
// ObservationBuilder shell containing DailyDateBuilder concrete class
// DailyObservationBuider has been superceded by ObservationBuilderDaily.
// To be deleted once removed from IMS
class DailyObservationBuilder: public ObservationBuilder {
public:
    static CClassConstSP const TYPE; 

private:
    DailyObservationBuilder():ObservationBuilder(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(DailyObservationBuilder, clazz);
        SUPERCLASS(ObservationBuilder);
        EMPTY_SHELL_METHOD(defaultDailyObservationBuilder);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultDailyObservationBuilder(){
        return new DailyObservationBuilder();
    }
};

CClassConstSP const DailyObservationBuilder::TYPE =
CClass::registerClassLoadMethod(
    "DailyObservationBuilder", typeid(DailyObservationBuilder), load);

// ObservationBuilderDaily builds a sample schedule to include all days except
// (optionally) weekends and a set of explicit excludeDates.
// Observation times are defaulted from ObservationTypes. 

/** get date list */
DateTimeArraySP ObservationBuilderDaily::dateList() const {
    static const string routine = "ObservationBuilderDaily::dates";
    try {
        DateTimeArraySP copyDates(copy(&dates));
        return copyDates;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** get obs type list */
ObservationTypeArraySP ObservationBuilderDaily::obsTypes() const {
    static const string routine = "ObservationBuilderDaily::obsTypes";
    try {
        ObservationTypeArraySP copyTypes(copy(&types));
        return copyTypes;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

void ObservationBuilderDaily::validatePop2Object() {
    static const string routine = "ObservationBuilderDaily::validatePop2Object";
    try {

        // Build a holiday object from excludeWeekends and excludeDates
        hols = HolidaySP(new Holiday("DateBuilder holidays",
                                        excludeDates,
                                        excludeWeekends));

        // Validate against singleton timeline and reverse order
        if(startDate >= endDate) {
            throw ModelException(routine, 
                "Start date " + startDate.toString() + 
                " has to be strictly before end date " + 
                endDate.toString());
        }

        // Check that start and end dates are not holidays
        if(hols->isHoliday(startDate)) {
            throw ModelException(routine, 
                "Start date " + startDate.toString() + 
                " is an excluded date");
        }

        if(hols->isHoliday(endDate)) {
            throw ModelException(routine, 
                "End date " + endDate.toString() + 
                " is an excluded date");
        }

        // Build types and dates
        ObservationTypeSP startType = ObservationType::make(startObsType);
        ObservationTypeSP intermType = ObservationType::make(obsType);
        ObservationTypeSP endType = ObservationType::make(endObsType);

        // Add start date
        dates.push_back(DateTime(startDate.getDate(), 
                    DateTime::timeConvert(startType->indicativeTime())));
        types.push_back(startType);

        // Add intermediate dates
        DateTime date(startDate.getDate() + 1, 
                    DateTime::timeConvert(intermType->indicativeTime()));

        while (date.getDate() < endDate.getDate()) {
            // Add only if not an excluded date
            if (hols->isBusinessDay(date)) {
                dates.push_back(date);
                types.push_back(intermType);
            }
            date = date.rollDate(1);
        }

        // Add last date
        dates.push_back(DateTime(endDate.getDate(),
                    DateTime::timeConvert(endType->indicativeTime())));
        types.push_back(endType);

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

ObservationBuilderDaily::ObservationBuilderDaily(const DateTime& startDate,
                        const DateTime& endDate,
                        bool excludeWeekends,
                        const DateTimeArray& excludeDates,
                        const string& startObsType,
                        const string& obsType,
                        const string& endObsType): CObject(TYPE), startDate(startDate),
                                endDate(endDate), excludeWeekends(excludeWeekends),
                                excludeDates(excludeDates), startObsType(startObsType),
                                obsType(obsType), endObsType(endObsType),
                                dates(0), types(0)
{
    validatePop2Object();
}

ObservationBuilderDaily::ObservationBuilderDaily(): 
                CObject(TYPE), 
                excludeWeekends(true), excludeDates(0),
                dates(0), types(0) {}

ObservationBuilderDaily::ObservationBuilderDaily(CClassConstSP clazz): 
         CObject(clazz), 
         excludeWeekends(true), excludeDates(0),
         dates(0), types(0)  {}

class ObservationBuilderDailyHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObservationBuilderDaily, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IObservationBuilder);
        EMPTY_SHELL_METHOD(defaultObservationBuilderDaily);
        FIELD(startDate, "Start date");
        FIELD(endDate, "End date");
        FIELD(excludeWeekends, "Exclude weekends");
        FIELD_MAKE_OPTIONAL(excludeWeekends);
        FIELD(excludeDates, "Explicit dates to exclude");
        FIELD_MAKE_OPTIONAL(excludeDates);
        FIELD(hols, "Holidays");
        FIELD_MAKE_TRANSIENT(hols);
        FIELD(startObsType, "Observation type for start date");
        FIELD(obsType, "Observation type for intermediate dates");
        FIELD(endObsType, "Observation type for end date");
        FIELD(dates, "Dates");
        FIELD_MAKE_TRANSIENT(dates);
        FIELD(types, "Observation types");
        FIELD_MAKE_TRANSIENT(types);
        clazz->setPublic(); // make visible to EAS/spreadsheet        
    }
    
    static IObject* defaultObservationBuilderDaily(){
        return new ObservationBuilderDaily();
    }
};

CClassConstSP const ObservationBuilderDaily::TYPE =
CClass::registerClassLoadMethod(
         "ObservationBuilderDaily", typeid(ObservationBuilderDaily), ObservationBuilderDailyHelper::load);

//////////////////////////////////////////////////////////////////:w
///////////////////

// <<IObservationBuilder>> interface. get date list.
DateTimeArraySP SimpleObservationSchedule::dateList() const {
    static const string routine = "SimpleObservationBuilder::dates";
        try {
            return dates->dates();
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
}

// <<IObservationBuilder>> interface. get obs type list.
ObservationTypeArraySP SimpleObservationSchedule::obsTypes() const {
    static const string routine = "SimpleObservationBuilder::dates";
        try {
            ObservationTypeArraySP types(new ObservationTypeArray(dates->size(), obsType));
            return types;
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
}

// Conversion function for DateTimeArray.
IObjectSP SimpleObservationSchedule::fromDateTimeArray(
                                       const IObjectSP& object, 
                                       CClassConstSP    requiredType) {
    // Copy and convert the DateTimeArray into an IDateBuilder. This is a
    // non-const action, so clone the original object and then convert it.
    IObjectSP myClone = IObjectSP(object->clone());
    CObject::checkType(myClone, IDateBuilder::TYPE);
    IDateBuilderSP dates = IDateBuilderSP(dynamic_cast<IDateBuilder*>(myClone.get()));
    // The default behaviour is to use the 'NotUsed' ObservationType.
    ObservationTypeSP notUsed = ObservationType::make("NotUsed");
    return IObjectSP(new SimpleObservationSchedule(dates,
                notUsed, "DEFAULT"));    
}

// <<IRevertTypeConvert>> interface. Revert object to type for interface.
IObjectSP SimpleObservationSchedule::revert(const string& interfaceType) const {
    static const string method = "SimpleObservationSchedule::revert";
    try {
        if (interfaceType != IRevertTypeConvert::PYRAMID) {
            throw ModelException( method, 
                                  "Cannot convert a SimpleObservationSchedule for"
                                  " the interface " + interfaceType);               
        }
        
        // SimpleObservationSchedule may be represented as a DateTimeArray
        IObjectSP dateArray(dates->dates());
        return dateArray;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Explicit constructor
SimpleObservationSchedule::SimpleObservationSchedule(
                IDateBuilderSP theDateBuilder,
                ObservationTypeSP theObsType,
                const string& theSource ) :
                CObject(TYPE),
                dates(theDateBuilder),
                obsType(theObsType),
                source(theSource) { }

// QLib object factory registration
void SimpleObservationSchedule::load(CClassSP& clazz) {
    REGISTER(SimpleObservationSchedule, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IObservationBuilder);
    IMPLEMENTS(IRevertTypeConvert);
    EMPTY_SHELL_METHOD(defaultSimpleObservationSchedule);
    FIELD(dates, "Observation date builder");
    FIELD(obsType, "Type of observation for this schedule");
    FIELD(source, "Name of the source");

    clazz->setPublic(); // make visible to EAS/spreadsheet

    // Register conversion magic to allow a raw DateTimeArray to be
    // transformed into an <<IObservationBuilder>> SimpleObservationSchedule
    // The *interface* type is used to make this type conversion available
    // to all methods requiring an IObservationBuilder.
    registerObjectFromArrayMethod( DateTimeArray::TYPE,
                                   IObservationBuilder::TYPE,
                                   &fromDateTimeArray );
}

IObject* SimpleObservationSchedule::defaultSimpleObservationSchedule() {
    return new SimpleObservationSchedule();
}

SimpleObservationSchedule::SimpleObservationSchedule() :
        CObject(TYPE) {}

SimpleObservationSchedule::SimpleObservationSchedule(CClassConstSP clazz) :
        CObject(clazz) {}

CClassConstSP const SimpleObservationSchedule::TYPE =
    CClass::registerClassLoadMethod("SimpleObservationSchedule",
                                    typeid(SimpleObservationSchedule),
                                    SimpleObservationSchedule::load);

DRLIB_END_NAMESPACE
