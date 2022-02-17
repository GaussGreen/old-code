#include "edginc/config.hpp"
#include "edginc/OptionSched.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/************************ OptionSchedInit ************************/

OptionSchedInit::OptionSchedInit() : CObject(TYPE), stubPos(StubPos::STUB_NONE) {}

void OptionSchedInit::validatePop2Object(void) {
    if (startDate>endDate) 
        throw ModelException("startDate>endDate");
    if (!frequency.get()) 
        throw ModelException("frequency needed");
}


void OptionSchedInit::modify(OptionSchedDates &sched) {
    try {
        DateTimeArray dates;

        RatesUtils::genNewEvenDateList(startDate, endDate, 
            frequency->approxAnnualFrequency(), stubPos==StubPos::STUB_END, dates);

        DateTimeArray &notifDates = getDates(sched.notifDate);
        notifDates.clear();
        notifDates.insert(notifDates.end(), dates.begin(), dates.end()); 
        DateTimeArray &exerciseDates = getDates(sched.exerciseDate);
        exerciseDates.clear();
        exerciseDates.insert(exerciseDates.end(), dates.begin(), dates.end()); 
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void OptionSchedInit::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(OptionSchedInit, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(startDate, "");
    FIELD(endDate, "");
    FIELD(frequency, "");
    FIELD(stubPos, "");   FIELD_MAKE_OPTIONAL(stubPos);

    Addin::registerConstructor(Addin::UTILITIES, OptionSchedInit::TYPE);
}

CClassConstSP const OptionSchedInit::TYPE = CClass::registerClassLoadMethod(
    "OptionSchedInit", typeid(OptionSchedInit), load);

/************************ OptionSchedDates ************************/

static void check(int a, int b, const string& s, bool zeroOk) {
    if (a==b) return;
    if (a==0 && zeroOk) return;
    throw ModelException(s+".size="+Format::toString(a)+" should be "+Format::toString(b));
}

void OptionSchedDates::validatePop2Object() {
    try {
        if (initializer.get()) {
            initializer->modify(*this);
            // re-apply the modifiers to adjust the dates
            notifDate.modifyDates();
            exerciseDate.modifyDates();
        }

        nbDates = notifDate.size();
        if (!nbDates) 
            throw ModelException("notifDate date list empty");
        check(exerciseDate.size(), nbDates, "exerciseDate", false);

        DateTime::ensureIncreasing(notifDate.getDates(), "notif", true);
        DateTime::ensureIncreasing(exerciseDate.getDates(), "exercise", true);

        // notification dates must be <= exercise date
        int i;
        for (i = 0; i < nbDates; i++){
            if (notifDate[i] > exerciseDate[i]) {
                throw ModelException("notifDate[" 
                    + Format::toString(i) + "] ("
                    + notifDate[i].toString() + ") must be <= exerciseDate["
                    + exerciseDate[i].toString()+"]");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void OptionSchedDates::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(OptionSchedDates, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(initializer, ""); FIELD_MAKE_OPTIONAL(initializer);
    FIELD(notifDate, ""); FIELD_MAKE_OPTIONAL(notifDate);
    FIELD(exerciseDate, "");   FIELD_MAKE_OPTIONAL(exerciseDate);

    FIELD(nbDates, ""); FIELD_MAKE_TRANSIENT(nbDates);

    Addin::registerConstructor(Addin::UTILITIES, OptionSchedDates::TYPE);
}

CClassConstSP const OptionSchedDates::TYPE = CClass::registerClassLoadMethod(
    "OptionSchedDates", typeid(OptionSchedDates), load);

/******************************** OptionSchedDisp ********************************/

class OptionSchedDisp : public CObject {
public:
    static CClassConstSP const TYPE;

    struct Dates {
        enum Enum {NOTIF, EXERCISE};
    };
    typedef BoxedEnum<Dates::Enum> DatesBoxedEnum;

    IObjectSP addin() {
        DateTimeArray const *a;
        switch (datesToDisplay) {
            case Dates::NOTIF: a = &val->notifDate.getDates(); break;
            case Dates::EXERCISE:   a = &val->exerciseDate.getDates(); break;
            default: ASSERT(0);
        }
        return IObjectSP(new DateTimeArray(*a));
    }

private:
    OptionSchedDatesSP val;
    Dates::Enum datesToDisplay;

    OptionSchedDisp() : CObject(TYPE)          {}
    static IObject* defaultConstructor(void) {return new OptionSchedDisp();}
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(OptionSchedDisp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);

        FIELD(val, "The OptionSched to display");    
        FIELD(datesToDisplay, "");    
        Addin::registerObjectMethod(
            "OptionSchedDisp", Addin::UTILITIES,
            "Displays the dates of a OptionSched object",
            false, Addin::expandSimple, &OptionSchedDisp::addin);
    }
};

CClassConstSP const OptionSchedDisp::TYPE = CClass::registerClassLoadMethod(
    "OptionSchedDisp", typeid(OptionSchedDisp), load);


START_PUBLIC_ENUM_DEFINITION(OptionSchedDisp::Dates::Enum, "");
ENUM_VALUE_AND_NAME(OptionSchedDisp::Dates::NOTIF, "NOTIF", "");
ENUM_VALUE_AND_NAME(OptionSchedDisp::Dates::EXERCISE, "EXERCISE", "");
END_ENUM_DEFINITION(OptionSchedDisp::Dates::Enum);


bool OptionSchedLoad() {return OptionSchedDates::TYPE !=0;}

DRLIB_END_NAMESPACE
