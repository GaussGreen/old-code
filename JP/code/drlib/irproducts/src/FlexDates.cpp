#include "edginc/config.hpp"
#include "edginc/FlexDates.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE


/************************ DatesModifier ************************/
void IDatesModifier::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IDatesModifier, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IDatesModifier::TYPE = CClass::registerInterfaceLoadMethod(
    "IDatesModifier", typeid(IDatesModifier), load);

/******************************** FlexDates ********************************/


void FlexDates::filterOnWeight(DateTimeArray &datesIn, DoubleArray &weights, DateTimeArray &datesOut) {
    try {
        if (datesIn.size()!=weights.size()) throw ModelException("datesIn.size()!=weights.size()");
        datesOut.clear();
        datesOut.reserve(datesIn.size());
        for (int i=0; i<datesIn.size(); ++i) {
            if (!Maths::isZero(weights[i])) {
                datesOut.push_back(datesIn[i]);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "FlexDates::filterOnWeight");
    }
}

void FlexDates::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FlexDates, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD       (datesExpl, ""); FIELD_MAKE_OPTIONAL(datesExpl); 
    FIELD       (modifier, "");  FIELD_MAKE_OPTIONAL(modifier); 
    FIELD(dates, "");     FIELD_MAKE_TRANSIENT(dates); 

    Addin::registerObjectMethod(
        "DatesModifier_view", Addin::UTILITIES,
        "Returns a DateTimeArray from a FlexDates",
        true, Addin::expandSimple, &FlexDates::addinDisplay);

    Addin::registerConstructor(Addin::UTILITIES, FlexDates::TYPE);
}

CClassConstSP const FlexDates::TYPE = CClass::registerClassLoadMethod(
    "FlexDates", typeid(FlexDates), load);

/******************************** DatesSched ********************************/

void DatesSched::modify(DateTimeArray &dates) {
    try {
        if (startDate>matDate) throw ModelException("startDate>matDate");
        if (!frequency.get()) throw ModelException("frequency needed");
        RatesUtils::genNewEvenDateList(startDate, matDate, 
            frequency->approxAnnualFrequency(), stubPos==StubPos::STUB_END, dates);
    }
    catch (exception& e) {
        throw ModelException(e, "DatesSched::modify");
    }
}

void DatesSched::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DatesSched, clazz);
    SUPERCLASS(CObject)
    IMPLEMENTS(IDatesModifier);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(startDate, "");    
    FIELD(matDate, "");      
    FIELD       (frequency, "");    
    FIELD(stubPos,"");     

    Addin::registerConstructor(Addin::UTILITIES, DatesSched::TYPE);
}

CClassConstSP const DatesSched::TYPE = CClass::registerClassLoadMethod(
    "DatesSched", typeid(DatesSched), load);

/******************************** DatesAdj ********************************/

void DatesAdj::modify(DateTimeArray &dates) {
    try {
        if (modifier.get()) 
            modifier->modify(dates);

        if (dates.size()==0) 
            return;

        switch (rangeBound) {
            case RangeBound::START: dates.pop_back(); break;
            case RangeBound::END:   { 
                DateTimeArray tmp = dates; 
                dates.clear(); 
                dates.insert(dates.end(), tmp.begin()+1, tmp.end()); 
                break; 
            }
            case RangeBound::NONE: break;
            default: ASSERT(0);
        }

        if (shiftMode!=ShiftMode::NONE) {
            if (!shiftLength.get())
                throw ModelException("shift requested but shiftLength not provided");

            for (int i=0; i<dates.size(); ++i) {
                dates[i] = RatesUtils::fwdByExpiry(dates[i], *shiftLength, shiftMode == ShiftMode::FWD);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "DatesAdj::modify");
    }
}

void DatesAdj::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DatesAdj, clazz);
    SUPERCLASS(CObject)
    IMPLEMENTS(IDatesModifier);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD       (modifier, "");   FIELD_MAKE_OPTIONAL(modifier);   
    FIELD(rangeBound, ""); FIELD_MAKE_OPTIONAL(rangeBound);
    FIELD       (shiftLength,""); FIELD_MAKE_OPTIONAL(shiftLength);
    FIELD(shiftMode,"");   FIELD_MAKE_OPTIONAL(shiftMode);

    Addin::registerConstructor(Addin::UTILITIES, DatesAdj::TYPE);
}

CClassConstSP const DatesAdj::TYPE = CClass::registerClassLoadMethod(
    "DatesAdj", typeid(DatesAdj), load);


START_PUBLIC_ENUM_DEFINITION(DatesAdj::RangeBound::Enum, "");
ENUM_VALUE_AND_NAME(DatesAdj::RangeBound::START, "START", "");
ENUM_VALUE_AND_NAME(DatesAdj::RangeBound::END, "END", "");
ENUM_VALUE_AND_NAME(DatesAdj::RangeBound::NONE, "NONE", "");
END_ENUM_DEFINITION(DatesAdj::RangeBound::Enum);

START_PUBLIC_ENUM_DEFINITION(DatesAdj::ShiftMode::Enum, "");
ENUM_VALUE_AND_NAME(DatesAdj::ShiftMode::NONE, "NONE", "");
ENUM_VALUE_AND_NAME(DatesAdj::ShiftMode::FWD, "FWD", "");
ENUM_VALUE_AND_NAME(DatesAdj::ShiftMode::BWD, "BWD", "");
END_ENUM_DEFINITION(DatesAdj::ShiftMode::Enum);

/******************************** FlexDatesDisp ********************************/

class FlexDatesDisp : public CObject,
                      virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;
    virtual IObjectSP run() {return IObjectSP(new DateTimeArray(val->getDates()));}
private:
    FlexDatesSP val;
    CMarketDataSP market;
    FlexDatesDisp() : CObject(TYPE)          {}
    static IObject* defaultConstructor(void) {return new FlexDatesDisp();}
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FlexDatesDisp, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);

        FIELD(val, "The FlexDates to display");    
        FIELD(market, "Market to use to compute the dates (if holidays needed). NOT IMPLEMENTED YET.");    
        FIELD_MAKE_OPTIONAL(market);
        Addin::registerObjectMethod(
            "FlexDatesDisp", Addin::UTILITIES,
            "Displays the dates of a FlexDates object",
            false, Addin::expandSimple, & FlexDatesDisp::run);
    }
};

CClassConstSP const FlexDatesDisp::TYPE = CClass::registerClassLoadMethod(
    "FlexDatesDisp", typeid(FlexDatesDisp), load);


bool FlexDatesLoad() {return FlexDates::TYPE !=0;}

DRLIB_END_NAMESPACE
