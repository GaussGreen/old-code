//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BarrierLevel.cpp
//
//   Description : Results data holder for barrier level today
//
//   Author      : Andrew J Swain
//
//   Date        : 10 June 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_BARRIERLEVEL_CPP
#include "edginc/BarrierLevel.hpp"

DRLIB_BEGIN_NAMESPACE

// how many days forward we report barriers
const int BarrierLevel::DATE_WINDOW = 7;

BarrierLevel::BarrierLevel(bool isUp, const DateTime& date, double level,
                           bool isContinuous): 
    CObject(TYPE), isUp(isUp), date(date), 
    level(level), isContinuous(isContinuous) {}

// given today's date, when should we report barriers up to ?
DateTime BarrierLevel::barrierWindow(const DateTime& today) {
    // we're too low level to have access to holidays at this point
    return today.rollDate(DATE_WINDOW);
}

BarrierLevel::BarrierLevel() : CObject(TYPE), isUp(false), 
    level(0.0), isContinuous(true){}

//// Routes through equalTo. Method added to support
//// instantiating array template
bool BarrierLevel::operator==(const BarrierLevel& rhs) const{
    return equalTo(&rhs);
}

class BarrierLevelHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BarrierLevel, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBarrierLevel);
        FIELD(isUp, "is an up barrier");
        FIELD(date, "date");
        FIELD(level, "level");
        FIELD(isContinuous, "is it a continuous barrier");
        FIELD_MAKE_OPTIONAL(isContinuous); // for now
    }
    static IObject* defaultBarrierLevel(){
        return new BarrierLevel();
    }
};

CClassConstSP const BarrierLevel::TYPE = CClass::registerClassLoadMethod(
    "BarrierLevel", typeid(BarrierLevel), BarrierLevelHelper::load);

DEFINE_TEMPLATE_TYPE(BarrierLevelArray);

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<BarrierLevel>::toIObject(
    const BarrierLevel& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<BarrierLevel>::toIObject(BarrierLevel& value){
    return BarrierLevelSP::attachToRef(&value);
}

/** Turns the IObjectSP into a BarrierLevel */
BarrierLevel arrayObjectCast<BarrierLevel>::fromIObject(IObjectSP& value){
    BarrierLevel *dtPtr = dynamic_cast<BarrierLevel *>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a BarrierLevel");
    }
    return *dtPtr;
}

// explicit clone for arrays of BarrierLevels - for performance
IObject* arrayClone<BarrierLevel>::clone(const CArray* arrayToClone){
    const BarrierLevelArray& theArray = 
        static_cast<const BarrierLevelArray&>(*arrayToClone);
    return new BarrierLevelArray(theArray);
}

DRLIB_END_NAMESPACE

