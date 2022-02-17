//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FwdLevel.cpp
//
//   Description : spot level scenario - set spot to supplied value
//                 at a future date
//
//   Author      : Andrew J Swain
//
//   Date        : 9 June 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FwdLevel.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor with explicit spot level & date */
FwdLevel::FwdLevel(double spot, const DateTime& fwdDate):
    SpotLevel(TYPE, spot), fwdDate(fwdDate){}

/** for reflection */
FwdLevel::FwdLevel(): 
    SpotLevel(TYPE, 0.0){
}

// given today's date, when is spot defined
DateTime FwdLevel::spotDate(const DateTime& today) const {
    return fwdDate;
}

class FwdLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FwdLevel, clazz);
        SUPERCLASS(SpotLevel);
        EMPTY_SHELL_METHOD(defaultFwdLevel);
        FIELD(fwdDate, "when fwd is defined for");
    }

    static IObject* defaultFwdLevel(){
        return new FwdLevel();
    }
};

CClassConstSP const FwdLevel::TYPE = CClass::registerClassLoadMethod(
    "FwdLevel", typeid(FwdLevel), FwdLevelHelper::load);

DRLIB_END_NAMESPACE
