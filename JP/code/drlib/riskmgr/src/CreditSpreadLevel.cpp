//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CreditSpreadsLevel.cpp
//
//   Description : credit spread level scenario - set CS to supplied value
//
//   Author      : Tycho von Rosenvinge and Andre Segger
//
//   Date        : 23 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditSpreadLevel.hpp"

DRLIB_BEGIN_NAMESPACE

CreditSpreadLevel::Shift::~Shift(){} // empty

/** constructor with explicit spread level */
CreditSpreadLevel::CreditSpreadLevel(double spread):
    ScalarPerturbation(TYPE, spread){}

/** for reflection */
CreditSpreadLevel::CreditSpreadLevel(): ScalarPerturbation(TYPE){ }

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditSpreadLevel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool CreditSpreadLevel::nameMatches(const OutputName&         name,
                                    IObjectConstSP          obj){
    // cast obj to CreditSpreadLevel::Shift and then invoke name method
    const Shift& creditSpreadLevelObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(creditSpreadLevelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void CreditSpreadLevel::appendName(OutputNameArray&          namesList,
                          IObjectConstSP          obj){
    // cast obj to CreditSpreadLevel::Shift and then invoke name method
    const Shift& creditSpreadLevelObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(
        new OutputName(creditSpreadLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditSpreadLevel::shift(IObjectSP obj) {
    // cast obj to CreditSpreadLevel::Shift and then invoke shift method
    Shift& creditSpreadLevelObj = 
        dynamic_cast<Shift&>(*obj);
    return creditSpreadLevelObj.sensShift(this);
}

class CreditSpreadLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditSpreadLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultCreditSpreadLevel);
    }

    static IObject* defaultCreditSpreadLevel(){
        return new CreditSpreadLevel();
    }
};

CClassConstSP const CreditSpreadLevel::TYPE = CClass::registerClassLoadMethod(
    "CreditSpreadLevel", typeid(CreditSpreadLevel), CreditSpreadLevelHelper::load);

CClassConstSP const CreditSpreadLevel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CreditSpreadLevel::Shift", typeid(CreditSpreadLevel::Shift), 0);

DRLIB_END_NAMESPACE
