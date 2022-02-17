//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadLevel.cpp
//
//   Description : CDS par spread level scenario - set spreads to supplied value
//
//   Author      : Andrew McCleery
//
//   Date        : 16 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ParSpreadLevel.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadLevel::IShift::~IShift(){} // empty

/** constructor with explicit spread level */
ParSpreadLevel::ParSpreadLevel(double spread): 
    ScalarPerturbation(TYPE, spread){}

/** for reflection */
ParSpreadLevel::ParSpreadLevel(): ScalarPerturbation(TYPE){
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ParSpreadLevel::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    IShift interface */
bool ParSpreadLevel::nameMatches(const OutputName&         name,
                                 IObjectConstSP            obj){
    // cast obj to ParSpreadLevel::IShift and then invoke name method
    const IShift& parSpreadLevelObj = 
        dynamic_cast<const IShift&>(*obj);
    return name.equals(parSpreadLevelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's IShift interface */
void ParSpreadLevel::appendName(OutputNameArray&          namesList,
                                IObjectConstSP            obj){
    // cast obj to ParSpreadLevel::IShift and then invoke name method
    const IShift& parSpreadLevelObj = 
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(
        new OutputName(parSpreadLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ParSpreadLevel::shift(IObjectSP obj) {
    // cast obj to ParSpreadLevel::IShift and then invoke shift method
    IShift& parSpreadLevelObj = 
        dynamic_cast<IShift&>(*obj);
    return parSpreadLevelObj.sensShift(this);
}

class ParSpreadLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ParSpreadLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultParSpreadLevel);
        // no fields
    }

    static IObject* defaultParSpreadLevel(){
        return new ParSpreadLevel();
    }
};

CClassConstSP const ParSpreadLevel::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadLevel", typeid(ParSpreadLevel), ParSpreadLevelHelper::load);

CClassConstSP const ParSpreadLevel::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ParSpreadLevel::IShift", typeid(ParSpreadLevel::IShift), 0);

DRLIB_END_NAMESPACE
