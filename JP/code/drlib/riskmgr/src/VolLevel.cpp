//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolLevel.hpp
//
//   Description : vol level scenario - set vol to supplied value
//
//   Author      : Andrew J Swain
//
//   Date        : 8 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolLevel.hpp"

DRLIB_BEGIN_NAMESPACE
VolLevel::Shift::Shift(){} // empty
VolLevel::Shift::~Shift(){} // empty

/** constructor with explicit vol level */
VolLevel::VolLevel(double vol): ScalarPerturbation(TYPE, vol){}

/** for reflection */
VolLevel::VolLevel(): ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VolLevel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VolLevel::nameMatches(const OutputName&         name,
                           IObjectConstSP          obj){
    // cast obj to VolLevel::Shift and then invoke name method
    const Shift& volLevelObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(volLevelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VolLevel::appendName(OutputNameArray&          namesList,
                          IObjectConstSP          obj){
    // cast obj to VolLevel::Shift and then invoke name method
    const Shift& volLevelObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(volLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VolLevel::shift(IObjectSP obj) {
    // cast obj to VolLevel::Shift and then invoke shift method
    Shift& volLevelObj =
        dynamic_cast<Shift&>(*obj);
    return volLevelObj.sensShift(this);
}

class VolLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultVolLevel);
    }

    static IObject* defaultVolLevel(){
        return new VolLevel();
    }
};

CClassConstSP const VolLevel::TYPE = CClass::registerClassLoadMethod(
    "VolLevel", typeid(VolLevel), VolLevelHelper::load);

CClassConstSP const VolLevel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VolLevel::Shift", typeid(VolLevel::Shift), 0);

DRLIB_END_NAMESPACE
