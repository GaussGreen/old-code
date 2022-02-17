//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotShift.hpp
//
//   Description : spot level scenario - set spot to supplied value
//
//   Author      : Mark A Robson
//
//   Date        : 16 Jul 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SpotShift.hpp"

DRLIB_BEGIN_NAMESPACE

SpotShift::Shift::~Shift(){} // empty

/** constructor with explicit % spot shift amount */
SpotShift::SpotShift(double spot): ScalarPerturbation(TYPE, spot){}

/** for reflection */
SpotShift::SpotShift(): ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP SpotShift::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool SpotShift::nameMatches(const OutputName&         name,
                            IObjectConstSP          obj){
    // cast obj to SpotShift::Shift and then invoke name method
    const Shift& SpotShiftObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(SpotShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void SpotShift::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to SpotShift::Shift and then invoke name method
    const Shift& SpotShiftObj = dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(SpotShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool SpotShift::shift(IObjectSP obj) {
    // cast obj to SpotShift::Shift and then invoke shift method
    Shift& SpotShiftObj = dynamic_cast<Shift&>(*obj);
    return SpotShiftObj.sensShift(this);
}

class SpotShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SpotShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultSpotShift);
    }

    static IObject* defaultSpotShift(){
        return new SpotShift();
    }
};

CClassConstSP const SpotShift::TYPE = CClass::registerClassLoadMethod(
    "SpotShift", typeid(SpotShift), SpotShiftHelper::load);

CClassConstSP const SpotShift::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "SpotShift::Shift", typeid(SpotShift::Shift), 0);


DRLIB_END_NAMESPACE
