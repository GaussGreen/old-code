//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolParallelShift.cpp
//
//   Description : vol shift scenario - add parallel shift to vol
//
//   Author      : Andrew J Swain
//
//   Date        : 5 June 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/VegaProxyParallel.hpp"

DRLIB_BEGIN_NAMESPACE
VolParallelShift::Shift::~Shift(){} // empty
VolParallelShift::Shift::Shift(){} // empty

const double VolParallelShift::MIN_SPOT_VOL = 0.05;
const double VolParallelShift::MIN_FWD_VOL  = 0.05;

/** constructor with explicit shift */
VolParallelShift::VolParallelShift(double shift):
    ScalarPerturbation(TYPE, shift), shiftFundProxies(true){}

/** for reflection */
VolParallelShift::VolParallelShift(): 
    ScalarPerturbation(TYPE), shiftFundProxies(true) {}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VolParallelShift::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VolParallelShift::nameMatches(const OutputName&         name,
                                   IObjectConstSP          obj){
    // cast obj to VolParallelShift::Shift and then invoke name method
    const Shift& volShiftObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(volShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VolParallelShift::appendName(OutputNameArray&          namesList,
                                  IObjectConstSP          obj){
    // cast obj to VolParallelShift::Shift and then invoke name method
    const Shift& volShiftObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VolParallelShift::shift(IObjectSP obj) {
    // cast obj to VolParallelShift::Shift and then invoke shift method
    Shift& volShiftObj =
        dynamic_cast<Shift&>(*obj);
    return volShiftObj.sensShift(this);
}

// allows SensControl to control its own shifting - default is via SensMgr
bool VolParallelShift::findAndShift(IObjectSP         objectToShift, 
                                    OutputNameConstSP name) {
    try {
        bool shifted = 
            ScalarPerturbation::findAndShift(objectToShift, name);// call parent
        if (shiftFundProxies) {
            // build up a proxy parallel shift and do that too
            shifted = VegaProxyParallel(shiftSize).
                findAndShift(objectToShift, name) || shifted;
        }
        return shifted;
    } catch (exception& e) {
        throw ModelException(e, "VolParallelShift::findAndShift");
    }
}

class VolParallelShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolParallelShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultVolParallelShift);
        FIELD(shiftFundProxies, "shift fund proxy vols");
        FIELD_MAKE_OPTIONAL(shiftFundProxies)
    }

    static IObject* defaultVolParallelShift(){
        return new VolParallelShift();
    }
};

CClassConstSP const VolParallelShift::TYPE = CClass::registerClassLoadMethod(
    "VolParallelShift", typeid(VolParallelShift), VolParallelShiftHelper::load);

CClassConstSP const VolParallelShift::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VolParallelShift::Shift", typeid(VolParallelShift::Shift), 0);

DRLIB_END_NAMESPACE
