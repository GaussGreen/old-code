//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolAbsoluteShift.cpp
//
//   Description : Scenario shift where an absolute shift (e.g. a 10% shift
//                 to 20% vol gives 30% vol) is applied for each benchmark.
//                 Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//                 1Y -> 5Y, 5Y ->  30Y etc).
//
//   Author      : Andrew J Swain
//
//   Date        : 9 May 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolAbsoluteShift.hpp"

DRLIB_BEGIN_NAMESPACE

VolAbsoluteShift::IShift::IShift(){} // empty
VolAbsoluteShift::IShift::~IShift(){} // empty

/** for reflection */
VolAbsoluteShift::VolAbsoluteShift(): MultiExpiryStepShift(TYPE) {}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VolAbsoluteShift::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VolAbsoluteShift::nameMatches(const OutputName&     name,
                                   IObjectConstSP        obj){
    // cast obj to VolAbsoluteShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(volShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VolAbsoluteShift::appendName(OutputNameArray&       namesList,
                                  IObjectConstSP         obj){
    // cast obj to VolAbsoluteShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VolAbsoluteShift::shift(IObjectSP obj) {
    // cast obj to VolAbsoluteShift::Shift and then invoke shift method
    IShift& volShiftObj = dynamic_cast<IShift&>(*obj);
    return volShiftObj.sensShift(this);
}

    
void VolAbsoluteShift::validatePop2Object() {
    try {
        MultiExpiryShift::validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e, "VolAbsoluteShift::validatePop2Object");
    }
}

class VolAbsoluteShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute shifts to vol by benchmark");
        REGISTER(VolAbsoluteShift, clazz);
        SUPERCLASS(MultiExpiryStepShift);
        EMPTY_SHELL_METHOD(defaultVolAbsoluteShift);
    }

    static IObject* defaultVolAbsoluteShift(){
        return new VolAbsoluteShift();
    }
};

CClassConstSP const VolAbsoluteShift::TYPE = CClass::registerClassLoadMethod(
    "VolAbsoluteShift", typeid(VolAbsoluteShift), VolAbsoluteShiftHelper::load);

CClassConstSP const VolAbsoluteShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VolAbsoluteShift::IShift", typeid(VolAbsoluteShift::IShift), 0);


DRLIB_END_NAMESPACE
