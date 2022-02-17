//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadPropShift.cpp
//
//   Description : CDS Par Spread Proportional Shift scenario 
//                 - apply proportional shift to CDS spread curve
//
//   Author      : Andrew McCleery
//
//   Date        : 24 February 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ParSpreadPropShift.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadPropShift::IShift::~IShift(){} // empty

/** constructor with explicit shift */
ParSpreadPropShift::ParSpreadPropShift(double shift):
    ScalarPerturbation(TYPE, shift){}

/** for reflection */
ParSpreadPropShift::ParSpreadPropShift():
    ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ParSpreadPropShift::shiftInterface() const{
    return IShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    IShift interface */
bool ParSpreadPropShift::nameMatches(const OutputName&         name,
                                     IObjectConstSP            obj){
    // cast obj to ParSpreadPropShift::IShift and then invoke name method
    const IShift& parSpreadPropShiftObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(parSpreadPropShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's IShift interface */
void ParSpreadPropShift::appendName(OutputNameArray&          namesList,
                                    IObjectConstSP            obj){
    // cast obj to ParSpreadPropShift::IShift and then invoke name method
    const IShift& parSpreadPropShiftObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(
        new OutputName(parSpreadPropShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ParSpreadPropShift::shift(IObjectSP obj) {
    // cast obj to ParSpreadPropShift::IShift and then invoke shift method
    IShift& parSpreadPropShiftObj =
        dynamic_cast<IShift&>(*obj);
    return parSpreadPropShiftObj.sensShift(this);
}

class ParSpreadPropShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ParSpreadPropShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultParSpreadPropShift);
        // no fields
    }

    static IObject* defaultParSpreadPropShift(){
        return new ParSpreadPropShift();
    }
};

CClassConstSP const ParSpreadPropShift::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadPropShift", typeid(ParSpreadPropShift), ParSpreadPropShiftHelper::load);

CClassConstSP const ParSpreadPropShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ParSpreadPropShift::IShift", typeid(ParSpreadPropShift::IShift), 0);

DRLIB_END_NAMESPACE
