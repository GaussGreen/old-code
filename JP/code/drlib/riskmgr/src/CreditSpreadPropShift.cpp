//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadPropShift.cpp
//
//   Description : Credit Spread Proportional Shift scenario - apply proportional shift to CS curve
//
//   Author      : Andrew McCleery
//
//   Date        : 24 February 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditSpreadPropShift.hpp"

DRLIB_BEGIN_NAMESPACE

CreditSpreadPropShift::IShift::~IShift(){} // empty

/** constructor with explicit shift */
CreditSpreadPropShift::CreditSpreadPropShift(double shift):
    ScalarPerturbation(TYPE, shift){
}

/** for reflection */
CreditSpreadPropShift::CreditSpreadPropShift(): ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditSpreadPropShift::shiftInterface() const{
    return IShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    IShift interface */
bool CreditSpreadPropShift::nameMatches(const OutputName&         name,
                                        IObjectConstSP            obj){
    // cast obj to CreditSpreadPropShift::IShift and then invoke name method
    const IShift& creditSpreadPropShiftObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(creditSpreadPropShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's IShift interface */
void CreditSpreadPropShift::appendName(OutputNameArray&          namesList,
                                       IObjectConstSP            obj){
    // cast obj to CreditSpreadPropShift::IShift and then invoke name method
    const IShift& creditSpreadPropShiftObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(
        new OutputName(creditSpreadPropShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditSpreadPropShift::shift(IObjectSP obj) {
    // cast obj to CreditSpreadPropShift::IShift and then invoke shift method
    IShift& creditSpreadPropShiftObj =
        dynamic_cast<IShift&>(*obj);
    return creditSpreadPropShiftObj.sensShift(this);
}


class CreditSpreadPropShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditSpreadPropShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultCreditSpreadPropShift);
        // no fields
    }

    static IObject* defaultCreditSpreadPropShift(){
        return new CreditSpreadPropShift();
    }
};

CClassConstSP const CreditSpreadPropShift::TYPE = CClass::registerClassLoadMethod(
    "CreditSpreadPropShift", typeid(CreditSpreadPropShift), CreditSpreadPropShiftHelper::load);

CClassConstSP const CreditSpreadPropShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CreditSpreadPropShift::IShift", typeid(CreditSpreadPropShift::IShift), 0);

DRLIB_END_NAMESPACE
