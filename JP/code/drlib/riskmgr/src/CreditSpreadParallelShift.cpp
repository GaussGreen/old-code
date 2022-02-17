//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadParallelShift.cpp
//
//   Description : Credit Spread Parallel Shift scenario
//    	           - apply parallel shift to credit spread curve
//
//   Author      : Andrew McCleery
//
//   Date        : 25 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditSpreadParallelShift.hpp"

DRLIB_BEGIN_NAMESPACE

CreditSpreadParallelShift::IShift::~IShift(){} // empty

// constructor with explicit shift size
CreditSpreadParallelShift::CreditSpreadParallelShift(double shiftSize):
    ScalarPerturbation(TYPE, shiftSize){}

// for reflection
CreditSpreadParallelShift::CreditSpreadParallelShift(): 
    ScalarPerturbation(TYPE){}

// returns the interface identifying what an object has to do in order
// to be support the tweak that this object represents
CClassConstSP CreditSpreadParallelShift::shiftInterface() const{
    return IShift::TYPE;
}

// Returns true if the supplied object matches the supplied name
// for this sensitivity. The object must implement the
// IShift interface
bool CreditSpreadParallelShift::nameMatches(const OutputName&         name,
                                            IObjectConstSP            obj) {
    // cast obj to CreditSpreadParallelShift::IShift and then invoke name method
    const IShift& CSParallelShiftObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(CSParallelShiftObj.sensName(this));
}

// Appends the name(s) of the supplied object with respect to this
// sensitivity to the supplied list. The object must implement the
// IShift interface
void CreditSpreadParallelShift::appendName(OutputNameArray&          namesList,
                                           IObjectConstSP            obj){
    // cast obj to CreditSpreadParallelShift::IShift and then invoke name method
    const IShift& CSParallelShiftObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(CSParallelShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditSpreadParallelShift::shift(IObjectSP obj) {
    // cast obj to CreditSpreadParallelShift::IShift and then invoke shift method
    IShift& CSParallelShiftObj =
        dynamic_cast<IShift&>(*obj);
    return CSParallelShiftObj.sensShift(this);
}


class CreditSpreadParallelShiftHelper {
public:
    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditSpreadParallelShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultCSParallelShift);
        // no fields
    }

    static IObject* defaultCSParallelShift(){
        return new CreditSpreadParallelShift();
    }
};

CClassConstSP const CreditSpreadParallelShift::TYPE = CClass::registerClassLoadMethod(
    "CreditSpreadParallelShift", typeid(CreditSpreadParallelShift), CreditSpreadParallelShiftHelper::load);

CClassConstSP const CreditSpreadParallelShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CreditSpreadParallelShift::IShift", typeid(CreditSpreadParallelShift::IShift), 0);

DRLIB_END_NAMESPACE
