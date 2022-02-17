//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadParallelShift.cpp
//
//   Description : CDS Par Spread Parallel Shift scenario
//                 - apply parallel shift to CDS par spread curve
//
//   Author      : Andrew McCleery
//
//   Date        : 25 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ParSpreadParallelShift.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadParallelShift::IShift::~IShift(){} // empty

// constructor with explicit shift size 
ParSpreadParallelShift::ParSpreadParallelShift(double shiftSize):
    ScalarPerturbation(TYPE, shiftSize){}

// for reflection 
ParSpreadParallelShift::ParSpreadParallelShift(): 
    ScalarPerturbation(TYPE){}

// returns the interface identifying what an object has to do in order
// to be support the tweak that this object represents 
CClassConstSP ParSpreadParallelShift::shiftInterface() const{
    return IShift::TYPE;
}
 
// Returns true if the supplied object matches the supplied name
// for this sensitivity. The object must implement the
// IShift interface 
bool ParSpreadParallelShift::nameMatches(const OutputName&         name,
                                         IObjectConstSP            obj) {
    // cast obj to ParSpreadParallelShift::IShift and then invoke name method
    const IShift& PSParallelShiftObj = 
        dynamic_cast<const IShift&>(*obj);
    return name.equals(PSParallelShiftObj.sensName(this));
}

// Appends the name(s) of the supplied object with respect to this
// sensitivity to the supplied list. The object must implement the
// IShift interface 
void ParSpreadParallelShift::appendName(OutputNameArray&          namesList,
                                        IObjectConstSP            obj){
    // cast obj to ParSpreadParallelShift::IShift and then invoke name method
    const IShift& PSParallelShiftObj = 
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(PSParallelShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ParSpreadParallelShift::shift(IObjectSP obj) {
    // cast obj to ParSpreadParallelShift::IShift and then invoke shift method
    IShift& PSParallelShiftObj = 
        dynamic_cast<IShift&>(*obj);
    return PSParallelShiftObj.sensShift(this);
}

class ParSpreadParallelShiftHelper {
public:
    // Invoked when Class is 'loaded' 
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ParSpreadParallelShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultPSParallelShift);
        // no fields
    }

    static IObject* defaultPSParallelShift(){
        return new ParSpreadParallelShift();
    }
};

CClassConstSP const ParSpreadParallelShift::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadParallelShift", typeid(ParSpreadParallelShift), ParSpreadParallelShiftHelper::load);

CClassConstSP const ParSpreadParallelShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ParSpreadParallelShift::IShift", typeid(ParSpreadParallelShift::IShift), 0);
DRLIB_END_NAMESPACE
