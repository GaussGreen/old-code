//----------------------------------------------------------------------------
//
//   Group       : EDG Quantitative Research
//
//   Filename    : ContangoTermShift.cpp
//
//   Description : Scenario shift where an absolute shift (e.g. a 10% shift
//                 to 20% rate gives 30% rate) is applied for each benchmark.
//                 Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//                 1Y -> 5Y, 5Y ->  30Y etc).
//
//   Author      : Andrew J Swain
//
//   Date        : 19 April 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ContangoTermShift.hpp"

DRLIB_BEGIN_NAMESPACE

ContangoTermShift::IShift::IShift(){}  // empty
ContangoTermShift::IShift::~IShift(){} // empty

/** for reflection */
ContangoTermShift::ContangoTermShift(): MultiExpiryStepShift(TYPE) {}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ContangoTermShift::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool ContangoTermShift::nameMatches(const OutputName&     name,
                                   IObjectConstSP        obj){
    // cast obj to ContangoTermShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(volShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void ContangoTermShift::appendName(OutputNameArray&       namesList,
                                  IObjectConstSP         obj){
    // cast obj to ContangoTermShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool ContangoTermShift::shift(IObjectSP obj) {
    // cast obj to ContangoTermShift::Shift and then invoke shift method
    IShift& volShiftObj = dynamic_cast<IShift&>(*obj);
    return volShiftObj.sensShift(this);
}

    
void ContangoTermShift::validatePop2Object() {
    try {
        MultiExpiryShift::validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e, "ContangoTermShift::validatePop2Object");
    }
}

class ContangoTermShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute shifts to contango curve by benchmark");
        REGISTER(ContangoTermShift, clazz);
        SUPERCLASS(MultiExpiryStepShift);
        EMPTY_SHELL_METHOD(defaultContangoTermShift);
    }

    static IObject* defaultContangoTermShift(){
        return new ContangoTermShift();
    }
};

CClassConstSP const ContangoTermShift::TYPE = CClass::registerClassLoadMethod(
    "ContangoTermShift", typeid(ContangoTermShift), ContangoTermShiftHelper::load);

CClassConstSP const ContangoTermShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ContangoTermShift::IShift", typeid(ContangoTermShift::IShift), 0);


DRLIB_END_NAMESPACE
