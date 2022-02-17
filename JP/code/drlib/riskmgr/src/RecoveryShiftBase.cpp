//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoveryShiftBase.cpp
//
//   Description : Sensitivity to Recovery
//
//   Author      : Stephen Hope
//
//   Date        : 19 August 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RecoveryShiftBase.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
RecoveryShiftBase::IShift::~IShift(){} // empty
RecoveryShiftBase::IRestorableShift::~IRestorableShift(){} // empty

const string RecoveryShiftBase::NAME = "RECOVERY_TWEAK";
const double RecoveryShiftBase::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
RecoveryShiftBase::RecoveryShiftBase(CClassConstSP clazz,
                                     const string& outputName,
                                     double        shiftSize):
    ScalarShift(clazz, outputName, shiftSize){}

/** constructor with no shift size */
RecoveryShiftBase::RecoveryShiftBase(CClassConstSP clazz,
                                     const string& outputName):
    ScalarShift(clazz, outputName){}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity. */
double RecoveryShiftBase::divisor() const{
    static const string routine("RecoveryShiftBase::divisor");
    double shiftSize;
    try{
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(routine, "Shift size is zero");
        }
        return (shiftSize/0.01);

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP RecoveryShiftBase::shiftInterface() const{
    return IShift::TYPE;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP RecoveryShiftBase::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this sensitivity's
    Shift interface */
bool RecoveryShiftBase::nameMatches(const OutputName&  name,
                                    IObjectConstSP     obj){
    // cast obj to RecoveryShiftBase::Shift and then invoke name method
    const IShift& recoveryTweakObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(recoveryTweakObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    sensitivity's Shift interface */
void RecoveryShiftBase::appendName(OutputNameArray&   namesList,
                                   IObjectConstSP     obj){
    // cast obj to RecoveryShiftBase::Shift and then invoke name method
    const IShift& RecoveryShiftBaseObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(
        new OutputName(RecoveryShiftBaseObj.sensName(this)));
    namesList.push_back(outputName);
}

bool RecoveryShiftBase::shift(IObjectSP obj) {
    // cast obj to RecoveryShiftBase::Shift and then invoke shift method
    IShift& RecoveryShiftBaseObj = dynamic_cast<IShift&>(*obj);
    return RecoveryShiftBaseObj.sensShift(this);
}

void RecoveryShiftBase::restore(IObjectSP obj) {
    // cast obj to RecoveryShiftBase::Shift and then invoke restore method
    IRestorableShift& RecoveryShiftBaseObj = 
        dynamic_cast<IRestorableShift&>(*obj);
    RecoveryShiftBaseObj.sensRestore(this);
}

void RecoveryShiftBase::negateShift()
{
    setShiftSize(-getShiftSize());
}

/** Invoked when class is 'loaded' */
void RecoveryShiftBase::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RecoveryShiftBase, clazz);
    SUPERCLASS(ScalarShift);
    IMPLEMENTS(Additive);
    // no fields
    // register how to build our sensitivity
};

CClassConstSP const RecoveryShiftBase::TYPE = CClass::registerClassLoadMethod(
    "RecoveryShiftBase", typeid(RecoveryShiftBase), load);

CClassConstSP const RecoveryShiftBase::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "RecoveryShiftBase::Shift", typeid(RecoveryShiftBase::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(RecoveryShiftBase::IRestorableShift, clazz);
    EXTENDS(RecoveryShiftBase::IShift);
}

CClassConstSP const RecoveryShiftBase::IRestorableShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "RecoveryShiftBase::RestorableShift", 
    typeid(RecoveryShiftBase::IRestorableShift),
    restorableShiftLoad);

DRLIB_END_NAMESPACE
