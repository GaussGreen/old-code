//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoveryShift.cpp
//
//   Description : Scenario for adjusting CDS Recovery
//
//   Author      : Mark A Robson
//
//   Date        : 26 July 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RecoveryShift.hpp"

DRLIB_BEGIN_NAMESPACE


RecoveryShift::~RecoveryShift()
{}


RecoveryShift::RecoveryShift(CClassConstSP clazz,
                             const string& outputName):
    RecoveryShiftBase(clazz, outputName), 
    relativeShift(false), 
    recoveryFloor(0.0), 
    recoveryCap(1.0)
{}

RecoveryShift::RecoveryShift(double shiftSize,
                             bool   relativeShift):
    RecoveryShiftBase(TYPE, NAME, shiftSize), 
    relativeShift(relativeShift), 
    recoveryFloor(0.0), 
    recoveryCap(1.0)
{}


/** for reflection */
RecoveryShift::RecoveryShift():
    RecoveryShiftBase(TYPE, NAME), 
    relativeShift(false), recoveryFloor(0.0), recoveryCap(1.0)
{}

IObject* RecoveryShift::defaultConstructor(){
    return new RecoveryShift();
}

/** check the inputs */
void RecoveryShift::validatePop2Object(){
    if (recoveryFloor > recoveryCap){
        throw ModelException("RecoveryShift::validatePop2Object",
                             "recoveryFloor > recoveryCap");
    }
}

/** Returns the new recovery level given the original one */
double RecoveryShift::applyShift(double unadjRecovery){
    setInitialValue(unadjRecovery);
    double shiftSize = getShiftSize();
    // make adjustment
    double recovery = relativeShift?
        unadjRecovery * (1.0+shiftSize): (unadjRecovery+shiftSize);
    // then apply cap/floor
    recovery = Maths::max(recoveryFloor, Maths::min(recoveryCap, recovery));
    return recovery;
}

/** Returns the original recovery level given the adjusted one */
double RecoveryShift::undoShift(double adjRecovery){
    return getInitialValue();
}

/** overridden: pretty meaningless for a scenario */
double RecoveryShift::divisor() const{
    return 1.0;
}

/** overridden: fails */
void RecoveryShift::calculate(TweakGroup* tweakGroup,
                              CResults*    results) {
    throw ModelException("RecoveryShift::calculate", "Not supported");
}
    

/** Invoked when SmallClass is 'loaded' */
void RecoveryShift::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RecoveryShift, clazz);
    SUPERCLASS(RecoveryShiftBase);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(relativeShift, "True: multiplicative shift of 1+shiftSize"
                 " else additive shift of shiftSize");
    FIELD_MAKE_OPTIONAL(relativeShift);
    FIELD(recoveryFloor, "Floor for recovery, default is 0.0");
    FIELD_MAKE_OPTIONAL(recoveryFloor);
    FIELD(recoveryCap, "Cap for recovery, default is 1.0");
    FIELD_MAKE_OPTIONAL(recoveryCap);
}

CClassConstSP const RecoveryShift::TYPE = CClass::registerClassLoadMethod(
    "RecoveryShift", typeid(RecoveryShift), load);
const string RecoveryShift::NAME = "RECOVERY_SHIFT";

/** Included in RiskMgrLib::linkInClasses() to force link to include this */
bool RecoveryShiftLinkIn() {
    return RecoveryShift::TYPE != NULL;
}


DRLIB_END_NAMESPACE
