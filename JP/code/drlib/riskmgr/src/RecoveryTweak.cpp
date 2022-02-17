//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoveryTweak.cpp
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
#include "edginc/RecoveryTweak.hpp"

DRLIB_BEGIN_NAMESPACE

//// Sensitivity to CDS Recovery

/** constructor with explicit shift size */
RecoveryTweak::RecoveryTweak(double shiftSize):
    RecoveryShiftBase(TYPE, NAME, shiftSize)
{}

/** Returns the new recovery level given the original one */
double RecoveryTweak::applyShift(double unadjRecovery) {
    double newRecovery = unadjRecovery + getShiftSize();
    //if the recovery has gone out of bounds
    //then negate the shift and re-apply
    if ((newRecovery < 0.0) || (newRecovery > 1.0))
    {
        negateShift();
        newRecovery = unadjRecovery + getShiftSize();
    }
    //may still be out of bounds but let the caller
    //deal with that...
    return newRecovery;
}

/** Returns the original recovery level given the adjusted one */
double RecoveryTweak::undoShift(double adjRecovery) {
    return (adjRecovery - getShiftSize());
}


Sensitivity* RecoveryTweak::Factory::createDefault() {
    return new RecoveryTweak(RecoveryTweak::DEFAULT_SHIFT);
}


Sensitivity* RecoveryTweak::Factory::createScalar(double shiftSize) {
    return new RecoveryTweak(shiftSize);
}


/** for reflection */
RecoveryTweak::RecoveryTweak():
    RecoveryShiftBase(TYPE, NAME, DEFAULT_SHIFT)
{}


IObject* RecoveryTweak::defaultRecoveryTweak() {
    return new RecoveryTweak();
}


/** Invoked when SmallClass is 'loaded' */
void RecoveryTweak::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RecoveryTweak, clazz);
    SUPERCLASS(RecoveryShiftBase);
    EMPTY_SHELL_METHOD(defaultRecoveryTweak);
    // no fields
    // register how to build our sensitivity
    SensitivityFactory::addSens(RecoveryTweak::NAME,
                                new Factory(),
                                new RecoveryTweak(DEFAULT_SHIFT),
                                IShift::TYPE);
}


CClassConstSP const RecoveryTweak::TYPE = CClass::registerClassLoadMethod(
    "RecoveryTweak", typeid(RecoveryTweak), load);
const string RecoveryTweak::NAME = "RECOVERY_TWEAK";
const double RecoveryTweak::DEFAULT_SHIFT = 0.01;

/** Included in RiskMgrLib::linkInClasses() to force link to include this */
bool RecoveryTweakLinkIn() {
    return RecoveryTweak::TYPE != NULL;
}


DRLIB_END_NAMESPACE
