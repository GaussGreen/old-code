/**
 * @file TweakOutcome.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TweakOutcome.hpp"

DRLIB_BEGIN_NAMESPACE

const double noOldValue = -1.234e56;

TweakOutcome::TweakOutcome(double distance, bool tweakMembers):
    CObject(TYPE),
    _oldValue(noOldValue),
    _newValue(distance),
    _tweakMembers(tweakMembers)
{}

TweakOutcome::TweakOutcome(double oldValue, double newValue, bool tweakMembers):
    CObject(TYPE),
    _oldValue(oldValue),
    _newValue(newValue),
    _tweakMembers(tweakMembers)
{}

TweakOutcome::~TweakOutcome() {}

bool TweakOutcome::hasOldValue() const {
    return _oldValue != noOldValue;
}

double TweakOutcome::oldValue() const {
    if (_oldValue == noOldValue)
        throw ModelException(
            "TweakOutcome::oldValue",
            "No meaningful 'value before tweaking' is available");
    return _oldValue;
}

double TweakOutcome::distance() const {
    return _oldValue == noOldValue ? _newValue : _newValue - _oldValue;
}

bool TweakOutcome::tweakMembers() const {
    return _tweakMembers;
}

void TweakOutcome::setTweakMembers(bool tweakMembers) {
    _tweakMembers = tweakMembers;
}

static IObject* defaultTweakOutcome() {
    return new TweakOutcome(0., false);
}

void TweakOutcome::load(CClassSP& clazz) {
    REGISTER(TweakOutcome, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultTweakOutcome);
    FIELD(_oldValue, "_oldValue");
    FIELD(_newValue, "_newValue");
    FIELD(_tweakMembers, "_tweakMembers");
}

CClassConstSP const TweakOutcome::TYPE = CClass::registerClassLoadMethod(
    "TweakOutcome", typeid(TweakOutcome), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
