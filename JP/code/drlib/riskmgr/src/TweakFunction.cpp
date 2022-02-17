/**
 * @file TweakFunction.cpp
 */

#include <float.h>
#include <string.h>
#include <stdlib.h>
#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/Maths.hpp"
#include "edginc/TweakFunction.hpp"

DRLIB_BEGIN_NAMESPACE

TweakFunction::TweakFunction(CClassConstSP type):
    CObject(type)
{}

TweakFunction::~TweakFunction() {}

double TweakFunction::operator ()(double argument, double x,
                                  const Range& range, bool clip) const {
    TRACE_METHOD;
    try {
        if (!Range::variableIsInRange(range, x)) {
            throw ModelException("Value is already outside its valid range "
                                 "before tweaking");
        }

        double y = value(argument, x, range);

        if (clip) {
            if (!range.getLower().isInfinite()) {
                double l = range.getLower().getValue();
                if (y < l) {
                    y = range.getLower().isClosedBracket() ? l : .5 * (x + l);
                }
            }
            if (!range.getUpper().isInfinite()) {
                double u = range.getUpper().getValue();
                if (y > u) {
                    y = range.getUpper().isClosedBracket() ? u : .5 * (x + u);
                }
            }
        }

        if (!Range::variableIsInRange(range, y)) {
            throw ModelException("Tweaked value " + Format::toString(y) +
                                 " would go out of its valid range");
        }

        TRACE("Old value was " << x << ", new value is " << y);

        return y;
    }
    catch (exception& e) {
        throw ModelException(e, "TweakFunction::operator ()()",
            "Tweaking value " + Format::toString(x) + " (valid range " +
            range.toString() + ")");
    }
}

string TweakFunction::repr() const {
    return getClass()->getName();
}

TweakFunctionConstSP TweakFunction::fromRepr(string repr) {
    try {
        return TweakFunctionSP::dynamicCast(
             IObjectSP(CClass::forName(repr)->newInstance()));
    }
    catch (exception& e) {
        throw ModelException(e, "TweakFunction::fromRepr()");
    }
}

void TweakFunction::load(CClassSP& clazz) {
    REGISTER(TweakFunction, clazz);
    SUPERCLASS(CObject);
}

CClassConstSP const TweakFunction::TYPE = CClass::registerClassLoadMethod(
    "TweakFunction", typeid(TweakFunction), load);

struct AdditiveTweakFunction: public TweakFunction {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(AdditiveTweakFunction, clazz);
        SUPERCLASS(TweakFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<AdditiveTweakFunction>::iObject);
    }

    AdditiveTweakFunction(): TweakFunction(TYPE) {}

    double value(double argument, double x, const Range&) const {
        TRACE_METHOD;
        return x + argument;
    }

    bool zeroIsNoop() const { return true; }

    string toString() const { return "+"; }
};

CClassConstSP const AdditiveTweakFunction::TYPE = CClass::registerClassLoadMethod(
    "AdditiveTweakFunction", typeid(AdditiveTweakFunction), load);

TweakFunctionConstSP TweakFunction::additive() {
    static TweakFunctionConstSP it(new AdditiveTweakFunction());
    return it;
}

struct MultiplicativeTweakFunction: public TweakFunction {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(MultiplicativeTweakFunction, clazz);
        SUPERCLASS(TweakFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<MultiplicativeTweakFunction>::iObject);
    }

    MultiplicativeTweakFunction(): TweakFunction(TYPE) {}

    double value(double argument, double x, const Range&) const {
        return x * (1 + argument);
    }

    bool zeroIsNoop() const { return true; }

    string toString() const { return "×"; }
};

CClassConstSP const MultiplicativeTweakFunction::TYPE = CClass::registerClassLoadMethod(
    "MultiplicativeTweakFunction", typeid(MultiplicativeTweakFunction), load);

TweakFunctionConstSP TweakFunction::multiplicative() {
    static TweakFunctionConstSP it(new MultiplicativeTweakFunction());
    return it;
}

struct ExponentialTweakFunction: public TweakFunction {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(ExponentialTweakFunction, clazz);
        SUPERCLASS(TweakFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<ExponentialTweakFunction>::iObject);
    }

    ExponentialTweakFunction(): TweakFunction(TYPE) {}

    double value(double argument, double x, const Range&) const {
        TRACE_METHOD;
        return x * exp(argument);
    }

    bool zeroIsNoop() const { return true; }

    string toString() const { return "exp"; }
};

CClassConstSP const ExponentialTweakFunction::TYPE = CClass::registerClassLoadMethod(
    "ExponentialTweakFunction", typeid(ExponentialTweakFunction), load);

TweakFunctionConstSP TweakFunction::exponential() {
    static TweakFunctionConstSP it(new ExponentialTweakFunction());
    return it;
}

struct AdaptiveTweakFunction: public TweakFunction {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(AdaptiveTweakFunction, clazz);
        SUPERCLASS(TweakFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<AdaptiveTweakFunction>::iObject);
    }

    AdaptiveTweakFunction(): TweakFunction(TYPE) {}

    double value(double argument, double x, const Range& range) const {
        TRACE_METHOD;

        if (range.getUpper().isInfinite()) {
            return range.getLower().isInfinite() ?
                x + argument :
                range.getLower().getValue() +
                    (x - range.getLower().getValue()) * exp(argument);
        }
        else {
            if (range.getLower().isInfinite()) {
                return range.getUpper().getValue() -
                         (range.getUpper().getValue() - x) * exp(-argument);
            }
            else {
                double cap = range.getUpper().getValue();
                double floor = range.getLower().getValue();

                double l = log((cap - floor) / (x - floor) - 1);
                if (!Maths::finite(l)) {
                    throw ModelException("Numerical error doing sigmoidal tweak");
                }

                return floor + (cap - floor) * (1 / (1 + exp(l - argument)));
            }
        }
    }

    bool zeroIsNoop() const { return true; }

    string toString() const { return "+/exp/sig"; }
};

CClassConstSP const AdaptiveTweakFunction::TYPE = CClass::registerClassLoadMethod(
    "AdaptiveTweakFunction", typeid(AdaptiveTweakFunction), load);

TweakFunctionConstSP TweakFunction::adaptive() {
    static TweakFunctionConstSP it(new AdaptiveTweakFunction());
    return it;
}

// This guy is deprecated --- it's only here to support any legacy data there might
// be lurking in Pyramid

struct AutomaticTweakFunction: public TweakFunction {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(AutomaticTweakFunction, clazz);
        SUPERCLASS(TweakFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<AutomaticTweakFunction>::iObject);
    }

    AutomaticTweakFunction(): TweakFunction(TYPE) {}

    double value(double argument, double x, const Range& range) const {
        TRACE_METHOD;

        if (range.getUpper().isInfinite()) {
            return range.getLower().isInfinite() ?
                x + argument :
                range.getLower().getValue() +
                   (x - range.getLower().getValue()) * (1 + argument);
        }
        else {
            if (range.getLower().isInfinite()) {
                return range.getUpper().getValue() -
                       (range.getUpper().getValue() - x) * (1 - argument);
            }
            else {
                double cap = range.getUpper().getValue();
                double floor = range.getLower().getValue();

                double l = log((cap - floor) / (x - floor) - 1);
                if (!Maths::finite(l)) {
                    throw ModelException("Numerical error doing sigmoidal tweak");
                }

                return floor + (cap - floor) * (1 / (1 + exp(l - argument)));
            }
        }
    }

    bool zeroIsNoop() const { return true; }

    string toString() const { return "+/×/sig"; }
};

CClassConstSP const AutomaticTweakFunction::TYPE = CClass::registerClassLoadMethod(
    "AutomaticTweakFunction", typeid(AutomaticTweakFunction), load);

struct SetterTweakFunction: public TweakFunction {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(SetterTweakFunction, clazz);
        SUPERCLASS(TweakFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<SetterTweakFunction>::iObject);
    }

    SetterTweakFunction(): TweakFunction(TYPE) {}

    double value(double argument, double x, const Range&) const {
        TRACE_METHOD;
        return argument;
    }

    bool zeroIsNoop() const { return false; }

    string toString() const { return "="; }
};

CClassConstSP const SetterTweakFunction::TYPE = CClass::registerClassLoadMethod(
    "SetterTweakFunction", typeid(SetterTweakFunction), load);

TweakFunctionConstSP TweakFunction::setter() {
    static TweakFunctionConstSP it(new SetterTweakFunction());
    return it;
}

DRLIB_END_NAMESPACE
