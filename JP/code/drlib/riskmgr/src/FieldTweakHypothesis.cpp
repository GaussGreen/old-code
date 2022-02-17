/**
 * @file FieldTweakHypothesis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/IFieldTweak.hpp"
#include "edginc/IAbstractRiskProperty.hpp"
#include "edginc/FieldRiskAxis.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/FieldTweakHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

FieldTweakHypothesis::FieldTweakHypothesis(FieldRiskAxisConstSP axis,
                                           double coefficient):
    AbstractPropertyTweakHypothesis(TYPE),
    _axis(axis),
    _coefficient(coefficient)
{}

FieldTweakHypothesis::FieldTweakHypothesis():
    AbstractPropertyTweakHypothesis(TYPE),
    _coefficient(0.)
{}

FieldTweakHypothesis::~FieldTweakHypothesis() {}

IFieldTweakConstSP FieldTweakHypothesis::tweak() const {
    if (!_tweak) {
        _tweak = _axis->tweak->scaled(_coefficient);
    }

    return _tweak;
}

IRiskAxisConstSP FieldTweakHypothesis::axis() const {
    return _axis;
}

TweakOutcome FieldTweakHypothesis::sensShift(IObjectSP root) const {
    try {
        double distance = tweak()->apply(root);
        return TweakOutcome(_axis->absoluteDistance ? distance : _coefficient,
                            false);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

CClassConstSP FieldTweakHypothesis::shiftableInterface() const {
    return _axis->property->subjectInterface();
}

CClassConstSP FieldTweakHypothesis::restorableInterface() const {
    return CClassConstSP();
}

void FieldTweakHypothesis::sensRestore(IObjectSP object) const {
    throw ModelException(__FUNCTION__, "not implemented");
}

OutputNameConstSP FieldTweakHypothesis::getMarketDataName() const {
    // should just be fieldAxis()->marketDataName() --- _marketDataName
    // override is for temporary compatibility with
    // Control::getCurrentSensitivity()

    return !_marketDataName ? _axis->marketDataName() : _marketDataName;
}

IObjectConstSP FieldTweakHypothesis::getQualifier() const {
    IArrayConstSP qs = tweak()->qualifiers();
    return !qs || (qs->getLength() != 1) ? IObjectConstSP()
                                         : qs->get(0);
}

double FieldTweakHypothesis::axisCoefficient() const {
    return _coefficient;
}

string FieldTweakHypothesis::sensName(IObjectConstSP object) const {
    try {
        const INamedObject* m =
            dynamic_cast<const INamedObject*>(object.get());
        return m ? m->getName() : !object ? "(null)" : object->toString();
    }
    catch (exception& e) {
        throw ModelException(e, "FieldTweakHypothesis::sensName()");
    }
}

void FieldTweakHypothesis::setAxisCoefficient(double c) {
    // compat mechanism for (?) ImpliedScalarShift
    _coefficient = c;
    _tweak.reset();
}

string FieldTweakHypothesis::toString() const {
    IObjectConstSP q = getQualifier();

    return Format::toString(axisCoefficient()) + " tweak of " +
           shiftableInterface()->getName() + " " +
           (!getMarketDataName() ? string()
                                 : getMarketDataName()->toString() + "'s ") +
           tweak()->toString();
}

void FieldTweakHypothesis::setMarketDataName(OutputNameConstSP n) {
    _marketDataName = n;
}

IObject* FieldTweakHypothesis::emptyShell() {
  return new FieldTweakHypothesis();
}

void FieldTweakHypothesis::load(CClassSP& clazz) {
  REGISTER(FieldTweakHypothesis, clazz);
  SUPERCLASS(AbstractPropertyTweakHypothesis);
  FIELD(_axis, "_axis");
  FIELD(_coefficient, "_coefficient");
  FIELD(_tweak, "_tweak");
  FIELD_MAKE_TRANSIENT(_tweak);
  FIELD(_marketDataName, "_marketDataName");
  FIELD_MAKE_TRANSIENT(_marketDataName);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const FieldTweakHypothesis::TYPE = CClass::registerClassLoadMethod(
  "FieldTweakHypothesis", typeid(FieldTweakHypothesis), FieldTweakHypothesis::load);

DRLIB_END_NAMESPACE
