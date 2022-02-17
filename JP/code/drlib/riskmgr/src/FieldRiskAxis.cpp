/**
 * @file FieldRiskAxis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/IAbstractRiskProperty.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/IFieldTweak.hpp"
#include "edginc/FieldRiskAxis.hpp"
#include "edginc/FieldTweakHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE


FieldRiskAxis::FieldRiskAxis(IAbstractRiskPropertyConstSP property,
                             OutputNameConstSP name,
                             IFieldTweakConstSP tweak,
                             bool absoluteDistance):
    CObject(TYPE),
    property(property),
    name(name),
    tweak(tweak),
    absoluteDistance(absoluteDistance)
{
    try {
        validate();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

FieldRiskAxis::FieldRiskAxis():
    CObject(TYPE),
    absoluteDistance(absoluteDistance)
{}

void FieldRiskAxis::validate() {
    if (!!name &&
          !MarketObject::TYPE->isAssignableFrom(property->subjectInterface())) {
        throw ModelException(
            "You can't specify a market-name-specific tweak to a class (" +
            property->subjectInterface()->getName() + ") "
            "which isn't a MarketObject");
    }
}

void FieldRiskAxis::validatePop2Object() {
    try {
        validate();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

FieldRiskAxis::~FieldRiskAxis() {}

IAbstractRiskPropertyConstSP FieldRiskAxis::abstractProperty() const {
    return property;
}

OutputNameConstSP FieldRiskAxis::marketDataName() const {
    return name;
}

IHypothesisConstSP FieldRiskAxis::hypothesis(double coeff) const {
    try {
        return coeff == 0 && tweak->zeroIsNoop() ?
            IHypothesis::noop() :
            IHypothesisConstSP(new FieldTweakHypothesis(
                FieldRiskAxisConstSP::attachToRef(this), coeff));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

RiskAxisConstSP FieldRiskAxis::frozen() const {
    return tweak->legacyRepresentation(property, name, absoluteDistance);
}

string FieldRiskAxis::toString() const {
    return property->subjectInterface()->getName() + " " +
           (!name ? "" : name->toString() + " ") +
           tweak->toString();
}



IObject* FieldRiskAxis::emptyShell() {
  return new FieldRiskAxis();
}

void FieldRiskAxis::load(CClassSP& clazz) {
  REGISTER(FieldRiskAxis, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskAxis);
  FIELD(property, "property");
  FIELD(name, "name");
  FIELD_MAKE_OPTIONAL(name);
  FIELD(tweak, "tweak");
  FIELD(absoluteDistance, "absoluteDistance");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const FieldRiskAxis::TYPE = CClass::registerClassLoadMethod(
  "FieldRiskAxis", typeid(FieldRiskAxis), FieldRiskAxis::load);

DRLIB_END_NAMESPACE
