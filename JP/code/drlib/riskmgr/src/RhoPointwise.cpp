/**
 * @file RhoPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

const string RhoPointwise::NAME = "RHO_POINTWISE";
const double RhoPointwise::DEFAULT_SHIFT = 0.0001;
const double RhoPointwise::SENSITIVITY_UNIT = 0.0001;

RhoPointwise::RhoPointwise(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, name)
{}

RhoPointwise::RhoPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv RhoPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<RatePointwise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

RhoPointwise::~RhoPointwise() {}

void RhoPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(RhoPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<RhoPointwise>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<RhoPointwise>(), 
                                new RhoPointwise(RhoPointwise::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<RatePointwise>::TYPE);
}

CClassConstSP const RhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "RhoPointwise", typeid(RhoPointwise), load);

bool RhoPointwiseLinkIn() {
    return RhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

