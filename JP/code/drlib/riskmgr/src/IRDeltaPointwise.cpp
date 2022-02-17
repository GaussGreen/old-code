/**
 * @file IRDeltaPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/IRDeltaPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IRRatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

const string IRDeltaPointwise::NAME = "IRDELTA_POINTWISE";
const double IRDeltaPointwise::DEFAULT_SHIFT = 0.0001;
const double IRDeltaPointwise::SENSITIVITY_UNIT = 0.0001;

IRDeltaPointwise::IRDeltaPointwise(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, name)
{}

IRDeltaPointwise::IRDeltaPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv IRDeltaPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<IRRatePointwise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

IRDeltaPointwise::~IRDeltaPointwise() {}

void IRDeltaPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(IRDeltaPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<IRDeltaPointwise>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<IRDeltaPointwise>(), 
                                new IRDeltaPointwise(IRDeltaPointwise::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<IRRatePointwise>::TYPE);
}

CClassConstSP const IRDeltaPointwise::TYPE = CClass::registerClassLoadMethod(
    "IRDeltaPointwise", typeid(IRDeltaPointwise), load);

bool IRDeltaPointwiseLinkIn() {
    return IRDeltaPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

