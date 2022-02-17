/**
 * @file VegaPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/VegaPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

const string VegaPointwise::NAME = "VEGA_POINTWISE";
const double VegaPointwise::DEFAULT_SHIFT = 0.005;
const double SENSITIVITY_UNIT = 0.01;

VegaPointwise::VegaPointwise(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, name)
{}

VegaPointwise::VegaPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv VegaPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<VolPointwise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

VegaPointwise::~VegaPointwise() {}

void VegaPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VegaPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VegaPointwise>::iObject);
    SensitivityFactory::addSens(VegaPointwise::NAME, 
                                new GenericSensitivityFactory<VegaPointwise>(), 
                                new VegaPointwise(VegaPointwise::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<VolPointwise>::TYPE);
}

CClassConstSP const VegaPointwise::TYPE = CClass::registerClassLoadMethod(
    "VegaPointwise", typeid(VegaPointwise), load);

bool VegaPointwiseLinkIn() {
    return VegaPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
