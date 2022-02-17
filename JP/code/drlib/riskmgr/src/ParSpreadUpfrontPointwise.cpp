/**
 * @file ParSpreadUpfrontPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/ParSpreadUpfronts.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/ParSpreadUpfrontPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

const string ParSpreadUpfrontPointwise::NAME = "PAR_SPREAD_UPFRONT_POINTWISE";
const double ParSpreadUpfrontPointwise::DEFAULT_SHIFT = 0.0001;
const double SENSITIVITY_UNIT = 0.0001;

ParSpreadUpfrontPointwise::ParSpreadUpfrontPointwise(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, name)
{}

ParSpreadUpfrontPointwise::ParSpreadUpfrontPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv ParSpreadUpfrontPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<ParSpreadUpfronts>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

ParSpreadUpfrontPointwise::~ParSpreadUpfrontPointwise() {}

void ParSpreadUpfrontPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(ParSpreadUpfrontPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<ParSpreadUpfrontPointwise>::iObject);
    SensitivityFactory::addSens(ParSpreadUpfrontPointwise::NAME, 
                                new GenericSensitivityFactory<ParSpreadUpfrontPointwise>(), 
                                new ParSpreadUpfrontPointwise(ParSpreadUpfrontPointwise::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<ParSpreadUpfronts>::TYPE);
}

CClassConstSP const ParSpreadUpfrontPointwise::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadUpfrontPointwise", typeid(ParSpreadUpfrontPointwise), load);

bool ParSpreadUpfrontPointwiseLinkIn() {
    return ParSpreadUpfrontPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
