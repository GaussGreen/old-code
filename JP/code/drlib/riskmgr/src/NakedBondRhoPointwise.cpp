/**
 * @file NakedBondRhoPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NakedBondRhoPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/IScalarDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

const string NakedBondRhoPointwise::NAME = "NAKED_BOND_RHO_POINTWISE";
const double NakedBondRhoPointwise::DEFAULT_SHIFT = RhoPointwise::DEFAULT_SHIFT;
const double NakedBondRhoPointwise::SENSITIVITY_UNIT =RhoPointwise::SENSITIVITY_UNIT;

NakedBondRhoPointwise::NakedBondRhoPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv NakedBondRhoPointwise::deriv() const {
    return Deriv(IResultsFunction::outputRequest(OutputRequest::NAKED_BOND_PRICE),
                 RiskProperty<RatePointwise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

NakedBondRhoPointwise::~NakedBondRhoPointwise() {}

void NakedBondRhoPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(NakedBondRhoPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<NakedBondRhoPointwise>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<NakedBondRhoPointwise>(), 
                                new NakedBondRhoPointwise(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<RatePointwise>::TYPE);
}

CClassConstSP const NakedBondRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "NakedBondRhoPointwise", typeid(NakedBondRhoPointwise), load);

bool NakedBondRhoPointwiseLinkIn() {
    return NakedBondRhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

