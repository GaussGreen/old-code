/**
 * @file NakedBondRhoParallel.cpp
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
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NakedBondRhoParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/RhoParallel.hpp"

DRLIB_BEGIN_NAMESPACE

const string NakedBondRhoParallel::NAME = "NAKED_BOND_RHO_PARALLEL";
const double NakedBondRhoParallel::DEFAULT_SHIFT = RhoParallel::DEFAULT_SHIFT;
const double NakedBondRhoParallel::SENSITIVITY_UNIT =RhoParallel::SENSITIVITY_UNIT;

NakedBondRhoParallel::NakedBondRhoParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv NakedBondRhoParallel::deriv() const {
    return Deriv(IResultsFunction::outputRequest(OutputRequest::NAKED_BOND_PRICE),
                 RiskProperty<RateParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

NakedBondRhoParallel::~NakedBondRhoParallel() {}

void NakedBondRhoParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(NakedBondRhoParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<NakedBondRhoParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<NakedBondRhoParallel>(), 
                                new NakedBondRhoParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<RateParallel>::TYPE);
}

CClassConstSP const NakedBondRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "NakedBondRhoParallel", typeid(NakedBondRhoParallel), load);

bool NakedBondRhoParallelLinkIn() {
    return NakedBondRhoParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

