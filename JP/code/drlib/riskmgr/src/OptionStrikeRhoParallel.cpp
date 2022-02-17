/**
 * @file OptionStrikeRhoParallel.cpp
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
#include "edginc/OptionStrikeRhoParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/RhoParallel.hpp"

DRLIB_BEGIN_NAMESPACE

const string OptionStrikeRhoParallel::NAME = "OPTION_STRIKE_RHO_PARALLEL";
const double OptionStrikeRhoParallel::DEFAULT_SHIFT = RhoParallel::DEFAULT_SHIFT;
const double OptionStrikeRhoParallel::SENSITIVITY_UNIT =RhoParallel::DEFAULT_SHIFT;

OptionStrikeRhoParallel::OptionStrikeRhoParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv OptionStrikeRhoParallel::deriv() const {
    return Deriv(IResultsFunction::outputRequest(OutputRequest::OPTION_ON_CONVERTIBLE_STRIKE),
                 RiskProperty<RateParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

OptionStrikeRhoParallel::~OptionStrikeRhoParallel() {}

void OptionStrikeRhoParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OptionStrikeRhoParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<OptionStrikeRhoParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<OptionStrikeRhoParallel>(), 
                                new OptionStrikeRhoParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<RateParallel>::TYPE);
}

CClassConstSP const OptionStrikeRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "OptionStrikeRhoParallel", typeid(OptionStrikeRhoParallel), load);

bool OptionStrikeRhoParallelLinkIn() {
    return OptionStrikeRhoParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

