/**
 * @file OptionStrikeRhoPointwise.cpp
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
#include "edginc/OptionStrikeRhoPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/IScalarDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

const string OptionStrikeRhoPointwise::NAME = "OPTION_STRIKE_RHO_POINTWISE";
const double OptionStrikeRhoPointwise::DEFAULT_SHIFT = RhoPointwise::DEFAULT_SHIFT;
const double OptionStrikeRhoPointwise::SENSITIVITY_UNIT =RhoPointwise::SENSITIVITY_UNIT;

OptionStrikeRhoPointwise::OptionStrikeRhoPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv OptionStrikeRhoPointwise::deriv() const {
    return Deriv(IResultsFunction::outputRequest(OutputRequest::OPTION_ON_CONVERTIBLE_STRIKE),
                 RiskProperty<RatePointwise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

OptionStrikeRhoPointwise::~OptionStrikeRhoPointwise() {}

void OptionStrikeRhoPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(OptionStrikeRhoPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<OptionStrikeRhoPointwise>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<OptionStrikeRhoPointwise>(), 
                                new OptionStrikeRhoPointwise(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<RatePointwise>::TYPE);
}

CClassConstSP const OptionStrikeRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "OptionStrikeRhoPointwise", typeid(OptionStrikeRhoPointwise), load);

bool OptionStrikeRhoPointwiseLinkIn() {
    return OptionStrikeRhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

