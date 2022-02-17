/**
 * @file LossMappingSens.cpp

	Author   : Juan Carlos Porras      Dec 06
	Definitions for std dev tweak to loss mapping
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/LossMappingSens.hpp"
#include "edginc/StdDevForLossMapping.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"

DRLIB_BEGIN_NAMESPACE

const string LossMappingSens::NAME = "LOSS_MAPPING_SENS";
const double ONE = 1.0;
const double LossMappingSens::DEFAULT_SHIFT = ONE;

LossMappingSens::LossMappingSens(double shiftSize,
                                 const string& packetName):
    ScalarRiskPropertySensitivity(TYPE, packetName, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv LossMappingSens::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<StdDevForLossMapping>::SP(),
                 IScalarDerivative::oneSided(),
                 ONE);
}

LossMappingSens::~LossMappingSens() {}

static IScalarRiskPropertyConstSP builtinProp(IObjectConstSP) {
    return RiskProperty<StdDevForLossMapping>::SP();
}

void LossMappingSens::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LossMappingSens, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<LossMappingSens>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<LossMappingSens>(), 
                                new LossMappingSens(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<StdDevForLossMapping>::TYPE);

    // We use this as a test case for the "generic greeks" internal override mechanism:
    // see testing/flexiblegreeksinp/internal-sensitivity.xml

    FieldSensitivityDefinition::registerBuiltin(
        "1-sided LossMappingSens",
        RQFFactory_1packet<LossMappingSens>::newOne);

    // Also make RiskProperty<LossMappingSens> available through the lower-level
    // property-override mechanism: see same test

    IFieldRiskPropertyDefinition::registerScalarBuiltin(
       "Loss Mapping", builtinProp);
}

CClassConstSP const LossMappingSens::TYPE = CClass::registerClassLoadMethod(
    "LossMappingSens", typeid(LossMappingSens), load);

bool LossMappingSensLinkIn() {
    return LossMappingSens::TYPE != NULL;
}

DRLIB_END_NAMESPACE


