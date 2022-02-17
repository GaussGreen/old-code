/**
 * @file ParSpreadRhoParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/ParSpreadParallel.hpp"
#include "edginc/ParSpreadRhoParallel.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"

DRLIB_BEGIN_NAMESPACE

const string ParSpreadRhoParallel::NAME = "PAR_SPREAD_RHO_PARALLEL";
const double ONE_BASIS_POINT = 0.0001;
const double ParSpreadRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

ParSpreadRhoParallel::ParSpreadRhoParallel(double shiftSize,
                                           const string& packetName):
    ScalarRiskPropertySensitivity(TYPE, packetName, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv ParSpreadRhoParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<ParSpreadParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 ONE_BASIS_POINT);
}

ParSpreadRhoParallel::~ParSpreadRhoParallel() {}

static IScalarRiskPropertyConstSP builtinProp(IObjectConstSP) {
    return RiskProperty<ParSpreadParallel>::SP();
}

void ParSpreadRhoParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(ParSpreadRhoParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<ParSpreadRhoParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<ParSpreadRhoParallel>(), 
                                new ParSpreadRhoParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<ParSpreadParallel>::TYPE);

    // We use this as a test case for the "generic greeks" internal override mechanism:
    // see testing/flexiblegreeksinp/internal-sensitivity.xml

    FieldSensitivityDefinition::registerBuiltin(
        "1-sided ParSpreadRhoParallel",
        RQFFactory_1packet<ParSpreadRhoParallel>::newOne);

    // Also make RiskProperty<ParSpreadParallel> available through the lower-level
    // property-override mechanism: see same test

    IFieldRiskPropertyDefinition::registerScalarBuiltin(
       "par spread parallel", builtinProp);
}

CClassConstSP const ParSpreadRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadRhoParallel", typeid(ParSpreadRhoParallel), load);

bool ParSpreadRhoParallelLinkIn() {
    return ParSpreadRhoParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

