/**
 * @file CorrSwapSamplingAdjTweak.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/CorrSwapSamplingAdjTweak.hpp"
#include "edginc/CorrSwapSamplingAdjAbsolute.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CorrSwapSamplingAdjTweak::NAME = "CORR_SWAP_SAMPLING_ADJ_TWEAK";
const double CorrSwapSamplingAdjTweak::DEFAULT_SHIFT = 0.001;
const double CorrSwapSamplingAdjTweak::SENSITIVITY_UNIT = 0.01;

CorrSwapSamplingAdjTweak::CorrSwapSamplingAdjTweak(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv CorrSwapSamplingAdjTweak::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CorrSwapSamplingAdjAbsolute>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CorrSwapSamplingAdjTweak::~CorrSwapSamplingAdjTweak() {}

void CorrSwapSamplingAdjTweak::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CorrSwapSamplingAdjTweak, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CorrSwapSamplingAdjTweak>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CorrSwapSamplingAdjTweak>(), 
                                new CorrSwapSamplingAdjTweak(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CorrSwapSamplingAdjAbsolute>::TYPE);
}

CClassConstSP const CorrSwapSamplingAdjTweak::TYPE = CClass::registerClassLoadMethod(
    "CorrSwapSamplingAdjTweak", typeid(CorrSwapSamplingAdjTweak), load);

bool CorrSwapSamplingAdjTweakLinkIn() {
    return CorrSwapSamplingAdjTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

