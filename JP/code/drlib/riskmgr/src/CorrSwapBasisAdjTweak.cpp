/**
 * @file CorrSwapBasisAdjTweak.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/CorrSwapBasisAdjTweak.hpp"
#include "edginc/CorrSwapBasisAdjAbsolute.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CorrSwapBasisAdjTweak::NAME = "CORR_SWAP_BASIS_ADJ_TWEAK";
const double CorrSwapBasisAdjTweak::DEFAULT_SHIFT = 0.001;
const double CorrSwapBasisAdjTweak::SENSITIVITY_UNIT = 0.01;

CorrSwapBasisAdjTweak::CorrSwapBasisAdjTweak(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv CorrSwapBasisAdjTweak::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CorrSwapBasisAdjAbsolute>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CorrSwapBasisAdjTweak::~CorrSwapBasisAdjTweak() {}

void CorrSwapBasisAdjTweak::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CorrSwapBasisAdjTweak, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CorrSwapBasisAdjTweak>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CorrSwapBasisAdjTweak>(), 
                                new CorrSwapBasisAdjTweak(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CorrSwapBasisAdjAbsolute>::TYPE);
}

CClassConstSP const CorrSwapBasisAdjTweak::TYPE = CClass::registerClassLoadMethod(
    "CorrSwapBasisAdjTweak", typeid(CorrSwapBasisAdjTweak), load);

bool CorrSwapBasisAdjTweakLinkIn() {
    return CorrSwapBasisAdjTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

