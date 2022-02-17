/**
 * @file CorrSwapBasisAdjParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/CorrSwapBasisAdjParallel.hpp"
#include "edginc/CorrSwapBasisAdjAbsolute.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CorrSwapBasisAdjParallel::NAME = "CORR_SWAP_BASIS_ADJ_PARALLEL";
const double CorrSwapBasisAdjParallel::DEFAULT_SHIFT = 0.001;
const double CorrSwapBasisAdjParallel::SENSITIVITY_UNIT = 0.01;

CorrSwapBasisAdjParallel::CorrSwapBasisAdjParallel(double shiftSize):
    AllNamesRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

AllNamesRiskPropertySensitivity::Deriv CorrSwapBasisAdjParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CorrSwapBasisAdjAbsolute>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CorrSwapBasisAdjParallel::~CorrSwapBasisAdjParallel() {}

void CorrSwapBasisAdjParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CorrSwapBasisAdjParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(AllNamesRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CorrSwapBasisAdjParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CorrSwapBasisAdjParallel>(), 
                                new CorrSwapBasisAdjParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CorrSwapBasisAdjAbsolute>::TYPE);
}

CClassConstSP const CorrSwapBasisAdjParallel::TYPE = CClass::registerClassLoadMethod(
    "CorrSwapBasisAdjParallel", typeid(CorrSwapBasisAdjParallel), load);

bool CorrSwapBasisAdjParallelLinkIn() {
    return CorrSwapBasisAdjParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
