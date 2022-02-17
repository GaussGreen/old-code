/**
 * @file CreditIndexBasisRhoParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/CreditIndexBasisParallel.hpp"
#include "edginc/CreditIndexBasisRhoParallel.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CreditIndexBasisRhoParallel::NAME = "CREDIT_INDEX_BASIS_RHO_PARALLEL";
const double ONE_BASIS_POINT = 0.0001;
const double CreditIndexBasisRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

CreditIndexBasisRhoParallel::CreditIndexBasisRhoParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv CreditIndexBasisRhoParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CreditIndexBasisParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 ONE_BASIS_POINT);
}

CreditIndexBasisRhoParallel::~CreditIndexBasisRhoParallel() {}

void CreditIndexBasisRhoParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CreditIndexBasisRhoParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CreditIndexBasisRhoParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CreditIndexBasisRhoParallel>(), 
                                new CreditIndexBasisRhoParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CreditIndexBasisParallel>::TYPE);
}

CClassConstSP const CreditIndexBasisRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "CreditIndexBasisRhoParallel", typeid(CreditIndexBasisRhoParallel), load);

bool CreditIndexBasisRhoParallelLinkIn() {
    return CreditIndexBasisRhoParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

