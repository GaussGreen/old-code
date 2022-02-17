/**
 * @file CreditIndexBasisRhoPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/CreditIndexBasisPointwise.hpp"
#include "edginc/CreditIndexBasisRhoPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"

DRLIB_BEGIN_NAMESPACE

const string CreditIndexBasisRhoPointwise::NAME = "CREDIT_INDEX_BASIS_RHO_POINTWISE";
const double ONE_BASIS_POINT = 0.0001;
const double CreditIndexBasisRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

CreditIndexBasisRhoPointwise::~CreditIndexBasisRhoPointwise() {}

void CreditIndexBasisRhoPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CreditIndexBasisRhoPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CreditIndexBasisRhoPointwise>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CreditIndexBasisRhoPointwise>(), 
                                new CreditIndexBasisRhoPointwise(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CreditIndexBasisPointwise>::TYPE);
}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv CreditIndexBasisRhoPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 ParSpreadMaxTweakSize::adapted(
                     RiskProperty<CreditIndexBasisPointwise>::SP()),
                 IScalarDerivative::oneSided(),
                 ONE_BASIS_POINT);
}

CreditIndexBasisRhoPointwise::CreditIndexBasisRhoPointwise(double shiftSize):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>* newCreditIndexBasisRhoPointwise(double shiftSize) {
    return new CreditIndexBasisRhoPointwise(shiftSize);
}

CClassConstSP const CreditIndexBasisRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "CreditIndexBasisRhoPointwise", typeid(CreditIndexBasisRhoPointwise), load);

bool CreditIndexBasisRhoPointwiseLinkIn() {
    return CreditIndexBasisRhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
