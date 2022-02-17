/**
 * @file CRVegaParallel.cpp
 * 
 * Group : Credit QR NY
 * 
 * File  : CRVegaParallel.cpp
 * 
 * Description : the tweak tag that describes a parallel tweak for spread volatility cubes. The greek is calculated by tweaking all market volatility smiles at all its
 * defined expiries simultaneously.
 * 
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
#include "edginc/CRVolParallel.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/CRVegaParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CRVegaParallel::NAME = "CR_VEGA_PARALLEL";
const double CRVegaParallel::DEFAULT_SHIFT = 0.001;
const double CRVegaParallel::SENSITIVITY_UNIT = 0.01;

CRVegaParallel::CRVegaParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv CRVegaParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CRVolParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CRVegaParallel::~CRVegaParallel() {}

void CRVegaParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CRVegaParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CRVegaParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CRVegaParallel>(), 
                                new CRVegaParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CRVolParallel>::TYPE);
}

CClassConstSP const CRVegaParallel::TYPE = CClass::registerClassLoadMethod(
    "CRVegaParallel", typeid(CRVegaParallel), load);

bool CRVegaParallelLinkIn() {
    return CRVegaParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

