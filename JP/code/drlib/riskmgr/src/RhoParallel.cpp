/**
 * @file RhoParallel.cpp
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
#include "edginc/RhoParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string RhoParallel::NAME = "RHO_PARALLEL";
const double RhoParallel::DEFAULT_SHIFT = 0.0001;
const double RhoParallel::SENSITIVITY_UNIT = 0.0001;

RhoParallel::RhoParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv RhoParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<RateParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

RhoParallel::~RhoParallel() {}

void RhoParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(RhoParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<RhoParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<RhoParallel>(), 
                                new RhoParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<RateParallel>::TYPE);
}

CClassConstSP const RhoParallel::TYPE = CClass::registerClassLoadMethod(
    "RhoParallel", typeid(RhoParallel), load);

bool RhoParallelLinkIn() {
    return RhoParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

