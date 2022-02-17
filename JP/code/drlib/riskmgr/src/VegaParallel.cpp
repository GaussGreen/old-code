/**
 * @file VegaParallel.cpp
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
#include "edginc/VolParallel.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string VegaParallel::NAME = "VEGA_PARALLEL";
const double VegaParallel::DEFAULT_SHIFT = 0.001;
const double VegaParallel::SENSITIVITY_UNIT = 0.01;

VegaParallel::VegaParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv VegaParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<VolParallel>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

VegaParallel::~VegaParallel() {}

void VegaParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VegaParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VegaParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<VegaParallel>(), 
                                new VegaParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<VolParallel>::TYPE);
}

CClassConstSP const VegaParallel::TYPE = CClass::registerClassLoadMethod(
    "VegaParallel", typeid(VegaParallel), load);

bool VegaParallelLinkIn() {
    return VegaParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

