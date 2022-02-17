//----------------------------------------------------------------------------
//   Group       : QR
//
//   Filename    : ParSpreadUpfrontParallel.cpp
//
//   Description : ParSpread curve Upfront parallel shift
//
//   Date        : Nov 2007
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/ParSpreadUpfrontParallelTP.hpp"
#include "edginc/ParSpreadUpfrontParallel.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string ParSpreadUpfrontParallel::NAME = "PAR_SPREAD_UPFRONT_PARALLEL";
const double ONE_CENT = 0.0001;
const double ParSpreadUpfrontParallel::DEFAULT_SHIFT = ONE_CENT;

ParSpreadUpfrontParallel::ParSpreadUpfrontParallel(double shiftSize,
                                                   const string& packetName):
    ScalarRiskPropertySensitivity(TYPE, packetName, shiftSize)
{}

ParSpreadUpfrontParallel::ParSpreadUpfrontParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv ParSpreadUpfrontParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<ParSpreadUpfrontParallelTP>::SP(),
                 IScalarDerivative::oneSided(), 
                 ONE_CENT);
}

ParSpreadUpfrontParallel::~ParSpreadUpfrontParallel() {}

void ParSpreadUpfrontParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(ParSpreadUpfrontParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<ParSpreadUpfrontParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<ParSpreadUpfrontParallel>(), 
                                new ParSpreadUpfrontParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<ParSpreadUpfrontParallelTP>::TYPE);
}

CClassConstSP const ParSpreadUpfrontParallel::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadUpfrontParallel", typeid(ParSpreadUpfrontParallel), load);

bool ParSpreadUpfrontParallelLinkIn() {
    return ParSpreadUpfrontParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
