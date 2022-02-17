//----------------------------------------------------------------------------
//   Group       : QR
//
//   Filename    : PrepayParallel.cpp
//
//   Description : ABCDS prepay curve horizontal parallel shift
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/PrepayParallelTP.hpp"
#include "edginc/PrepayParallel.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string PrepayParallel::NAME = "PREPAY_PARALLEL";
const double THREE_MONTH = 0.25;
const double PrepayParallel::DEFAULT_SHIFT = THREE_MONTH;
const double SENSITIVITY_UNIT = 1;

PrepayParallel::PrepayParallel(double shiftSize,
                               const string& packetName):
    ScalarRiskPropertySensitivity(TYPE, packetName, shiftSize)
{}

PrepayParallel::PrepayParallel(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv PrepayParallel::deriv() const {
    double sensUnit = getShiftSize();
    return Deriv(IResultsFunction::price(),
                 RiskProperty<PrepayParallelTP>::SP(),
                 IScalarDerivative::oneSided(),
                 sensUnit);
}

PrepayParallel::~PrepayParallel() {}

void PrepayParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(PrepayParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<PrepayParallel>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<PrepayParallel>(), 
                                new PrepayParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<PrepayParallelTP>::TYPE);
}

CClassConstSP const PrepayParallel::TYPE = CClass::registerClassLoadMethod(
    "PrepayParallel", typeid(PrepayParallel), load);

bool PrepayParallelLinkIn() {
    return PrepayParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

