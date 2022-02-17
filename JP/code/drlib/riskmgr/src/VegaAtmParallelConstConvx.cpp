/**
 * @file VegaAtmParallelConstConvx.cpp
 */
#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VolAtmParallelConstConvx.hpp"
#include "edginc/VegaAtmParallelConstConvx.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"



DRLIB_BEGIN_NAMESPACE

const string VegaAtmParallelConstConvx::NAME = "VEGA_ATM_PARALLEL_CONST_CONVX";
const double VegaAtmParallelConstConvx::DEFAULT_SHIFT = 0.001;
const double VegaAtmParallelConstConvx::SENSITIVITY_UNIT = 0.01;

VegaAtmParallelConstConvx::VegaAtmParallelConstConvx(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

ScalarRiskPropertySensitivity::Deriv VegaAtmParallelConstConvx::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<VolAtmParallelConstConvx>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

VegaAtmParallelConstConvx::~VegaAtmParallelConstConvx() {}

void VegaAtmParallelConstConvx::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VegaAtmParallelConstConvx, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VegaAtmParallelConstConvx>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<VegaAtmParallelConstConvx>(), 
                                new VegaAtmParallelConstConvx(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<VolAtmParallelConstConvx>::TYPE);
}

CClassConstSP const VegaAtmParallelConstConvx::TYPE = CClass::registerClassLoadMethod(
    "VegaAtmParallelConstConvx", typeid(VegaAtmParallelConstConvx), load);

bool VegaAtmParallelConstConvxLinkIn() {
    return VegaAtmParallelConstConvx::TYPE != NULL;
}

DRLIB_END_NAMESPACE
