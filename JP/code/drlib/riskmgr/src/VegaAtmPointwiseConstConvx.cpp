/**
 * @file VegaAtmPointwiseConstConvx.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VolAtmPointwiseConstConvx.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/VegaAtmPointwiseConstConvx.hpp"

DRLIB_BEGIN_NAMESPACE

const string VegaAtmPointwiseConstConvx::NAME = "VEGA_ATM_POINTWISE_CONST_CONVX";
const double VegaAtmPointwiseConstConvx::DEFAULT_SHIFT = 0.001;
const double SENSITIVITY_UNIT = 0.01;

VegaAtmPointwiseConstConvx::VegaAtmPointwiseConstConvx(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, name)
{}

VegaAtmPointwiseConstConvx::VegaAtmPointwiseConstConvx(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv VegaAtmPointwiseConstConvx::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<VolAtmPointwiseConstConvx>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

VegaAtmPointwiseConstConvx::~VegaAtmPointwiseConstConvx() {}

void VegaAtmPointwiseConstConvx::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VegaAtmPointwiseConstConvx, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VegaAtmPointwiseConstConvx>::iObject);
    SensitivityFactory::addSens(VegaAtmPointwiseConstConvx::NAME, 
                                new GenericSensitivityFactory<VegaAtmPointwiseConstConvx>(), 
                                new VegaAtmPointwiseConstConvx(VegaAtmPointwiseConstConvx::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<VolAtmPointwiseConstConvx>::TYPE);
}

CClassConstSP const VegaAtmPointwiseConstConvx::TYPE = CClass::registerClassLoadMethod(
    "VegaAtmPointwiseConstConvx", typeid(VegaAtmPointwiseConstConvx), load);

bool VegaAtmPointwiseConstConvxLinkIn() {
    return VegaAtmPointwiseConstConvx::TYPE != NULL;
}

DRLIB_END_NAMESPACE
