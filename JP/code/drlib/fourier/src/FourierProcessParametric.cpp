/**
 * @file FourierProcessParametric.cpp
 */

#include "edginc/config.hpp"
#include "edginc/FourierProcessParametric.hpp"

DRLIB_BEGIN_NAMESPACE

FourierProcessParametric::FourierProcessParametric(CClassConstSP type):
    FourierProcess(type),
    allowRiskMapping(false)
{}

FourierProcessParametric::~FourierProcessParametric() {}

IModel::WantsRiskMapping FourierProcessParametric::wantsRiskMapping() const {
    return allowRiskMapping ? IModel::riskMappingAllowed :
                              IModel::riskMappingDisallowed;
}

void FourierProcessParametric::load(CClassSP& clazz) {
    REGISTER(FourierProcessParametric, clazz);
    SUPERCLASS(FourierProcess);
    FIELD(allowRiskMapping,
                 "Enable Black-Scholes greeks if suitable risk mapping matrix "
                   "is available in the market data");
    FIELD_MAKE_OPTIONAL(allowRiskMapping);
}

CClassConstSP const FourierProcessParametric::TYPE = CClass::registerClassLoadMethod(
    "FourierProcessParametric", typeid(FourierProcessParametric), load);

DRLIB_END_NAMESPACE
