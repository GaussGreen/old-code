/**
 * @file MCPathConfigParametric.cpp
 */

#include "edginc/config.hpp"
#include "edginc/MCPathConfigParametric.hpp"
#include "edginc/Dependence.hpp"

DRLIB_BEGIN_NAMESPACE

MCPathConfigParametric::MCPathConfigParametric(CClassConstSP type):
    MCPathConfig(type),
    allowRiskMapping(false)
{}

MCPathConfigParametric::MCPathConfigParametric(CClassConstSP type,
                                               DependenceMakerSP dependenceMaker):
    MCPathConfig(type, dependenceMaker),
    allowRiskMapping(false)
{}

MCPathConfigParametric::~MCPathConfigParametric() {}

IModel::WantsRiskMapping MCPathConfigParametric::wantsRiskMapping() const {
    return allowRiskMapping ? IModel::riskMappingAllowed :
                              IModel::riskMappingDisallowed;
}

void MCPathConfigParametric::load(CClassSP& clazz) {
    REGISTER(MCPathConfigParametric, clazz);
    SUPERCLASS(MCPathConfig);
    FIELD(allowRiskMapping,
                 "Enable Black-Scholes greeks if suitable risk mapping matrix "
                   "is available in the market data");
    FIELD_MAKE_OPTIONAL(allowRiskMapping);
}

CClassConstSP const MCPathConfigParametric::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigParametric", typeid(MCPathConfigParametric), load);

DRLIB_END_NAMESPACE
