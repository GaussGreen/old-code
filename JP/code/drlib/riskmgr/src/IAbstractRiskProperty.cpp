/**
 * @file IAbstractRiskProperty.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IAbstractRiskProperty.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/RiskMappingMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

IAbstractRiskProperty::IAbstractRiskProperty()
{}

IAbstractRiskProperty::~IAbstractRiskProperty() {}

RiskMappingMatrixConstSP IAbstractRiskProperty::riskMappingMatrix(
                                                    IObjectConstSP world,
                                                    OutputNameConstSP name) const {
    return RiskMappingMatrixSP();
}

CClassConstSP const IAbstractRiskProperty::TYPE = CClass::registerInterfaceLoadMethod(
    "IAbstractRiskProperty", typeid(IAbstractRiskProperty), 0);

DEFINE_TEMPLATE_TYPE(IAbstractRiskPropertyArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
