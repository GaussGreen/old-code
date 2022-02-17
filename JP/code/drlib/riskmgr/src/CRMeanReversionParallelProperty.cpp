/**
 * @file CRMeanReversionParallelProperty.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CRMeanReversionParallelProperty.hpp"

DRLIB_BEGIN_NAMESPACE

CRMeanReversionParallelProperty::CRMeanReversionParallelProperty(): CObject(TYPE) {}
CRMeanReversionParallelProperty::~CRMeanReversionParallelProperty() {}

static void CRMeanReversionParallelProperty_load(CClassSP& clazz) {
    REGISTER(CRMeanReversionParallelProperty, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CRMeanReversionParallelProperty>::iObject);
}

CClassConstSP const CRMeanReversionParallelProperty::TYPE = CClass::registerClassLoadMethod("CRMeanReversionParallelProperty", typeid(CRMeanReversionParallelProperty), CRMeanReversionParallelProperty_load);

RiskProperty_TYPES(CRMeanReversionParallelProperty)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
