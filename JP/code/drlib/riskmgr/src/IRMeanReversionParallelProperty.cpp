/**
 * @file IRMeanReversionParallelProperty.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/IRMeanReversionParallelProperty.hpp"

DRLIB_BEGIN_NAMESPACE

IRMeanReversionParallelProperty::IRMeanReversionParallelProperty(): CObject(TYPE) {}
IRMeanReversionParallelProperty::~IRMeanReversionParallelProperty() {}

static void IRMeanReversionParallelProperty_load(CClassSP& clazz) {
    REGISTER(IRMeanReversionParallelProperty, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<IRMeanReversionParallelProperty>::iObject);
}

CClassConstSP const IRMeanReversionParallelProperty::TYPE = CClass::registerClassLoadMethod("IRMeanReversionParallelProperty", typeid(IRMeanReversionParallelProperty), IRMeanReversionParallelProperty_load);

RiskProperty_TYPES(IRMeanReversionParallelProperty)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
