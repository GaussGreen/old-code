/**
 * @file StdDevForLossMapping.cpp
	Author   : Juan Carlos Porras      Dec 06
	Definitions for std dev tweak to loss mapping
 */

#include "edginc/config.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/StdDevForLossMapping.hpp"
#include "edginc/DefaultConstructor.hpp"

DRLIB_BEGIN_NAMESPACE


StdDevForLossMapping::StdDevForLossMapping(): CObject(TYPE) {}
StdDevForLossMapping::~StdDevForLossMapping() {}

static void StdDevForLossMapping_load(CClassSP& clazz) {
    REGISTER(StdDevForLossMapping, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<StdDevForLossMapping>::iObject);
}

CClassConstSP const StdDevForLossMapping::TYPE = 
	CClass::registerClassLoadMethod(
		"StdDevForLossMapping", 
		typeid(StdDevForLossMapping), 
		StdDevForLossMapping_load);



RiskProperty_TYPES(StdDevForLossMapping)

DRLIB_END_NAMESPACE
