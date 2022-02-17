/**
 * @file IDynamicsParameter.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IDynamicsParameter.hpp"

DRLIB_BEGIN_NAMESPACE

IDynamicsParameter::IDynamicsParameter()
{}

IDynamicsParameter::~IDynamicsParameter() {}

CClassConstSP const IDynamicsParameter::TYPE = CClass::registerInterfaceLoadMethod(
    "IDynamicsParameter", typeid(IDynamicsParameter), 0);

DRLIB_END_NAMESPACE
