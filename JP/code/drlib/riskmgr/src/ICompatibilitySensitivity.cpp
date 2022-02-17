/**
 * @file ICompatibilitySensitivity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/ICompatibilitySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

ICompatibilitySensitivity::ICompatibilitySensitivity() {}
ICompatibilitySensitivity::~ICompatibilitySensitivity() {}

CClassConstSP const ICompatibilitySensitivity::TYPE = CClass::registerInterfaceLoadMethod(
    "ICompatibilitySensitivity", typeid(ICompatibilitySensitivity), 0);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
