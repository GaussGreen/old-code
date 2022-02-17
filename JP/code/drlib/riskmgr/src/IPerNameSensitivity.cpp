/**
 * @file IPerNameSensitivity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IPerNameSensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

IPerNameSensitivity::IPerNameSensitivity()
{}

IPerNameSensitivity::~IPerNameSensitivity() {}

CClassConstSP const IPerNameSensitivity::TYPE = CClass::registerInterfaceLoadMethod(
    "IPerNameSensitivity", typeid(IPerNameSensitivity), 0);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
