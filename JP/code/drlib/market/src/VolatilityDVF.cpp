//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolatilityDVF.cpp
//
//   Description : Interface for those vol's that return a CVolProcessedDVF
//
//   Author      : JNJ
//
//   Date        : 06 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Class.hpp"
DRLIB_BEGIN_NAMESPACE

CClassConstSP const IVolatilityDVF::TYPE = CClass::registerInterfaceLoadMethod(
    "IVolatilityDVF", typeid(IVolatilityDVF), 0);

DRLIB_END_NAMESPACE
