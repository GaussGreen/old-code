//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolatilityBS.cpp
//
//   Description : Interface for those vol's that return a CVolProcessedBS
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/Class.hpp"
DRLIB_BEGIN_NAMESPACE

CClassConstSP const IVolatilityBS::TYPE = CClass::registerInterfaceLoadMethod(
    "IVolatilityBS", typeid(IVolatilityBS), 0);

DRLIB_END_NAMESPACE
