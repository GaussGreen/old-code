//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/NextStrike.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const INextStrike::TYPE = 
CClass::registerInterfaceLoadMethod("INextStrike", typeid(INextStrike), 0);

DRLIB_END_NAMESPACE
