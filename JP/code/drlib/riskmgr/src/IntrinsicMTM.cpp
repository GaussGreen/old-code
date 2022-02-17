//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/IntrinsicMTM.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IIntrinsicMTM::TYPE = 
CClass::registerInterfaceLoadMethod("IIntrinsicMTM", typeid(IIntrinsicMTM), 0);

DRLIB_END_NAMESPACE
