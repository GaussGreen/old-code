//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/ScaleOutputs.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IScaleOutputs::TYPE = 
CClass::registerInterfaceLoadMethod("IScaleOutputs", typeid(IScaleOutputs), 0);

DRLIB_END_NAMESPACE
