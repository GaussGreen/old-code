//----------------------------------------------------------------------------
//
//   Filename    : ObservationOverride.cpp
//
//   Description : Interface for classes which provide overrides for
//                 observations for IMarketObservables 
//                 e.g. PastValues provides such overrides
//
//   Author      : Ian Stares
//
//   Date        : July 4 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ObservationOverride.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IObservationOverride::TYPE = CClass::registerInterfaceLoadMethod(
    "IObservationOverride", typeid(IObservationOverride), load);

void IObservationOverride::load(CClassSP& clazz){
    REGISTER_INTERFACE(IObservationOverride, clazz);
    EXTENDS(IObject);
}

DRLIB_END_NAMESPACE




