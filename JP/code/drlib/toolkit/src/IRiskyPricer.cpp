//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : IRiskyPricer.cpp
//
//   Description : Interface that must be implemented if an Instrument needs to use
//                 a 'risky price' (i.e., a price assuming that equity is subject to 
//                 the risk of default) for a given sensitivity.
//                 Adapted from IRiskyPricer.cpp
//
//   Author      : Milan Kovacevic
//
//   Date        : 19th November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRiskyPricer.hpp"

DRLIB_BEGIN_NAMESPACE

IRiskyPricer::IRiskyPricer() {}

IRiskyPricer::~IRiskyPricer() {}

CClassConstSP const IRiskyPricer::TYPE = CClass::registerInterfaceLoadMethod(
    "IRiskyPricer", typeid(IRiskyPricer), 0);


DRLIB_END_NAMESPACE






