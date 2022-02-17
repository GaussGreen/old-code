//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Additive.cpp
//
//   Description : Marker interface for sensitivities to show if they can
//                 be added together
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

Additive::Additive() {}

Additive::~Additive() {}

CClassConstSP const Additive::TYPE = CClass::registerInterfaceLoadMethod(
    "Additive", typeid(Additive), 0);

DRLIB_END_NAMESPACE
