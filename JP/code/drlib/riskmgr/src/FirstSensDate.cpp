//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : FirstSensDate.cpp
//
//   Description : Interface for instruments that know when to begin
//                 tweaking pointwise greeks. Modeled after LastSensDate.cpp
//
//   Author      : Sean Chen
//
//   Date        : 5 Aug 2005
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/FirstSensDate.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

FirstSensDate::FirstSensDate() {}

FirstSensDate::~FirstSensDate() {}

CClassConstSP const FirstSensDate::TYPE = CClass::registerInterfaceLoadMethod(
    "FirstSensDate", typeid(FirstSensDate), 0);

DRLIB_END_NAMESPACE
