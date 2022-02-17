//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : E2CModel.cpp
//
//   Description : Interface that must be implemented if an Instrument needs to use
//
//   Author      : Andre Segger
//
//   Date        : 07 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/E2CModel.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

IE2CModel::IE2CModel() {}

CClassConstSP const IE2CModel::TYPE = CClass::registerInterfaceLoadMethod(
    "IE2CModel", typeid(IE2CModel), 0);


DRLIB_END_NAMESPACE






