//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IXLInterfaceMap.cpp
//
//   Description : Allow objects to have a different Interface 
//                 on the spreadsheet to the internal class model
//
//   Author      : Stephen Hope
//
//   Date        : 5th March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IXLInterfaceMap.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

IXLInterfaceMap::IXLInterfaceMap() {}

IXLInterfaceMap::~IXLInterfaceMap() {}

CClassConstSP const IXLInterfaceMap::TYPE = CClass::registerInterfaceLoadMethod(
    "IXLInterfaceMap", typeid(IXLInterfaceMap), 0);


DRLIB_END_NAMESPACE
