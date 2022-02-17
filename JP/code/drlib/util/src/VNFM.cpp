//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VNFM.cpp
//
//   Description : An interface for VNFM approximation
//
//   Author      : Keith Jia
//
//   Date        : August 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VNFM.hpp"

DRLIB_BEGIN_NAMESPACE

void VNFM::load(CClassSP& clazz) {
    REGISTER_INTERFACE(VNFM, clazz);
    EXTENDS(IObject);
}

VNFM::~VNFM() {}

VNFM::VNFM(CClassConstSP clazz)
{}

   
CClassConstSP const VNFM::TYPE = CClass::registerInterfaceLoadMethod(
    "VNFM", typeid(VNFM), load);


DRLIB_END_NAMESPACE
