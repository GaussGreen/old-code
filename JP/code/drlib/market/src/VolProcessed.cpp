//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessed.cpp
//
//   Description : Abstract Processed Vol Interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLPROCESSED_CPP
#include "edginc/VolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

IVolProcessed::~IVolProcessed(){}

IVolProcessed::IVolProcessed() {}

static void myLoad(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IVolProcessed, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IVolProcessed::TYPE = CClass::registerClassLoadMethod(
    "VolProcessed", typeid(IVolProcessed), myLoad);

DRLIB_END_NAMESPACE
