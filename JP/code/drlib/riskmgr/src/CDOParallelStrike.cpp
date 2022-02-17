//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CDOParallelStrike.cpp
//
//   Description : Tag class for parallel shift in CDO strikes 
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CDOParallelStrike.hpp"

DRLIB_BEGIN_NAMESPACE

CDOParallelStrike::CDOParallelStrike(): CObject(TYPE) {}
CDOParallelStrike::~CDOParallelStrike() {}

static void CDOParallelStrike_load(CClassSP& clazz) {
    REGISTER(CDOParallelStrike, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CDOParallelStrike>::iObject);
}

CClassConstSP const CDOParallelStrike::TYPE = CClass::registerClassLoadMethod("CDOParallelStrike", typeid(CDOParallelStrike), CDOParallelStrike_load);

RiskProperty_TYPES(CDOParallelStrike)

DRLIB_END_NAMESPACE
