//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CombinableResult.cpp
//
//   Description : Interface for results that can be added together
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CombinableMixedResult.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

CombinableResult::CombinableResult(){}
CombinableResult::~CombinableResult(){}

CClassConstSP const CombinableResult::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CombinableResult", typeid(CombinableResult), 0);

CombinableMixedResult::CombinableMixedResult() {}

static void loadCombinableMixedResult(CClassSP& clazz){
    REGISTER_INTERFACE(CombinableMixedResult, clazz);
    EXTENDS(CombinableResult);
}

CClassConstSP const CombinableMixedResult::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CombinableMixedResult", typeid(CombinableMixedResult), 
    loadCombinableMixedResult);

DRLIB_END_NAMESPACE
