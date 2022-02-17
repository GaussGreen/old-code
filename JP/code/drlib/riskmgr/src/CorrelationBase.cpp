//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationBase.cpp
//
//   Description : CorrelationBase  - allows models to manipulate correlations
//
//   Author      : Mark A Robson
//
//   Date        : 19 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CORRELATIONBASE_CPP
#include "edginc/CorrelationBase.hpp"

DRLIB_BEGIN_NAMESPACE
CorrelationBase::~CorrelationBase(){}

CorrelationBase::CorrelationBase(CClassConstSP clazz): MarketObject(clazz){}

void CorrelationBase::load(CClassSP& clazz){
    REGISTER(CorrelationBase, clazz);
    SUPERCLASS(MarketObject);
    clazz->enableCloneOptimisations(); // for derived types
}

CClassConstSP const CorrelationBase::TYPE = CClass::registerClassLoadMethod(
    "CorrelationBase", typeid(CorrelationBase), load);
    
DRLIB_END_NAMESPACE

