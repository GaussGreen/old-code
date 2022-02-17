//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : IExpectedLossInterpolator.cpp
//
//   Description : Interface for expected loss interpolators
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/IExpectedLossInterpolator.hpp"

DRLIB_BEGIN_NAMESPACE


IExpectedLossInterpolator::~IExpectedLossInterpolator()
{}


void IExpectedLossInterpolator::load (CClassSP& clazz) {
    REGISTER_INTERFACE(IExpectedLossInterpolator, clazz);
    EXTENDS(IObject);
}


CClassConstSP const IExpectedLossInterpolator::TYPE = 
    CClass::registerInterfaceLoadMethod("IExpectedLossInterpolator", 
                                        typeid(IExpectedLossInterpolator), 
                                        load);

DRLIB_END_NAMESPACE
