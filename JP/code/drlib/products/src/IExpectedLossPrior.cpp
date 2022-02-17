//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : IExpectedLossPrior.cpp
//
//   Description : Base class for expected loss priors
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IExpectedLossPrior.hpp"

DRLIB_BEGIN_NAMESPACE


IExpectedLossPrior::~IExpectedLossPrior()
{}


void IExpectedLossPrior::load (CClassSP& clazz) {
    REGISTER_INTERFACE(IExpectedLossPrior, clazz);
    EXTENDS(IObject);
}


CClassConstSP const IExpectedLossPrior::TYPE = 
    CClass::registerInterfaceLoadMethod("IExpectedLossPrior", 
                                        typeid(IExpectedLossPrior), 
                                        load);

DRLIB_END_NAMESPACE
