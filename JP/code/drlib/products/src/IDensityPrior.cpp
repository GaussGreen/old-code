//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : IDensityPrior.cpp
//
//   Description : Base class for probability density priors
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IDensityPrior.hpp"

DRLIB_BEGIN_NAMESPACE


IDensityPrior::~IDensityPrior()
{}


void IDensityPrior::load (CClassSP& clazz) {
    REGISTER_INTERFACE(IDensityPrior, clazz);
    EXTENDS(IObject);
}


CClassConstSP const IDensityPrior::TYPE = 
    CClass::registerInterfaceLoadMethod("IDensityPrior", 
                                        typeid(IDensityPrior), 
                                        load);

DRLIB_END_NAMESPACE

