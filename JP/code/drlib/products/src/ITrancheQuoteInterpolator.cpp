//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ITrancheQuoteInterpolator.cpp
//
//   Description : Interface for classes taht can interpolate tranche quotes and produce an ExpectedLossSurface
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ITrancheQuoteInterpolator.hpp"

DRLIB_BEGIN_NAMESPACE


ITrancheQuoteInterpolator::~ITrancheQuoteInterpolator()
{}


void ITrancheQuoteInterpolator::load (CClassSP& clazz) {
    REGISTER_INTERFACE(ITrancheQuoteInterpolator, clazz);
    EXTENDS(IObject);
}


CClassConstSP const ITrancheQuoteInterpolator::TYPE = 
    CClass::registerInterfaceLoadMethod("ITrancheQuoteInterpolator", 
                                        typeid(ITrancheQuoteInterpolator), 
                                        load);

DRLIB_END_NAMESPACE
