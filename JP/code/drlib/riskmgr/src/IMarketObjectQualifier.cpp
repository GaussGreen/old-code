//------------------------------------------------------------------------------
//
//   Group       : QR - Core Analytics 
//
//   Description : General qualifier for distinguishing multiple instances of
//                 market data with the same (root) name and type
//
//   Author      : Andrew Greene 
//
//   Date        : 31 August 2006
//
//------------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/IMarketObjectQualifier.hpp"

#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

IMarketObjectQualifier::IMarketObjectQualifier()
{
}

IMarketObjectQualifier::~IMarketObjectQualifier()
{
}

void IMarketObjectQualifier::load(CClassSP& clazz)
{
    REGISTER_INTERFACE(IMarketObjectQualifier, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IMarketObjectQualifier::TYPE =
    CClass::registerInterfaceLoadMethod("IMarketObjectQualifier",
                                        typeid(IMarketObjectQualifier),
                                        load);

DRLIB_END_NAMESPACE
