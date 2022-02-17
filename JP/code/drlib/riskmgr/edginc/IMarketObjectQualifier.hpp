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

#ifndef IMARKET_OBJECT_QUALIFIER_HPP
#define IMARKET_OBJECT_QUALIFIER_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL IMarketObjectQualifier : virtual public IObject
{
public:
    static CClassConstSP const TYPE;

    IMarketObjectQualifier();
    virtual ~IMarketObjectQualifier();

    virtual bool equals(const IMarketObjectQualifier* imdq) = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<IMarketObjectQualifier> IMarketObjectQualifierConstSP;
typedef smartPtr<IMarketObjectQualifier> IMarketObjectQualifierSP;

DRLIB_END_NAMESPACE

#endif // IMARKET_OBJECT_QUALIFIER_HPP
