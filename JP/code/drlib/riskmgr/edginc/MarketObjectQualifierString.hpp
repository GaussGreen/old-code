//------------------------------------------------------------------------------
//
//   Group       : QR - Core Analytics 
//
//   Description : A simple market object qualifier just comprising a string
//
//   Author      : Andrew Greene 
//
//   Date        : 14 September 2006
//
//------------------------------------------------------------------------------

#ifndef MOQ_MARKET_OBJECT_QUALIFIER_STRING_HPP
#define MOQ_MARKET_OBJECT_QUALIFIER_STRING_HPP

#include "edginc/IMarketObjectQualifier.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL MarketObjectQualifierString : public CObject,
                                                public virtual IMarketObjectQualifier
{
public:
    static CClassConstSP const TYPE;

    MarketObjectQualifierString();
    MarketObjectQualifierString(string qualifierString);
    bool equals(const IMarketObjectQualifier* imdq);
    static IObject* createFromString(CClassConstSP requiredType, const string& data);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultMarketObjectQualifierString();

    string qualifierString;
};

DECLARE(MarketObjectQualifierString);

DRLIB_END_NAMESPACE

#endif // MOQ_MARKET_OBJECT_QUALIFIER_STRING_HPP
