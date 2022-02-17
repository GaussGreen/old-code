//------------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : A market object qualifier comprising string array matching each element
//
//   Date        : 14 September 2006
//
//------------------------------------------------------------------------------

#ifndef MOQ_MARKET_OBJECT_QUALIFIERS_HPP
#define MOQ_MARKET_OBJECT_QUALIFIERS_HPP

#include "edginc/IMarketObjectQualifier.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL MarketObjectQualifiers : public CObject,
                                           public virtual IMarketObjectQualifier
{
public:
    static CClassConstSP const TYPE;

    MarketObjectQualifiers();
    MarketObjectQualifiers(CStringArrayConstSP qualifiers);

    // returns true only if each element of string array matches
    bool equals(const IMarketObjectQualifier* imdq);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultMarketObjectQualifiers();

    CStringArrayConstSP qualifiers;
};

DECLARE(MarketObjectQualifiers);

DRLIB_END_NAMESPACE

#endif // MOQ_MARKET_OBJECT_QUALIFIER_ARRAY_HPP
