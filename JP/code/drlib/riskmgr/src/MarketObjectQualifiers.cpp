//------------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : A market object qualifier comprising string array matching each element
//
//   Date        : 14 September 2006
//
//------------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/Addin.hpp"
#include "edginc/MarketObjectQualifiers.hpp"

DRLIB_BEGIN_NAMESPACE

MarketObjectQualifiers::MarketObjectQualifiers():  CObject(TYPE){}

MarketObjectQualifiers::MarketObjectQualifiers(StringArrayConstSP qualifiers):
    CObject(TYPE), qualifiers(qualifiers){}

// returns true only if each element of string array matches
bool MarketObjectQualifiers::equals(const IMarketObjectQualifier* imdq)
{
    const MarketObjectQualifiers* moqs;

    try {
        moqs = DYNAMIC_CONST_CAST(MarketObjectQualifiers, imdq);
    }
    catch (ModelException) {
        return false;
    }

    if (qualifiers->size() != moqs->qualifiers->size())
        return false;

    bool isEqual = true;
    for (int i=0; i<qualifiers->size(); ++i){
        isEqual &= ((*qualifiers)[i] == (*moqs->qualifiers)[i]);
        if (!isEqual)
            return false;
    }

    return true;
}

CClassConstSP const MarketObjectQualifiers::TYPE =
    CClass::registerClassLoadMethod("MarketObjectQualifiers",
                                    typeid(MarketObjectQualifiers),
                                    load);

void MarketObjectQualifiers::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MarketObjectQualifiers, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMarketObjectQualifier);
    EMPTY_SHELL_METHOD(defaultMarketObjectQualifiers);
    FIELD(qualifiers, "Qualifier string array");
    Addin::registerConstructor(Addin::UTILITIES, MarketObjectQualifiers::TYPE);
}

IObject* MarketObjectQualifiers::defaultMarketObjectQualifiers()
{
    return new MarketObjectQualifiers;
}

bool MarketObjectQualifiersLinkIn()
{
    return MarketObjectQualifiers::TYPE != NULL;
}

DRLIB_END_NAMESPACE
