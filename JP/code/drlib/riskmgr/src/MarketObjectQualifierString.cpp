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

#include "edginc/config.hpp"

#include "edginc/MarketObjectQualifierString.hpp"

#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

MarketObjectQualifierString::MarketObjectQualifierString():
    CObject(TYPE)
{}

MarketObjectQualifierString::MarketObjectQualifierString(string qualifierString):
    CObject(TYPE),
    qualifierString(qualifierString)
{}

bool MarketObjectQualifierString::equals(const IMarketObjectQualifier* imdq)
{
    const MarketObjectQualifierString* moqs;

    try {
        moqs = DYNAMIC_CONST_CAST(MarketObjectQualifierString, imdq);
    }
    catch (ModelException) {
        return false;
    }

    return moqs->qualifierString == qualifierString;
}

IObject* MarketObjectQualifierString::createFromString(CClassConstSP requiredType,
                                                       const string& data)
{
    if (!TYPE->isAssignableFrom(requiredType)) {
        throw ModelException(requiredType->getName() +
                             " not derived from MarketObjectQualifierString");
    }

    return new MarketObjectQualifierString(data);
}

CClassConstSP const MarketObjectQualifierString::TYPE =
    CClass::registerClassLoadMethod("MarketObjectQualifierString",
                                    typeid(MarketObjectQualifierString),
                                    load);

void MarketObjectQualifierString::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MarketObjectQualifierString, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMarketObjectQualifier);
    EMPTY_SHELL_METHOD(MarketObjectQualifierString::defaultMarketObjectQualifierString);
    FIELD(qualifierString, "Qualifier String");
    CString::registerObjectFromStringMethod(TYPE, createFromString);
}

IObject* MarketObjectQualifierString::defaultMarketObjectQualifierString()
{
    return new MarketObjectQualifierString;
}

DRLIB_END_NAMESPACE
