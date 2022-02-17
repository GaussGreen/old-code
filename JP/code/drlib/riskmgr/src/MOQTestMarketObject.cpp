//------------------------------------------------------------------------------
//
//   Group       : QR - Core Analytics 
//
//   Description : Market object for testing retrieval of market data using
//                 market data qualifiers
//
//   Author      : Andrew Greene 
//
//   Date        : 14 September 2006
//
//------------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/MOQTestMarketObject.hpp"

DRLIB_BEGIN_NAMESPACE

MOQTestMarketObject::MOQTestMarketObject():
    MarketObject(TYPE)
{}

/** Returns the name by which this piece of market data is referenced by
    other market data or by instruments/models etc. For example a
    collection of different type of parameterized equity volatilities for
    the FTSE would all have the same name eg “FTSE”. Similarly, base and
    swaption interest rate volatilities for the same yield curve would
    have the same root name. The default implementation is to return getName() */
string MOQTestMarketObject::getRootName() const
{
    return rootName;
}

/** When more than one market object of appropriate type has the same root
    name an additional qualifier is needed in order to determine which
    instance should be used. It is proposed that an object of type
    IMarketObjectQualifier is used to qualify MarketObjects which have the
    same "root name". This type would contain only one method, namely a
    pure virtual equals(IMarketObjectQualifier*) method. One obvious
    implementation of such an interface is a class that just wraps a single
    string. Following the examples in getRootName() this could be, for
    example, "FLOW" or "EXOTIC" for equity volatilies or "BASE" or
    "SWAPTION" for interest rate volatilities. (These strings are only
    examples - there are no predefined values.). The default implementation
    returns a NULL smart pointer. */
IMarketObjectQualifierSP MOQTestMarketObject::getMarketObjectQualifier() const
{
    return IMarketObjectQualifierSP(
        DYNAMIC_CAST(IMarketObjectQualifier, qualifier.clone()));
}

/** Returns the identifier for this object as supplied by the client. It
    is the name typically used to identify sensitivities. Note that it is
    not the name used when retrieving market data (see getRootName()). As
    such it is not a ‘magic’ string, i.e., its value does not have to
    appear anywhere else in the input data. */
string MOQTestMarketObject::getName() const
{
    return name;
}

CClassConstSP const MOQTestMarketObject::TYPE =
    CClass::registerClassLoadMethod("MOQTestMarketObject",
                                    typeid(MOQTestMarketObject),
                                    load);

void MOQTestMarketObject::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MOQTestMarketObject, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(MOQTestMarketObject::defaultMOQTestMarketObject);

    FIELD(rootName , "Root name"       );
    FIELD(name     , "Name"            );
    FIELD(qualifier, "String qualifier");
    FIELD(data     , "Integer data"    );
}

IObject* MOQTestMarketObject::defaultMOQTestMarketObject()
{
    return new MOQTestMarketObject;
}

bool MOQTestMarketObjectLoad()
{
    return MOQTestMarketObject::TYPE != NULL;
}

DRLIB_END_NAMESPACE
