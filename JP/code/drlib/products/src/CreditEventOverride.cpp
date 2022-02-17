//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CreditEventOverride.cpp
//
//   Description : Simple class used by credit instruments (e.g., CIS) in order 
//                 to override (at instrument level) a name's default-related 
//                 parameters.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CreditEventOverride.hpp"
#include "edginc/CreditEventOverrideName.hpp"

DRLIB_BEGIN_NAMESPACE


CreditEventOverride::CreditEventOverride() : CObject(TYPE)
{}


CreditEventOverride::~CreditEventOverride()
{}


/** Called immediately after object constructed */
void CreditEventOverride::validatePop2Object() {
    static const string method = "CreditEventOverride::validatePop2Object";
    
    // Verify that there is no more than one entry in the array for each name
    int numOverrides = overrideArray->size();
    for (int i=0; i<numOverrides; ++i) {
        const string name = (*overrideArray)[i]->getName();

        for (int j=i+1; j<numOverrides; ++j) {
            if ((*overrideArray)[j]->getName() == name) {
                throw ModelException(method,
                                     "There are at least two overrides for "
                                     "name " + name);
            }
        }
    }
}

/** GetMarket implementation*/
void CreditEventOverride::getMarket(const IModel* model, 
                                    const MarketData* market)
{
    for (int i=0; i < overrideArray->size(); ++i) {
        (*overrideArray)[i]->getMarket(model, market);
    }
}


/** Returns the override information for the supplied name. If no override
 * information is available returns ICreditEventOverrideNameSP() */
ICreditEventOverrideNameSP CreditEventOverride::getOverrideForName(
   const string& name) const 
{
    for (int i=0; i < overrideArray->size(); ++i) {
        if ((*overrideArray)[i]->getName() == name) {
            return (*overrideArray)[i];
        }
    }
    return ICreditEventOverrideNameSP();
}


void CreditEventOverride::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditEventOverride, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICreditEventOverride);
    EMPTY_SHELL_METHOD(defaultCreditEventOverride);

    FIELD(overrideArray, "Array of ICreditEventOverrideNames");
}


IObject* CreditEventOverride::defaultCreditEventOverride() {
    return new CreditEventOverride();
}

CClassConstSP const CreditEventOverride::TYPE = 
    CClass::registerClassLoadMethod("CreditEventOverride", 
                                    typeid(CreditEventOverride), 
                                    load);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool CreditEventOverrideLoad() {
    return (CreditEventOverride::TYPE != 0);
}


DRLIB_END_NAMESPACE
