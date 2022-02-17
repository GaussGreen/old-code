//----------------------------------------------------------------------------
//
//   Group       : EDR
//
//   Filename    : CreditDefaultSwap.cpp
//
//   Description : Credit default swap - 'flow' version (rather than 'generic')
//
//   Author      : Mark A Robson
//
//   Date        : January 12, 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CredDefSwap.hpp"
#include "edginc/ICreditEventOverrideName.hpp"

DRLIB_BEGIN_NAMESPACE

/** Flow version of a Credit Default Swap. Done as a separate type in case
    we want to completely separate it from the generic version. For
    speedy development just is a 'copy' of the generic CredDefSwap */
class CreditDefaultSwap: public CredDefSwap {
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        try{
            // we really want the dcc field in CredDefSwap to be mandatory
            // but for backward compatibility do it here instead
            validateDCCSuppliedAndNoBDC();
            CredDefSwap::validatePop2Object();
        } catch (exception& e){
            throw ModelException(e, "CreditDefaultSwap::validatePop2Object");
        }
    }
    
    /** instrument specific market data handling */
    virtual void getMarket(const IModel* model, const MarketData* market) {
        CredDefSwap::getMarket(model,market);
    }   
    
private:
    CreditDefaultSwap(): CredDefSwap(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditDefaultSwap, clazz);
        SUPERCLASS(CredDefSwap);
        EMPTY_SHELL_METHOD(defaultConstructor);
    }

    static IObject* defaultConstructor(){
        return new CreditDefaultSwap();
    }
};

CClassConstSP const CreditDefaultSwap::TYPE = CClass::registerClassLoadMethod(
    "CreditDefaultSwap", typeid(CreditDefaultSwap), load);

/* to ensure class is linked in */
bool CreditDefaultSwapLoad(){
    return true;
}

DRLIB_END_NAMESPACE

