//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketFactor.cpp
//
//   Description : MarketFactor interface + other registration code
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/YieldAsset.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE
//// SenSensitiveStrikeDescriptor
SensitiveStrikeDescriptor::SensitiveStrikeDescriptor(): forwardOnly(false) {}

////// IMarketFactor /////////
IMarketFactor::IMarketFactor(){}
IMarketFactor::~IMarketFactor(){}
void IMarketFactor::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IMarketFactor, clazz);
    EXTENDS(IGetMarket);
}

CClassConstSP const IMarketFactor::TYPE = CClass::registerInterfaceLoadMethod(
    "IMarketFactor", typeid(IMarketFactor), load);

DEFINE_TEMPLATE_TYPE(IMarketFactorArray);

/** If supplied market factor wrapper is using the market data cache, then
    retrieves, and makes if necessary, a market factor composed of the
    underlying asset together with the requested currency treatment.
    The different ccyTreatments are as per CAsset. Currently only
    CAssets can have anything other 'N' currency treatment */
void IMarketFactor::getMarketData(const IModel*                 model, 
                                  const MarketData*             market,
                                  const string&                 ccyTreatment,
                                  const string&                 payOutYCName,
                                  MarketWrapper<IMarketFactor>& factor){
    static const string method("IMarketFactor::getMarketData");
    bool checkCcyTreatment = false;
    try{
        const string& name = factor.getName();
        if (factor.usingCache()){
            // find out the type of the object in the cache
            CClassConstSP factorType = market->getDataType(name, TYPE);
            if (CAsset::TYPE->isAssignableFrom(factorType)){
                // delegate to CAsset
                factor.setObject(CAssetSP(CAsset::makeAssetUsingCache(
                                              model, market, ccyTreatment,
                                              payOutYCName, name)));
            } else {
                // retrieve object - check ccy treatment
                MarketObjectSP mo(market->GetData(model, name, TYPE));
                factor.setObject(mo);
                checkCcyTreatment = true;
            }
        } else {
            if (CAsset::TYPE->isInstance(factor.get())){
                // defer to Asset
                CAssetWrapper assetWrapper(
                    CAssetSP::dynamicCast(factor.getSP()));
                CAsset::getAssetMarketData(model, market, ccyTreatment, 
                                           payOutYCName, assetWrapper);
                factor.setObject(assetWrapper.getSP());
            } else {
                // check ccy treatment
                checkCcyTreatment = true;
                // ensure factor has all its market data
                factor->getMarket(model, market);
            }
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
    if (checkCcyTreatment && 
        ccyTreatment != CAsset::CCY_TREATMENT_NONE &&
        ccyTreatment != CAsset::CCY_TREATMENT_VANILLA){
        throw ModelException(method, "Only "
                             "currency treatment (N)one supported for factors "
                             " that are not assets");
    }
}
#if 0
class IMarketFactor::Action: virtual public ObjectIteration::IAction{
public:
    string           name;
    IMarketFactorSP  marketFactor;
    ObjectIteration* objectIteration; // reference to the iteration
    Action(const string& name): name(name){}

    /** method invoked by recurse routine to record name of object */
    bool invoke(const IObjectSP& obj){
        IMarketFactorSP mf(IMarketFactorSP::dynamicCast(obj));
        if (mf->getName() == name){
            marketFactor = mf; // save
            objectIteration->quitRecursion(true); // exit immediately
            return false; // redundant
        }
        return true;
    }
};

/** Given the name and type of a market factor find the first instance
    of it in the supplied 'container' */
IMarketFactorSP IMarketFactor::find(IObjectConstSP container,
                                    const string&  name,
                                    CClassConstSP  clazz,
                                    bool           failIfNotFound){
    static const string method("IMarketFactor::find");
    Action action(name);
    if (!IMarketFactor::TYPE->isAssignableFrom(clazz)){
        throw ModelException(method, "Type "+clazz->getName()+
                             " is not a market factor");
    }
    // create the ObjectIteration class
    ObjectIteration iteration(action, clazz);
    action.objectIteration = &iteration; // allow access from Action class
    // then recurse over all the relevant objects
    iteration.recurse(IObjectSP::constCast(container));
    IMarketFactorSP mf(action.marketFactor);
    if (failIfNotFound && !mf){
        throw ModelException(method, "No object of type "+clazz->getName()+
                             " with name "+name+" found");
    }
    return mf;
}
#endif
///// support for wrapper class /////

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(IMarketFactorWrapper);

/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<IMarketFactorWrapper>::toIObject(
    const IMarketFactorWrapper& value){
    return IObjectConstSP::attachToRef(&value);
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<IMarketFactorWrapper>::toIObject(
    IMarketFactorWrapper& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DateTime */
IMarketFactorWrapper arrayObjectCast<IMarketFactorWrapper>::fromIObject(
    IObjectSP& value){
    return (*DYNAMIC_CAST(IMarketFactorWrapper, value.get()));
}

DEFINE_TEMPLATE_TYPE(IMarketFactorWrapperArray);

////// IGeneralAsset /////////
IGeneralAsset::IGeneralAsset(){}
IGeneralAsset::~IGeneralAsset(){}
void IGeneralAsset::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IGeneralAsset, clazz);
    EXTENDS(IMarketFactor);
}

CClassConstSP const IGeneralAsset::TYPE = CClass::registerInterfaceLoadMethod(
    "IGeneralAsset", typeid(IGeneralAsset), load);

////// IPriceAsset /////////
IPriceAsset::IPriceAsset(){}
IPriceAsset::~IPriceAsset(){}
void IPriceAsset::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IPriceAsset, clazz);
    EXTENDS(IGeneralAsset);
}

CClassConstSP const IPriceAsset::TYPE = CClass::registerInterfaceLoadMethod(
    "IPriceAsset", typeid(IPriceAsset), load);

////// IYieldAsset /////////
IYieldAsset::~IYieldAsset(){}
void IYieldAsset::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IYieldAsset, clazz);
    EXTENDS(IGeneralAsset);
}

CClassConstSP const IYieldAsset::TYPE = CClass::registerInterfaceLoadMethod(
    "IYieldAsset", typeid(IYieldAsset), load);

DRLIB_END_NAMESPACE
