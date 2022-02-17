//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetCcyCollector.cpp
//
//   Description : asset currency collector class
//
//   Author      : André Segger
//
//   Date        : 13 Jul 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/Asset.hpp"
#include "edginc/ObjectIteration.hpp"

DRLIB_BEGIN_NAMESPACE
/** Invoked when TestCollect is 'loaded' */
void AssetCcyCollector::load(CClassSP& clazz){
    REGISTER(AssetCcyCollector, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICollector);
}

CClassConstSP const AssetCcyCollector::TYPE = CClass::registerClassLoadMethod(
    "AssetCcyCollector", typeid(AssetCcyCollector), load);

AssetCcyCollector::AssetCcyCollector(const string& imntIsoCode) : 
        CObject(TYPE), imntIsoCode(imntIsoCode), empty(false) {}

AssetCcyCollector::AssetCcyCollector() : 
        CObject(TYPE), empty(true) {}

void AssetCcyCollector::currencyValidate(const string& isoCode,
                                         const string& assetName){
    if (!empty) {
        if (isoCode != imntIsoCode) {
            throw ModelException(
                "AssetCcyCollector::currencyValidate",
                "The currency of asset " + assetName + 
                " (" + isoCode + ") is different from \n"
                "the payoff currency (" + 
                imntIsoCode + ")");
        }
    } else {
        // stash the incoming ccy
        empty = false;
        imntIsoCode = isoCode;
    }
}

void AssetCcyCollector::validateAllCurrencies(IObjectConstSP       obj,
                                              const YieldCurve*    discCcy){
    if (discCcy) {
        validateAllCurrencies(obj, discCcy->getCcy());
    }
}

void AssetCcyCollector::validateAllCurrencies(IObjectConstSP obj, 
                                              const string&  discCcy)
{
    // class for handling call back
    class Action: public ObjectIteration::IActionConst{
    public:
        AssetCcyCollector* collector;
        // called by ObjectIteration
        bool invoke(const ObjectIteration::State& state, IObjectConstSP obj) {
            obj->accept(collector);
            return false;
        }
    };
        
    AssetCcyCollector collector(discCcy);
    // create the action object
    Action action;
    action.collector = &collector;
    // build an instance of the class that drives the recursion
    ObjectIteration iteration(Asset::TYPE);
    // go
    iteration.recurse(action, obj);
}

const string& AssetCcyCollector::ccyCode() const {
    if (empty) {
        throw ModelException("AssetCcyCollector::ccyCode",
                             "no ccy has been found");
    }
    return imntIsoCode;
}

DRLIB_END_NAMESPACE

