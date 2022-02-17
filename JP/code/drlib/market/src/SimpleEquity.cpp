//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimpleEquity.hpp
//
//   Description : Single factor, no currency treatment, equity asset
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SimpleEquity.hpp"


DRLIB_BEGIN_NAMESPACE

SimpleEquity::~SimpleEquity(){}

/** returns the asset name */
string SimpleEquity::getName() const {
    // default to equity name unless we have our own
    if (name.empty()) {
        return equity->getName();
    }
    else {
        return name;
    }
}

/** Allows asset to be currency struck.
    Combines market and instrument data together to give a
    Processed Vol. Here the processed volatility is a processed
    struck volatility ie it reflects the combination of this
    asset together with the supplied FX asset and the
    correlation between this CVolBase and the vol of the
    FX. Note that the struckAsset is indeed the struckAsset cf
    'this' which is the non struck asset */
CVolProcessed* SimpleEquity::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      struckAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    return vol->getProcessedVol(volRequest, struckAsset,
                                fxAsset, eqFXCorr);
}

/** adds the credit spread to the asset's growth curve */
void SimpleEquity::makeRisky(ICreditCurveSP  creditSpreads,
                             const  DateTime *maturityDate)
{
       equity->makeRiskyEquity(creditSpreads,maturityDate);
}

// for given trade date, return settle date and name of the asset
void SimpleEquity::delivery(
    const DateTime&               tradeDate, 
    double                        quantity, 
    double                        price,
    PhysicalDeliveryByAssetArray* byAsset) const {
    try {
        PhysicalDeliverySP deliv(new PhysicalDelivery(quantity,
                                                      price,
                                                      tradeDate,
                                                      settleDate(tradeDate)));
        PhysicalDelivery::ByAsset delivByAsset(getName(), deliv.get()); 

        byAsset->push_back(delivByAsset);
    }
    catch (exception& e) {
        throw ModelException(e, "SimpleEquity::delivery");
    }
}

/** Populates the object with the market data that this object
needs.  In particular this will just call the equivalen function
in EquityBase and then retrieve the AssetHistory(s)*/
void SimpleEquity::getMarket(const IModel* model, const MarketData* market) {
    // call parent's getMarket
    EquityBase::getMarket(model, market);
}

/* for reflection */
SimpleEquity::SimpleEquity(): EquityBase(TYPE){}

class SimpleEquityHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SimpleEquity, clazz);
        SUPERCLASS(EquityBase);
        IMPLEMENTS(ICanBeRisky);
        IMPLEMENTS(Asset::IStruckable);
        IMPLEMENTS(ICanPhysicallySettle);
        EMPTY_SHELL_METHOD(defaultSimpleEquity);
        FIELD(equity, "Stock");
        FIELD(vol, "Volatility");
        FIELD(name, "Name");
        FIELD_MAKE_OPTIONAL(name);
    }

    static IObject* defaultSimpleEquity(){
        return new SimpleEquity();
    }
};

/** constructor  */
SimpleEquity::SimpleEquity(const Equity*   equity,
                           const CVolBase* vol,
                           CClassConstSP const type): EquityBase(type, equity, vol){}

CClassConstSP const SimpleEquity::TYPE = CClass::registerClassLoadMethod(
    "SimpleEquity", typeid(SimpleEquity), SimpleEquityHelper::load);


DRLIB_END_NAMESPACE

