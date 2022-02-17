
#include "edginc/config.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/MarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class is useful for retrieving market data when you don't have a model
    to hand. (Recall that the Model drives how the market data is retrieved
    from the cache.) */

NonPricingModel::~NonPricingModel(){}

/** Does what the name of this class suggests - throws an exception */
void NonPricingModel::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    throw ModelException("NonPricingModel::Price",
                         "This model does not price!");
}


/** Override default createMDF in order to return the local MDF */
MarketDataFetcherSP NonPricingModel::createMDF() const {
    return mdf;
}

IModel::WantsRiskMapping NonPricingModel::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

void NonPricingModel::load(CClassSP& clazz){
    REGISTER(NonPricingModel, clazz);
    SUPERCLASS(CModel);
}

NonPricingModel::NonPricingModel(): CModel(TYPE),
      mdf(new MarketDataFetcher()){}

/** uses supplied MarketDataFetcher to fetch data */
NonPricingModel::NonPricingModel(
      const MarketDataFetcherSP& mdf): CModel(TYPE), mdf(mdf){}


CClassConstSP const NonPricingModel::TYPE = CClass::registerClassLoadMethod(
    "NonPricingModel", typeid(NonPricingModel), load);

DRLIB_END_NAMESPACE
