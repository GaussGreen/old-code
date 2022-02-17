
#ifndef NON_PRICING_MODEL_HPP
#define NON_PRICING_MODEL_HPP

#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class is useful for retrieving market data when you don't have a model
    to hand. (Recall that the Model drives how the market data is retrieved
    from the cache.) */
class RISKMGR_DLL NonPricingModel: public CModel{
public:
    static CClassConstSP const TYPE;
    ~NonPricingModel();

    /** Does what the name of this class suggests - throws an exception */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    /** uses MarketDataFetcher to fetch data */
    NonPricingModel();

    /** uses supplied MarketDataFetcher to fetch data */
    NonPricingModel(const MarketDataFetcherSP& mdf);

    /** Override default createMDF in order to return the local MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Just returns riskMappingIrrelevant.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

private:
    MarketDataFetcherSP   mdf; // $unregistered
    static void load(CClassSP& clazz);
};

DECLARE(NonPricingModel);

DRLIB_END_NAMESPACE
#endif
