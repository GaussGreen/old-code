//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IModel.hpp
//
//   Description : Interface class for Models
//
//   Author      : Jose Hilera
//
//   Date        : 28 June 2005
//
//----------------------------------------------------------------------------


#ifndef EDG_I_MODEL_H
#define EDG_I_MODEL_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IInstrumentCollection);
FORWARD_DECLARE(MarketObject);
FORWARD_DECLARE(Scenario);
FORWARD_DECLARE(CInstrument);
FORWARD_DECLARE(IInstrumentCollection);
FORWARD_DECLARE(CorrelationBase);
FORWARD_DECLARE(MarketDataFetcher);
FORWARD_DECLARE(MarketData);
FORWARD_DECLARE(IModel);

class PriceCounter;
class ExposureHighlighter;
class SensControl;

/** Interface class for Models */
class RISKMGR_DLL IModel: public virtual IObject {
public:
    static CClassConstSP const TYPE;

    IModel();

    virtual ~IModel();

    /** Does everything, namely retrieve market data, applies scenario 
        (if non null), calculates price and greeks. Default implementation
        here just calls static go method below. For models which want to
        delegate *everything* (including market data retrieval etc) to another
        model, then this method and the other 'go' method can be overridden. */
    virtual CResultsArraySP go(IInstrumentCollectionSP instruments,
                               ScenarioSP              scenario, // optional
                               CControlSP              control,
                               MarketDataSP            market) = 0;

    /** Same as above go method above but for on a single instrument. The two
        methods are separate to allow eg a different write to file. This 
        isn't ideal but keeps the ability to create a regression file which
        actually mirrors what was actually done. Alternatives would be a method
        on IInstrumentCollection or some ugly switch on type of instruments */
    virtual CResultsSP go(CInstrumentSP instrument,
                          ScenarioSP    scenario, // optional
                          CControlSP    control,
                          MarketDataSP  market) = 0;

    /** main control - calculates price and sensitivities */
    virtual CResults* Run(CInstrument* instrument, 
                          CControl*    control) = 0;
    
    /** main control - calculates price and sensitivities */
    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments,
                                     CControl* control) = 0;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results) = 0;

    /** calculate prices and store results */
    virtual void PriceMulti(IInstrumentCollectionSP instruments, 
                            CControl* control, 
                            CResultsArraySP results) = 0;

    /** Retrieves the market object with the specififed name and type 
     * (and retrieves its market data) */
    virtual MarketObjectSP getMarketObject(
        const MarketData* market,
        const string&     name,
        CClassConstSP     type,
        const string&     domesticYCName) const = 0;

    virtual void setDomesticYCName (string discountYieldCurveName) const = 0;

    /** Obtains the instrument and model market data. It is here, 
     * in the model, because the domestic yield curve name needs to be set
     * BEFORE getting the instrument's market data */
    virtual void getInstrumentAndModelMarket (const MarketData*  market,
                                              CInstrument* inst) = 0;

    /** Method to obtain the instruments' and model's market data. It is here, 
     * in the model, because the domestic yield curve name needs to be set
     * in the mdf attribute of the model BEFORE getting the instrument's
     * market data */
    virtual void getInstrumentsAndModelMarket (MarketDataConstSP market,
                                               IInstrumentCollectionSP insts) = 0;

    /** Returns a [deep] copy of the market data with supplied name
        and type from the given market data cache but does not retrieve its
        market data. This gives the
        model a chance to choose a specific type of market data rather
        than just a general instance. For example, the method could
        request a Black-Scholes Vol rather than just any old vol */
    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const = 0;

    /** Having retrieved a market object from the cache (or if presupplied) it
        is necessary to ensure that it has all its market data. Rather than
        call getMarket on the MarketObject, this method should be called 
        instead as it allows the model to have context information 
        regarding any calls to GetMarket */
    virtual void getComponentMarketData(const MarketData*    market,
                                        MarketObjectSP       mo) const = 0;

    /** Same as above, but tells the model what the 'domestic' currency is.
        Here domestic means what currency we want things to be in eg the
        instrument's discount currency */
    virtual void getComponentMarketData(const MarketData* market,
                                        MarketObjectSP    mo,
                                        const string&     domesticYCName) const = 0;

    /** Retrieve a correlation from the cache with the supplied name and
     *  type specified by corrType. type1 and type2 are the types of the
     *  two objects with respect to which this is a correlation. The
     *  return correlation will be correcly configured for the
     *  appropriate sensitivity (eg PHI, FX_PHI) etc. Note that a model
     *  may decide that it does not use this type of correlation (eg IR v EQ) 
     *  and return a 'dummy' version if say the object does not live in the
     *  cache. In a similar manner, the returned object may be sensitive to
     *  no sensitivities if the model is not going to use the correlation */
    virtual CorrelationBaseSP getCorrelation(const string&     corrName,
                                             CClassConstSP     type1,
                                             CClassConstSP     type2,
                                             CClassConstSP     corrType,
                                             const MarketData* market) const = 0;
                                
    /** Invoked for each piece of market data (whether already inline in
     *  instrument or pulled from cache). It complements GetMarket in that 
     *  this method works for instruments which have the market data inside 
     *  them already */
    virtual MarketObjectSP modifyMarketData(
        const MarketData*     market,
        const CClassConstSP&  clazz,         // what type was originally requested
        const MarketObjectSP& mo) const = 0; /* what GetMarket returned or what was
                                              * "inline" already */

    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required */
    virtual void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments) = 0;

    /** Return the name of the 'domestic' yield curve - only available (at best)
        when fetching market data. Returns an empty string if not available.
        Some points to note:
        1. This is really a property of the MarketDataFetcher but currently
        we have a Model around it. 
        2. Here 'domestic' means, for example, the instrument currency. It 
        reflects the currency we want to transform, eg an asset, into. 
        It is therefore not necessarily the instrument currency (eg stuck asset
        inside a protected XCB).
        3. The motivation for this method is for CDSParSpreads (eg CDO) where
        they can be buried many layers down making it difficult to pass the
        additional data needed down */
    virtual const string& getDomesticYCName() const = 0;

    /** get valueDate */
    virtual DateTime getValueDate() const = 0;

    /** utility wrapper around Price */
    virtual double calcPrice(CInstrument* instrument, CControl* control) = 0;

    /** override a control shift (eg for delta on trees)
        returns true if new control is constructed else returns 0 */
    virtual SensControl* AlterControl(
         const SensControl* currSensControl) const = 0;

    /** possible return values for wantsRiskMapping() */
    enum WantsRiskMapping {
        riskMappingIrrelevant, riskMappingDisallowed, riskMappingAllowed
    };

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * RiskMapping enables us to get (an approximation to) Black-Scholes Delta,
     * VegaParallel etc. for instruments priced with parametric vol (VolSVJ,
     * SRMEQ::Vol, ...: anything carrying the IDynamicsParameter tag
     * interface).
     *
     * Most models are defined directly in Black-Scholes terms (using
     * VolSurface or whatever) and should return riskMappingIrrelevant.
     * Example: MonteCarloImpliedDefault::wantsRiskMapping().
     *
     * Models which may use an IDynamicsParameter parametric vol should
     * have a flag <tt>allowRiskMapping</tt>, exposed in the public
     * interface, and return
     *
     *    -  riskMappingIrrelevant if they're not using parametric
     *       vol in this particular instance;
     *    -  riskMappingDisallowed if they are but the flag is off;
     *    -  riskMappingAllowed if they are and the flag is on.
     *
     * Example: FourierModel::wantsRiskMapping().
     *
     * Compound models built up from several sub-models should return
     *
     *    -  riskMappingDisallowed if any of the sub-models do;
     *    -  otherwise, riskMappingAllowed if any of the sub-models do;
     *    -  otherwise, riskMappingIrrelevant.
     *
     * Example: IAggregateModel::wantsRiskMapping(),
     * VarCapCVModel::wantsRiskMapping().
     */

    virtual WantsRiskMapping wantsRiskMapping() const = 0;

    /** called after a series of pricings to indicate that the object
        will not be used again to calculate any
        sensitivities. Basically, state information can be stored
        inside the Model - some of which might be expensive in terms
        of memory. Since the model is constructed by clients, we do
        not control when the Model object is freed. This gives a
        chance for Models to free any expensive caches */
    virtual void flush() = 0;

    /** Returns the market data fetcher */
    virtual MarketDataFetcherSP getMDF() const = 0;

    /** returns a PriceCounter - a "model" that does everything except
    actually price, so you can get a count of the number of times it was asked to price */
    virtual PriceCounter* priceCounter() = 0;

    /** returns an ExposureHighlighter - a "model" that does everything except
        actually price, so you get to see what market data it uses 
        Default implementation supplied */
    virtual ExposureHighlighter* exposureHighlighter() = 0;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
