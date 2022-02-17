//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : MultiModel.hpp
//
//   Description : A model which contains an array of other models
//
//   Author      : Linus Thand
//
//   Date        : 24 May 2006
//
//----------------------------------------------------------------------------



#ifndef EDG_MULTI_MODEL_H
#define EDG_MULTI_MODEL_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IModel.hpp"

DRLIB_BEGIN_NAMESPACE

/** MultiModel contains an array of models, and delegates the pricing of 
    each instrument in an InstrumentCollection to the appropriate model,
    appropriate meaning Model X prices Instrument X, Model X+1 prices 
    Instrument X+1 et.c. 
    The model will not run unless the number of instruments is the same as the 
    number of wrapped models. This is used to price an InstrumentCollection 
    where not all Instruments should be priced by the same models, e.g. when 
    solving for reinvestment in a cashflow after a name substitution on a 
    CDO portfolio (i.e. a managed trade).
*/

class MultiModel: public CObject,
                  public virtual IModel
{
 public:
    static CClassConstSP const TYPE;

    virtual PriceCounter* priceCounter();
 
    virtual ~MultiModel();

    IModelSP getModel(int n) const;

    /** Does everything, namely retrieve market data, applies scenario 
        (if non null), calculates price and greeks. Default implementation
        here just calls static go method below. For models which want to
        delegate *everything* (including market data retrieval etc) to another
        model, then this method and the other 'go' method can be overridden. */
    virtual CResultsArraySP go(IInstrumentCollectionSP instruments,
                               ScenarioSP              scenario, // optional
                               CControlSP              control,
                               MarketDataSP            market);

    /** Same as above go method above but for on a single instrument. The two
        methods are separate to allow eg a different write to file. This 
        isn't ideal but keeps the ability to create a regression file which
        actually mirrors what was actually done. Alternatives would be a method
        on IInstrumentCollection or some ugly switch on type of instruments */

    virtual CResultsSP go(CInstrumentSP instrument,
                          ScenarioSP    scenario, // optional
                          CControlSP    control,
                          MarketDataSP  market);

    virtual void flush();

    virtual void PriceMulti(IInstrumentCollectionSP instruments, 
                            CControl*               control, 
                            CResultsArraySP         results);

    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments, 
                                     CControl*               control);

    /** The following public methods should never be used, they should already
    have been delegated to the wrapped models.
    */

    virtual CResults* Run(CInstrument* instrument, CControl* control);
 
    virtual MarketObjectSP getMarketObject(const MarketData* market,
                                           const string&     name,
                                           CClassConstSP     type,
                                           const string&     domesticYCName) 
                                           const;

    virtual void setDomesticYCName (string discountYieldCurveName) const;

    virtual void getInstrumentAndModelMarket (const MarketData* market,
                                              CInstrument*      inst);

    virtual void getInstrumentsAndModelMarket (MarketDataConstSP       market,
                                               IInstrumentCollectionSP insts);

    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const;

    virtual void getComponentMarketData(const MarketData* market,
                                        MarketObjectSP    mo) const;

    virtual void getComponentMarketData(const MarketData* market,
                                        MarketObjectSP    mo,
                                        const string&     domesticYCName) const;

    virtual CorrelationBaseSP getCorrelation(const string&     corrName,
                                             CClassConstSP     type1,
                                             CClassConstSP     type2,
                                             CClassConstSP     corrType,
                                             const MarketData* market) const;
                                
    virtual MarketObjectSP modifyMarketData(const MarketData*     market,
                                            const CClassConstSP&  clazz,
                                            const MarketObjectSP& mo) const; 

    virtual void getMarket(const MarketData*       market,
                           IInstrumentCollectionSP instruments);

    virtual const string& getDomesticYCName() const;

    virtual DateTime getValueDate() const;

    double calcPrice(CInstrument* instrument, CControl* control);

    virtual SensControl* AlterControl(const SensControl* currSensControl) const;

    MarketDataFetcherSP getMDF() const; 

    virtual void Price(CInstrument* instrument, 
                       CControl*    control, 
                       CResults*    results);
                                           
    virtual WantsRiskMapping wantsRiskMapping() const;

    virtual ExposureHighlighter* exposureHighlighter();

 private:
    MultiModel(CClassConstSP clazz);
    MultiModel();
    MultiModel(const MultiModel& rhs); 
    MultiModel& operator=(const MultiModel& rhs);
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz); 

    CResultsArraySP miniGo(IInstrumentCollectionSP instruments,
                           ScenarioSP              scenario, // optional
                           CControlSP              control,
                           MarketDataSP            market);

    /* Fields */
    IModelArraySP models; // An array of wrapped models
};

FORWARD_DECLARE(MultiModel);


DRLIB_END_NAMESPACE
#endif
