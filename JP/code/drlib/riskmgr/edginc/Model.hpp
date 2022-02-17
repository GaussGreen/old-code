
#ifndef EDG_MODEL_H
#define EDG_MODEL_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/IModel.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

class SensControl;
class ExposureHighlighter;
class CashflowInfo;

/** Base class for all Models. Models direct the pricing of an instrument.
    Each instance of a model must specify an interface that an instrument
    must implement in order to be priced using that model */
class RISKMGR_DLL CModel: public CObject,
              public virtual IModel,
              public virtual Theta::IShift 
{
public:
    friend class ModelHelper;
    static CClassConstSP const TYPE;

    virtual ~CModel();

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

    /** Does everything, namely retrieve market data, applies scenario (if
        non null), calculates price and greeks. Static method to allow
        other implementations of IModel to use this if they want */
    static CResultsArraySP go(IModelSP                model,
                              IInstrumentCollectionSP instruments,
                              ScenarioSP              scenario, // optional
                              CControlSP              control,
                              MarketDataSP            market);

    /** Does everything, namely retrieve market data, applies scenario (if
        non null), calculates price and greeks as well as save to file of
        inputs and outputs. Static method to allow
        other implementations of IModel to use this if they want */
    static CResultsSP go(IModelSP      model,
                         CInstrumentSP instrument,
                         ScenarioSP    scenario, // optional
                         CControlSP    control,
                         MarketDataSP  market);

    /** wrapper round RunMulti */
    virtual CResults* Run(CInstrument* instrument, 
                          CControl*    control);

    /** main control - calculates prices and sensitivities. This method is
        implemented and, in general, should not be overridden */
    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments, 
                                     CControl* control);

    /** calculate prices and store results */
    virtual void PriceMulti(IInstrumentCollectionSP instruments, 
                            CControl* control, 
                            CResultsArraySP results);
    
    /** Essentially a wrapper around MarketData::GetData(model, name, type).
        Retrieves the market object with
        the specififed name and type (and retrieves its market data). The
        domesticYCName flags what the current 'domestic' currency is (see
        getDomesticYCName() for more information). */
    virtual MarketObjectSP getMarketObject(
        const MarketData* market,
        const string&     name,
        CClassConstSP     type,
        const string&     domesticYCName) const;


    /** Method to set the domestic yield curve name in the mdf  */
    virtual void setDomesticYCName (string discountYieldCurveName) const;


    /** Method to obtain the instrument and model market data. It is here, 
     * in the model, because the domestic yield curve name needs to be set
     * in the mdf attribute of the model BEFORE getting the instrument's
     * market data */
    virtual void getInstrumentAndModelMarket (const MarketData*  market,
                                              CInstrument* inst);

    /** Method to obtain the instruments' and model's market data. It is here, 
     * in the model, because the domestic yield curve name needs to be set
     * in the mdf attribute of the model BEFORE getting the instrument's
     * market data */
    virtual void getInstrumentsAndModelMarket (MarketDataConstSP market,
                                               IInstrumentCollectionSP insts);

    /** Returns a [deep] copy of the market data with supplied name
        and type from the given market data cache but does not retrieve its
        market data. This gives the
        model a chance to choose a specific type of market data rather
        than just a general instance. For example, the method could
        request a Black-Scholes Vol rather than just any old vol. The
        default implementation provided by CModel just asks the market
        data fetcher for the object of the given type */
    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const;

    /** Having retrieved a market object from the cache (or if presupplied) it
        is necessary to ensure that it has all its market data. Rather than
        call getMarket on the MarketObject, this method should be called 
        instead as it allows the model to have context information 
        regarding any calls to GetMarket. The default implementation,
        of course, just calls getMarket on the MarketObject. */
    virtual void getComponentMarketData(const MarketData*    market,
                                        MarketObjectSP       mo) const;

    /** Same as above, but tells the model what the 'domestic' currency is.
        Here domestic means what currency we want things to be in eg the
        instrument's discount currency */
    virtual void getComponentMarketData(const MarketData* market,
                                        MarketObjectSP    mo,
                                        const string&     domesticYCName) const;

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
                                             const MarketData* market) const;
                                
    /** Invoked for each piece of market data (whether already inline in
     *  instrument or pulled from cache). The default implementation just
     *  returns mo. Derived classes can use this to replace market data
     *  objects with other instances. It compliments GetMarket in that this
     *  method works for instruments which have the market data inside them
     *  already */
    virtual MarketObjectSP modifyMarketData(
        const MarketData*     market,
        const CClassConstSP&  clazz,     // what type was originally requested
        const MarketObjectSP& mo) const; /* what GetMarket returned or what was
                                          * "inline" already */

    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required. Default implementation does nothing */
    virtual void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments);

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
    virtual const string& getDomesticYCName() const;

    /** get valueDate */
    virtual DateTime getValueDate() const;

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** utility wrapper around Price */
    double calcPrice(CInstrument* instrument, 
                     CControl*    control);

    /** override a control shift (eg for delta on trees)
        returns true if new control is constructed else returns 0 */
    virtual SensControl* AlterControl(
        const SensControl* currSensControl) const;

    /** called after a series of pricings to indicate that the object
        will not be used again to calculate any
        sensitivities. Basically, state information can be stored
        inside the Model - some of which might be expensive in terms
        of memory. Since the model is constructed by clients, we do
        not control when the Model object is freed. This gives a
        chance for Models to free any expensive caches. The default
        implementation does nothing */
    virtual void flush();

    /** Marker interface for IIntoProduct */
    class RISKMGR_DLL IModelIntoProduct: virtual public IObject {
    public:
        static CClassConstSP const TYPE;
        virtual ~IModelIntoProduct(){}
    protected:
        IModelIntoProduct(){}
    };

    /** Marker interface for IProdCreator 
        this is to filter out */
    class RISKMGR_DLL IProdCreator: virtual public IObject {
    public:
        static CClassConstSP const TYPE;
        // Returns the value of the component
        // (historical value for past dates, deterministic estimation for future dates)
        // If it cannot be done (model needed), return 0. cfi->amountType updated.
        virtual double getValue(DateTime date, CashflowInfo &cfi ) const;
        // Specifies that the value will be needed at resetDates[i] (future and past fixings)
        virtual void addResetDates(const DateTimeArray &resetDates) {}
        // calls addResetDates() with one date
        void addResetDate(DateTime resetDate);

        // Similar to GetMarket but at the component level. 
        // Must be also used to call setup() on the underling IProdCreators
        virtual void setup(const IModel* model, const MarketData* market);
        virtual ~IProdCreator(){}
    protected:
        IProdCreator(){}
    private:
    	static void load(CClassSP& clazz);
    };
	DECLARE(IProdCreator);

    /** Returns the market data fetcher */
    MarketDataFetcherSP getMDF() const; 

    /** returns a PriceCounter - a "model" that does everything except
    actually price, so you can get a count of times it was asked to price */
    virtual PriceCounter* priceCounter();

    /** returns an ExposureHighlighter - a "model" that does everything except
        actually price, so you get to see what market data it uses 
        Default implementation supplied */
    virtual ExposureHighlighter* exposureHighlighter();

    /** Releases a previously created MDF */
    virtual void releaseMDF() const;

protected:
    CModel(CClassConstSP clazz);
    //// normally 'today' is populated when retrieving market data, however
    //// if this step has been missed then the date needs to be explicitly
    //// specified
    CModel(CClassConstSP clazz, const DateTime& today);


    /** Creates a MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /** Populates valueData field in this object from MarketData - this is
        is done in getInstrumentsAndModelMarket but if you override that method
        this method needs to be called */
    void fetchToday(MarketDataConstSP market);
private:
    // these 2 methods needed to avoid compiler trying to create these functions
    // (due to export) although it does not have full class info for components.
    CModel(const CModel &rhs);
    CModel& operator=(const CModel& rhs);

    /* Method to check if there is a mdf in place, or create one otherwise */
    void checkMDF() const;

    static CResultsArraySP miniGo(
        IModelSP                model,
        IInstrumentCollectionSP instruments,
        ScenarioSP              scenario, // optional
        CControlSP              control,
        MarketDataSP            market);

    /// fields //////////////
    mutable MarketDataFetcherSP mdf;  // $unregistered
    DateTime valueDate; // transient
};

typedef CModel Model;
typedef CModel::IProdCreator IProdCreator; // for convenience
DECLARE(IProdCreator);

// typedef for smart pointers to Model
typedef smartConstPtr<CModel> CModelConstSP;
typedef smartPtr<CModel> CModelSP;
typedef array<IModelSP, CModel> CModelArray;
typedef smartPtr<CModelArray> CModelArraySP;
typedef smartConstPtr<CModelArray> CModelArrayConstSP;
#ifndef QLIB_MODEL_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CModel>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CModel>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<IModelSP _COMMA_ CModel>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CModelArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CModelArray>);
EXTERN_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetSmartPtr<IModelSP>(IModelSP* t));
EXTERN_TEMPLATE(void RISKMGR_DLL FieldSetSmartPtr<IModelSP>(IModelSP* t, IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CModel>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CModel>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<IModelSP _COMMA_ CModel>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CModelArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CModelArray>);
INSTANTIATE_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetSmartPtr<IModelSP>(IModelSP* t));
INSTANTIATE_TEMPLATE(void RISKMGR_DLL FieldSetSmartPtr<IModelSP>(IModelSP* t, IObjectSP o));
#endif

//----------------------------------------------------
// Allow participating Models to be given a collective name
// Used for VolPreferred per-model-family
class RISKMGR_DLL IModelFamily {
public:
    // string* to allow easy backwards compatability via 0
    virtual const string* getFamilyName() const = 0;

    virtual ~IModelFamily() {};
};

DRLIB_END_NAMESPACE
#endif
