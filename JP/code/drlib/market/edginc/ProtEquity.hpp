#ifndef PROT_EQUITY_HPP
#define PROT_EQUITY_HPP

#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/FXVolBase.hpp"
#include "edginc/SpotLevelProbability.hpp"
#include "edginc/HaveEquity.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/IEqVolNamePair.hpp"

DRLIB_BEGIN_NAMESPACE

/** Simple implementation of volatility where the volatility is captured by
    a single scalar. */
class MARKET_DLL ProtEquity: public CAsset,
                  virtual public Asset::IQuanto,
                  virtual public VegaSkewParallel::IShift,
                  virtual public VegaSkewPointwise::IShift,
                  virtual public Theta::IShift,
                  virtual public INextStrike,
                  virtual public SpotLevelProbability::Shift,
                  virtual public IHaveEquity,
                  virtual public DeltaSurface::IShift,
                  virtual public VolRelativeShift::IShift,
                  virtual public IEqVolNamePair,
                  virtual public ObjectIteration::IOverride {
public:
    static CClassConstSP const TYPE;
    friend class ProtEquityHelper;

    /** Pull out the base date, vol, fx vol and correlation from the
        market data */
    void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** returns the equity's name */
    virtual string getTrueName() const;

    YieldCurveWrapper  getYC() const;

    /** get a (const) plain asset */
    CAssetConstSP  getPlainAsset() const;

    // access the FX asset
    const Asset* getFXAsset() const;

    /** could replace this by a collector */
    string getVolName() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** <<Asset::IQuanto>> implementation.
     *  Calculates the expected spot price of the underlying asset at
     *  the given date ie no currency protection adjustment is made */
    virtual double unadjustedFwdValue(const DateTime& date) const;

    /** <<Asset::IQuanto>> implementation.
     *  Calculates the expected spot price of the underlying asset at
     *  the given dates ie no currency protection adjustment is made */
    virtual void unadjustedFwdValue(const DateTimeArray& dates,
                                    CDoubleArray&        result) const;

    /** <<Asset::IQuanto>> implementation.
     *  return the Eq/FX correlation */
    virtual const Correlation* getCorrelation() const;

    /** <<Asset::IQuanto>> implementation.
     *  Returns a processed vol - which combines the vol market data with
     *  the instrument data in the volRequest */
    virtual IVolProcessed* getProcessedFXVol( const CVolRequest* volRequest ) const;

    /** Calculates the expected spot price of the asset at each of the
        given dates */
    virtual void fwdValue(const DateTimeArray& dates,
                          CDoubleArray&        result) const;

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const;

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const;

    /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;
    

    /** Constructor needed for case when instrument specifies
        underlying asset and currency treatment. */
    ProtEquity(const string&        name,
               const Equity*        equity,
               const string&        eqVolName,
               const string&        fxVolName,
               const string&        protCcyCode,
               bool                 homoGreeks);


    ProtEquity(const string&        name,
               const Equity*        equity,
               const CVolBaseWrapper  vol,
               const FXVolBaseWrapper fxVol,
               const Correlation*         corrEqFX,
               const DateTime&             baseDate,
               const string&           protCcyIsoCode,
               const string&        protCcyCode,
               const double         scale,
               bool                 homoGreeks);

    /** validation code to be called after object construction */
    virtual void validatePop2Object();

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    virtual string sensName(VegaSkewParallel* shift) const;
    /** Override the name matching mechanism to deal with fx vol */
    bool sensNameMatches(VegaSkewParallel* shift, const OutputName& name) const;
    /** Override the name adding mechanism to deal with fx vol */
    void sensAppendName(VegaSkewParallel* shift, OutputNameArray& namesList) const;
    /** Shifts the object for VegaSkewParallel */
    virtual bool sensShift(VegaSkewParallel* shift);

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    virtual string sensName(VegaSkewPointwise* shift) const;
    /** Override the name matching mechanism to deal with fx vol */
    bool sensNameMatches(VegaSkewPointwise* shift, const OutputName& name) const;
    /** Override the name adding mechanism to deal with fx vol */
    void sensAppendName(VegaSkewPointwise* shift, OutputNameArray& namesList) const;
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    /** Shifts the object for VegaSkewPointwise */
    virtual bool sensShift(VegaSkewPointwise* shift);

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /// implementation of DeltaSurface::IShift interface
    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

    // implementation of VolRelativeShift::IShift interface
    virtual string sensName(VolRelativeShift* shift) const;
    virtual bool sensShift(VolRelativeShift* shift);

    /** returns sensitive strikes for a given vol request */
    virtual void getSensitiveStrikes(
        const CVolRequest* volRequest,
        OutputNameConstSP outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP sensitiveStrikes) const;

    /** given a current spot level, get the next strike on the vol surface where
        the slope is non-differentiable */
    virtual double getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const;

    /** Returns the name of the stock/asset - used to determine
        whether to shift the object */
    virtual string sensName(SpotLevelProbability* shift) const;

    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(SpotLevelProbability* shift);

    /** returns dividend list */
    DividendListConstSP  getDivList() const;

    /** return the VolBaseWrapper */
    const CVolBaseWrapper& getVol()const;
    
    /** return the FXVolBaseWrapper */
    const FXVolBaseWrapper& getFXVol()const;

    /** Is the Eq/FX correlation in the cache ? */
    bool useCorrFromCache()const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** returns a (copied) smart pointer to the equity */
    virtual EquitySP getEquity() const { return equity; }

    /* for EqVolNamePair */
    virtual bool getNamePairs(string& eqName, string& volName) const;

    //// Returns the ccy treatment for the asset.  
    virtual string getCcyTreatment() const;

    /** Can this asset physically settle? */
    virtual bool canPhysicallySettle() const;

    virtual bool recurse(const CFieldConstSP& field,
                         const CClassConstSP& targetClass) const;

    // the IMarketObservable interface for retrieving a single sample
    virtual double pastValue(const DateTime&             sampleDate,
                             const ObservationType*      obsType,
                             const ObservationSource*    source,
                             const FixingType*           fixType,
                             const IObservationOverride* overrides,
                             const SamplingConvention*   sampleRule) const;

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;

    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const;

    // the IMarketObservable interface for 
    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const;
private:
    ProtEquity();
    ProtEquity(const ProtEquity& rhs);
    ProtEquity& operator=(const ProtEquity& rhs);
    /** calculates protected adjustment for a set of dates */
    CDoubleArraySP protAdjustment(const DateTimeArray& dates,
                                  const double         spotPrice,
                                  const CDoubleArray&  unprotectedFwds) const;

    static void acceptValueDateCollector(const ProtEquity* asset, 
                                         CValueDateCollector* collector);
    static void acceptDeltaShift(const ProtEquity* asset,
                                 ShiftSizeCollector* collector);
    static void acceptImntCcy(const ProtEquity* asset,
                              AssetCcyCollector* collector);
    static void acceptHoliday(const ProtEquity* asset,
                              HolidayCollector* collector);
    static void acceptCriticalDateCollector(const ProtEquity* asset,
                                            CriticalDateCollector* collector);


    /////////// fields /////////////////
    string           name;     // asset's name
    EquitySP         equity;
    CVolBaseWrapper  vol;
    FXVolBaseWrapper fxVol;         // vol of FX rate 
    CorrelationSP    corrEqFX;      // equity fx correlation 
    DateTime         baseDate;      // the value date
    string           protCcyCode;   // what we're protecting into 
    // transient
    bool             getCorrFromCache; // true: get correlation from cache
    double           scale;            // for scaling equity fwd price by ccy adjustment
    string           protCcyIsoCode;   // the ISO code of the protected ccy
    bool             coherentGreeks;   // whether to use homogeneous Greeks e.g. Delta vs FXDelta 
    FXAssetWrapper   fx;               // the FX for the FX vol
};

typedef smartPtr<ProtEquity> ProtEquitySP;
typedef smartConstPtr<ProtEquity> ProtEquityConstSP;

DRLIB_END_NAMESPACE
#endif // PROT_EQUITY_HPP
