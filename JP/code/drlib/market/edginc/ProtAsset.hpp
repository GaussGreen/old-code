//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ProtAsset.hpp
//
//   Description : Generalized protected asset
//
//   Author      : Andrew J Swain
//
//   Date        : 2 October 2001
//
//
//----------------------------------------------------------------------------

#ifndef PROTASSET_HPP
#define PROTASSET_HPP

#include "edginc/Asset.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/FXVolBase.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

/** Generalized protected asset */
class MARKET_DLL ProtAsset: public CAsset,
                 virtual public Asset::IQuanto,
                 virtual public Theta::IShift,
                 virtual public VegaSkewParallel::IShift,
                 virtual public VegaSkewPointwise::IShift,
                 virtual public INextStrike,
                 virtual public ObjectIteration::IOverride {
public:
    static CClassConstSP const TYPE;
    friend class ProtAssetHelper;

    /** overrides CObject version to allow for easy default */
    virtual bool accept(ICollector* collector) const;
    
    /** Pull out the base date, vol, fx vol and correlation from the
        market data */
    void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** returns the 'true' name of the asset that is being protected */
    virtual string getTrueName() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const;

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

   /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;
    

    /** Constructor needed for case when instrument specifies
        underlying asset and currency treatment. */
    ProtAsset(const string& name,
              const Asset*  asset,
              const string& fxVolName,
              const string& protCcyCode,
              bool homoGreeks);


    /** validation code to be called after object construction */
    virtual void validatePop2Object();

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** Returns the name of the vol - used to determine whether to tweak
        the object. This is used to facilitate correct spot values on the FX vol */
    virtual string sensName(VegaSkewParallel* shift) const;
    /** Override the name matching mechanism to deal with fx vol */
    bool sensNameMatches(VegaSkewParallel* shift, const OutputName& name) const;
    /** Override the name adding mechanism to deal with fx vol */
    void sensAppendName(VegaSkewParallel* shift, OutputNameArray& namesList) const;
    /** Shifts the object for VegaSkewParallel */
    virtual bool sensShift(VegaSkewParallel* shift);

    /** Returns the name of the vol - used to determine whether to tweak
        the object This is used to facilitate correct spot values on the FX vol  */
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

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                CResults*       results,
                                const DateTime& maturityDate) const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** returns dividend list - to be retired - it's not the right way to 
        do this */
    DividendListSP  getAllDivsBetweenDates(const DateTime& start,
                                           const DateTime& end) const;

    //// Returns the ccy treatment for the asset.  
    virtual string getCcyTreatment() const;

    /** Can this asset physically settle? */
    virtual bool canPhysicallySettle() const;

    virtual bool recurse(const CFieldConstSP& field,
                         const CClassConstSP& targetClass) const;

    /** get a (const) plain asset */
    CAssetConstSP  getPlainAsset() const;

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
    ProtAsset();
    ProtAsset(const ProtAsset& rhs);
    ProtAsset& operator=(const ProtAsset& rhs);
    /** calculates protected adjustment for a set of dates */
    CDoubleArraySP protAdjustment(const DateTimeArray& dates,
                                  const double         spotPrice,
                                  const CDoubleArray&  unprotectedFwds) const;

    static void acceptValueDateCollector(const ProtAsset*     asset, 
                                         CValueDateCollector* collector);
    static void acceptImntCcy(const ProtAsset*   asset,
                              AssetCcyCollector* collector);
    static void acceptNameCollector(const ProtAsset*    asset, 
                                    AssetNameCollector* collector);

    /////////// fields /////////////////
    string            name;          // asset's name
    CAssetWrapper     asset;
    FXVolBaseWrapper  fxVol;         // vol of FX rate 
    CorrelationSP     corrAssetFX;   // asset fx correlation 
    DateTime          baseDate;      // the value date
    YieldCurveWrapper protYC;        /* what we're protecting into - note that
                                        we don't populate this wrapper to 
                                        avoid generating rho sensitivities */
    // transient
    string            protYCCode;       // iso code of protYC
    bool              getCorrFromCache; // true: get correlation from cache
    double            scale;            // for scaling asset fwd price by ccy adjustment
    bool              coherentGreeks;   // whether to use homogeneous Greeks e.g. Delta vs FXDelta 
    FXAssetWrapper    fx;               // the FX for the FX vol
};

typedef smartPtr<ProtAsset> ProtAssetSP;
typedef smartConstPtr<ProtAsset> ProtAssetConstSP;

DRLIB_END_NAMESPACE
#endif
