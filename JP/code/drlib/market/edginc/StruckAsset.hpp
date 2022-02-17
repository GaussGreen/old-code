//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StruckAsset.hpp
//
//   Description : Implement asset for ccy struck assets 
//
//   Author      : Mark A Robson
//
//   Date        : 16 Oct 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_STRUCK_ASSET_HPP
#define EDR_STRUCK_ASSET_HPP
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/AssetFairValue.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of CAsset for ccy struck assets - we need the asset we're
    protecting to implement an additional interface in order to handle 
    the vol */
class MARKET_DLL StruckAsset: public CAsset,
                   public virtual CAsset::IStruck,
                   public virtual INextStrike,
                   public virtual Theta::IShift,
                   public virtual IAssetFairValue {
public:
    static CClassConstSP const TYPE;
    friend class StruckAssetHelper;

    virtual ~StruckAsset();

    /** Constructor needed for case when instrument specifies
        underlying asset and currency treatment */
    StruckAsset(const string& name,         // name for struck asset
                const string& assetName,    // asset to make struck
                const string& fxAssetName); // name for fx asset

    /** Validation */
    virtual void validatePop2Object();

    /** Pull out the asset date and correlation from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** returns the equity's name */
    virtual string getTrueName() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;


    /** Calculates the expected spot price of the asset at the given dates */
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
    virtual DateTime settleDate(const DateTime& tradeDate) const;

    /** (IStruck interface) Returns the fx spot */
    virtual double getFXSpot() const;

    /** (IStruck interface) Returns the fx forward value */
    virtual double fxFwdValue(const DateTime& date) const;

    /** Shifts the object using given shift. */
    bool sensShift(Theta* shift);

    /** overrides CObject version to allow for easy default */
    virtual bool accept(ICollector* collector) const;

   /** returns sensitive strikes for a given vol request */
    virtual void getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const;

    /** given a current spot level, get the next strike on the vol
        surface where the slope is non-differentiable */
    virtual double getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** returns dividend list - to be retired - it's not the right way to 
        do this */
    DividendListSP  getAllDivsBetweenDates(const DateTime& start,
                                           const DateTime& end) const;

    /** Returns fair value of asset */
    virtual double fairValue() const;

    // Tactical method to get (const) underlying (non-struck) asset 
    CAssetConstSP  getPlainAsset() const;

    //// Returns the ccy treatment for the asset.  
    virtual string getCcyTreatment() const;

    void checkStruckSourceAndObsType(const ObservationSource* source,
                                     const ObservationType* type) const;

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
    virtual bool isHoliday(const DateTime&            sampleDate,
                           const ObservationSource*   source) const;
private:
    StruckAsset();
    StruckAsset(const StruckAsset& rhs);
    StruckAsset& operator=(const StruckAsset& rhs);

    static void acceptValueDateCollector(const StruckAsset*   asset, 
                                         CValueDateCollector* collector);

    static void acceptImntCcy(const StruckAsset* asset,
                              AssetCcyCollector* collector);

    static void acceptNameCollector(const StruckAsset*  asset, 
                                    AssetNameCollector* collector);

    string             name;     // struck asset's name
    CAssetWrapper      asset;    // non struck asset
    FXAssetWrapper     fx;       // what its struck into
    CorrelationSP      corrEqFx; // correlation between equity & fx
    DateTime           baseDate; // the value date
    // transient
    bool               getCorrFromCache; // true: get correlation from cache
};

typedef smartPtr<StruckAsset> StruckAssetSP;
typedef smartConstPtr<StruckAsset> StruckAssetConstSP;

DRLIB_END_NAMESPACE
#endif
