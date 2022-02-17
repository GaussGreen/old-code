//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StruckEquity.hpp
//
//   Description : Implement asset for simple ccy struck equity 
//
//   Author      : Mark A Robson
//
//   Date        : 8 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_STRUCK_EQUITY_HPP
#define EDG_STRUCK_EQUITY_HPP
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/SpotLevelProbability.hpp"
#include "edginc/CanBeRisky.hpp"
#include "edginc/HaveEquity.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/IEqVolNamePair.hpp"

DRLIB_BEGIN_NAMESPACE

/* You could derive this class from SimpleEquity (or something similar) but
   given its very small the extra clarity of a separate implementation both
   for this and SimpleEquity seems worth it */

/** Implementation of CAsset for a simply ccy struck equity. */
class MARKET_DLL StruckEquity: public CAsset,
                    virtual public CAsset::IStruck,
                    virtual public VegaSkewParallel::IShift,
                    virtual public VegaSkewPointwise::IShift,
                    virtual public INextStrike,
                    virtual public SpotLevelProbability::Shift,
                    virtual public ICanBeRisky,
                    virtual public IHaveEquity,
                    virtual public DeltaSurface::IShift,
                    virtual public IAssetFairValue,
                    virtual public VolRelativeShift::IShift,
                    virtual public IEqVolNamePair {
public:
    static CClassConstSP const TYPE;
    friend class StruckEquityHelper;

    ~StruckEquity();

    /** Constructor needed for case when instrument specifies
        underlying asset and currency treatment */
    StruckEquity(const string& name, // name for struck asset
                 const Equity* equity,
                 const string& volName, // name for vol
                 const string& fxAssetName); // name for fx asset

    StruckEquity(const string& name, // name for struck asset
                 const Equity* equity,
                 const CVolBaseWrapper    vol, // name for vol
                 const FXAssetWrapper     fx,
				 const Correlation*       corrEqFx); // name for fx asset 

    /** Pull out the vol, fx asset and correlation from the
        market data */
    void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** returns the equity's name */
    virtual string getTrueName() const;

    /** Returns fair value of asset */
    virtual double fairValue() const;

    /** get a (const) plain asset */
    CAssetConstSP  getPlainAsset() const;

    /** could replace this by a collector */
    string getVolName() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;


    /** Calculates the expected spot price of the asset at the given dates */
    void fwdValue(const DateTimeArray& dates,
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

    /** (IStruck interface) Returns the fx spot */
    virtual double getFXSpot() const;

    /** (IStruck interface) Returns the fx forward value */
    virtual double fxFwdValue(const DateTime& date) const;
    
    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    string sensName(VegaSkewParallel* shift) const;

    /** Shifts the object for VegaSkewParallel */
    bool sensShift(VegaSkewParallel* shift);

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    string sensName(VegaSkewPointwise* shift) const;
    
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;

    /** Shifts the object for VegaSkewPointwise */
    bool sensShift(VegaSkewPointwise* shift);

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

    FXAssetConstSP getFX() const;
	
    /** return the VolBaseWrapper */
    const CVolBaseWrapper& getVol()const;

    /** Is the Eq/FX correlation in the cache ? */
    bool useCorrFromCache()const;
    
    /** return the Eq/FX correlation */
    const Correlation* getCorrelation()const;

    /** return the FXAssetWrapper */
    const FXAssetWrapper& getFXAsset()const;
    
    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** adds the credit spread to the asset's growth curve */
    virtual void makeRisky(ICreditCurveSP creditSpreads,
        const  DateTime *maturityDate=0);

    /** returns a (copied) smart pointer to the equity */
    virtual EquitySP getEquity() const { return equity; }

    /* for IEqVolNamePair */
    virtual bool getNamePairs(string& eqName, string& volName) const;

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
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const;
private:
    StruckEquity();
    StruckEquity(const StruckEquity& rhs);
    StruckEquity& operator=(const StruckEquity& rhs);
    void validatePop2Object();

    static void acceptValueDateCollector(const StruckEquity*  asset, 
                                         CValueDateCollector* collector);

    static void acceptDeltaShift(const StruckEquity* asset, 
                                 ShiftSizeCollector* collector);

    static void acceptImntCcy(const StruckEquity* asset,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const StruckEquity* asset,
                              HolidayCollector*   collector);


    static void acceptCriticalDateCollector(const StruckEquity*    asset,
                                            CriticalDateCollector* collector);

    string             name;     // asset's name
    EquitySP           equity;   // how to grow stock forward
    CVolBaseWrapper    vol;      // the stock's vol
    FXAssetWrapper     fx;       // what its struck into
    CorrelationSP      corrEqFx; // correlation between equity & fx
    // transient
    bool               getCorrFromCache; // true: get correlation from cache
};

typedef smartPtr<StruckEquity> StruckEquitySP;
typedef smartConstPtr<StruckEquity> StruckEquityConstSP;

DRLIB_END_NAMESPACE
#endif
