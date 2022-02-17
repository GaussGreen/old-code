//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PercXCBWithVol.hpp
//
//   Description : Implement asset for percentage weighted cross 
//                 currency basket with vol override
//
//   Author      : Mark A Robson
//
//   Date        : 9 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_PERC_XCB_WITH_VOL_HPP
#define EDG_PERC_XCB_WITH_VOL_HPP
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/StartDateCollector.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/FutureExpiryCollector.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/SpotShift.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/BasketSpot.hpp"
#include "edginc/IEqVolNamePair.hpp"

DRLIB_BEGIN_NAMESPACE

/* Separate implementation from other xcb's for clarity */

/** Implementation of CAsset for percentage weighted cross currency
    basket with vol override. */
class MARKET_DLL PercXCBWithVol: public CAsset,
    // public DeltaBasket::IShift, (need to decide what we're doing)
                      virtual public IAssetFairValue,
                      virtual public VegaSkewParallel::IShift,
                      virtual public VegaSkewPointwise::IShift,
                      virtual public ITweakableWithRespectTo<BasketSpot>,
                      virtual public Theta::IShift,
                      virtual public DeltaSurface::IShift,
                      virtual public VolRelativeShift::IShift,
                      virtual public SpotShift::Shift,
                      virtual public IEqVolNamePair {
public:
    static CClassConstSP const TYPE;
    friend class PercXCBWithVolHelper;

    /** Validation */
    virtual void validatePop2Object();

    /** Pull out the component assets from the market data */
    void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

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

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

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

//    /** Calculates the expected spot price of the asset at the given date if
//        the spot price had the given value spot on spotDate */
//    virtual double fwdFwd(const DateTime& spotDate,
//                          double          spot, 
//                          const DateTime& fwdDate) const;

   /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;

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

    /** Shifts the object using given shift (see Theta::Shift)*/
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

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                CResults*       results,
                                const DateTime& maturityDate) const;

    /** Returns the name of the XCB for Delta */
    virtual string sensName(const BasketSpot*) const;
    /** Shifts all components in the XCB for Delta */
    virtual TweakOutcome sensShift(const PropertyTweak<BasketSpot>& shift);

    /** Returns fair value of XCB */
    virtual double fairValue() const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** Returns name identifying this object for SpotShift */
    virtual string sensName(SpotShift* shift) const;
    /** Shifts the object using given shift (see SpotShift::Shift)*/
    virtual bool sensShift(SpotShift* shift);

    /* for IEqVolNamePair */
    virtual bool getNamePairs(string& eqName, string& volName) const;

    /** Can this asset physically settle? */
    virtual bool canPhysicallySettle() const;

private:
    PercXCBWithVol();
    PercXCBWithVol(const PercXCBWithVol& rhs);
    PercXCBWithVol& operator=(const PercXCBWithVol& rhs);

    /** Returns a reference to the internal weights to use for
        combining assets */
    const CDoubleArray& weightsRef() const;

    static void acceptCollector(const PercXCBWithVol* asset, 
                                StartDateCollector* collector);

    static void acceptNameCollector(const PercXCBWithVol* asset, 
                                    AssetNameCollector* collector);

    static void acceptFutureCollector(const PercXCBWithVol* asset, 
                                      FutureExpiryCollector* collector);

    static void acceptValueDateCollector(const PercXCBWithVol* asset, 
                                         CValueDateCollector* collector);

    static void acceptDeltaShift(const PercXCBWithVol* asset, 
                                 ShiftSizeCollector* collector);

    static void acceptImntCcy(const PercXCBWithVol* asset,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const PercXCBWithVol* asset,
                              HolidayCollector* collector);

    static void acceptWrapperNameCollector(const PercXCBWithVol* asset, 
                                           WrapperNameCollector* collector);

    static void acceptCriticalDateCollector(const PercXCBWithVol* asset, 
                                            CriticalDateCollector* collector);
    
    /** Helper for centralised sampling/isda adjustment for XCBs
    gets a list of dates for which the given sample will finally be known for 
    each component. Results are put into obsDates array and the returned bool
    is false if the date is omitted (note can't omit some components and not others)*/
    bool getObsDates(const DateTime&           sampleDate,
                     const SamplingConvention* sampleRule,
                     DateTimeArray&            obsDates,
                     DateTime&                 finalObsDate) const;

    /////////// fields /////////////
    string             name;        // asset's name
    CVolBaseWrapper    basketVol;   // explicit basket vol
    CAssetWrapperArray assets;      // array of assets
    CStringArray       ccyTreatments; // None, Struck or Prot for each asset
    string             basketYCName;  // name of basket currency
    CDoubleArray       pubWeights;  // percentage weight of each asset
    HolidayWrapper     marketHols;  // market holidays
    DateTime           baseDate;    // today - used to determine fwd starting
    DateTime           startDate;   // when basket starts
    CDoubleArray       spotsAtStart; /* spots of assets at start date. List
                                        always exists but values are not 
                                        valid if fwd starting */
    // fields for AssetHistory sampling
    // XCB sampling source/observation type is hard wired int the XCB definition
    StringArray             sources;
    ObservationTypeArray    obsTypes;
    StringArray             fxSources;
    ObservationTypeArray    fxObsTypes;

    // transient fields built from the pieces above
    ObservationSourceArray  obsSourceObjs;
    ObservationTypeArray    obsTypeObjs;

    // not in registration
    mutable CDoubleArray  intWeights;  /* do not access. Always go through // $unregistered
                                          weightsRef which returns a
                                          ref to internal weights */
    // transient fields
    string             basketCcyCode; 
};

typedef smartPtr<PercXCBWithVol> PercXCBWithVolSP;
typedef smartConstPtr<PercXCBWithVol> PercXCBWithVolConstSP;

DRLIB_END_NAMESPACE
#endif
