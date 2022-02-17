//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XCB.hpp
//
//   Description : Implement asset for percentage or unit weighted cross 
//                 currency basket with derived composite volatility
//
//   Author      : Mark A Robson
//
//   Date        : 12 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_XCB_HPP
#define EDG_XCB_HPP
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/Delta.hpp"
#include "edginc/StartDateCollector.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/FutureExpiryCollector.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/CorrelationSkew.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/CorrelationTerm.hpp"
#include "edginc/SpotShift.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/BasketSpot.hpp"
#include "edginc/IPDFBoundaryProb.hpp"

DRLIB_BEGIN_NAMESPACE
class CVolRequestDVF;

/* Separate implementation from other xcb's for clarity */

/** Implementation of CAsset for percentage weighted cross currency
    basket with vol override. */
class MARKET_DLL XCB: public CAsset,
           virtual public IAssetFairValue,
           virtual public ITweakableWithRespectTo<BasketSpot>,
           virtual public Theta::IShift,
           virtual public DeltaSurface::IShift,
           virtual public IPDFBoundaryProb,
           virtual public SpotShift::Shift {
public:
    static CClassConstSP const TYPE;
    friend class XCBHelper;
    friend class XLGetCorrelationMatrixAddin;

    static const string NO_SMILE;
    static const string FIXED_STRIKE_SMILE;
    static const string FLOAT_STRIKE_SMILE;

    ~XCB();

    /** Validation */
    virtual void validatePop2Object();

    /** Pull out the component assets & correlations from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market);

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
                             const SamplingConvention*  sampleRule) const;

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

   /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;

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

    /** Returns the name of the XCB for (Basket-)Delta */
    virtual string sensName(const BasketSpot*) const;

    /** Shifts all components in the XCB for Delta */
    virtual TweakOutcome sensShift(const PropertyTweak<BasketSpot>& shift);

    /** Shifts the object using given shift (see Theta::Shift)*/
    virtual bool sensShift(Theta* shift);

    /// implementation of DeltaSurface::IShift interface
    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

    /** Returns fair value of XCB */
    virtual double fairValue() const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    virtual void getCorrelationSkewParameters(const IModel* model, const MarketData* market);

    /** Returns name identifying this object for SpotShift */
    virtual string sensName(SpotShift* shift) const;
    /** Shifts the object using given shift (see SpotShift::Shift)*/
    virtual bool sensShift(SpotShift* shift);

    /** Can this asset physically settle? */
    virtual bool canPhysicallySettle() const;

    /** Returns pdfBoundaryProb */
    virtual double getPDFBoundaryProb() const;

private:
    XCB();
    XCB(const XCB& rhs);
    XCB& operator=(const XCB& rhs);

    /** Do we use smile when interpolating the vol */
    bool useSmileVolInterp() const;

    /** Returns true if basket needs initial spots */
    bool needInitialSpots() const;

    /** calculate the vol interp scale factor for smile */
    double volInterpSmileScale(
        double        basketSpot,    /* (I) the current basket spot price */
        int           i) const;      /* (I) index of asset in compAsset */

    /** Returns a reference to the internal weights to use for
        combining assets */
    const CDoubleArray& weightsRef() const;

    static void acceptStartDateCollector(const XCB*                asset, 
                                         StartDateCollector* collector);

    static void acceptNameCollector(const XCB* asset, 
                                    AssetNameCollector* collector);

    static void acceptFutureCollector(const XCB* asset, 
                                      FutureExpiryCollector* collector);

    static void acceptValueDateCollector(const XCB* asset, 
                                         CValueDateCollector* collector);

    static void acceptDeltaShift(const XCB* asset,
                                 ShiftSizeCollector* collector);

    static void acceptHoliday(const XCB* asset,
                              HolidayCollector* collector);

    static void acceptImntCcy(const XCB* asset,
                              AssetCcyCollector* collector);

    static void acceptWrapperNameCollector(const XCB* asset,
                                           WrapperNameCollector* collector);

    static void acceptCriticalDateCollector(const XCB* asset, 
                                            CriticalDateCollector* collector);

    void alterDeltaShiftSize(ShiftSizeCollector* collector,
                             const Delta*        sensControl) const;

    CVolRequestLNArraySP getComponentVolRequests(
        const CVolRequestLN* volRequest) const;

    double calculateMaxShiftSize(int tweakedAssetIndex,
                                 int i,
                                 const double& strike,
                                 const double& volSurfaceStrike,
                                 const double& weight,
                                 const double& assetSpot,
                                 const double& basketSpot) const;

    VolSurfaceSP defaultVolSurface() const;
    void computeCorrelationSkewAndPower(double&  correlationSkew,
                                        double&  correlationSkewPower) const;

    CVolProcessed * getProcessedVolLN(const CVolRequestLN* volRequest) const;
    CVolProcessed * getProcessedVolDVF(const CVolRequestDVF* volRequest) const;

    /** Helper for centralised sampling/isda adjustment for XCBs
    gets a list of dates for which the given sample will finally be known for 
    each component. Results are put into obsDates array and the returned bool
    is false if the date is omitted (note can't omit some components and not others)*/
    bool getObsDates(const DateTime&           sampleDate,
                     const SamplingConvention* sampleRule,
                     DateTimeArray&            obsDates,
                     DateTime&                 finalObsDate) const;

    // smile types
    static int const noSmile;     /* (C) no smile - all at the money */
    static int const fixedStrike; /* (D) fixed strike smile */
    static int const floatStrike; /* (E) floating strike smile */

    //  fields
    string             name;        // asset's name
    CAssetWrapperArray assets;      // array of assets
    CStringArray       ccyTreatments; // None, Struck or Prot for each asset
    string             basketYCName;  // name of basket currency
    bool               unitWeights; // true: unit weighted, false: % weighted
    CDoubleArray       pubWeights;  /* number of units or percentage of
                                       each asset (as initially supplied) */
    HolidayWrapper     marketHols;  // market holidays
    DateTime           baseDate;    // today - used to determine fwd starting
    DateTime           startDate;   // when basket starts
    CDoubleArray       spotsAtStart; /* spots of assets at start date. List
                                        always exists but values are not 
                                        valid if fwd starting */
    string             smileType;   // type of smile
    DoubleMatrix       correlations; /* between assets - this should be 
                                        initialised in the getMarketData 
                                        method using  Correlations */
    TimeMetricSP       timeMetric;  // for composite vol
    CorrelationSkewWrapperArray 
                       corrSkewNames; // names identifying corr skew
    CDoubleArray       corrSkewWeights; // weights for the above

    // use term structure of correlation when computing composite vol
    bool                    useCorrelationTerm;    

    // transient but tweakable fields
    CorrelationCommonArray  corrObjects;    // array of correlation objects -- transient
    CorrelationTermArray    corrTermArray;  // array of correlation term objects -- transient
    
    // fields for AssetHistory sampling
    // XCB sampling source/observation type is hard wired int the XCB definition
    StringArray             sources;
    ObservationTypeArray    obsTypes;
    StringArray             fxSources;
    ObservationTypeArray    fxObsTypes;

    // transient fields built from the pieces above
    ObservationSourceArray  obsSourceObjs;
    ObservationTypeArray    obsTypeObjs;

    int              smile;         // smile */
    mutable double   basketAtStart; /* derived basket spot at start - only
                                       valid if needInitialSpots() is true */
    mutable CDoubleArray  intWeights;  /* do not access. Always go through
                                          weightsRef which returns a
                                          ref to internal weights */
    string              basketCcyCode; 
    mutable DoubleArray pdfStrikes;    /* cached strikes for pdf calculator/
                                          local vol*/
    double pdfBoundaryProb;            /* Used to limit strikes returned from CVolProcessedBS (transient) */ 
    // not used, but too painful to remove
    CorrelationTermArray    corrTermArrayNotUsed; // $unregistered
};

typedef smartPtr<XCB> XCBSP;
typedef smartConstPtr<XCB> XCBConstSP;

DRLIB_END_NAMESPACE
#endif
