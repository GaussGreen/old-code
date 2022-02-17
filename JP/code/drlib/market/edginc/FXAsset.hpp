//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXAsset.hpp
//
//   Description : FX rate as an asset
//
//
//----------------------------------------------------------------------------
#ifndef FX_ASSET_HPP
#define FX_ASSET_HPP
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE
class FXAsset;
#ifndef QLIB_FXASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<FXAsset>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<FXAsset>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<FXAsset>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<FXAsset>);
#endif
DRLIB_END_NAMESPACE

#include "edginc/YieldCurve.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/Spot.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/FXDelta.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/FXRateCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/FXVolBase.hpp"
#include "edginc/SpotLevel.hpp"
#include "edginc/SpotShift.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/AssetHistoryContainer.hpp"
#include "edginc/IEqVolNamePair.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

class FXAsset;
typedef MarketWrapper<FXAsset> FXAssetWrapper;
#ifndef QLIB_FXASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<FXAsset>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<FXAsset>);
#endif

/** Implementation of CAsset for a single stock with no currency treatment */
class MARKET_DLL FXAsset: public CAsset, 
               virtual public ITweakableWithRespectTo<Spot>,
               virtual public FXDelta::RestorableShift,
               virtual public SpotLevel::Shift,
               virtual public SpotShift::Shift,
               virtual public Theta::IShift,
               virtual public ObjectIteration::IOverride,
               virtual public DeltaSurface::IShift,
               virtual public IEqVolNamePair {
public:
    static CClassConstSP const TYPE;
    friend class FXAssetHelper;
    friend class IrConverter;

    /** Pull out the two yield curves, value date and fx vol from the
        market data */
    virtual void getMarket(const IModel* model, const MarketData* market);

    static void getMarketForStruck(FXAssetWrapper&      wrapper,
                                   const IModel*        model, 
                                   const MarketData*    market);

    /** Initialises this piece of market data - records the pair of names
        idenitfying the correlation */
    virtual void initialise(MarketData* market);


    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** Returns name of (fx) vol inside asset - needed when creating a prot
        equity on the fly */
    string getVolName() const;

    /** return the FXVolBaseWrapper */
    const FXVolBaseWrapper& getFXVol()const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates the expected spot price of the asset at a given set of dates */
    virtual void fwdValue(const DateTimeArray& dateList,
                                CDoubleArray&  result) const;

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

    /** Returns 'base' ccy iso code (eg for FTSE struck into dollars,
        this would be USD) */
    string getBaseCcyIsoCode() const;
     
    /** Returns 'risk' ccy iso code (eg for FTSE struck into dollars,
        this would be GBP) */
    string getRiskCcyIsoCode() const;

    /** validation code to be called after object construction */
    virtual void validatePop2Object();

    /** Returns name identifying this object for FXDelta */
    virtual string sensName(FXDelta* shift) const;

    /** Shifts the object using given shift (see FXDelta::Shift)*/
    virtual bool sensShift(FXDelta* shift);

    /** Restores the object to its original form */
    virtual void sensRestore(FXDelta* shift);

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** Implements SpotLevel scenario */
    /** Returns name identifying this object for SpotLevel */
    virtual string sensName(SpotLevel* shift) const;
    /** Shifts the object using given shift (see SpotLevel::Shift)*/
    virtual bool sensShift(SpotLevel* shift);

    /** Returns name identifying this object */
    virtual string sensName(SpotShift* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(SpotShift* shift);

    /** Returns name for delta */
    virtual string sensName(const Spot*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<Spot>& shift);

    /// implementation of DeltaSurface::IShift interface
    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

    /** return a pdf calculator provided our vol supports it */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /* for IEqVolNamePair */
    virtual bool getNamePairs(string& eqName, string& volName) const;

    /** Accessor methods - here until we work out the general framework needed
        for SRM3 for dealing with the different assets */
    YieldCurveConstSP getBaseCcy() const; // aka domestic
    YieldCurveConstSP getRiskCcy() const; // aka foreign

    /** used to override vega style greeks - implements
        ObjectIteration::IOverride interface */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    /** Returns today's date */
    const DateTime& getToday() const;

    // Returns homogeneous greeks flag
    bool useCoherentGreeks() const;

    virtual void getSensitiveStrikes(const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const;

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

    /** returns the whole history. Used by IndexSpecFX */
    AssetHistoryContainerConstSP getAssetHistoryContainer() const;

private:
    FXAsset();
    FXAsset(const FXAsset& rhs);
    FXAsset& operator=(const FXAsset& rhs);

    static void acceptValueDateCollector(const FXAsset* asset, 
                                         CValueDateCollector* collector);

    static void acceptFXRateCollector(const FXAsset* asset, 
                                      CFXRateCollector* collector);

    static void acceptImntCcy(const FXAsset* asset,
                              AssetCcyCollector* collector);


    string            name;
    DateTime          today;       /* value date */
    HolidayWrapper    holidays;    /* NOT NEEDED */
    YieldCurveWrapper riskCcy;     /* aka foreign */
    YieldCurveWrapper baseCcy;     /* aka domestic */
    double            spotFX;      /* base ccy/risk ccy */
    FXVolBaseWrapper  fxVol;       /* cross fx vol */
    AssetHistoryContainerSP history;
        
    // derived (cached) data - needs to be cleared for Theta calculations
    DateTime     spotDate;

    /* whether to use homogeneous Greeks e.g. Delta vs FXDelta */
    bool coherentGreeks; 
};

typedef smartPtr<FXAsset> FXAssetSP;
typedef smartConstPtr<FXAsset> FXAssetConstSP;

DRLIB_END_NAMESPACE
#endif
