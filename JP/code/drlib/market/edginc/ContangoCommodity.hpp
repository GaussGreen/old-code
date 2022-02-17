//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ContangoCommodity.hpp
//
//   Description : Commodity asset where forward prices are expressed in terms of a contango curve
//
//   Author      : Andrew McCleery
//
//   Date        : 03 October 2005
//
//
//----------------------------------------------------------------------------

#ifndef _CONTANGOCOMMODITY_HPP
#define _CONTANGOCOMMODITY_HPP

#include "edginc/Commodity.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/Spot.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/Theta.hpp"
#include "edginc/SpotLevelProbability.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/ContangoRhoParallel.hpp"
#include "edginc/ContangoRhoPointwise.hpp"
#include "edginc/ForwardParallel.hpp"
#include "edginc/ForwardPointwise.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/SpotLevel.hpp"
#include "edginc/SpotShift.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/AssetHistoryContainer.hpp"
#include "edginc/ContangoTermShift.hpp"
#include "edginc/IEqVolNamePair.hpp"

DRLIB_BEGIN_NAMESPACE

// Commodity is just a marker to tell them apart from Equity & FX
class MARKET_DLL ContangoCommodity: public Commodity,
                         virtual public CAsset::IStruckable,
                         virtual public IAssetFairValue,
                         virtual public ITweakableWithRespectTo<Spot>,
                         virtual public Theta::IShift,
                         virtual public VegaSkewParallel::IShift,
                         virtual public VegaSkewPointwise::IShift,
                         virtual public SpotLevelProbability::Shift, 
                         virtual public DeltaSurface::IShift,
                         virtual public IEqVolNamePair,
                         virtual public ContangoRhoParallel::IRestorableShift,
                         virtual public ContangoRhoPointwise::IRestorableShift,
                         virtual public ForwardParallel::IRestorableShift,
                         virtual public ForwardPointwise::IRestorableShift,
                         virtual public VegaMatrix::IShift,
                         virtual public SpotLevel::Shift,
                         virtual public SpotShift::Shift,
                         virtual public VolRelativeShift::IShift,
                         virtual public ContangoTermShift::IShift {
public:
    static CClassConstSP const TYPE;

    // general object methods
    /** Validation */
    virtual void validatePop2Object();

    // IGeneralAsset methods
    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const;

    /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const;

    /** returns the spot price */
    virtual double getSpot() const;

    /** Pull out the component assets & correlations from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market);

    // IPriceAsset methods
    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    // asset methods
    /** returns the asset name */
    virtual string getName() const;

    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const;

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

    // IStruckable methods
    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        asset together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. Note that the struckAsset is indeed the struckAsset cf
        'this' which is the non struck asset */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      struckAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    // IAssetFairValue methods
    /** Returns fair value of commodity */
    virtual double fairValue() const;

    // sensitivity & scenario interfaces
 
    /** Returns name for delta */
    virtual string sensName(const Spot*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<Spot>& shift);

    /** Shifts the object using given shift (see Theta::Shift)*/
    virtual bool sensShift(Theta* shift);

    // implementation of VegaSkewParallel::IShift interface
    virtual string sensName(VegaSkewParallel* shift) const;
    virtual bool sensShift(VegaSkewParallel* shift);

    // implementation of VegaSkewPointwise::IShift interface
    virtual string sensName(VegaSkewPointwise* shift) const;
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    virtual bool sensShift(VegaSkewPointwise* shift);


    /** Returns the name of the stock/asset - used to determine
        whether to shift the object */
    virtual string sensName(SpotLevelProbability* shift) const;

    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(SpotLevelProbability* shift);

    /// implementation of DeltaSurface::IShift interface
    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

    /* for IEqVolNamePair */
    virtual bool getNamePairs(string& eqName, string& volName) const;

    /** Returns the name of the contango curve - used to determine
        whether to tweak the object */
    virtual string sensName(ContangoRhoParallel* shift)const;

    /** Shifts the object using given shift */    
    virtual bool sensShift(ContangoRhoParallel* shift);

    /** Restores the object to its original form */
    virtual void sensRestore(ContangoRhoParallel* shift);

    /** Returns the name of the contango curve - used to determine whether 
        to tweak the object */
    virtual string sensName(ContangoRhoPointwise* shift)const;

    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this contango curve */
    virtual ExpiryArrayConstSP sensExpiries(ContangoRhoPointwise* shift)const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(ContangoRhoPointwise* shift);

    /** Restores the object to its original form */
    virtual void sensRestore(ContangoRhoPointwise* shift);  

    /** Returns name identifying vol for vega matrix */
    virtual string sensName(VegaMatrix* shift) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaMatrix* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(VegaMatrix* shift);

    /** Implements SpotLevel scenario */
    /** Returns name identifying this object for SpotLevel */
    virtual string sensName(SpotLevel* shift) const;
    /** Shifts the object using given shift (see SpotLevel::Shift)*/
    virtual bool sensShift(SpotLevel* shift);

    /** Returns name identifying this object for Delta */
    virtual string sensName(SpotShift* shift) const;
    /** Shifts the object using given shift (see Delta::Shift)*/
    virtual bool sensShift(SpotShift* shift);

    /** Returns name identifying this object for Delta */
    virtual string sensName(VolRelativeShift* shift) const;
    /** Shifts the object using given shift (see Delta::Shift)*/
    virtual bool sensShift(VolRelativeShift* shift);

    /** Returns name identifying this object for ContangoTermShift */
    virtual string sensName(ContangoTermShift* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(ContangoTermShift* shift);

    virtual void getSensitiveStrikes(const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const;

    // Backs out a contango curve from the supplied forward points (offsets to spot)
    virtual CashFlowArraySP contangoRates(ExpiryArrayConstSP expiries, 
                                         ExpiryArrayConstSP outputExpiries,
                                         DoubleArrayConstSP fwdPts, 
                                         string bdc) const;

    // simon's new sensitivies - we tweak the fwd prices as opposed to the forward rates
    virtual string sensName(ForwardPointwise* shift)const;
    virtual ExpiryArrayConstSP sensExpiries(ForwardPointwise* shift) const;
    virtual bool sensShift(ForwardPointwise* shift);
    virtual void sensRestore(ForwardPointwise* shift);

    virtual string sensName(ForwardParallel* shift)const;
    virtual bool sensShift(ForwardParallel* shift);
    virtual void sensRestore(ForwardParallel* shift);

protected:
    ContangoCommodity(CClassConstSP clazz); 

private:
    friend class ContangoCommodityHelper;

    ContangoCommodity();
    ContangoCommodity(const ContangoCommodity& rhs);
    ContangoCommodity& operator=(const ContangoCommodity& rhs);

    // minimum set of collectors to aid validation
    static void acceptValueDateCollector(const ContangoCommodity*       asset, 
                                         CValueDateCollector* collector);

    static void acceptDeltaShift(const ContangoCommodity*      asset,
                                 ShiftSizeCollector* collector);

    static void acceptImntCcy(const ContangoCommodity*     asset,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const ContangoCommodity*    asset,
                              HolidayCollector* collector);

    // Gets the contango rate for date given contangoDates and contangoRates
    double contangoRate(const DateTime& date,
                                   const DateTimeArray& contangoDates,
                                   const DoubleArray& contangoRates) const;

    // Factored out tweaking functionality as it is repeated in several places
    double tweakFwdValue(const double shiftSize,
		                 const int bucket) const;

    // check to see if fwd dates have been built and cached if not build and cache
    void checkCache() const;
	
    // build and cache fwd dates
    void buildFwdDates(const DateTime&) const;

    //  fields
    DateTime           today;
    string             name;
    double             spot;
    SettlementSP       settles;

    // For expressing the forward
    ExpiryArray        expiries;    // these pertain to the spotDate rather than today (offset according to settles)
    DoubleArray        contango;
    mutable DateTimeArray      fwdDates;    // array of fwd dates;
    mutable bool       gotCache;    // have we built our cache yet

    CVolBaseWrapper    vol;

    // need a USD curve too (as that's the quoting ccy)
    YieldCurveWrapper  ccy;

    // the asset history
    AssetHistoryContainerSP history;
};

typedef smartConstPtr<ContangoCommodity> ContangoCommodityConstSP;
typedef smartPtr<ContangoCommodity> ContangoCommoditySP;

DRLIB_END_NAMESPACE
#endif
