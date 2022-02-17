//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : ProxyVol.hpp
//
//   Description : Volatility generated from a set of proxies (a composite vol)
//                 This used to be inside the Fund class but was factored out
//                 for wider use (i.e. with assets like ContangoCommodity)
//
//                  Vol surface is built from proxies as if we're doing an XCB
//                  with an optional spread applied followed by an optional skew
//                  Finally a minVol is used to floor the resultant vols
//
//   Author      : Ian S Stares
//
//   Date        : 27 September 2006
//
//
//----------------------------------------------------------------------------

#ifndef PROXY_VOL_HPP
#define PROXY_VOL_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/StartDateCollector.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/FutureExpiryCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/DeltaProxy.hpp"
#include "edginc/VegaProxyParallel.hpp"
#include "edginc/VegaProxyPointwise.hpp"
#include "edginc/VegaProxyMatrix.hpp"
#include "edginc/VegaSkewProxyParallel.hpp"
#include "edginc/VegaSkewProxyPointwise.hpp"
#include "edginc/RootTimeVegaProxy.hpp"
#include "edginc/DeltaSurfaceProxy.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/CorrelationTerm.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/VolTypeSensitiveStrikes.hpp"
#include "edginc/VolDeltaShiftSize.hpp"
#include "edginc/IPDFBoundaryProb.hpp"

DRLIB_BEGIN_NAMESPACE

class ProxyVol;
typedef smartPtr<ProxyVol> ProxyVolSP;
typedef smartConstPtr<ProxyVol> ProxyVolConstSP;

// ProxyVol of assets **
class MARKET_DLL ProxyVol: public CVolBase,
            virtual public ObjectIteration::IOverride,
            virtual public IVolatilityBS,
            virtual public IPDFCalculator,
            virtual public IRestorableWithRespectTo<VolParallel>,
            virtual public IRestorableWithRespectTo<VolPointwise>,
            virtual public VegaMatrix::IShift,
            virtual public RootTimeVega::IRestorableShift,
            virtual public VegaSkewParallel::IShift,
            virtual public VegaSkewPointwise::IShift,
            virtual public VolParallelShift::Shift,
            virtual public VolBenchmarkShift::Shift,
            virtual public PowerVega::Shift,
            virtual public ITweakableWithRespectTo<RateParallel>,
            virtual public ITweakableWithRespectTo<RatePointwise>,
            virtual public DeltaProxy::Shift,
            virtual public VegaProxyParallel::Shift,
            virtual public VegaProxyPointwise::IShift,
            virtual public VegaProxyMatrix::IShift,
            virtual public VegaSkewProxyParallel::IShift,
            virtual public VegaSkewProxyPointwise::IShift,
            virtual public RootTimeVegaProxy::IShift,
            virtual public Theta::IShift,
            virtual public IVolDeltaShiftSize,
            virtual public DeltaSurfaceProxy::IShift,
            virtual public VolRelativeShift::IShift, 
            virtual public VolAbsoluteShift::IShift,
            virtual public IPDFBoundaryProb,
            virtual public IVolTypeSensitiveStrikes {
public:
    static CClassConstSP const TYPE;

    // so Fund can build one inside itself (for now)
    ProxyVol(string&             name,
             YieldCurveWrapper&  yc,
             CAssetWrapperArray& components,
             CStringArray&       ccyTreatments,
             DoubleArray&        weights,
             HolidayWrapper&     marketHols,
             DateTime&           baseDate,
             DoubleMatrix&       correlations,
             TimeMetricSP&       timeMetric,
             ExpiryArraySP&      spreadDates,
             DoubleArray&        volSpread,
             DoubleArray&        volSkew,
             double              minVol,
             bool                correlationTerm = false);

    /** Validation */
    virtual void validatePop2Object();

    // Pull out the component assets & correlations from the market data **
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns name of vol */
    virtual string getName() const;

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    /** handles the delta shift size adjustment */
    virtual void adjustDeltaShiftSize(ShiftSizeCollector* collector,
                                      const string assetName,
                                      double spot) const;

    /** returns sensitive strikes for a given vol request and asset */
    virtual void getSensitiveStrikes(const CAsset* asset,
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const;

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                CResults*       results,
                                const DateTime& maturityDate) const;

    /** Shifts the object using given shift (see Theta::Shift)*/
    virtual bool sensShift(Theta* shift);

    /** Returns name identifying vol for vega parallel */
    virtual string sensName(const VolParallel*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak);
    /** Restores the object to its original form */
    virtual void sensRestore(const PropertyTweak<VolParallel>& tweak);

    /** Returns name identifying vol for vega pointwise */
    virtual string sensName(const VolPointwise*) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryWindowArrayConstSP sensQualifiers(const VolPointwise*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<VolPointwise>&);
    /** Restores the object to its original form */
    virtual void sensRestore(const PropertyTweak<VolPointwise>&);

    /** Returns name identifying vol for vega matrix */
    virtual string sensName(VegaMatrix* shift) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaMatrix* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(VegaMatrix* shift);

    /** Returns name identifying vol for root time vega */
    virtual string sensName(RootTimeVega* shift) const;
    /** Shifts the object using given shift for root time vega */
    virtual bool sensShift(RootTimeVega* shift);
    /** Restores the object to its original form */
    virtual void sensRestore(RootTimeVega* shift);

    /// implementation of VegaSkewParallel::IShift interface
    virtual string sensName(VegaSkewParallel* shift) const;
    virtual bool sensShift(VegaSkewParallel* shift);

    /// implementation of VegaSkewPointwise::IShift interface
    virtual string sensName(VegaSkewPointwise* shift) const;
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    virtual bool sensShift(VegaSkewPointwise* shift);

    // Implements VolParallelShift scenario 
    /** Returns name identifying this object for VolParallelShift */
    virtual string sensName(VolParallelShift* shift) const;
    /** Shifts the object using given shift (see VolParallelShift::Shift)*/
    virtual bool sensShift(VolParallelShift* shift);
  
    // Implements VolBenchmarkShift scenario 
    /** Returns name identifying this object for VolBenchmarkShift */
    virtual string sensName(VolBenchmarkShift* shift) const;
    /** Shifts the object using given shift (see VolBenchmarkShift::Shift)*/
    virtual bool sensShift(VolBenchmarkShift* shift);

    // Implements PowerVega scenario
    /** Returns name identifying this object for PowerVega */
    virtual string sensName(PowerVega* shift) const;
    /** Shifts the object using given shift (see PowerVega::Shift)*/
    virtual bool sensShift(PowerVega* shift);

    /** Returns name identifying yield curve for rho parallel */
    virtual string sensName(const RateParallel* shift) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<RateParallel>& shift);

    /** Returns the name of the yield curve - used to determine whether 
        to tweak the object */
    virtual string sensName(const RatePointwise* shift) const;   
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this  yield curve */
    virtual ExpiryWindowArrayConstSP sensQualifiers(const RatePointwise* shift) const;
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
   virtual TweakOutcome sensShift(const PropertyTweak<RatePointwise>& shift);

    /** Returns true if this object matches the supplied name with
        respect to the DeltaProxy sensitivity */
    virtual bool sensNameMatches(DeltaProxy*       shift,
                                 const OutputName& name) const;
    
    /** Appends the name(s) of this object with respect to
        the DeltaProxy sensitivity to the supplied list. */
    virtual void sensAppendName(DeltaProxy*      shift,
                                OutputNameArray& namesList) const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(DeltaProxy* shift);

    /** Returns true if this object matches the supplied name with 
        respect to the VegaProxyParallel sensitivity */
    virtual bool sensNameMatches(VegaProxyParallel* shift,
                                 const OutputName&  name) const;
    
    /** Appends the name(s) of this object with respect to
        the VegaProxyParallel sensitivity to the supplied list. */
    virtual void sensAppendName(VegaProxyParallel* shift,
                                OutputNameArray&   namesList) const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(VegaProxyParallel* shift);

    /** Returns true if this object matches the supplied name with 
        respect to the VegaProxyPointwise sensitivity */
    virtual bool sensNameMatches(VegaProxyPointwise* shift,
                                 const OutputName&   name) const;
    
    /** Appends the name(s) of this object with respect to
        the VegaProxyPointwise sensitivity to the supplied list. */
    virtual void sensAppendName(VegaProxyPointwise* shift,
                                OutputNameArray&    namesList) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaProxyPointwise* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(VegaProxyPointwise* shift);

    /** Returns true if this object matches the supplied name with 
        respect to the sensitivity */
    virtual bool sensNameMatches(VegaProxyMatrix*  shift,
                                 const OutputName& name) const;
    
    /** Appends the name(s) of this object with respect to
        the sensitivity to the supplied list. */
    virtual void sensAppendName(VegaProxyMatrix* shift,
                                OutputNameArray& namesList) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaProxyMatrix* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(VegaProxyMatrix* shift);

    /** Returns true if this object matches the supplied name with 
        respect to the sensitivity */
    virtual bool sensNameMatches(VegaSkewProxyParallel* shift,
                                 const OutputName&      name) const;
    
    /** Appends the name(s) of this object with respect to
        the sensitivity to the supplied list. */
    virtual void sensAppendName(VegaSkewProxyParallel* shift,
                                OutputNameArray&       namesList) const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(VegaSkewProxyParallel* shift);

    /** Returns true if this object matches the supplied name with 
        respect to the sensitivity */
    virtual bool sensNameMatches(VegaSkewProxyPointwise* shift,
                                 const OutputName&       name) const;
    
    /** Appends the name(s) of this object with respect to
        the sensitivity to the supplied list. */
    virtual void sensAppendName(VegaSkewProxyPointwise* shift,
                                OutputNameArray&        namesList) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewProxyPointwise* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(VegaSkewProxyPointwise* shift);

    /** Returns true if this object matches the supplied name with 
        respect to the sensitivity */
    virtual bool sensNameMatches(RootTimeVegaProxy* shift,
                                 const OutputName&  name) const;
    
    /** Appends the name(s) of this object with respect to
        the sensitivity to the supplied list. */
    virtual void sensAppendName(RootTimeVegaProxy* shift,
                                OutputNameArray&   namesList) const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(RootTimeVegaProxy* shift);

    /** implementation of DeltaSurfaceProxy::IShift interface*/
    virtual bool sensNameMatches(DeltaSurfaceProxy* shift,
                                 const OutputName&  name) const;
    
    virtual void sensAppendName(DeltaSurfaceProxy* shift,
                                OutputNameArray&   namesList) const;
    
    virtual bool sensShift(DeltaSurfaceProxy* shift);

    // implementation of VolRelativeShift::IShift interface
    virtual string sensName(VolRelativeShift* shift) const;
    virtual bool sensShift(VolRelativeShift* shift);

    // implementation of VolAbsoluteShift::IShift interface
    virtual string sensName(VolAbsoluteShift* shift) const;
    virtual bool sensShift(VolAbsoluteShift* shift);

    /** used to turn off proxy greeks */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

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

    /** Returns pdfBoundaryProb */
    virtual double getPDFBoundaryProb() const;

    // return a pdf calculator */
    virtual PDFCalculator* getPDFCalculator(const PDFRequest* request,
                                            const CAsset*     asset) const;

private:
    friend class ProxyVolHelper;
    friend class FundVol;

    ProxyVol();
    ProxyVol(const ProxyVol& rhs);
    ProxyVol& operator=(const ProxyVol& rhs);

    // calculate the vol interp scale factor for smile **
    double volInterpSmileScale(
        double assetSpot,    /* (I) the current asset spot price */
        int    i) const;    /* (I) index of asset in compAsset */

    /** Returns a reference to the internal weights to use for
        combining assets */
    const CDoubleArray& weightsRef(const CAsset* asset) const;

    VolSurfaceSP buildVegaMatrixVol(const CAsset* asset) const;

    static void acceptCollector(const ProxyVol* proxyVol, 
                                StartDateCollector* collector);

    static void acceptNameCollector(const ProxyVol* proxyVol, 
                                    AssetNameCollector* collector);

    static void acceptFutureCollector(const ProxyVol* proxyVol, 
                                      FutureExpiryCollector* collector);

    static void acceptValueDateCollector(const ProxyVol* proxyVol, 
                                         CValueDateCollector* collector);

    static void acceptImntCcy(const ProxyVol* proxyVol,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const ProxyVol* proxyVol,
                              HolidayCollector* collector);

    CVolRequestLNArraySP getComponentVolRequests(
        const CVolRequest* volRequest,
        const CAsset* asset) const;

    bool needToSkew() const;
    double measureSkew(Expiry* expiry, VolSurface* vol, const CAsset* asset) const;

    VolSurfaceSP volSkewATM(const CAsset* asset) const;
    
    VolSurfaceSP defaultVolSurface(const CAsset* asset) const;

    void applyVegaSkew(VolSurface* vol) const;

    DoubleArraySP interpLevels(const CVolRequest* volRequest,
                               const CAsset* asset) const;
    VolSurfaceSP  volSurface(const CVolRequest* volRequest,
                             const CAsset* asset) const;
    VolSurfaceSP  volSurface(const DoubleArray& strikes,
                             const CAsset* asset) const;
    DoubleArray interpVolCurve(const CVolRequest* volRequest,
                               CompositeVol*      compVol,
                               const CAsset* asset) const;

    CDoubleMatrixSP interpVolMatrix(
        const ExpiryArraySP benchmarks,
        const DoubleArray&  strikes,
        const CAsset* asset) const;  

    // Appends the market data sensNames inside supplied proxy to the passed names array
    void proxySensNames(SensControlPerName* sens, const Asset* proxy, OutputNameArray& sensNames) const;

    // Returns the expiries of the target market data present in the proxies.
    ExpiryArrayConstSP proxySensExpiries(VectorShift *sens) const;

    //  fields
    string             name;        
    YieldCurveWrapper  yc;          // yc for the asset's ccy
    CAssetWrapperArray proxies;       // array of proxy assets
    CStringArray       ccyTreatments; // None, Struck or Prot for each asset
    DoubleArray        weights;       // percentage of each asset 
    HolidayWrapper     marketHols;    // market holidays
    TimeMetricSP       timeMetric;  
    ExpiryArraySP      spreadDates; // for spread adjustment
    DoubleArray        volSpread;   // ditto
    DoubleArray        volSkew;     // for skew adjustment
    double             minVol;      // a floor to the final vol
    bool               useCorrelationTerm;    

    // optional fields
    DateTime           baseDate;
    DoubleMatrix       correlations; /* between assets - this should be 
                                        initialised in the getMarketData 
                                        method using Correlations */

    //transient fields
    mutable DoubleArray adjWeights;   /* do not access. Always go through
                                         weightsRef which returns a
                                        ref to internal weights */
    DoubleArray        sensSpread;    // for generated (i.e. top) vol level sensitivities

    bool                useVegaMatrix;
    VegaMatrixSP        vegaMatrix;
    bool                useSkewParallel;
    VegaSkewParallelSP  skewParallel;
    bool                useSkewPointwise;
    VegaSkewPointwiseSP skewPointwise;

    mutable DoubleArray pdfStrikes;    /* cached strikes for pdf calculator */

    double pdfBoundaryProb;            /* Used to limit strikes returned from CVolProcessedBS (transient) */ 

    // transient but tweakable fields
    CorrelationCommonArray        corrObjects;    // array of correlation objects -- transient
    CorrelationTermArray    corrTermArray;  // array of correlation term objects -- transient*/
};

DRLIB_END_NAMESPACE
#endif
