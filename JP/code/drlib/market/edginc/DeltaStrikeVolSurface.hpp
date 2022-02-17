//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaStrikeVolSurface.hpp
//
//   Description : A vol surface where deltas form the x axis rather than
//                 strikes - what you might expect for FX or gold
//
//   Author      : Andrew J Swain
//
//   Date        : 18 October 2004
//
//
//----------------------------------------------------------------------------

#ifndef _DELTASTRIKEVOLSURFACE_HPP
#define _DELTASTRIKEVOLSURFACE_HPP
#include "edginc/RollingSettlement.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/IPDFBoundaryProb.hpp"


//#include "edginc/VegaMatrix.hpp"
//#include "edginc/DeltaSurface.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL DeltaStrikeVolSurface: public CVolBase,
                             virtual public IPDFCalculator,
                             virtual public IVolatilityBS,
                             virtual public Theta::IShift,
                             virtual public IRestorableWithRespectTo<VolParallel>,
                             virtual public ITweakableWithRespectTo<VolPointwise>,
                             virtual public VolLevel::Shift,
                             virtual public VolParallelShift::Shift,
                             virtual public VolBenchmarkShift::Shift,
                             virtual public VolRelativeShift::IShift,
                             virtual public VolAbsoluteShift::IShift,
                             virtual public RootTimeVega::IShift,
                             virtual public IPDFBoundaryProb,
                             virtual public PowerVega::Shift {
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object();

    /** Returns name of vol */
    virtual string getName() const;

    /** Returns pdfBoundaryProb */
    virtual double getPDFBoundaryProb() const;

    // VolBase methods
    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    /* simple linear (non forward starting) interpolation at absolute strike */
    CVolProcessed* getProcessedVol(double absoluteStrike,
                                     const CAsset* asset) const;

    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        CVolBase together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    void checkCache() const;

    // Checks the passed asset has the same name and TYPE as the asset on the vol surface
    void checkAsset(const CAsset* incoming) const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    // IPDFCalculator
    virtual PDFCalculator* getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const;

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** returns the base date within the vol - needed as there is a tendency
    to use the VolSurface as a data holder into the library. Probably
    better to be able to get a new VolCurve class from a vol surface */
    const DateTime& getBaseDate() const;

    /** Returns the array of expiries in the vol suurface. See comments
        under getBaseDate. To review in conjunction with parameterised vols */
    ExpiryArrayConstSP getExpiries() const;

    TimeMetricSP getTimeMetric() const;

    /** Returns the array of expiries in the vol surface converted to
        actual dates */
    const DateTimeArray& getDates() const;

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

    /** Implements VolLevel scenario */
    /** Returns name identifying this object for VolLevel */
    virtual string sensName(VolLevel* shift) const;
    /** Shifts the object using given shift (see VolLevel::Shift)*/
    virtual bool sensShift(VolLevel* shift);

    /** Implements VolParallelShift scenario */
    /** Returns name identifying this object for VolParallelShift */
    virtual string sensName(VolParallelShift* shift) const;
    /** Shifts the object using given shift (see VolParallelShift::Shift)*/
    virtual bool sensShift(VolParallelShift* shift);

    /** Implements VolBenchmarkShift scenario */
    /** Returns name identifying this object for VolBenchmarkShift */
    virtual string sensName(VolBenchmarkShift* shift) const;
    /** Shifts the object using given shift (see VolBenchmarkShift::Shift)*/
    virtual bool sensShift(VolBenchmarkShift* shift);

    /** Implements VolRelativeShift scenario */
    /** Returns name identifying this object for VolRelativeShift */
    virtual string sensName(VolRelativeShift* shift) const;
    /** Shifts the object using given shift (see VolRelativeShift::IShift)*/
    virtual bool sensShift(VolRelativeShift* shift);

    /** Implements VolAbsoluteShift scenario */
    /** Returns name identifying this object for VolAbsoluteShift */
    virtual string sensName(VolAbsoluteShift* shift) const;
    /** Shifts the object using given shift (see VolAbsoluteShift::IShift)*/
    virtual bool sensShift(VolAbsoluteShift* shift);

    /** Returns name identifying vol for root time vega */
    virtual string sensName(RootTimeVega* shift) const;
    /** Shifts the object using given shift for root time vega */
    virtual bool sensShift(RootTimeVega* shift);

    /** Implements PowerVega scenario */
    /** Returns name identifying this object for PowerVega */
    virtual string sensName(PowerVega* shift) const;
    /** Shifts the object using given shift (see PowerVega::Shift)*/
    virtual bool sensShift(PowerVega* shift);

    /** Returns name identifying vol for vega pointwise */
    //virtual string sensName(VegaMatrix* shift) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    //virtual ExpiryArrayConstSP sensExpiries(VegaMatrix* shift) const;
    /** Shifts the object using given shift */
    //virtual bool sensShift(VegaMatrix* shift);

    // Translate the delta based surface into a fixed strike surface with the supplied
    // strikes and expiries
    VolSurfaceSP fixedStrikeSurface(const DoubleArray& strikes, 
                                    ExpiryArrayConstSP expiries, 
                                    const CAsset* asset) const;

    void skewShift(VegaSkewPointwise* shift, const CAsset* asset) const;
    void skewShift(VegaSkewParallel* shift, const CAsset* asset) const;

    VolSurfaceSP deltaSurfaceShift(DeltaSurface* shift, const CAsset* asset) const;

    double getBoundaryProb() const;

private:

    friend class FixedStrikeVol;
    friend class StrikeForDelta;
    friend class DeltaStrikeVolSurfaceHelper;
    DeltaStrikeVolSurface();
    DeltaStrikeVolSurface(
         const DeltaStrikeVolSurface*           surface,
         const DoubleArray&                     vols,
         const ExpiryArraySP&                   expiries,
         const DateTimeArraySP&                 dates,
         const DoubleArraySP&                   tradTime);

    DeltaStrikeVolSurface(
        const DeltaStrikeVolSurface*            surface,
        const DoubleMatrix&                     vol,
        const ExpiryArraySP                     expiries);

    static const double FWD_START_MIN_FWD_VOL;

    static void acceptValueDateCollector(DeltaStrikeVolSurface* volSurface,
                                         CValueDateCollector*   collector);

    void buildYearFracCache() const;
    void zapPastBM(bool useTradTime);

    //void skewShift(int expiryIndex, VegaSkewParallel* shift);
    //void skewShift(int expiryIndex, VegaSkewPointwise* shift);

    // Returns the strike corresponding to the deltas at timeIdx
    DoubleArraySP strikeForDelta(const CAsset* asset, int timeIdx) const;

    // Backs out the strike from a Black-Scholes delta
    static double strikeForDelta(double F, double T, double discFactor,
                          double sigma, double delta, bool isFwdDelta,
                          string optionType, double S);

    // Create fixed strike vol matrix with the supplied
    // strikes and expiries 
    CDoubleMatrixSP fixedStrikeVols(const DoubleArray& strikes, 
                                   ExpiryArrayConstSP expiries, 
                                   const CAsset* asset) const;

    // Return the forward price of the asset corresponding to the supplied expiry
    double fwdValue(const CAsset* asset, const DateTime& expiryDate) const;

    class Interp;
    friend class Interp;
    typedef smartPtr<Interp> InterpSP;

    // fields
    string             name;    // name of the vol
    DateTime           baseDate;
    TimeMetricSP       metric;
    ExpiryArraySP      expiries;
    BoolArray          isFwdDelta;   // spot or forward
    DoubleArray        deltas;
    CDoubleMatrixSP    vol;
    StringArray        deltaType;      // P (put), C (call), S (straddle), ATMF, ATM
    YieldCurveWrapper  discount;        // discount curve of vol surface Numeraire (for PreciousMetals this is 
                                        // the asset base/quoting currency: USD). For FX we may need to cater for the 
                                        // possibility that this is either the asset base or risk currency 
                                        // (-> may need to switch from base to risk during the interpolation)
    RollingSettlementSP  settlement;

    double pdfBoundaryProb;            /* Used to limit strikes returned from CVolProcessedBS (transient) */ 

    
    // internals
    mutable bool                 gotCache; // have we built our cache yet
    mutable DoubleArrayConstSP   tradYears; /* trading time, in years,
                                               between bm dates */
    mutable DateTimeArrayConstSP dates; /* cached - internally derived
                                           from expiries */

};

typedef smartPtr<DeltaStrikeVolSurface> DeltaStrikeVolSurfaceSP;
typedef smartConstPtr<DeltaStrikeVolSurface> DeltaStrikeVolSurfaceConstSP;
#ifndef QLIB_DELTASTRIKEVOLSURFACE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<DeltaStrikeVolSurface>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<DeltaStrikeVolSurface>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<DeltaStrikeVolSurface>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<DeltaStrikeVolSurface>);
#endif

// support for wrapper class
typedef MarketWrapper<DeltaStrikeVolSurface> DeltaStrikeVolSurfaceWrapper;
#ifndef QLIB_DELTASTRIKEVOLSURFACE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<DeltaStrikeVolSurface>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<DeltaStrikeVolSurface>);
#endif

DRLIB_END_NAMESPACE
#endif
