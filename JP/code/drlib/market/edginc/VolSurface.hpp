//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSurface.hpp
//
//   Description : Surface Based Implementation of BS Vol Interface
//
//   Author      : Mark A Robson
//
//   Date        : 21 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_VOLSURFACE_HPP
#define EDG_VOLSURFACE_HPP

#include "edginc/VolBase.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/Theta.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/AllStrikes.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/VolDeltaShiftSize.hpp"

DRLIB_BEGIN_NAMESPACE
class CVolProcessedBS;

/** Implementation of volatility where the volatility is captured by a
    volatility surface. The surface has strikes given by absolute
    amounts together with an ExpiryArray which determines the
    benchmark dates. The ExpiryArray can be made up of a combination
    of fixed and relative dates. */
class MARKET_DLL VolSurface: public CVolBase, 
                  virtual public IPDFCalculator,
                  virtual public IVolatilityBS,
                  virtual public IRestorableWithRespectTo<VolParallel>,
                  virtual public ITweakableWithRespectTo<VolPointwise>,
                  virtual public VegaMatrix::IShift,
                  virtual public RootTimeVega::IShift,
                  virtual public VolLevel::Shift,
                  virtual public VolParallelShift::Shift,
                  virtual public VolBenchmarkShift::Shift,
                  virtual public PowerVega::Shift,
                  virtual public VegaSkewParallel::IShift,
                  virtual public VegaSkewPointwise::IShift,
                  virtual public INextStrike,
                  virtual public Theta::IShift,
                  virtual public IAllStrikes,
                  virtual public IVolDeltaShiftSize,
                  virtual public DeltaSurface::IShift,
                  virtual public VolRelativeShift::IShift,
                  virtual public VolAbsoluteShift::IShift                   
{
public:
    static CClassConstSP const TYPE;
    friend class VolSurfaceHelper;
    friend class VolSurfaceAddin;

    /** Returns name of vol */
    virtual string getName() const;

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

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

    /** overrides default */
    virtual void validatePop2Object();
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

    /** Returns name identifying vol for vega pointwise */
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
    
    /// implementation of VegaSkewParallel::IShift interface
    virtual string sensName(VegaSkewParallel* shift) const;
    virtual bool sensShift(VegaSkewParallel* shift);

    /// implementation of VegaSkewPointwise::IShift interface
    virtual string sensName(VegaSkewPointwise* shift) const;
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    virtual bool sensShift(VegaSkewPointwise* shift);

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

    /** Implements PowerVega scenario */
    /** Returns name identifying this object for PowerVega */
    virtual string sensName(PowerVega* shift) const;
    /** Shifts the object using given shift (see PowerVega::Shift)*/
    virtual bool sensShift(PowerVega* shift);

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /// implementation of DeltaSurface::IShift interface
    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

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

    /** Constructor - in general though object is created via data dict.
        Note that the metric must already be populated with its market
        data (ie its holidays) */
    VolSurface(const string&        volName,
               const TimeMetric*    metric,
               const DoubleArray&   strikes,
               const DoubleMatrix&  vol,
               const ExpiryArray*   expiries,
               const DateTime&      baseDate);

    //// Like above constructor but takes variances rather than vols
    static VolSurface* createFromVars(const string&        volName,
                                      const TimeMetric*    metric,
                                      const DoubleArray&   strikes,
                                      const DoubleMatrix&  vars,
                                      const ExpiryArray*   expiries,
                                      const DateTime&      baseDate);

    /** Creates a vol surface from an existing surface but uses the
        supplid strikes and vols - this is optimised for performance */
    VolSurface(const VolSurface*    surface,
               const CDoubleArray&  strikes,
               const CDoubleMatrix& matrix);

    /** Updates this vol surface with the supplied strikes and vols. Note
        that no copy is made of the vols. (Also note that the vol surface
        may at any point overwrite or alter them.) The return value is the
        new VolSurface to use - it may or may not be the same as 'this'.
        Use smart pointers to handle memory management.
        If validate is true then the strikes are validated to be increasing
        and distinct and the vols are validated to be positive. Only set
        to this false if you're confident that the data will always be
        good */
    VolSurface* update(const CDoubleArray&    strikes,
                       const CDoubleMatrixSP& vols,
                       bool                   validate);

    /* simple linear (non forward starting) interpolation at absolute strike */
    CVolProcessedBS* getProcessedVol(double absoluteStrike) const;

    /** Same as getProcessedVol(CVolRequest, CAsset) but takes smart pointer
        to surface and solves memory ownership problems */
    static CVolProcessed* getProcessedVol(
    const smartConstPtr<VolSurface>&  volSurface,
    const CVolRequest*                volRequest,
    const CAsset*                     asset);

    /** given a current spot level, get the next strike on the vol surface */
    virtual double getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const;

    /** returns a double list of all strikes on the vol surface. */
    virtual DoubleArraySP getAllStrikes() const;

    /** handles the delta shift size adjustment */
    virtual void adjustDeltaShiftSize(ShiftSizeCollector* collector,
                                      const string assetName,
                                      double spot) const;

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

    /** Returns the array of strikes in the vol surface */
    const DoubleArray& getStrikes() const;

    /** Returns the vol matrix */
    const DoubleMatrix& getVolMatrix() const;

    /** minimum fwd vol for fwd starting interpolation */
    static const double FWD_START_MIN_FWD_VOL;


	/** minimum fwd variance for fwd starting interpolation, set to 2 bps */
	static const double FWD_START_MIN_FWD_VARIANCE;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual PDFCalculator* getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const;
   
    static const string ABSOLUTE;
    static const string SPOTMONEYNESS;
    static const string FWDMONEYNESS;

    ~VolSurface();

private:
    static void acceptValueDateCollector(const VolSurface*    volSurface,
                                         CValueDateCollector* collector);

    VolSurface();
    VolSurface(const VolSurface &rhs);
    VolSurface& operator=(const VolSurface& rhs);
    void skewShift(int expiryIndex, VegaSkewParallel* shift);
    void skewShift(int expiryIndex, VegaSkewPointwise* shift);
    void buildYearFracCache() const;
    void safeFwdVol();
    void checkCache() const;
    void zapPastBM(bool useTradTime);

    class Interp;
    friend class Interp;
    typedef smartPtr<Interp> InterpSP;

    /** Builds vol surface using existing values from surface but with
        new single strike, vols and bm dates. Note no clones of incoming
        data made */
    VolSurface(const VolSurface*         surface,
               double                    strike,
               const DoubleArray&        vols,
               const ExpiryArraySP&      expiries,
               const DateTimeArraySP&    dates,     // must match expiries
               const DoubleArraySP&      tradTime); // must match dates

    // fields ////
    string             name;    // name of the vol
    TimeMetricSP       metric;  /* ideally would be TimeMetricConstSP
                                   but call to getMarket makes it tricky */
    DoubleArray        strikes;
    CDoubleMatrixSP    vol;       
    ExpiryArraySP      expiries;
    DateTime           baseDate;

    /* transient fields (won't appear in dd interface) */
    mutable bool                 gotCache; // have we built our cache yet
    mutable DoubleArrayConstSP   tradYears; /* trading time, in years, 
                                               between bm dates */
    mutable DateTimeArrayConstSP dates; /* cached - internally derived
                                           from expiries */
    // MAR: not used currently - decision pending
    mutable InterpSP             procVol; // cached copy of last one created $unregistered

    string             strikeMode;  // Strike Representation: Absolute, FwdMoneyness etc.
};

typedef smartConstPtr<VolSurface> VolSurfaceConstSP;
typedef smartPtr<VolSurface> VolSurfaceSP;
#ifndef QLIB_VOLSURFACE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<VolSurface>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<VolSurface>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<VolSurface>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<VolSurface>);
#endif

DRLIB_END_NAMESPACE
#endif

