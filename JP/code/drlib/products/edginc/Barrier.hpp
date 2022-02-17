//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Barrier.hpp
//
//   Description : Methods that a Barrier classes must implement
//
//   Date        : Feb 03
//
//
//----------------------------------------------------------------------------

#ifndef EDR_BARRIER_HPP
#define EDR_BARRIER_HPP

#include "edginc/Schedule.hpp"
#include "edginc/BarrierLevel.hpp"
#include "edginc/SVGenBarrier.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/BarrierUtil.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

class IBarrierBendMaker;
class IBarrierUtil;
typedef refCountPtr<IBarrierUtil> IBarrierUtilSP;
class IBarrierMgr;
class SVGenBarrierHV;
typedef refCountPtr<SVGenBarrierHV> SVGenBarrierHVSP;

// Barrier Base Class
class PRODUCTS_DLL Barrier : public CObject,
                virtual public LegalTerms::Shift,
                virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE;
    static const string EUROPEAN_MONITORING;
    static const string DAILY_MONITORING;
    static const string CONTINUOUS_MONITORING;
    static const double MINIMUM_SPREAD;
    static const double ADJUST_CONSTANT;
    static const double BUS_DAYS_IN_YEAR;
    static const double SQRT_ONE_DAY_YEAR_FRAC;
    static const double BUS_DAYS_IN_YEAR_BA;
    static const double SQRT_ONE_DAY_YEAR_FRAC_BA;

    /** Conversion method for state variables */
    /** Creates a SVGenBarrierHV from a set of data */
    virtual SVGenBarrierHVSP convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset) const = 0;

    /** Creates a SVGenBarrierHVT from a set of data */
    virtual SVGenBarrierHVTSP convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset, int iDummy) const = 0;

    // barrier adjustement for continuous monitoring
    static void BarrierAdjustment(double vol, bool isUp, double& barrierLevel);

    // for barrier (tree) products....
    /**  set up barrier levels for each step */
    // return the barrier map array (See DateTime::createMapping)
    static IntArray setStepLevel(ScheduleSP in, double scale, const DateTimeArray& stepDates, CDoubleArray& out);

    virtual ~Barrier();

    /* Methods a barrier must implement */
    virtual void validate(const MonteCarlo* model, 
                          int nbAssets) = 0;

    virtual void validate(int nbAssets)const = 0;

    virtual void amendIsHit(const IMultiFactors* mFactors,
                              IRefLevel::IStateVarSP refLevel,
                              const DateTimeArrayArray& monitoringDates,
                              const DateTime& today) = 0;

    // Create the interpolated barrier levels from sample dates
    virtual void createInterpBarrier(const DateTime& valueDate,
                                     const DateTimeArray& samples)const = 0;

    virtual void overrideEventBarrier(const DateTime& valueDate)const = 0;

    // Allow barrier to influence mon dates in inst if Daily
    virtual DateTimeArray getFutureMonitoringDates(const DateTime&         today,
                                                   const DateTimeArray&    instMonitorDates, 
                                                   smartPtr<IMCPathConfig> pathConfig) = 0;

    virtual bool createBarrierUtil(const DateTime& valueDate) const = 0;

    virtual refCountPtr<IBarrierUtil> getBarrierUtil() const = 0;

    virtual void pathUpdated(const IMCPathGenerator* pathGen) = 0;

    //! Obtain the schedule for the barrier
    /** This assumes that the barrier contains the dates and levels and
        interpolation style and is able to comminicate that via a 
        schedule */
    virtual ScheduleSP getBarrierSchedule() const = 0;

    // Get the interp level from the barrier
    virtual double getInterpLevel(const IMCPathGenerator* pathGen,
                                  const IMCProduct* product,
                                  int iAsset)const = 0;
                                  
    // For barrier adjustment N-Factor
    virtual void adjustLN(const IMCPathGenerator* pathGen,
                          const IMCProduct*        product,
                          const IMultiFactors*    assets,
                          int                     pathIdx)= 0;

    // For barrier adjustment N-Factor
    virtual void adjustImplied(const IMCPathGenerator*  pathGen,
                               const IMCProduct*         product,
                               const IMultiFactors*     assets,
                               const YieldCurve*        discount,
                               InstrumentSettlement*    instSettle,
                               int                      pathIdx)= 0;

    /** Another way to look at things... */
    // Still not terribly clear how the per-asset is to operate : should we have 
    // barrier class instance per asset, or one barrier which deals with several assets?
    // For now have 2 - one which returns values purely for a single asset
    // - and one which gives an "overall" value for all assets.
    // Don't make them const - may allow "on-demand" type operation?
    virtual double hitValue(const  IMCPathGenerator*  pathGen,
                            int    assetIdx) = 0;
    virtual double hitValue(const  IMCPathGenerator*  pathGen) = 0;

    // When does any hit condition get met?
    // Return is MC-specific and in context of the provided pathGen : the stepIdx at which 
    // hit occured. Supported for 'E' type only
    virtual void hitValueAndTime(const   IMCPathGenerator*  pathGen,
                                 int     assetIdx,
                                 double& value,
                                 bool&   isConditionMet,
                                 int&    metAtStep) = 0;
    virtual void hitValueAndTime(const  IMCPathGenerator*  pathGen,
                                 double& value,
                                 bool&   isConditionMet,
                                 int&    metAtStep) = 0;

    virtual void hitValueAndTimeAsBasket(const  IMCPathGenerator*  pathGen,
                                         const  IAggregateMakerSP& bsk,
                                         double& value,
                                         bool&   isConditionMet,
                                         int&    metAtStep) = 0;

    // 3 routines, somewhat hacky for multi-period monitoring. Pending redesign with StateVars
    virtual void hitValueAndTimeAtStep(const IMCPathGenerator*  pathGen,
                                       int                      iStep,
                                       int                      iAsset,
                                       double                   refLevel,
                                       double&                  value,           
                                       bool&                    isConditionMet,
                                       int&                     metAtStep) = 0;
    virtual void hitValueAndTimeAtMat(const IMCPathGenerator*  pathGen,
                                      int                      iStep,
                                      int                      iAsset,
                                      double                   refLevel,
                                      bool                     isConditionMetSoFar,
                                      double&                  value,           
                                      bool&                    isConditionMet,
                                      int&                     metAtStep) = 0;
    virtual void multiHelper(DoubleArray&       hitValues,
                             const BoolArray&   isConditionMet,
                             const IntArray&    metAtStep,
                             double&            overallHitValue,
                             bool&              overallIsConditionMet,
                             int&               overallMetAtStep) = 0;

    virtual void multiHelper2(DoubleArray&       hitValues,
                              const BoolArray&   isConditionMet,
                              double&            overallHitValue,
                              double&            overallBinaryValue) = 0;
    
    // Help client code report barrier levels
    virtual BarrierLevelArraySP reportLevels(const DateTime& valueDate,
                                             double          refLevel,
                                             int             assetIdx) const = 0;

//    virtual BarrierLevelArraySP reportLevelsAsBasket(const  IAggregateMakerSP& bsk,
//                                                     const DateTime& valueDate,
//                                                     const DoubleArray assetPerf,
//                                                     double          refLevel,
//                                                     int             assetIdx) const = 0;

	//// roll through time replacing risk barrier in past by legal one
	virtual bool sensShift(Theta* theta) = 0;

	/** Satisfy LegalTerms::Shift interface */
    virtual bool sensShift(LegalTerms* shift) = 0;

    // gets barrier breach events for last date in pathGen
    virtual void getEvents(const  IMCPathGenerator*  pathGen, 
                            EventResults* events,
                            bool useIsHitFlags,
                            const string barrierName,
                            StringArray* assNames) = 0;

    virtual void getEvents(EventResults* events,
                            const string barrierName,
                            int   thisPeriodIdx,
                            const IntArrayArray& metAtStep,
                            const StringArray& assNames,
                            const DoubleMatrix& assLevels,
                            const DoubleMatrix& refLevels) = 0;

    virtual bool isMonitoringPoint(const DateTimeArray& monitoringDates,
                                    const DateTime& today) = 0;

protected:

    // Calculate fwd vols between monitor  dates - for use in closed form barrier adjustment 
    static void calculateFwdVols(const  DoubleArray&   interpBarrier,
                          const  CAsset&        asset,
                          double                refLevel,
                          const  DateTime&      startDate, 
                          const  DateTimeArray& monDates,
                          const  DateTime&      today, 
                          DoubleArray&          vols,
                          DoubleArray&          days,
                          DoubleArray&          spotVols,
                          DoubleArray&          spotDays);

    struct PRODUCTS_DLL DigitalDiffFunc{
        DateTime          valueDate;
        DateTime          startDate;
        DateTime          matDate;
        bool              isCall;
        bool              fwdStarting;
        bool              oneContract;
        double            notional;
        mutable double    highStrike;
        double            initialSpot;
        const InstrumentSettlement* instSettle;
        const CAsset*     asset;
        const YieldCurve* discount;
        double            target;

        double operator()(double lStrike) const;

    };

    /* Find the adjusted barrier level that satisfies :
       N(Xadj) = Digital(Badj(t))
       where Xadj is a closed form adjusted random variable and B(t) is the unadjusted barrier level */
    static double digitalRoot(DigitalDiffFunc& digital,
                       double initGuess);

//    class BskPerfFunc;

protected:
    friend class BarrierHelper;
    Barrier(); // not implemented
    Barrier(CClassConstSP clazz);
private:
    Barrier(const Barrier& rhs);
    Barrier& operator=(const Barrier& rhs);
};

typedef smartConstPtr<Barrier> BarrierConstSP;
typedef smartPtr<Barrier> BarrierSP;
#ifndef QLIB_BARRIER_CPP
EXTERN_TEMPLATE(class PRODUCTS_DLL_SP smartPtr<Barrier>);
#else
INSTANTIATE_TEMPLATE(class PRODUCTS_DLL smartPtr<Barrier>);
#endif

// The simplest barrier class
class PRODUCTS_DLL BarrierSchedule : public Barrier{
public:
    
    /** Conversion to state variable representation */
    virtual SVGenBarrierHVSP convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset) const;

    /** Creates a SVGenBarrierHVT from a set of data */
    virtual SVGenBarrierHVTSP convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset, int iDummy) const;

    static CClassConstSP const TYPE;
    friend class BarrierScheduleHelper;

    // constructor to allow other product class to have own internal Barrier Class.
    BarrierSchedule(ScheduleSP   levels,
                                 bool         isOut,
                                 bool         isUp,
                                 ScheduleSP   economicLevels,
                                 int          numHits,
                                 IntArraySP   isHit,
                                 DateTime     hitDate,
                                 string       monitorType,
                                 bool         smoothing,
                                 double       lowSpread,
                                 double       highSpread,
                                 bool         isSimultaneousBreach = false,
                                 bool         isInternalMonitoringDates = false);

    virtual ~BarrierSchedule();

    virtual void validate(const MonteCarlo* model, 
                          int nbAssets);
    
    void validate(int nbAssets)const;

    string monTypeString() const;

    bool isKO() const { return isOut; }

    virtual void amendIsHit(const IMultiFactors* mFactors,
                              IRefLevel::IStateVarSP refLevel,
                              const DateTimeArrayArray& monitoringDates,
                              const DateTime& today);

    virtual void createInterpBarrier(const DateTime& valueDate,
                                     const DateTimeArray& samples)const;

    virtual void overrideEventBarrier(const DateTime& valueDate)const;

    // Allow barrier to influence mon dates in inst if Daily
    virtual DateTimeArray getFutureMonitoringDates(const DateTime&         today,
                                                   const DateTimeArray&    instMonitorDates, 
                                                   smartPtr<IMCPathConfig> pathConfig);

    // instance of IBarrierPay creation
    virtual bool createBarrierUtil(const DateTime& valueDate) const;

    virtual refCountPtr<IBarrierUtil> getBarrierUtil() const;

    virtual void pathUpdated(const IMCPathGenerator* pathGen);

    virtual ScheduleSP getBarrierSchedule() const;

    virtual double getInterpLevel(const IMCPathGenerator* pathGen,
                                  const IMCProduct* product,
                                  int iAsset)const;

    virtual void adjustLN(const IMCPathGenerator*  pathGen,
                          const IMCProduct*         product,
                          const IMultiFactors*     assets,
                          int                      pathIdx);

    virtual void adjustImplied(const IMCPathGenerator* pathGen,
                               const IMCProduct*        product,
                               const IMultiFactors*    assets,
                               const YieldCurve*       discount,
                               InstrumentSettlement*   instSettle,
                               int                     pathIdx);

    virtual double hitValue(const  IMCPathGenerator*  pathGen,
                            int    assetIdx);
    virtual double hitValue(const  IMCPathGenerator*  pathGen);

    virtual void hitValueAndTime(const   IMCPathGenerator*  pathGen,
                                 int     assetIdx,
                                 double& value,
                                 bool&   isHit,
                                 int&    hitTime);
    virtual void hitValueAndTime(const  IMCPathGenerator*  pathGen,
                                 double& value,
                                 bool&   isHit,
                                 int&    hitTime);

    virtual void hitValueAndTimeAsBasket(const  IMCPathGenerator*  pathGen,
                                         const  IAggregateMakerSP& bsk,
                                         double& value,
                                         bool&   isConditionMet,
                                         int&    metAtStep);

    // Several routines, somewhat hacky for multi-period monitoring. Pending redesign with StateVars
    virtual void hitValueAndTimeAtStep(const IMCPathGenerator*  pathGen,
                                       int                      iStep,
                                       int                      iAsset,
                                       double                   refLevel,
                                       double&                  value,           
                                       bool&                    isConditionMet,
                                       int&                     metAtStep);
    virtual void hitValueAndTimeAtMat(const IMCPathGenerator*  pathGen,
                                      int                      iStep,
                                      int                      iAsset,
                                      double                   refLevel,
                                      bool                     isConditionMetSoFar,
                                      double&                  value,           
                                      bool&                    isConditionMet,
                                      int&                     metAtStep);

    virtual void multiHelper(DoubleArray&       hitValues,
                             const BoolArray&   isConditionMet,
                             const IntArray&    metAtStep,
                             double&            overallHitValue,
                             bool&              overallIsConditionMet,
                             int&               overallMetAtStep);

    virtual void multiHelper2(DoubleArray&       hitValues,
                              const BoolArray&   isConditionMet,
                              double&            overallHitValue,
                              double&            overallBinaryValue); // unsmoothed

    virtual BarrierLevelArraySP reportLevels(const DateTime& valueDate,
                                             double          refLevel,
                                             int             assetIdx) const;

//    virtual BarrierLevelArraySP reportLevelsAsBasket(const  IAggregateMakerSP& bsk,
//                                                     const DateTime& valueDate,
//                                                     const DoubleArray assetPerf,
//                                                     double          refLevel,
//                                                     int             assetIdx) const;

	//// roll through time replacing risk barrier in past by legal one
	virtual bool sensShift(Theta* theta);

    /** Satisfy LegalTerms::Shift interface */
    virtual bool sensShift(LegalTerms* shift);

    virtual void getEvents(const  IMCPathGenerator*  pathGen, 
                            EventResults* events,
                            bool useIsHitFlags,
                            const string barrierName,
                            StringArray* assNames);

    // gets barrier breach events by looking at state after MC has run past
    virtual void getEvents(EventResults* events,
                            const string barrierName,
                            int   thisPeriodIdx,
                            const IntArrayArray& metAtStep,
                            const StringArray& assNames,
                            const DoubleMatrix& assLevels,
                            const DoubleMatrix& refLevels);

    // gets barrier breach event for RRK and (possibly) TriggerECO, ie assets where a basket
    // is monitored for breach
    void getEvents(const IMCPathGenerator* pathGen,
        EventResults* events,
        const string barrierName,
        const IAggregateMakerSP& bsk);

    virtual bool isMonitoringPoint(const DateTimeArray& monitoringDates,
                                    const DateTime& today);

private:
    BarrierSchedule();
    BarrierSchedule(const BarrierSchedule& rhs);
    BarrierSchedule& operator=(const BarrierSchedule& rhs);
    virtual void validatePop2Object();

protected:
    BarrierSchedule(CClassConstSP clazz);

private:
    int                    numHits;           // How many hits trigger a breach
    mutable IntArraySP     isHit;             // has the barrier already been breached for each asset
    DateTime               hitDate;           // when it was breached
    string                 monitorType;       /* "E" - European : monitor on given dates only
                                                 "D" - Daily monitoring between given dates
                                                 "C" - Continuous monitoring between given dates */    
    // For smoothing
    bool                   smoothing;         // Is smoothing to be used ?
    double                 lowSpread;
    double                 highSpread;

    ScheduleSP                  levels;            // the schedule of dates and levels
    bool                        isOut;
    bool                        isUp;
    bool                        isSimultaneousBreach;
    string                      closedFormDates;   //!< Specifies which dates to use for closed form
    string                      closedFormMethod;  //!< Specifies the type of closed form to be used
    string                      maxNumDivsPerYear; //!< Maximum number of divs per year to be used in BB: either numeric or "All"

    // Optional
    ScheduleSP                  economicLevels;
    smartPtr<IBarrierBendMaker> bendMaker;

    // Transient
    mutable int            numAssets;         // already have from isHit->size() but convenient
    mutable DoubleArraySP  smoothHitValues;   // avoid alloc lower down
    mutable CBoolArraySP   isHitPerAsset;     // working area recording hit/not for each asset
    mutable CBoolArraySP   isHitPerAssetSoFar; // past storage version
    mutable bool           useAdjBarrier;     // use the adjusted Barrier levels
    mutable bool           forceDailySteps;   // for LV - insert extra dates to allow "D"/"C" monitoring $unregistered

    // for multi-asset case ... I hate this design!
    mutable CIntArray      pastHitTime;       // may have had assets hit in the past
    mutable int            numRealHitsSoFar;      // past of count of hits (for multi-asset case)
    mutable bool           isConditionMetSoFar; // for any overall condition
    mutable int            globalPastHitTime; // if overall hit condition met in past - this is the step when. 
                                              // Goodness knows how we'd get a DateTime -> hitDate!?

    // Some transient fields for state variables
    CBoolArraySP           amendedIsHit;          //!< IsHit && current spot <> current barrier
    bool                   doneAmendedIsHit;      //!< Whether amendedIsHit has been set or not

        // Transient
    mutable DateTimeArraySP     barrierDates;      // actually all the sim dates
    mutable DoubleArraySP       interpBarrier;     // interpolated barrier to match barrier dates
    mutable DateTime			origDate;		// value date when interpolated barrier was built
    mutable DoubleArrayArraySP  adjBarLevels;      // closed form adjusted levels per asset
    mutable IntArray            areMet;            // useful for multiHelper

    mutable refCountPtr<IBarrierUtil> iBarrier;          // cache the IBarrier $unregistered
    mutable vector<IBarrierUtilSP> ibsPerAsset;  // cache the IBarriers per asset, used for barrier adjustmentz $unregistered

    // following three variable and related functionality maybe better to move to BarrierUtil.
    mutable bool                isInternalMonitoringDates;  // true to use instrument's monitor date  $unregistered
                                                            // rather than barrier class's monitor.
    mutable vector<bool>        isMonitorStep;              // is barrier minitor step? $unregistered
    mutable IntArray            barMap; // $unregistered

    /** For sort algorithm - allows an int array to be sorted based upon
    the value of a corresponding int array */
    class PRODUCTS_DLL IntSortHelper{
        IntArray ints;
    public:
        IntSortHelper(IntArray ints);

        //// we want best first (ie largest to smallest)
        bool operator()(int i1, int i2);
    };
};


typedef smartConstPtr<BarrierSchedule> BarrierScheduleConstSP;
typedef smartPtr<BarrierSchedule> BarrierScheduleSP;


/******************************************************************/
// DoubleBarrier class
/******************************************************************/

class PRODUCTS_DLL DoubleBarrier : public Barrier{    
private:
    BarrierSP       barrierA;
    BarrierSP       barrierB;
    // isAnd is used to indicate:
    //      if A and B same type, one KI/O or need both(not simultaneous) to KI/KO 
    //      if A and B are KI and KO respectively, false=KI_CANCEL_KO, true=KI_KEEP_KO
    bool            isAnd;
    bool            isSimultaneousBreach;

    // transient fields
    mutable refCountPtr<IBarrierUtil>   iBarrierA; // $unregistered
    mutable refCountPtr<IBarrierUtil>   iBarrierB; // $unregistered
    mutable refCountPtr<IBarrierUtil>   iBarrier; // $unregistered

public:
    static CClassConstSP const TYPE;

    DoubleBarrier(bool isAnd,
        bool           isSimultaneousBreach,
        BarrierSP      barrierA, 
        BarrierSP      barrierB);
    ~DoubleBarrier();

    virtual SVGenBarrierHVSP convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset) const;

    /** Creates a SVGenBarrierHVT from a set of data */
    virtual SVGenBarrierHVTSP convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset, int iDummy) const;
    virtual void validatePop2Object();

    virtual void validate(const MonteCarlo* model, 
                          int nbAssets);

    virtual void validate(int nbAssets) const;

    virtual void amendIsHit(const IMultiFactors* mFactors,
                              IRefLevel::IStateVarSP refLevel,
                              const DateTimeArrayArray& monitoringDates,
                              const DateTime& today);

    virtual void createInterpBarrier(const DateTime& valueDate,
                                     const DateTimeArray& samples)const;

    virtual void overrideEventBarrier(const DateTime& valueDate)const;

    // Allow barrier to influence mon dates in inst if Daily
    virtual DateTimeArray getFutureMonitoringDates(const DateTime&         today,
                                                   const DateTimeArray&    instMonitorDates, 
                                                   smartPtr<IMCPathConfig> pathConfig);

    virtual bool createBarrierUtil(const DateTime& valueDate) const;

    virtual refCountPtr<IBarrierUtil> getBarrierUtil() const;
    
    virtual void pathUpdated(const IMCPathGenerator* pathGen);
        
    virtual ScheduleSP getBarrierSchedule() const;

    virtual double getInterpLevel(const IMCPathGenerator* pathGen,
                                  const IMCProduct* product,
                                  int iAsset)const;

    virtual void adjustLN(const IMCPathGenerator* pathGen,
                          const IMCProduct*        product,
                          const IMultiFactors*    assets,
                          int                     pathIdx);

    virtual void adjustImplied(const IMCPathGenerator*  pathGen,
                               const IMCProduct*         product,
                               const IMultiFactors*     assets,
                               const YieldCurve*        discount,
                               InstrumentSettlement*    instSettle,
                               int                      pathIdx);

    virtual double hitValue(const  IMCPathGenerator*  pathGen,
                            int    assetIdx);
    virtual double hitValue(const  IMCPathGenerator*  pathGen);

    virtual void hitValueAndTime(const   IMCPathGenerator*  pathGen,
                                 int     assetIdx,
                                 double& value,
                                 bool&   isConditionMet,
                                 int&    metAtStep);
    virtual void hitValueAndTime(const  IMCPathGenerator*  pathGen,
                                 double& value,
                                 bool&   isConditionMet,
                                 int&    metAtStep);

    virtual void hitValueAndTimeAsBasket(const  IMCPathGenerator*  pathGen,
                                         const  IAggregateMakerSP& bsk,
                                         double& value,
                                         bool&   isConditionMet,
                                         int&    metAtStep);

    virtual void hitValueAndTimeAtStep(const IMCPathGenerator*  pathGen,
                                       int                      iStep,
                                       int                      iAsset,
                                       double                   refLevel,
                                       double&                  value,           
                                       bool&                    isConditionMet,
                                       int&                     metAtStep);
    virtual void hitValueAndTimeAtMat(const IMCPathGenerator*  pathGen,
                                      int                      iStep,
                                      int                      iAsset,
                                      double                   refLevel,
                                      bool                     isConditionMetSoFar,
                                      double&                  value,           
                                      bool&                    isConditionMet,
                                      int&                     metAtStep);
    virtual void multiHelper(DoubleArray&       hitValues,
                             const BoolArray&   isConditionMet,
                             const IntArray&    metAtStep,
                             double&            overallHitValue,
                             bool&              overallIsConditionMet,
                             int&               overallMetAtStep);

    virtual void multiHelper2(DoubleArray&       hitValues,
                              const BoolArray&   isConditionMet,
                              double&            overallHitValue,
                              double&            overallBinaryValue);
    
    virtual BarrierLevelArraySP reportLevels(const DateTime& valueDate,
                                             double          refLevel,
                                             int             assetIdx) const;

//    virtual BarrierLevelArraySP reportLevelsAsBasket(const  IAggregateMakerSP& bsk,
//                                                     const DateTime& valueDate,
//                                                     const DoubleArray assetPerf,
//                                                     double          refLevel,
//                                                     int             assetIdx) const;

    //// roll through time replacing risk barrier in past by legal one
	virtual bool sensShift(Theta* theta);

    /** Satisfy LegalTerms::Shift interface */
    virtual bool sensShift(LegalTerms* shift);

    // used by NFB and RKO to get events
    virtual void getEvents(const  IMCPathGenerator*  pathGen, 
                            EventResults* events,
                            bool useIsHitFlags,
                            const string barrierName,
                            StringArray* assNames);

    // gets barrier breach events by looking at state after MC has run past
    // used by CKD
    virtual void getEvents(EventResults* events,
                            const string barrierName,
                            int   thisPeriodIdx,
                            const IntArrayArray& metAtStep,
                            const StringArray& assNames,
                            const DoubleMatrix& assLevels,
                            const DoubleMatrix& refLevels);

    virtual bool isMonitoringPoint(const DateTimeArray& monitoringDates,
                                    const DateTime& today);

protected:
    virtual void hitValueAndTime(const IMCPathGenerator*  pathGen,
                                 int        beginStep,
                                 int        endStep,
                                 const double *bendAmts,
                                 double&    value, 
                                 bool&      isConditionMet, 
                                 int&       metAtStep,
                                 int&       numRealHits,
                                 CBoolArray *isHitPerAsset,
                                 bool&      isConditionMetSoFar, 
                                 CBoolArray *isHitPerAssetSoFar, 
                                 int&       numRealHitsSoFar,
                                 int&       globalPastHitTime);
private:
    DoubleBarrier();

    friend class DoubleBarrierHelper;
};

typedef smartConstPtr<DoubleBarrier> DoubleBarrierConstSP;
typedef smartPtr<DoubleBarrier> DoubleBarrierSP;

// The union of all Barrier types to allow IMS to handle abstract types
class PRODUCTS_DLL BarrierUnion : public CObject {
public:
    static CClassConstSP const TYPE;
    ~BarrierUnion();
    friend class BarrierUnionHelper;
    virtual void validatePop2Object();
    
    Barrier* getBarrier(){
        return realBarrier.get();
    }

protected:
    BarrierUnion(CClassConstSP clazz);

private:
    BarrierUnion();
    BarrierUnion(const BarrierUnion& rhs);
    BarrierUnion& operator=(const BarrierUnion& rhs);


    string                              barrierType;     // Simple, Schedule, MultiSchedule
    BarrierScheduleSP                   barrierSchedule;
    smartPtr<DoubleBarrier>             doubleBarrier;
    
    // Transient
    BarrierSP  realBarrier;

};

typedef smartPtr<BarrierUnion> BarrierUnionSP;

DRLIB_END_NAMESPACE

#endif

