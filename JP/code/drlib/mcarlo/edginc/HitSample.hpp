//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HitSample.hpp
//
//   Description : Class for computing probabilities and sampling of hitting times etc.
//
//   Date        : February 2004
//
//----------------------------------------------------------------------------
#ifndef HITSAMPLE_HPP
#define HITSAMPLE_HPP


#include "edginc/Class.hpp"
#include "edginc/FunctionWrapper.hpp"
#include "edginc/Dividend.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/** Interface for HitSample object */
class MCARLO_DLL HitSample: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    /** Virtual destructor */
    virtual ~HitSample() {};

    /** A sample for 1(tau<T) given a uniform u */
    virtual int hitNoHitSample(double u) const = 0;

    /** A sample for tau given a uniform u */
    virtual double hitTimeSample(double u) const = 0;

    /** A sample for tau|tau<T given a uniform u */
    virtual double hitTimeGivenHitSample(double u) const = 0;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP &clazz);
};

typedef smartPtr<HitSample> HitSampleSP;
typedef smartConstPtr<HitSample> HitSampleConstSP;
typedef vector<HitSampleSP> HitSampleArray;

////////////////////////////////////////////////////////////////////////

/** Computes the Hit Sample distributions */
class MCARLO_DLL HitDistributionBM: public CObject,
                         virtual public HitSample {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    HitDistributionBM(double valueStart, double valueEnd, double barrierStart,
        double barrierEnd, double vol, double length, bool isUp);

    /** Partial constructor assuming that updateBridge will be called */
    HitDistributionBM(double barrierStart, double barrierEnd,
        double vol, double length, bool isUp);

    /** The cumulative distribution function for 0<tau<t<T */
    double hitTimeProb(double t) const;

    /** The cumulative distribution function for tau|0<tau<t<T */
    double hitTimeGivenHitProb(double t) const;

    /** The hitting probability within the BB */
    double getHittingProb() const;

    /** A sample for 1(tau<T) given a uniform u */
    int hitNoHitSample(double u) const;

    /** A sample for tau given a uniform u */
    double hitTimeSample(double u) const;

    /** A sample for tau|tau<T given a uniform u */
    double hitTimeGivenHitSample(double u) const;

    /** Updates the bridge for a different path */
    void updateBridge(double valueStart, double valueEnd);

private:
    /** Default constructor */
    HitDistributionBM();

    /** Returns if bridge hits at start/end */
    bool hitEndpoint(bool start) const;

    /** Validates and populates transient fields */
    virtual void validatePop2Object();

    /** Updates transient fields w.r.t. changes of the BB endpoints */
    void update();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultHitDistributionBM();

    // Fields
    double valueStart;          //!< Spot at start of BB
    double valueEnd;            //!< Spot at end of BB
    double barrierStart;        //!< Barrier at start of BB
    double barrierEnd;          //!< Barrier at end of BB
    double vol;                 //!< Integrated variance from start to end
    double length;              //!< BB length "T"
    bool   isUp;                //!< Flag for Up/Down barrier

    // Transient fields
    double hitProb;             //!< The hitting probability withing the interval
    double commonNum;
    double num1;
    double num2;
    bool   hitStart;            //!< Whether the bridge has hit the barrier at start
    bool   hitEnd;              //!< Whether the bridge has hit the barrier at the end
    bool   zeroVol;             //!< Whether there is zero vol in this BB
    bool   zeroLength;          //!< Whether there is zero trading time in this BB

    /** Wrapper class for inverting cumulative distributions */
    typedef const HitDistributionBM* Ptr;
    typedef double (HitDistributionBM::* Func1DConstPtr)(double x) const;
    typedef MemFuncWrapper<Ptr, Func1DConstPtr> HitDistWrapper;
    typedef FuncShiftWrapper<HitDistWrapper>  HitDistSolveWrapper;
    DECLARE_REF_COUNT(HitDistSolveWrapper);

    // Mutable in order to reset the target within a const sample method
    mutable HitDistSolveWrapperSP probTime; // $unregistered
    mutable HitDistSolveWrapperSP probCondTime; // $unregistered

    static const double INVERSION_TOLERANCE;   //!< Tolerance used in inverting distribution
    static const double HUGE_HITTING_TIME;     //!< Infinite hitting time i.e. no hit

    /** Helper function for hitting time computation */
    double hitTimeSampleHelper(HitDistSolveWrapperSP& wrapper,
        double tLow, double tHigh, double u) const;
};

typedef smartPtr<HitDistributionBM> HitDistributionBMSP;
typedef smartConstPtr<HitDistributionBM> HitDistributionBMConstSP;
typedef vector<HitDistributionBMSP> HitDistributionBMArray;

////////////////////////////////////////////////////////////////////////

class MCARLO_DLL HitNoHitBB {
public:
    /** Partial constructor assuming that updateBridge will be called */
    HitNoHitBB(double barrierStart, double barrierEnd,
               double var, double length, bool isUp, const DividendArray& divArray, const DoubleArray& divTimes);

    /** Partial constructor assuming that updateBridge will be called */
    HitNoHitBB(double var, double length, bool isUp, const DividendArray& divArray, const DoubleArray& divTimes);

    /** Destructor */
    ~HitNoHitBB() {}

    /** Returns a sample for the Brownian Bridge*/
    int updateAndSample(double valueStart, double ValueEnd, double uniform);

    /** Returns a sample for the Brownian Bridge*/
    int updateAndSample(double barrierStart, double barrierEnd, double valueStart, double ValueEnd, double uniform);

private:

    /** Computes the dividend weight for the forward brownian bridge */
    void computeAdjustments(const DividendArray& divArray, const DoubleArray& divTimes);

    /** Computes the dividend weights */
    double computeWeight(double t, double CT, int Wtype);

    /** Validates and populates the fields */
    void validate() ;

    /** Returns the hitting probability for the Brownian Bridge with one dividend */
    void updateBridgeOneDiv();

    /** Updates the bridge for a different path */
    void updateBridge(double valueStart, double valueEnd);

    /** Updates the bridge for a different path and barrier */
    void updateBridge(double barrierStart, double barrierEnd, double valueStart, double valueEnd);

    /** A sample for 1(tau<T) given a uniform u */
    int hitNoHitSample(double u) const;

    /** Returns whether the endpoints are hit */
    bool areEndpointsHit() const;

    /** Returns the hitting probability for the Brownian Bridge */
    double getHitProb() const;

#if 0

    /** Computes the dividend weight for the forward brownian bridge */
    double fwdBBDivAdjustment(const DividendArray& divArray, const DoubleArray& divTimes) const;

    /** Computes the dividend weight for the forward brownian bridge */
    double revBBDivAdjustment(const DividendArray& divArray, const DoubleArray& divTimes) const;

    /** Produces a hitting time sample conditional on a hit event */
    static double hitTimeProb(
        double x, double y,
        double bStart, double bEnd,
        double sqrtVar, double length, bool isUp,
        double u);

#endif

    // Fields
    double valueStart;          //!< Spot at start of BB
    double valueEnd;            //!< Spot at end of BB
    double barrierStart;        //!< Barrier at start of BB
    double barrierEnd;          //!< Barrier at end of BB
    double var;                 //!< Integrated variance from start to end
    double length;              //!< Time length of Brownian Bridge (time metric dependent)
    bool   isUp;                //!< Flag for Up/Down barrier
    double forAdj;              //!< Forward dividend adjustment
    double revAdj;              //!< Reverse dividend adjustment
    double intAdj;              //!< Intermediate dividend adjustment
    double CTime ;              //!< Dividend centered Time

    // Transient fields
    double hitProb;             //!< The hitting probability withing the interval
    bool   zeroVar;             //!< Whether there is zero variance in this BB
    double sqrtVar;             //!< Square root of variance
    bool   isHit;               //!< Whether endpoints hit
};

DECLARE_REF_COUNT(HitNoHitBB);

DRLIB_END_NAMESPACE

#endif
