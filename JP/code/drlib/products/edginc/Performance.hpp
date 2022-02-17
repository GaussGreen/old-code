//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Performance.hpp
//
//   Description : Measures performance of assets
//
//   Author      : Mark A Robson
//
//   Date        : 24 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PERFORMANCE_HPP
#define EDR_PERFORMANCE_HPP

#include "edginc/DateTime.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"

#define STATE_VARIABLES

#ifdef STATE_VARIABLES
#include "edginc/StateVariable.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/SVGenSpot.hpp"

#endif

DRLIB_BEGIN_NAMESPACE

#ifdef STATE_VARIABLES
class IPerformanceSV;
typedef refCountPtr<IPerformanceSV> IPerformanceSVSP;
typedef vector<IPerformanceSVSP> IPerformanceSVArray;
typedef refCountPtr<IPerformanceSVArray> IPerformanceSVArraySP;

///////////////////////////////////////////////////////////////////////////

/** Provides ability to calculate a 'performance' of a single asset */
class PRODUCTS_DLL IPerformanceSV: virtual public IStateVariableGen,
                      virtual public IStateVariableClient {
public:
    /** Defines an interface geared for use by the monte carlo - it contains
        'state' information and provides a clear distinction between user
        interface data and state variables as needed by the Monte Carlo */
    class PRODUCTS_DLL IStateVar: public virtual IStateVariable {
    public:
        /** Calculates the performance of ith asset (i'th asset in the
            path generators view of the world) */
        virtual void calcPerf(DoubleArray& perf) const = 0;

        /** Provides a LN Vol request for this performance measurement */
        virtual CVolRequestLN* getVolInterp(const MCProductClient* mcProduct, 
                                            int iAsset) const = 0;

        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level */
        virtual double maxDeltaScalingFactor(int iAsset) const = 0;

        virtual ~IStateVar() {}
    };
    typedef smartPtr<IStateVar> IStateVarSP;

    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries, int iAsset) const = 0;

    /** Returns the Performance state var object for this Performance object */
    virtual IStateVarSP getPerfStateVar(IStateVarSP               oldStateVar,
                                        IStateVariableGen::IStateGen* pathGen) const = 0;
    
    virtual ~IPerformanceSV() {}
};

///////////////////////////////////////////////////////////////////////////

class PerformanceSV;
typedef refCountPtr<PerformanceSV> PerformanceSVSP;
typedef vector<PerformanceSVSP> PerformanceSVArray;
typedef refCountPtr<PerformanceSVArray> PerformanceSVArraySP;


/** Simple performance where perfType, strike, participation, and 
    sampleDates are the same across all assets (a better name anyone?) */
class PRODUCTS_DLL PerformanceSV: virtual public IPerformanceSV {
private:
    // Fields
    StringArray               perfTypePerAsset;         //!< Pertormance type
    DoubleArray               strikePerAsset;           //!< Strikes
    DoubleArray               participationPerAsset;    //!< Participations
    SVGenSpotSP                  spotPathGen;              //!< Spot path generator
    IRefLevel::IStateVarGenSP refLevelGen;              //!< Ref Level generator

public:
    class StateVar;
    typedef smartPtr<StateVar> StateVarSP;
    friend class StateVar;

    /** Simple performance state variable */
    class PRODUCTS_DLL StateVar: virtual public IPerformanceSV::IStateVar {
    private:
        // Fields
        const PerformanceSV*       performance;
        int                        nbAssets;
        mutable DoubleArray        sumOut;
        bool                       isPast;         //!< Whether we process the past or future
        SVGenSpotSP                   spotPathGen;    //!< Spot path generator
        IRefLevel::IStateVarGenSP  refLevelGen;    //!< Ref Level generator

        // Transient fields
        SVGenSpot::IStateVarSP        spotPath;       //!< Spot path state variable
        IRefLevel::IStateVarSP     refLevel;       //!< Reference level for barrier

    public:
        virtual bool doingPast() const;

        /** Returns the number of performances that this object will be
            calculating ie the size of the perf parameter in calcPerf() */
        virtual int numPerformances() const;

        /** Calculates collectively the performances for all assets */
        virtual void calcPerf(DoubleArray& perf) const;
        
        // for the LogNormal path generator
        virtual CVolRequestLN* getVolInterp(const MCProductClient*  mcProduct,
                                            int                     iAsset) const;
        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level  */
        virtual double maxDeltaScalingFactor(int iAsset) const;

        /** Allows the barrier class to switch from past to future */
        void update(IStateVariableGen::IStateGen* pathGen);

        StateVar(const PerformanceSV*             performance, 
                 const IRefLevel::IStateVarGenSP& refLevelGen,
                 const SVGenSpotSP&                  spotPathGen);

    protected:
        /** Calculates the performance of ith asset */
        double calcPerf(int iAsset) const;
    };

public:
    /** Full constructor */
    PerformanceSV(const StringArray&        perfTypePerAsset,
                  const DoubleArray&        strikePerAsset,
                  const DoubleArray&        participationPerAsset,
                  SVGenSpotSP                  spotPathGen,
                  IRefLevel::IStateVarGenSP refLevelGen);

    /** Implementation of IStateVariableClient */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Returns the Performance state var object for this Performance object */
    virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                 IStateVariableGen::IStateGen* pathGen) const;

    /** Returns the Performance state var object for this Performance object */
    virtual IStateVarSP getPerfStateVar(IStateVarSP               oldStateVar,
                                        IStateVariableGen::IStateGen* pathGen) const;

    virtual void recordDates(SimSeries* simSeries, int iAsset) const;

    static PerformanceSVSP mergePerfs(const PerformanceSVArray& perfs,
                                      IRefLevel::IStateVarGenSP refLevelGen);
};

typedef refCountPtr<IPerformanceSVArray> IPerformanceSVArraySP;

#endif 

/** Provides ability to calculate a 'performance' of a single asset */
class PRODUCTS_DLL IAssetPerformance: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    /** Defines an interface geared for use by the monte carlo - it contains
        'state' information and provides a clear distinction between user
        interface data and state variables as needed by the Monte Carlo */
    class PRODUCTS_DLL IMCPerf{
    public:     
        /** Calculates the performance of ith asset (i'th asset in the
            path generators view of the world) */
        virtual double calcPerf(const MCPathGenerator*  pathGen,
                                int                      iPath,
                                int                      iAsset) const = 0;      

        /** Provides a LN Vol request for this performance measurement */
        virtual CVolRequestLN* getVolInterp(
            const IMCProduct*        mcProduct,
            const MCPathGenerator* pathGen,
            int                     iAsset) const = 0;

        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level */
        virtual double maxDeltaScalingFactor(
            const MCPathGenerator* futurePathGen, 
            int                     iAsset) const = 0;

        virtual ~IMCPerf(){}
    };
    typedef refCountPtr<IMCPerf> IMCPerfSP;

    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries, int iAsset) const = 0;

    /** Returns the IMCPerf object for this Performance object */
    virtual IMCPerf* getMCPerf(const DateTime&    today,
                               const SimSeries*   allDates,
                               int                iAsset) const = 0;

    /** Converts performance object to state variable representation */
    virtual PerformanceSVSP getStateVarGen(
        IRefLevel::IStateVarGenSP refLevelGen, int nbAssets) = 0;

    ~IAssetPerformance(){}
private:
    static void load(CClassSP& clazz);

};

typedef smartPtr<IAssetPerformance> IAssetPerformanceSP;
typedef array<IAssetPerformanceSP, IAssetPerformance> IAssetPerformanceArray;
typedef smartPtr<IAssetPerformanceArray> IAssetPerformanceArraySP;


/** Provides ability to calculate a 'performance' (which
    means an array of doubles (one per asset) */
class PRODUCTS_DLL IPerformance: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    // ideally these would be static const char's but MSVC doesn't let you
    // define them in the header file, and if you put them in the source file
    // you can't then switch (as in the keyword) on them - at least not in
    // another file
#define PERF_TYPE_FORWARD   'F'
#define PERF_TYPE_CALL      'C'
#define PERF_TYPE_PUT       'P'
#define PERF_TYPE_STRADDLE  'S'

    /** Defines an interface geared for use by the monte carlo - it contains
        'state' information and provides a clear distinction between user
        interface data and state variables as needed by the Monte Carlo */
    class PRODUCTS_DLL IMCPerf{
    public:
        /** Calculates the performance of each asset returning the results
            in perf (which must of the required size) */
        virtual void calcPerf(const MCPathGenerator*  pathGen,
                              int                      iPath,
                              DoubleArray&             perf) const = 0;
             
        /** Returns the number of performances that this object will be
            calculating ie the size of the perf parameter in calcPerf() */
        virtual int numPerformances() const = 0;
       
        /** Provides a LN Vol request for this performance measurement */
        virtual CVolRequestLN* getVolInterp(
            const IMCProduct*        mcProduct,
            const MCPathGenerator* pathGen,
            int                     iAsset) const = 0;
       
        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. It assumes the refLevel()
            function contains valid estimates of the ref level  */
        virtual double maxDeltaScalingFactor(
            const MCPathGenerator* futurePathGen, 
            int                     iAsset) const = 0;

        virtual ~IMCPerf(){}
    };
    typedef refCountPtr<IMCPerf> IMCPerfSP;

    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries) const = 0;

    /** Returns the IMCPerf object for this Performance object */
    virtual IMCPerf* getMCPerf(const DateTime&    today,
                               const SimSeries*   allDates) const = 0;

    /** Converts performance object to state variable representation */
    virtual PerformanceSVSP getStateVarGen(
        IRefLevel::IStateVarGenSP refLevelGen, int nbAssets) = 0;

    ~IPerformance(){}

    /** Contains useful [static] methods */
    class PRODUCTS_DLL Util{
    public:
        /** Validate the perfType flag (must be 'F', 'C', or 'P'), checks
            the strike isn't negative, and checks the sampleDates are
            increasing */
        static void validatePerfFlags(const string&        perfType,
                                      double               strike,
                                      const DateTimeArray& sampleDates);

        /** switching on perfType:
            (not recognised) - returns originalPerf
            (PERF_TYPE_CALL) - returns Maths::max(0.0, perf[iAsset])
            (PERF_TYPE_PUT) - returns Maths::max(0.0, -perf[iAsset])
            [make inline in case performance an issue] */
        static double calcPerf(const string&        perfType,
                               double               originalPerf){
            switch (perfType[0]) {
                /*case 'F':
                  Nothing extra here
                  break;*/
            case PERF_TYPE_CALL:
                return Maths::max(0.0, originalPerf);
            case PERF_TYPE_PUT:
                return Maths::max(0.0,-originalPerf);
            case PERF_TYPE_STRADDLE:
                return fabs(originalPerf);
            }
            return originalPerf; /* for performance (clients should check flags
                                    before hand) */
        }
    };
private:
    static void load(CClassSP& clazz);

};

typedef smartPtr<IPerformance> IPerformanceSP;
typedef array<IPerformanceSP, IPerformance> IPerformanceArray;
typedef smartPtr<IPerformanceArray> IPerformanceArraySP;

DRLIB_END_NAMESPACE

#endif

