//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenBarrier.hpp
//
//   Description : Barrier state variable
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#ifndef SVGenBarrier_HPP
#define SVGenBarrier_HPP

#include "edginc/StateVariable.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/RefLevel.hpp"
#include "edginc/MCPathGenerator.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/StateVariableGen.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Results.hpp"
#include "edginc/Asset.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_REF_COUNT(SVGenBarrierHV);
//FORWARD_DECLARE_REF_COUNT(SVGenBarrierHVBB);
//FORWARD_DECLARE_REF_COUNT(SVGenBarrierHVBarAdj);


// Data structures

/** Barrier data per asset */
class MCARLO_DLL BarrierPerAssetData {
public:
    /** Full constructor */
    BarrierPerAssetData(const string& monitorType,
                        ScheduleSP levels,
                        ScheduleSP economicLevels,
                        DateTimeArraySP monitoringDates,
                        const DateTime& smoothDate,
                        const DateTime& valueDate,
                        bool isOut,
                        bool isUp,
                        bool isHit);

    ~BarrierPerAssetData() {}

    // Fields
    ScheduleSP       levels;                //!< Schedule of dates and levels
    ScheduleSP       economicLevels;        //!< Economic schedule of dates and levels
    DateTimeArraySP  monitoringDates;       //!< Monitoring dates
    bool             isOut;                 //!< In or Out
    bool             isUp;                  //!< Up or Down
    bool             isHit;                 //!< Whether asset has hit already for "C" or "D"

    // Transient fields
    DoubleArray      barrierLevels;         //!< The barrier levels at monitoring dates
    double           barrierSmoothDate;     //!< Barrier level at smooth date
    IntArray         barMap;                // map of Barrier Date in monitoringDates.

    /** Returns the first future monitoring date for "D" type monitoring
        in a DateTimeArray of size 1. If the barrier is in the past or if
        monitorType is different from "D" an empty DateTimeArray is returned */
    DateTimeArray getFirstFutureDailyMonDate(
        const CAsset& asset, const DateTime& valueDate, const string& monitorType) const;

private:
    BarrierPerAssetData();                                          //!< Disabled empty constructor
    BarrierPerAssetData(const BarrierPerAssetData& rhs);            //!< Disabled copy constructor
    BarrierPerAssetData& operator=(const BarrierPerAssetData& rhs); //!< Disabled assignment operator
};

DECLARE_REF_COUNT(BarrierPerAssetData);

/** Data specifying smoothing for barrier monitoring */
class MCARLO_DLL SmoothingData {
public:
    /** Full constructor */
    SmoothingData(bool smoothing, double lowSpread, double highSpread,
        const DateTime& smoothDate);

    /** Destructor */
    ~SmoothingData() {}

    // Fields
    bool        smoothing;        //!< Use smoothing
    double      lowSpread;        //!< Low strike for call spread
    double      highSpread;       //!< High strike for call spread
    DateTime    smoothDate;       //!< Date at which to apply smoothing (maturity)

private:
    SmoothingData();                                    //!< Disabled empty constructor
    SmoothingData(const SmoothingData& rhs);            //!< Disabled copy constructor
    SmoothingData& operator=(const SmoothingData& rhs); //!< Disabled assignment operator
};

DECLARE_REF_COUNT(SmoothingData);

/** Barrier data */
class MCARLO_DLL BarrierData {
public:
    /** Full constructor */
    BarrierData(const string& monitorType,
                const string& closedFormDates,
                const string& closedFormMethod,
                const string& maxNumDivsPerYearString,
                const vector<BarrierPerAssetDataSP>& barrierPerAsset,
                int numHits,
                SmoothingDataSP smoothingData,
                const DateTime& hitDate,
                const DateTime& valueDate);

    /** Destructor */
    ~BarrierData() {}

    /** Indicates if the barrier is european or CF */
    bool european() const;

    // Fields
    // Barrier definition
    string                        monitorType;          //!< "E", "D" or "C"
    string                        closedFormDates;      //!< "INPUT", "ENDPOINTS", "BENCHMARK", "AUTO"
    string                        closedFormMethod;     //!< "BROWNIAN_BRIDGE", "BARRIER_ADJUSTMENT"
    int                           numHits;              //!< How many hits trigger a breach
    DateTime                      hitDate;              //!< Date structure hit for "C" or "D"
    DateTime                      valueDate;            //!< Value date

    // Transient fields
    int                           maxNumDivsPerYear;    //!< Maximum number of divs per year
    int                           numAssets;            //!< Number of assets
    bool                          isBarrierInPast;      //!< All dates are in the past
    bool                          allAssetsHit;         //!< All assets hit in past for CF

    /** Updates the isHit flag for iAsset. Used when past gets processed */
    void amendIsHit(bool isAssetHit, int iAsset);

    // Monitoring flags
    static const string EUROPEAN_MONITORING;        //!< "E"
    static const string DAILY_MONITORING;           //!< "D"
    static const string CONTINUOUS_MONITORING;      //!< "C"


    static const string DEFAULT;                    //!< Default flag used in various places

    // Specify how to what constitutes the monitoring timeline for closed form
    static const string INPUT;                      //!< "C", "D" monitoring between input dates
    static const string BENCHMARK;                  //!< "C", "D" add benchmark between start and end dates
    static const string ENDPOINTS;                  //!< "C", "D" monitor at beginning and end only
    static const string AUTO;                       //!< "C", "D" monitor at beginning and end only

    // Closed form methodology
    static const string BROWNIAN_BRIDGE;            //!< BB
    static const string BARRIER_ADJUSTMENT;         //!< ADJ

    // Dividend treatment for Brownian Bridges
    static const string DEFAULT_MAX_DIVS_PER_YEAR;  //!< Default max divs per year
    static const string USE_ALL_DIVS;               //!< Use all dividends in Brownian Bridge
    static const string NO_DIVS;                    //!< Don't use dividends at all
    static const int MIN_DIVS_PER_PERIOD;           //!< Minimum number of divs per period irrespective of length

    /** Validates the combination of strings */
    static void validateMethodology(const string& monitorType, bool useStateVars,
                                    string& method, string& dates);

    /** Accessor for asset barriers. Returns const ConstSP& in case client
        needs to be mega optimal (these are refCounts and copying costs) */
    /// But nothing costs more than bugs
    BarrierPerAssetDataSP getAssetBarrier(int iAsset) const;

    /** Accessor for smoothing data. Returns const ConstSP& in case client
        needs to be mega optimal (these are refCounts and copying costs) */
    SmoothingDataSP getSmoothingData() const;

private:
    vector<BarrierPerAssetDataSP> barrierPerAsset;  //!< Barrier for each asset
    SmoothingDataSP               smoothingData;    //!< Smoothing specification

    /** Populates allAssetsHit field */
    void amendAllAssetsHit();

    BarrierData();                                  //!< Disabled empty constructor
    BarrierData(const BarrierData& rhs);            //!< Disabled copy constructor
    BarrierData& operator=(const BarrierData& rhs); //!< Disabled assignment operator
};

DECLARE_REF_COUNT(BarrierData);

//////////////////////////////////////////////////////////////////////


/** SVGenBarrierHV class can generate a HitValue state variable. */
class MCARLO_DLL SVGenBarrierHV: virtual public IStateVariableGen {
public:
    /** Interface for Barrier Hit Value state variable */
    class MCARLO_DLL IStateVar: public virtual IStateVariable {
    public:
        /** Returns a hit / no hit value per asset. It's a double because
            we might use smoothing */
        virtual double hitValue(int iAsset) const = 0;

        /** Records BARRIER_LEVEL output requests in results */
        virtual void recordBarrierLevels(CControl* control,
                                         Results*  results,
                                         const IMultiFactors* mAsset) const = 0;
    };

    class MCBarrierHVStateVarTrivial;   //!< Implementation for past closed form

    typedef smartPtr<IStateVar> IStateVarSP;
    typedef smartConstPtr<IStateVar> IStateVarConstSP;
    typedef vector<IStateVarSP> IStateVarArray;

    /** Destructor */
    virtual ~SVGenBarrierHV() {}

    /** Method that returns a hit value state variable */
    virtual IStateVarSP getHitValueSV(IStateVarSP               oldStateVar,
                                      IStateVariableGen::IStateGen* pathGen) const = 0;

    /** Allows access to data */
    virtual BarrierDataConstSP getBarrierData() const = 0;

    /** Helper function to deduce the smoothed hitting value for
        a particular asset */
    static double hitValuePayoff(const SmoothingDataConstSP& smoothingData,
                                 const BarrierPerAssetDataConstSP& assetBarrier,
                                 const SVPath& spotSmoothPath,
                                 double refLevel,
                                 bool isAssetHit,
                                 bool isAssetHitInPast,
                                 bool doingPast);
};


//////////////////////////////////////////////////////////////////////


/** The class can generate a European type HitValue state variable.
    It is a client of state variables. */
class MCARLO_DLL SVGenBarrierHVEur: virtual public SVGenBarrierHV,
                         virtual public IStateVariableClient {
public:
    /** European barrier implementation */
    class MCARLO_DLL StateVar: virtual public SVGenBarrierHV::IStateVar {
    public:
        /** Constructor */
        StateVar(BarrierDataConstSP data,
                 IRefLevel::IStateVarGenSP refLevelGen,
                 SVGenSpotSP spotMonGen,
                 SVGenSpotSP spotSmoothGen);

        /** Allows the barrier class to switch from past to future */
        void update(IStateVariableGen::IStateGen* pathGen);

        /** Hit value per asset */
        double hitValue(int iAsset) const;

        /** Indicates past Vs future */
        bool doingPast() const;

        /** Records BARRIER_LEVEL output requests in results */
        virtual void recordBarrierLevels(CControl* control,
                                         Results*  results,
                                         const IMultiFactors* mAsset) const;

    private:
        // Fields
        BarrierDataConstSP         data;              //!< Data
        IRefLevel::IStateVarGenSP  refLevelGen;       //!< Ref Level generator
        SVGenSpotSP                   spotMonGen;        //!< Spot generator for E type monitoring
        SVGenSpotSP                   spotSmoothGen;     //!< Spot generator for smoothing

        // Transient fields
        SVGenSpot::IStateVarSP        monitoringPath;    //!< Discrete path to monitor
        SVGenSpot::IStateVarSP        smoothPath;        //!< Spot path for smoothing
        IRefLevel::IStateVarSP     refLevel;          //!< Reference level for barrier

        bool                       isPast;            //!< Whether we process the past or future

        mutable BoolArray          hitNoHit;          //!< Working area
        mutable BoolArray          hitNoHitInPast;    //!< Whether asset hit in past
        mutable BoolArray          hitNoHitSoFar;     //!< Preservation area
    };

    typedef smartPtr<StateVar> StateVarSP;
    typedef smartConstPtr<StateVar> StateVarConstSP;

    /** Constuctor */
    SVGenBarrierHVEur(BarrierDataSP data,
                      IRefLevel::IStateVarGenSP refLevelGen);

    /** Implementation of IStateVariableClient */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Implementation of IStateVariableGen */
    IStateVarSP getHitValueSV(IStateVarSP               oldStateVar,
                              IStateVariableGen::IStateGen* pathGen) const;

    /** Part of IStateVariableGen */
    IStateVariableSP create(IStateVariableSP             oldStateVar,
                         IStateVariableGen::IStateGen* pathGen) const;

    /** Allows access to data */
    virtual BarrierDataConstSP getBarrierData() const;

    /** Allows access to ref level */
    IRefLevel::IStateVarGenSP  getRefLevel() const;

    /** Allows access to smooth state var */
    SVGenSpotSP getSmoothPath() const;

private:
    /** Disabled default constructor */
    SVGenBarrierHVEur();

    // Fields
    BarrierDataSP             data;            //!< The barrier data
    IRefLevel::IStateVarGenSP refLevelGen;     //!< Ref Level generator

    // Transient fields
    SVGenSpotSP                  spotMonGen;      //!< Spot generator for E type monitoring
    SVGenSpotSP                  spotSmoothGen;   //!< Spot generator for smoothing
};


//////////////////////////////////////////////////////////////////////


/** The class can generate a barrier adjusted HitValue state variable
    by delegating the creation to the path generator. It is both a
    client of state variables and an elementary state variable. */
class MCARLO_DLL SVGenBarrierHVBarAdj: virtual public SVGenBarrierHV,
                            virtual public IElemStateVariableGen,
                            virtual public IStateVariableClient {
public:
    /** Constuctor */
    SVGenBarrierHVBarAdj(BarrierDataSP data,
                         IRefLevel::IStateVarGenSP refLevelGen,
                         const IMultiFactors* mAsset);

    /** Implementation of IStateVariableClient */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Implementation of IStateVariableGen */
    IStateVarSP getHitValueSV(IStateVarSP               oldStateVar,
                              IStateVariableGen::IStateGen* pathGen) const;

    /** Part of IStateVariableGen */
    IStateVariableSP create(IStateVariableSP             oldStateVar,
                         IStateVariableGen::IStateGen* pathGen) const;

    /** Allows access to data */
    virtual BarrierDataConstSP getBarrierData() const;

    /** Allows access to ref level */
    IRefLevel::IStateVarGenSP  getRefLevel() const;

    /** Allows access to spot path state var */
    SVGenSpotSP getSpotPath() const;

    /** Allows access to smooth state var */
    SVGenSpotSP getSmoothPath() const;

    void attachSVGen(IElemStateVariableGenVisitor*) const;

private:
    /** Disabled default constructor */
    SVGenBarrierHVBarAdj();

    // Fields
    BarrierDataSP             data;            //!< The barrier data
    IRefLevel::IStateVarGenSP refLevelGen;     //!< Ref Level generator
    const IMultiFactors*      mAsset;          //!< Multi factors

    // Transient fields
    SVGenSpotSP                  spotMonGen;      //!< Spot generator for E type monitoring
    SVGenSpotSP                  spotSmoothGen;   //!< Spot generator for smoothing
};

//DECLARE_REF_COUNT(SVGenBarrierHVBarAdj);
typedef vector<const SVGenBarrierHVBarAdj*> SVGenBarrierHVBarAdjArray;
typedef refCountPtr<SVGenBarrierHVBarAdj> SVGenBarrierHVBarAdjSP;
typedef refCountPtr<const SVGenBarrierHVBarAdj> SVGenBarrierHVBarAdjConstSP;


//////////////////////////////////////////////////////////////////////

class MCProductClient;

/** The class can generate a Brownian Bridge HitValue state variable
    by delegating the creation to the path generator. It is both a
    client of state variables and an elementary state variable. */
class MCARLO_DLL SVGenBarrierHVBB: virtual public SVGenBarrierHV,
                        virtual public IElemStateVariableGen,
                        virtual public IStateVariableClient {
public:
    /** Brownian bridge Hit/NoHit state variable. Deduces hit no hit value by
        reading from the simulated Brownan Bridges hit no hit value path. */
    class MCARLO_DLL StateVar: virtual public SVGenBarrierHV::IStateVar {
    public:
        /** Full constructor */
        StateVar(IStateVariableGen::IStateGen*   pastPathGen,
                 IStateVariableGen::IStateGen*   futPathGen,
                 BarrierDataConstSP          data,
                 IRefLevel::IStateVarGenSP   refLevelGen,
                 SVGenSpotSP                    spotSmoothGen,
                 SVGenSpotSP                    spotAtMonStartGen,
                 SVGenSpotSP                    spotAtMonEndGen,
                 const vector<IntArray>&     offsets,
                 const vector<IntArray>&     hitNoHitPath);

        /** Returns a hit / no hit value per asset. It's a double because
            we might use smoothing. Part of the SVGenBarrierHV::IStateVar. */
        virtual double hitValue(int iAsset) const;

        /** Part of the IStateVariable. */
        virtual bool doingPast() const;

        /** Records BARRIER_LEVEL output requests in results */
        virtual void recordBarrierLevels(CControl* control,
                                         Results*  results,
                                         const IMultiFactors* mAsset) const;

    private:
        StateVar();                                 //!< Disabled default constructor
        StateVar(const StateVar& rhs);              //!< Disabled copy constructor
        StateVar& operator=(const StateVar& rhs);   //!< Disabled assignment operator

        // Fields
        vector<int>                  monitorAtStart;     //!< Monitor explicitly first date for "D" monitoring
        vector<int>                  monitorAtEnd;       //!< Monitor explicitly last date for "S" barriers
        DoubleArray                  barrierPctStart;    //!< Barrier at start
        DoubleArray                  barrierPctEnd;      //!< Barrier at end
        BarrierDataConstSP           data;               //!< Barrier data
        vector<IntArray>             offsets;            //!< Offsets for reading from path
        const vector<IntArray>&      hitNoHitPath;       //!< Hit No Hit paths

        IRefLevel::IStateVarGenSP    refLevelGen;        //!< Ref Level generator
        SVGenSpotSP                     spotSmoothGen;      //!< Spot generator for smoothing
        SVGenSpotSP                     spotAtMonStartGen;  //!< Spot generator for monitoring start
        SVGenSpotSP                     spotAtMonEndGen;    //!< Spot generator for monitoring end

        IRefLevel::IStateVarSP       refLevel;           //!< RefLevel state variable
        SVGenSpot::IStateVarSP          smoothPath;         //!< SmoothPath state variable
        SVGenSpot::IStateVarSP          spotAtMonStartPath; //!< Spot at monitoring start state variable
        SVGenSpot::IStateVarSP          spotAtMonEndPath;   //!< Spot at monitoring end state variable
    };
    typedef smartPtr<StateVar> StateVarSP;


    /** Constuctor */
    SVGenBarrierHVBB(BarrierDataSP data,
                     IRefLevel::IStateVarGenSP refLevelGen,
                     const IMultiFactors* mAsset);

    /** Implementation of IStateVariableClient */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Implementation of IStateVariableGen */
    IStateVarSP getHitValueSV(IStateVarSP               oldStateVar,
                              IStateVariableGen::IStateGen* pathGen) const;

    /** Part of IStateVariableGen */
    IStateVariableSP create(IStateVariableSP             oldStateVar,
                         IStateVariableGen::IStateGen* pathGen) const;

    /** Allows access to data */
    virtual BarrierDataConstSP getBarrierData() const;

    /** Allows access to ref level */
    IRefLevel::IStateVarGenSP  getRefLevelGen() const;

    /** Allows access to smooth state var */
    SVGenSpotSP getSpotSmoothGen() const;

    /** Allows access to monitoring at start */
    SVGenSpotSP getSpotAtMonStartGen() const;

    /** Allows access to monitoring at end */
    SVGenSpotSP getSpotAtMonEndGen() const;

    void attachSVGen(IElemStateVariableGenVisitor*) const;

    /** Provides certain helper utilities for creating a barrier
        timeline based on various rules */
    class MCARLO_DLL BarrierDatesHelper {
    public:
        /** Helper function that creates barrier timeline from input data
            for a set of Brownian Bridge barriers based on a string. The
            resulting dates will be one of i) endpoints only or ii) input
            dates only or iii) endpoints with benchmark dates inserted. */
        static DateTimeArraySP getBarrierTimeline(
            const MCProductClient* prod,
            const vector<const SVGenBarrierHVBB*>& barrierGens,
            DateTimeArray& simDates);

       /** Creates an array of dates that are equal to startDate offset
           by each MaturityPeriod and have the input time */
        static DateTimeArray createDates(const DateTime& startDate,
                                         const MaturityPeriodArray& tenors,
                                         int time);
    private:
        /** Used to initialize static tenors field */
        static MaturityPeriodArray createTenors();

        static const MaturityPeriodArray tenors;    //!< Default VolSurface tenors
    };

private:
    /** Disabled default constructor */
    SVGenBarrierHVBB();

    // Fields
    BarrierDataSP             data;               //!< The barrier data
    IRefLevel::IStateVarGenSP refLevelGen;        //!< Ref Level generator
    const IMultiFactors*      mAsset;             //!< Multi factors

    // Transient fields
    SVGenSpotSP                  spotSmoothGen;      //!< Spot generator for smoothing
    SVGenSpotSP                  spotAtMonStartGen;  //!< Spot generator at start    (used for 1st day monitoring in "D")
    SVGenSpotSP                  spotAtMonEndGen;    //!< Spot generator at maturity (used for "S"tairs barriers)
};

//DECLARE_REF_COUNT(SVGenBarrierHVBB);
typedef vector<const SVGenBarrierHVBB*> SVGenBarrierHVBBArray;
typedef refCountPtr<SVGenBarrierHVBB> SVGenBarrierHVBBSP;
typedef refCountPtr<const SVGenBarrierHVBB> SVGenBarrierHVBBConstSP;

//////////////////////////////////////////////////////////////////////


/** SVGenBarrierHVStruc class is both a user of state variables
    (eg european barrier) or a generator of state variables */
class MCARLO_DLL SVGenBarrierHVStruct: public virtual IStateVariableGen,
                            public virtual IStateVariableClient {
public:
    /** Aggregates hits accross assets */
    class MCARLO_DLL StateVar: virtual public IStateVariable {
    public:
        /** Constructor */
        StateVar(SVGenBarrierHVSP barrierGen, bool isStructOut);

        /** Hit value method operates on all assets */
        double hitValue() const;

        /** Allows the barrier class to switch from past to future */
        void update(IStateVariableGen::IStateGen* pathGen);

        /** Indicates past Vs future */
        virtual bool doingPast() const;

        /** Records BARRIER_LEVEL output requests in results */
        void recordBarrierLevels(CControl* control,
                                 Results*  results,
                                 const IMultiFactors* mAsset) const;

    private:
        // Fields
        SVGenBarrierHVSP            barrierGen;     //!< The barrierHV generator

        // Transient fields
        SVGenBarrierHV::IStateVarSP barrierSV;      //!< Barrier state variable

        mutable DoubleArray         hitValues;      //!< Hit array for assets
        int                         numAssets;      //!< Number of assets
        int                         numHits;        //!< Number of hits that trigger event
        bool                        isStructOut;    //!< Out or In barrier
        bool                        isPast;         //!< Whether we process the past or future
    };

    typedef smartPtr<StateVar> StateVarSP;
    typedef smartConstPtr<StateVar> StateVarConstSP;

    /** Constuctor from barrier HitValue generator. The structure isStructOut
        flag is set to the common barrier isOut flag */
    SVGenBarrierHVStruct(SVGenBarrierHVSP barrierGen);


    /** Constuctor where isStructOut flag is explicitly specified */
    SVGenBarrierHVStruct(SVGenBarrierHVSP barrierGen, bool isStructOut);

    /** Destructor */
    virtual ~SVGenBarrierHVStruct() {}

    /** Implementation of IStateVariableClient */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Implementation of IStateVariableGen */
    StateVarSP getHitValueStructSV(StateVarSP                oldStateVar,
                                   IStateVariableGen::IStateGen* pathGen) const;

    /** Part of IStateVariableGen */
    IStateVariableSP create(IStateVariableSP             oldStateVar,
                         IStateVariableGen::IStateGen* pathGen) const;

    /** Allows access to data */
    virtual BarrierDataConstSP getBarrierData() const;

private:
    /** Disabled default constructor */
    SVGenBarrierHVStruct();

    // Fields
    SVGenBarrierHVSP  barrierGen;   //!< The barrierHV generator
    bool              isStructOut;  //!< Whether it is an In or Out structure
};
DECLARE_REF_COUNT(SVGenBarrierHVStruct);

//////////////////////////////////////////////////////////////////////

// The following classes should be nested inside the
// corresponding generators

class BarrierT;

/** Interface for Barrier Hit Time state variable */
class MCARLO_DLL IMCBarrierHTStateVar: public virtual IStateVariable {
public:
    virtual int hitTime(int iAsset) const = 0;
};

typedef smartPtr<IMCBarrierHTStateVar> IMCBarrierHTStateVarSP;
typedef smartConstPtr<IMCBarrierHTStateVar> IMCBarrierHTStateVarConstSP;

/** Interface that Path Generators implement if they support a
    barrier hit time state variable */
class MCARLO_DLL IMCBarrierHTGen {
public:
    /** Destructor */
    virtual ~IMCBarrierHTGen() {}

    /** create a barrier hit time state variable. Note, could have
        several methods here for different barriers */
    virtual IMCBarrierHTStateVar* createHitTimeSV(const BarrierT* barrier) = 0;
};


//////////////////////////////////////////////////////////////////////

    /** Interface for Barrier Hit Value and Time state variables */
class MCARLO_DLL IMCBarrierHVTimeStateVarVat: public virtual SVGenBarrierHV::IStateVar,
                                   public virtual IMCBarrierHTStateVar {
};

typedef smartPtr<IMCBarrierHVTimeStateVarVat> IMCBarrierHVTimeStateVarVatSP;
typedef smartConstPtr<IMCBarrierHVTimeStateVarVat> IMCBarrierHVTimeStateVarVatConstSP;

////////////////////////////////////////////////////////////////////////
// SVGenBarrierHVT class can generate a HitValue state variable and Time
////////////////////////////////////////////////////////////////////////
//-------------------//
// Interface class
//-------------------//
class SVGenBarrierHVT: virtual public IStateVariableGen{
public:
    /** European barrier implementation */
    class IStateVar: virtual public IStateVariable {
    public:
        /** return hit or not.  Also return hitValue and hit Time Step Index*/
        virtual bool hitValueTime(double &hitValue, int &iStep) const = 0;

        /** Records BARRIER_LEVEL output requests in results */
        virtual void recordBarrierLevels(CControl* control,
                                         Results*  results,
                                         const IMultiFactors* mAsset) const = 0;
    };

    typedef smartPtr<IStateVar> IStateVarSP;
    typedef smartConstPtr<IStateVar> IStateVarConstSP;

    /** Deconstuctor */
    virtual ~SVGenBarrierHVT(){};

    /** Implementation of IStateVariableGen */
    virtual IStateVarSP getHitValueTimeSV(IStateVarSP   oldStateVar,
                                  IStateVariableGen::IStateGen* pathGen) const = 0;

    virtual BarrierDataConstSP getBarrierData() const = 0;
};

DECLARE_REF_COUNT(SVGenBarrierHVT);

//---------------------------//
// only for simHit Euro Case.
//---------------------------//
class MCARLO_DLL SVGenBarrierHVTSimEur: virtual public SVGenBarrierHVT,
                                        virtual public IStateVariableClient {
public:
    /** European barrier implementation */
    class StateVar: virtual public SVGenBarrierHVT::IStateVar {
    public:
        /** Constructor */
        StateVar(BarrierDataConstSP data,
                 IRefLevel::IStateVarGenSP refLevelGen,
                 SVGenSpotSP spotMonGen);

        /** Allows the barrier class to switch from past to future */
        void update(IStateVariableGen::IStateGen* pathGen);

        /** return hit or not.  Also return hitValue and hit Time Step Index*/
        bool hitValueTime(double &hitValue, int &iStep) const;

        /** Indicates past Vs future */
        bool doingPast() const;

        /** Records BARRIER_LEVEL output requests in results */
        virtual void recordBarrierLevels(CControl* control,
                                         Results*  results,
                                         const IMultiFactors* mAsset) const;

    private:
        // Fields
        BarrierDataConstSP         data;              //!< Data
        IRefLevel::IStateVarGenSP  refLevelGen;       //!< Ref Level generator
        SVGenSpotSP                   spotMonGen;        //!< Spot generator for E type monitoring

        // Transient fields
        SVGenSpot::IStateVarSP        monitoringPath;    //!< Discrete path to monitor
        IRefLevel::IStateVarSP     refLevel;          //!< Reference level for barrier

        bool                       isPast;            //!< Whether we process the past or future

        mutable BoolArray          hitNoHit;          //!< Working area
        mutable DoubleArray        hitValues;      //!< Hit array for assets
        
        bool isStructOut;
        mutable bool isPastHit;     // for past management;
        mutable int  metAtStepPast; // for past managemetn;
    };

    typedef smartPtr<StateVar> StateVarSP;
    typedef smartConstPtr<StateVar> StateVarConstSP;

    /** Constuctor */
    SVGenBarrierHVTSimEur(BarrierDataSP data,
                    IRefLevel::IStateVarGenSP refLevelGen);

    /** Implementation of IStateVariableClient */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;
    
    /** Implementation of IStateVariableGen */
    IStateVarSP getHitValueTimeSV(IStateVarSP               oldStateVar,
                                  IStateVariableGen::IStateGen* pathGen) const;
    
    /** Part of IStateVariableGen */
    IStateVariableSP create(IStateVariableSP             oldStateVar,
                            IStateVariableGen::IStateGen* pathGen) const;

    /** Allows access to data */
    virtual BarrierDataConstSP getBarrierData() const;

    /** Allows access to ref level */
    IRefLevel::IStateVarGenSP  getRefLevel() const;

private:
    /** Disabled default constructor */
    SVGenBarrierHVTSimEur();

    // Fields
    BarrierDataSP             data;            //!< The barrier data
    IRefLevel::IStateVarGenSP refLevelGen;     //!< Ref Level generator

    // Transient fields
    SVGenSpotSP                  spotMonGen;      //!< Spot generator for E type monitoring
};

DECLARE_REF_COUNT(SVGenBarrierHVTSimEur);

DRLIB_END_NAMESPACE

#endif

