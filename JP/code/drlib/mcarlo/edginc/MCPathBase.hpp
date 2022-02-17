//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathGenerator.cpp
//
//   Description : Monte Carlo base classes
//                 for future path generators
//
//   Date        : March 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_MCPATHBASE_HPP
#define EDR_MCPATHBASE_HPP
#include "edginc/config.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/RefLevel.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

// Whether to use a macro for mapping the path or no
#define MAPPATH_MACRO


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_REF_COUNT(MCPathBase);

/** Produces a simulation timeline given a SimSeries and a RefLevel.
    Also produces information on past/future dates etc. */
class MCARLO_DLL MCProductTimeline :
                public virtual VirtualDestructorBase
{
public:

    /** Destructor */
    ~MCProductTimeline () {}

    const SimSeries*      simSeries;      //!< Holds simulation dates (not RefLevel dates)
    const IRefLevel*      refLevelObj;    //!< RefLevel object

    /** future steps in RefLevel are relative to simStart date and not today
        This is mainly for the situation when you have simStart date =
        fwd start date and we want to avoid calculating random numbers for
        a zero time interval (or, more to the point, match the EDG library) */
    DateTime today;              //!< Today
    DateTime simStartDate;       //!< Simulation start date
    bool     simHasPast;         //!< Indicates whether product has past

    int      numFutRefLevel;     //!< Number of future dates in ref level
    int      numFutSteps;        //!< Includes future dates in refLevel
    int      totalNumSteps;      //!< numPastDates + numFutSteps
    int      numPastDates;       //!< Excludes past dates in refLevel
    int      totalNumPastDates;  //!< Includes historic ref levl
    int      offsetToClientPath; //!< Typically numFutRefLevel unless any overlap

    DateTimeArray futureDates;   //!< Set of future dates

private:
    friend class MCPathBase;

    /** Disable default constructor */
    MCProductTimeline();

    /** Constructor */
    MCProductTimeline(const IMCProduct* prod,
                      const MCPathGeneratorSP& pastPathGenerator);

};

typedef smartConstPtr<MCProductTimeline> MCProductTimelineConstSP;
typedef smartPtr<MCProductTimeline> MCProductTimelineSP;


/** Creates elementary reference level data for Monte Carlo path generators
    e.g. fwds, fwdAtSimStart, reference level etc. */
class MCARLO_DLL RefLevelData : public virtual VirtualDestructorBase
{
public:

    /** Destructor */
    ~RefLevelData() {}

    vector<DoubleArray>   refLevels;      //!< [iAsset][iPath]
    DoubleArray           fwdsAtSimStart; //!< By asset
    IRefLevel::IMCPathSP  refLevelPath;   //!< Works out refLevels

private:
    friend class MCPathBase;

    /** Disabled default constructor */
    RefLevelData();

    /** Constructor where fwds are produced at product dates only */
    RefLevelData(const MCProductTimelineSP& timeline,
                 const IntArray& nbPaths,
                 const IMultiFactors* mAsset,
                 const IMCProduct* prod,
                 const MCPathGeneratorSP& pastPathGenerator);

};

typedef smartConstPtr<RefLevelData> RefLevelDataConstSP;
typedef smartPtr<RefLevelData> RefLevelDataSP;


/** Interface class that is wrapped by the 2 path generators. */
class MCARLO_DLL MCPathBase{
public:
    virtual ~MCPathBase() {};

    class MCARLO_DLL Generator: virtual public MCPathGenerator{
    public:
        ~Generator();
        /** Returns the path base with which the Generator was built */
        virtual MCPathBaseSP getPathBase() = 0;
    };
    typedef refCountPtr<Generator> GeneratorSP;

    /** Creates a [future] path generator using supplied MCPathBase.
        Essentially a wrapper mapping requirements of MCPathGenerator onto
        properties of MCPathBase (eg handles different dates per
        asset etc) */
    static GeneratorSP createPathGenerator(
        const MCPathGeneratorSP& pastPathGenerator,
        MCPathConfig*            pathConfig,
        const IMCProduct*         prod,
        MCPathBaseSP             basePathGen);

    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0.  */
    virtual double maxDriftProduct(int iAsset) const = 0;

    /** Obtains timeline object from base */
    virtual MCProductTimelineConstSP getTimeline() const = 0;

    /** Returns number of assets */
    virtual int NbSimAssets() const = 0;

    /** Returns the reference level for iAsset, iPath (write access) */
    virtual double& refLevel(int iAsset, int iPath) = 0;

    /** Returns the reflevel path */
    virtual IRefLevel::IMCPathSP& refLevelPath() const = 0;

    /** Returns the paths for iAsset, iPath (write access) */
    virtual double* Path(int iAsset, int iPath) = 0;

    /** Returns the number of paths per asset */
    virtual int nbPaths(int iAsset) const = 0;

    /** Returns if it is a single path generator or not */
    virtual bool isSinglePath() const = 0;

    /** Generate path for specified path in simulation (pathIdx), for
        specified asset (iAsset) and specified vol interp (iPath) */
    virtual void generatePath(int pathIdx, int iAsset, int iPath) = 0;

    /** Generates random numbers. Antithetics is handled in here. Called
        once on each run. */
    virtual void drawRandomNumbers(int pathIdx) = 0;

#ifdef MAPPATH_MACRO
    class MCSimpleFuturePathGenerator;
    class MCFuturePathGenerator;
#endif

protected:
    /** Produces the timeline information. Constructor of
        MCProductTimeline is private and only this method can access it
        via friendship */
    static MCProductTimelineSP getProductTimeline(const IMCProduct* prod,
        const MCPathGeneratorSP& pastPathGenerator);

    /** Create market data */
    static RefLevelDataSP getRefData(
        const MCProductTimelineSP& timeline,
        const IntArray& nbPaths,
        const IMultiFactors* mAsset,
        const IMCProduct* prod,
        const MCPathGeneratorSP& pastPathGenerator);

#ifndef MAPPATH_MACRO
private:
    class MCSimpleFuturePathGenerator;
    class MCFuturePathGenerator;
#endif
};

class MCARLO_DLL MCGammaAdjustedPathGenerator {
public:
	virtual const DoubleMatrix& getVols(int iAsset) const = 0;
	virtual const DoubleMatrix& getDrifts(int iAsset) const = 0;
	virtual const DoubleMatrix& getRandoms(int iAsset) const = 0;
};


DRLIB_END_NAMESPACE
#endif
