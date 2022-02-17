//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IMCPathConfig.hpp
//
//   Description : Interface for Monte Carlo path generator factory
//
//   Date        : June 2001
//
//
//----------------------------------------------------------------------------
#ifndef EDR_IMCPATHGENERATOR_HPP
#define EDR_IMCPATHGENERATOR_HPP
#include "edginc/Class.hpp"
#include "edginc/MCPathGenerator.hpp"
#include "edginc/Random.hpp"
#include "edginc/LRGenerator.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Theta.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class IMCProduct;
class IQMCGenDiffusibleAsset;
class IPastValues;

/** This interface is what MonteCarlo requires of a PathConfig. 
    Specifically it requires it to be able to build an IPathGenerator - 
    One for dealing with historic dates, and one for future dates. The
    interface uses smart pointers to allow complete flexibility to classes
    implementing this interface. For example an IMCPathConfig could
    implement the IPathGenerator interface as well or the pastPathGenerator
    could be the same object as the futurePathGenerator */
class MCARLO_DLL IMCPathConfig : public virtual IObject {
public:
    /** These ints are used in a bitwise manner in futurePathGenerator to
        indicate what caching is needed */
    static const int PATH_RANDOM_ACCESS;
    static const int PATH_CACHE;

    static CClassConstSP const TYPE;
    
    virtual ~IMCPathConfig();

    /** Returns the amount of storage space the path generator needs to
        save paths (ie if PATH_CACHE bit is set in the call to 
        futurePathGenerator()). It does not need to be exact. */
    virtual int storagePerPath(IMCProduct* product) const = 0;

    /** Creates a PathGenerator used for historic dates */
    virtual MCPathGeneratorSP pastPathGenerator(
        const IMCProduct*) = 0;

    /** Creates a PathGenerator used for future dates. The cachingMode
        parameter indicates what should be cached between tweaks. For example,
        if a path generator is
        created using with "cachingMode & PATH_RANDOM_ACCESS = true", then 
        subsequently another path generator can be created (with
        "cachingMode & PATH_RANDOM_ACCESS = false") which will
        support the pathIdx parameter in the generatePath() method, ie the
        paths can be generated in any order. */
    virtual MCPathGeneratorSP futurePathGenerator(
        int                                cachingMode,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates ) = 0; 


    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
     *  for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const = 0;

    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required. */
    virtual void getMarket(const IModel*      model,
                           const MarketData*  market,
                           IInstrumentCollectionSP instrument) = 0;

    /** For some vol models this may not make sense */
    virtual bool vegaMatrixSupported() const = 0;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const = 0;

    /** Saves [internally] the current state of the random number
        generator.  This state will be used to initialise the
        generator every time a [future] path generator is created. If
        no state has been saved init() should be used to initialise
        the generator.  (It might be cleaner if the state was stored
        somewhere else but this is ok for now) */
    virtual void saveRandomGeneratorState() = 0;

    /** Creates a LRGenerator for the specified sensitivity . The
        pathGen supplied must be that returned from the
        futurePathGenerator() method. If a LRGenerator is not
        supported for this PathConfig then null should be returned */
    virtual LRGenerator* createLRGenerator(
        const MCPathGeneratorSP& pathGen,
        int                      nbIter,
        int                      nbSubSamples) = 0;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

typedef smartPtr<IMCPathConfig> IMCPathConfigSP;
class MCARLO_DLL MCPathConfig: public CObject,
                    public virtual IMCPathConfig,
                    public virtual Theta::Shift {
public:
    static CClassConstSP const TYPE;
    /** Utility class for caching paths between tweaks. For use by derived
        classes */
    class MCARLO_DLL PathCache{
    private:
        /* for performance, cache stored as one big array of doubles
           See write for memory layout */
        double*                          paths;
        int                              pathOffset; // num dbls per path
        IntArray                         assetOffsets; // offsets for each asset
        BoolArray                        haveRefLevels; // per asset
        BoolArray                        valid; // [asset]
        bool                             allowWrite; // allow updates
        // friend class MCPathConfig;
    public:
        ~PathCache();
        /** Validates that the cache is the right size. Very much a
            defensive piece of programming but it is only called once per
            construction of path generator so very cheap. Derived classes of
            MCPathConfig should call this with the appropriate data: length is
            number of points in cache excluding the refLevel. */
        void validateSize(int iAsset, int length, bool storeRefLevel) const;
        
        /** Path stored in cache is written to dest. Additinally stored
         * refLevel is returned if stored */
        void read(int iAsset, int pathIdx, double* dest, double& refLevel)const;
        /** Supplied path is written into cache with specified refLevel if
            cache has been configured to store refLevels */
        void write(int iAsset, int pathIdx, const double* src, double refLevel);

        /** Returns true if the cache has a valid path for this asset */
        bool isValid(int iAsset) const;
        //// returns true if the cache should be updated if isCacheValid = false
        bool updateAllowed(int iAsset) const;
    
        PathCache(int              numPaths, 
                  const IntArray&  numDatesPerAsset,
                  const BoolArray& storeRefLevels);
        PathCache();
        const double* internalRead(int iAsset, int pathIdx, double* dest) const;
        void configureCacheValidity(const IntArray& assets);
        void invalidateAllPaths();
        bool isDummy() const;
    };
    typedef refCountPtr<PathCache> PathCacheSP;
    
    class MCARLO_DLL RandomCache{
    private:
        // friend class MCPathConfig;
        double*   correlatedRands; // in one long line
        double*   uncorrelatedRands; // in one long line
        int       numRandsPerPath;
        int       numRandPaths;
        int       numAssets;
        int       numDatesPerAsset;
        bool      validCorrelatedRands;
        bool      validUncorrelatedRands;
        bool      allowWrite; // allow updates
        bool      storeCorrelatedRands;
        bool      storeUncorrelatedRands;
        bool      numRandsPerPathSet;
    public:
        ~RandomCache();
        /** Informs cache about number of assets and number of dates ie store
            numAssets * numDatesPerAsset for each simulated path */
        void configure(int numAssets, int numDates);
        /** Populates supplied DoubleMatrix with cached randoms for
            specified path.  Either correlated or uncorrelated
            randoms depending on 'correlated'. Results
            undefined if cache not written to or randPathIdx has value
            outside that specified in constructor or randoms not of
            sufficient size (as specified in configure) */
        void read(int randPathIdx, bool correlated, DoubleMatrix& randoms)const;
        /** Copies supplied DoubleMatrix into cache for specified path.
            Randoms should either be correlated or uncorrelated randoms
            depending on 'correlated'. configure() must be invoked
            before write can be invoked. Results undefined if
            randPathIdx has value outside that specified in constructor or
            randoms not of sufficient size (as specified in configure) */
        void write(int randPathIdx,bool correlated,const DoubleMatrix& randoms);
        /** Low level write function - used MCRandom::blockFillByDate which
            generates randoms in SRM3 style (ie outer most loops across 
            dates) */
        void write(int randPathIdx, int factorIdx, int dateIdx, bool correlated,
                   double random); 
        /** Are the contents of the cache valid */
        bool isValid(bool correlated) const;
        /** Are we allowed to write correlatedRands (true) /
            uncorrelated (false) randoms to the cache */
        bool updateAllowed(bool correlated) const;
        /** Returns the number of paths for which this cache stores random
            numbers */
        int getNumRandPaths() const;

        RandomCache(int  numRandPaths, 
                    bool storeNonCorrelatedRandoms, 
                    bool storeUncorrelatedRands);
        RandomCache();
        void configure(bool correlatedRandomsValid);
    };
    typedef refCountPtr<RandomCache> RandomCacheSP;

    ~MCPathConfig();
    MCPathConfig(CClassConstSP clazz);
    MCPathConfig(CClassConstSP clazz, DependenceMakerSP dependenceMaker);

    virtual void validatePop2Object();
    
    /** Returns the amount of storage space the path generator needs
        to save paths (ie if PATH_CACHE bit is set in the call to
        futurePathGenerator()). It does not need to be
        exact. Implementation here assumes two doubles per date per
        asset */
    virtual int storagePerPath(IMCProduct* product) const;

    /** Override default clone to do shallow copy of 'randStates' */
    IObject* clone() const;

    /** Creates a PathGenerator used for historic dates */
    virtual MCPathGeneratorSP pastPathGenerator(
        const IMCProduct*);

    /** See comment on pure virtual declaration */
    virtual MCPathGeneratorSP futurePathGenerator(
        int                                cachingMode,
        int                                numPaths,
        const MCPathGeneratorSP& pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates);

    /** Invoked after instrument has got its market data. Allows model to
        get any extra data required. */
    virtual void getMarket(const IModel* model,
                           const MarketData* market,
                           IInstrumentCollectionSP instrument);

    /** Saves [internally] the current state of the random number generator.
        This state will be used to initialise the generator every time
        a [future] path generator is created. If no state has been saved
        init() should be used to initialise the generator. */
    virtual void saveRandomGeneratorState();

    /** Returns the random generator pre-initialised and ready for use.
        Motivation: to split MC pricing into blocks need to initialise the
        generator's state correctly for each block.
        Note this returns a reference to the one inside */
    IRandomSP getRandomGenerator();

    /** Returns the cache of paths (which is preserved across tweaks). Only
        paths that are marked as valid can be used. Will not return null */
    PathCacheSP getPathCache();

    /** Returns the cache of asset random numbers (which is preserved across tweaks). 
        Clients must respect return value of isValid(). Will not return null */
    RandomCacheSP getRandomCache();         // asset random cache

    /** Returns the cache of idiosyn factors resp market factor random numbers (which is preserved across tweaks). 
        Clients must respect return value of isValid(). Will not return null */
    RandomCacheSP getIdiosynFactorRandomCache();
    RandomCacheSP getMarketFactorRandomCache();

    /** Returns the number of bytes used for random number storage per path.
        Implementation here assumes one double per asset per date. Do not
        invoke if there are no sim dates in future */
    virtual int randomStoragePerPath(IMCProduct* product) const;

    /** whether we're doing 'isCarefulRandoms' ie consistent random numbers
        when doing theta */
    virtual bool carefulRandoms() const = 0;

    /** Clears out appropriate caches under a theta shift */
    virtual bool sensShift(Theta* shift);

    /** See comment on pure virtual declaration. This default implementation
        returns null */
    virtual LRGenerator* createLRGenerator(
        const MCPathGeneratorSP& pathGen,
        int                      nbIter,
        int                      nbSubSamples);

    bool getMomentMatchingFlag() const {return momentMatching;}

    virtual void createUtil( 
        IQMCGenDiffusibleAsset* gen, 
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    virtual DateTimeArray getCriticalDates(
        IQMCGenDiffusibleAsset* gen,
        const DateTime& start,       // likely to be "today"
        const DateTime& finish); // the latest requested date

    virtual void setDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        IQMCGenDiffusibleAsset* gen,       // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);    // factorized correlations

protected:
    /** Essentially a pass through for futurePathGenerator except that the
        relevant caches are created/updated etc */
    virtual MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates) = 0;

public:
    DependenceMakerSP getDependenceMaker() const { return dependenceMaker; }
    SkewMakerSP       getSkewMaker() const { return skewMaker; }

protected:
    IRandomSP                     rand;
    DependenceMakerSP             dependenceMaker;
    SkewMakerSP                  skewMaker;
    IRandom::StateSP              initRandState; // transient, not registered $unregistered
    PathCacheSP                   pathCache; // $unregistered
    RandomCacheSP                 randomCache; // $unregistered
    RandomCacheSP                 idiosynFactorRandomCache; // $unregistered
    RandomCacheSP                 marketFactorRandomCache; // $unregistered
    bool                          momentMatching;
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<MCPathConfig> MCPathConfigConstSP;
typedef smartPtr<MCPathConfig> MCPathConfigSP;

DRLIB_END_NAMESPACE
#endif
