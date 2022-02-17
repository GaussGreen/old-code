//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfig.cpp
//
//   Description : 
//
//   Date        : July 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/LRGenerator.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/DependenceLocalCorr.hpp"
#include "edginc/QMCGenDiffusibleAsset.hpp"
#include "edginc/PastValues.hpp"

DRLIB_BEGIN_NAMESPACE
IMCPathConfig::~IMCPathConfig(){}
void IMCPathConfig::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IMCPathConfig, clazz);
    EXTENDS(IObject);
}
const int IMCPathConfig::PATH_RANDOM_ACCESS = 1;
const int IMCPathConfig::PATH_CACHE = 2;

CClassConstSP const IMCPathConfig::TYPE = CClass::registerInterfaceLoadMethod(
    "IMCPathConfig", typeid(IMCPathConfig), load);

MCPathConfig::~MCPathConfig(){}

MCPathConfig::MCPathConfig(CClassConstSP clazz) : 
    CObject(clazz), dependenceMaker(new DependenceMakerGauss()), skewMaker(0), momentMatching(false) { }

MCPathConfig::MCPathConfig(CClassConstSP clazz, DependenceMakerSP dependenceMaker) : 
    CObject(clazz), dependenceMaker(dependenceMaker), skewMaker(0), momentMatching(false) { }

    void MCPathConfig::validatePop2Object() {}

/** For shallow copy of randStates */
IObject* MCPathConfig::clone() const{
    MCPathConfig& myCopy = 
        dynamic_cast<MCPathConfig&>(*(CObject::clone()));
    myCopy.initRandState = initRandState; // shallow copy
    myCopy.pathCache = pathCache; // shallow copy
    myCopy.randomCache = randomCache; // shallow copy
    myCopy.idiosynFactorRandomCache = idiosynFactorRandomCache; // shallow copy
    myCopy.marketFactorRandomCache = marketFactorRandomCache; // shallow copy
    myCopy.momentMatching = momentMatching;
    return &myCopy;
}

void MCPathConfig::createUtil( 
    IQMCGenDiffusibleAsset* gen,
    const DateTime&        today,
    DateTimeArrayConstSP   simDates)
{
    // do nothing by default
}
DateTimeArray MCPathConfig::getCriticalDates(
    IQMCGenDiffusibleAsset* gen,
    const DateTime& start,       // likely to be "today"
    const DateTime& finish)      // the latest requested date
{
    // do nothing by default
    return DateTimeArray();
}

void MCPathConfig::setDiffusibleAsset(
    const DateTime&     today,              // base date of the run
    IQMCGenDiffusibleAsset* gen,       // to get max/min boundaries, etc
    const IPastValues*  pastValues,         // historic values
    DependenceSP        dependence)         // factorized correlations
{
    // do nothing by default
}

/** Creates a past path generator */
MCPathGeneratorSP MCPathConfig::pastPathGenerator(
    const IMCProduct* prod) {
    return MCPathGeneratorSP(new PastPathGenerator(prod));
}

/** Creates a PathGenerator used for future dates which supports
    'random access' to the paths. That is if a path generator is
    created using savePaths = true, then subsequently another path
    generator can be created (with savePaths = false) which will
    support the pathIdx parameter in the generatePath() method, ie the
    paths can be generated in any order. */
MCPathGeneratorSP MCPathConfig::futurePathGenerator(
    int                      cachingMode,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control, 
    Results*                 results,
    DateTimeArray&           simDates ) {
    bool isPricing = control->isPricing();
    int  myNumStates = (numPaths+1)/2;
    bool cachePaths = (cachingMode & PATH_CACHE)? true: false;
    // cachePaths => cacheRandoms, so no need to save states
    bool mySaveStates = !cachePaths && (cachingMode & PATH_RANDOM_ACCESS);
    /* no point saving paths if no greeks */
    if (mySaveStates){
        throw ModelException("MCPathConfig::futurePathGenerator", 
                             "QuickGreeks are not supported without path caching");
    }

    if (cachePaths){
        // caching paths: use new cache - find out number of future dates
        int numAssets = prod->getNumAssets();
        const DateTime& today = prod->getToday();
        const DateTime& simStartDate = prod->getEffectiveSimStartDate();
        const IRefLevel* refLevel = prod->getRefLevel();
        const SimSeries* simSeries = prod->getSimSeries();
        BoolArray storeRefLevels(numAssets); // only cache overall refLevel
        IntArray numDatesPerAsset(numAssets);
        for (int i = 0; i < numAssets; i++){
            numDatesPerAsset[i] = simSeries->numFutureDates(today, i);
            // note simulation starts from simStartDate so 1st refLevel
            // might be 'known'
            storeRefLevels[i] = refLevel->numFutureDates(simStartDate, i) > 0;
        }
        pathCache = PathCacheSP(
            new PathCache(numPaths, numDatesPerAsset, storeRefLevels));
        // then set up random cache
        randomCache = RandomCacheSP(new RandomCache(myNumStates, true, true));
        if (skewMaker.get()) {
            randomCache = RandomCacheSP(new RandomCache(myNumStates*2, true, false));
            idiosynFactorRandomCache = RandomCacheSP(new RandomCache(myNumStates, false, true));
            marketFactorRandomCache = RandomCacheSP(new RandomCache(myNumStates, false, true));
        } else {
            idiosynFactorRandomCache = RandomCacheSP(new RandomCache()); // dummy
            marketFactorRandomCache = RandomCacheSP(new RandomCache()); // dummy
        }
    } else if (!pathCache || isPricing){
        pathCache = PathCacheSP(new PathCache());
        randomCache = RandomCacheSP(new RandomCache());
        idiosynFactorRandomCache = RandomCacheSP(new RandomCache()); // dummy
        marketFactorRandomCache = RandomCacheSP(new RandomCache()); // dummy
    } else if (!pathCache->isDummy()){
        // invalidate paths for appropriate assets
        // Work out which assets which are being tweaked
        // Start by gettting hold of MultiFactors
        const IMultiFactors* multiFactor = prod->getMultiFactors();
        // then which greek we're doing
        SensitivitySP sens(control->getCurrentSensitivity());
        // getSensitiveAssets() does the real work
        IntArray sensitivePhiAssets(multiFactor->getSensitiveAssets(
                                        sens.get(), true)); // include phi etc
        if (sensitivePhiAssets.empty()){
            IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                                         sens.get(), false)); // exclude phi etc
            pathCache->configureCacheValidity(sensitiveAssets);
            randomCache->configure(true); // correlated randoms are ok
            // dont have correlated RN in idiosynFactorRandomCache or marketFactorRandomCache
        } else {
            // mark all paths as invalid
            pathCache->invalidateAllPaths();
            randomCache->configure(false); // correlated randoms are not ok
            // dont have correlated RN in idiosynFactorRandomCache or marketFactorRandomCache            
        }
    }
    // Bit comparison to deduce whether user has asked for caching
    bool cachingRequested = (cachingMode & IMCPathConfig::PATH_CACHE) ? true : false;
    
    return makePathGenerator(cachingRequested, numPaths, pastPathGenerator, 
                             prod, control, results, simDates);
}

/** Returns the cache of paths (which is preserved across tweaks). Only
    paths that are marked as valid can be used. Will not return null */
MCPathConfig::PathCacheSP MCPathConfig::getPathCache(){
    if (!pathCache){
        throw ModelException("MCPathConfig::getPathCache",
                             "Internal error: "
                             "futurePathGenerator not yet called");
    }
    return pathCache;
}

/** Returns the cache of asset random numbers (which is preserved across tweaks). 
    Clients must respect return value of isValid(). Will not return null */
MCPathConfig::RandomCacheSP MCPathConfig::getRandomCache(){
    if (!randomCache){
        throw ModelException("MCPathConfig::getRandomCache",
                             "Internal error: "
                             "futurePathGenerator not yet called");
    }
    return randomCache;
}

/** Returns the cache of idiosyn factors resp market factor random numbers (which is preserved across tweaks). 
        Clients must respect return value of isValid(). Will not return null */
MCPathConfig::RandomCacheSP MCPathConfig::getIdiosynFactorRandomCache(){
    if (!idiosynFactorRandomCache){
        throw ModelException("MCPathConfig::getIdiosynFactorRandomCache",
                             "Internal error: "
                             "futurePathGenerator not yet called");
    }
    return idiosynFactorRandomCache;
}

MCPathConfig::RandomCacheSP MCPathConfig::getMarketFactorRandomCache(){
    if (!marketFactorRandomCache){
        throw ModelException("MCPathConfig::getMarketFactorRandomCache",
                             "Internal error: "
                             "futurePathGenerator not yet called");
    }
    return marketFactorRandomCache;
}

/** Invoked after instrument has got its market data. Allows model to
    get any extra data required. Nothing by default */
void MCPathConfig::getMarket(const IModel*      model,
                             const MarketData*  market,
                             IInstrumentCollectionSP instrument){
    if (DependenceMakerGaussTerm::TYPE->isInstance(dependenceMaker)) {
        MarketObjectArraySP mo = market->GetAllDataWithType(CorrelationCategory::TYPE); 
        if (mo->size()==0) {
            // switch back since there are no correlation category objects
            dependenceMaker = DependenceMakerSP(new DependenceMakerGauss());
        }
    }
    if (skewMaker.get()) {
        MarketObjectArraySP mo = market->GetAllDataWithType(CorrelationCategory::TYPE); 
        if (mo->size()==0) {
            // switch back since there are no correlation category objects            
            skewMaker = SkewMakerSP(0);
        }
    }
}
    
/** Saves [internally] the current state of the random number generator.
    This state will be used to initialise the generator every time
    a [future] path generator is created. If no state has been saved
    init() should be used to initialise the generator. */
void MCPathConfig::saveRandomGeneratorState(){
    initRandState = IRandom::StateSP(rand->getState());
}

/** Returns the random generator pre-initialised and ready for use.
    Motivation: to split MC pricing into blocks need to initialise the
    generator's state correctly for each block */
IRandomSP MCPathConfig::getRandomGenerator(){
    if (!initRandState){
        rand->init();
    } else {
        rand->setState(initRandState.get());
    }
    return rand;
}

/** Returns the number of bytes used for random number storage per path.
    Implementation here assumes one double per asset per date. Do not
    invoke if there are no sim dates in future */
int MCPathConfig::randomStoragePerPath(IMCProduct* product) const{
    int numAssets = product->getNumAssets();
    const DateTime& today = product->getToday();
    const IRefLevel* refLevel = product->getRefLevel();
    const SimSeries* simSeries = product->getSimSeries();
    int numFutureRefLevelDates = refLevel->getFutureDates(today).size();
    int numFutureSimDates = simSeries->getFutureDates(today).size();

    if (skewMaker.get()) {   
        // onlye rough approximation -- at least one simdate per month, on average 2 factors
        int monthsUntilMaturity = (double)(simSeries->getLastDate().daysDiff(today)) / 30.0;        
        int addSimDatesUsed = Maths::max(monthsUntilMaturity, numFutureRefLevelDates + numFutureSimDates);

        int nbMarketFactors = numAssets > 10 ? 3 : 2;
        
        // we store both correlated and uncorrelated numbers but only for
        // every other path. Hence no times by 2.
        return (sizeof(double) * (
                numAssets * (numFutureRefLevelDates + numFutureSimDates) 
                + nbMarketFactors * addSimDatesUsed));

    } else {
        // we store both correlated and uncorrelated numbers but only for
        // every other path. Hence no times by 2.
        return (sizeof(double) *
            numAssets * (numFutureRefLevelDates + numFutureSimDates));
    }
}   

/** Returns the amount of storage space the path generator needs to
    save paths. It does not need to be exact. Implementation here
    assumes two doubles per date per asset */
int MCPathConfig::storagePerPath(IMCProduct* product) const{
    if (!product->hasFuture()){
        return 0; // nothing to simulate
    }
        
    int numAssets = product->getNumAssets();
    const DateTime& today = product->getToday();
    // First handle path cache
    const IRefLevel* refLevel = product->getRefLevel();
    const SimSeries* simSeries = product->getSimSeries();
    // need to include avg in + usual sim dates
    int numFutureDates = 0;
    for (int i = 0; i < numAssets; i++){
        if (refLevel->numFutureDates(today, i) > 0){
            numFutureDates++;
        }
        numFutureDates += simSeries->numFutureDates(today, i);
    }
    int numBytes = (sizeof(double) * numFutureDates);
    
    numBytes += randomStoragePerPath(product); // For Random Cache
    return numBytes;
}

/** Clears out appropriate caches under a theta shift */
bool MCPathConfig::sensShift(Theta* shift){
    pathCache = PathCacheSP();
    randomCache = RandomCacheSP();
    return true;
}

/** See comment on pure virtual declaration. This default implementation
    returns null */
LRGenerator* MCPathConfig::createLRGenerator(
        const MCPathGeneratorSP& pathGen,
        int                      nbIter,
        int                      nbSubSamples){
    return 0;
}

/** Invoked when Class is 'loaded' */
void MCPathConfig::load(CClassSP& clazz) {
    REGISTER(MCPathConfig, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMCPathConfig);
    IMPLEMENTS(Theta::Shift);
    FIELD(rand, "Random number generator");
    FIELD(dependenceMaker, "dependence maker");
    FIELD(skewMaker, "skew maker");
    FIELD(momentMatching, "moment matching");
    FIELD_MAKE_OPTIONAL(dependenceMaker);
    FIELD_MAKE_OPTIONAL(momentMatching);
    FIELD_MAKE_OPTIONAL(skewMaker);
}

CClassConstSP const MCPathConfig::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfig", typeid(MCPathConfig), load);

/** Creates cache ready to store paths */
MCPathConfig::PathCache::PathCache(int              numPaths, 
                                   const IntArray&  numDatesPerAsset,
                                   const BoolArray& storeRefLevels):
    paths(0), assetOffsets(numDatesPerAsset.size()+1),
    haveRefLevels(storeRefLevels),
    valid(numDatesPerAsset.size(), false), allowWrite(true){
    int numAssets = numDatesPerAsset.size();
    int totalNumDoubles = 0;
    for (int i = 0; i < numAssets; i++){
        totalNumDoubles += numDatesPerAsset[i];
        if (haveRefLevels[i]){
            totalNumDoubles++;
        }
    }
    pathOffset = totalNumDoubles; // offset for each path
    totalNumDoubles *= numPaths;
    if (totalNumDoubles > 0){
        paths = new double[totalNumDoubles];
        assetOffsets[0] = 0;
        for (int j = 1; j < numAssets+1; j++){
            assetOffsets[j] = assetOffsets[j-1] + numDatesPerAsset[j-1] + 
                (haveRefLevels[j-1]? 1: 0);
        }
    }
}

/** Creates empty cache. isCacheValid() and updateAllowed() will return false */
MCPathConfig::PathCache::PathCache(): paths(0), allowWrite(false) {}

//// returns true if the cache should be updated if isCacheValid = false
bool MCPathConfig::PathCache::updateAllowed(int iAsset) const{
    return allowWrite; // don't support per asset [yet]
}

/** Returns true if the cache has a valid path for this asset */
bool MCPathConfig::PathCache::isValid(int iAsset) const{
    return (!valid.empty() && valid[iAsset]);
}

/** Validates that the cache is the right size. Very much a defensive piece
    of programming but it is only called once per construction of path
    generator so very cheap. Derived classes of MCPathConfig should call
    this with the appropriate data: length is number of points in cache 
    excluding the refLevel. */
void MCPathConfig::PathCache::validateSize(
    int iAsset, int length, bool storeRefLevel) const{
    static const string method("MCPathConfig::PathCache::validateSize");
    if (isValid(iAsset) || updateAllowed(iAsset)){
        if (iAsset < 0 || iAsset >= assetOffsets.size()-1){
            // this really shouldn't happen
            throw ModelException(method, "Asset index out of bounds");
        }
        int assetOffset = assetOffsets[iAsset];
        int cacheLength = assetOffsets[iAsset+1] - assetOffset;
        int requestedLength = length + (storeRefLevel? 1: 0);
        if (requestedLength != cacheLength){
            // hopefully this won't happen
            throw ModelException(method, "Internal error: cache misconfigured");
        }
    }
}

/** Path stored in cache is written to dest. If cache has been
    configured to store refLevels for specified asset then refLevel is
    populated */
void MCPathConfig::PathCache::read(
    int iAsset, int pathIdx, double* dest, double& refLevel) const{
    int assetOffset = assetOffsets[iAsset];
    int length = assetOffsets[iAsset+1] - assetOffset;
    const double* path = paths + (pathOffset * pathIdx) + assetOffset;
    if (haveRefLevels[iAsset]){
        refLevel = *path;
        path++; // move over refLevel
        length--; // one less double to copy over
    }
    memcpy(dest, path, sizeof(double) * length);
}

/** Supplied path is written into cache with specified refLevel if cache has
    been configured to store refLevels for specified asset */
void MCPathConfig::PathCache::write(int iAsset, int pathIdx, 
                                    const double* src, double refLevel){
    int assetOffset = assetOffsets[iAsset];
    int length = assetOffsets[iAsset+1] - assetOffset;
    double* path = paths + (pathOffset * pathIdx) + assetOffset;
    if (haveRefLevels[iAsset]){
        *path = refLevel;
        path++; // move over refLevel
        length--; // one less double to copy over
    }
    memcpy(path, src, sizeof(double) * length);
}

/** The paths for each asset in the specified list are marked as invalid. Those
    that are not are marked as valid. */
void MCPathConfig::PathCache::configureCacheValidity(const IntArray& assets){
    // mark everything as being valid again
    valid = BoolArray(valid.size(), true);
    // and then switch off specified assets
    for (int i = 0; i < assets.size(); i++){
        valid[assets[i]] = false;
    }
    allowWrite = false;
}

/** Marks paths for all assets as invalid */
void MCPathConfig::PathCache::invalidateAllPaths(){
    valid = BoolArray(valid.size(), false);
    allowWrite = false;
}

/** Returns true if cache is just a dummy ie constructed using empty 
    constructor */
bool MCPathConfig::PathCache::isDummy() const{
    return (valid.empty());
}

MCPathConfig::PathCache::~PathCache(){
    delete[] paths;
}

MCPathConfig::RandomCache::~RandomCache(){
    delete[] correlatedRands;
    delete[] uncorrelatedRands;
}

/** Informs cache about number of assets and number of dates ie store
    numAssets * numDatesPerAsset for each simulated path */
void MCPathConfig::RandomCache::configure(
    int numAssets, int numDatesPerAsset){
    static const string method("MCPathConfig::RandomCache::configure");
    // only do anything if we actually caching anything
    if (storeCorrelatedRands || storeUncorrelatedRands){
        if (!numRandsPerPathSet){
            // first time in
            this->numAssets = numAssets;
            this->numDatesPerAsset = numDatesPerAsset;
            numRandsPerPath = numAssets * numDatesPerAsset;
            numRandsPerPathSet = true;
            if (numRandPaths > 0){
                int numDbls = numRandsPerPath * numRandPaths;
                QLIB_VERIFY (numDbls < 0x10000000, 
                    "Too many random numbers needed - try using less sample paths or switching off matchLegacySRM3 if using SRM");
                try{
                    if (!correlatedRands && storeCorrelatedRands){
                        correlatedRands = new double[numDbls];
                    }
                    if (!uncorrelatedRands && storeUncorrelatedRands){
                        uncorrelatedRands = new double[numDbls];
                    }
                } catch (exception& e){
                    throw ModelException(e, method, "Failed when trying to "
                                         "allocate space for "+
                                         Format::toString(numDbls)+" doubles "
                                         "for random number cache");
                }
            }
        } else {
            if (this->numAssets != numAssets ||
                this->numDatesPerAsset != numDatesPerAsset){
                throw ModelException(method,
                                     "Number of assets/dates has changed "
                                     "between greeks");
            }
        }
    }
}

/** Populates supplied DoubleMatrix with cached randoms for specified path.
    Either correlated or uncorrelated randoms depending on 
    correlatedRandsStored(). Results undefined if cache not written to
    or randPathIdx has value outside that specified in constructor or
    randoms not of sufficient size (as specified in configure) */
void MCPathConfig::RandomCache::read(
    int randPathIdx, bool correlated, DoubleMatrix& randoms) const{
    const double* rand = (correlated? correlatedRands: uncorrelatedRands) + 
        numRandsPerPath * randPathIdx;
    for (int i = 0; i < numAssets; i++, rand += numDatesPerAsset){
        memcpy(randoms[i], rand, sizeof(double) * numDatesPerAsset);
    }
}

/** Copies supplied DoubleMatrix into cache for specified path.
    Randoms should either be correlated or uncorrelated randoms
    depending on correlatedRandsStored(). configure() must be invoked
    before write can be invoked. Results undefined if
    randPathIdx has value outside that specified in constructor or
    randoms not of sufficient size (as specified in configure) */
void MCPathConfig::RandomCache::write(
    int randPathIdx, bool correlated, const DoubleMatrix& randoms){
    double* rand = (correlated? correlatedRands: uncorrelatedRands) + 
        numRandsPerPath * randPathIdx;
    for (int i = 0; i < numAssets; i++, rand += numDatesPerAsset){
        memcpy(rand, randoms[i], sizeof(double) * numDatesPerAsset);
    }
}
/** Low level write function - used MCRandom::blockFillByDate which
    generates randoms in SRM3 style (ie outer most loops across 
    dates) */
void MCPathConfig::RandomCache::write(
    int randPathIdx, int factorIdx, int dateIdx, 
    bool correlated, double random){
    double* rand = (correlated? correlatedRands: uncorrelatedRands) + 
        numRandsPerPath * randPathIdx;
    rand += numDatesPerAsset * factorIdx;
    rand[dateIdx] = random;
}

/** Are the contents of the cache valid */
bool MCPathConfig::RandomCache::isValid(bool correlated) const{
    return (correlated? validCorrelatedRands: validUncorrelatedRands);
}

/** Are we allowed to write correlatedRands (true) /
    uncorrelated (false) randoms to the cache */
bool MCPathConfig::RandomCache::updateAllowed(bool correlated) const{
    return (allowWrite && (correlated? 
                           storeCorrelatedRands: storeUncorrelatedRands));
}

/** Returns the number of paths for which this cache stores random numbers */
int MCPathConfig::RandomCache::getNumRandPaths() const{
    return numRandPaths;
}

/** Create dummy cache that will not store anything */
MCPathConfig::RandomCache::RandomCache(): 
    correlatedRands(0), uncorrelatedRands(0), 
    numRandsPerPath(0), numRandPaths(0), numAssets(0), numDatesPerAsset(0),
    validCorrelatedRands(false), validUncorrelatedRands(false), 
    allowWrite(false), storeCorrelatedRands(false),
    storeUncorrelatedRands(false), numRandsPerPathSet(false){}

/** Create cache for storing either correlated or non correlated randoms*/
MCPathConfig::RandomCache::RandomCache(
    int numRandPaths, bool storeCorrelatedRands, bool storeUncorrelatedRands):
    correlatedRands(0), uncorrelatedRands(0), 
    numRandsPerPath(0), numRandPaths(numRandPaths), numAssets(0),
    numDatesPerAsset(0), 
    validCorrelatedRands(false), validUncorrelatedRands(false), 
    allowWrite(true), storeCorrelatedRands(storeCorrelatedRands),
    storeUncorrelatedRands(storeUncorrelatedRands), numRandsPerPathSet(false){}

/** Switches allowWrite off, sets validUncorrelatedRands if stored and
    sets validCorrelatedRands to supplied value */
void MCPathConfig::RandomCache::configure(bool correlatedRandomsValid)
{
    allowWrite = false;
    validUncorrelatedRands = storeUncorrelatedRands;
    validCorrelatedRands = correlatedRandomsValid;
}

DRLIB_END_NAMESPACE


