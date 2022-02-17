//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathBase.cpp
//
//   Description : MonteCarlo::IPathGenerator Methods for
//                 MCSimpleFuturePathGenerator and
//                 MCFuturePathGenerator
//
//
//   Date        : June 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/DECLARE.hpp"

#define CAREFUL_RAND_DEBUG 0
#if CAREFUL_RAND_DEBUG
#endif


#ifdef MAPPATH_MACRO

#define MAPPATH(isSinglePath, iPath) \
    isSinglePath ? 0 : iPath;

#endif

DRLIB_BEGIN_NAMESPACE


MCProductTimeline::MCProductTimeline(const IMCProduct* prod,
                                     const MCPathGeneratorSP& pastPathGenerator):
simSeries(prod->getSimSeries()), refLevelObj(prod->getRefLevel()) {
    static const string routine = "MCProductTimeline::MCProductTimeline";
    try {
        today = prod->getToday();
        simStartDate = prod->getEffectiveSimStartDate();

        simHasPast = pastPathGenerator->hasPast();

        // set up the variables in this structure
        DateTimeArray futureRefLevelDates =
            refLevelObj->getFutureDates(simStartDate /* NB not today */);
        DateTimeArray futurePathDates = simSeries->getFutureDates(today);
        int coincidence = 0;
        if (!futureRefLevelDates.empty()) {
            if (!futurePathDates.empty()) {
                if (futureRefLevelDates.back().isGreater(futurePathDates.front())) {
                    throw ModelException(routine, "Overlapping reference level and"
                                         " simulation dates");
                } else if (futureRefLevelDates.back()==futurePathDates.front()) {
                    coincidence = 1;
                }
            }
        }

        // Ensure that there are no dates between today and simStartDate
        int iStep;
        for(iStep = 0; iStep < futureRefLevelDates.size(); iStep++) {
            const DateTime& thisDate = futureRefLevelDates[iStep];
            if(today <= thisDate && thisDate < simStartDate) {
                throw ModelException(
                    "Future reference level date " +
                    thisDate.toString() +
                    " is before simulation start date " +
                    simStartDate.toString());
            }
        }
        for(iStep = 0; iStep < futurePathDates.size(); iStep++) {
            const DateTime& thisDate = futurePathDates[iStep];
            if(today <= thisDate && thisDate < simStartDate) {
                throw ModelException(
                    "Future simulation date " +
                    thisDate.toString() +
                    " is before simulation start date " +
                    simStartDate.toString());
            }
        }

        // then combine lists of dates together with simulation start
        // date
        numFutRefLevel = futureRefLevelDates.size();
        numFutSteps = numFutRefLevel + futurePathDates.size() - coincidence;

        totalNumSteps = simSeries->getAllDates().size() + numFutRefLevel - coincidence;
        numPastDates = totalNumSteps - numFutSteps; // excludes ref level
        // totalNumPastDates includes historic ref level dates
        const DateTimeArray& allRefLevelDates = refLevelObj->getAllDates();
        totalNumPastDates = allRefLevelDates.size() - 
            numFutRefLevel + numPastDates;
#if 0
        // This is #'d out since the effect is benign and we wish to avoid a massive (within-STE
        // jump in many test files). The code is here to allow easy confirmation of the cause 
        // of STE-sized diffs for relevant cases (e.g. spi/Coupons09).
        // Also account for overlap in past dates - deal only with the case caught above
        // since otherwise inst would not have survived to have past dates
        if (numPastDates>0) {
            const DateTimeArray& allSimDates = simSeries->getAllDates();
            if (allRefLevelDates.back() == allSimDates.front()) {
                // we have a simple overlap, so don't count it. I think this
                // affects random number order only.
                totalNumPastDates--;
            }
        }
#endif
        offsetToClientPath = numFutRefLevel - coincidence; // don't double count

        futureDates = DateTimeArray(1, simStartDate);
        futureDates.reserve(1+ numFutSteps);
        futureDates.insert(futureDates.end(),
                           futureRefLevelDates.begin(),
                           futureRefLevelDates.end());
        if (coincidence) {
            // don't double count any overlap - remove the first futurePathDates which
            // is already in futureDates via futureRefLevelDates.end()
            futurePathDates.erase(futurePathDates.begin());
            // used to do this but the fwd start case only works because there IS a duplication
            // at the start. It's all wrong since that is strictly a PAST date, but trying to do
            // minimal change here...
            //        DateTime::removeDuplicates(futureDates,
            //                                   false);
        }
        futureDates.insert(futureDates.end(),
                           futurePathDates.begin(),
                           futurePathDates.end());
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


////////////////////////////////////////////////////////////////

RefLevelData::RefLevelData(const MCProductTimelineSP& timeline,
                           const IntArray& nbPaths,
                           const IMultiFactors* mAsset,
                           const IMCProduct* prod,
                           const MCPathGeneratorSP& pastPathGenerator)
{
    // Create dimensions
    int numAssets = nbPaths.size();
    refLevels.resize(numAssets);
    fwdsAtSimStart.resize(numAssets);

    // Obtain fwds
    int iAsset;
    for (iAsset=0; iAsset < numAssets; iAsset++) {
        refLevels[iAsset] = DoubleArray(nbPaths[iAsset]);

        // save fwd at sim start
        fwdsAtSimStart[iAsset] = mAsset->factorFwdValue(iAsset, timeline->simStartDate);
    }

    refLevelPath = IRefLevel::IMCPathSP(
        timeline->refLevelObj->createMCPath(timeline->today,
                                            fwdsAtSimStart,
                                            prod->getMCPastValues()));

    // populate refLevel using pastPathGenerator. Note that generatePath
    // will correctly populate future average in dates - however part
    // of the specification is that until generatePath is called we should
    // return the same values as the past path generator
    for (iAsset=0; iAsset < numAssets; iAsset++) {
        /* next line should make you think. The issue is that the
           numFutRefLevel is wrt simulation start date. So for fwd starting
           the reference level is 'known' - it is the fwd (since the
           simulation starts at sim date). However, the past path generator
           uses the spot for the reference level. In generatePath() the
           refLevels are only populated if numFutRefLevel > 0 */
        double level = timeline->numFutRefLevel == 0?
            refLevelPath->refLevel(iAsset, (double *)0):
            pastPathGenerator->refLevel(iAsset, 0);
        for (int iPath = 0; iPath < nbPaths[iAsset]; iPath++){
            refLevels[iAsset][iPath] = level;
        }
    }
}


////////////////////////////////////////////////////////////////


/*******************************
MCPathBase
********************************/

MCProductTimelineSP MCPathBase::getProductTimeline(const IMCProduct* prod,
    const MCPathGeneratorSP& pastPathGenerator){
    return MCProductTimelineSP(new MCProductTimeline(
        prod, pastPathGenerator));
}


RefLevelDataSP MCPathBase::getRefData(
        const MCProductTimelineSP& timeline,
        const IntArray& nbPaths,
        const IMultiFactors* mAsset,
        const IMCProduct* prod,
        const MCPathGeneratorSP& pastPathGenerator) {

    return RefLevelDataSP (new RefLevelData(timeline, nbPaths, mAsset,
        prod, pastPathGenerator));
}


MCPathBase::Generator::~Generator(){}

/** Non virtual implementation of MCPathBase. Caches pointers
    to MCPathBase member fields to avoid some virtual calls.
    This speeds up all accessor functions. Only necessary
    virtual calls are propagated. The class is used by
    MCSimpleFuturePathGenerator and MCFuturePathGenerator.

    This is only a temporary solution around the problem that
    MCPathBase is not itself a MCPathGenerator. */
class MCPathBaseWrapper {
public:
    /** Constructor */
    MCPathBaseWrapper(const MCPathBaseSP& basePathGen):
    basePathGen(basePathGen), numAssets(basePathGen->NbSimAssets()),
    isSinglePath(basePathGen->isSinglePath()), timeline(basePathGen->getTimeline()),
    refLevelPath(basePathGen->refLevelPath()) {
        // Initialize nbPaths, paths and reflevels
        nbPaths   = IntArray(numAssets);
        paths     = vector<vector<double*> >(numAssets);
        refLevels = vector<vector<double*> >(numAssets);

        for (int iAsset=0; iAsset < numAssets; iAsset++) {
            nbPaths[iAsset]   = basePathGen->nbPaths(iAsset);
            paths[iAsset]     = vector<double*>(nbPaths[iAsset]);
            refLevels[iAsset] = vector<double*>(nbPaths[iAsset]);
            for (int iPath = 0; iPath < nbPaths[iAsset]; iPath++) {
                // Fill in the pointers
                paths[iAsset][iPath] = basePathGen->Path(iAsset, iPath);
                refLevels[iAsset][iPath] = &basePathGen->refLevel(iAsset, iPath);
            }
        }
    }

    /** Destructor */
    ~MCPathBaseWrapper() {}

#ifndef MAPPATH_MACRO
    // Helper function
    inline int mapPath(int iPath) const { return isSinglePath ? 0 : iPath; }
#endif

    // Accessors
    inline MCPathBaseSP getPathBase() const { return basePathGen; }
    inline int NbSimAssets() const { return numAssets; }
    inline const MCProductTimelineConstSP& getTimeline() const { return timeline; }
    inline const IRefLevel::IMCPathSP& getRefLevelPath() const { return refLevelPath; }
    inline int getNbPaths(int iAsset) const { return nbPaths[iAsset]; }

#ifdef MAPPATH_MACRO
    friend class MCPathBase::MCFuturePathGenerator;

    inline double* path(int iAsset, int iPath) const {
        int mappedPath = MAPPATH(isSinglePath, iPath);
        return paths[iAsset][mappedPath];
    }
    inline double& refLevel(int iAsset, int iPath) const {
        int mappedPath = MAPPATH(isSinglePath, iPath);
        return *refLevels[iAsset][mappedPath];
    }
#else
    inline double* path(int iAsset, int iPath) const {
        return paths[iAsset][mapPath(iPath)];
    }
    inline double& refLevel(int iAsset, int iPath) const {
        return *refLevels[iAsset][mapPath(iPath)];
    }
#endif

    // Propagating functions
    inline double maxDriftProduct(int iAsset) const {
        return basePathGen->maxDriftProduct(iAsset);
    }
    inline void generatePath(int pathIdx, int iAsset, int iPath) {
        basePathGen->generatePath(pathIdx, iAsset, iPath);
    }
    inline void drawRandomNumbers(int pathIdx) {
        basePathGen->drawRandomNumbers(pathIdx);
    };
private:
    MCPathBaseSP             basePathGen;   //!< Base PathGen

    int                      numAssets;     //!< Number of assets
    bool                     isSinglePath;  //!< isSinglePath generator or not
    MCProductTimelineConstSP timeline;      //!< Product timeline
    IRefLevel::IMCPathSP     refLevelPath;  //!< RefLevel path
    IntArray                 nbPaths;       //!< Number of paths [iAsset]
    vector<vector<double*> > paths;         //!< Pointers to path[iAsset][iPath][iStep]
    vector<vector<double*> > refLevels;     //!< Pointers to refLevel[iAsset][iPath]
};
DECLARE_REF_COUNT(MCPathBaseWrapper);

/** Class for generating future paths when all assets have the same set
    of simulation dates */
class MCPathBase::MCSimpleFuturePathGenerator: public Generator{
private:
    int                       beginIdx;   /* index into paths */
    int                       endIdx;     /* index into paths */
    int                       nowPathIdx; // Current path idx
    MCPathConfig::PathCacheSP pathCache;  // cache of generated paths
    MCPathBaseWrapper         wrapper;    // Quick access to MCPathBase fields
public:
    /** Returns the path base with which the Generator was built */
    virtual MCPathBaseSP getPathBase(){
        return wrapper.getPathBase();
    }

    virtual int NbSimAssets() const {
        return wrapper.NbSimAssets();
    }

    virtual double refLevel(int iAsset, int iPath) const {
        return wrapper.refLevel(iAsset, iPath);
    }

    /** Returns the sum of the absolute value of the 'drifts'
        between each of the simulation dates being simulated (ie
        between begin() and end()) (where the first ratio is between
        the first simulation date and the simulation start date). Note
        past path generator will return 0 */
    virtual double maxDriftProduct(int iAsset) const{
        return wrapper.maxDriftProduct(iAsset);
    }

    virtual bool hasPast() const {
        return wrapper.getTimeline()->simHasPast;
    }

    virtual bool doingPast() const {
        return false;
    }

    virtual const double* Path(int iAsset, int iPath) const{
        /* client expects array to be indexed such that the ref level dates
           are excluded */
        int offset = wrapper.getTimeline()->offsetToClientPath;
        return wrapper.path(iAsset, iPath) + offset;
    }

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    virtual int begin(int iAsset) const{
        return beginIdx; // independent of iAsset
    }

    //// see  begin(int iAsset)
    virtual int end(int iAsset) const{
        return endIdx;  // independent of iAsset
    }

    virtual void generatePath(int pathIdx){
        const MCProductTimelineConstSP& timeline = wrapper.getTimeline();

        // cache the pathIdx for getPathIndex()
        nowPathIdx = pathIdx;

        // get hold of the correlated random numbers
        wrapper.drawRandomNumbers(pathIdx);
        // for ease, read some variables
        bool futureRefLevels   = timeline->numFutRefLevel > 0;
        int offsetToFuturePath = timeline->numPastDates +
                                 timeline->offsetToClientPath;
        // then populate paths field in conjunction with cache
        for (int iAsset = 0; iAsset < wrapper.NbSimAssets(); iAsset++) {
            bool isCacheValid = pathCache->isValid(iAsset);
            int  numVolReq = wrapper.getNbPaths(iAsset);
            for (int iPath = 0; iPath < numVolReq; iPath++) {
                double& refLevel = wrapper.refLevel(iAsset, iPath);
                double* path = wrapper.path(iAsset, iPath);
                double* futPath = path+offsetToFuturePath;
                if (!isCacheValid || iPath > 0 ){
                    // need to generate path
                    wrapper.generatePath(pathIdx, iAsset, iPath);
                    // calculate refLevel
                    if (futureRefLevels){
                        refLevel = wrapper.getRefLevelPath()->refLevel(iAsset, path);
                    }
                    // then save it (only for first vol interp)
                    if (iPath == 0 && pathCache->updateAllowed(iAsset)){
                        pathCache->write(iAsset, pathIdx, futPath,
                                         refLevel);
                    }
                } else {
                    // can read path (and refLevel if in future) from cache
                    pathCache->read(iAsset, pathIdx, futPath, refLevel);
                }
            }
        }
    }

    virtual int getPathIndex() const {
        return nowPathIdx;
    }

    MCSimpleFuturePathGenerator(
        const MCPathGeneratorSP& pastPathGenerator,
        MCPathConfig*            pathConfig,
        const IMCProduct*         prod,
        MCPathBaseSP             basePathGen):
    // parent base class handles most things
    nowPathIdx(-1), pathCache(pathConfig->getPathCache()), wrapper(basePathGen) {

        static const string routine("MCSimpleFuturePathGenerator");
        try{
            // just need to set up our own variables
            beginIdx = pastPathGenerator->end(0); // start at the previous end
            endIdx = wrapper.getTimeline()->simSeries->getAllDates().size();
            // and to populate our path
            for (int iAsset=0; iAsset < wrapper.NbSimAssets(); iAsset++) {
                /* copy historic path information over (if there is past
                   then this implies numFutRefLevel = 0) */
                for (int i = 0; i < wrapper.getNbPaths(iAsset); i++){
                    const double* path = pastPathGenerator->Path(iAsset, i);
                    for (int j = 0; j < beginIdx; j++) {
                        wrapper.path(iAsset, i)[j] = path[j];
                    }
                }
            }
            // ensure path cache is configured properly
            int cacheLength = endIdx - beginIdx;
            bool futureRefLevels = wrapper.getTimeline()->numFutRefLevel > 0;
            for (int i = 0; i < wrapper.NbSimAssets(); i++){
                pathCache->validateSize(i, cacheLength, futureRefLevels);
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    IStateVariableSP create(const IStateVariableGen* svGen) {
        throw ModelException(
            "MCSimpleFuturePathGenerator::create",
            "Path generator does not support state variables");
    }
};


/** Class for generating future paths when each asset has a different set
    of simulation dates */
class MCPathBase::MCFuturePathGenerator: public Generator {
private:
    // fields ///
    vector<DoubleMatrix>      assetPaths; /* paths for each asset's dates */
    IntArray                  beginIdx;   /* index, per asset, into paths */
    IntArray                  endIdx;     /* index, per asset, into paths */
    IntArray                  assetNumFutRefLevel; /* number of future dates in
                                                   ref level per asset */
    IntArray                  assetPathIdx; // scratch pad for some calcs
    MCPathConfig::PathCacheSP pathCache; // cache of generated paths
    int                       nowPathIdx; // current path idx

    MCPathBaseWrapper         wrapper;    // Quick access to MCPathBase fields
public:
    /** Returns the path base with which the Generator was built */
    virtual MCPathBaseSP getPathBase(){
        return wrapper.getPathBase();
    }

    virtual int NbSimAssets() const {
        return wrapper.NbSimAssets();
    }

    virtual double refLevel(int iAsset, int iPath) const {
        return wrapper.refLevel(iAsset, iPath);
    }

    /** Returns the sum of the absolute value of the 'drifts'
        between each of the simulation dates being simulated (ie
        between begin() and end()) (where the first ratio is between
        the first simulation date and the simulation start date). Note
        past path generator will return 0 */
    virtual double maxDriftProduct(int iAsset) const{
        return wrapper.maxDriftProduct(iAsset);
    }

    virtual bool hasPast() const {
        return wrapper.getTimeline()->simHasPast;
    }

    virtual bool doingPast() const {
        return false;
    }

#ifdef MAPPATH_MACRO
    virtual const double* Path(int iAsset, int iPath) const{
        /* client expects array to be indexed such that the ref level dates
           are excluded */
        int mappedPath = MAPPATH(wrapper.isSinglePath, iPath);
        return (assetPaths[iAsset][mappedPath] +
                assetNumFutRefLevel[iAsset]);
    }
#else
    virtual const double* Path(int iAsset, int iPath) const{
        /* client expects array to be indexed such that the ref level dates
           are excluded */
        return (assetPaths[iAsset][wrapper.mapPath(iPath)] +
                assetNumFutRefLevel[iAsset]);
    }
#endif

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    virtual int begin(int iAsset) const{
        return beginIdx[iAsset];
    }

    //// see  begin(int iAsset)
    virtual int end(int iAsset) const{
        return endIdx[iAsset];
    }

    virtual void generatePath(int pathIdx) {
        const MCProductTimelineConstSP& timeline = wrapper.getTimeline();

        // cache the pathIdx for getPathIndex()
        nowPathIdx = pathIdx;
        // get hold of the correlated random numbers
        wrapper.drawRandomNumbers(pathIdx);
        // for ease, read some variables
        // start by calculating basePathGen->paths where cache is not valid
        int iAsset;
        for (iAsset = 0; iAsset < wrapper.NbSimAssets(); iAsset++) {
            bool isCacheValid = pathCache->isValid(iAsset);
            int numPaths = wrapper.getNbPaths(iAsset);
            for (int iPath = 0; iPath < numPaths; iPath++) {
                if (!isCacheValid || iPath > 0 ){
                    // need to generate path
                    wrapper.generatePath(pathIdx, iAsset, iPath);
                }
            }
        }

        // need to copy over relevant values for each asset for which the cache
        // is invalid
        assetPathIdx = beginIdx; // need to keep track of where we are
        int totalNumSteps = timeline->totalNumSteps;
        int numFutRefLevel = timeline->numFutRefLevel;
        for (int iStep = timeline->numPastDates;
             iStep < totalNumSteps; iStep++){
            // need to switch between doing RefLevel and SimSeries
            bool  doingRefLevel = iStep < timeline->numFutRefLevel;
            const IntArray& assetsOnDate = doingRefLevel?
                wrapper.getRefLevelPath()->assetsOnDate(iStep):
                timeline->simSeries->assetsOnDate(iStep - numFutRefLevel);
            for (int i = 0; i < assetsOnDate.size(); i++){
                int iAsset = assetsOnDate[i];
                int numPaths = wrapper.getNbPaths(iAsset);
                DoubleMatrix& assetPath = assetPaths[iAsset];
                bool isCacheValid = pathCache->isValid(iAsset);
                for (int iPath = 0; iPath < numPaths; iPath++){
                    double* allAssetsPath = wrapper.path(iAsset, iPath);
                    if (!isCacheValid || iPath > 0){
                        assetPath[iPath][assetPathIdx[iAsset]] =
                            allAssetsPath[iStep];
                    }
                }
                assetPathIdx[iAsset]++;
            }
        }
        // then read from or write to cache
        for (iAsset = 0; iAsset < wrapper.NbSimAssets(); iAsset++) {
            bool isCacheValid = pathCache->isValid(iAsset);
            DoubleMatrix& assetPath = assetPaths[iAsset];
            int numPaths = wrapper.getNbPaths(iAsset);
            for (int iPath = 0; iPath < numPaths; iPath++) {
                double* path = wrapper.path(iAsset, iPath);
                double& refLevel = wrapper.refLevel(iAsset, iPath);
                // for cache, want to jump over ref level dates and historic
                // sim dates
                double* futPath = assetPath[iPath] +
                    assetNumFutRefLevel[iAsset] + beginIdx[iAsset];
                if (!isCacheValid || iPath > 0 ){
                    // calculate refLevel
                    if (assetNumFutRefLevel[iAsset] > 0){ // performance boost
                        refLevel = wrapper.getRefLevelPath()->refLevel(iAsset, path);
                    }
                    // then save it (only for first vol interp)
                    if (iPath == 0 && pathCache->updateAllowed(iAsset)){
                        pathCache->write(iAsset, pathIdx,
                                         futPath, refLevel);
                    }
                } else {
                    // can read path (and refLevel if in future) from cache
                    pathCache->read(iAsset, pathIdx, futPath, refLevel);
                }
            }
        }
    }

    virtual int getPathIndex() const {
        return nowPathIdx;
    }

    MCFuturePathGenerator(
        const MCPathGeneratorSP& pastPathGenerator,
        MCPathConfig*            pathConfig,
        const IMCProduct*         prod,
        MCPathBaseSP             basePathGen):
        // parent base class handles most things
        assetPaths(basePathGen->NbSimAssets()),
        beginIdx(basePathGen->NbSimAssets()),
        endIdx(basePathGen->NbSimAssets()),
        assetNumFutRefLevel(basePathGen->NbSimAssets()),
        assetPathIdx(basePathGen->NbSimAssets()),
        pathCache(pathConfig->getPathCache()),
        nowPathIdx(-1), wrapper(basePathGen) {

        static const string routine("MCFuturePathGenerator");
        try{
            const MCProductTimelineConstSP& timeline = wrapper.getTimeline();

            for (int iAsset=0; iAsset < wrapper.NbSimAssets(); iAsset++) {
                assetNumFutRefLevel[iAsset] = timeline->refLevelObj->numFutureDates(
                        prod->getEffectiveSimStartDate(), iAsset);
                int nbInterps = wrapper.getNbPaths(iAsset);
                // start at the previous end
                beginIdx[iAsset] = pastPathGenerator->end(iAsset);
                endIdx[iAsset] = timeline->simSeries->numDates(iAsset);
                // total number includes reference level
                int numSteps = endIdx[iAsset] + assetNumFutRefLevel[iAsset];
                // allocate space
                assetPaths[iAsset] = DoubleMatrix(nbInterps, numSteps);
                /* copy historic path information over (if there is past
                   then this implies numFutRefLevel = 0) */
                for (int i = 0; i < wrapper.getNbPaths(iAsset); i++){
                    const double* path = pastPathGenerator->Path(iAsset, i);
                    for (int j = 0; j < beginIdx[iAsset]; j++){
                        assetPaths[iAsset][i][j] = path[j];
                    }
                }
            }
            // ensure path cache is configured properly
            for (int i = 0; i < wrapper.NbSimAssets(); i++){
                pathCache->validateSize(i, endIdx[i] - beginIdx[i],
                                        assetNumFutRefLevel[i] > 0);
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    IStateVariableSP create(const IStateVariableGen* svGen) {
        throw ModelException(
            "MCSimpleFuturePathGenerator::create",
            "Path generator does not support state variables");
    }
};

/** Creates a [future] path generator using supplied MCPathBase.
    Essentially a wrapper mapping requirements of MCPathGenerator onto
    properties of MCPathBase (eg handles different dates per
    asset etc) */
MCPathBase::GeneratorSP MCPathBase::createPathGenerator(
    const MCPathGeneratorSP& pastPathGenerator,
    MCPathConfig*            pathConfig,
    const IMCProduct*         prod,
    MCPathBaseSP             basePathGen){
    Generator* pathGen;
    /* switch between general case and specific case where each asset
       has the same simulation dates */
    if (prod->getRefLevel()->sameDatesPerAsset() &&
        prod->getSimSeries()->sameDatesPerAsset()){
        pathGen = new MCSimpleFuturePathGenerator(pastPathGenerator,
                                                  pathConfig,
                                                  prod,
                                                  basePathGen);
    } else {
        pathGen = new MCFuturePathGenerator(pastPathGenerator,
                                            pathConfig,
                                            prod,
                                            basePathGen);
    }
    return GeneratorSP(pathGen);
}


DRLIB_END_NAMESPACE
