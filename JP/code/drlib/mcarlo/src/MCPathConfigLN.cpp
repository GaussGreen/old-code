//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathGeneratorLN.cpp
//
//   Description : Monte Carlo path generator
//                 Since random number gen is here - so too must all
//                 the work that is expected to be done with the same
//                 set of randoms. So, "iPath", "iAsset", "iStep"
//                 and at some stage this might also mean "iTweak" for
//                 internal sensitivities
//
//
//   Date        : June 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Delta.hpp"
#include "edginc/CrossGamma.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductEngineClient.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/DependenceLocalCorr.hpp"

#define STATE_VARIABLES
#define STATEVAR_CACHING

#ifdef STATE_VARIABLES
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenBarrier.hpp"
#include "edginc/HitSample.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/MCPricing.hpp"  // Could remain in .cpp directory?
#include <algorithm>
#endif

DRLIB_BEGIN_NAMESPACE


/*****************************************
MCPathConfigLN
*****************************************/
/** Class splits its work into three smaller classes. This is possibly more
    code but it makes it much simpler and easier to develop */
class MCPathConfigLN: public MCPathConfig {
private:
    // fields ////
    string                        volType;
    bool                          isCarefulRandoms;
    string                        dependenceType;

    // transient
    class Generator;
    friend class Generator;
    class LRDelta;
    class LRCrossGamma;
    class MyLRGenerator;

    MCCacheManager cacheMgr;    //!< Cache manager $unregistered

public:
    class LRDeltaCommon;
    static CClassConstSP const TYPE;

    IObject* clone() const{
        MCPathConfigLN& myCopy = dynamic_cast<MCPathConfigLN&>(*(MCPathConfig::clone()));
        myCopy.cacheMgr = cacheMgr; // shallow copy
        return &myCopy;
    }

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const{
        MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));        
        dependenceMaker->modifyMarketDataFetcher(mdf);
        if(skewMaker.get()) {
            skewMaker->modifyMarketDataFetcher(mdf);
        }
        return mdf;        
    }


    /** See comment on pure virtual declaration. Defined below */
    virtual LRGenerator* createLRGenerator(
        const MCPathGeneratorSP& pathGen,
        int                      nbIter,
        int                      nbSubSamples);

#ifdef STATE_VARIABLES
    class Gen;              //!< Referee path generator class
    class PathGenSpot;      //!< Spot PathGen
    class PathGenHVBB;      //!< Brownian Bridge PathGen for hitting values
    class BarrAdj;          //!< Barrier adjustment class

    friend class Gen;
    friend class PathGenSpot;
    friend class PathGenHVBB;

    typedef refCountPtr<PathGenSpot> PathGenSpotSP;
    typedef refCountPtr<PathGenHVBB> PathGenHVBBSP;
    typedef refCountPtr<BarrAdj>     BarrAdjSP;
#endif


protected:
    MCPathConfigLN(CClassConstSP clazz,
                   const string& volType,
                   const string& dependenceType):
    MCPathConfig(clazz, DependenceMakerGaussTermSP(new DependenceMakerGaussTerm())),
    volType(volType), isCarefulRandoms(false), dependenceType("not used") {}

private:
    MCPathConfigLN(): MCPathConfig(TYPE, DependenceMakerGaussTermSP(new DependenceMakerGaussTerm())),
                      volType(IVolatilityBS::TYPE->getName()),
                      isCarefulRandoms(false),
                      dependenceType("not used") {}

    virtual bool vegaMatrixSupported() const { return true; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * Irrelevant since it's B-S not parametric.  See
     * IModel::wantsRiskMapping().
     */

    IModel::WantsRiskMapping wantsRiskMapping() const {
        return IModel::riskMappingIrrelevant;
    }

protected:
    /** Creates a past path generator */
    MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod);

    /** Creates a future path generator */
    MCPathGeneratorSP futurePathGenerator(
        int                      cachingMode,
        int                      numPaths,
        const MCPathGeneratorSP& pastPathGenerator,
        const IMCProduct*         prod,
        Control*                 control,
        Results*                 results,
        DateTimeArray&           simDates );

    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control,
        Results*                           results,
        DateTimeArray&                     simDates ); // defined below

    int storagePerPath(IMCProduct* product) const;

    bool sensShift(Theta* shift);

    virtual bool carefulRandoms() const { return isCarefulRandoms; }
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigLN, clazz);
        SUPERCLASS(MCPathConfig);
        EMPTY_SHELL_METHOD(defaultMCPathConfigLN);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);  // for now
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);
        FIELD(isCarefulRandoms, "isCarefulRandoms");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigLN(){
        return new MCPathConfigLN();
    }
};

CClassConstSP const MCPathConfigLN::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigLN", typeid(MCPathConfigLN), load);

/*****************************************
MCFuturePathGenBaseLN
*****************************************/
class MCPathConfigLN::Generator: virtual public MCPathBase,
                                 public DependenceMakerGauss::Support,
								 public MCGammaAdjustedPathGenerator,
                                 public DependenceMakerGaussTerm::Support,
                                 public DependenceMakerLocalCorr::Support {
    friend class MCPathConfigLN::LRDeltaCommon;
    friend class MCPathConfigLN;
    const static double MAX_DRIFT_RAND_NUMBER; /* what bounds to use when
                                                  calculating max drift */
    vector<CVolRequestLNArray> volRequests;    // [numAssets]
    vector<DoubleMatrix>       drifts;     /* x[numAssets] each with
                                              NbVolReq(iAsset) cols x
                                              NbSteps rows */
    vector<DoubleMatrix>       vols;
    DoubleArray                maxDrifts; // product of MAX(1, drifts)

    IMCRandomSP                randomGen;  // Random number generator
    DependenceSP               dependence; // Dependence object
    MCProductTimelineSP        timeline;   // Product timeline
    RefLevelDataSP             refData;    // RefLevel data

    IntArray              numPaths;     // Number of paths per asset
    vector<DoubleArray>   fwds;
    const IMultiFactors*  mAsset;       // Multi asset interface
    int                   numAssets;    // number of assets
    vector<DoubleMatrix>  productPaths; // [iAsset][iPath][iStep]

    // Additional fields used for LocalCorr dependence maker
    DependenceMakerSP           dependenceMaker;
    IRandomSP                   rand;
    MCPathConfig::RandomCacheSP idiosynFactorRandomCache;
    MCPathConfig::RandomCacheSP marketFactorRandomCache;
    bool                        isCarefulRandoms;

    // simulates random numbers
    virtual void drawRandomNumbers(int pathIdx) {
        randomGen->generate(pathIdx);
    }

public:
    Generator(MCPathConfigLN*          mcPathConfigLN,
              const MCPathGeneratorSP& pastPathGenerator,
              const IMCProduct*         prod):
        volRequests(prod->getNumAssets()),
        drifts(prod->getNumAssets()),
        vols(prod->getNumAssets()),
        maxDrifts(prod->getNumAssets(), 1.0),
        numPaths(getNumPaths(prod, pastPathGenerator, prod->getNumAssets())),
        fwds(prod->getNumAssets()),
        mAsset(prod->getMultiFactors()),
        numAssets(prod->getNumAssets()),
        productPaths(prod->getNumAssets()),
        // need to pass on the following three fields to LocalCorr
        dependenceMaker(mcPathConfigLN->getDependenceMaker()),
        rand(mcPathConfigLN->getRandomGenerator()), 
        idiosynFactorRandomCache(mcPathConfigLN->getIdiosynFactorRandomCache()),
        marketFactorRandomCache(mcPathConfigLN->getMarketFactorRandomCache()), 
        isCarefulRandoms(mcPathConfigLN->carefulRandoms()) {

        static const string routine("MCPathConfigLN::Generator");
        try{
            const IMCProductLN* prodLN =
                dynamic_cast<const IMCProductLN*>(prod);
            if (!prodLN){
                throw ModelException(routine, "Product does not implement "
                                     "IMCProductLN interface");
            }

            // Obtain product timeline
            timeline = getProductTimeline(prod, pastPathGenerator);

            // Obtain market data
            refData = getRefData(timeline, numPaths, mAsset, prod, pastPathGenerator);

            // Paths
            int iAsset;
            for (iAsset=0; iAsset < numAssets; iAsset++) {
                fwds[iAsset] = DoubleArray(timeline->futureDates.size(), 0.0);
                productPaths[iAsset] = DoubleMatrix(numPaths[iAsset], timeline->totalNumSteps);
                mAsset->factorFwdValues(iAsset, timeline->futureDates, fwds[iAsset]);
            }

            for (iAsset=0; iAsset < numAssets; iAsset++) {
                // for each asset have an array of vol requests
                volRequests[iAsset] =
                    prodLN->getVolInterp(pastPathGenerator.get(), iAsset);
            }
            prodLN->initialiseLN(pastPathGenerator.get());

            /* now we know how many paths per asset we can allocate
               space for them. Also note that we can't ask for
               the vol interps until the past has been computed */
            for (iAsset=0; iAsset < numAssets; iAsset++) {
                int nbInterps = volRequests[iAsset].size();
                drifts[iAsset]   = DoubleMatrix(nbInterps, timeline->numFutSteps);
                vols[iAsset]     = DoubleMatrix(nbInterps, timeline->numFutSteps);
            }

            // precompute drift and vols could possibly have a class
            // "process" which did this - for now leave as local functions
            generateDriftAndVols(timeline->futureDates);

            // set up Dependence            

            if (mcPathConfigLN->skewMaker.get() && (numAssets>1)) {
                DependenceMakerLocalCorr* dpm = 
                    dynamic_cast<DependenceMakerLocalCorr*>(mcPathConfigLN->skewMaker.get());
                if (!dpm) {
                    throw ModelException(routine, "Internal Error");
                }
                dependence = dpm->createDependence(this);
                // a special one IRandomSP
                randomGen = IMCRandomSP(new MCRandom(
                    0,                                  // path generator
                    dependence,                         // dependence
                    IRandomSP(),                        // random number generator
                    mcPathConfigLN->getRandomCache(),   // randomCache
                    true,                               // not used (isCarefulRandoms)
                    timeline->numFutSteps,              // numDates
                    numAssets,                          // numFactors
                    timeline->totalNumPastDates));      // numPastDates
            } else {
                dependence = mcPathConfigLN->dependenceMaker->createDependence(this);
                // the classic one         
                randomGen = IMCRandomSP(new MCRandom(
                    0,
                    dependence,
                    mcPathConfigLN->getRandomGenerator(),
                    mcPathConfigLN->getRandomCache(),
                    mcPathConfigLN->carefulRandoms(),
                    timeline->numFutSteps,
                    numAssets,
                    timeline->totalNumPastDates));
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigLN::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // implement methods since derivation from DependenceMakerGaussTerm::Support
    virtual DateTimeArray getSimDates() const {
        return timeline->futureDates;
    }
    virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const {
        // bool interpolateAtmFwd not used for LN
        // non SV approach: do not pass today to Driver
        int nbAssets = vols.size();
        int nbDates = timeline->futureDates.size() - 1;
        DoubleMatrix fwdVarAtDates(nbAssets, nbDates);
        for (int iAsset=0; iAsset<nbAssets; iAsset++) {
            for (int iDate=0; iDate<nbDates; iDate++) {
                fwdVarAtDates[iAsset][iDate] =
                    vols[iAsset][0][iDate]*vols[iAsset][0][iDate]; // take 1st column
            }
        }
        return fwdVarAtDates;
    }
    virtual const IMultiFactors* getMultiFactors() const {
        return mAsset;
    }

    // implement methods since derivation from DependenceMakerLocalCorr::Support
    virtual DependenceMakerSP getDependenceMaker() const {
        return dependenceMaker;
    }
    virtual int getNbPastDates() const {
        return timeline->totalNumPastDates;
    }
    virtual const IRandomSP getRandomGenerator() const {
        return rand;
    }
    virtual const MCPathConfig::RandomCacheSP getIdiosynFactorRandomCache() const {
        return idiosynFactorRandomCache;
    }
    virtual const MCPathConfig::RandomCacheSP getMarketFactorRandomCache() const {
        return marketFactorRandomCache;
    }
    virtual bool carefulRandoms() const {
        return isCarefulRandoms;
    }

protected:
    /** Generate path for specified path in simulation (pathIdx), for
        specified asset (iAsset) and specified vol interp (iPath) */
    virtual void generatePath(int pathIdx, int iAsset, int iPath){
        const DoubleMatrix& randoms = randomGen->getRandomNumbers();

        // populate paths field
        DoubleMatrix& assetPath = productPaths[iAsset]; // for ease/speed
        const DoubleMatrix& assetVols = vols[iAsset];
        const double* assetRandoms = randoms[iAsset];
        const DoubleMatrix& assetDrifts = drifts[iAsset];
        double* thisPath = assetPath[iPath];
        const double* pathVols = assetVols[iPath];
        const double* thisDrift = assetDrifts[iPath];
        for (int iStep = 0, modStep = timeline->numPastDates;
             iStep < timeline->numFutSteps; iStep++, modStep++) {

            double sigmaDtDz = pathVols[iStep] *assetRandoms[iStep];

            thisPath[modStep] = iStep == 0?
                refData->fwdsAtSimStart[iAsset]: thisPath[modStep-1];

            thisPath[modStep] *= thisDrift[iStep];

            thisPath[modStep] *= exp(sigmaDtDz);
        }
    }

	/** 3 functions implemented for gamma adjusted vol */
	virtual const DoubleMatrix& getVols(int iAsset) const {
		return vols[iAsset];
	}

	virtual const DoubleMatrix& getDrifts(int iAsset) const {
		return drifts[iAsset];
	}


	virtual const DoubleMatrix& getRandoms(int iAsset) const {
		return randomGen->getRandomNumbers();
	}

    /** Returns an array of integers containing the number of paths
        per asset. This information is inside the VolInterp */
    IntArray getNumPaths(const IMCProduct* prod,
                         const MCPathGeneratorSP& pastPathGenerator,
                         const int numAssets) const{

        static const string routine("MCPathConfigLN::Generator::getNumPaths");

        const IMCProductLN* prodLN = dynamic_cast<const IMCProductLN*>(prod);
        if (!prodLN){
            throw ModelException(routine, "Product does not implement "
                                 "IMCProductLN interface");
        }

        IntArray numPaths(numAssets);

        int iAsset; // MSVC broken
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            // for each asset have an array of vol requests
            numPaths[iAsset] =
                prodLN->getVolInterp(pastPathGenerator.get(), iAsset).size();
        }

        return numPaths;
    }

    void generateDriftAndVols(const DateTimeArray& futurePathDates) {
        static const string routine("MCFuturePathGenLN::"
                                    "generateDriftAndVols");
        try{
            CDoubleArray fwdVar(timeline->numFutSteps+1); // should be -1?
            for(int iAsset=0;iAsset<numAssets;iAsset++) {
                double theDrift = 1.0;
                for(int iPath=0; iPath < volRequests[iAsset].size(); iPath++)
                {
                    // interpolate the vol using our LN request
                    // XXX factor vs asset
                    CVolProcessedSP vol(
                        mAsset->factorGetProcessedVol(iAsset,
                                                      volRequests[iAsset][iPath].get()));
                    CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));
                    volBS->CalcVar(futurePathDates, volBS->forward, fwdVar);

                    // Referenced to the simulated dates (not the simStartDate)
                    // drifts[iStep] is drift FROM prev point TO iStep
                    // vols[iStep] is vol FROM prev point TO iStep
                    double totalVar = 0.0; // used for max drift calc
                    for (int iStep = 0; iStep < timeline->numFutSteps; iStep ++)
                    {
                        double drift = fwds[iAsset][iStep+1]/
                            fwds[iAsset][iStep];
                        drift *= exp(-0.5 * fwdVar[iStep]);
                        drifts[iAsset][iPath][iStep] = drift;
                        totalVar += fwdVar[iStep];
                        theDrift *= Maths::max(1.0, drift);
                        double vol = sqrt(fwdVar[iStep]);
                        vols[iAsset][iPath][iStep] = vol;
                    }
                    // here we're simulating what happens in generatePath
                    // want theDrift to represent [reasonable] worst
                    // case drift. A value of 2 for MAX_DRIFT_RAND_NUMBER
                    // means that we catch >96% of paths
                    theDrift *= exp(sqrt(totalVar)*MAX_DRIFT_RAND_NUMBER);
                    if (theDrift > maxDrifts[iAsset]){
                        maxDrifts[iAsset] = theDrift;
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }
    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0 */
    double maxDriftProduct(int iAsset) const{
        return maxDrifts[iAsset];
    }

    /** Obtains timeline object from base */
    MCProductTimelineConstSP getTimeline() const {
        return timeline;
    }

    /** Returns number of assets */
    int NbSimAssets() const {
        return numAssets;
    }

    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath) {
        return refData->refLevels[iAsset][iPath];
    }

    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const {
        return refData->refLevelPath;
    }

    /** Returns the paths */
    double* Path(int iAsset, int iPath) {
        return productPaths[iAsset][iPath];
    }

    /** Returns the number of paths per asset */
    int nbPaths(int iAsset) const {
        return numPaths[iAsset];
    }

    /** Returns if it is a single path generator or not */
    bool isSinglePath() const {
        return false;
    }
};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


#ifdef STATE_VARIABLES


/////////////////////////////////////////////////////////////////////////////

/** Referee class that distributes simulation to components
    e.g. spot, hitTime etc. */
class MCPathConfigLN::Gen: virtual public MCPathGenerator, // For backward compatibility
                           virtual public MCStatelessPathGen,
                           virtual public IStateVariableGen::IStateGen {
public:
    /** Constructor */
    Gen(MCPathConfigLN*          mcPathConfigLN,
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient,
        MCCacheManager&          cacheMgr,
        bool                     cachingRequested,
        DateTimeArray&           simDates);

    /** MCpathGenerator methods */
    // Deprecated methods
    virtual int NbSimAssets() const;

    virtual const double* Path(int iAsset, int iPath) const;

    virtual double refLevel(int iAsset, int iPath) const;

    virtual double maxDriftProduct(int iAsset) const;

    virtual int begin(int iAsset) const;

    virtual int end(int iAsset) const;

    // Live methods
    virtual bool hasPast() const;

    virtual bool doingPast() const;

    virtual void generatePath(int pathIdx);

    virtual void advance();

    virtual void reset();

    virtual int getPathIndex() const;

    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen);

private:
    StateVarDBase                    svDBase;        //!< Collection of Generators + statevars
    MCPathConfigLN::PathGenSpotSP    pathGenSpot;    //!< Spot generator
    MCPathConfigLN::PathGenHVBBSP    pathGenHVBB;    //!< Hit value Brownian Bridge generator
    MCPathConfigLN::BarrAdjSP        barrAdj;        //!< Barrier adjustment
    int                              nowPathIdx;     //!< Current path idx
};


/////////////////////////////////////////////////////////////////////////////

/** Spot path generator for LogNormal process using state variable approach */
class MCPathConfigLN::PathGenSpot: virtual public MCPathGen,
                                   public DependenceMakerGauss::Support,
                                   public DependenceMakerGaussTerm::Support {
private:
    // friend class MCPathConfigLN::LRDeltaCommon;
    // friend class MCPathConfigLN;
    // friend class MCPathConfigLN::Gen;

    int                        numAssets;       //!< Number of assets
    DateTimeArray              simTimeline;     //!< Today + strictly future merged dates
    bool                       simHasPast;      //!< Whether product has past

    vector<CVolRequestLNSP>    volRequests;     //!< Volatility request [iAsset]
    DoubleArray                spotAtStart;     //!< Spot at simulation start
    DoubleMatrix               drifts;          //!< Log spot drift  at simTimeline     [iAsset][iStep]
    DoubleMatrix               vols;            //!< Sqrt of fwd var at simTimeline     [iAsset][iStep]
    DoubleMatrix               fwds;            //!< Asset Forwards  at simTimeline     [iAsset][iStep]
    DoubleMatrix               driverPaths;     //!< Driver paths    at simTimeline     [iAsset][iStep]
    vector<BoolArray>          doSpot;          //!< Simulate spot   at simTimeline     [iAsset][iStep]
    vector<DoubleArray>        spotProdPaths;   //!< Spot price at asset specific dates [iAsset][iStep]
    vector<int>                spotProdOffset;  //!< Starting point for future path per asset
    vector<double>             maxDrifts;       //!< product of MAX(1, drifts)

    MCRandomGenCacheSP         randomGen;       //!< Random number generator for normals
    DependenceSP               dependence;      //!< Dependence object e.g. gaussian

    const IMultiFactors*       mAsset;          //!< Multi asset interface
    SVGenSpotPathCacheSP          spotPathCache;   //!< Cache of generated paths

    static const MCCache::KeyUtils::Key keyPathCache;   //!< Key for spot path cache
    static const MCCache::KeyUtils::Key keyRandCache;   //!< Key for spot path random numbers cache
    static const double MAX_DRIFT_RAND_NUMBER;          //!< Bound to use for max drift

    vector<IAdvanceableStateVariable*> advanceSvSet; // for stateless advance
public:
    /** Constructor */
    PathGenSpot(MCPathConfigLN*          mcPathConfigLN,
                const PastPathGenSpotSP& pastPathGenSpot,
                const MCProductClient*   prodClient,
                const SVGenSpotArray&       spotGenArray,
                const DateTimeArray&     driverDates,
                MCCacheManager&          cacheMgr,
                bool                     cachingRequested,
                StateVarDBase&           svDBase,
                DateTimeArray&           simDates):
    numAssets(prodClient->getNumAssets()),
    simHasPast(pastPathGenSpot->hasPast()),
    volRequests(numAssets),
    spotAtStart(numAssets),
    doSpot(numAssets),
    spotProdPaths(numAssets),
    spotProdOffset(numAssets, 0),
    maxDrifts(numAssets, 1.0),
    mAsset(prodClient->getMultiFactors()) {

        static const string routine("MCPathConfigLN::PathGenSpot::PathGenSpot");
        try{
            const DateTime& today = prodClient->getToday();
            const IPastValues* pastValues = prodClient->getMCPastValues();

            // Get driver and spot dates
            for(int iDate = 0; iDate < driverDates.size(); iDate++) {
                if(today > driverDates[iDate]) {
                    throw ModelException("Driver dates must be strictly in the future.");
                }
            }
            DateTimeArraySP spotGenDates = MCPath::getAllDates(spotGenArray);

            // Enrich driver dates by spot dates and drop past dates
            DateTimeArray allDriverDates = DateTime::merge(driverDates, *spotGenDates);
            DateTimeArray futDriverDates = today.getFutureDates(allDriverDates);

            // Create simulation timeline
            simTimeline = futDriverDates;
            simTimeline.insert(simTimeline.begin(), today);
            int nbRandoms = simTimeline.size() - 1;

            if (dynamic_cast<const IMCStatelessProductClient*>(prodClient))
                simDates = simTimeline;

            // Allocate memory for matrices now that timeline is known
            drifts      = DoubleMatrix(numAssets, simTimeline.size());
            vols        = DoubleMatrix(numAssets, simTimeline.size());
            fwds        = DoubleMatrix(numAssets, simTimeline.size());
            driverPaths = DoubleMatrix(numAssets, simTimeline.size());

            // At this stage, future spot dates are a subset of driver dates.
            // So we will always simulate the driver and sometimes the spot.

            // Create spot paths and driver paths per asset
            DateTimeArrayArray spotDatesPerAsset(numAssets);
            vector<const double*> spotPtrs(numAssets);
            vector<int> spotBeginInd(numAssets);
            vector<int> spotEndInd(numAssets);

            int iAsset;
            for (iAsset=0; iAsset < numAssets; iAsset++) {
                // 1) SPOT PATHS
                // Populate past values
                DateTimeArraySP spotAssetDates =
                    MCPath::getAllDates(spotGenArray, iAsset);
                spotProdPaths[iAsset]   = DoubleArray(spotAssetDates->size());
                DateTimeArray assetPastDates = today.getPastDates(*spotAssetDates);
                DoubleArray assetPastValues =
                    pastValues->getPastValues(assetPastDates, iAsset, today);
                spotProdOffset[iAsset] = assetPastValues.size();
                int iStep;
                for(iStep = 0; iStep < assetPastValues.size(); iStep++) {
                    spotProdPaths[iAsset][iStep] = assetPastValues[iStep];
                }

                // Create spot mappings
                spotDatesPerAsset[iAsset] = *spotAssetDates;
                spotPtrs[iAsset] = &spotProdPaths[iAsset][0];
                spotBeginInd[iAsset] = assetPastDates.size();
                spotEndInd[iAsset]   = spotAssetDates->size();

                // 2) SPOT AND DRIVER FLAGS
                BoolArray& doSpotPerAsset = doSpot[iAsset];
                doSpotPerAsset = BoolArray(simTimeline.size(), false);
                // Be careful not to use iStep = 0 because today might be
                // an asset date but it has been dealt with in the past
                for(iStep = 1; iStep < simTimeline.size(); iStep++) {
                    int occurances = count(spotAssetDates->begin(),
                                           spotAssetDates->end(),
                                           simTimeline[iStep]);
                    doSpotPerAsset[iStep] = occurances > 0 ? true: false;
                }
            }

#ifdef STATEVAR_CACHING
            // Caching of spot paths
            IntArray numDatesPerAsset(numAssets, 0);
            for (iAsset=0; iAsset < numAssets; iAsset++) {
                numDatesPerAsset[iAsset] = spotEndInd[iAsset] - spotBeginInd[iAsset];
            }
            spotPathCache = SVGenSpotPathCache::createCache(
                cacheMgr, keyPathCache, numDatesPerAsset, cachingRequested);

            // Caching of ramdom numbers
            MCRandomCacheSP randCache = MCRandomCache::createCache(cacheMgr,
                keyRandCache, numAssets, nbRandoms, cachingRequested);
#endif

            // Create SVGenSpot::IStateVars and put them in database
            MCPath::IStateVarArray spotSVArray(MCPath::createPaths(
                false,
                spotGenArray,
                spotDatesPerAsset,
                spotBeginInd,
                spotEndInd,
                spotPtrs,
                maxDrifts));

            unsigned int iVar;
            for(iVar = 0; iVar < spotGenArray.size(); iVar++) {
                svDBase.append(spotGenArray[iVar], spotSVArray[iVar]);
            }

            // populate advanceables
            svDBase.filterAdvanceables(advanceSvSet);

            // Create Forwards at simulation timeline
            // Paths
            for (iAsset=0; iAsset < numAssets; iAsset++) {
                double* assetFwds = fwds[iAsset];
                for(int iStep = 0; iStep < simTimeline.size(); iStep++) {
                    assetFwds[iStep] = mAsset->factorFwdValue(iAsset, simTimeline[iStep]);
                }
                spotAtStart[iAsset]   = fwds[iAsset][0];
            }

            const IMCProductLN* prodLN =
                dynamic_cast<const IMCProductLN*>(prodClient);
            if (!prodLN){
                throw ModelException(routine, "Product does not implement the "
                                     "IMCProductLN interface");
            }

            for (iAsset=0; iAsset < numAssets; iAsset++) {
                // for each asset have an array of vol requests
                CVolRequestLNArray assetRequests =
                    prodLN->getVolInterp(0, iAsset);    // pass null past path gen
                if(assetRequests.size() > 1) {
                    throw ModelException("Multiple paths per asset not supported.");
                }
                volRequests[iAsset] = assetRequests.front();
            }

            // precompute drift and vols could possibly have a class
            // "process" which did this - for now leave as local functions
            generateDriftAndVols(simTimeline);

            // set up Dependence at driver dates
            dependence = mcPathConfigLN->dependenceMaker->createDependence(this);

            // Initialize random number generator
            randomGen = MCRandomGenCacheSP(new MCRandomGenCache(
                0,
                dependence,
                mcPathConfigLN->getRandomGenerator(),
                randCache,
                mcPathConfigLN->carefulRandoms(),
                nbRandoms,
                numAssets,
                allDriverDates.size() - futDriverDates.size()));

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** Simulates paths for all assets. Part of the MCPathGen IFace. */
    virtual void generatePath(int pathIdx) {
        // Draw random numbers
        randomGen->generate(pathIdx);

        for (int iAsset = 0; iAsset < numAssets; iAsset++) {
#ifdef STATEVAR_CACHING
            int offset = spotProdOffset[iAsset];
            double* futPath = &spotProdPaths[iAsset][offset];
            if (spotPathCache->isValid(iAsset)){
                // Read paths from cache
                spotPathCache->read(iAsset, pathIdx, futPath);
            } else {
                // Diffuse paths and write them to cache
                generatePath(pathIdx, iAsset);
                if (spotPathCache->updateAllowed(iAsset)) {
                    spotPathCache->write(iAsset, pathIdx, futPath);
                }
            }
#else
            generatePath(pathIdx, iAsset);
#endif
        }

        reset();
    }

    /** Part of the MCPathGen IFace. */
    virtual bool doingPast() const {
        return false;
    }

    bool hasPast() const {
        return simHasPast;
    }

    double maxDriftProduct(int iAsset) const {
        return maxDrifts[iAsset];
    }

    // for stateless payoff
    virtual void reset() {
        for (size_t i = 0; i < advanceSvSet.size(); ++i)
            advanceSvSet[i]->reset();
    }

    virtual void advance() {
        for (size_t i = 0; i < advanceSvSet.size(); ++i)
            advanceSvSet[i]->advance();
    }

    /** Information passed to driver users i.e. driver
        values and vols for all simulation dates */
    class Driver {
    public:
        Driver(const DateTimeArray& simTimeline,
               const DoubleMatrix&  driverPath,
               const DoubleMatrix&  driverVols):
        simTimeline(simTimeline), driverPath(driverPath), driverVols(driverVols) {}

        ~Driver() {}

        inline static double spotFromDriver(double driver) {
            return MCPathConfigLN::PathGenSpot::spotFromDriver(driver);
        }

        inline static double driverFromSpot(double spot) {
            return MCPathConfigLN::PathGenSpot::driverFromSpot(spot);
        }

        const DateTimeArray& simTimeline;     //!< Timeline of driver
        const DoubleMatrix&  driverPath;      //!< Driver paths at simTimeline
        const DoubleMatrix&  driverVols;      //!< Driver forward vols at simTimeline
    };
    DECLARE_REF_COUNT(Driver);

    /** Returns a driver object */
    DriverConstSP getDriver() const {
        return DriverConstSP(new Driver(simTimeline, driverPaths, vols));
    }

    /** Mapping function from spot to driver */
    inline static double spotFromDriver(double driver) {
        return exp(driver);
    }

    /** Mapping function from driver to spot */
    inline static double driverFromSpot(double spot) {
        return log(spot);
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigLN::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // implement methods since derivation from DependenceMakerGaussTerm::Support
    virtual DateTimeArray getSimDates() const {
        return simTimeline;
    }
    virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const {
        // bool interpolateAtmFwd not used for LN
        // SV approach is different from nonSV approach: today needs to be passed to Driver
        int nbDates = simTimeline.size()-1;
        int nbAssets = vols.numCols();
        DoubleMatrix fwdVarAtDates(nbAssets,nbDates);
        for (int iDate=0; iDate<nbDates; iDate++) {
            for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                fwdVarAtDates[iAsset][iDate] =
                    vols[iAsset][iDate+1]*vols[iAsset][iDate+1];
            }
        }
        return fwdVarAtDates;
    }
    virtual const IMultiFactors* getMultiFactors() const {
        return mAsset;
    }

protected:
    /** Generate path for specified path in simulation (pathIdx), for
        specified asset (iAsset) */
    void generatePath(int pathIdx, int iAsset) {
        const DoubleMatrix& randoms = randomGen->getRandomNumbers();

        // populate driver and spot paths
        double* spotPath         = &spotProdPaths[iAsset][0];
        double* driverPath       = driverPaths[iAsset];
        const double* vol        = vols[iAsset];
        const double* drift      = drifts[iAsset];
        const double* stdNormals = randoms[iAsset];
        // WARNING: DO NOT USE const bool* myBoolPtr = &doSpot[iAsset][0]
        // as the implementation of array<bool, bool> is using T = unsigned int
        // and the pointer arithmetic does not work
        const BoolArray& doSpotPerAsset = doSpot[iAsset];

        driverPath[0] = drift[0];
        for (int iStep = 0, modStep = spotProdOffset[iAsset];
             iStep < simTimeline.size() - 1; iStep++) {
            // Diffuse driver according to normal SDE
            int iNextStep = iStep + 1;
            driverPath[iNextStep] = driverPath[iStep] + drift[iNextStep] +
                                    vol[iNextStep] * stdNormals[iStep];

            // Evaluate spot if necessary
            if(doSpotPerAsset[iNextStep]) {
                // Cheat here to make computation faster:
                // driverPath[iStep] is generally big as driverPath[0] = log(spot)
                // and exp(driverPath[iStep]) is slow.
                // Instead, take exponent of the increment and multiply by spot
                spotPath[modStep] = spotAtStart[iAsset] *
                    spotFromDriver(driverPath[iNextStep] - driverPath[0]);
                modStep++;
            }
        }
    }

    /** Precomputes drifts and vols for LN path generator */
    void generateDriftAndVols(const DateTimeArray& simTimeline) {
        static const string routine("MCPathConfigLN::PathGenSpot::generateDriftAndVols");
        try{
            // Get forward variance between dates hence N-1
            CDoubleArray fwdVar(simTimeline.size() - 1);
            for(int iAsset=0;iAsset<numAssets;iAsset++) {
                CVolProcessedSP vol(
                    mAsset->factorGetProcessedVol(iAsset,
                                                  volRequests[iAsset].get()));
                CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));
                volBS->CalcVar(simTimeline, volBS->forward, fwdVar);

                // Referenced to the simulated dates (not the simStartDate)
                // drifts[iStep] is drift FROM prev point TO iStep
                // vols[iStep] is vol FROM prev point TO iStep
                double theDrift = 1.0;
                double totalVar = 0.0; // used for max drift calc

                drifts[iAsset][0] = driverFromSpot(spotAtStart[iAsset]);
                vols[iAsset][0]   = 0.0;
                for (int iStep = 0; iStep < simTimeline.size() - 1; iStep ++) {
                    int iNextStep = iStep + 1;

                    // Lognormal drift and volatility
                    double var   = fwdVar[iStep];
                    double drift = log(fwds[iAsset][iNextStep]/fwds[iAsset][iStep]) -
                                   0.5 * var;

                    drifts[iAsset][iNextStep] = drift;
                    vols[iAsset][iNextStep]   = sqrt(var);

                    // Used for quick greeks
                    totalVar += var;
                    theDrift *= Maths::max(1.0, exp(drift));
                }
                // here we're simulating what happens in generatePath
                // want theDrift to represent [reasonable] worst
                // case drift. A value of 2 for MAX_DRIFT_RAND_NUMBER
                // means that we catch >96% of paths
                theDrift *= exp(sqrt(totalVar)*MAX_DRIFT_RAND_NUMBER);
                if (theDrift > maxDrifts[iAsset]){
                    maxDrifts[iAsset] = theDrift;
                }
            }

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }
};

// Initialize cache keys
const MCCache::KeyUtils::Key MCPathConfigLN::PathGenSpot::keyPathCache = MCCache::KeyUtils::getNewKey();
const MCCache::KeyUtils::Key MCPathConfigLN::PathGenSpot::keyRandCache = MCCache::KeyUtils::getNewKey();

 /** Bound to use when calculating max drift */
const double MCPathConfigLN::PathGenSpot::MAX_DRIFT_RAND_NUMBER = 3.0;


/////////////////////////////////////////////////////////////////////////////


/** Hitting values path generator for LogNormal process
    using state variable approach */
class MCPathConfigLN::PathGenHVBB: virtual public MCPathGen,
                                   public DependenceMakerGauss::Support,
                                   public DependenceMakerGaussTerm::Support {
public:
    /** Constructor */
    PathGenHVBB(MCPathConfigLN*              mcPathConfigLN,
                const MCProductClient*       prodClient,
                IStateVariableGen::IStateGen*    pastPathGen,
                IStateVariableGen::IStateGen*    futPathGen,
                const PathGenSpotSP&         pathGenSpot,
                const SVGenBarrierHVBBArray& barrierHVBBGens,
                MCCacheManager&              cacheMgr,
                bool                         cachingRequested,
                StateVarDBase&               svDBase):
    numBarriers(barrierHVBBGens.size()),
    numAssets(prodClient->getNumAssets()),
    mAsset(prodClient->getMultiFactors()),
    hitNoHitPath(numAssets),
    doHitNoHit(numAssets) {

        static const string routine = "MCPathConfigLN::PathGenHVBB::PathGenHVBB";

        try {
            if(!barrierHVBBGens.size()) {
                throw ModelException("No barriers to simulate. Internal error.");
            }

            // Construct a LinearInterpolant for the standard normal distribution
            int stdDevs = 6;
            int numNormals = 1000;
            LinearInterpolatorSP interpolator = LinearInterpolatorSP(new LinearInterpolator());
            tabulatedNorm = EquidistantLinearInterpolantNonVirtualSP::constCast(
                LinearInterpolator::computeEquidInterpNV(
                *interpolator, -stdDevs, stdDevs, numNormals, N1));

            // Validate the object and assume for the moment that barriers are disjoint
            validate(barrierHVBBGens);

            // Obtain the driver from the path gen
            driver = pathGenSpot->getDriver();

            // Join all barriers together in a common timeline for BBs.
            joinBarriers(pastPathGen, futPathGen, barrierHVBBGens, driver, svDBase);

#ifdef STATEVAR_CACHING
            // Caching of HitNoHit values
            IntArray numDatesPerAsset(numAssets, 0);
            for (int iAsset=0; iAsset < numAssets; iAsset++) {
                numDatesPerAsset[iAsset] = hitNoHitPath[iAsset].size();
            }
            hitValuePathCache = MCHitValuePathCache::createCache(
                cacheMgr, keyHitValueCache, numDatesPerAsset, cachingRequested);

            // Caching of ramdom numbers
            MCRandomCacheSP randCache = MCRandomCache::createCache(cacheMgr,
                keyRandCache, numAssets, numBBs, cachingRequested);
#endif

            // set up Dependence at driver dates
            dependence = mcPathConfigLN->dependenceMaker->createDependence(this);

            // Get random number generator
            IRandomSP rand = mcPathConfigLN->getRandomGenerator();
            randomGen = MCRandomGenCacheSP(new MCRandomGenCache(
                0,
                dependence,
                rand,
                randCache,
                mcPathConfigLN->carefulRandoms(),
                numBBs,
                numAssets,
                0));

        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Simulates brownian bridges for all assets. Part of the MCPathGen IFace. */
    virtual void generatePath(int pathIdx) {
        // Draw random numbers for all assets
        randomGen->generate(pathIdx);

        // Generate path per asset. Need to read / write from cache here
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
#ifdef STATEVAR_CACHING
            int* futPath = &hitNoHitPath[iAsset][0];
            if (hitValuePathCache->isValid(iAsset)){
                // Read paths from cache
                hitValuePathCache->read(iAsset, pathIdx, futPath);
            } else {
                // Diffuse paths and write them to cache
                generatePath(pathIdx, iAsset);
                if (hitValuePathCache->updateAllowed(iAsset)) {
                    hitValuePathCache->write(iAsset, pathIdx, futPath);
                }
            }
#else
            generatePath(pathIdx, iAsset);
#endif
        }
    }

    /** Part of the MCPathGen IFace. */
    virtual bool doingPast() const {
        return false;
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigLN::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // implement getGaussTermData since derivation from DependenceMakerGaussTerm::Support
    virtual DateTimeArray getSimDates() const {
        return driver->simTimeline;
    }
    virtual DoubleMatrix getFwdVarAtDates(bool interpolateAtmFwd) const {
        // bool interpolateAtmFwd not used for LN
        int nbDates = driver->simTimeline.size()-1;
        int nbAssets = driver->driverVols.numCols();
        DoubleMatrix fwdVarAtDates(nbAssets,nbDates);
        for (int iDate=0; iDate<nbDates; iDate++) {
            for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                fwdVarAtDates[iAsset][iDate] =
                    driver->driverVols[iAsset][iDate+1]*driver->driverVols[iAsset][iDate+1];
            }
        }
        return fwdVarAtDates;
    }
    virtual const IMultiFactors* getMultiFactors() const {
        return mAsset;
    }
private:
    /** Simulates brownian bridges for specified asset. */
    void generatePath(int pathIdx, int iAsset) {
        // Obtain random numbers
        const DoubleMatrix& randoms = randomGen->getRandomNumbers();

        const double* assetRandoms = randoms[iAsset];
        const double* assetDriver  = driver->driverPath[iAsset];
        int* assetHitNoHit         = &hitNoHitPath[iAsset][0];
        BoolArray& doHitNoHitAsset = doHitNoHit[iAsset];

        // Loop on the path and do BBs
        for(int iStep = 0, modStep = 0; iStep < maxNumBBs; iStep++) {
            if(doHitNoHitAsset[iStep]) {
                const BridgeData::AssetBridgeDataSP& assetBBData = bridgeData[iStep]->assetData[iAsset];
                double ref = assetBBData->refLevelSV->refLevel(iAsset);
                double driverRef = driver->driverFromSpot(ref);

                // Update BB, obtain uniform and get a sample
                HitNoHitBB& bb = *assetBBData->bbSample;
                double driverStart = assetDriver[iStep] - driverRef;
                double driverEnd   = assetDriver[iStep + 1] - driverRef;
                double uniform = tabulatedNorm->value(assetRandoms[modStep]);
                assetHitNoHit[iStep] = bb.updateAndSample(driverStart,driverEnd,uniform);
                modStep++;
            }
        }
    }

    /** Basic validation of the object */
    void validate(const SVGenBarrierHVBBArray& barrierHVBBGens) {
        static const string method = "MCPathConfigLN::PathGenHVBB::validate";

        try {
            // Get the barrier data
            vector<BarrierDataConstSP> barriers;
            unsigned int iBarrier;
            for(iBarrier = 0; iBarrier < barrierHVBBGens.size(); iBarrier++) {
                barriers.push_back(barrierHVBBGens[iBarrier]->getBarrierData());
            }

            // First implementation: only allow disjoint regions for each asset
            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                DateTimeArray startDates, endDates;
                for(iBarrier = 0; iBarrier < barriers.size(); iBarrier++) {
                    const DateTimeArray& monDates = *barriers[iBarrier]->getAssetBarrier(iAsset)->monitoringDates;
                    // Start and end dates for this barrier for the corresponding asset
                    startDates.push_back(monDates.front());
                    endDates.push_back(monDates.back());
                }

                // Validate that startDates and EndDates of other barriers are mutually disjoint
                // for this asset
                for(iBarrier = 0; iBarrier < barriers.size(); iBarrier++) {
                    for(unsigned int jBarrier = iBarrier + 1; jBarrier < barriers.size(); jBarrier++) {
                        bool disjoint =
                            (endDates[iBarrier] < startDates[jBarrier]) ||
                            (endDates[jBarrier] < startDates[iBarrier]);
                        if(! disjoint) {
                            throw ModelException("Barrier objects must be disjoint.");
                        }
                    }
                }
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Create common timeline from N barriers */
    void joinBarriers(IStateVariableGen::IStateGen*    pastPathGen,
                      IStateVariableGen::IStateGen*    futPathGen,
                      const SVGenBarrierHVBBArray& barrierHVBBGens,
                      PathGenSpot::DriverConstSP   driver,
                      StateVarDBase&               svDBase) {

        static const string method = "MCPathConfigLN::PathGenHVBB::joinBarriers";

        try {
            int iAsset, iStep, iBarrier;

            // Offsets for reading hitNoHitPath in each SV
            vector<vector<IntArray> >  offsets(numBarriers);
            for(iBarrier = 0; iBarrier < numBarriers; iBarrier++) {
                offsets[iBarrier] = vector<IntArray>(numAssets);
                for(iAsset = 0; iAsset < numAssets; iAsset++) {
                    offsets[iBarrier][iAsset] = IntArray(0);
                }
            }

            // Flag dates at which we need samples
            maxNumBBs = driver->simTimeline.size() - 1;
            if(maxNumBBs < 1) {
                throw ModelException("Need at least 2 dates.");
            }

            // Decide on number of divs to be used in the simulation
            int maxNumDivsPerYear = barrierHVBBGens[0]->getBarrierData()->maxNumDivsPerYear;
            const DateTime& valueDate = barrierHVBBGens[0]->getBarrierData()->valueDate;
            for(iBarrier = 1; iBarrier < numBarriers; iBarrier++) {
                if(maxNumDivsPerYear != barrierHVBBGens[iBarrier]->getBarrierData()->maxNumDivsPerYear) {
                    throw ModelException(
                        "Maximum number of divs per years has to be common between barriers."
                        "Check barrier objects number " +
                        Format::toString(iBarrier + 1) + " and " +
                        Format::toString(iBarrier));
                }
            }

            const DateTimeArray& simTimeline = driver->simTimeline;
            const DateTime& simStart = simTimeline.front();
            const DateTime& simEnd   = simTimeline.back();
            int totalNbDivs = maxNumDivsPerYear == -1 ? maxNumDivsPerYear :
                (int)(maxNumDivsPerYear * simStart.yearFrac(simEnd));
            if(maxNumDivsPerYear != -1 && maxNumDivsPerYear != 0) {
                // Use a minimum number of divs
                totalNbDivs = Maths::max(totalNbDivs, BarrierData::MIN_DIVS_PER_PERIOD);
            }

            // Brownian Bridge data
            bridgeData = vector<BridgeDataSP>(maxNumBBs);
            vector<TimeMetricConstSP> assetTimeMetrics(numAssets);
            vector<DividendListSP> driverDivs(numAssets);
            for(iAsset = 0; iAsset < numAssets; iAsset++) {
                doHitNoHit[iAsset]   = BoolArray(maxNumBBs);
                hitNoHitPath[iAsset] = IntArray(maxNumBBs);

                // Get asset time metric
                ATMVolRequestSP  volRequest(new ATMVolRequest());
                CVolProcessedSP vol(mAsset->factorGetProcessedVol(iAsset, volRequest.get()));
                assetTimeMetrics[iAsset] = vol->GetTimeMetric();

                // Get asset dividends
                driverDivs[iAsset] = AssetUtil::getDiscreteDivs(
                    &mAsset->getAsset(iAsset), valueDate, simStart, simEnd,
                    totalNbDivs, DividendCollector::DOLLAR_TO_YIELD);
            }

            // Flag if date is useful or not
            for(iStep = 0; iStep < maxNumBBs; iStep++) {
                int iNextStep = iStep + 1;
                const DateTime& lastDate = simTimeline[iStep];
                const DateTime& thisDate = simTimeline[iStep + 1];
                BridgeDataSP bbData(new BridgeData(lastDate, thisDate, numAssets));

                for(iAsset = 0; iAsset < numAssets; iAsset++) {
                    for(iBarrier = 0; iBarrier < numBarriers; iBarrier++) {
                        BarrierDataConstSP data = barrierHVBBGens[iBarrier]->getBarrierData();
                        BarrierPerAssetDataConstSP assetBarrier = data->getAssetBarrier(iAsset);
                        DateTimeArraySP monitoringDates = assetBarrier->monitoringDates;

                        // Define the start and end date for this asset's barrier
                        DateTime startDate = monitoringDates->front();
                        const DateTime& endDate   = monitoringDates->back();
                        DateTimeArray firstFutureDailyMonDate = assetBarrier->getFirstFutureDailyMonDate(
                            mAsset->getAsset(iAsset) , valueDate, data->monitorType);
                        if(firstFutureDailyMonDate.size()) {
                            // This will make sure we will not do BBs before firstFutureDailyMonDate
                            startDate = firstFutureDailyMonDate[0];
                        }

                        if(startDate < thisDate && thisDate <= endDate) {
                            doHitNoHit[iAsset][iStep] = true; // sample the asset
                            // Get the reference level from the appropriate barrier object
                            IRefLevel::IStateVarSP oldStateVar = barrierHVBBGens[iBarrier]->
                                getRefLevelGen()->getRefLevelSV(IRefLevel::IStateVarSP(   ),
                                              pastPathGen);
                            IRefLevel::IStateVarSP refLevelSV = barrierHVBBGens[iBarrier]->getRefLevelGen()->
                                getRefLevelSV(oldStateVar, futPathGen);

                            // Create driver's barrier at beginning and end
                            ScheduleSP schedule = assetBarrier->levels;
                            double barrierStart = driver->driverFromSpot(schedule->interpolate(lastDate));
                            double barrierEnd;
                            if(CString::equalsIgnoreCase(schedule->getInterp(), Schedule::INTERP_STAIRS)) {
                                // Step interpolation: use end barrier
                                barrierEnd = barrierStart;
                            } else {
                                // Other i.e. None of Linear: use interpolation at start of barrier
                                barrierEnd = driver->driverFromSpot(schedule->interpolate(thisDate));
                            }
                            double driverVol = driver->driverVols[iAsset][iNextStep];
                            double driverVar = Maths::square(driverVol);

                            // See if we need to do a barrier shift for daily monitoring
                            bool adjust = CString::equalsIgnoreCase(
                                data->monitorType, BarrierData::DAILY_MONITORING);
                            if(adjust) {
                                // Get business days between dates using holidays
                                HolidayConstSP hols = AssetUtil::getHoliday(&mAsset->getAsset(iAsset));
                                int numBusDays = hols->businessDaysDiff(lastDate, thisDate);
                                double barrierShift;
                                if(numBusDays) {
                                    int isUp = assetBarrier->isUp ? 1 : -1;
                                    barrierShift = isUp * ADJUST_CONSTANT * driverVol / sqrt((double)(numBusDays));
                                } else {
                                    barrierShift = 0.0;
                                }

                                barrierStart += barrierShift;
                                barrierEnd   += barrierShift;
                            }

                            double bridgeLength = assetTimeMetrics[iAsset]->yearFrac(lastDate, thisDate);

                            // Obtain dividends for this brownian bridge
                            DividendListSP assetBridgeDivs(driverDivs[iAsset]->getAllDivsBetweenDates(lastDate, thisDate));
                            DoubleArraySP divTimes(new DoubleArray(assetBridgeDivs->getArray().size()));
                            DateTimeArrayConstSP exDivDates = assetBridgeDivs->getExDivDates();
                            assetTimeMetrics[iAsset]->yearFrac(lastDate, *exDivDates, *divTimes);

                            HitNoHitBBSP bb(new HitNoHitBB(barrierStart,
                                                           barrierEnd,
                                                           driverVar,
                                                           bridgeLength,
                                                           assetBarrier->isUp,
                                                           assetBridgeDivs->getArray(),
                                                           *divTimes));

                            bbData->assetData[iAsset] = BridgeData::AssetBridgeDataSP(new
                                BridgeData::AssetBridgeData(iAsset,
                                                            barrierStart,
                                                            barrierEnd,
                                                            refLevelSV,
                                                            bb,
                                                            bridgeLength,
                                                            assetBridgeDivs,
                                                            divTimes));

                            // Create offsets and BB object
                            offsets[iBarrier][iAsset].push_back(iStep);
                        }
                    }
                }
                bridgeData[iStep] = bbData;
            }

            // Count on how many dates we need to do something
            numBBs = 0;
            for(iStep = 0; iStep < maxNumBBs; iStep++) {
                for(iAsset = 0; iAsset < numAssets; iAsset++) {
                    if(doHitNoHit[iAsset][iStep]) {
                        numBBs++;
                        break;
                    }
                }
            }

            // Populate the map of stateVarGens, StateVars
            for(iBarrier = 0; iBarrier < numBarriers; iBarrier++) {
                const SVGenBarrierHVBB* barrierBBGen = barrierHVBBGens[iBarrier];
                svDBase.append(barrierBBGen, SVGenBarrierHVBB::StateVarSP(new
                    SVGenBarrierHVBB::StateVar(pastPathGen,
                                               futPathGen,
                                               barrierBBGen->getBarrierData(),
                                               barrierBBGen->getRefLevelGen(),
                                               barrierBBGen->getSpotSmoothGen(),
                                               barrierBBGen->getSpotAtMonStartGen(),
                                               barrierBBGen->getSpotAtMonEndGen(),
                                               offsets[iBarrier],
                                               hitNoHitPath)));
            }
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Contains N-factor data for BB between 2 dates */
    class BridgeData {
    public:
        /** Contains 1-factor data for BB between 2 dates */
        class AssetBridgeData {
        public:
            /** Full constructor */
            AssetBridgeData(int iAsset,
                            double barrierStartPct,
                            double barrierEndPct,
                            IRefLevel::IStateVarSP  refLevelSV,
                            HitNoHitBBSP bbSample,
                            double bridgeLength,
                            DividendListSP divList,
                            DoubleArraySP divTimes):
            iAsset(iAsset), barrierStartPct(barrierStartPct), barrierEndPct(barrierEndPct),
            refLevelSV(refLevelSV), bbSample(bbSample), bridgeLength(bridgeLength),
            divList(divList), divTimes(divTimes) {
                // Compute weight for forward problem
                const DividendArray& divArray = divList->getArray();
                double meanTime = 0.0;
                double sumDivs  = 0.0;
                for(int iDiv = 0; iDiv < divArray.size(); iDiv++) {
                    double tDiv = (*divTimes)[iDiv];
                    sumDivs  += divArray[iDiv].getDivAmount();
                    meanTime += divArray[iDiv].getDivAmount() * tDiv;
                }
                if(!Maths::isZero(sumDivs) && bridgeLength > 0.0) {
                    meanTime /= sumDivs;
                    fwdWeight = 1.0 - meanTime / bridgeLength;
                } else {
                    // Just use the forward problem
                    fwdWeight = 1.0;
                }
            }

            /** Destructor */
            virtual ~AssetBridgeData() {}

            int                     iAsset;             //!< Asset index
            double                  barrierStartPct;    //!< Barrier % at start of BB
            double                  barrierEndPct;      //!< Barrier % at end of BB
            IRefLevel::IStateVarSP  refLevelSV;         //!< Reference level state variable
            HitNoHitBBSP            bbSample;           //!< Reverse Brownian Bridge
            double                  bridgeLength;       //!< Length of bridge
            DividendListSP          divList;            //!< Asset dividends for BB
            DoubleArraySP           divTimes;           //!< Dividend time relative to start of bridge
            double                  fwdWeight;          //!< Weight for forward problem

        private:
            /** Disabled default constructor */
            AssetBridgeData();
        };

        typedef refCountPtr<AssetBridgeData> AssetBridgeDataSP;


        /** Full constructor */
        BridgeData(const DateTime& startDate,
                   const DateTime& endDate,
                   int numAssets):
        startDate(startDate), endDate(endDate), assetData(numAssets) {}

        /** Destructor */
        virtual ~BridgeData() {}

        DateTime                  startDate;          //!< Start date for the BB
        DateTime                  endDate;            //!< End date for the BB
        vector<AssetBridgeDataSP> assetData;          //!< Brownian Bridge Per asset data

    private:
        /** Disabled default constructor */
        BridgeData();
    };

    typedef refCountPtr<BridgeData> BridgeDataSP;

    // Fields
    int                            numBarriers;        //!< Number of barriers
    int                            numAssets;          //!< Number of assets
    const IMultiFactors*           mAsset;             //!< Multi asset interface
    int                            maxNumBBs;          //!< Maximum nb of BBs = driver dates - 1
    int                            numBBs;             //!< Number of actual BBs
    DependenceSP                   dependence;         //!< Dependence object e.g. gaussian
    MCRandomGenCacheSP             randomGen;          //!< Random number generator for normals

    MCHitValuePathCacheSP          hitValuePathCache;  //!< Hit value path cache
    MCRandomCacheSP                randCache;          //!< Random number cache

    vector<IntArray>               hitNoHitPath;       //!< HitNoHit flag for BB [iAsset][iStep]
    vector<BoolArray>              doHitNoHit;         //!< Whether HitNoHit is required [iAsset][iStep]

    PathGenSpot::DriverConstSP     driver;             //!< Driver of pathGenSpot
    vector<BridgeDataSP>           bridgeData;         //!< Data for Brownian Bridges [iStep]

    // tabulated normal distribution. Use fast non virtual inteprolant
    EquidistantLinearInterpolantNonVirtualSP tabulatedNorm; //!< Tabulated normal distribution

    static const MCCache::KeyUtils::Key keyHitValueCache;   //!< Key for hit no hit path cache
    static const MCCache::KeyUtils::Key keyRandCache;       //!< Key for hit no hit path random numbers cache
    static const double ADJUST_CONSTANT;                    //!< Barrier adjustment constant
};

// Initialize cache keys
const MCCache::KeyUtils::Key MCPathConfigLN::PathGenHVBB::keyHitValueCache = MCCache::KeyUtils::getNewKey();
const MCCache::KeyUtils::Key MCPathConfigLN::PathGenHVBB::keyRandCache = MCCache::KeyUtils::getNewKey();
const double MCPathConfigLN::PathGenHVBB::ADJUST_CONSTANT = 0.5826;


/////////////////////////////////////////////////////////////////////////////


/** Class responsible for barrier adjustment */
class MCPathConfigLN::BarrAdj {
public:
    static const double ADJUST_CONSTANT;
    static const double BUS_DAYS_IN_YEAR;
    static const double SQRT_ONE_DAY_YEAR_FRAC;

    BarrAdj(const SVGenBarrierHVBarAdjArray& barrierHVBarAdjGenArray,
            IStateVariableGen::IStateGen*        pastPathGen,
            IStateVariableGen::IStateGen*        futPathGen,
            const MCProductClient*           prodClient,
            StateVarDBase&                   svDBase) {

        static const string routine = "MCPathConfigLN::BarrAdj::BarrAdj";

        const IMultiFactors* mAsset = prodClient->getMultiFactors();
        const DateTime& today = prodClient->getToday();

        // Obtain driver from path gen spot
        // PathGenSpot::DriverConstSP driver = pathGenSpot->getDriver();
        // const DateTimeArray driverDates = driver->simTimeline;

        // Loop over all barrier objects
        unsigned int iVar;
        HolidaySP hols(Holiday::weekendsOnly());
        for(iVar = 0; iVar < barrierHVBarAdjGenArray.size(); iVar++) {
            // Extract the barrier data and reference level
            const SVGenBarrierHVBarAdj* barrierGen = barrierHVBarAdjGenArray[iVar];
            BarrierDataConstSP data = barrierGen->getBarrierData();

            // Obtain reference level state variable
            IRefLevel::IStateVarGenSP refLevelGen = barrierGen->getRefLevel();
            IRefLevel::IStateVarSP oldStateVar = refLevelGen->getRefLevelSV(
                IRefLevel::IStateVarSP(   ), pastPathGen);
            IRefLevel::IStateVarSP refLevelSV =
                refLevelGen->getRefLevelSV(oldStateVar, futPathGen);

            // Deduce fwd starting flags
            const DateTime& startDate = refLevelGen->getAllDates().front();
            bool fwdStarting = startDate.isGreater(today);

            const DateTime& smoothDate = data->getSmoothingData()->smoothDate;

            // Loop over assets
            for(int iAsset = 0; iAsset < data->numAssets; iAsset++) {
                // Deconst asset barrier to adjust its value
                BarrierPerAssetDataSP assetBarrier =
                    CONST_POINTER_CAST<BarrierPerAssetData>(data->getAssetBarrier(iAsset));
                const DateTimeArray& monDates = *assetBarrier->monitoringDates;
                DoubleArray& barrierLevels = assetBarrier->barrierLevels;
                int UoD = assetBarrier->isUp? -1 : 1;
                double refLevel = refLevelSV->refLevel(iAsset);

                // Let's make sure that smoothing date is one of the monitoring dates
                // otherwise it's not clear what is meant by using a shifted barrier
                // at the smoothing date
                // Find smoothing date and replace smoothing barrier
                int iStep, smoothIdx;
                int nbDates = monDates.size();
                for(iStep = 0; iStep < nbDates; iStep++) {
                    if(smoothDate == monDates[iStep]) {
                        smoothIdx = iStep;
                        break;
                    }
                }

                if(iStep == nbDates) {
                    throw ModelException(
                        routine,
                        "Smoothing date must be one of the monitoring dates.");
                }

                // const DoubleArray& driverVols = driver->driverVols[iAsset];
                // Loop over monitoring dates and shift the barrier
                DoubleArray vols(nbDates, 0.0);
                DoubleArray days(nbDates, 0.0);

                // Create vol request
                double interpLevel = 1.0;
                LinearStrikeTSVolRequestSP volRequest(new LinearStrikeTSVolRequest(
                    interpLevel, startDate, monDates.back(), fwdStarting));

                for(iStep = 0; iStep < nbDates - 1; iStep++) {
                    // Compute the vol by integrating the driver's variance
                    // between the desired dates
                    if (today >= monDates[iStep + 1]) {
                        vols[iStep + 1] = 0.0;
                    } else {
                        interpLevel = barrierLevels[iStep + 1];
                        if (!fwdStarting) {
                            interpLevel *= refLevel;
                        }
                        volRequest->setStrike(interpLevel);
                        CVolProcessedSP vol(mAsset->factorGetProcessedVol(iAsset, volRequest.get()));
                        CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));

                        TimeMetricConstSP metric = vol->GetTimeMetric();
                        double vol2 = volBS->CalcVol(today, monDates[iStep + 1]);
                        double t2 = metric->yearFrac(today, monDates[iStep + 1]);
                        vols[iStep + 1]  = vol2;

                        if (today < monDates[iStep]) {
                            double t1 = metric->yearFrac(today, monDates[iStep]);
                            double vol1 = volBS->CalcVol(today, monDates[iStep]);

                            // Compute vol if there is trading time otherwise set to 0.0
                            double dt = t2 - t1;
                            if(Maths::isZero(dt)) {
                                vols[iStep + 1]  = 0.0;
                            } else {
                                vols[iStep + 1]  = sqrt((vol2*vol2*t2 - vol1*vol1*t1) / dt);
                            }
                        }
                    }

                    DateTime fromDate = ((today > monDates[iStep]))? today: monDates[iStep];
                    DateTime yesterday = DateTime(today.getDate()-1, today.getTime());

                    days[iStep + 1] = hols->businessDaysDiff(
                        (today > monDates[iStep])? yesterday: fromDate, monDates[iStep+1]);

                }

                for(iStep = 0; iStep < nbDates; iStep++) {
                    barrierLevels[iStep] *= exp(UoD * ADJUST_CONSTANT * vols[iStep] *
                        sqrt(days[iStep] / BUS_DAYS_IN_YEAR));

                    if (CString::equalsIgnoreCase(data->monitorType,
                        BarrierData::DAILY_MONITORING)) {

                        // for daily adjust back from adjusted-cts to daily interval
                        barrierLevels[iStep] *= exp(-UoD * ADJUST_CONSTANT * vols[iStep] *
                            SQRT_ONE_DAY_YEAR_FRAC);
                    }
                }
                // Adjust barrier at smoothing date
                assetBarrier->barrierSmoothDate = barrierLevels[smoothIdx];
            }

            SVGenBarrierHVEur::StateVarSP barrierEuropean(
                new SVGenBarrierHVEur::StateVar(
                data, barrierGen->getRefLevel(), barrierGen->getSpotPath(),
                barrierGen->getSmoothPath()));
            barrierEuropean->update(pastPathGen);
            barrierEuropean->update(futPathGen);
            svDBase.append(barrierGen, barrierEuropean);
        }
    }
};

const double MCPathConfigLN::BarrAdj::ADJUST_CONSTANT = 0.5826;
const double MCPathConfigLN::BarrAdj::BUS_DAYS_IN_YEAR = 262.0;
const double MCPathConfigLN::BarrAdj::SQRT_ONE_DAY_YEAR_FRAC = sqrt(1.0 / BUS_DAYS_IN_YEAR);


/////////////////////////////////////////////////////////////////////////////


MCPathConfigLN::Gen::Gen(MCPathConfigLN*          mcPathConfigLN,
                         const MCPathGeneratorSP& pastPathGenerator,
                         const MCProductClient*   prodClient,
                         MCCacheManager&          cacheMgr,
                         bool                     cachingRequested,
                         DateTimeArray&           simDates):
nowPathIdx(0) {
    static const string routine = "MCPathConfigLN::Gen::Gen";

    try {
        // Collect state variables from product and categorize them
        StateVariableCollectorSP svCollector(new StateVariableCollector());
        prodClient->collectStateVars(svCollector);
        IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

        // 1) Spot requests
        SVGenSpotArray spotGenArray = filterStateVars<SVGenSpot>(stateVarGenArray);
        if(!spotGenArray.size()) {
            throw ModelException("No spot paths specified.");
        }

        // 2a) Brownian Bridge Barrier hit value requests
        SVGenBarrierHVBBArray barrierHVBBGens =
            filterStateVars<SVGenBarrierHVBB>(stateVarGenArray);

        // 2b) Barrier Adjustment hit value requests
        SVGenBarrierHVBarAdjArray barrierHVBarAdjGenArray =
            filterStateVars<SVGenBarrierHVBarAdj>(stateVarGenArray);

        if(barrierHVBBGens.size() && barrierHVBarAdjGenArray.size()) {
            throw ModelException(
                "Cannot support simultaneously Brownian "
                "Bridges and Barrier Adjustment methodologies. ");
        }

        // 3) Create discount factor state variables
        SVGenDiscFactorArray discFactors =
            filterStateVars<SVGenDiscFactor>(stateVarGenArray);
        unsigned int iVar;
        for(iVar = 0; iVar < discFactors.size(); iVar++) {
            IStateVariableSP sv(discFactors[iVar]->
                             determinsticSV(false /* not doing past */));
            svDBase.append(discFactors[iVar], sv);
        }

        // 4) Create expected discount factor (ZCBs) state variables
        vector<const SVGenExpectedDiscFactor*>  expDiscFactors(
            filterStateVars<SVGenExpectedDiscFactor>(stateVarGenArray));
        const DateTime& today = prodClient->getToday();
        for(iVar = 0; iVar < expDiscFactors.size(); iVar++) {
            IStateVariableSP sv(expDiscFactors[iVar]->determinsticSV(
                                 today, false /* not doing past */));
            svDBase.append(expDiscFactors[iVar], sv);
        }

        // Nothing else is supported
        if(!stateVarGenArray.empty()) {
            throw ModelException("Unable to recognize all state variable types.");
        }

        // Create driver generators for BarrierHVs
        DateTimeArraySP allDriverDates =
            SVGenBarrierHVBB::BarrierDatesHelper::getBarrierTimeline(
            prodClient, barrierHVBBGens, simDates);

        // Create a spot path generator for the SVGenSpot and MCDriver
        PastPathGen* pastPathGen = dynamic_cast<PastPathGen*>(pastPathGenerator.get());
        if(!pastPathGen) {
            throw ModelException("Past path generator is not of PastPathGen type.");
        }
        pathGenSpot = MCPathConfigLN::PathGenSpotSP(new
            MCPathConfigLN::PathGenSpot(mcPathConfigLN,
            pastPathGen->getPathGenSpot(), prodClient, spotGenArray,
            *allDriverDates, cacheMgr, cachingRequested, svDBase, simDates));

        // Create a hit value path generator. Pass the driverGens and the
        // pathGenSpot so that the HV simulator can obtain the driver SVs
        if(barrierHVBBGens.size()) {
            pathGenHVBB = MCPathConfigLN::PathGenHVBBSP(new
            MCPathConfigLN::PathGenHVBB(mcPathConfigLN, prodClient,
            pastPathGen, this, pathGenSpot, barrierHVBBGens,
            cacheMgr, cachingRequested, svDBase));
        }

        // Create barrier adjustment state variables.
        if(barrierHVBarAdjGenArray.size()) {
            barrAdj = MCPathConfigLN::BarrAdjSP(new
            MCPathConfigLN::BarrAdj(barrierHVBarAdjGenArray,
            pastPathGen, this, prodClient, svDBase));
        }

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


int MCPathConfigLN::Gen::NbSimAssets() const {
    static const string routine = "MCPathConfigLN::Gen::NbSimAssets";
    throw ModelException(routine, "Method is retired for StateVars");
}

const double* MCPathConfigLN::Gen::Path(int iAsset, int iPath) const {
    static const string routine = "MCPathConfigLN::Gen::Path";
    throw ModelException(routine, "Method is retired for StateVars");
};

double MCPathConfigLN::Gen::refLevel(int iAsset, int iPath) const {
    static const string routine = "MCPathConfigLN::Gen::refLevel";
    throw ModelException(routine, "Method is retired for StateVars");
}

double MCPathConfigLN::Gen::maxDriftProduct(int iAsset) const {
    static const string routine = "MCPathConfigLN::Gen::maxDriftProduct";
    try {
        return pathGenSpot->maxDriftProduct(iAsset);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

int MCPathConfigLN::Gen::begin(int iAsset) const {
    static const string routine = "MCPathConfigLN::Gen::begin";
    throw ModelException(routine, "Method is retired for StateVars");
}

int MCPathConfigLN::Gen::end(int iAsset) const{
    static const string routine = "MCPathConfigLN::Gen::end";
    throw ModelException(routine, "Method is retired for StateVars");
}

// Live methods
bool MCPathConfigLN::Gen::hasPast() const {
    return pathGenSpot->hasPast();
}

bool MCPathConfigLN::Gen::doingPast() const {
    return false;
}

void MCPathConfigLN::Gen::generatePath(int pathIdx) {
    nowPathIdx = pathIdx;
    pathGenSpot->generatePath(pathIdx);
    if(pathGenHVBB.get()) {
        pathGenHVBB->generatePath(pathIdx);
    }
}

void MCPathConfigLN::Gen::advance() {
    pathGenSpot->advance();
}

void MCPathConfigLN::Gen::reset() {
    pathGenSpot->reset();
}

int MCPathConfigLN::Gen::getPathIndex() const {
    return nowPathIdx;
}

/** Returns the state variable corresponding to generator.
    Part of the IStateVariableGen::IStateGen IFace */
IStateVariableSP MCPathConfigLN::Gen::create(const IStateVariableGen* svGen) {
    static const string routine = "MCPathConfigLN::Gen::create";

    try {
        return svDBase.find(svGen);
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


/** Creates a past path generator */
MCPathGeneratorSP MCPathConfigLN::pastPathGenerator(const IMCProduct* prod) {
    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
    if(!prodClient){
        // Old approach
        return MCPathConfig::pastPathGenerator(prod);
    } else {
        // State variables approach
        return MCPathGeneratorSP(new PastPathGen(prodClient));
    }
}


/** Creates a PathGenerator used for future dates which supports
    'random access' to the paths. That is if a path generator is
    created using savePaths = true, then subsequently another path
    generator can be created (with savePaths = false) which will
    support the pathIdx parameter in the generatePath() method, ie the
    paths can be generated in any order. */
MCPathGeneratorSP MCPathConfigLN::futurePathGenerator(
    int                      cachingMode,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control,
    Results*                 results,
    DateTimeArray&           simDates ) {

    static const string routine = "MCPathConfigLN::futurePathGenerator";

    try {
        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
        if(!prodClient) {
            // Old approach
            return MCPathConfig::futurePathGenerator(cachingMode, numPaths,
                pastPathGenerator, prod, control, results, simDates);
        } else {
            pathCache = PathCacheSP(new PathCache());
            randomCache = RandomCacheSP(new RandomCache());

            // State variables approach
            bool isPricing = control->isPricing();
            bool cachePaths = (cachingMode & PATH_CACHE)? true: false;
            // cachePaths => cacheRandoms, so no need to save states
            bool mySaveStates = !cachePaths && (cachingMode & PATH_RANDOM_ACCESS);
            /* no point saving paths if no greeks */
            if (mySaveStates){
                throw ModelException("MCPathConfigLN::futurePathGenerator",
                                     "QuickGreeks are not supported without path caching");
            }

            if (cachePaths || isPricing) {
                // Kill them
                cacheMgr = MCCacheManager();
            } else {

                const IMultiFactors* multiFactor = prod->getMultiFactors();
                int numAssets = prod->getNumAssets();
                SensitivitySP sens(control->getCurrentSensitivity());
                // getSensitiveAssets() does the real work
                IntArray sensitivePhiAssets(multiFactor->getSensitiveAssets(
                                            sens.get(), true)); // include phi etc
                if (sensitivePhiAssets.empty()){
                    IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                                             sens.get(), false)); // exclude phi etc
                    cacheMgr.disallowAssets(sensitiveAssets);
                } else {
                    // mark all paths as invalid
                    IntArray sensitiveAssets(numAssets);
                    for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                        sensitiveAssets[iAsset] = iAsset;
                    }
                    cacheMgr.disallowAssets(sensitiveAssets);
                }
            }

            // Bit comparison to deduce whether user has asked for caching
            bool cachingRequested = (cachingMode & IMCPathConfig::PATH_CACHE) ? true : false;

            MCPathGeneratorSP gen = makePathGenerator(
                cachingRequested, numPaths, pastPathGenerator, prod, control, results, simDates);

            if (cachePaths || isPricing) {
                // Caches have been built. Configure them now
                cacheMgr.configureCache(numPaths);
            }

            return gen;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigLN::makePathGenerator(
    bool                     cachingRequested,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control,
    Results*                 results,
    DateTimeArray&           simDates ){

    static const string routine = "MCPathConfigLN::makePathGenerator";

    try {
        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
        if(!prodClient){
            // Old approach
            MCPathBaseSP futurePathGenBaseLN(
                new Generator(this, pastPathGenerator, prod));
            return MCPathBase::createPathGenerator(pastPathGenerator,
                                                   this,
                                                   prod,
                                                   futurePathGenBaseLN);
        } else {
            // State variables approach
            return MCPathGeneratorSP(new MCPathConfigLN::Gen(
                this,
                pastPathGenerator,
                prodClient,
                cacheMgr,
                cachingRequested,
                simDates));
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


int MCPathConfigLN::storagePerPath(IMCProduct* product) const{
    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(product);
    if(!prodClient){
        // Call parent
        return MCPathConfig::storagePerPath(product);
    } else {
        // State variables approach

        if (!product->hasFuture()) {
            return 0; // nothing to simulate
        }

        // Let's just create one on the fly...
        smartPtr<MCPathConfigLN> pathConfig(copy(this)); // copy to avoid const problems
        MCPathGeneratorSP pastPathGenerator(pathConfig->pastPathGenerator(product));
        product->pathGenUpdated(pastPathGenerator.get());
        SensitivityArrayConstSP    sens;
        OutputRequestArrayConstSP  request;
        Control control(sens, request, false, "");
        Results results;
        DateTimeArray simDates;
        MCPathGeneratorSP pathGen(pathConfig->futurePathGenerator(2, // cache
                                                                  2, // num paths
                                                                  pastPathGenerator,
                                                                  product,
                                                                  &control,
                                                                  &results,
                                                                  simDates ));

        int memRequirement = pathConfig->cacheMgr.storagePerRun();
        return memRequirement;
    }
}

/** Does special things for Theta-type tweaks */
bool MCPathConfigLN::sensShift(Theta* shift) {
    MCPathConfig::sensShift(shift); // call parent's method
    cacheMgr = MCCacheManager();    // Blank cache manager now

    return true;
}


#endif

// 3 problems with reporting standard error for greeks
// (i) How to request the output
// (ii) Where to store the output
// (iii) How to ensure it is 'aggregated' properly ie if sim is split into
// blocks you will get the wrong numbers
// Current status: Standard error for greeks always calculated and stored in
// a debug output packet. MonteCarlo::Run sorts out aggregation across blocks

#define xLR_TEST /* slow way of doing same calculation - just in code once
                    so it's in RCS for reference */
class MCPathConfigLN::LRDeltaCommon{
private:
    Generator*             gen;
    DoubleArray            cachedDeltaRatios;
    int                    pathIdx; // initially set to -1
    DoubleMatrix           corrInverse; /* inverse of correlation matrix */
    DoubleArray            spotsAtSimStart;
    DoubleArray            sqrtVars; // square root of vars to 1st time point
public:
    static LRDeltaCommon* create(const MCPathGeneratorSP& pathGen){
        static const string method("MCPathConfigLN::LRDeltaCommon");
        // get hold of our generator
        MCPathBaseSP pathBase =
            dynamic_cast<MCPathBase::Generator&>(*pathGen).getPathBase();
        Generator* gen = &dynamic_cast<MCPathConfigLN::Generator&>(*pathBase);
        if (gen->timeline->numFutRefLevel > 0){
            throw ModelException(method,
                                 "Likelihood Ratio not support for forward"
                                 " starting/averaging in");
        }
        // futureDates always contains today/sim start date. This needs
        // reviewing/fixing if numFutRefLevel > 0
        if (gen->timeline->futureDates.size() < gen->timeline->numFutRefLevel+2){
            throw ModelException(method, "Internal error - too few dates");
        }
        if (gen->timeline->today.yearFrac(gen->timeline->futureDates[gen->timeline->numFutRefLevel+1]) <
            (double)(DateTime::VOL_TIME_IN_DAY - 1)/
            (double)DateTime::VOL_TIME_IN_YEAR){
            // This is somewhat arbitrary. We could just fail but that is
            // awkward in prod. cf Fwd Start case above. Also want clear cut
            // case when algorithm switches
            return 0; // variance too small - flag that we can't do it
        }

        DoubleArray  sqrtVars(gen->numAssets);
        // just pick information needed at first sim date in future
        for (int i = 0; i < sqrtVars.size(); i++){
            const DoubleMatrix& assetVols = gen->vols[i];
            // just using 1st vol request (per asset)
            sqrtVars[i] = assetVols[0][0];
            if (Maths::isZero(sqrtVars[i])){
                // variance really is too small - flag that we can't do it
                return 0;
            }
        }

        return new LRDeltaCommon(gen, sqrtVars);
    }
private:
    // constructor
    LRDeltaCommon(Generator*                gen,
                  const DoubleArray&        sqrtVars):
        gen(gen), pathIdx(-1), sqrtVars(sqrtVars){
        static const string method("MCPathConfigLN::LRDeltaCommon");
        int numAssets = gen->numAssets;
        cachedDeltaRatios = DoubleArray(numAssets);
        spotsAtSimStart = gen->refData->fwdsAtSimStart;
        DoubleMatrix correlations(gen->dependence->getCorrelations(0));
        corrInverse = correlations.computeInverse();
    }
public:
    /** Must be invoked before any deltaRatio/gammaRatio etc methods are
        called. Allows common calculations to be cached */
    void moveToNextPath(int newPathIdx){
        const DoubleMatrix& randoms = gen->randomGen->getRandomNumbers();
        if (newPathIdx != pathIdx){
            pathIdx = newPathIdx; // update
            // Then calculate likelihood ratio for delta. The
            // derivation of the formula is very complex (see docs)
            // but amazingly collapses to that given below
            int numAssets = spotsAtSimStart.size(); // for ease
            for (int idx = 0; idx < numAssets; idx++){
                double ratio = 0.0;
                const double* matrixCol = corrInverse[idx];
                for (int iAsset = 0; iAsset < numAssets; iAsset++){
                    ratio += matrixCol[iAsset] * /*gen->*/randoms[iAsset][0];
                }
                ratio /= spotsAtSimStart[idx] * sqrtVars[idx];
                cachedDeltaRatios[idx] = ratio;
            }
        }
    }
    // how many assets are being simulated
    int numAssets() const {
        return gen->numAssets;
    }

    /** How to turn delta for asset into delta for underlying whose name
        is returned from getOutputName() */
    double assetDeltaScalingFactor(int assetIdx) const{
        const CAsset& asset = gen->mAsset->getAsset(assetIdx);
        if (CAsset::IStruck::TYPE->isInstance(&asset)){
            const CAsset::IStruck& struckAsset =
                dynamic_cast<const CAsset::IStruck&>(asset);
            return struckAsset.getFXSpot();
        }
        return 1.0;
    }
    //// returns new OutputName(assetGetTrueName()) on selected asset
    OutputNameSP getOutputName(int assetIdx) const{
        return OutputNameSP(
            new OutputName(gen->mAsset->assetGetTrueName(assetIdx)));
    }

    //// computes the likelihood ratio for delta wrt specified asset
    double deltaRatio(int idx) const{
        return cachedDeltaRatios[idx];
    }
    double gammaRatio(int idx) const{
        double deltaR = deltaRatio(idx);
        double gammaR = deltaR * (deltaR - 1.0/spotsAtSimStart[idx]) -
            corrInverse[idx][idx]/(Maths::square(spotsAtSimStart[idx] *
                                                 sqrtVars[idx]));
        return gammaR;
    }

    double crossGammaRatio(int idx1, int idx2) const{
        double deltaR1 = deltaRatio(idx1);
        double deltaR2 = deltaRatio(idx2);
        double xGammaR = deltaR1 * deltaR2 - corrInverse[idx1][idx2] /
            (spotsAtSimStart[idx1] * sqrtVars[idx1] *
             spotsAtSimStart[idx2] * sqrtVars[idx2]);
        return xGammaR;
    }
};

/** Handles computation of delta & gamma via LR */
class  MCPathConfigLN::LRDelta: public LRGenerator::GreekCalculator{
private:
    refCountPtr<LRDeltaCommon>       deltaCommon;
    vector<int>                      assets; // which assets to do
    vector<IMCPricesSimpleSP>        prices;
    double                           origShiftSize;
public:
    virtual ~LRDelta(){}
    LRDelta(const refCountPtr<LRDeltaCommon>& deltaCommon,
            vector<int>&                 assets, // what the deriv is wrt
            vector<IMCPricesSimpleSP>& prices,// preallocated: 2*assets.size
            double                       origShiftSize):
        deltaCommon(deltaCommon), assets(assets),
        prices(prices), origShiftSize(origShiftSize) {}

    /** Calculate contrubution for greek from this path */
    virtual void processPath(int pathIdx, double price){
        bool zeroPrice = Maths::isZero(price);
        if (!zeroPrice){
            deltaCommon->moveToNextPath(pathIdx);
        }
        int numAssets = assets.size();
        for (int i = 0; i < numAssets; i++){
            int assetIdx = assets[i];
            prices[i]->add(zeroPrice?
                           0.0: price * deltaCommon->deltaRatio(assetIdx));
            prices[numAssets+i]->add(zeroPrice?
                                     0.0:
                                     price * deltaCommon->gammaRatio(assetIdx));
        }
    }
    /** Store the result for the greek corresponding to this LRGenerator */
    virtual void storeResult(Results*   results,
                             double     pvFactor){
        results->removePacket(Delta::NAME); // needed if Turbo model
        results->removePacket(Delta::SECOND_ORDER_NAME);// needed if Turbo model
        for (unsigned int i = 0; i < assets.size(); i++){
            double deltaValue, deltaStdErr;
            // do delta
            prices[i]->getResult(&deltaValue, &deltaStdErr);
            double deltaScalingFactor = deltaCommon->assetDeltaScalingFactor(i);
            deltaValue *= pvFactor * deltaScalingFactor;
            deltaStdErr *= pvFactor * deltaScalingFactor;
            OutputNameSP name(deltaCommon->getOutputName(assets[i]));
            if (results->exists(Delta::NAME, name)){
                throw ModelException("MCPathConfigLN::LRDelta::storeResult",
                                     "Delta for "+name->toString()+" already "
                                     "calculated. Duplicate assets "
                                     "in simulation?");
            }
            results->storeScalarGreek(deltaValue, Delta::NAME, name);
            /* This is a bit questionable, but if delta/gamma is calculated
               using LR, but x gamma from external tweaking then we must
               store a 'shift size' in order for that algorithm to work */
            // store divisor used (needed for cross derivatives)
            results->storeScalarGreek(origShiftSize,
                                      Delta::NAME+Results::SHIFT_SIZE_POSTFIX,
                                      name);
             // do gamma
            double gammaValue, gammaStdErr;
            prices[assets.size()+i]->getResult(&gammaValue, &gammaStdErr);
            gammaValue *= pvFactor * Maths::square(deltaScalingFactor);
            gammaStdErr *= pvFactor * Maths::square(deltaScalingFactor);
            results->storeScalarGreek(gammaValue, Delta::SECOND_ORDER_NAME,
                                      name);
            results->storeScalarGreek(deltaStdErr,
                                      MonteCarlo::DEBUG_GREEK_STD_ERR+
                                      Delta::NAME,
                                      name);
            results->storeScalarGreek(gammaStdErr,
                                      MonteCarlo::DEBUG_GREEK_STD_ERR+
                                      Delta::SECOND_ORDER_NAME,
                                      name);
        }
    }
};

class  MCPathConfigLN::LRCrossGamma: public LRGenerator::GreekCalculator{
private:
    refCountPtr<LRDeltaCommon>      deltaCommon;
    vector<int>                     assets1; // which assets to do
    vector<int>                     assets2; // which assets to do
    vector<IMCPricesSimpleSP>     prices;
public:
    virtual ~LRCrossGamma(){}
    //// Note: require assets1.size() = assets2.size()
    LRCrossGamma(const refCountPtr<LRDeltaCommon>&  deltaCommon,
                 vector<int>&                 assets1, // what the deriv is wrt
                 vector<int>&                 assets2, // what the deriv is wrt
                 vector<IMCPricesSimpleSP>& prices):// preallocated:
    deltaCommon(deltaCommon),
    assets1(assets1), assets2(assets2), prices(prices) {}

    /** Calculate contrubution for greek from this path */
    virtual void processPath(int pathIdx, double price){
        bool zeroPrice = Maths::isZero(price);
        if (!zeroPrice){
            deltaCommon->moveToNextPath(pathIdx);
        }
        int numAssets = assets1.size();
        for (int i = 0; i < numAssets; i++){
            prices[i]->add(zeroPrice?
                          0.0:
                          price * deltaCommon->crossGammaRatio(assets1[i],
                                                               assets2[i]));
        }
    }
    /** Store the result for the greek corresponding to this LRGenerator */
    virtual void storeResult(Results*   results,
                             double     pvFactor){
        results->removePacket(CrossGamma::NAME); // needed if Turbo model
        for (unsigned int j = 0; j < assets1.size(); j++){
            OutputNameSP name1(deltaCommon->getOutputName(assets1[j]));
            OutputNameSP name2(deltaCommon->getOutputName(assets2[j]));
            OutputNameSP name(new OutputName(name1.get(), name2.get()));
            double value, stdErr;
            prices[j]->getResult(&value, &stdErr);
            double deltaScalingFactor1 =
                deltaCommon->assetDeltaScalingFactor(assets1[j]);
            double deltaScalingFactor2 =
                deltaCommon->assetDeltaScalingFactor(assets2[j]);
            value *= pvFactor * deltaScalingFactor1 * deltaScalingFactor2;
            stdErr *= pvFactor * deltaScalingFactor1 * deltaScalingFactor2;
            for (int i = 0; i < 2; i++){
                results->storeScalarGreek(value, CrossGamma::NAME, name);
                results->storeScalarGreek(stdErr,
                                          MonteCarlo::DEBUG_GREEK_STD_ERR+
                                          CrossGamma::NAME,
                                          name);
                name = OutputNameSP(new OutputName(name2.get(), name1.get()));
            }
       }
    }
};

class MCPathConfigLN::MyLRGenerator: public LRGenerator{
    // fields
    refCountPtr<LRDeltaCommon> deltaCommon;
    int                        nbIter;
    int                        nbSubSamples;

    vector<int> getAssetIdxs(const Sensitivity* sens){
        // to do - handle names listed in sensitivity
        vector<int> assets(deltaCommon->numAssets());
        for (unsigned int i = 0; i < assets.size(); i++){
            assets[i] = i;
        }
        return assets;
    }
    vector<IMCPricesSimpleSP> createPrices(int numPrices){
        vector<IMCPricesSimpleSP> prices(numPrices);
        for (unsigned int i = 0; i < prices.size(); i++){
            prices[i] = IMCPricesSimpleSP(
                new MCPricesSimple(nbIter, nbSubSamples));
        }
        return prices;
    }
public:
    MyLRGenerator(LRDeltaCommon*           deltaCommon,
                  int                      nbIter,
                  int                      nbSubSamples):
        deltaCommon(deltaCommon), nbIter(nbIter), nbSubSamples(nbSubSamples) {}

    /** Returns object which handles greek calculation for specified
        sensitivity */
    GreekCalculator* greekCalculator(const Sensitivity* sens,
                                     bool               calcStdErr){
        if (!deltaCommon){
            return 0; // too close to sample date - flag that we can't do it
        }
        if (calcStdErr){
            // see comments above for HACK_FOR_STD_ERR for reason why
            throw ModelException("MCPathConfigLN::greekCalculator", "Standard "
                                 "error for greeks not supported");
        }
        CClassConstSP clazz = sens->getClass();
        if (clazz == Delta::TYPE){
            vector<int> assets(getAssetIdxs(sens));
            vector<IMCPricesSimpleSP> prices(createPrices(2 * assets.size()));
            const Delta& delta = dynamic_cast<const Delta&>(*sens);
            return new LRDelta(deltaCommon, assets,
                               prices, delta.getShiftSize());
        }
        if (clazz == CrossGamma::TYPE){
            // to do - handle names listed in sensitivity
            vector<int> assets(getAssetIdxs(sens));
            vector<int> assets1;
            vector<int> assets2;
            for (unsigned int i = 0; i < assets.size(); i++){
                for (unsigned int j = i+1; j < assets.size(); j++){
                    assets1.push_back(i);
                    assets2.push_back(j);
                }
            }
            vector<IMCPricesSimpleSP> prices(createPrices(assets1.size()));
            return new LRCrossGamma(deltaCommon, assets1, assets2, prices);
        }
        return 0;
    }
};

/** See comment on pure virtual declaration. */
LRGenerator* MCPathConfigLN::createLRGenerator(
    const MCPathGeneratorSP& pathGen,
    int                      nbIter,
    int                      nbSubSamples){
    LRDeltaCommon* deltaCommon = LRDeltaCommon::create(pathGen);
    // deltaCommon may be null
    return new MyLRGenerator(deltaCommon, nbIter, nbSubSamples);
}

/** "bread & butter" path config that can be captured in Pyramid using
    current IMS */
class MCPathConfigLNDefault: public MCPathConfigLN {
public:
    static CClassConstSP const TYPE;

private:
    MCPathConfigLNDefault():
        MCPathConfigLN(
            TYPE,
            IVolatilityBS::TYPE->getName(),
            "not used"){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigLNDefault, clazz);
        SUPERCLASS(MCPathConfigLN);
        EMPTY_SHELL_METHOD(defaultMCPathConfigLNDefault);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigLNDefault(){
        return new MCPathConfigLNDefault();
    }
};

CClassConstSP const MCPathConfigLNDefault::TYPE =
CClass::registerClassLoadMethod(
    "MCPathConfigLNDefault", typeid(MCPathConfigLNDefault), load);

/* what bounds to use when calculating max drift */
const double MCPathConfigLN::Generator::MAX_DRIFT_RAND_NUMBER = 3.0;

// external symbol to allow class to be forced to be linked in
bool MCPathConfigLNLoad(){
    return (MCPathConfigLN::TYPE != 0 && MCPathConfigLNDefault::TYPE != 0);
}

DRLIB_END_NAMESPACE
