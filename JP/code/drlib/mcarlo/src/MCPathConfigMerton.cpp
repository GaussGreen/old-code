//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigMerton.cpp
//
//   Author      : Oliver Brockhaus
//
//   Description : Monte Carlo path generator
//                 Since random number gen is here - so too must all 
//                 the work that is expected to be done with the same
//                 set of randoms. So, "iPath", "iAsset", "iStep"
//                 and at some stage this might also mean "iTweak" for
//                 internal sensitivities
//
//
//   Date        : April 2003
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/VolMerton.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/DependenceGauss.hpp"

#define CAREFUL_RAND_DEBUG 0
#if CAREFUL_RAND_DEBUG
#endif

DRLIB_BEGIN_NAMESPACE



/***************************************** 
MCPathConfigMerton
*****************************************/
/** Class splits its work into three smaller classes. This is possibly more
    code but it makes it much simpler and easier to develop */
class MCPathConfigMerton : public MCPathConfig {
private:
    // fields ////
    string                        volType;
    bool                          isCarefulRandoms;
    string                        dependenceType;
    bool                          isSimJumpDates;

    class Generator;
    friend class Generator;
public:
    static CClassConstSP const TYPE;

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const {
        return MarketDataFetcherSP(
            new MDFAssetVol(volType, VolSurface::TYPE->getName()));
    }

protected:
    MCPathConfigMerton(CClassConstSP clazz, const string& volType, const string& dependenceType):
        MCPathConfig(clazz), volType(volType), isCarefulRandoms(false), isSimJumpDates(false) {
    }

private:
    MCPathConfigMerton(): MCPathConfig(TYPE), volType("VolMerton"),
        isCarefulRandoms(false), dependenceType("not used"),
        isSimJumpDates(false)
    {}

    virtual bool vegaMatrixSupported() const { return true; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * Not sure if we will ever want this?  See IModel::wantsRiskMapping().
     */

    IModel::WantsRiskMapping wantsRiskMapping() const {
        return IModel::riskMappingIrrelevant;
    }

protected:
    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                     cachingRequested,
        int                      numPaths,
        const MCPathGeneratorSP& pastPathGenerator,
        const IMCProduct*         prod,
        Control*                 control, 
        Results*                 results,
        DateTimeArray&           simDates ); // defined below
    /*{

        MCFuturePathGenMertonSP futurePathGenBaseMerton(
            new MCPathConfigMerton::Generator(saveStates, numPaths, this, pastPathGenerator, prod));
        return MCPathBase::createPathGenerator(pastPathGenerator,
                                               this,
                                               prod,
                                               futurePathGenBaseMerton);        
    }*/
    
    virtual bool carefulRandoms() const { return isCarefulRandoms; }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigMerton, clazz);
        SUPERCLASS(MCPathConfig);
        EMPTY_SHELL_METHOD(defaultMCPathConfigMerton);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);  // for now 
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);
        FIELD(isCarefulRandoms, "isCarefulRandoms");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms); 
        FIELD(isSimJumpDates, "isSimJumpDates");
        FIELD_MAKE_OPTIONAL(isSimJumpDates);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigMerton(){
        return new MCPathConfigMerton();
    }
};

CClassConstSP const MCPathConfigMerton::TYPE = 
CClass::registerClassLoadMethod("MCPathConfigMerton", typeid(MCPathConfigMerton), load);

/** "bread & butter" path config that can be captured in Pyramid using 
     current IMS */

class MCPathConfigMertonDefault: public MCPathConfigMerton {
public:
    static CClassConstSP const TYPE;
    
private:
    MCPathConfigMertonDefault():
        MCPathConfigMerton(
            TYPE,
            ""/*IVolatilityBS::TYPE->getName()*/,
            "not used"){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigMertonDefault, clazz);
        SUPERCLASS(MCPathConfigMerton);
        EMPTY_SHELL_METHOD(defaultMCPathConfigMertonDefault);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigMertonDefault(){
        return new MCPathConfigMertonDefault();
    }
};

CClassConstSP const MCPathConfigMertonDefault::TYPE = 
CClass::registerClassLoadMethod("MCPathConfigMertonDefault", typeid(MCPathConfigMertonDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigMertonLoad(){
    return (MCPathConfigMerton::TYPE != 0 && MCPathConfigMertonDefault::TYPE != 0);
}

///////////////////////////////////////////////////////////////////////////////////
// Generator
///////////////////////////////////////////////////////////////////////////////////

//typedef smartConstPtr<Poisson> PoissonConstSP;
//typedef smartPtr<Poisson> PoissonSP;
//typedef smartConstPtr<VolProcessedMerton> VolProcessedMertonConstSP;
//typedef smartPtr<VolProcessedMerton> VolProcessedMertonSP;

/** Class containing the basic variable and methods that (hopefully) will 
be used by a lot of path generators. */
class MCPathConfigMerton::Generator: virtual public MCPathBase,
                                     virtual public IMCRandom::Callbacks,
                                     public DependenceMakerGauss::Support {
    const static double MAX_DRIFT_RAND_NUMBER; /* what bounds to use when
                                                  calculating max drift */
    vector<CDoubleArray>      drifts;         // x[numAssets] each with NbSteps rows
    vector<CDoubleArray>      vols;
    DoubleArray               maxDrifts;      // product of MAX(1, drifts)
public:
//    friend class MCSimpleFuturePathGeneratorMerton;
//    friend class MCFuturePathGeneratorMerton;
    Generator(int                      numSimPaths,
              MCPathConfigMerton*      mcPathConfigMerton, 
              const MCPathGeneratorSP& pastPathGenerator,
              const IMCProduct*         prod);
    
    /** Obtains timeline object from base */
    MCProductTimelineConstSP getTimeline() const {
        return timeline;
    }

    /** Returns number of assets */
    int NbSimAssets() const;

    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath);

    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const;

    /** Returns the path for iAsset, iPath */
    double* Path(int iAsset, int iPath);

    /** Returns the number of paths per asset */
    int nbPaths(int iAsset) const {
        return 1;
    }

    /** Returns if it is a single path generator or not */
    bool isSinglePath() const {
        return true;
    }

    /** Configures pathGen for antithetics */
    virtual void configureAntithetics();

    /** Configures pathGen for nonAntithetics. We need a better name */
    virtual void configureNonAntithetics();

    /** Draws random numbers */
    virtual void drawRandomNumbers(int pathIdx);

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const;

protected:
    // fields

    // Merton
    int                       maxNumJumps;
    double                    quantile;       // probability of having more than numJumps jumps
    IntArray                  numJumps;
    DoubleMatrix              randomJumps;    // numAssets cols x NbSteps rows
    DoubleMatrix              randomJumpTimes;// 1 row, maxNumJumps cols
    DoubleArray               dt;
    DoubleMatrix              correlationJumps;
    DependenceSP              dependenceJumps;
    PoissonSP                 dependenceJumpTimes;
    bool                      isSimJumpDates;

    // transient
    mutable TimeMetricSP     timeMetric;

    IMCRandomSP              randomGen;      // Random number generator for paths
    IMCRandomSP              randomGenJumps; //Random number generator for jump sizes
    DependenceSP             dependence;   // Dependence object
    MCProductTimelineSP      timeline;     // Product timeline
    RefLevelDataSP           refData;      // RefLevel data

    vector<DoubleArray>      fwds;         // fwds [iAsset][iStep]
    const IMultiFactors*     mAsset;       // Multi asset interface
    int                      numAssets;    // number of assets
    vector<DoubleArray>      productPaths; // [iAsset][iStep]

    IRandomSP                rand; // Keep it here to be used in # jumps
    bool                     isCarefulRandoms;

public:
    void generateDriftAndVols(const DateTimeArray& futurePathDates);
    
    vector<VolProcessedMertonSP> volProcessed;

protected:
    void generatePath(int pathIdx, int iAsset, int iPath);
    
    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0 */
    double maxDriftProduct(int iAsset) const{
        return maxDrifts[iAsset];
    }
};

/** Split the cases of same simulation dates per asset and different 
    simulation dates per asset. This is slightly more code but a lot
    simpler. */

///////////////////////////////////////////////////////////////////////////////////
// Generator, MCSimpleFuturePathGeneratorMerton, MCFuturePathGeneratorMerton
///////////////////////////////////////////////////////////////////////////////////
/***************************************** 
Generator
*****************************************/
MCPathConfigMerton::Generator::Generator(int                      numSimPaths,
                                         MCPathConfigMerton*      mcPathConfigMerton, 
                                         const MCPathGeneratorSP& pastPathGenerator,
                                         const IMCProduct*         prod):
drifts(numAssets), 
vols(numAssets), 
maxDrifts(numAssets, 1.0),
fwds(prod->getNumAssets()),
mAsset(prod->getMultiFactors()), 
numAssets(prod->getNumAssets()),
productPaths(prod->getNumAssets()),
rand(mcPathConfigMerton->getRandomGenerator()) {    
    static const string routine("Generator");
    try{
        rand->init();
        
        // Obtain product timeline
        timeline = getProductTimeline(prod, pastPathGenerator);

        // Obtain market data
        IntArray nbPaths(numAssets, 1);
        refData = getRefData(timeline, nbPaths, mAsset, prod, pastPathGenerator);

        // Paths
        int iAsset;
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            productPaths[iAsset] = DoubleArray(timeline->totalNumSteps);
        }


        //vector<VolProcessedMertonSP> volProcessed;
        volProcessed.resize(numAssets);
        VolRequestMertonSP volRequest (new VolRequestMerton());
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            // for each asset have an array of vol requests
            volProcessed[iAsset] 
                = VolProcessedMertonSP(dynamic_cast<VolProcessedMerton*>
                    (mAsset->factorGetProcessedVol(iAsset,volRequest.get()))); 
        }

        // Merton
        if( isSimJumpDates ) {
            // 99.99% quantile of number of jumps
            double quantile;
            volProcessed[0]->Quantile( 
                timeline->today,
                timeline->futureDates[timeline->futureDates.size()-1], 
                0.00001, 
                1000,
                &quantile,&maxNumJumps );
        }
        else {
            maxNumJumps = timeline->numFutSteps;
        }

        if( maxNumJumps>0 ) {
            randomJumps = DoubleMatrix(numAssets, Maths::min(maxNumJumps,timeline->numFutSteps));
        }
        if( maxNumJumps>0 ) {
            randomJumpTimes = DoubleMatrix(1,maxNumJumps);
        }
        numJumps = IntArray(timeline->numFutSteps);
            
        /* now we know how many paths per asset we can allocate
           space for them. Also note that we can't ask for
           the vol interps until the past has been computed */
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            drifts[iAsset]  = DoubleArray(timeline->numFutSteps);
            vols[iAsset]    = DoubleArray(timeline->numFutSteps);
        }
        correlationJumps = DoubleMatrix(numAssets, numAssets);

        // precompute drift and vols could possibly have a class
        // "process" which did this - for now leave as local functions
        generateDriftAndVols(timeline->futureDates);

        // set up Dependence
        if( (maxNumJumps>0) /*&& (theVol->crashSizeUncertainty>0)*/ ) {
            randomJumps = DoubleMatrix(numAssets, Maths::min(maxNumJumps,timeline->numFutSteps));
        }
        dependence = mcPathConfigMerton->dependenceMaker->createDependence(this);
        dependenceJumps     = DependenceSP(new Gauss(correlationJumps));
        dependenceJumpTimes = PoissonSP(new Poisson(volProcessed[0]->getCrashRate()));

        if( !isSimJumpDates ) {
            // preprocess...
            DoubleArray timeInYears(timeline->numFutSteps);
            int iStep;
            for( iStep=0; iStep<timeline->numFutSteps-1; iStep++ ) {
                timeInYears[iStep] = dt[iStep]-dt[iStep+1];
            }
            timeInYears[iStep] = dt[iStep];
            dependenceJumpTimes->numJumpsPreprocess( timeInYears, 100 );
        }

        // Initialize random number generator for paths. No callback to pathGen
        randomGen = IMCRandomSP(new MCRandomNoCache(
            0,
            dependence,
            rand,
            mcPathConfigMerton->carefulRandoms(),
            timeline->numFutSteps,
            numAssets,
            timeline->totalNumPastDates));

        // Initialize random number generator for jumps. Callback to pathGen
        randomGenJumps = IMCRandomSP(new MCRandomNoCache(
            this,
            dependenceJumps,
            rand,
            mcPathConfigMerton->carefulRandoms(),
            Maths::min(maxNumJumps,timeline->numFutSteps),
            numAssets,
            0));

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


/** Returns number of assets */
int MCPathConfigMerton::Generator::NbSimAssets() const {
    return numAssets;
}

/** Returns the reference level for iAsset, iPath */
double& MCPathConfigMerton::Generator::refLevel(int iAsset, int iPath) {
    return refData->refLevels[iAsset][0];
}
    

/** Returns the reflevel path */
IRefLevel::IMCPathSP& MCPathConfigMerton::Generator::refLevelPath() const {
    return refData->refLevelPath;
}


/** Returns the paths */
double* MCPathConfigMerton::Generator::Path(int iAsset, int iPath) {
    return &productPaths[iAsset][0];
}


/** simulates random numbers */
void MCPathConfigMerton::Generator::drawRandomNumbers(int pathIdx) {
    // WE NEED TO REVISIT THIS PIECE OF CODE
    
    // We need to find a robust way of aligning the different 
    // path generators for antithetics, careful randoms etc.
    // The notion of careful randoms is particularly tricky here
    // as we are wasting a random number of random numbers (dpending on
    // number of jumps).
    
    // Draw random numbers for paths and jump sizes
    randomGen->generate(pathIdx);
    randomGenJumps->generate(pathIdx);
    return;
}

// implement getGaussData since derivation from DependenceMakerGauss::Support
CDoubleMatrixConstSP MCPathConfigMerton::Generator::getGaussData() const {
    static const string method("MCPathConfigMerton::Generator::getGaussData");
    try {
        return mAsset->factorsCorrelationMatrix();
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Configures pathGen for antithetics */
void MCPathConfigMerton::Generator::configureAntithetics() {
    // All done by MCRandom
    return;
}


/** Configures pathGen for nonAntithetics. We need a better name */
void MCPathConfigMerton::Generator::configureNonAntithetics() {
    // Simulate number of jumps only at nonAntithetic paths

    // maximum # of jumps
    int numSteps = randomJumpTimes.numRows();
    int iStep;
    if(isCarefulRandoms) {
        // Fetch backwards: no real reason to do this given that 
        // there is a random number of random numbers and that we
        // use the same random number generator as the paths.
        for(iStep=numSteps-1; iStep>=0; iStep--){
			rand->fetch(1, &randomJumpTimes[0][iStep]);
        }
    } else {
        // Fetch forward
        rand->fetch(numSteps, &randomJumpTimes[0][0]);
    }

    if( isSimJumpDates ) {
        dependenceJumpTimes->correlateSeries(randomJumpTimes, 0); // dummy pathIdx
    }
    else {
        dependenceJumpTimes->numJumps(randomJumpTimes);
    }
    return;
}


void MCPathConfigMerton::Generator::generatePath(int pathIdx, int iAsset, int iPath){
    const DoubleMatrix& randoms = randomGen->getRandomNumbers();
    const DoubleMatrix& randomJumps = randomGenJumps->getRandomNumbers();
    
    // populate paths field
    const double* assetRandoms  = randoms[iAsset];
    const double* jumpRandoms   = randomJumps[iAsset];
    double*       thisPath      = &productPaths[iAsset][0]; // for ease/speed
    DoubleArray&  thisVol       = vols[iAsset];
    DoubleArray&  thisDrift     = drifts[iAsset];
    for (int iStep = 0, modStep = timeline->numPastDates, jStep=0;
         iStep < timeline->numFutSteps; iStep++, modStep++)
    {
        double sigmaDtDz  = thisVol[iStep]  * assetRandoms[iStep];
        thisPath[modStep] = iStep==0 ? refData->fwdsAtSimStart[iAsset] : thisPath[modStep-1];
        thisPath[modStep] *= thisDrift[iStep];            
        thisPath[modStep] *= exp(sigmaDtDz);
        if( randomJumpTimes[0][iStep]/*numJumps[iStep]*/>0 )
        {
            thisPath[modStep] *= volProcessed[iAsset]->CalcJump(jumpRandoms[jStep],numJumps[iStep]);
            jStep++;
        }
    }
}

void MCPathConfigMerton::Generator::generateDriftAndVols(const DateTimeArray& futurePathDates) {
    static const string routine("MCFuturePathGenMertonBase::"
                                "generateDriftAndVols");
    try{

        int iStep, iAsset, jAsset;

        dt.resize(timeline->numFutSteps);
        for (iStep = 0; iStep < timeline->numFutSteps; iStep ++)
        {
            // futurePathDates.size() = numFutSteps+1, first date is today!
            dt[iStep] = volProcessed[0]->GetTimeMetric()->yearFrac(
                futurePathDates[iStep],
                futurePathDates[timeline->numFutSteps]);
        }

        for(iAsset=0;iAsset<numAssets;iAsset++) {
            correlationJumps[iAsset][iAsset] = 1;
            double beta_i = volProcessed[iAsset]->getBeta();
            for(jAsset=iAsset+1;jAsset<numAssets;jAsset++) {
                double beta_j = volProcessed[jAsset]->getBeta();
                correlationJumps[iAsset][jAsset] = beta_i*beta_j;
                correlationJumps[jAsset][iAsset] = beta_i*beta_j;
            }
        }

        for(iAsset=0;iAsset<numAssets;iAsset++) {

            /* Compute terms to use in path generation.
             *   dS = S.(r-q-k.lambda).dt 
             *           + S.sqrt(sqVol).dW + S.dq
             *      =  'MyDrift'.S 
             *           + 'MyVols'.S.gaussian + S.Y.dq
             */
            // Referenced to the simulated dates (not the simStartDate)
            // drifts[iStep] is drift FROM prev point TO iStep
            // vols[iStep] is sqrt(var) FROM prev point TO iStep
            volProcessed[iAsset]->CalcDrift(    futurePathDates,
                                                VolProcessedMerton::forward,
                                                drifts[iAsset]);
            volProcessed[iAsset]->CalcVar(      futurePathDates,
                                                VolProcessedMerton::forward,
                                                vols[iAsset]);
            for (iStep = 0; iStep < timeline->numFutSteps; iStep ++)
            {
                vols[iAsset][iStep]   = sqrt(vols[iAsset][iStep]);
                drifts[iAsset][iStep] *= fwds[iAsset][iStep+1] / fwds[iAsset][iStep];
            }
            // here we're simulating what happens in generatePath
            // want theDrift to represent [reasonable] worst
            // case drift. A value of 2 for MAX_DRIFT_RAND_NUMBER
            // means that we catch >96% of paths
            //theDrift *= exp(sqrt(sqMeanVol*tau)*MAX_DRIFT_RAND_NUMBER); // ???
            //if (theDrift > maxDrifts[iAsset]){
            //    maxDrifts[iAsset] = theDrift;
            //}
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigMerton::makePathGenerator(
    bool                     cachingRequested,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control, 
    Results*                 results,
    DateTimeArray&           simDates ) {

    static const string routine = "MCPathConfigMerton::makePathGenerator";

    if(cachingRequested) {
        throw ModelException(routine, "Paths caching is not supported in MCPathConfigMerton");
    }

    MCPathBaseSP futurePathGenBaseMerton(
        new MCPathConfigMerton::Generator(numPaths, this, pastPathGenerator, prod));
    return MCPathBase::createPathGenerator(pastPathGenerator,
                                           this,
                                           prod,
                                           futurePathGenBaseMerton);        
}

/* what bounds to use when calculating max drift */
const double MCPathConfigMerton::Generator::MAX_DRIFT_RAND_NUMBER = 3.0;

DRLIB_END_NAMESPACE
