//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigLVMerton.cpp
//
//   Description : Monte Carlo path generator using Merton implied LV
//
//   Date        : 16 Feb 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/LocalVolGrid.hpp"
#include "edginc/VolMertonLVProcessed.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/ProtAsset.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/DependenceGauss.hpp"

DRLIB_BEGIN_NAMESPACE

#define MAX_STEP_SIZE 100		// max time step size 100 yrs
// make sure that max step size is MAX_STEP_SIZE
#define CAP_DT(dt)	((dt)=(dt)>(MAX_STEP_SIZE)?(MAX_STEP_SIZE):(dt))

// #define DEBUG_FILE_NAME 
#ifdef DEBUG_FILE_NAME
FILE * file=fopen("C:\\temp\\MertonLV.txt","w");
#endif

/***************************************** 
MCPathConfigLVMerton class
*****************************************/
/** Class splits its work into three smaller classes. This is possibly more
    code but it makes it much simpler and easier to develop */
class MCPathConfigLVMerton: public MCPathConfig {
private:

    string          volType;
    string          dependenceType;
    bool            isCarefulRandoms;


    double          skewStepScaling;
    int             numVolGrid;
    double          stdevGridRange;

    bool            useTweakingForTimeDerivs;
    double          tweakStrikeUnscaled;
    double          tweakTimeUnscaled;
    double          probDensRatioMin;
    bool            useMidPoint;

    // transient fields
    DateTimeArraySP               simDatesLV; // cached between tweaks
    // for reflection
    MCPathConfigLVMerton();

    class Generator;
    friend class Generator;
public:
    static CClassConstSP const TYPE;


    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const{
        return MarketDataFetcherSP(new MarketDataFetcherLNSpline(volType));
    }

    virtual int randomStoragePerPath(IMCProduct* product) const;

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

    virtual bool carefulRandoms() const { return isCarefulRandoms; }

    /** Throws away cached sim dates for Theta-type tweaks */
    virtual bool sensShift(Theta* shift) {
        MCPathConfig::sensShift(shift); // call parents method
        simDatesLV.reset();
        return true;
    }
protected:
	MCPathConfigLVMerton(CClassConstSP clazz, const string& volType, const string& dependenceType):
         MCPathConfig(clazz), volType(volType), isCarefulRandoms(false){}

    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,    
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates );

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigLVMerton, clazz);
        SUPERCLASS(MCPathConfig);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        EMPTY_SHELL_METHOD(defaultMCPathConfigLVMerton);

        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);  // for now 
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);

        FIELD(skewStepScaling, "scaling factor for num of steps, higher number give more steps. -1=one step per trading day");
        FIELD_MAKE_OPTIONAL(skewStepScaling);        
        FIELD(numVolGrid, "num of pts for pre-computed vol grid (default=41)");
        FIELD_MAKE_OPTIONAL(numVolGrid);        
        FIELD(stdevGridRange, "stdev of vol grid range (default = 6.0)");
        FIELD_MAKE_OPTIONAL(stdevGridRange);        

        FIELD(useTweakingForTimeDerivs, "useTweakingForTimeDerivs");
        FIELD_MAKE_OPTIONAL(useTweakingForTimeDerivs);
        FIELD(probDensRatioMin, "probDensRatioMin");
        FIELD_MAKE_OPTIONAL(probDensRatioMin);
        FIELD(tweakStrikeUnscaled, "tweakStrikeUnscaled");
        FIELD_MAKE_OPTIONAL(tweakStrikeUnscaled);
        FIELD(tweakTimeUnscaled, "tweakTimeUnscaled");
        FIELD_MAKE_OPTIONAL(tweakTimeUnscaled);
        FIELD(useMidPoint, "useMidPoint");
        FIELD_MAKE_OPTIONAL(useMidPoint);        

        FIELD_NO_DESC(simDatesLV);
        FIELD_MAKE_TRANSIENT(simDatesLV); 

        FIELD(isCarefulRandoms, "isCarefulRandoms");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms);
        clazz->setPublic(); // make visible to EAS/spreadsheet

	}

    static IObject* defaultMCPathConfigLVMerton(){
        return new MCPathConfigLVMerton();
    }
};

CClassConstSP const MCPathConfigLVMerton::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigLVMerton", typeid(MCPathConfigLVMerton), load);

///////////////////////////////////////////////////////////////////////////////////
// Generator
///////////////////////////////////////////////////////////////////////////////////

typedef smartConstPtr<Poisson> PoissonConstSP;
typedef smartPtr<Poisson> PoissonSP;
typedef smartConstPtr<VolMertonLVProcessed> VolMertonLVProcessedConstSP;
typedef smartPtr<VolMertonLVProcessed> VolMertonLVProcessedSP;

/** Class containing the basic variable and methods that (hopefully) will 
be used by a lot of path generators. */
class MCPathConfigLVMerton::Generator: virtual public MCPathBase,
                                       virtual public IMCRandom::Callbacks,
                                       public DependenceMakerGauss::Support {

public:
    Generator(int                      numSimPaths,
              MCPathConfigLVMerton*    mcPathConfigLVMerton, 
              const MCPathGeneratorSP& pastPathGenerator,
              const IMCProduct*         prod);
   
	// future steps for simulation, this is futureDates and LV dates combined
    DateTimeArraySP simDatesLV; 
    // at each step, true=is a futureDate, fasle=not a future sample date
    BoolArray isSampleDate; 
    BoolArray isQuanto;

     //vector<VolMertonLVProcessedSP> locVols;    // [numAssets] 

    // local vol grid set up at the begining of simulation for spline interpolation
    // each asset has a lvGrid which contain lv grid for every path step
    vector<ILocalVolGridSP>     locVolGrid;

    DoubleArray            maxDrifts; // product of MAX(1, drifts)

	DoubleArray		JumpCorrection;


//    FXAssetConstSP fxAsset;
    double              eqFXCorr;
    vector<DoubleArray> fxVar;

    void generateDriftAndVols(const DateTimeArray& futurePathDates);
       
    vector<VolMertonLVProcessedSP> volProcessed;

    /** Draws random numbers */
    virtual void drawRandomNumbers(int pathIdx);
    
    virtual void generatePath(int pathIdx, int iAsset, int iPath);

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
    int NbSimAssets() const;

    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath) ;
    
    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const;

    /** Returns the paths */
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

    /** Configures pathGen for nonAntithetics */
    virtual void configureNonAntithetics();
    
    /** generate simulation dates. First using vol skew to generate steps and then merge with futureDates */
    void genSimDates(const DateTime& startDate, double skewStepScaling);

    /** pre-compute fwds arrays, get quanto data here for now*/
    void preCompute(int numVolGrid, double stdevGridRange, vector<DateTimeArraySP> cumDates);

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const;

protected:
    // fields
    // Merton
    int                       maxNumJumps;
    double                    quantile;       // probability of having more than numJumps jumps
    IntArray                  numJumps;
    DoubleMatrix              randomJumpTimes;// 1 row, maxNumJumps cols
    DoubleArray               dt;
    DoubleMatrix              correlationJumps;
    DependenceSP              dependenceJumps;
    PoissonSP                 dependenceJumpTimes;
	long idum;

    // transient
    mutable TimeMetricSP      timeMetric;

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
};

///////////////////////////////////////////////////////////////////////////////////
// Generator, MCSimpleFuturePathGeneratorMerton, MCFuturePathGeneratorMerton
///////////////////////////////////////////////////////////////////////////////////
/***************************************** 
Generator
*****************************************/
MCPathConfigLVMerton::Generator::Generator(int                      numSimPaths,
                                           MCPathConfigLVMerton*    mcPathConfigLVMerton, 
                                           const MCPathGeneratorSP& pastPathGenerator,
                                           const IMCProduct*         prod):
    simDatesLV(mcPathConfigLVMerton->simDatesLV), 
    maxDrifts(prod->getNumAssets(), 1.0),
    fwds(prod->getNumAssets()),
    mAsset(prod->getMultiFactors()), 
    numAssets(prod->getNumAssets()),
    productPaths(prod->getNumAssets()),
    rand(mcPathConfigLVMerton->getRandomGenerator()),
    isCarefulRandoms(mcPathConfigLVMerton->carefulRandoms()) {
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

        // sizes
        locVolGrid.resize(numAssets);
        volProcessed.resize(numAssets);
        fxVar.resize(numAssets);
		JumpCorrection.resize(numAssets);

		// obtain and sort cum div dates. no local vol averaging across div dates
		vector<DateTimeArraySP> cumDates(numAssets);
        EventAssetMove divEvent;
		const DateTime& startDate = prod->getEffectiveSimStartDate();
		const DateTime& mat = timeline->futureDates[timeline->futureDates.size()-1];
        int numDivs = (int)(4*startDate.yearFrac(mat))+1; // 4 divs per year selected as critical dates
		for (int i=0; i<numAssets; i++)
        {// each asset selects largestN per year div dates to go into critical dates
            if (AssetUtil::getJumpEvents(&mAsset->getAsset(i), 
                                         startDate, 
                                         mat,
                                         numDivs,
                                         divEvent))
            {	
                cumDates[i] = divEvent.getCritDate(0, true);
				sort(cumDates[i]->begin(), cumDates[i]->end());
			}
			else
				cumDates[i]  = DateTimeArraySP(new DateTimeArray()); // create empty array in case of no div dates
        }
        // init vols
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            // local vol request for each asset
            LocVolRequest volRequest(timeline->today, // always start simulate today
                                     prod->getEffectiveSimStartDate() > timeline->today,
                                     false,
                                     mcPathConfigLVMerton->useTweakingForTimeDerivs,
                                     mcPathConfigLVMerton->useMidPoint,
                                     mcPathConfigLVMerton->tweakStrikeUnscaled,
                                     mcPathConfigLVMerton->tweakTimeUnscaled,
                                     mcPathConfigLVMerton->probDensRatioMin);

            volProcessed[iAsset] = VolMertonLVProcessedSP(dynamic_cast<VolMertonLVProcessed*>(
                                                     mAsset->factorGetProcessedVol(iAsset, &volRequest))); 

			JumpCorrection[iAsset] =volProcessed[iAsset]->getMertonParam("JumpRate")
						*(exp(volProcessed[iAsset]->getMertonParam("JumpMean")+
						0.5*volProcessed[iAsset]->getMertonParam("JumpWidth")
						*volProcessed[iAsset]->getMertonParam("JumpWidth"))-1.0);
        }
        // reuse sim dates between tweaks for numerical stability
        if (simDatesLV->empty()){
            /** generate simulation dates. First using vol skew to
                generate steps and then merge with futureDates */
            genSimDates(prod->getEffectiveSimStartDate(),
                        mcPathConfigLVMerton->skewStepScaling);
        }
        
        // set future dates flag
        isSampleDate.resize(simDatesLV->size(), false);
        isQuanto.resize(numAssets, false);

        int iSearch = 0;
        for (int j=0; j<simDatesLV->size(); j++)
        {
            if (timeline->futureDates[iSearch] == (*simDatesLV)[j])
            {
                isSampleDate[j] = true;
                iSearch++;
                if (iSearch == timeline->futureDates.size())
                    break;
            }
        }

        // pre compute fwds, vol grid and quanto stuff
        preCompute(mcPathConfigLVMerton->numVolGrid, mcPathConfigLVMerton->stdevGridRange, cumDates);

        // Merton
        // 99.99% quantile of number of jumps
        double quantile;
        volProcessed[0]->Quantile( 
            timeline->today,(*simDatesLV)[simDatesLV->size()-1], 0.00001, 1000,
            &quantile,&maxNumJumps );

        if( maxNumJumps>0 ) {
            randomJumpTimes = DoubleMatrix(1,maxNumJumps);
        }
        numJumps = IntArray(simDatesLV->size()-1);
            
        /* now we know how many paths per asset we can allocate
           space for them. Also note that we can't ask for
           the vol interps until the past has been computed */
        correlationJumps = DoubleMatrix(numAssets, numAssets);

        // precompute drift and vols could possibly have a class
        // "process" which did this - for now leave as local functions
        generateDriftAndVols(*simDatesLV);

		// set up Dependence
        dependence = mcPathConfigLVMerton->dependenceMaker->createDependence(this);
        dependenceJumps     = DependenceSP(new Gauss(correlationJumps));
        dependenceJumpTimes = PoissonSP(new Poisson(volProcessed[0]-> getMertonParam("JumpRate")));

        // Initialize random number generator for paths. No callback to pathGen
        randomGen = IMCRandomSP(new MCRandomNoCache(
            0,
            dependence,
            rand,
            mcPathConfigLVMerton->carefulRandoms(),
            simDatesLV->size()-1,
            numAssets,
            timeline->totalNumPastDates));

        // Initialize random number generator for jumps. Callback to pathGen
        randomGenJumps = IMCRandomSP(new MCRandomNoCache(
            this,
            dependenceJumps,
            rand,
            mcPathConfigLVMerton->carefulRandoms(),
            Maths::min(maxNumJumps,simDatesLV->size()-1),
            numAssets,
            0));

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


/** Returns number of assets */
int MCPathConfigLVMerton::Generator::NbSimAssets() const {
    return numAssets;
}

/** Returns the reference level for iAsset, iPath */
double& MCPathConfigLVMerton::Generator::refLevel(int iAsset, int iPath) {
    return refData->refLevels[iAsset][0];
}
    
/** Returns the reflevel path */
IRefLevel::IMCPathSP& MCPathConfigLVMerton::Generator::refLevelPath() const {
    return refData->refLevelPath;
}

/** Returns the paths */
double* MCPathConfigLVMerton::Generator::Path(int iAsset, int iPath) {
    return &productPaths[iAsset][0];
}


/** simulates random numbers */
void MCPathConfigLVMerton::Generator::drawRandomNumbers(int pathIdx) {
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


void MCPathConfigLVMerton::Generator::generatePath(int pathIdx, int iAsset, int iPath){
    const DoubleMatrix& randoms = randomGen->getRandomNumbers();
    const DoubleMatrix& randomJumps = randomGenJumps->getRandomNumbers();
    
    int numPathSteps = simDatesLV->size();
    // populate paths field
    const double* assetRandoms  = randoms[iAsset];
    const double* jumpRandoms   = randomJumps[iAsset];
    double*       thisPath      = &productPaths[iAsset][0]; // for ease/speed
    const DoubleArray& assetFwds = fwds[iAsset];
    ILocalVolGrid& pathVols = *locVolGrid[iAsset];
	double logS=log(refData->fwdsAtSimStart[iAsset]/assetFwds[0]);
	double vol_dt;
    for (int iStep = 0, modStep = timeline->numPastDates, jStep=0;
         iStep < numPathSteps-1; iStep++)
    {
        vol_dt = pathVols.interpLocVar(iStep, logS);

		// add quanto adjustment
        if (isQuanto[iAsset])
			logS += -eqFXCorr * fxVar[iAsset][iStep] * vol_dt;

		logS  += vol_dt * (assetRandoms[iStep]-0.5*vol_dt)
			-(dt[iStep]-dt[iStep+1])*JumpCorrection[iAsset];
        if( numJumps[iStep]>0 )
        {
            logS += volProcessed[iAsset]->CalcJump(jumpRandoms[jStep],numJumps[iStep]);
            jStep++;
        }
        if(isSampleDate[iStep+1]){
			thisPath[modStep] = assetFwds[iStep+1] * exp(logS);
            modStep++;
            if (modStep > timeline->numFutSteps + timeline->numPastDates) // can't be
                throw ModelException("MCPathBaseLVMerton::generatePath",
                                     "internal error");
		}
    }
}

///////////////////bad name for the function, but it is the one used by Olivier Brockhaus.///used differently by him, of course.
void MCPathConfigLVMerton::Generator::generateDriftAndVols(const DateTimeArray& futurePathDates) {
    static const string routine("MCFuturePathGenMertonBase::"
                                "generateDriftAndVols");
    try{

        int iStep, iAsset, jAsset;

        dt.resize(simDatesLV->size());
        for (iStep = 0; iStep < simDatesLV->size(); iStep ++)
        {
            // futurePathDates.size() = numFutSteps+1, first date is today!
            dt[iStep] = volProcessed[0]->GetTimeMetric()->yearFrac(
                (*simDatesLV)[iStep],
                (*simDatesLV)[simDatesLV->size()-1]);
        }

        for(iAsset=0;iAsset<numAssets;iAsset++) {
            correlationJumps[iAsset][iAsset] = 1;
            for(jAsset=iAsset+1;jAsset<numAssets;jAsset++) {
                correlationJumps[iAsset][jAsset] = 0.0;////can add correlation as Olivier Brockhaus did in his Merton Model
                correlationJumps[jAsset][iAsset] = 0.0;
            }
        }

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


/** Configures pathGen for antithetics */
void MCPathConfigLVMerton::Generator::configureAntithetics() {
    // All done by the MCRandom class
    return;
}


/** Configures pathGen for nonAntithetics */
void MCPathConfigLVMerton::Generator::configureNonAntithetics() {
    //Draw random jump times when not doing antithetics
    
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
	dependenceJumpTimes->correlateSeries(randomJumpTimes, 0); // dummy pathIdx

	/////////////////////////////different from Oliver Brockhaus//////////////
	///put the randomJumpTimes in the corresponding intervals of simDatesLV////
    int i, j;
	for(i=0;i<numJumps.size();i++) numJumps[i]=0;
	for(i=0;i<randomJumpTimes.numRows();i++) {
		for (j=0;j<simDatesLV->size()-1;j++){
			if ((dt[j]>=randomJumpTimes[0][i])&&(dt[j+1]<randomJumpTimes[0][i]))
				numJumps[j]+=1;
		}
	}
    return;
}


/////////////////////////////
/** pre-compute fwds arrays, get quanto data here for now */
void MCPathConfigLVMerton::Generator::preCompute(int numVolGrid, double stdevGridRange, vector<DateTimeArraySP> cumDates)
{
	idum=-1;
    double sqr_dt;
    int i, j;

    // for each asset
    for (i=0; i<numAssets; i++)
    {
        // just for getting dt for each asset
        TimeMetricConstSP metric = volProcessed[i]->GetTimeMetric();
        // need to cache the forward at all simulation dates
        fwds[i] = DoubleArray(simDatesLV->size());

        // current way of dealing with quanto
        CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(i));

        if (ProtAsset::TYPE->isInstance(asset)) {
            throw ModelException("MCPathConfigLVMerton::Generator::preCompute", 
                                 "Protected assets (" + asset->getTrueName() + 
                                 ") of type " + asset->getClass()->getName() + 
                                 " are not supported yet");
        }    

        isQuanto[i] = ProtEquity::TYPE->isInstance(asset);
        if (isQuanto[i])
        {
            fxVar[i] = DoubleArray(simDatesLV->size()-1);
            ProtEquityConstSP eq = ProtEquityConstSP::dynamicCast(asset);
            // compute plain fwds
            CAssetConstSP plainAsset = eq->getPlainAsset();
            plainAsset->fwdValue(*simDatesLV, fwds[i]);

            // get fx vol and corr
            FXVolBaseWrapper fxVol = eq->getFXVol();
            eqFXCorr = eq->getCorrelation()->getCorrelation();;
            // interpolate the fx vols atm 
            ATMVolRequestSP fxVolRequest(new ATMVolRequest());
            // interpolate the vol
            const Asset* fx = eq->getFXAsset();
            CVolProcessedSP  fxVoltmp(fxVol->getProcessedVol(fxVolRequest.get(), fx));
            // cast to the type of vol we're expecting
            CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(fxVoltmp);
            volFXBS->CalcVol(*simDatesLV, CVolProcessedBS::forward, fxVar[i]);

            for (j=0; j<fxVar[i].size(); j++)
            {
				sqr_dt = sqrt(metric->yearFrac((*simDatesLV)[j], (*simDatesLV)[j+1]));
                fxVar[i][j] *= sqr_dt; // using asset dt
            }
		}
        else{
            mAsset->factorFwdValues(i, *simDatesLV, fwds[i]);
        }

        // generate local vol grid for all sample points
        
        locVolGrid[i] = LocalVolGrid::createLocalVolGrid(fwds[i], volProcessed[i], *simDatesLV, 
            *cumDates[i], numVolGrid, stdevGridRange);
		maxDrifts[i] = locVolGrid[i]->maxDrift();
    }
}

// implement getGaussData since derivation from DependenceMakerGauss::Support
CDoubleMatrixConstSP MCPathConfigLVMerton::Generator::getGaussData() const {
        static const string method("MCPathConfigLVMerton::Generator::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

/** generate simulation dates. First using vol skew to generate steps and then merge with futureDates */
void MCPathConfigLVMerton::Generator::genSimDates(const DateTime& startDate, double skewStepScaling)
{
    static const string method("MCPathBaseLVMerton::genSimDates");

    try{
        const DateTime& mat = timeline->futureDates[timeline->futureDates.size()-1];
        const double one_week = 1.0/52.0;
        const double adjFactor = 2.0; // term structure scaling adjustment, an empirical factor here

        CSliceDouble k(2), var(2);
        double tmp, maxVolDiff = 0.0;
        DoubleArray spots(numAssets);
        DateTimeArray t(2), t_next(2);

        // just use the one time metric
        TimeMetricConstSP metric = volProcessed[0]->GetTimeMetric();

        // use 90-100 skew (vol*sqrt(deltat)) and term structure to determine step size
        // compute starting maxVolDiff
        t[0] = metric->impliedTime(startDate, 1.0/365, tmp); // QuadExp does not give skew  for start=today
        t[1] = metric->impliedTime(t[0], one_week, tmp); // tmp not used, one week loc vol*sqrt(deltat)
        for (int iAsset=0; iAsset<numAssets; iAsset++)
        {
            // compute 90-100 fwd local vol spread
            spots[iAsset] = mAsset->factorGetSpot(iAsset);
            k[0] = 0.9*spots[iAsset];
            k[1] = spots[iAsset];

            volProcessed[iAsset]->CalcLocVar(&k, t, &var, false); // no intradayinterp
            tmp = fabs(sqrt(var[0]) - sqrt(var[1]));
            if (tmp > maxVolDiff)
                maxVolDiff = tmp;

            // also take into account of term structure
            var[0] = volProcessed[iAsset]->computeImpVol(t[1], k[1]);
            var[1] = volProcessed[iAsset]->computeImpVol(mat, k[1]);
            tmp = adjFactor*fabs(var[1] - var[0])*sqrt(one_week);
            if (tmp > maxVolDiff)
                maxVolDiff = tmp;
        }

        // now space time axis according to maxVolDiff
        simDatesLV->clear();
        if (skewStepScaling == -1.0) // daily steps requested
        {
            double deltat = 1.0/metric->volDays(startDate, startDate.rollDate(365));
            DateTime nextDate(metric->impliedTime(startDate, deltat, tmp)); // tmp not used
            while (nextDate < mat)
            {
                simDatesLV->push_back(nextDate);
                nextDate = metric->impliedTime(nextDate, deltat, tmp);
            }
        }
        else if (maxVolDiff>0.0 && skewStepScaling > 0.0) // if starting skew is 0 then we assume zero skew
        {
			CSliceDouble var_next(2);
            double deltat = 1.0/100.0/maxVolDiff/skewStepScaling; // const does not matter
			CAP_DT(deltat); // cap deltat at MAX_STEP_SIZE. 

			DateTime nextDate(metric->impliedTime(startDate, deltat, tmp)); // tmp not used
            while (nextDate < mat)
            {
                simDatesLV->push_back(nextDate);

                // use 90-100 skew (vol*sqrt(deltat)) and term structure to determine step size
                // compute starting maxVolDiff
                t[0] = nextDate;
                t[1] = metric->impliedTime(t[0], deltat, tmp); // tmp not used
                t_next[0] = t[1];
                t_next[1] = metric->impliedTime(t_next[0], deltat, tmp);
                maxVolDiff = 1.0e-20; // set to almost 0
                for (int i=0; i<numAssets; i++)
                {
                    // compute 90-100 fwd local vol spread
                    k[0] = 0.9*spots[i];
                    k[1] = spots[i];

                    volProcessed[i]->CalcLocVar(&k, t, &var, false); // no intradayinterp

                    tmp = fabs(sqrt(var[0]) - sqrt(var[1]))*sqrt(one_week/deltat);
                    if (tmp > maxVolDiff)
                        maxVolDiff = tmp;

                    // also take into account of term structure
					// this makes the algorithm more stable. otherwise, vol surface noise may give tmp~0 when it shouldn't be 
                    volProcessed[i]->CalcLocVar(&k, t_next, &var_next, false); // no intradayinterp
                    tmp = fabs(sqrt(var_next[0]) - sqrt(var_next[1]))*sqrt(one_week/deltat);
                    if (tmp > maxVolDiff)
                        maxVolDiff = tmp;
                    tmp = adjFactor*fabs(sqrt(var_next[1]) - sqrt(var[1]))*sqrt(one_week/deltat);
                    if (tmp > maxVolDiff)
                        maxVolDiff = tmp;
                    tmp = adjFactor*fabs(sqrt(var_next[0]) - sqrt(var[0]))*sqrt(one_week/deltat);
                    if (tmp > maxVolDiff)
                        maxVolDiff = tmp;
                }

                deltat = 1.0/100.0/maxVolDiff/skewStepScaling;
				CAP_DT(deltat); // cap deltat at MAX_STEP_SIZE. 
                nextDate = metric->impliedTime(nextDate, deltat, tmp);
            }
        }

        // now, add futureDates and sort
        MCProductTimelineSP timelineTemp(MCProductTimelineSP::constCast(timeline));
        simDatesLV->insert(simDatesLV->end(), 
                           timelineTemp->futureDates.begin(), 
                           timelineTemp->futureDates.end());

        sort(simDatesLV->begin(), simDatesLV->end());
        // remove duplicate dates
        for (vector<DateTime>::iterator iter = simDatesLV->begin()+1; iter != simDatesLV->end(); iter++)
        {
            if (*iter == *(iter-1))
            {
                simDatesLV->erase(iter);
                iter --;
            }
        }

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

////////////////////////////////////////////////
// MCPathConfigLVMerton class member functions
////////////////////////////////////////////////
MCPathConfigLVMerton::MCPathConfigLVMerton():MCPathConfig(TYPE),
                                            volType("VolMertonLV"),
                                            dependenceType("not used"),
                                            isCarefulRandoms(false),
                                            skewStepScaling(1.0),
                                            numVolGrid(101),
                                            stdevGridRange(6.0),
                                            useTweakingForTimeDerivs(true),
                                            tweakStrikeUnscaled(0.01),
                                            tweakTimeUnscaled(0.001),
                                            probDensRatioMin(0.01),
                                            useMidPoint(false){}

/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigLVMerton::makePathGenerator(
    bool                               cachingRequested,
    int                                numPaths,
    const MCPathGeneratorSP&           pastPathGenerator,
    const IMCProduct*                   prod,
    Control*                           control, 
    Results*                           results,
    DateTimeArray&                     simDates ){
    static const string routine = "MCPathConfigLVMerton::makePathGenerator";

    if(cachingRequested) {
        throw ModelException(routine, "Paths caching is not supported in MCPathConfigLVMerton");
    }

    // create empty sim dates array if simDatesLV is null or on new pricing run
    if (!simDatesLV || control->isPricing()){
        simDatesLV = DateTimeArraySP(new DateTimeArray());
    }
    refCountPtr<MCPathConfigLVMerton::Generator> futurePathGenBaseMertonLV(
        new MCPathConfigLVMerton::Generator(numPaths,
                                            this,
                                            pastPathGenerator,
                                            prod));
    
    // store num of sim dates
    if (control && control->isPricing()){
        int simSteps = futurePathGenBaseMertonLV->simDatesLV->size();
        results->storeScalarGreek(simSteps-1, Results::DEBUG_PACKET, 
                                  OutputNameSP(
                                      new OutputName("MertonLV_SIM_STEPS_USED")));
    }

    // Construct the MCPathGenerator
    return MCPathBase::createPathGenerator(pastPathGenerator,
                                           this,
                                           prod,
                                           futurePathGenBaseMertonLV);

}

/** Returns the number of bytes used for random number storage per
    path. Do not invoke if there are no sim dates in future */
int MCPathConfigLVMerton::randomStoragePerPath(IMCProduct* product) const{
    // Otherwise have to worry about how many random numbers we use
    // This is a bit of a pain as currently there is no easy way to get hold
    // of the number of dates unless we build the entire path generator
    smartPtr<MCPathConfigLVMerton>
        pathConfig(copy(this)); // copy to avoid const problems
    MCPathGeneratorSP pastPathGenerator(pathConfig->pastPathGenerator(product));
    SensitivityArrayConstSP    sens;
    OutputRequestArrayConstSP  request;
    Control control(sens, request, false, "");
    Results results;
    DateTimeArray simDates;
    MCPathGeneratorSP pathGen(pathConfig->futurePathGenerator(0, // no caching
                                                              2, // num paths 
                                                              pastPathGenerator,
                                                              product,
                                                              &control,
                                                              &results,
                                                              simDates ));
    if (!pathConfig->simDatesLV){
        throw ModelException("MCPathConfigLVMerton::storagePerPath", "Internal "
                             "error - no simulation dates");
    }
    // we store both correlated and uncorrelated numbers but only for
    // every other path. Hence no times by 2.
    return (sizeof(double) * pathConfig->simDatesLV->size() * 
            product->getNumAssets());
}

/** MC LV model */
class MCPathConfigLVMertonDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MCPathConfigLVMertonDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigLVMertonDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloLV);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloLV(){
        return new MCPathConfigLVMertonDefault();
    }
};

CClassConstSP const MCPathConfigLVMertonDefault::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigLVMertonDefault", typeid(MCPathConfigLVMertonDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigLVMertonLoad(){
    return (MCPathConfigLVMerton::TYPE != 0 && MCPathConfigLVMertonDefault::TYPE !=0);
}

DRLIB_END_NAMESPACE
