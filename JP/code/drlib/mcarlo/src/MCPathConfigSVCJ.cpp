//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSVCJ.cpp
//
//   Description : Monte Carlo path generator for VolSVCJ
//
//   Date        : 20 Oct 2004
//
//
//   $Log: MCPathConfigSVCJ.cpp,v $
//   Revision 1.13  2005/03/10 23:02:05  rgaragnon
//   changed state variable part to jump arrival simulation
//
//   Revision 1.12  2005/02/03 19:20:21  nshen
//   name change for state variable.
//
//   Revision 1.11  2005/01/11 16:37:35  nshen
//   meanVar dealt with in VolSVCJ now.
//
//   Revision 1.10  2005/01/06 17:14:57  nshen
//   added meanVar time dependency; fixed mem leak.
//
//   Revision 1.9  2004/12/14 17:20:12  nshen
//   checked in the right file.
//
//   Revision 1.8  2004/12/09 22:58:19  nshen
//   corrected RandomGamma (and just directly use imsl). InitVol simulation works correctly.
//
//   Revision 1.7  2004/11/30 21:49:41  rgaragnon
//   corrected a bug in CDFMapping object construction
//
//   Revision 1.6  2004/11/30 18:04:23  nshen
//   changed getMarket().
//
//   Revision 1.5  2004/11/30 16:05:26  rguichard
//   Added convert method to RandomGamma
//
//   Revision 1.4  2004/11/24 17:02:15  nshen
//   needs to use CDFMapping asset for market pdf.
//
//   Revision 1.3  2004/11/23 22:09:06  nshen
//   added CDFMapping support; started generating initialVol for stationary variance distribution.
//
//   Revision 1.2  2004/11/11 16:46:41  nshen
//   pass iterator by value to fix unix build.
//
//   Revision 1.1  2004/11/10 21:57:48  nshen
//   simulate jump arrivals.
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/ProtAsset.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MCPathConfigParametric.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/VolSVCJ.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/DependenceGauss.hpp"

// for Implied Mapping
#include "edginc/CDFMapping.hpp"
#include "edginc/PDFCalculatorMaker.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"

#define STATE_VARIABLES

#ifdef STATE_VARIABLES
#include "edginc/mathlib.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#endif

DRLIB_BEGIN_NAMESPACE

/////////////////////////////////////////////////////////////////////////
//-----------------  Gamma distribution generator ----------------
// a simple wrap of imsl gamma generator, left here until we have a generic one
/////////////////////////////////////////////////////////////////////////
class MCRandGamma{
public:
    MCRandGamma(double alpha, int seed):
    alpha(alpha),seed(seed){
        imsl_random_seed_set(seed);
    }

    void fetch(int numToFetch, double* rands) {
        imsl_d_random_gamma(numToFetch, alpha, IMSL_RETURN_USER, rands);
    }

private:
    double  alpha;
    int     seed;
};
typedef refCountPtr<MCRandGamma> MCRandGammaSP;

// Structure used to store parameters for zbrent root finder
typedef struct _JumpCorrelationData
{
    double lowJumpProba; // lower jump marginal probability
    double lowJumpIntensity; // lower jump marginal intensity
    double highJumpProba; // higher jump marginal probability
    double highJumpIntensity; // higher jump marginal intensity
    double normInvLowProba; // lowJumpProba quantile of N1
    double normInvHighProba; // highJumpProba quantile of N1
    double correlation; // correlation between the corresponding asset
    double targetProba; // target conditional jump probability participation

} JumpCorrelationData;

// Computes the conditional jump participation probability for two assets
double objectivFunc(double jumpCorrelation, void *objectivFuncDataMem)
{
    JumpCorrelationData* JumpCorrelationStruct = (JumpCorrelationData*)objectivFuncDataMem;

    double lowJumpProba = JumpCorrelationStruct->lowJumpProba; // lower jump marginal probability
    double normInvLowProba = JumpCorrelationStruct->normInvLowProba; // lowJumpProba quantile of N1
    double normInvHighProba = JumpCorrelationStruct->normInvHighProba; // highJumpProba quantile of N1
    double targetProba = JumpCorrelationStruct->targetProba; // target conditional jump probability participation

    // we want to solve func = 0 in the variable jumpCorrelation
    double func = ((1 / lowJumpProba) * N2(normInvLowProba, normInvHighProba, jumpCorrelation)) - targetProba;
    return func;
}

///////////////////////////////////////////////////////////////////////////////////
// MCPathConfig
///////////////////////////////////////////////////////////////////////////////////
class MCPathConfigSVCJ: public MCPathConfigParametric {
public:
    static CClassConstSP const TYPE;

    /* create MDF class so we can pass info to CDFMapping - may want to
       put this class on CDFMapping itself - unclear at present */
    class MDF: public MDFAssetVol{
        CDFMappingSP cdfMapping; // reference

    public:
        static CClassConstSP const TYPE;
        mutable bool insideMDFFetchGetData; /* Need to remember if we are inside the
                                             * fetch, because otherwise there could be
                                             * an infinite recursion when calling
                                             * market->GetData */

        MDF(const string& volType, CDFMappingSP cdfMapping):
            MDFAssetVol(volType, VolSurface::TYPE->getName()),
            cdfMapping(cdfMapping), insideMDFFetchGetData(false) {}

        virtual MarketObjectSP fetch(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type,
                                     const IModel*        model) const {
            /* if it is an asset, instruct the MDF to behave like a
             * MarketDataFetcherLN, ie, set collapseToVolSurface to true -
             * and then grab the asset and set it in the CDFMaping */
            if (!insideMDFFetchGetData &&
                cdfMapping.get() &&
                CAsset::TYPE->isAssignableFrom(type))
            {
                MDF* ncThis = const_cast<MDF*>(this); // dodgy...

                // Set the MDF to behave like a MDFLN, ie, set the
                // collapseToVolSurface to true and the volType to
                // "VolPreferred"
                bool originalCollapseToVolSurface =
                    ncThis->setCollapseToVolSurface(true);
                string originalVolType = ncThis->setVolType("VolPreferred");

                insideMDFFetchGetData = true;
                MarketObjectSP asset(market->GetData(model, name, type));
                insideMDFFetchGetData = false;

                cdfMapping->setAsset(asset);

                // Restore the original values
                ncThis->setCollapseToVolSurface(originalCollapseToVolSurface);
                ncThis->setVolType(originalVolType);
            }
            return MDFAssetVol::fetch(market, name, type, model);
        }
    };

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retreiving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const {
        return MarketDataFetcherSP(new MDF(VolSVCJ::TYPE->getName(),
                                           cdfMappingParam));
    }

protected:
    // for reflection
    MCPathConfigSVCJ(CClassConstSP clazz):
    MCPathConfigParametric(clazz),
    isCarefulRandoms(false),
    dependenceType("not used"),
    nbStepsPerYear(12),
    spotDiscreteScheme(VolSVJ_SpotDSTypeHelper::getDefaultName()),
    varDiscreteScheme(VolSVJ_VarDSTypeHelper::getDefaultName()),
    spotDiscreteSchemeIndex(-1),
    varDiscreteSchemeIndex(-1){}

    virtual void validatePop2Object(){
        static const string routine("MCPathConfigSVCJ::validatePop2Object");
        try {
            Maths::checkPositive(nbStepsPerYear, "nbStepsPerYear");
            spotDiscreteSchemeIndex = VolSVJ_SpotDSTypeHelper::getIndex(spotDiscreteScheme);
            varDiscreteSchemeIndex = VolSVJ_VarDSTypeHelper::getIndex(varDiscreteScheme);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Returns the number of bytes used for random number storage per
        path. Do not invoke if there are no sim dates in future */
    virtual int randomStoragePerPath(IMCProduct* product) const{
        // Otherwise have to worry about how many random numbers we use
        // This is a bit of a pain as currently there is no easy way to get hold
        // of the number of dates unless we build the entire path generator
        smartPtr<MCPathConfigSVCJ> pathConfig(copy(this)); // copy to avoid const problems
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
        if (!pathConfig->simDates){
            throw ModelException("MCPathConfigSVCJ::storagePerPath", "Internal "
                                 "error - no simulation dates");
        }
        // we store both correlated and uncorrelated numbers but only for
        // every other path. Hence no times by 2.
        return (sizeof(double) * pathConfig->simDates->size() *
                product->getNumAssets());
    }

    virtual bool vegaMatrixSupported() const {
        return false;
    }

    virtual bool carefulRandoms() const {
        return isCarefulRandoms;
    }

    /** Throws away cached sim dates for Theta-type tweaks */
    virtual bool sensShift(Theta* shift) {
        MCPathConfig::sensShift(shift); // call parents method
        simDates.reset();
        return true;
    }

    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control,
        Results*                           results,
        DateTimeArray&                     simDates );

    /** Creates a past path generator */
    MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod);

    /** Helper function used by both Generator and PathGenFineProcess.
        Creates a time line that is uniformally spaced between each
        fixing date aka 'futureDate' */
    static DateTimeArraySP genSimDates(const DateTimeArray& futureDates,
                                       const TimeMetric&    timeMetric,
                                       int                  nbStepsPerYear,
                                       IntArray&            fixDateOffsets) {
        static const string method("MCPathConfigSVCJ::Generator::genSimDates");
        try{
            // calculate approx number of sim dates
            int nbFutureDates = futureDates.size();
            ASSERT(nbFutureDates > 1);
            DateTime startDate = futureDates.front();
            DateTime horizonDate = futureDates.back();
            double horizonYearFrac = timeMetric.yearFrac(startDate, horizonDate);
            int approxNbSteps = (int)(ceil(horizonYearFrac * nbStepsPerYear));
            DateTimeArraySP simDates(new DateTimeArray());
            simDates->reserve(approxNbSteps + 1);
            // loop over future dates and insert additional points as needed
            // calculate fixing date offsets simultaneously
            fixDateOffsets.resize(nbFutureDates - 1);
            DateTime currFutureDate = startDate;
            simDates->push_back(currFutureDate);
            for (int iDate = 1, iOffset = 0; iDate < nbFutureDates; iOffset = iDate++){
                // calculate mesh size
                DateTime nextFutureDate = futureDates[iDate];
                double yearFrac = timeMetric.yearFrac(currFutureDate, nextFutureDate);
                int nbSteps = (int)(ceil(yearFrac * nbStepsPerYear));
                double delta = yearFrac / nbSteps;
                // create sim dates
                DateTime currDate = currFutureDate;
                fixDateOffsets[iOffset] = 0;
                for (int iStep = 1; iStep < nbSteps; ++iStep){
                    double notused;
                    currDate = timeMetric.impliedTime(currDate,
                                                      delta,
                                                      notused);
                    simDates->push_back(currDate);
                    fixDateOffsets[iOffset]++;
                }
                simDates->push_back(nextFutureDate);
                fixDateOffsets[iOffset]++;
                currFutureDate = nextFutureDate;
            }
            return simDates;
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Helper function used by both Generator and PathGenFineProcess.
        Creates a time line that is uniformally spaced between each
        fixing date aka 'futureDate'
        setting stepTypes in this version */
    static DateTimeArraySP genSimDates(const DateTimeArray& futureDates,
                                       const TimeMetric&    timeMetric,
                                       int                  nbStepsPerYear,
                                       list<VolSVCJ::STEP_TYPE>&     stepTypes) {
        static const string method("MCPathConfigSVCJ::Generator::genSimDates");
        try{
            // calculate approx number of sim dates
            int nbFutureDates = futureDates.size();
            ASSERT(nbFutureDates > 1);
            DateTime startDate = futureDates.front();
            DateTime horizonDate = futureDates.back();
            double horizonYearFrac = timeMetric.yearFrac(startDate, horizonDate);
            int approxNbSteps = (int)(ceil(horizonYearFrac * nbStepsPerYear));
            DateTimeArraySP simDates(new DateTimeArray());
            simDates->reserve(approxNbSteps + 1);
            // loop over future dates and insert additional points as needed
            DateTime currFutureDate = startDate;
            simDates->push_back(currFutureDate);
            stepTypes.push_back(VolSVCJ::START_DATE);
            for (int iDate = 1 ; iDate < nbFutureDates; iDate++){
                // calculate mesh size
                DateTime nextFutureDate = futureDates[iDate];
                double yearFrac = timeMetric.yearFrac(currFutureDate, nextFutureDate);
                int nbSteps = (int)(ceil(yearFrac * nbStepsPerYear));
                double delta = yearFrac / nbSteps;
                // create sim dates
                DateTime currDate = currFutureDate;
                for (int iStep = 1; iStep < nbSteps; ++iStep){
                    double notused;
                    currDate = timeMetric.impliedTime(currDate,
                                                      delta,
                                                      notused);
                    simDates->push_back(currDate);
                    stepTypes.push_back(VolSVCJ::DIFF_DATE);
                }
                simDates->push_back(nextFutureDate);
                stepTypes.push_back(VolSVCJ::SAMPLE_DATE);
                currFutureDate = nextFutureDate;
            }
            return simDates;
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MCPathConfigSVCJ, clazz);
        SUPERCLASS(MCPathConfigParametric);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        EMPTY_SHELL_METHOD(defaultMCPathConfigSVCJ);
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);
        FIELD(nbStepsPerYear, "Nb of steps per year");
        FIELD_MAKE_OPTIONAL(nbStepsPerYear);

        FIELD(cdfMappingParam, "CDF Implied Mapping parameters")
        FIELD_MAKE_OPTIONAL(cdfMappingParam);

        FIELD(spotDiscreteScheme, "Spot discretization scheme ("
                                         + VolSVJ_SpotDSTypeHelper::getNameList()
                                         + ")");
        FIELD_MAKE_OPTIONAL(spotDiscreteScheme);
        FIELD(spotDiscreteSchemeIndex, "");
        FIELD_MAKE_TRANSIENT(spotDiscreteSchemeIndex);
        FIELD(varDiscreteScheme, "Variance discretization scheme ("
                                         + VolSVJ_VarDSTypeHelper::getNameList()
                                         + ")");
        FIELD_MAKE_OPTIONAL(varDiscreteScheme);
        FIELD(varDiscreteSchemeIndex, "");
        FIELD_MAKE_TRANSIENT(varDiscreteSchemeIndex);
        FIELD_NO_DESC(simDates);
        FIELD_MAKE_TRANSIENT(simDates);

        FIELD_NO_DESC(cdfMapping);
        FIELD_MAKE_TRANSIENT(cdfMapping);

        FIELD(isCarefulRandoms, "Use careful randoms if true; don't otherwise");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms);
    }

    static IObject* defaultMCPathConfigSVCJ(){
        return new MCPathConfigSVCJ();
    }

    // for reflection
    MCPathConfigSVCJ():
    MCPathConfigParametric(TYPE),
    isCarefulRandoms(false),
    dependenceType("not used"),
    nbStepsPerYear(12),
    spotDiscreteScheme(VolSVJ_SpotDSTypeHelper::getDefaultName()),
    varDiscreteScheme(VolSVJ_VarDSTypeHelper::getDefaultName()),
    spotDiscreteSchemeIndex(-1),
    varDiscreteSchemeIndex(-1){}

    class Generator;
    friend class Generator;

#ifdef STATE_VARIABLES
    // Referee class
    class Gen;
    friend class Gen;
    // PathGen for spot
    class PathGenFineProcess;
    friend class PathGenFineProcess;
    typedef refCountPtr<PathGenFineProcess> PathGenFineProcessSP;
#endif

    // visible fields
    bool            isCarefulRandoms;
    string          dependenceType;
    int             nbStepsPerYear;
    CDFMappingSP    cdfMappingParam;
    string          spotDiscreteScheme;
    string          varDiscreteScheme;

    // transient fields
    int                spotDiscreteSchemeIndex;
    int                varDiscreteSchemeIndex;
    DateTimeArraySP    simDates; // cached between tweaks

    CDFMappingArray    cdfMapping;
};

CClassConstSP const MCPathConfigSVCJ::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigSVCJ", typeid(MCPathConfigSVCJ), load);


/** Creates elementary reference level data for Monte Carlo path generators
    e.g. fwds, fwdAtSimStart, reference level etc. */
class RefLevelDataSVCJ :    public virtual VirtualDestructorBase
{
public:
	/** Constructor where fwds are produced at product dates only */
    RefLevelDataSVCJ(const MCProductTimelineSP& timeline,
                 const IntArray& nbPaths,
                 const IMultiFactors* mAsset,
                 const IMCProduct* prod,
                 const MCPathGeneratorSP& pastPathGenerator);

    /** Destructor */
    ~RefLevelDataSVCJ() {}

    vector<DoubleArray>   refLevels;      //!< [iAsset][iPath]
    DoubleArray           fwdsAtSimStart; //!< By asset
    IRefLevel::IMCPathSP  refLevelPath;   //!< Works out refLevels

private:

    friend class MCPathBase;

	/** Disabled default constructor */
    RefLevelDataSVCJ();
};

typedef smartConstPtr<RefLevelDataSVCJ> RefLevelDataSVCJConstSP;
typedef smartPtr<RefLevelDataSVCJ> RefLevelDataSVCJSP;

RefLevelDataSVCJ::RefLevelDataSVCJ(const MCProductTimelineSP& timeline,
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

		// compute plain fwds at sim start
		CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));
		if (ProtEquity::TYPE->isInstance(asset) && !(AssetUtil::isBasket(asset))) {
			ProtEquityConstSP eq = ProtEquityConstSP::dynamicCast(asset);
			CAssetConstSP plainAsset = eq->getPlainAsset();
			fwdsAtSimStart[iAsset] = plainAsset->fwdValue(timeline->simStartDate);
		}
		else {
			fwdsAtSimStart[iAsset] = mAsset->factorFwdValue(iAsset, timeline->simStartDate);
		}
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



///////////////////////////////////////////////////////////////////////////////////
// Generator
///////////////////////////////////////////////////////////////////////////////////
class MCPathConfigSVCJ::Generator: public MCPathBase,
                                   virtual public IMCRandom::Callbacks,
                                   public DependenceMakerGauss::Support {
public:
    Generator(int                      numSimPaths,
              MCPathConfigSVCJ*         pathConfig,
              const MCPathGeneratorSP& pastPathGenerator,
              const IMCProduct*         prod):
    pathConfig(pathConfig),
    mAsset(prod->getMultiFactors()),
    nbAssets(prod->getNumAssets()),
    simDates(pathConfig->simDates),
    fwds(nbAssets),
    vols(nbAssets),
    timeMetrics(nbAssets),
    productPaths(nbAssets),
    maxDrifts(nbAssets, 1.0),
    isCarefulRandoms(pathConfig->carefulRandoms()),
    tradYrs(nbAssets),
    instVars(nbAssets),
    integratedVars(nbAssets),
	isQuanto(nbAssets, false),
	eqFXCorr(nbAssets, 0.0),
	fxVar(nbAssets) {
        static const string method("MCPathConfigSVCJ::Generator::Generator");
        try{
            // Obtain product timeline
            timeline = getProductTimeline(prod, pastPathGenerator);

            // Some consistency checks
            if (prod->getToday() != prod->getEffectiveSimStartDate()){
                throw ModelException(method, "fwd start style not supported");
            }
            ASSERT(timeline->futureDates[0] == prod->getEffectiveSimStartDate());

            // Obtain market data
            IntArray nbPaths(nbAssets, 1);
            refData = RefLevelDataSVCJSP(new RefLevelDataSVCJ(timeline, nbPaths, mAsset, prod, pastPathGenerator));

            // Initialize vols and time metrics etc
            int iAsset;
            simulateInitialVol = false;
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                // copy the vols that are in the multiasset
                vols[iAsset] = VolSVCJSP::dynamicCast(
                    VolRequestRaw::copyVolBase(*mAsset, iAsset));
                // copy time metric
                timeMetrics[iAsset] = TimeMetricSP(copy(&vols[iAsset]->getTimeMetric()));
                // and set flag for simulating initial vol
                simulateInitialVol |= vols[iAsset]->isRandomInitialVolatility();
            }

            // reuse sim dates between tweaks for numerical stability
            if (simDates->empty()){
                simDates = genSimDates(timeline->futureDates,
                                       *timeMetrics[0],       // pick a time metric !
                                       pathConfig->nbStepsPerYear,
                                       stepTypes);
            }

            // Initialize paths, i.e. logSpot, inst var and integrated var paths
            // Also, initialize productPaths.
            // NB: the latter contains the past (if any), the fomer don't
            logSpots.resize(timeline->futureDates.size());
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                productPaths[iAsset].resize(timeline->totalNumSteps);
                instVars[iAsset].resize(simDates->size()) ;
                integratedVars[iAsset].resize(simDates->size()) ;
            }

            // compute trading years
            DateTime startDate = (*simDates)[0];
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                tradYrs[iAsset].resize(simDates->size());
                list<double>::iterator t = tradYrs[iAsset].begin();
                for (int iDate = 0; iDate < simDates->size(); ++iDate, t++) {
                    *t = timeMetrics[iAsset]->yearFrac(startDate, (*simDates)[iDate]);

                }
            }

			// pre compute fwds, vol grid and quanto stuff
			preComputeFwds();

            // market factor jump intensity (sum of asset jump intensities)
            double jumpRate = 0.0;
            for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
            {
                jumpRate += vols[iAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE) *
                    (tradYrs[iAsset].back() / tradYrs[0].back()); // intensities adjusted to account for different
                                                                  // trading time conventions between assets
            }

            // arrival times
            dependenceJumpTimes.resize(1); // jump time process for market factor
            dependenceJumpTimes[0] = PoissonSP(new Poisson(jumpRate));

            const double QUANTILE_ERR = 0.00001; // 99.99% quantile of number of jumps
            double quantile;
            maxNumJumps = 0; // maximum number of market factor jumps

            int num = Quantile(tradYrs[0].back(),
                               jumpRate,
                               QUANTILE_ERR, (int)(1.0/QUANTILE_ERR), &quantile);
            if (maxNumJumps < num)
            {
                maxNumJumps = num;
            }

            // set up Dependence and DependenceJump
            dependence = pathConfig->dependenceMaker->createDependence(this);

            DependenceMakerGauss dependenceMakerJump;
            CDoubleMatrixSP correlationsJump = computeJumpCorrelations();
            //DependenceSP dependenceJump = dependenceMakerJump.createDependence(correlationsJump, this);
            DependenceSP dependenceJump(new Gauss(*correlationsJump));
            CDoubleMatrixSP correlationsJumpTimesMF = CDoubleMatrixSP(new CDoubleMatrix(1, 1));
            (*correlationsJumpTimesMF)[0][0] = 1.0;
            //DependenceSP dependenceJumpTimesMF =
            //    dependenceMakerJump.createDependence(correlationsJumpTimesMF, this);
            DependenceSP dependenceJumpTimesMF(new Gauss(*correlationsJumpTimesMF));

            // Initialize random number generator for spots
            // NB the size of the random deviates depends on the scheme
            // that is requested
            int nbSpotDeviates;
            if (pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                nbSpotDeviates = simDates->size() - 1;
            }
            else if(pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                nbSpotDeviates = logSpots.size() - 1;
            }
            else{
                throw ModelException(method,
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
            }
            int nbVarDeviates;
            if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else{
                throw ModelException(method,
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            // maxNumJumps added to num of random numbers to be generated.
            // could save random num generation if variable num fetching is possible
            randomGenSpot = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbSpotDeviates + maxNumJumps,
                nbAssets,
                timeline->totalNumPastDates));
            randomGenVar = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbVarDeviates + maxNumJumps,
                nbAssets,
                0));

            randomGenUniJumpTimesMF = IMCRandomSP(new MCRandomNoCache(
                0,
                dependenceJumpTimesMF,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                maxNumJumps,
                1, // market factor
                0));

            if (nbAssets > 1)
            {
                // no need to simulate these variables in the single-asset framework
                // since the asset participates to all market jumps
                randomGenUniJump = IMCRandomSP(new MCRandomNoCache(
                    0,
                    dependenceJump,
                    pathConfig->getRandomGenerator(),
                    pathConfig->carefulRandoms(),
                    maxNumJumps,
                    nbAssets,
                    0));
            }
            randomGenNormalJumpVar = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                maxNumJumps,
                nbAssets,
                0));
            randomGenNormalJumpSpot = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                maxNumJumps,
                nbAssets,
                0));
            // if (at least one) initialVol is to be simulated
            if (simulateInitialVol){
                randomGenInitVol.resize(nbAssets);
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    double shape = vols[iAsset]->gammaShape();
                    randomGenInitVol[iAsset] = MCRandGammaSP (new MCRandGamma(shape, 97));
                }
            }

            // marginal jump participation probabilities
            jumpParticipationProbs = computeJumpParticipationProbs();

            // do Implied Mapping if requested
            if (!!pathConfig->cdfMappingParam && timeline->futureDates.size()>0){
                // get market pdf calculator for each asset
                // create vol request first
                pathConfig->cdfMapping.resize(nbAssets);
                DateTimeArray futDates(timeline->futureDates);
                futDates.erase(futDates.begin());

                VolRequestLNStrikeSP volRequest(new LinearStrikeVolRequest(
                                               100.0, // any strike will do
                                               timeline->today,
                                               futDates.back(),
                                               false));

                PDFRequestLNStrike pdfLNStrkRequest(volRequest.get());

                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // get market pdf calculator for each asset
                    PDFDefaultLNStrikeSP pdfMarket(
                                          new PDFDefaultLNStrike(timeline->today,
															&(mAsset->getAsset(iAsset)),
															&pdfLNStrkRequest));

                    // get vol, manage memory
                    CVolProcessedSP vol(vols[iAsset]->getProcessedVol(volRequest.get(),
                                                                   &(mAsset->getAsset(iAsset))));
                    // get model pdf calculator for each asset
                    PDFCalculatorSP pdfModel(
                                        PDFCalculatorMaker::makePDFCalculator("PDFFourier",
                                                                         timeline->today,
                                                                         pathConfig->cdfMappingParam->getModel().get(),
                                                                         &(mAsset->getAsset(iAsset)),
                                                                         vol.get()));

                    pathConfig->cdfMapping[iAsset] = CDFMappingSP(
                            new CDFMapping(futDates,
                                pathConfig->cdfMappingParam.get()));
                    // creates model and market spots grids
                    pathConfig->cdfMapping[iAsset]->setSpotsGrids(futDates,
												 //CAssetConstSP(&(mAsset->getAsset(iAsset))),
												 pathConfig->cdfMappingParam->getAsset(
														mAsset->getAsset(iAsset).getName()),
												 volRequest,
                                                 timeline->today,
                                                 timeline->today,
                                                 pdfMarket,
                                                 pdfModel);

                    // build the spline interpolant from model spots to market spots
                    // cdfMapping[iAsset]->buildSplineInterpolant();

                    // build the linear interpolant from model spots to market spots
                    pathConfig->cdfMapping[iAsset]->buildLinearInterpolant();
                }
            }

        } catch (exception& e){
            throw ModelException(e, method);
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

    /** Computes the jump correlation matrix used to correlate the marginal jump participation probabilities */
    CDoubleMatrixSP computeJumpCorrelations() const {
        DoubleArraySP jumpParticipationProbs = computeJumpParticipationProbs();
        CDoubleMatrixConstSP correlations = mAsset->factorsCorrelationMatrix();
        CDoubleMatrixSP jumpCorrelations = CDoubleMatrixSP(new CDoubleMatrix(nbAssets, nbAssets));

        // computes jump correlations
        // this matrix is used in the Gaussian copula meant to correlate participation to market factor jumps
        int iAsset;
        for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
        {
            int jAsset;
            for (jAsset = iAsset ; jAsset < nbAssets ; jAsset++)
            {
                // iAsset = jAsset, correlation = 1.0
                if (iAsset == jAsset)
                {
                    (*jumpCorrelations)[iAsset][jAsset] = 1.0;
                }
                // if iAsset != jAsset, compute jump correlation based on asset correlation
                else
                {
                    JumpCorrelationData JumpCorrelationStruct;

                    // jump intensities for assets i anf j
                    double intensityIAsset = vols[iAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE);
                    double intensityJAsset = vols[jAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE);

                    // set lowJumpIntensityAsset and highJumpIntensityAsset
                    if (intensityIAsset <= intensityJAsset)
                    {
                        // lower jump marginal probability
                        JumpCorrelationStruct.lowJumpProba = (*jumpParticipationProbs)[iAsset];
                        // lower jump marginal intensity
                        JumpCorrelationStruct.lowJumpIntensity = intensityIAsset;
                        // higher jump marginal probability
                        JumpCorrelationStruct.highJumpProba = (*jumpParticipationProbs)[jAsset];
                        // higher jump marginal intensity
                        JumpCorrelationStruct.highJumpIntensity = intensityJAsset;
                    }
                    else
                    {
                        // lower jump marginal probability
                        JumpCorrelationStruct.lowJumpProba = (*jumpParticipationProbs)[jAsset];
                        // lower jump marginal intensity
                        JumpCorrelationStruct.lowJumpIntensity = intensityJAsset;
                        // higher jump marginal probability
                        JumpCorrelationStruct.highJumpProba = (*jumpParticipationProbs)[iAsset];
                        // higher jump marginal intensity
                        JumpCorrelationStruct.highJumpIntensity = intensityIAsset;
                    }

                    // lowJumpProba quantile of N1
                    JumpCorrelationStruct.normInvLowProba = N1Inverse(JumpCorrelationStruct.lowJumpProba);
                    // highJumpProba quantile of N1
                    JumpCorrelationStruct.normInvHighProba = N1Inverse(JumpCorrelationStruct.highJumpProba);
                    // correlation between the corresponding asset
                    JumpCorrelationStruct.correlation = (*correlations)[iAsset][jAsset];

                    // target conditional jump probability participation
                    double ratio = JumpCorrelationStruct.highJumpIntensity / JumpCorrelationStruct.lowJumpIntensity;
                    double correlation = JumpCorrelationStruct.correlation;
                    JumpCorrelationStruct.targetProba = Maths::min(Maths::max(ratio * correlation, 0.0), 1.0);
                    // JumpCorrelationStruct.targetProba = Maths::max(correlation, 0.0);

                    // jump correlation is equal to 1.0 if target probability is 1.0
                    if (Maths::isZero(JumpCorrelationStruct.targetProba - 1.0))
                    {
                        (*jumpCorrelations)[iAsset][jAsset] = 1.0;
                    }
                    // jump correlation is equal to 1.0 if target probability is -1.0
                    else if (Maths::isZero(JumpCorrelationStruct.targetProba + 1.0))
                    {
                        (*jumpCorrelations)[iAsset][jAsset] = -1.0;
                    }
                    else
                    {
                        (*jumpCorrelations)[iAsset][jAsset] = zbrentUseful(&objectivFunc,          /* (I) The function to find the root of */
                                                                           &JumpCorrelationStruct, /* (I) Parameter block */
                                                                           -1.0,                   /* (I) Low value for x */
                                                                           1.0,                    /* (I) High value for x */
                                                                           0.001);                 /* (I) Tolerance (0.1% correlation) */
                    }
                    (*jumpCorrelations)[iAsset][jAsset] = max(-0.99, min(0.99, (*jumpCorrelations)[iAsset][jAsset]));
                    (*jumpCorrelations)[jAsset][iAsset] = (*jumpCorrelations)[iAsset][jAsset];
                }
            }
        }

        // transform jumpCorrelations into a true correlation matrix if needed
        double squareErr;
        CDoubleMatrix correl(jumpCorrelations->symmToCorrel(&squareErr, 0.001));
        CDoubleMatrixSP jumpCorrs(copy(&correl)); // eigenValueFloor set to 0.0
        return jumpCorrs;
    }

    /** Computes the jump participation marginal probabilities */
    DoubleArraySP computeJumpParticipationProbs() const {
        DoubleArraySP jumpParticipationProbs = DoubleArraySP(new DoubleArray(nbAssets));

        // market factor jump intensity (sum of asset jump intensities)
        double jumpRateMF = 0.0;
        DoubleArray jumpRates(nbAssets);

        int iAsset;
        for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
        {
            // intensities adjusted to account for different trading time conventions between assets
            jumpRates[iAsset] = vols[iAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE) *
                (tradYrs[iAsset].back() / tradYrs[0].back());
            jumpRateMF += jumpRates[iAsset];
        }

        // marginal jump participation marginal probabilities
        for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
        {
            (*jumpParticipationProbs)[iAsset] = jumpRates[iAsset] / jumpRateMF;
        }

        return jumpParticipationProbs;
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigSVCJ::Generator::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

protected:

	friend class MCPathConfigSVCJ;

    /** Configures pathGen for antithetics */
    void configureAntithetics() {
    };

    /** Configures pathGen for nonAntithetics. We need a better name */
    void configureNonAntithetics() {
    };

    /** Draws random numbers */
    virtual void drawRandomNumbers(int pathIdx){
        // Draw random numbers for paths and jump sizes
        randomGenSpot->generate(pathIdx);
        randomGenVar->generate(pathIdx);
        randomGenUniJumpTimesMF->generate(pathIdx);
        if (nbAssets > 1)
        {
            randomGenUniJump->generate(pathIdx);
        }
        randomGenNormalJumpVar->generate(pathIdx);
        randomGenNormalJumpSpot->generate(pathIdx);
    }

    virtual void generatePath(int pathIdx, int iAsset, int iPath) {
        static const string method = "MCPathConfigSVCJ::Generator::generatePath";
        try{
            const DoubleMatrix& randomSpot = randomGenSpot->getRandomNumbers();
            const DoubleMatrix& randomVol = randomGenVar->getRandomNumbers();
            const DoubleMatrix& randomUniJumpTimesMFMatrix = randomGenUniJumpTimesMF->getRandomNumbers();
            const DoubleMatrix& randomNormalJumpVar = randomGenNormalJumpVar->getRandomNumbers();
            const DoubleMatrix& randomNormalJumpSpot = randomGenNormalJumpSpot->getRandomNumbers();

            double randomInitVar = 0.0;
            if (vols[iAsset]->isRandomInitialVolatility()){
                randomGenInitVol[iAsset]->fetch(1, &randomInitVar);
            }

            int varDSType = pathConfig->varDiscreteSchemeIndex;
            int spotDSType = pathConfig->spotDiscreteSchemeIndex;
            int numJumps = 0;

            DoubleMatrix randomJumpTimes(1,maxNumJumps); // 1 column, maxNumJumps rows
             // simulate jump arrival times
            if (maxNumJumps >0){
                int iJump;
                for (iJump = 0; iJump < maxNumJumps; ++iJump)
                    randomJumpTimes[0][iJump] = randomUniJumpTimesMFMatrix[0][iJump];
                // arrival times
                dependenceJumpTimes[0]->correlateSeries(randomJumpTimes, pathIdx);
                // truncate jump dates to last sim date
                for (iJump = 0; iJump < maxNumJumps; ++iJump){
                    if (randomJumpTimes[0][iJump] > tradYrs[0].back()){
                        numJumps = iJump;
                        break;
                    }
                }
            }

            BoolArray jumpParticipation(maxNumJumps); // maxNumJumps elements
            if (nbAssets > 1)
            {
                const DoubleMatrix& randomUniJumpMatrix = randomGenUniJump->getRandomNumbers();
                int iJump;
                for (iJump = 0 ; iJump < maxNumJumps ; ++iJump)
                {
                    double uniformParticipation = Maths::max(N1(randomUniJumpMatrix[iAsset][iJump]), 0.00000001);
                    jumpParticipation[iJump] =
                        (uniformParticipation < (*jumpParticipationProbs)[iAsset]) ? true : false;
                }
            }
            else
            {
                int iJump;
                for (iJump = 0 ; iJump < maxNumJumps ; ++iJump)
                {
                    jumpParticipation[iJump] = true;
                }
            }

            // insert jumps: jumpTypesAdded and
            vector<VolSVCJ::STEP_TYPE_ITER> jumpTypesAdded;
            // jumpsAdded[0] -> tradYrs, [1] -> instVars, [2] -> integratedVars, [3] -> logSpots */
            vector<vector<VolSVCJ::LIST_ITER> >jumpsAdded(4);
            vector<double> varJumps; // keeps values of instVar jumps
            vols[iAsset]->addJumps(tradYrs[iAsset],
                                   randomJumpTimes[0],
                                   numJumps,
                                   jumpParticipation,
                                   randomNormalJumpVar[iAsset],
                                   randomNormalJumpSpot[iAsset],
                                   stepTypes,
                                   instVars[iAsset],
                                   integratedVars[iAsset],
                                   logSpots,
                                   varJumps,
                                   jumpTypesAdded,
                                   jumpsAdded,
                                   randomInitVar);

            // generate the var paths
            if (varDSType == VolSVJ_VarDSTypeHelper::EnumList::EULER)
                vols[iAsset]->generateVarPathsEuler(tradYrs[iAsset].begin(),
                                                    tradYrs[iAsset].size(),
                                                    randomVol[iAsset],
                                                    instVars[iAsset].begin(),
                                                    integratedVars[iAsset].begin());
            else if (varDSType == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER)
                vols[iAsset]->generateVarPathsTransEuler(tradYrs[iAsset].begin(),
                                                         tradYrs[iAsset].size(),
                                                         randomVol[iAsset],
                                                         instVars[iAsset].begin(),
                                                         integratedVars[iAsset].begin());
            else
                throw ModelException("VolSVCJ::MCGenerateVarPaths",
                                     "unexpected variance discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(varDSType));

            // then generate the spot path
            if (spotDSType == VolSVJ_SpotDSTypeHelper::EnumList::EULER)
                vols[iAsset]->generateSpotPathsEuler(tradYrs[iAsset].begin(),
                                                     tradYrs[iAsset].size(),
                                                     stepTypes.begin(),
                                                     randomSpot[iAsset],
                                                     randomVol[iAsset],
                                                     logSpots.begin(),
                                                     instVars[iAsset].begin(),
                                                     integratedVars[iAsset].begin());
            else if (spotDSType == VolSVJ_SpotDSTypeHelper::EnumList::EXACT)
                vols[iAsset]->generateSpotPathsExact(tradYrs[iAsset].begin(),
                                                     tradYrs[iAsset].size(),
                                                     stepTypes.begin(),
                                                     randomSpot[iAsset],
                                                     varJumps,
                                                     logSpots.begin(),
                                                     instVars[iAsset].begin(),
                                                     integratedVars[iAsset].begin());
            else
                throw ModelException("VolSVCJ::MCGenerateExJumpPaths",
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(spotDSType));

            // remove inserted jump points
            if (numJumps>0){
                int nStepsAdded = jumpsAdded[0].size();

                for (int j = 0; j < nStepsAdded; j++){
                    stepTypes.erase(jumpTypesAdded[j]);
                    tradYrs[iAsset].erase(jumpsAdded[0][j]);
                    instVars[iAsset].erase(jumpsAdded[1][j]);
                    integratedVars[iAsset].erase(jumpsAdded[2][j]);
                    logSpots.erase(jumpsAdded[3][j]);
                }
            }

			// Compute Quanto adjustment
			DoubleArray adjQuanto(logSpots.size(), 0);
			if (isQuanto[iAsset]) {
				int nbSimDates = stepTypes.size();
				VolSVCJ::LIST_ITER instV = instVars[iAsset].begin();
				VolSVCJ::STEP_TYPE_ITER stepTypeIter = stepTypes.begin();
				int idxAdj = 0;
				instV++;
				stepTypeIter++;
				for (int iSimDates = 1; iSimDates < nbSimDates; iSimDates++, instV++, stepTypeIter++) {
					adjQuanto[idxAdj] -= eqFXCorr[iAsset] * fxVar[iAsset][iSimDates - 1] * sqrt(Maths::max((*instV),0.0));
					if (*stepTypeIter == VolSVCJ::SAMPLE_DATE) {
						idxAdj++;
						adjQuanto[idxAdj] = adjQuanto[idxAdj - 1];
					}
				}
			}

            // then, populate spot path
            const DoubleArray& fwds = this->fwds[iAsset];
            DoubleArray& spots = productPaths[iAsset];
            int nbFutFixDates = logSpots.size();
            // NB: spots[timeline->numPastDates] corresponds to futureDates[1], i.e.
            // the first future product/fixing date -- not today
            VolSVCJ::LIST_ITER logS = logSpots.begin();
            logS++;
            for (int iFutFixDates = 1, iFixDates = timeline->numPastDates;
                 iFutFixDates < nbFutFixDates; ++iFutFixDates, ++iFixDates, logS++) {
					spots[iFixDates] = fwds[iFutFixDates] * exp(*logS) * exp(adjQuanto[iFutFixDates - 1]);
	        }

            // do Implied Mapping if requested
            if (!!pathConfig->cdfMappingParam && timeline->futureDates.size()>0){
                double* s = &*(spots.begin() + timeline->numPastDates);

                pathConfig->cdfMapping[iAsset]->impliedMCSpot(s,
                                                            s,
                                                            0,
                                                            timeline->numFutSteps);
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Obtains timeline object from base */
    MCProductTimelineConstSP getTimeline() const {
        return timeline;
    }

    /** Returns number of assets */
    int NbSimAssets() const {
        return nbAssets;
    }

    /** Returns the reference level for iAsset, iPath */
    const double& refLevel(int iAsset, int iPath) const {
        return refData->refLevels[iAsset][0];
    }

    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath) {
        return refData->refLevels[iAsset][0];
    }

    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const {
        return refData->refLevelPath;
    }

    /** Returns the path for iAsset, iPath */
    const double* Path(int iAsset, int iPath) const {
        return &productPaths[iAsset][0];
    }

    /** Returns the path for iAsset, iPath */
    double* Path(int iAsset, int iPath)  {
        return &productPaths[iAsset][0];
    }

    /** Returns the number of paths per asset */
    int nbPaths(int iAsset) const {
        return 1;
    }

    /** Returns if it is a single path generator or not */
    bool isSinglePath() const {
        return true;
    }

	/** pre-compute fwds arrays, get quanto data here for now */
    void preComputeFwds() {
        // for each asset
        for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
            // need to cache the forward at all simulation dates
            fwds[iAsset].resize(timeline->numFutSteps + 1);
            // current way of dealing with quanto
            CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));

            if (ProtAsset::TYPE->isInstance(asset)) {
                throw ModelException("preComputeFwds",
                                     "Protected assets (" + asset->getTrueName() +
                                     ") of type " + asset->getClass()->getName() +
                                         " are not supported yet");
            }

            isQuanto[iAsset] = ProtEquity::TYPE->isInstance(asset);
            if (isQuanto[iAsset]) {
                fxVar[iAsset] = DoubleArray(simDates->size()-1);
                ProtEquityConstSP eq = ProtEquityConstSP::dynamicCast(asset);
                // compute plain fwds
                CAssetConstSP plainAsset = eq->getPlainAsset();
                plainAsset->fwdValue(timeline->futureDates, fwds[iAsset]);
                // get fx vol and corr
                FXVolBaseWrapper fxVol = eq->getFXVol();
                eqFXCorr[iAsset] = eq->getCorrelation()->getCorrelation();
                // interpolate the fx vols atm
                ATMVolRequestSP fxVolRequest(new ATMVolRequest());
                // interpolate the vol
                const Asset* fx = eq->getFXAsset();
                CVolProcessedSP  fxVoltmp(fxVol->getProcessedVol(fxVolRequest.get(), fx));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(fxVoltmp);
                volFXBS->CalcVol(*simDates, CVolProcessedBS::forward, fxVar[iAsset]);
                // calculate fx vars
				const TimeMetric& metric = vols[iAsset]->getTimeMetric();
                for (int iStep = 0; iStep < fxVar[iAsset].size(); iStep++) {
				    double dt = metric.yearFrac((*simDates)[iStep], (*simDates)[iStep+1]);
                    fxVar[iAsset][iStep] *= dt; // using asset dt
                }
		    }
            else{
                mAsset->factorFwdValues(iAsset, timeline->futureDates, fwds[iAsset]);
            }
        }
    }

    // fields
    MCPathConfigSVCJ*       pathConfig;
    const IMultiFactors*    mAsset;            // multi asset interface
    int                     nbAssets;          // number of assets
    DateTimeArraySP         simDates;          // fine grid of simulation dates
    list<VolSVCJ::STEP_TYPE>stepTypes;         // keeps type of each sim date
    vector<DoubleArray>     fwds;              // fwds [iAsset][iStep]
    VolSVCJArray            vols;              // array of VolSVCJs
    TimeMetricArray         timeMetrics;       // array of time metrics
    vector<DoubleArray>     productPaths;      // spot path [iAsset][iStep]
    IMCRandomSP             randomGenSpot;     // Random number generator for spots
    IMCRandomSP             randomGenVar;      // Random number generator for vars
    IMCRandomSP             randomGenUniJumpTimesMF; // Random number generator for jump times uniform deviates
    IMCRandomSP             randomGenUniJump; // Random number generator for jump participations uniform deviates
    DoubleArraySP           jumpParticipationProbs; // marginal jump participation probability
    IMCRandomSP             randomGenNormalJumpVar;  // Random number generator for jump normal deviates
    IMCRandomSP             randomGenNormalJumpSpot;  // Random number generator for jump normal deviates
    bool                    simulateInitialVol;     // if init vol is simulated (stationary distribution)
    vector<MCRandGammaSP>   randomGenInitVol;  // Random number generator for init vol - stationary gamma distribution
    DependenceSP            dependence;     // Dependence object
    MCProductTimelineSP     timeline;          // Product timeline
    RefLevelDataSVCJSP      refData;           // RefLevel data
    DoubleArray             maxDrifts;         // XXX product of MAX(1, drifts)
    bool                    isCarefulRandoms;  // XXX

    int                     maxNumJumps;
    vector<list<double> >   tradYrs;		   // fine grid of trad years
    list<double>            logSpots;          // spot path [iAsset][iStep]
    vector<list<double> >   instVars;          // instantaneous variance path [iAsset][iStep]
    vector<list<double> >   integratedVars;    // integrated variance path [iAsset][iStep]
    vector<PoissonSP>       dependenceJumpTimes;

    BoolArray               isQuanto;
    DoubleArray             eqFXCorr;
    vector<DoubleArray>     fxVar;
};

#ifdef STATE_VARIABLES
/** Spot path generator using state variable approach */
class MCPathConfigSVCJ::PathGenFineProcess: virtual public MCPathGen,
                                            virtual public IMCRandom::Callbacks,
                                            public DependenceMakerGauss::Support {
public:
    PathGenFineProcess(int                             numSimPaths,
                       MCPathConfigSVCJ*                 pathConfig,
                       const PastPathGenSpotSP&        pastPathGenSpot,
                       const MCProductClient*          prodClient,
                       const SVGenSpotArray&              spotGenArray,
                       const MCQuadVarArray&           quadVarGenArray,
                       const MCSqrtAnnualQuadVarArray& sqrtAnnualQuadVarGenArray,
                       StateVarDBase&                  svDBase):
    pathConfig(pathConfig),
    mAsset(prodClient->getMultiFactors()),
    nbAssets(prodClient->getNumAssets()),
    simDates(pathConfig->simDates),
    fwds(nbAssets),
    vols(nbAssets),
    timeMetrics(nbAssets),
    tradYears(nbAssets),
    instVars(nbAssets),
    maxDrifts(nbAssets, 1.0),
    integratedVars(nbAssets),
	isQuanto(nbAssets, false),
	eqFXCorr(nbAssets, 0.0),
	fxVar(nbAssets),
    simHasPast(pastPathGenSpot->hasPast()),
    spotProdPaths(nbAssets),
    spotProdOffset(nbAssets, 0),
    quadVarProdPaths(nbAssets),
    quadVarProdOffset(nbAssets, 0),
    sqrtAnnualQuadVarProdPaths(nbAssets),
    sqrtAnnualQuadVarProdOffset(nbAssets, 0) {

        static const string method("MCPathConfigSVCJ::PathGenFineProcess::Generator");
        try{
            const DateTime& today = prodClient->getToday();

            // Create simulation timeline
            DateTimeArraySP spotGenDates = MCPath::getAllDates(spotGenArray);
            DateTimeArraySP quadVarGenDates = MCPath::getAllDates(quadVarGenArray);
            DateTimeArraySP sqrtAnnualQuadVarGenDates = MCPath::getAllDates(sqrtAnnualQuadVarGenArray);
            DateTimeArray allDates = DateTime::merge(*quadVarGenDates, *spotGenDates);
            allDates = DateTime::merge(*sqrtAnnualQuadVarGenDates, allDates);
            simTimeline = today.getFutureDates(allDates);
            numFutSteps = simTimeline.size();
            if (!numFutSteps){
                throw ModelException(method,
                                     "there is no future date to simulate");
            }
            simTimeline.insert(simTimeline.begin(), today);

            // Obtain market data
            IntArray nbPaths(nbAssets, 1);

            // Initialize vols and time metrics etc
            int iAsset;
            simulateInitialVol = false;
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                // copy the vols that are in the multiasset
                vols[iAsset] = VolSVCJSP::dynamicCast(
                    VolRequestRaw::copyVolBase(*mAsset, iAsset));
                // copy time metric
                timeMetrics[iAsset] = TimeMetricSP(copy(&vols[iAsset]->getTimeMetric()));
                // and set flag for simulating initial vol
                simulateInitialVol |= vols[iAsset]->isRandomInitialVolatility();
            }

            // reuse sim dates between tweaks for numerical stability
            if (simDates->empty()){
                simDates = genSimDates(simTimeline,
                                       *timeMetrics[0],       // pick a time metric !
                                       pathConfig->nbStepsPerYear,
                                       fixDateOffsets);
            }

            // Initialize paths, i.e. logSpot, inst var and integrated var paths
            logSpots.resize(simTimeline.size());
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                instVars[iAsset].resize(simDates->size()) ;
                integratedVars[iAsset].resize(simDates->size()) ;
            }

            if (spotGenArray.size()){
                // Create spot paths per asset
                DateTimeArrayArray spotDatesPerAsset(nbAssets);
                vector<const double*> spotPtrs(nbAssets);
                vector<int> spotBeginInd(nbAssets);
                vector<int> spotEndInd(nbAssets);

                const IPastValues* pastValues = prodClient->getMCPastValues();
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // Populate spot past values
                    DateTimeArraySP spotProdDates = MCPath::getAllDates(spotGenArray, iAsset);
                    spotProdPaths[iAsset] = DoubleArray(spotProdDates->size());
                    DateTimeArray spotPastProdDates = today.getPastDates(*spotProdDates);
                    DoubleArray spotPastValues = pastValues->getPastValues(spotPastProdDates, iAsset, today);
                    spotProdOffset[iAsset] = spotPastValues.size();
                    int iStep;
                    for(iStep = 0; iStep < spotPastValues.size(); iStep++) {
                        spotProdPaths[iAsset][iStep] = spotPastValues[iStep];
                    }

                    // Create spot mappings
                    spotDatesPerAsset[iAsset] = *spotProdDates;
                    spotPtrs[iAsset] = &spotProdPaths[iAsset][0];
                    spotBeginInd[iAsset] = spotPastProdDates.size();
                    spotEndInd[iAsset]   = spotProdDates->size();
                }

                // Create SVGenSpot::IStateVars and put them in database
                MCPath::IStateVarArray spotSVCJArray(MCPath::createPaths(
                    false,
                    spotGenArray,
                    spotDatesPerAsset,
                    spotBeginInd,
                    spotEndInd,
                    spotPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < spotGenArray.size(); iVar++) {
                    svDBase.append(spotGenArray[iVar], spotSVCJArray[iVar]);
                }
            }

            if (quadVarGenArray.size()){
                // Create quad var paths per asset
                DateTimeArrayArray quadVarDatesPerAsset(nbAssets);
                vector<const double*> quadVarPtrs(nbAssets);
                vector<int> quadVarBeginInd(nbAssets);
                vector<int> quadVarEndInd(nbAssets);

                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // populate quad var past values
                    DateTimeArraySP quadVarProdDates = MCPath::getAllDates(quadVarGenArray, iAsset);
                    quadVarProdPaths[iAsset] = DoubleArray(quadVarProdDates->size());
                    DateTimeArray quadVarPastProdDates = today.getPastDates(*quadVarProdDates);
                    if (quadVarPastProdDates.size()){
                        throw ModelException(method,
                                             "past quadratic variation samplings not currently supported");
                    }
                    quadVarProdOffset[iAsset] = 0;

                    // Create quadVar mappings
                    quadVarDatesPerAsset[iAsset] = *quadVarProdDates;
                    quadVarPtrs[iAsset] = &quadVarProdPaths[iAsset][0];
                    quadVarBeginInd[iAsset] = quadVarPastProdDates.size();
                    quadVarEndInd[iAsset]   = quadVarProdDates->size();
                }

                // Create MCQuadVar::IStateVars and put them in database
                MCPath::IStateVarArray quadVarSVCJArray(MCPath::createPaths(
                    false,
                    quadVarGenArray,
                    quadVarDatesPerAsset,
                    quadVarBeginInd,
                    quadVarEndInd,
                    quadVarPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < quadVarGenArray.size(); iVar++) {
                    svDBase.append(quadVarGenArray[iVar], quadVarSVCJArray[iVar]);
                }
            }

            if (sqrtAnnualQuadVarGenArray.size()){
                // Create quad var paths per asset
                DateTimeArrayArray sqrtAnnualQuadVarDatesPerAsset(nbAssets);
                vector<const double*> sqrtAnnualQuadVarPtrs(nbAssets);
                vector<int> sqrtAnnualQuadVarBeginInd(nbAssets);
                vector<int> sqrtAnnualQuadVarEndInd(nbAssets);

                // some historical values for the past
                pastAnnualQuadVars.resize(nbAssets);
                pastWeights.resize(nbAssets);

                const IPastValues* pastValues = prodClient->getMCSqrtAnnualQuadVarPastValues();
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // populate quad var past values
                    DateTimeArraySP sqrtAnnualQuadVarProdDates = MCPath::getAllDates(sqrtAnnualQuadVarGenArray, iAsset);
                    sqrtAnnualQuadVarProdPaths[iAsset] = DoubleArray(sqrtAnnualQuadVarProdDates->size());
                    DateTimeArray sqrtAnnualQuadVarPastProdDates = today.getPastDates(*sqrtAnnualQuadVarProdDates);
                    DoubleArray sqrtAnnualQuadVarPastValues
                        = pastValues->getPastValues(sqrtAnnualQuadVarPastProdDates, iAsset, today);
                    sqrtAnnualQuadVarProdOffset[iAsset] = sqrtAnnualQuadVarPastValues.size();
                    int iStep;
                    for(iStep = 0; iStep < sqrtAnnualQuadVarPastValues.size(); iStep++) {
                        sqrtAnnualQuadVarProdPaths[iAsset][iStep] = sqrtAnnualQuadVarPastValues[iStep];
                    }

                    // some precomputation need to be done to evaluate total variance
                    // when there is a past
                    if (sqrtAnnualQuadVarProdOffset[iAsset] > 0){
                        // get the realized var at valueDate
                        DateTimeArray valueDate(1, (*simDates)[0]);
#if 1
                        DoubleArray pastVols = pastValues->getPastValues(valueDate, iAsset, valueDate[0]);
                        pastAnnualQuadVars[iAsset] = Maths::square(pastVols[0]);
#else
                        DoubleArray pastAnnualQuadVar = pastValues->getPastValues(valueDate, iAsset, valueDate[0]);
                        pastAnnualQuadVars[iAsset] = pastAnnualQuadVar[0];
#endif
                        // get nb of past values up to and including valueDate
                        double nbPastValues = pastValues->getNbPastValues(valueDate[0], iAsset);
                        // get total nb of past values
                        double nbValues = pastValues->getNbPastValues(iAsset);
                        pastWeights[iAsset] = (nbPastValues - 1.0) / (nbValues - 1.0);
                    }

                    // Create sqrtAnnualQuadVar mappings
                    sqrtAnnualQuadVarDatesPerAsset[iAsset] = *sqrtAnnualQuadVarProdDates;
                    sqrtAnnualQuadVarPtrs[iAsset] = &sqrtAnnualQuadVarProdPaths[iAsset][0];
                    sqrtAnnualQuadVarBeginInd[iAsset] = sqrtAnnualQuadVarPastProdDates.size();
                    sqrtAnnualQuadVarEndInd[iAsset] = sqrtAnnualQuadVarProdDates->size();
                }

                // Create MCQuadVar::IStateVars and put them in database
                MCPath::IStateVarArray sqrtAnnualQuadVarSVCJArray(MCPath::createPaths(
                    false,
                    sqrtAnnualQuadVarGenArray,
                    sqrtAnnualQuadVarDatesPerAsset,
                    sqrtAnnualQuadVarBeginInd,
                    sqrtAnnualQuadVarEndInd,
                    sqrtAnnualQuadVarPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < sqrtAnnualQuadVarGenArray.size(); iVar++) {
                    svDBase.append(sqrtAnnualQuadVarGenArray[iVar], sqrtAnnualQuadVarSVCJArray[iVar]);
                }
            }

            // compute trading years
            DateTime startDate = (*simDates)[0];
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                tradYears[iAsset].resize(simDates->size());
                for (int iDate = 0; iDate < simDates->size(); ++iDate) {
                    tradYears[iAsset][iDate]
                        = timeMetrics[iAsset]->yearFrac(startDate,
                                                        (*simDates)[iDate]);
                }
            }

            // pre compute fwds, vol grid and quanto stuff
			preComputeFwds();

            // market factor jump intensity (sum of asset jump intensities)
            double jumpRate = 0.0;
            for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
            {
                jumpRate += vols[iAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE) *
                    (tradYrs[iAsset].back() / tradYrs[0].back()); // intensities adjusted to account for different
                                                                  // trading time conventions between assets
            }

            // arrival times
            dependenceJumpTimes.resize(1); // jump time process for market factor
            dependenceJumpTimes[0] = PoissonSP(new Poisson(jumpRate));

            const double QUANTILE_ERR = 0.00001; // 99.99% quantile of number of jumps
            double quantile;
            maxNumJumps = 0; // maximum number of market factor jumps

            int num = Quantile(tradYrs[0].back(),
                               jumpRate,
                               QUANTILE_ERR, (int)(1.0/QUANTILE_ERR), &quantile);
            if (maxNumJumps < num)
            {
                maxNumJumps = num;
            }

            // set up Dependence and DependenceJump
            dependence = pathConfig->dependenceMaker->createDependence(this);

            DependenceMakerGauss dependenceMakerJump;
            CDoubleMatrixSP correlationsJump = computeJumpCorrelations();
            //DependenceSP dependenceJump = dependenceMakerJump.createDependence(correlationsJump, this);
            DependenceSP dependenceJump(new Gauss(*correlationsJump));
            CDoubleMatrixSP correlationsJumpTimesMF = CDoubleMatrixSP(new CDoubleMatrix(1, 1));
            (*correlationsJumpTimesMF)[0][0] = 1.0;
            //DependenceSP dependenceJumpTimesMF =
            //    dependenceMakerJump.createDependence(correlationsJumpTimesMF, this);
            DependenceSP dependenceJumpTimesMF(new Gauss(*correlationsJumpTimesMF));

            // Initialize random number generator for spots
            // NB the size of the random deviates depends on the scheme
            // that is requested
            int nbSpotDeviates;
            if (pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EULER){
                nbSpotDeviates = simDates->size() - 1;
            }
            else if(pathConfig->spotDiscreteSchemeIndex == VolSVJ_SpotDSTypeHelper::EnumList::EXACT){
                nbSpotDeviates = logSpots.size() - 1;
            }
            else{
                throw ModelException(method,
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
            }
            int nbVarDeviates;
            if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else{
                throw ModelException(method,
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }

            // maxNumJumps added to num of random numbers to be genereated.
            // could save random num generation if variable num fetching is poosible
            randomGenSpot = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbSpotDeviates + maxNumJumps,
                nbAssets,
                spotGenDates->size() - numFutSteps));
            randomGenVar = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbVarDeviates + maxNumJumps,
                nbAssets,
                0));

            randomGenUniJumpTimesMF = IMCRandomSP(new MCRandomNoCache(
                0,
                dependenceJumpTimesMF,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                maxNumJumps,
                1, // market factor
                0));
            if (nbAssets > 1)
            {
                // no need to simulate these variables in the single-asset framework
                // since the asset participates to all market jumps
                randomGenUniJump = IMCRandomSP(new MCRandomNoCache(
                    0,
                    dependenceJump,
                    pathConfig->getRandomGenerator(),
                    pathConfig->carefulRandoms(),
                    maxNumJumps,
                    nbAssets,
                    0));
            }
            randomGenNormalJumpVar = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                maxNumJumps,
                nbAssets,
                0));
            randomGenNormalJumpSpot = IMCRandomSP(new MCRandomNoCache(
                0,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                maxNumJumps,
                nbAssets,
                0));
            // if (at least one) initialVol is to be simulated
            if (simulateInitialVol){
                randomGenInitVol.resize(nbAssets);
                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    double shape = vols[iAsset]->gammaShape();
                    randomGenInitVol[iAsset] = MCRandGammaSP (new MCRandGamma(shape, 97));
                }
            }

            // marginal jump participation probabilities
            jumpParticipationProbs = computeJumpParticipationProbs();

            // do Implied Mapping if requested
            if (!!pathConfig->cdfMappingParam && numFutSteps>0){
                // get market pdf calculator for each asset
                // create vol request first
                pathConfig->cdfMapping.resize(nbAssets);
                DateTimeArray  futDates(simTimeline);
                futDates.erase(futDates.begin());

                VolRequestLNStrikeSP volRequest(new LinearStrikeVolRequest(
                                               100.0, // any strike will do
                                               today,
                                               futDates.back(),
                                               false));

                PDFRequestLNStrike pdfLNStrkRequest(volRequest.get());

                for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                    // get market pdf calculator for each asset
                    PDFDefaultLNStrikeSP pdfMarket(
                                          new PDFDefaultLNStrike(today,
															&(mAsset->getAsset(iAsset)),
															&pdfLNStrkRequest));

                    // get vol, manage memory
                    CVolProcessedSP vol(vols[iAsset]->getProcessedVol(volRequest.get(),
                                                                   &(mAsset->getAsset(iAsset))));
                    // get model pdf calculator for each asset
                    PDFCalculatorSP pdfModel(
										PDFCalculatorMaker::makePDFCalculator("PDFFourier",
																		 today,
																		 pathConfig->cdfMappingParam->getModel().get(),
																		 &(mAsset->getAsset(iAsset)),
																		 vol.get()));

                    pathConfig->cdfMapping[iAsset] = CDFMappingSP(
                            new CDFMapping(futDates,
                                pathConfig->cdfMappingParam.get()));
                    // creates model and market spots grids
                    pathConfig->cdfMapping[iAsset]->setSpotsGrids(futDates,
                                                 //CAssetConstSP(&(mAsset->getAsset(iAsset))),
                                                 pathConfig->cdfMappingParam->getAsset(
                                                        mAsset->getAsset(iAsset).getName()),
                                                 volRequest,
                                                 today,
                                                 today,
                                                 pdfMarket,
                                                 pdfModel);

					// build the spline interpolant from model spots to market spots
					// cdfMapping[iAsset]->buildSplineInterpolant();

					// build the linear interpolant from model spots to market spots
					pathConfig->cdfMapping[iAsset]->buildLinearInterpolant();
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
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

    /** Computes the jump correlation matrix used to correlate the marginal jump participation probabilities */
    CDoubleMatrixSP computeJumpCorrelations() const {
        DoubleArraySP jumpParticipationProbs = computeJumpParticipationProbs();
        CDoubleMatrixConstSP correlations = mAsset->factorsCorrelationMatrix();
        CDoubleMatrixSP jumpCorrelations = CDoubleMatrixSP(new CDoubleMatrix(nbAssets, nbAssets));

        // computes jump correlations
        // this matrix is used in the Gaussian copula meant to correlate participation to market factor jumps
        int iAsset;
        for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
        {
            int jAsset;
            for (jAsset = iAsset ; jAsset < nbAssets ; jAsset++)
            {
                // iAsset = jAsset, correlation = 1.0
                if (iAsset == jAsset)
                {
                    (*jumpCorrelations)[iAsset][jAsset] = 1.0;
                }
                // if iAsset != jAsset, compute jump correlation based on asset correlation
                else
                {
                    JumpCorrelationData JumpCorrelationStruct;

                    // jump intensities for assets i anf j
                    double intensityIAsset = vols[iAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE);
                    double intensityJAsset = vols[jAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE);

                    // set lowJumpIntensityAsset and highJumpIntensityAsset
                    if (intensityIAsset <= intensityJAsset)
                    {
                        // lower jump marginal probability
                        JumpCorrelationStruct.lowJumpProba = (*jumpParticipationProbs)[iAsset];
                        // lower jump marginal intensity
                        JumpCorrelationStruct.lowJumpIntensity = intensityIAsset;
                        // higher jump marginal probability
                        JumpCorrelationStruct.highJumpProba = (*jumpParticipationProbs)[jAsset];
                        // higher jump marginal intensity
                        JumpCorrelationStruct.highJumpIntensity = intensityJAsset;
                    }
                    else
                    {
                        // lower jump marginal probability
                        JumpCorrelationStruct.lowJumpProba = (*jumpParticipationProbs)[jAsset];
                        // lower jump marginal intensity
                        JumpCorrelationStruct.lowJumpIntensity = intensityJAsset;
                        // higher jump marginal probability
                        JumpCorrelationStruct.highJumpProba = (*jumpParticipationProbs)[iAsset];
                        // higher jump marginal intensity
                        JumpCorrelationStruct.highJumpIntensity = intensityIAsset;
                    }

                    // lowJumpProba quantile of N1
                    JumpCorrelationStruct.normInvLowProba = N1Inverse(JumpCorrelationStruct.lowJumpProba);
                    // highJumpProba quantile of N1
                    JumpCorrelationStruct.normInvHighProba = N1Inverse(JumpCorrelationStruct.highJumpProba);
                    // correlation between the corresponding asset
                    JumpCorrelationStruct.correlation = (*correlations)[iAsset][jAsset];

                    // target conditional jump probability participation
                    double ratio = JumpCorrelationStruct.highJumpIntensity / JumpCorrelationStruct.lowJumpIntensity;
                    double correlation = JumpCorrelationStruct.correlation;
                    JumpCorrelationStruct.targetProba = Maths::min(Maths::max(ratio * correlation, 0.0), 1.0);
                    // JumpCorrelationStruct.targetProba = Maths::max(correlation, 0.0);

                    // jump correlation is equal to 1.0 if target probability is 1.0
                    if (Maths::isZero(JumpCorrelationStruct.targetProba - 1.0))
                    {
                        (*jumpCorrelations)[iAsset][jAsset] = 1.0;
                    }
                    // jump correlation is equal to 1.0 if target probability is -1.0
                    else if (Maths::isZero(JumpCorrelationStruct.targetProba + 1.0))
                    {
                        (*jumpCorrelations)[iAsset][jAsset] = -1.0;
                    }
                    else
                    {
                        (*jumpCorrelations)[iAsset][jAsset] = zbrentUseful(&objectivFunc,          /* (I) The function to find the root of */
                                                                           &JumpCorrelationStruct, /* (I) Parameter block */
                                                                           -1.0,                   /* (I) Low value for x */
                                                                           1.0,                    /* (I) High value for x */
                                                                           0.001);                 /* (I) Tolerance (0.1% correlation) */
                    }
                    (*jumpCorrelations)[iAsset][jAsset] = max(-0.99, min(0.99, (*jumpCorrelations)[iAsset][jAsset]));
                    (*jumpCorrelations)[jAsset][iAsset] = (*jumpCorrelations)[iAsset][jAsset];
                }
            }
        }

        // transform jumpCorrelations into a true correlation matrix if needed
        double squareErr;
        CDoubleMatrix correl(jumpCorrelations->symmToCorrel(&squareErr, 0.001));
        CDoubleMatrixSP jumpCorrs(copy(&correl)); // eigenValueFloor set to 0.0
        return jumpCorrs;
    }

    /** Computes the jump participation marginal probabilities */
    DoubleArraySP computeJumpParticipationProbs() const {
        DoubleArraySP jumpParticipationProbs = DoubleArraySP(new DoubleArray(nbAssets));

        // market factor jump intensity (sum of asset jump intensities)
        double jumpRateMF = 0.0;
        DoubleArray jumpRates(nbAssets);

        int iAsset;
        for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
        {
            // intensities adjusted to account for different trading time conventions between assets
            jumpRates[iAsset] = vols[iAsset]->getSVCJParam(VolSVCJ::COMMON_CRASH_RATE) *
                (tradYrs[iAsset].back() / tradYrs[0].back());
            jumpRateMF += jumpRates[iAsset];
        }

        // marginal jump participation marginal probabilities
        for (iAsset = 0 ; iAsset < nbAssets ; iAsset++)
        {
            (*jumpParticipationProbs)[iAsset] = jumpRates[iAsset] / jumpRateMF;
        }

        return jumpParticipationProbs;
    }

    bool hasPast() const {
        return simHasPast;
    }

    /** MCPathGen method */
    void generatePath(int pathIdx) {
        // Draw random numbers
        drawRandomNumbers(pathIdx);
        // loop over assets
        for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
            generatePath(pathIdx, iAsset);
        }
    }

    bool doingPast() const {
        return false;
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigSVCJ::PathGenFineProcess::getGaussData");
        try {
            return mAsset->factorsCorrelationMatrix();
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

private:

    /** Configures pathGen for antithetics */
    void configureAntithetics() {
    };

    /** Configures pathGen for nonAntithetics. We need a better name */
    void configureNonAntithetics() {
    };

    /** Draws random numbers */
    virtual void drawRandomNumbers(int pathIdx){
        // Draw random numbers for paths and jump sizes
        randomGenSpot->generate(pathIdx);
        randomGenVar->generate(pathIdx);
        randomGenUniJumpTimesMF->generate(pathIdx);
        if (nbAssets > 1)
        {
            randomGenUniJump->generate(pathIdx);
        }
        randomGenNormalJumpVar->generate(pathIdx);
        randomGenNormalJumpSpot->generate(pathIdx);
    }

    virtual void generatePath(int pathIdx, int iAsset) {
        static const string method = "MCPathConfigSVCJ::Generator::generatePath";
        try{
            const DoubleMatrix& randomSpot = randomGenSpot->getRandomNumbers();
            const DoubleMatrix& randomVol = randomGenVar->getRandomNumbers();
            const DoubleMatrix& randomUniJumpTimesMFMatrix = randomGenUniJumpTimesMF->getRandomNumbers();
            const DoubleMatrix& randomNormalJumpVar = randomGenNormalJumpVar->getRandomNumbers();
            const DoubleMatrix& randomNormalJumpSpot = randomGenNormalJumpSpot->getRandomNumbers();

            double randomInitVar = 0.0;
            if (vols[iAsset]->isRandomInitialVolatility()){
                randomGenInitVol[iAsset]->fetch(1, &randomInitVar);
            }

            int varDSType = pathConfig->varDiscreteSchemeIndex;
            int spotDSType = pathConfig->spotDiscreteSchemeIndex;
            int numJumps = 0;

            DoubleMatrix randomJumpTimes(1,maxNumJumps);// 1 column, maxNumJumps rows
             // simulate jump arrival times
            if (maxNumJumps >0){
                int iJump;
                for (iJump = 0; iJump < maxNumJumps; ++iJump)
                    randomJumpTimes[0][iJump] = randomUniJumpTimesMFMatrix[0][iJump];
                // arrival times
                dependenceJumpTimes[0]->correlateSeries(randomJumpTimes, pathIdx);
                // truncate jump dates to last sim date
                for (iJump = 0; iJump < maxNumJumps; ++iJump){
                    if (randomJumpTimes[0][iJump] > tradYrs[0].back()){
                        numJumps = iJump;
                        break;
                    }
                }
            }

            BoolArray jumpParticipation(maxNumJumps); // maxNumJumps elements
            if (nbAssets > 1)
            {
                const DoubleMatrix& randomUniJumpMatrix = randomGenUniJump->getRandomNumbers();
                int iJump;
                for (iJump = 0 ; iJump < maxNumJumps ; ++iJump)
                {
                    double uniformParticipation = Maths::max(N1(randomUniJumpMatrix[iAsset][iJump]), 0.00000001);
                    jumpParticipation[iJump] =
                        (uniformParticipation < (*jumpParticipationProbs)[iAsset]) ? true : false;
                }
            }
            else
            {
                int iJump;
                for (iJump = 0 ; iJump < maxNumJumps ; ++iJump)
                {
                    jumpParticipation[iJump] = true;
                }
            }

            // insert jumps: jumpTypesAdded and
            vector<VolSVCJ::STEP_TYPE_ITER> jumpTypesAdded;
            // jumpsAdded[0] -> tradYrs, [1] -> instVars, [2] -> integratedVars, [3] -> logSpots */
            vector<vector<VolSVCJ::LIST_ITER> >jumpsAdded(4);
            vector<double> varJumps; // keeps values of instVar jumps
            vols[iAsset]->addJumps(tradYrs[iAsset],
                                   randomJumpTimes[0],
                                   numJumps,
                                   jumpParticipation,
                                   randomNormalJumpVar[iAsset],
                                   randomNormalJumpSpot[iAsset],
                                   stepTypes,
                                   instVars[iAsset],
                                   integratedVars[iAsset],
                                   logSpots,
                                   varJumps,
                                   jumpTypesAdded,
                                   jumpsAdded,
                                   randomInitVar);

            // generate the var paths
            if (varDSType == VolSVJ_VarDSTypeHelper::EnumList::EULER)
                vols[iAsset]->generateVarPathsEuler(tradYrs[iAsset].begin(),
                                                    tradYrs[iAsset].size(),
                                                    randomVol[iAsset],
                                                    instVars[iAsset].begin(),
                                                    integratedVars[iAsset].begin());
            else if (varDSType == VolSVJ_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER)
                vols[iAsset]->generateVarPathsTransEuler(tradYrs[iAsset].begin(),
                                                         tradYrs[iAsset].size(),
                                                         randomVol[iAsset],
                                                         instVars[iAsset].begin(),
                                                         integratedVars[iAsset].begin());
            else
                throw ModelException("VolSVCJ::MCGenerateVarPaths",
                                     "unexpected variance discretization scheme of type "
                                     + VolSVJ_VarDSTypeHelper::getName(varDSType));

            // then generate the spot path
            if (spotDSType == VolSVJ_SpotDSTypeHelper::EnumList::EULER)
                vols[iAsset]->generateSpotPathsEuler(tradYrs[iAsset].begin(),
                                                     tradYrs[iAsset].size(),
                                                     stepTypes.begin(),
                                                     randomSpot[iAsset],
                                                     randomVol[iAsset],
                                                     logSpots.begin(),
                                                     instVars[iAsset].begin(),
                                                     integratedVars[iAsset].begin());
            else if (spotDSType == VolSVJ_SpotDSTypeHelper::EnumList::EXACT)
                vols[iAsset]->generateSpotPathsExact(tradYrs[iAsset].begin(),
                                                     tradYrs[iAsset].size(),
                                                     stepTypes.begin(),
                                                     randomSpot[iAsset],
                                                     varJumps,
                                                     logSpots.begin(),
                                                     instVars[iAsset].begin(),
                                                     integratedVars[iAsset].begin());
            else
                throw ModelException("VolSVCJ::MCGenerateExJumpPaths",
                                     "unexpected spot discretization scheme of type "
                                     + VolSVJ_SpotDSTypeHelper::getName(spotDSType));

            // remove inserted jump points
            if (numJumps>0){
                int nStepsAdded = jumpsAdded[0].size();

                for (int j = 0; j < nStepsAdded; j++){
                    stepTypes.erase(jumpTypesAdded[j]);
                    tradYrs[iAsset].erase(jumpsAdded[0][j]);
                    instVars[iAsset].erase(jumpsAdded[1][j]);
                    integratedVars[iAsset].erase(jumpsAdded[2][j]);
                    logSpots.erase(jumpsAdded[3][j]);
                }
            }

			// Compute Quanto adjustment
			DoubleArray adjQuanto(logSpots.size(), 0);
			if (isQuanto[iAsset]) {
				int nbSimDates = stepTypes.size();
				VolSVCJ::LIST_ITER instV = instVars[iAsset].begin();
				VolSVCJ::STEP_TYPE_ITER stepTypeIter = stepTypes.begin();
				int idxAdj = 0;
				instV++;
				stepTypeIter++;
				for (int iSimDates = 1; iSimDates < nbSimDates; iSimDates++, instV++, stepTypeIter++) {
					adjQuanto[idxAdj] -= eqFXCorr[iAsset] * fxVar[iAsset][iSimDates - 1] * sqrt(Maths::max((*instV),0.0));
					if (*stepTypeIter == VolSVCJ::SAMPLE_DATE) {
						idxAdj++;
						adjQuanto[idxAdj] = adjQuanto[idxAdj - 1];
					}
				}
			}

            // then, populate spot path
            const DoubleArray& fwds = this->fwds[iAsset];
            DoubleArray& spots = productPaths[iAsset];
            int nbFutFixDates = logSpots.size();
            // NB: spots[timeline->numPastDates] corresponds to futureDates[1], i.e.
            // the first future product/fixing date -- not today
            VolSVCJ::LIST_ITER logS = logSpots.begin();
            logS++;
            for (int iFutFixDates = 1, iFixDates = spotProdOffset[iAsset];
                 iFutFixDates < nbFutFixDates; ++iFutFixDates, ++iFixDates, logS++) {
                    spots[iFixDates] = fwds[iFutFixDates] * exp(*logS) * exp(adjQuanto[iFutFixDates - 1]);
            }

            // do Implied Mapping if requested
            if (!!pathConfig->cdfMappingParam && numFutSteps>0){
                double* s = &*(spots.begin() + spotProdOffset[iAsset]);

                pathConfig->cdfMapping[iAsset]->impliedMCSpot(s,
                                                            s,
                                                            0,
                                                            numFutSteps);
            }

            // followed by quad var path
            DoubleArray& quadVars = quadVarProdPaths[iAsset];
            VolSVCJ::LIST_ITER instV = instVars[iAsset].begin();
            VolSVCJ::LIST_ITER integratedV = integratedVars[iAsset].begin();
            // copy everything, lasy way, may be a bit slow
            vector<double> instVar(instVars[iAsset].size());
            vector<double> integratedVar(integratedVars[iAsset].size());
            int i;
            for (i=0; i<(int)instVar.size(); i++, instV++)
                instVar[i] = *instV;

            for (i=0; i<(int)integratedVar.size(); i++, integratedV++)
                integratedVar[i] = *integratedV;

            if (quadVars.size()){
                for (int iFutFixDate = 1, iQuadVarProdDate = quadVarProdOffset[iAsset];
                     iFutFixDate < nbFutFixDates; ++iFutFixDate, ++iQuadVarProdDate) {
                        int iIntegratedVar = fixDateOffsets[iFutFixDate - 1];
                        quadVars[iQuadVarProdDate] = integratedVar[iIntegratedVar];
                }
            }
            // and sqrt annualized quad var path
            DoubleArray& sqrtAnnualQuadVars = sqrtAnnualQuadVarProdPaths[iAsset];

            if (sqrtAnnualQuadVars.size()){
                bool hasPast = sqrtAnnualQuadVarProdOffset[iAsset] > 0;
                for (int iFutFixDate = 1, iSqrtAnnualQuadVarProdDate = sqrtAnnualQuadVarProdOffset[iAsset],
                     iIntegratedVar = fixDateOffsets[0];
                     iFutFixDate < nbFutFixDates;
                     iIntegratedVar += fixDateOffsets[iFutFixDate++],
                     ++iSqrtAnnualQuadVarProdDate) {
                    if (hasPast){
                        // future annualized var from today till current sim date
                        double tradYear = tradYears[iAsset][iIntegratedVar];
                        double futAnnualQuadVar = Maths::isZero(tradYear) ?
                                                  instVar[iIntegratedVar] :
                                                  integratedVar[iIntegratedVar] / tradYear;
                        // total annualized var from first sampling date till current sim date
                        double annualQuadVar = pastWeights[iAsset] * pastAnnualQuadVars[iAsset]
                                               + (1.0 - pastWeights[iAsset]) * futAnnualQuadVar;
                        // total annualized vol
                        sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate]
                            = sqrt(annualQuadVar);
                    }
                    else{
                        int integratedVarStartIdx = fixDateOffsets[0];
                        if (iIntegratedVar == integratedVarStartIdx){
                            sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate] = 0.0;
                        }
                        else{
                            // quad var increment between first sampling date and current sim date
                            double quadVarInc = integratedVar[iIntegratedVar]
                                                - integratedVar[integratedVarStartIdx];
                            // time elapsed since first sampling date
                            double timeDiff = tradYears[iAsset][iIntegratedVar]
                                              - tradYears[iAsset][integratedVarStartIdx];
                            // total annualized vol
                            sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate]
                                = sqrt(quadVarInc / timeDiff);
                        }
                    }
                }
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** pre-compute fwds arrays, get quanto data here for now */
    void preComputeFwds() {
        // for each asset
        for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
            // need to cache the forward at all simulation dates
            fwds[iAsset].resize(numFutSteps + 1);
            // current way of dealing with quanto
            CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));

            if (ProtAsset::TYPE->isInstance(asset)) {
                throw ModelException("preComputeFwds", "Protected assets (" + asset->getTrueName() +
                                     ") of type " + asset->getClass()->getName() +
                                     " are not supported yet");
            }

            isQuanto[iAsset] = ProtEquity::TYPE->isInstance(asset);
            if (isQuanto[iAsset]) {
                fxVar[iAsset] = DoubleArray(simDates->size()-1);
                ProtEquityConstSP eq = ProtEquityConstSP::dynamicCast(asset);
                // compute plain fwds
                CAssetConstSP plainAsset = eq->getPlainAsset();
                plainAsset->fwdValue(simTimeline, fwds[iAsset]);
                // get fx vol and corr
                FXVolBaseWrapper fxVol = eq->getFXVol();
                eqFXCorr[iAsset] = eq->getCorrelation()->getCorrelation();;
                // interpolate the fx vols atm
                ATMVolRequestSP fxVolRequest(new ATMVolRequest());
                // interpolate the vol
                const Asset* fx = eq->getFXAsset();
                CVolProcessedSP  fxVoltmp(fxVol->getProcessedVol(fxVolRequest.get(), fx));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(fxVoltmp);
                volFXBS->CalcVol(*simDates, CVolProcessedBS::forward, fxVar[iAsset]);
                // calculate fx vars
				const TimeMetric& metric = vols[iAsset]->getTimeMetric();
                for (int iStep = 0; iStep < fxVar[iAsset].size(); iStep++) {
                    double dt = metric.yearFrac((*simDates)[iStep], (*simDates)[iStep+1]);
                    fxVar[iAsset][iStep] *= dt; // using asset dt
                }
            }
            else {
                mAsset->factorFwdValues(iAsset, simTimeline, fwds[iAsset]);
            }
        }
    }

    // fields
    MCPathConfigSVCJ*       pathConfig;
    const IMultiFactors*    mAsset;            // multi asset interface
    int                     nbAssets;          // number of assets
    DateTimeArraySP         simDates;          // fine grid of simulation dates
    list<VolSVCJ::STEP_TYPE>stepTypes;         // keeps type of each sim date
    vector<DoubleArray>     fwds;              // fwds [iAsset][iStep]
    VolSVCJArray            vols;              // array of VolSVCJs
    TimeMetricArray         timeMetrics;       // array of time metrics
    vector<DoubleArray>     productPaths;      // spot path [iAsset][iStep]
    vector<DoubleArray>     tradYears;         // fine grid of trad years
    IntArray                fixDateOffsets;    // offsets for fixing dates
    IMCRandomSP             randomGenSpot;     // Random number generator for spots
    IMCRandomSP             randomGenVar;      // Random number generator for vars
    IMCRandomSP             randomGenUniJumpTimesMF; // Random number generator for jump times uniform deviates
    IMCRandomSP             randomGenUniJump; // Random number generator for jump participations uniform deviates
    DoubleArraySP           jumpParticipationProbs; // marginal jump participation probability
    IMCRandomSP             randomGenNormalJumpVar;  // Random number generator for jump normal deviates
    IMCRandomSP             randomGenNormalJumpSpot;  // Random number generator for jump normal deviates
    bool                    simulateInitialVol;     // if init vol is simulated (stationary distribution)
    vector<MCRandGammaSP>   randomGenInitVol;  // Random number generator for init vol - stationary gamma distribution
    DependenceSP            dependence;     // Dependence object
    vector<double>          pastAnnualQuadVars; // past realized vars at value date [iAsset]
    vector<double>          pastWeights;       // past weights at value date [iAsset]
    vector<double>          maxDrifts;         // XXX product of MAX(1, drifts)

    int                     maxNumJumps;
    vector<list<double> >   tradYrs;         // fine grid of trad years
    list<double>            logSpots;          // spot path [iAsset][iStep]
    vector<list<double> >   instVars;          // instantaneous variance path [iAsset][iStep]
    vector<list<double> >   integratedVars;    // integrated variance path [iAsset][iStep]
    vector<PoissonSP>       dependenceJumpTimes;

    BoolArray               isQuanto;
    DoubleArray             eqFXCorr;
    vector<DoubleArray>     fxVar;

    // Timeline fields
    DateTimeArray           simTimeline;     //!< Today + strictly future merged dates
    int                     numFutSteps;     //!< Number of strictly future merged dates
    bool                    simHasPast;      //!< Whether product has past

    vector<DoubleArray>     spotProdPaths;      // spot path [iAsset][iStep]
    vector<int>             spotProdOffset;  //!< Starting point for future path per asset

    vector<DoubleArray>     quadVarProdPaths;      // quadVar path [iAsset][iStep]
    vector<int>             quadVarProdOffset;  //!< Starting point for future path per asset

    vector<DoubleArray>     sqrtAnnualQuadVarProdPaths;      // sqrtAnnualQuadVar path [iAsset][iStep]
    vector<int>             sqrtAnnualQuadVarProdOffset;  //!< Starting point for future path per asset
};


/** Referee class that distributes simulation to components
    e.g. spot, hitTime etc. */
class MCPathConfigSVCJ::Gen: virtual public MCPathGenerator, // For backward compatibility
                           virtual public MCPathGen,
                           virtual public IStateVariableGen::IStateGen {
public:
    /** Constructor */
    Gen(int                      numSimPaths,
        MCPathConfigSVCJ*          pathConfig,
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient):
    nowPathIdx(0) {
        static const string routine = "MCPathConfigSVCJ::Gen::Gen";

        try {
            // Collect state variables from product and categorize them
            StateVariableCollectorSP svCollector(new StateVariableCollector());
            prodClient->collectStateVars(svCollector);
            IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

            // Spot requests
            SVGenSpotArray spotGenArray = filterStateVars<SVGenSpot>(stateVarGenArray);

            // QuadVar requests
            MCQuadVarArray quadVarGenArray = filterStateVars<MCQuadVar>(stateVarGenArray);

            // SqrtAnnualQuadVar requests
            MCSqrtAnnualQuadVarArray sqrtAnnualQuadVarGenArray = filterStateVars<MCSqrtAnnualQuadVar>(stateVarGenArray);

            // Create discount factor state variables
            SVGenDiscFactorArray discFactors = filterStateVars<SVGenDiscFactor>(stateVarGenArray);
            for(unsigned int iVar = 0; iVar < discFactors.size(); iVar++) {
                IStateVariableSP sv(discFactors[iVar]->determinsticSV(false /* not doing past */));
                svDBase.append(discFactors[iVar], sv);
            }
            // Nothing else is supported
            if(stateVarGenArray.size()) {
                throw ModelException("Unable to recognize all state variable types.");
            }

            // XXX should pass in a PastPathGen in the first place, really
            PastPathGen* pastPathGen = dynamic_cast<PastPathGen*>(pastPathGenerator.get());
            if(!pastPathGen) {
                throw ModelException("Past path generator is not of PastPathGen type.");
            }

            // Create a (spot, quad var, sqrt annualized quad var) path generator
            pathGenFineProcess = MCPathConfigSVCJ::PathGenFineProcessSP(
                new MCPathConfigSVCJ::PathGenFineProcess(numSimPaths,
                                                       pathConfig,
                                                       pastPathGen->getPathGenSpot(),
                                                       prodClient,
                                                       spotGenArray,
                                                       quadVarGenArray,
                                                       sqrtAnnualQuadVarGenArray,
                                                       svDBase));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** MCpathGenerator methods */
    // Deprecated methods
    int NbSimAssets() const {
        static const string routine = "MCPathConfigSVCJ::Gen::NbSimAssets";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    const double* Path(int iAsset, int iPath) const {
        static const string routine = "MCPathConfigSVCJ::Gen::Path";
        throw ModelException(routine, "Method is retired for StateVars");
    };

    double refLevel(int iAsset, int iPath) const {
        static const string routine = "MCPathConfigSVCJ::Gen::refLevel";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    double maxDriftProduct(int iAsset) const {
        static const string routine = "MCPathConfigSVCJ::Gen::maxDriftProduct";
        try {
            return pathGenFineProcess->maxDriftProduct(iAsset);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    int begin(int iAsset) const {
        static const string routine = "MCPathConfigSVCJ::Gen::begin";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    int end(int iAsset) const{
        static const string routine = "MCPathConfigSVCJ::Gen::end";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    // Live methods
    bool hasPast() const {
        return pathGenFineProcess->hasPast();
    }

    bool doingPast() const {
        return false;
    }

    void generatePath(int pathIdx) {
        nowPathIdx = pathIdx;
        pathGenFineProcess->generatePath(pathIdx);
    }

    int getPathIndex() const {
        return nowPathIdx;
    }

    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    IStateVariableSP create(const IStateVariableGen* svGen) {
        static const string routine = "MCPathConfigSVCJ::Gen::create";

        try {
            return svDBase.find(svGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    StateVarDBase                        svDBase;               //!< Collection of Generators + statevars
    MCPathConfigSVCJ::PathGenFineProcessSP pathGenFineProcess;    //!< Spot generator
    int                                  nowPathIdx;            //!< Current path idx
};
#endif

/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigSVCJ::makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control,
        Results*                           results,
        DateTimeArray&                     mainSimDates ){
    static const string routine = "MCPathConfigSVCJ::makePathGenerator";

    if(cachingRequested) {
        throw ModelException(routine, "Paths caching is not supported in MCPathConfigSVCJ");
    }

    // create empty sim dates array if simDates is null or on new pricing run
    if (!simDates || control->isPricing()){
        simDates = DateTimeArraySP(new DateTimeArray());
    }

    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
    if(!prodClient){
        refCountPtr<MCPathConfigSVCJ::Generator> futurePathGenBase(
            new MCPathConfigSVCJ::Generator(numPaths,
                                          this,
                                          pastPathGenerator,
                                          prod));

		OutputRequest* request;
		int numAssets = futurePathGenBase->nbAssets;
		if	((control->isPricing()) &&
			(request = control->requestsOutput(OutputRequest::FWD_AT_MAT))) {
			for (int i = 0; i < numAssets ; ++i) {
				CAssetConstSP asset = CAssetConstSP::attachToRef(&(futurePathGenBase->mAsset)->getAsset(i));
				results->storeRequestResult(request, futurePathGenBase->fwds[i].back(),
					asset.get()->getName());
			}
		}

        // store num of sim dates
        if (control && control->isPricing()){
            int nbSimDates = simDates->size();
            results->storeScalarGreek(nbSimDates - 1,
                                      Results::DEBUG_PACKET,
                                      OutputNameSP(
                                          new OutputName("MertonLV_SIM_STEPS_USED")));
        }

        // Construct the MCPathGenerator
        return MCPathBase::createPathGenerator(pastPathGenerator,
                                               this,
                                               prod,
                                               futurePathGenBase);
    } else {
        // State variables approach
        return MCPathGeneratorSP(new MCPathConfigSVCJ::Gen(
            numPaths,
            this,
            pastPathGenerator,
            prodClient));
    }
}

/** Creates a past path generator */
MCPathGeneratorSP MCPathConfigSVCJ::pastPathGenerator(const IMCProduct* prod) {
    static const string method = "MCPathConfigSVCJ::pastPathGenerator";
    try{
        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
        if(!prodClient){
            // Old approach
            return MCPathConfig::pastPathGenerator(prod);
        } else {
            // State variables approach
            return MCPathGeneratorSP(new PastPathGen(prodClient));
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** MC LV model */
class MonteCarloSVCJDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;

private:
    MonteCarloSVCJDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloSVCJDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloSVCJ);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloSVCJ(){
        return new MonteCarloSVCJDefault();
    }
};

CClassConstSP const MonteCarloSVCJDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloSVCJDefault", typeid(MonteCarloSVCJDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigSVCJLoad(){
    return (MCPathConfigSVCJ::TYPE != 0 && MonteCarloSVCJDefault::TYPE !=0);
}


/** MC LV model */
class MCPathConfigSVCJImplied: public MCPathConfigSVCJ {
public:
    static CClassConstSP const TYPE;

private:
    MCPathConfigSVCJImplied():MCPathConfigSVCJ(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MCPathConfigSVCJImplied, clazz);
        SUPERCLASS(MCPathConfigSVCJ);
        EMPTY_SHELL_METHOD(defaultMCPathConfigSVCJImplied);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigSVCJImplied(){
        return new MCPathConfigSVCJImplied();
    }
};

CClassConstSP const MCPathConfigSVCJImplied::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigSVCJImplied", typeid(MCPathConfigSVCJImplied), load);


class MonteCarloSVCJImpliedDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;

private:
    MonteCarloSVCJImpliedDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloSVCJImpliedDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloSVCJImplied);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloSVCJImplied(){
        return new MonteCarloSVCJImpliedDefault();
    }
};

CClassConstSP const MonteCarloSVCJImpliedDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloSVCJImpliedDefault", typeid(MonteCarloSVCJImpliedDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigSVCJImpliedLoad(){
    return (MCPathConfigSVCJImplied::TYPE != 0 && MonteCarloSVCJImpliedDefault::TYPE !=0);
}



/************************/
/*** SVCJ + LV scheme ***/
/************************/

class MCPathConfigSVCJLV : public MCPathConfigSVCJ {
public:
    static CClassConstSP const TYPE;

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retreiving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const {
        return MarketDataFetcherSP(
            new MDF( VolSVCJLV::TYPE->getName(), cdfMappingParam ) );
    }

protected:
    // for reflection
    MCPathConfigSVCJLV(CClassConstSP clazz) :
        MCPathConfigSVCJ(clazz) {}

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MCPathConfigSVCJLV, clazz);
        SUPERCLASS(MCPathConfigSVCJ);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        EMPTY_SHELL_METHOD(defaultMCPathConfigSVCJLV);
    }

    static IObject* defaultMCPathConfigSVCJLV(){
        return new MCPathConfigSVCJLV();
    }

    // for reflection
    MCPathConfigSVCJLV() :
        MCPathConfigSVCJ(TYPE) {}
};

/** MC LV model */
class MonteCarloSVCJLVDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;

private:
    MonteCarloSVCJLVDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloSVCJDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloSVCJLV);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloSVCJLV(){
        return new MonteCarloSVCJLVDefault();
    }
};

CClassConstSP const MCPathConfigSVCJLV::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigSVCJLV", typeid(MCPathConfigSVCJLV), load);

CClassConstSP const MonteCarloSVCJLVDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloSVCJLVDefault", typeid(MonteCarloSVCJLVDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigSVCJLVLoad(){
    return (MCPathConfigSVCJLV::TYPE != 0 && MonteCarloSVCJLVDefault::TYPE !=0);
}

DRLIB_END_NAMESPACE
