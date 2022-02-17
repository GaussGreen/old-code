//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSV.cpp
//
//   Description : Monte Carlo path generator for VolSV
//
//   Date        : 04 March 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include <algorithm>  
#include "edginc/mathlib.hpp"
#include "edginc/MCPathConfigParametric.hpp"
#include "edginc/MCPathGenerator.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/Random.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/StruckAsset.hpp"
#include "edginc/DependenceGauss.hpp"

#define STATE_VARIABLES

#ifdef STATE_VARIABLES
#include "edginc/SVGenSpot.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MonteCarlo.hpp"
#include <algorithm>
#endif



DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////////
// MCPathConfig
///////////////////////////////////////////////////////////////////////////////////
class MCPathConfigSV: public MCPathConfigParametric {
public: 
    static CClassConstSP const TYPE;

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const {
        return MarketDataFetcherSP(new MDFAssetVol(VolSV::TYPE->getName(), 
                                                   VolSurface::TYPE->getName()));
    }

private: 
    virtual void validatePop2Object(){
        static const string routine("MCPathConfigSV::validatePop2Object");
        try {
            Maths::checkPositive(nbStepsPerYear, "nbStepsPerYear");
            spotDiscreteSchemeIndex = VolSV_SpotDSTypeHelper::getIndex(spotDiscreteScheme);
            varDiscreteSchemeIndex = VolSV_VarDSTypeHelper::getIndex(varDiscreteScheme);
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
        smartPtr<MCPathConfigSV> pathConfig(copy(this)); // copy to avoid const problems
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
            throw ModelException("MCPathConfigSV::storagePerPath", "Internal "
                                 "error - no simulation dates");
        }
        // we store both correlated and uncorrelated numbers but only for
        // every other path. Hence no times by 2.
        return (sizeof(double) * pathConfig->simDates->size() * 
                product->getNumAssets());
    }


    virtual bool vegaMatrixSupported() const {
        return true;
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
        static const string method("MCPathConfigSV::Generator::genSimDates");
        try{
            // calculate approx number of sim dates
            int nbFutureDates = futureDates.size();
            ASSERT(nbFutureDates > 1);
            DateTime startDate = futureDates.front();
            DateTime horizonDate = futureDates.back();
            double horizonYearFrac = timeMetric.yearFrac(startDate, horizonDate);
            int approxNbSteps = static_cast<int>(ceil(horizonYearFrac * nbStepsPerYear));
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
                int nbSteps = static_cast<int>(ceil(yearFrac * nbStepsPerYear));
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

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MCPathConfigSV, clazz);
        SUPERCLASS(MCPathConfigParametric);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        EMPTY_SHELL_METHOD(defaultMCPathConfigSV);
        FIELD(dependenceType, "not used");
        FIELD_MAKE_OPTIONAL(dependenceType);
        FIELD(nbStepsPerYear, "Nb of steps per year");
        FIELD_MAKE_OPTIONAL(nbStepsPerYear);
        FIELD(spotDiscreteScheme, "Spot discretization scheme (" 
                                         + VolSV_SpotDSTypeHelper::getNameList()
                                         + ")");
        FIELD_MAKE_OPTIONAL(spotDiscreteScheme);
        FIELD(spotDiscreteSchemeIndex, "");
        FIELD_MAKE_TRANSIENT(spotDiscreteSchemeIndex);
        FIELD(varDiscreteScheme, "Variance discretization scheme (" 
                                         + VolSV_VarDSTypeHelper::getNameList()
                                         + ")");
        FIELD_MAKE_OPTIONAL(varDiscreteScheme);
        FIELD(varDiscreteSchemeIndex, "");
        FIELD_MAKE_TRANSIENT(varDiscreteSchemeIndex);        
        FIELD_NO_DESC(simDates);
        FIELD_MAKE_TRANSIENT(simDates); 
        FIELD(isCarefulRandoms, "Use careful randoms if true; don't otherwise");
        FIELD_MAKE_OPTIONAL(isCarefulRandoms);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMCPathConfigSV(){
        return new MCPathConfigSV();
    }

    // for reflection
    MCPathConfigSV():
    MCPathConfigParametric(TYPE),
    isCarefulRandoms(false),
    dependenceType("not used"),
    nbStepsPerYear(12),
    spotDiscreteScheme(VolSV_SpotDSTypeHelper::getDefaultName()),
    varDiscreteScheme(VolSV_VarDSTypeHelper::getDefaultName()),
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
    string          spotDiscreteScheme;
    string          varDiscreteScheme;

    // transient fields
    int                spotDiscreteSchemeIndex;
    int                varDiscreteSchemeIndex;
    DateTimeArraySP    simDates; // cached between tweaks
};

CClassConstSP const MCPathConfigSV::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigSV", typeid(MCPathConfigSV), load);

///////////////////////////////////////////////////////////////////////////////////
// Generator
///////////////////////////////////////////////////////////////////////////////////
class MCPathConfigSV::Generator: public MCPathBase,
                                 virtual public IMCRandom::Callbacks,
                                 public DependenceMakerGauss::Support {
public:
    Generator(int                      numSimPaths,
              MCPathConfigSV*          pathConfig, 
              const MCPathGeneratorSP& pastPathGenerator,
              const IMCProduct*         prod):
    pathConfig(pathConfig),
    mAsset(prod->getMultiFactors()), 
    nbAssets(prod->getNumAssets()),
    simDates(pathConfig->simDates),
    fwds(nbAssets),
    quanto(nbAssets), //quanto parameters
    vols(nbAssets),
    timeMetrics(nbAssets),
    tradYears(nbAssets),
    productPaths(nbAssets),
    instVars(nbAssets),
    integratedVars(nbAssets),
    maxDrifts(nbAssets, 1.0),
    isCarefulRandoms(pathConfig->carefulRandoms()) {
        static const string method("MCPathConfigSV::Generator::Generator");
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
            refData = getRefData(timeline, nbPaths, mAsset, 
                                 prod, pastPathGenerator);

            // Initialize vols and time metrics
            int iAsset;
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                // copy the vols that are in the multiasset
                vols[iAsset] = VolSVSP::dynamicCast(
                    VolRequestRaw::copyVolBase(*mAsset, iAsset));
                // copy time metric
                timeMetrics[iAsset] = TimeMetricSP(copy(&vols[iAsset]->getTimeMetric()));
            }

            // reuse sim dates between tweaks for numerical stability
            if (simDates->empty()){
                simDates = genSimDates(timeline->futureDates,
                                       *timeMetrics[0],       // pick a time metric !
                                       pathConfig->nbStepsPerYear,
                                       fixDateOffsets);
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
                tradYears[iAsset].resize(simDates->size());
                for (int iDate = 0; iDate < simDates->size(); ++iDate) {
                    tradYears[iAsset][iDate] 
                        = timeMetrics[iAsset]->yearFrac(startDate,
                                                        (*simDates)[iDate]);
                }
            }

            // pre compute fwds
            preComputeFwds();

		    // set up Dependence
            dependence = pathConfig->dependenceMaker->createDependence(this);

            // Initialize random number generator for spots
            // NB the size of the random deviates depends on the scheme 
            // that is requested
            int nbSpotDeviates;
            if (pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EULER){
                nbSpotDeviates = simDates->size() - 1;
            }
            else if(pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EXACT){
                nbSpotDeviates = logSpots.size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSV_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
            }
            int nbVarDeviates;
            if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSV_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            randomGenSpot = IMCRandomSP(new MCRandomNoCache(
                this,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbSpotDeviates,
                nbAssets,
                timeline->totalNumPastDates));
            randomGenVar = IMCRandomSP(new MCRandomNoCache(
                this,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbVarDeviates,
                nbAssets,
                0));
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

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigSV::Generator::getGaussData");
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
    }
    
    virtual void generatePath(int pathIdx, int iAsset, int iPath) {
        static const string method = "MCPathConfigSV::Generator::generatePath";
        try{
            const DoubleMatrix& randomSpot = randomGenSpot->getRandomNumbers();
            const DoubleMatrix& randomVol = randomGenVar->getRandomNumbers();    
            // XXX this list of nested if statements is obviously not viable
            if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::EULER){
                if (pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EULER){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EULER>(),
                                                  Int2Type<VolSV_VarDSType::EULER>(),
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else if(pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EXACT){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EXACT>(),
                                                  Int2Type<VolSV_VarDSType::EULER>(), 
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else{
                    throw ModelException(method, 
                                         "unexpected spot discretization scheme of type "
                                         + VolSV_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
                }
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                if (pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EULER){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EULER>(),
                                                  Int2Type<VolSV_VarDSType::VAR_TRANSFORM_EULER>(),
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else if(pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EXACT){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EXACT>(),
                                                  Int2Type<VolSV_VarDSType::VAR_TRANSFORM_EULER>(), 
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else{
                    throw ModelException(method, 
                                         "unexpected spot discretization scheme of type "
                                         + VolSV_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
                }
            }
            else{
                throw ModelException(method, 
                                     "unexpected variance discretization scheme of type "
                                     + VolSV_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            // then, populate spot path
            const DoubleArray& fwds = this->fwds[iAsset];
            DoubleArray& spots = productPaths[iAsset];
            int nbFutFixDates = logSpots.size();
            // NB: spots[timeline->numPastDates] corresponds to futureDates[1], i.e.
            // the first future product/fixing date -- not today
            for (int iFutFixDates = 1, iFixDates = timeline->numPastDates;
                 iFutFixDates < nbFutFixDates; ++iFutFixDates, ++iFixDates) {
                    // insert here drift adjustment for currency protection
                    spots[iFixDates] = fwds[iFutFixDates] * exp(logSpots[iFutFixDates]);
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
    // Reminder about timelines:
    // No state variables: coarse path (timeline->futureDates); fine path (simDates);
    // State variables: coarse path (simTimeline); fine path (simDates);
    void preComputeFwds() {

        static const string routine = "MCPathConfigSV::preComputeFwds";

        try {

            // for each asset
            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                // need to cache the forward at all simulation dates
                fwds[iAsset].resize(timeline->numFutSteps + 1);
                
                CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));

                if (CAsset::IStruck::TYPE->isInstance(asset)) {
                    throw ModelException(routine, "Struck assets (" + asset->getTrueName() + 
                                        ") are not supported yet");
                } 

                //using quanto of type QuantoParamSV
                bool isQuanto = Asset::IQuanto::TYPE->isInstance(asset);
                DoubleArray  fxSqrtVar = DoubleArray(simDates->size() - 1);
                double eqFXCorr = 0.0;
                if (isQuanto) {
                    const Asset::IQuanto& prot = dynamic_cast<const Asset::IQuanto&>(*asset.get());

                    // Get the unprotected forward values
                    prot.unadjustedFwdValue(timeline->futureDates, fwds[iAsset]);
                    
                    // get fx corr
                    eqFXCorr = prot.getCorrelation()->getCorrelation();
                    
                    // get fx vols atm
                    ATMVolRequestSP fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP volFX(prot.getProcessedFXVol( fxVolRequest.get() ));
                    // cast to the type of vol we're expecting
                    CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(volFX);
                    volFXBS->CalcVol(*simDates, CVolProcessedBS::forward, fxSqrtVar);
                    
                    // calculate fx vars
                    for (int iStep = 0; iStep < fxSqrtVar.size(); iStep++) {
                        double dt = tradYears[iAsset][iStep+1]-tradYears[iAsset][iStep];
                        fxSqrtVar[iStep] *= dt; // using asset dt
                    }
                    quanto[iAsset] = QuantoParamSVSP(new QuantoParamSV(isQuanto,eqFXCorr,fxSqrtVar));
		        } else{
                    quanto[iAsset] = QuantoParamSVSP(new QuantoParamSV(isQuanto,eqFXCorr,fxSqrtVar));
                    mAsset->factorFwdValues(iAsset, timeline->futureDates, fwds[iAsset]);
                }
                
            }
        } catch (exception& e) {
        throw ModelException(e, routine);
        }
    }


    // fields
    MCPathConfigSV*         pathConfig;
    const IMultiFactors*    mAsset;            // multi asset interface
    int                     nbAssets;          // number of assets
    DateTimeArraySP         simDates;          // fine grid of simulation dates
    vector<DoubleArray>     fwds;              // fwds [iAsset][iStep]
    VolSVArray              vols;              // array of VolSVs
    TimeMetricArray         timeMetrics;       // array of time metrics
    vector<DoubleArray>     tradYears;         // fine grid of trad years
    vector<DoubleArray>     productPaths;      // spot path [iAsset][iStep]
    DoubleArray             logSpots;          // spot path [iAsset][iStep]
    vector<DoubleArray>     instVars;          // instantaneous variance path [iAsset][iStep]
    vector<DoubleArray>     integratedVars;    // integrated variance path [iAsset][iStep]
    IntArray                fixDateOffsets;    // offsets for fixing dates
    IMCRandomSP             randomGenSpot;     // Random number generator for spots
    IMCRandomSP             randomGenVar;      // Random number generator for vars
    DependenceSP            dependence;        // Dependence object
    MCProductTimelineSP     timeline;          // Product timeline
    RefLevelDataSP          refData;           // RefLevel data
    DoubleArray             maxDrifts;         // XXX product of MAX(1, drifts)
    bool                    isCarefulRandoms;  // XXX
    vector<QuantoParamSVSP>     quanto;            // quanto parameters
};

#ifdef STATE_VARIABLES
/** Spot path generator using state variable approach */
class MCPathConfigSV::PathGenFineProcess: virtual public MCPathGen,
                                          virtual public IMCRandom::Callbacks,
                                          public DependenceMakerGauss::Support {
public:
    PathGenFineProcess(int                             numSimPaths, 
                       MCPathConfigSV*                 pathConfig, 
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
    quanto(nbAssets),
    vols(nbAssets),
    timeMetrics(nbAssets),
    tradYears(nbAssets),
    instVars(nbAssets),
    integratedVars(nbAssets),
    maxDrifts(nbAssets, 1.0),
    simHasPast(pastPathGenSpot->hasPast()),
    spotProdPaths(nbAssets),
    spotProdOffset(nbAssets, 0),
    quadVarProdPaths(nbAssets),
    quadVarProdOffset(nbAssets, 0),
    sqrtAnnualQuadVarProdPaths(nbAssets),
    sqrtAnnualQuadVarProdOffset(nbAssets, 0) {

        static const string method("MCPathConfigSV::Generator::Generator");
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

            // Initialize vols and time metrics
            int iAsset;
            for (iAsset = 0; iAsset < nbAssets; iAsset++) {
                // copy the vols that are in the multiasset
                vols[iAsset] = VolSVSP::dynamicCast(
                    VolRequestRaw::copyVolBase(*mAsset, iAsset));
                // copy time metric
                timeMetrics[iAsset] = TimeMetricSP(copy(&vols[iAsset]->getTimeMetric()));
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
                MCPath::IStateVarArray spotSVArray(MCPath::createPaths(
                    false,
                    spotGenArray,
                    spotDatesPerAsset,
                    spotBeginInd,
                    spotEndInd,
                    spotPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < spotGenArray.size(); iVar++) {
                    svDBase.append(spotGenArray[iVar], spotSVArray[iVar]);
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
                MCPath::IStateVarArray quadVarSVArray(MCPath::createPaths(
                    false,
                    quadVarGenArray,
                    quadVarDatesPerAsset,
                    quadVarBeginInd,
                    quadVarEndInd,
                    quadVarPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < quadVarGenArray.size(); iVar++) {
                    svDBase.append(quadVarGenArray[iVar], quadVarSVArray[iVar]);
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
                    sqrtAnnualQuadVarEndInd[iAsset]   = sqrtAnnualQuadVarProdDates->size();
                }

                // Create MCQuadVar::IStateVars and put them in database
                MCPath::IStateVarArray sqrtAnnualQuadVarSVArray(MCPath::createPaths(
                    false,
                    sqrtAnnualQuadVarGenArray,
                    sqrtAnnualQuadVarDatesPerAsset,
                    sqrtAnnualQuadVarBeginInd,
                    sqrtAnnualQuadVarEndInd,
                    sqrtAnnualQuadVarPtrs,
                    maxDrifts));

                for(unsigned int iVar = 0; iVar < sqrtAnnualQuadVarGenArray.size(); iVar++) {
                    svDBase.append(sqrtAnnualQuadVarGenArray[iVar], sqrtAnnualQuadVarSVArray[iVar]);
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

            // pre compute fwds
             preComputeFwds();
          
		    // set up Dependence
            dependence = pathConfig->dependenceMaker->createDependence(this);

            // Initialize random number generator for spots
            // NB the size of the random deviates depends on the scheme 
            // that is requested
            int nbSpotDeviates;
            if (pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EULER){
                nbSpotDeviates = simDates->size() - 1;
            }
            else if(pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EXACT){
                nbSpotDeviates = logSpots.size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSV_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
            }
            int nbVarDeviates;
            if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                nbVarDeviates = simDates->size() - 1;
            }
            else{
                throw ModelException(method, 
                                     "unexpected spot discretization scheme of type "
                                     + VolSV_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            randomGenSpot = IMCRandomSP(new MCRandomNoCache(
                this,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbSpotDeviates,
                nbAssets,
                spotGenDates->size() - numFutSteps));
            randomGenVar = IMCRandomSP(new MCRandomNoCache(
                this,
                dependence,
                pathConfig->getRandomGenerator(),
                pathConfig->carefulRandoms(),
                nbVarDeviates,
                nbAssets,
                0));
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
        static const string method("MCPathConfigSV::PathGenFineProcess::getGaussData");
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
    }
   
    virtual void generatePath(int pathIdx, int iAsset) {
        static const string method = "MCPathConfigSV::Generator::generatePath";
        try{
            const DoubleMatrix& randomSpot = randomGenSpot->getRandomNumbers();
            const DoubleMatrix& randomVol = randomGenVar->getRandomNumbers();    
            // XXX this list of nested if statements is obviously not viable
            if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::EULER){
                if (pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EULER){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EULER>(),
                                                  Int2Type<VolSV_VarDSType::EULER>(),
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else if(pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EXACT){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EXACT>(),
                                                  Int2Type<VolSV_VarDSType::EULER>(), 
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else{
                    throw ModelException(method, 
                                         "unexpected spot discretization scheme of type "
                                         + VolSV_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
                }
            }
            else if (pathConfig->varDiscreteSchemeIndex == VolSV_VarDSTypeHelper::EnumList::VAR_TRANSFORM_EULER){
                if (pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EULER){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EULER>(),
                                                  Int2Type<VolSV_VarDSType::VAR_TRANSFORM_EULER>(),
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else if(pathConfig->spotDiscreteSchemeIndex == VolSV_SpotDSTypeHelper::EnumList::EXACT){
                    vols[iAsset]->MCGeneratePaths(Int2Type<VolSV_SpotDSType::EXACT>(),
                                                  Int2Type<VolSV_VarDSType::VAR_TRANSFORM_EULER>(), 
                                                  tradYears[iAsset],
                                                  fixDateOffsets,
                                                  randomSpot[iAsset],
                                                  randomVol[iAsset],
                                                  quanto[iAsset],
                                                  logSpots,
                                                  instVars[iAsset],
                                                  integratedVars[iAsset]);
                }
                else{
                    throw ModelException(method, 
                                         "unexpected spot discretization scheme of type "
                                         + VolSV_SpotDSTypeHelper::getName(pathConfig->spotDiscreteSchemeIndex));
                }
            }
            else{
                throw ModelException(method, 
                                     "unexpected variance discretization scheme of type "
                                     + VolSV_VarDSTypeHelper::getName(pathConfig->varDiscreteSchemeIndex));
            }
            // then, populate spot path
            const DoubleArray& fwds = this->fwds[iAsset];
            int nbFutFixDates = logSpots.size();
            DoubleArray& spots = spotProdPaths[iAsset];
            if (spots.size()){
                for (int iFutFixDate = 1, iSpotProdDate = spotProdOffset[iAsset];
                     iFutFixDate < nbFutFixDates; ++iFutFixDate, ++iSpotProdDate) {
                        spots[iSpotProdDate] = fwds[iFutFixDate] * exp(logSpots[iFutFixDate]);
                        // insert here drift adjustment for currency protection
                        // drift already adjusted
	            }
            }
            // followed by quad var path
            DoubleArray& quadVars = quadVarProdPaths[iAsset];
            if (quadVars.size()){
                for (int iFutFixDate = 1, iQuadVarProdDate = quadVarProdOffset[iAsset];
                     iFutFixDate < nbFutFixDates; ++iFutFixDate, ++iQuadVarProdDate) {
                        int iIntegratedVar = fixDateOffsets[iFutFixDate - 1];
                        quadVars[iQuadVarProdDate] = integratedVars[iAsset][iIntegratedVar];
	            }
            }
            // and sqrt annualized quad var path
            DoubleArray& sqrtAnnualQuadVars = sqrtAnnualQuadVarProdPaths[iAsset];
            if (sqrtAnnualQuadVars.size()){
                bool hasPast = sqrtAnnualQuadVarProdOffset[iAsset] > 0;
                for (int iFutFixDate = 1, iSqrtAnnualQuadVarProdDate = sqrtAnnualQuadVarProdOffset[iAsset],
                     iIntegratedVar = fixDateOffsets[0];
                     iFutFixDate < nbFutFixDates; 
                     ++iFutFixDate,
                     ++iSqrtAnnualQuadVarProdDate) {
                    if (hasPast){
                        // future annualized var from today till current sim date
						double tradYear = tradYears[iAsset][iIntegratedVar];
                        double futAnnualQuadVar = Maths::isZero(tradYear) ? 
                                                  instVars[iAsset][iIntegratedVar] :
												  integratedVars[iAsset][iIntegratedVar] / tradYear;
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
                            double quadVarInc = integratedVars[iAsset][iIntegratedVar]
                                                - integratedVars[iAsset][integratedVarStartIdx];
                            // time elapsed since first sampling date
                            double timeDiff = tradYears[iAsset][iIntegratedVar]
                                              - tradYears[iAsset][integratedVarStartIdx];
                            // total annualized vol
                            sqrtAnnualQuadVars[iSqrtAnnualQuadVarProdDate] 
                                = sqrt(quadVarInc / timeDiff);
                        }
                    }
                    if (iFutFixDate < fixDateOffsets.size()){
                        iIntegratedVar += fixDateOffsets[iFutFixDate];
                    }
	            }
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }


    /** pre-compute fwds arrays, get quanto data here for now */
    // Reminder about timelines:
    // No state variables: coarse path (timeline->futureDates); fine path (simDates);
    // State variables: coarse path (simTimeline); fine path (simDates);
     void preComputeFwds() {

        static const string routine = "MCPathConfigSV::preComputeFwds";

        try {

            // for each asset
            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                // need to cache the forward at all simulation dates
                fwds[iAsset].resize(simTimeline.size());

                CAssetConstSP asset = CAssetConstSP::attachToRef(&mAsset->getAsset(iAsset));
                
                if (CAsset::IStruck::TYPE->isInstance(asset)) {
                    throw ModelException(routine, "Struck assets (" + asset->getTrueName() + 
                                        ") are not supported yet");
                } 

                //using quanto of type QuantoParamSV
                bool isQuanto = Asset::IQuanto::TYPE->isInstance(asset);
                double eqFXCorr = 0.0;
                DoubleArray  fxSqrtVar = DoubleArray(simDates->size() - 1);

                if (isQuanto) {
                    const Asset::IQuanto& prot = dynamic_cast<const Asset::IQuanto&>(*asset.get());
                    
                    // Get the unprotected forward values
                    prot.unadjustedFwdValue(simTimeline, fwds[iAsset]);
                    
                    // get fx corr
                    eqFXCorr = prot.getCorrelation()->getCorrelation();
                    
                    // get fx vols atm
                    ATMVolRequestSP fxVolRequest(new ATMVolRequest());
                    CVolProcessedSP volFX(prot.getProcessedFXVol( fxVolRequest.get() ));
                    // cast to the type of vol we're expecting
                    CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(volFX);
                    volFXBS->CalcVol(*simDates, CVolProcessedBS::forward, fxSqrtVar);
                    
                    // calculate fx vars
                    for (int iStep = 0; iStep < fxSqrtVar.size(); iStep++) {
                        double dt = tradYears[iAsset][iStep+1]-tradYears[iAsset][iStep];
                        fxSqrtVar[iStep] *= dt; // using asset dt
                    }
                    quanto[iAsset] = QuantoParamSVSP(new QuantoParamSV(isQuanto,eqFXCorr,fxSqrtVar));
                } else{
                    quanto[iAsset] = QuantoParamSVSP(new QuantoParamSV(isQuanto,eqFXCorr,fxSqrtVar));
                    mAsset->factorFwdValues(iAsset, simTimeline, fwds[iAsset]);
                }
           
            }

        } catch (exception& e) {
        throw ModelException(e, routine);
        }
    }


    // fields
    MCPathConfigSV*         pathConfig;
    const IMultiFactors*    mAsset;            // multi asset interface
    int                     nbAssets;          // number of assets
    DateTimeArraySP         simDates;          // fine grid of simulation dates
    vector<DoubleArray>     fwds;              // fwds [iAsset][iStep]
    VolSVArray              vols;              // array of VolSVs
    TimeMetricArray         timeMetrics;       // array of time metrics
    vector<DoubleArray>     tradYears;         // fine grid of trad years
    DoubleArray             logSpots;          // spot path [iAsset][iStep]
    vector<DoubleArray>     instVars;          // instantaneous variance path [iAsset][iStep]
    vector<DoubleArray>     integratedVars;    // integrated variance path [iAsset][iStep]
    IntArray                fixDateOffsets;    // offsets for fixing dates
    IMCRandomSP             randomGenSpot;     // Random number generator for spots
    IMCRandomSP             randomGenVar;      // Random number generator for vars
    DependenceSP            dependence;        // Dependence object
    vector<double>          pastAnnualQuadVars; // past realized vars at value date [iAsset]
    vector<double>          pastWeights;       // past weights at value date [iAsset]
    vector<double>          maxDrifts;         // XXX product of MAX(1, drifts)
    vector<QuantoParamSVSP>     quanto;         //quanto parameters

    
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
class MCPathConfigSV::Gen: virtual public MCPathGenerator, // For backward compatibility
                           virtual public MCPathGen,
                           virtual public IStateVariableGen::IStateGen {
public:
    /** Constructor */
    Gen(int                      numSimPaths,
        MCPathConfigSV*          pathConfig, 
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient):
    nowPathIdx(0) {
        static const string routine = "MCPathConfigSV::Gen::Gen";
        
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
            pathGenFineProcess = MCPathConfigSV::PathGenFineProcessSP(
                new MCPathConfigSV::PathGenFineProcess(numSimPaths, 
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
        static const string routine = "MCPathConfigSV::Gen::NbSimAssets";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    const double* Path(int iAsset, int iPath) const {
        static const string routine = "MCPathConfigSV::Gen::Path";
        throw ModelException(routine, "Method is retired for StateVars");
    }; 

    double refLevel(int iAsset, int iPath) const {
        static const string routine = "MCPathConfigSV::Gen::refLevel";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    double maxDriftProduct(int iAsset) const {
        static const string routine = "MCPathConfigSV::Gen::maxDriftProduct";
        try {
            return pathGenFineProcess->maxDriftProduct(iAsset);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    int begin(int iAsset) const {
        static const string routine = "MCPathConfigSV::Gen::begin";
        throw ModelException(routine, "Method is retired for StateVars");
    }

    int end(int iAsset) const{
        static const string routine = "MCPathConfigSV::Gen::end";
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
        static const string routine = "MCPathConfigSV::Gen::create";

        try {
            return svDBase.find(svGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }
  
private:
    StateVarDBase                        svDBase;               //!< Collection of Generators + statevars
    MCPathConfigSV::PathGenFineProcessSP pathGenFineProcess;    //!< Spot generator
    int                                  nowPathIdx;            //!< Current path idx
};
#endif

/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigSV::makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                   prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     mainSimDates ){
    static const string routine = "MCPathConfigSV::makePathGenerator";

    if(cachingRequested) {
        throw ModelException(routine, "Paths caching is not supported in MCPathConfigSV");
    }

    // create empty sim dates array if simDates is null or on new pricing run
    if (!simDates || control->isPricing()){
        simDates = DateTimeArraySP(new DateTimeArray());
    }

    const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
    if(!prodClient){
        refCountPtr<MCPathConfigSV::Generator> futurePathGenBase(
            new MCPathConfigSV::Generator(numPaths,
                                          this,
                                          pastPathGenerator,
                                          prod));
    
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
        return MCPathGeneratorSP(new MCPathConfigSV::Gen(
            numPaths,
            this, 
            pastPathGenerator,
            prodClient));
    }
}

/** Creates a past path generator */
MCPathGeneratorSP MCPathConfigSV::pastPathGenerator(const IMCProduct* prod) {
    static const string method = "MCPathConfigSV::pastPathGenerator";
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
class MonteCarloSVDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloSVDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloSVDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloSV);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloSV(){
        return new MonteCarloSVDefault();
    }
};

CClassConstSP const MonteCarloSVDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloSVDefault", typeid(MonteCarloSVDefault), load);

// external symbol to allow class to be forced to be linked in
bool MCPathConfigSVLoad(){
    return (MCPathConfigSV::TYPE != 0 && MonteCarloSVDefault::TYPE !=0);
}



DRLIB_END_NAMESPACE
