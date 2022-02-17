//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathGeneratorLV.cpp
//
//   Description : Monte Carlo path generator using LV
//
//   Date        : 6 Jan 2003
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/depnet.hpp"
#include "edginc/LocalVolGrid.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"
#include "edginc/DependenceLocalCorr.hpp"

DRLIB_BEGIN_NAMESPACE

#define P refCountPtr

using namespace depnet;

#define MAX_STEP_SIZE 100               // max time step size 100 yrs
// make sure that max step size is MAX_STEP_SIZE
#define CAP_DT(dt)      ((dt)=(dt)>(MAX_STEP_SIZE)?(MAX_STEP_SIZE):(dt))

struct LocalVolCache;

class MCPathConfigLV: public MCPathConfig {
public:
    static CClassConstSP const TYPE;

    /** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
        for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const{
        MarketDataFetcherSP mdf(new MarketDataFetcherLNSpline(volType));        
        dependenceMaker->modifyMarketDataFetcher(mdf);
        if(skewMaker.get()) {
            skewMaker->modifyMarketDataFetcher(mdf);
        }
        return mdf;        
    }

    virtual int randomStoragePerPath(IMCProduct* product) const;
    
    string          volType;
    string          dependenceType;

    double          skewStepScaling;
    int             numVolGrid;
    double          stdevGridRange;

    bool            useTweakingForTimeDerivs;
    double          tweakStrikeUnscaled;
    double          tweakTimeUnscaled;
    double          probDensRatioMin;
    bool            useMidPoint;

    // Possibilities for Local Vol tweaks
    enum LVTweaks {
        TWEAK_LOCAL_VOL = 0,
        FIX_LOCAL_VOL_EXCEPT_THETA = 1,
        FIX_LOCAL_VOL = 2
    };
    
    LVTweaks        tweakLocalVols;
    bool            cacheLocalVols;

private:

    P<LocalVolCache> cache; // $unregistered

    static IObject* emptyShell() { return new MCPathConfigLV(); }
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(MCPathConfigLV, clazz);
        SUPERCLASS(MCPathConfig);
        EMPTY_SHELL_METHOD(emptyShell);

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
        FIELD(tweakLocalVols, "Possibilities for LocalVol tweaks");
        FIELD_MAKE_OPTIONAL(tweakLocalVols);
        FIELD(cacheLocalVols, "Cache local vols across tweaks (on by default)");
        FIELD_MAKE_OPTIONAL(cacheLocalVols);
    }

    MCPathConfigLV():
        MCPathConfig(TYPE,
                     DependenceMakerGaussTermSP(new DependenceMakerGaussTerm())),
        volType("VolPreferred"),
        dependenceType("not used"),
        skewStepScaling(200.0),
        numVolGrid(101),
        stdevGridRange(6.0),
        useTweakingForTimeDerivs(true),
        tweakStrikeUnscaled(0.01),
        tweakTimeUnscaled(0.001),
        probDensRatioMin(0.01),
        useMidPoint(false),
        tweakLocalVols(TWEAK_LOCAL_VOL),
        cacheLocalVols(true)
    {}

    class Gen;              //!< Referee path generator class
    class PathGenSpot;      //!< Spot PathGen
    
    friend class Gen;
    friend class PathGenSpot;
    
    typedef refCountPtr<PathGenSpot> PathGenSpotSP;

public:

    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                  prod,
        Control*                           control, 
        Results*                           results,
        DateTimeArray&                     simDates );

    /** MCPathConfig overriden to switch state of RefLevel */
    MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod) {
        const MCProductClient* prodClient =
            dynamic_cast<const MCProductClient*>(prod);
        if (!prodClient) {
            // Old framework
            const_cast<IMCProduct*>(prod)->setSimStartDateToFirstRefDate();
            return MCPathConfig::pastPathGenerator(prod);
        }
        else {
            // State variables framework
            return MCPathGeneratorSP(new PastPathGen(prodClient));
        }
    }

    virtual bool vegaMatrixSupported() const { return false; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * Irrelevant since MCLV recalibrates the local vols on the fly.  See
     * IModel::wantsRiskMapping().
     */

    IModel::WantsRiskMapping wantsRiskMapping() const {
        return IModel::riskMappingIrrelevant;
    }

    virtual bool carefulRandoms() const { return false; }

    /** Throws away cached sim dates for Theta-type tweaks */
    virtual bool sensShift(Theta* shift) {
        MCPathConfig::sensShift(shift);
        cache.reset();
        return true;
    }
};

CClassConstSP const MCPathConfigLV::TYPE = CClass::registerClassLoadMethod(
    "MCPathConfigLV", typeid(MCPathConfigLV), load);

// Enum registration
START_PUBLIC_ENUM_DEFINITION(MCPathConfigLV::LVTweaks, "Values for LVTweaks");
ENUM_VALUE_AND_NAME(MCPathConfigLV::TWEAK_LOCAL_VOL, "TWEAK_LOCAL_VOL", 
                    "Tweak Local Vols i.e. Black Scholes Greeks");
ENUM_VALUE_AND_NAME(MCPathConfigLV::FIX_LOCAL_VOL_EXCEPT_THETA, "FIX_LOCAL_VOL_EXCEPT_THETA", 
                    "Fix Local Vol between all tweaks except Theta");
ENUM_VALUE_AND_NAME(MCPathConfigLV::FIX_LOCAL_VOL, "FIX_LOCAL_VOL", 
                    "Fix Local Vol between all tweaks");
END_ENUM_DEFINITION(MCPathConfigLV::LVTweaks);

// 
// ***************
//  LocalVolCache
// ***************
// 

static CVolProcessedDVFConstSP localVol(
        const IMultiFactors* factors, int factor,
        DateTime today, DateTime volStartDate,
        const MCPathConfigLV* pathConfig) {
    TRACE_METHOD;
    LocVolRequest volRequest(volStartDate,
                             volStartDate > today, // fwd start flag
                             false,
                             pathConfig->useTweakingForTimeDerivs,
                             pathConfig->useMidPoint,
                             pathConfig->tweakStrikeUnscaled,
                             pathConfig->tweakTimeUnscaled,
                             pathConfig->probDensRatioMin);

    return CVolProcessedDVFSP::dynamicCast(
      CVolProcessedSP(factors->factorGetProcessedVol(factor, &volRequest)));
}

static void computeSimDates(const DateTimeArray& futureObservationDates,
                            const vector<CVolProcessedDVFConstSP>& localVols,
                            const DateTime& startDate,
                            const vector<double>& spots,
                            const double& skewStepScaling,
                            DateTimeArray& simDates) {
    TRACE_METHOD;

    const int numAssets = localVols.size();
    const DateTime& mat = futureObservationDates.back();
    const double one_week = 1.0/52.0;
    const double adjFactor = 2.0; // term structure scaling adjustment, an empirical factor here

    CSliceDouble k(2), var(2);
    double tmp, maxVolDiff = 0.0;
    DateTimeArray t(2), t_next(2);

    // just use the one time metric
    TimeMetricConstSP metric = localVols[0]->GetTimeMetric();

    // use 90-100 skew (vol*sqrt(dt)) and term structure to determine step size
    // compute starting maxVolDiff
    t[0] = metric->impliedTime(startDate, 1.0/365, tmp); // QuadExp does not give skew  for start=today
    t[1] = metric->impliedTime(t[0], one_week, tmp); // tmp not used, one week loc vol*sqrt(dt)
    for (int iAsset=0; iAsset<numAssets; iAsset++)
    {
        // compute 90-100 fwd local vol spread
        // use fwd at start as strike in vol interp to get the skew
        // equivalent to today's ATM
        k[0] = 0.9*spots[iAsset];
        k[1] = spots[iAsset];

        localVols[iAsset]->CalcLocVar(&k, t, &var, false); // no intradayinterp
        tmp = fabs(sqrt(var[0]) - sqrt(var[1]));
        if (tmp > maxVolDiff)
            maxVolDiff = tmp;

        // also take into account of term structure
        var[0] = localVols[iAsset]->computeImpVol(t[1], k[1]);
        var[1] = localVols[iAsset]->computeImpVol(mat, k[1]);
        tmp = adjFactor*fabs(var[1] - var[0])*sqrt(one_week);
        if (tmp > maxVolDiff)
            maxVolDiff = tmp;
    }

    // now space time axis according to maxVolDiff
    simDates.clear();
    if (skewStepScaling == -1.0) // daily steps requested
    {
        double dt = 1.0/metric->volDays(startDate, startDate.rollDate(365));
        DateTime nextDate(metric->impliedTime(startDate, dt, tmp)); // tmp not used
        while (nextDate < mat)
        {
            simDates.push_back(nextDate);
            nextDate = metric->impliedTime(nextDate, dt, tmp);
        }
    }
    else if (skewStepScaling < -1) {
        double perYear = -skewStepScaling;

        for (int o = 0; o < futureObservationDates.size(); ++o) {
            DateTime last = o ? futureObservationDates[o-1] : startDate;
            double dt = metric->yearFrac(last, futureObservationDates[o]);
            int steps = int(ceil(perYear * dt));
            for (int e = 1; e < steps; ++e) {
                simDates.push_back(metric->impliedTime(last, e*dt/steps, tmp));
            }
            simDates.push_back(futureObservationDates[o]);
        }
    }
    else if (maxVolDiff>0.0 && skewStepScaling > 0.0) // if starting skew is 0 then we assume zero skew
    {
        CSliceDouble var_next(2);
        double dt = 1.0/100.0/maxVolDiff/skewStepScaling; // const does not matter
        CAP_DT(dt); // cap dt at MAX_STEP_SIZE. 

        DateTime nextDate(metric->impliedTime(startDate, dt, tmp)); // tmp not used
        while (nextDate < mat)
        {
            simDates.push_back(nextDate);

            // use 90-100 skew (vol*sqrt(dt)) and term structure to determine step size
            // compute starting maxVolDiff
            t[0] = nextDate;
            t[1] = metric->impliedTime(t[0], one_week, tmp); // tmp not used
            t_next[0] = metric->impliedTime(t[0], dt, tmp);
            t_next[1] = metric->impliedTime(t_next[0], one_week, tmp);
            maxVolDiff = 1.0e-20; // set to almost 0
            for (int i=0; i<numAssets; i++)
            {
                // compute 90-100 fwd local vol spread
                k[0] = 0.9*spots[i];
                k[1] = spots[i];

                localVols[i]->CalcLocVar(&k, t, &var, false); // no intradayinterp

                tmp = fabs(sqrt(var[0]) - sqrt(var[1]));
                if (tmp > maxVolDiff)
                    maxVolDiff = tmp;

                // also take into account of term structure
                // this makes the algorithm more stable. otherwise, vol surface noise may give tmp~0 when it shouldn't be 
                localVols[i]->CalcLocVar(&k, t_next, &var_next, false); // no intradayinterp
                tmp = fabs(sqrt(var_next[0]) - sqrt(var_next[1]));
                if (tmp > maxVolDiff)
                    maxVolDiff = tmp;
                tmp = adjFactor*fabs(sqrt(var_next[1]) - sqrt(var[1]));
                if (tmp > maxVolDiff)
                    maxVolDiff = tmp;
                tmp = adjFactor*fabs(sqrt(var_next[0]) - sqrt(var[0]));
                if (tmp > maxVolDiff)
                    maxVolDiff = tmp;
            }

            dt = 1.0/100.0/maxVolDiff/skewStepScaling;
            CAP_DT(dt); // cap dt at MAX_STEP_SIZE. 
            nextDate = metric->impliedTime(nextDate, dt, tmp);
        }
    }

    // merge in the instrument monitoring dates

    if (skewStepScaling >= -1) {
        simDates.insert(simDates.end(), 
                        futureObservationDates.begin(), 
                        futureObservationDates.end());
    }

    DateTime::doSortUniq(simDates);
}

static DateTimeArray cumDates(const IMultiFactors* factors, int f,
                              DateTime from, DateTime to) {
    TRACE_METHOD;

    // obtain and sort cum div dates. no local vol averaging across div dates

    EventAssetMove divEvent;
    int numDivs = (int)(4*from.yearFrac(to))+1; // 4 divs per year selected as critical dates
    DateTimeArray cumDates;

    if (AssetUtil::getJumpEvents(&factors->getAsset(f), from, to,
                                 numDivs, divEvent)) {
        cumDates = *divEvent.getCritDate(0, true);
        sort(cumDates.begin(), cumDates.end());
    }

    return cumDates;
}

struct Drift {
    DoubleArray forward;
    DoubleArray logForward;
    double eqFXCorr;
    DoubleArray fxSqrtVar;

    void set(const IMultiFactors* factors, int f,
             const DateTimeArray& simulationDates,
             CVolProcessedDVFConstSP localVol);
};

void Drift::set(const IMultiFactors* factors, int f,
                const DateTimeArray& simulationDates,
                CVolProcessedDVFConstSP localVol) {
    TRACE_METHOD;

    CAssetConstSP asset = CAssetConstSP::attachToRef(&factors->getAsset(f));

    if (CAsset::IStruck::TYPE->isInstance(asset)) {
        throw ModelException(__FUNCTION__,
                             "Struck assets (" + asset->getTrueName() +
                                 ") are not supported yet");
    }

    forward.resize(simulationDates.size());
    logForward.resize(simulationDates.size());

    if (Asset::IQuanto::TYPE->isInstance(asset)) {
        // current way of dealing with quanto

        const Asset::IQuanto& prot = dynamic_cast<const Asset::IQuanto&>(*asset.get());

        // Get the unprotected forward values
        prot.unadjustedFwdValue(simulationDates, forward);

        // get fx corr
        eqFXCorr = prot.getCorrelation()->getCorrelation();

        // get fx vols atm
        ATMVolRequestSP fxVolRequest(new ATMVolRequest());
        CVolProcessedSP volFX(prot.getProcessedFXVol( fxVolRequest.get() ));
        // cast to the type of vol we're expecting
        CVolProcessedBSSP volFXBS = CVolProcessedBSSP::dynamicCast(volFX);

        fxSqrtVar.resize(simulationDates.size() - 1);
        volFXBS->CalcVol(simulationDates, CVolProcessedBS::forward, fxSqrtVar);

        // just for getting dt for each asset
        TimeMetricConstSP metric = localVol->GetTimeMetric();
        for (int j = 0; j < fxSqrtVar.size(); ++j) {
            fxSqrtVar[j] *= sqrt(metric->yearFrac(simulationDates[j],
                                                  simulationDates[j+1]));
        }
    }
    else {
        factors->factorFwdValues(f, simulationDates, forward);
        fxSqrtVar.clear();
        eqFXCorr = 0.;
    }

    // Compute log forwards
    for(int s = 0; s < forward.size(); s++) {
        logForward[s] = log(forward[s]);
    }
}

static DependenceSP gaussDependence(DependenceMakerConstSP maker,
                                    CDoubleMatrixConstSP corrs,
                                    const DateTimeArray& simDates,
                                    const IMultiFactors* mfs,
                                    int nbPastDates,
                                    MCPathConfig* pathConfig){
    struct Corrs: DependenceMakerGauss::Support,
                  DependenceMakerGaussTerm::Support,
                  DependenceMakerLocalCorr::Support {
        CDoubleMatrixConstSP corrs;
        DateTimeArray simDates;
        const IMultiFactors* mfs;
        int nbPastDates;
        MCPathConfig* pathConfig;

        // implement getGaussData since derivation from DependenceMakerGauss::Support

        CDoubleMatrixConstSP getGaussData() const { return corrs; }

        // implement methods since derivation from DependenceMakerGaussTerm::Support

        DateTimeArray getSimDates() const {
            return simDates;
        }

        DoubleMatrix getFwdVarAtDates(bool) const {
            DateTimeArray toDates = simDates;
            DateTime valueDate = toDates.front();
            toDates.erase(toDates.begin());
            DoubleArrayArray var(mfs->NbAssets());
            for(int iAsset=0; iAsset<mfs->NbAssets(); iAsset++ ) {
                var[iAsset] = DoubleArray(toDates.size());
                const CAsset& thisAsset = mfs->getAsset(iAsset); 
                ATMVolRequestSP volRequest(new ATMVolRequest());
                volRequest->allowNegativeFwdVar(false); 
                CVolProcessedBSSP volBS(thisAsset.getProcessedVol(volRequest.get()));
                volBS->CalcVar(valueDate, 
                               toDates, 
                               volBS->forward,
                               var[iAsset]);
            }
            return DoubleMatrix(var);
        }

        const IMultiFactors* getMultiFactors() const {
            return mfs;
        }

        // implement methods since derivation from DependenceMakerLocalCorr::Support
        DependenceMakerSP getDependenceMaker() const {
            return pathConfig->getDependenceMaker();
        }
        int getNbPastDates() const {
            return nbPastDates;
        }
        const IRandomSP getRandomGenerator() const {
            return pathConfig->getRandomGenerator();
        }        
        const MCPathConfig::RandomCacheSP getIdiosynFactorRandomCache() const {
            return pathConfig->getIdiosynFactorRandomCache();
        }        
        const MCPathConfig::RandomCacheSP getMarketFactorRandomCache() const {
            return pathConfig->getMarketFactorRandomCache();
        }
        bool carefulRandoms() const {
            return pathConfig->carefulRandoms();                      
        }        
    };

    Corrs corrrs;
    corrrs.corrs = corrs;
    corrrs.simDates = simDates;
    corrrs.mfs = mfs;
    corrrs.nbPastDates = nbPastDates;
    corrrs.pathConfig = pathConfig;
    //pathConfig->getRandomGenerator()->init();
    if (pathConfig->getSkewMaker().get()) {
        DependenceMakerLocalCorr* dpm = 
            dynamic_cast<DependenceMakerLocalCorr*>(pathConfig->getSkewMaker().get());
        if (!dpm) {
            throw ModelException("method", "Internal Error");
        }
        return dpm->createDependence(&corrrs);       
    } else {
        return maker->createDependence(&corrrs);
    }
}

struct NoopMCRandomCallbacks: IMCRandom::Callbacks {
    void configureAntithetics() {}
    void configureNonAntithetics() {}
};

struct LocalVolCache: CObject {

    bool tweaked;
    P<Var<MCPathConfigLV*> > pathConfig;

    // what is the relation between these two?
    P<Var<DateTimeArray> > futureObservationDates; // put here so we can ASSERT below
    P<Var<int> > numPastDates;
    P<Var<DateTime> > volStart;
    P<Var<DateTime> > today;
    P<Var<const IMultiFactors*> > multiFactors;

    vector<P<Var<int> > > factors;
    P<Var<CDoubleMatrixConstSP> > correls;
    P<Var<DateTimeArray> > simulationDates;
    vector<P<Var<double> > > spotsAtStart;
    vector<P<Var<double> > > spotsAtVolStart;
    vector<P<Var<Drift> > > drifts;
    vector<P<Var<ILocalVolGridSP> > > localVolGrids;
    P<Var<DependenceSP> > dependence;

    MCRandomSP randomGen;
    vector<BoolArray> observed;

    LocalVolCache(int numFactors):
        CObject(TYPE),
        tweaked(false),
        pathConfig(new Var<MCPathConfigLV*>()),
        futureObservationDates(withNoChangeTest(new Var<DateTimeArray>())),
        numPastDates(withNoChangeTest(new Var<int>())),
        volStart(withNoChangeTest(new Var<DateTime>())),
        today(withNoChangeTest(new Var<DateTime>())),
        multiFactors(new Var<const IMultiFactors*>()),
        factors(numFactors),
        spotsAtStart(numFactors),
        spotsAtVolStart(numFactors),
        drifts(numFactors),
        localVolGrids(numFactors),
        observed(max(numFactors, 1))
    {
#define NAME(X) X->name = #X
        NAME(pathConfig); NAME(futureObservationDates); NAME(numPastDates); NAME(volStart); NAME(today); NAME(multiFactors);
#undef NAME        

        P<Var<DateTime> > maturity = mappedGet(futureObservationDates, -1);

        vector<P<Var<CVolProcessedDVFConstSP> > > localVols(numFactors);

        for (int f = 0; f < numFactors; ++f) {
            factors[f] = root(f);
            spotsAtVolStart[f] = mapped(&IMultiFactors::factorFwdValue,
                                        multiFactors, factors[f], volStart);
            spotsAtStart[f] = mapped(&IMultiFactors::factorFwdValue,
                                     multiFactors, factors[f],
                                     mappedGet(futureObservationDates, 0));
            localVols[f] = mapped(localVol, multiFactors, factors[f], today,
                                  volStart, pathConfig);
        }

        simulationDates = mapped(computeSimDates, futureObservationDates,
                                 grouped(localVols), volStart,
                                 grouped(spotsAtVolStart),
                                 mapped_member(&MCPathConfigLV::skewStepScaling,
                                               pathConfig));
        for (int f = 0; f < numFactors; ++f) {
            drifts[f] = setWith(&Drift::set, multiFactors, factors[f],
                                simulationDates, localVols[f]);
        }

        for (int f = 0; f < numFactors; ++f) {
            localVolGrids[f] = mapped(
                LocalVolGrid::createLocalVolGrid,
                mapped_member(&Drift::forward, drifts[f]),
                localVols[f],
                simulationDates,
                mapped(cumDates, multiFactors, factors[f], volStart, maturity),
                mapped_member(&MCPathConfigLV::numVolGrid, pathConfig),
                mapped_member(&MCPathConfigLV::stdevGridRange, pathConfig));
        }

        correls = mapped(&IMultiFactors::factorsCorrelationMatrix,
                         multiFactors);

        dependence = mapped(
            gaussDependence,
            mapped(&MCPathConfig::getDependenceMaker, pathConfig),
            correls,
            simulationDates,
            multiFactors,
            numPastDates,
            pathConfig);
    }

    MCRandomSP newMCRandom() const {
        static NoopMCRandomCallbacks callbacks;

        if (pathConfig->value()->getSkewMaker().get() && (numFactors()>1)) {
            // a special one IRandomSP
            return MCRandomSP(new MCRandom(
                &callbacks,                             // path generator
                dependence->value(),                    // dependence
                IRandomSP(),                            // random number generator
                pathConfig->value()->getRandomCache(),  // randomCache
                true,                                   // not used (isCarefulRandoms)
                simulationDates->value().size() - 1,    // numDates
                numFactors(),                           // numFactors
                numPastDates->value() - futureObservationDates->value().size()));
        } else {
            // the classic one         
            return MCRandomSP(new MCRandom(
                &callbacks, dependence->value(),
                pathConfig->value()->getRandomGenerator(),
                pathConfig->value()->getRandomCache(),
                pathConfig->value()->carefulRandoms(),
                simulationDates->value().size() - 1,
                numFactors(),
                numPastDates->value() - futureObservationDates->value().size()));
        }
    }

    int numFactors() const { return factors.size(); }

    void generatePath(int pathIndex, int asset, double *path) const;
};

/** Generate path for specified path in simulation (pathIdx), for
    specified asset */

void LocalVolCache::generatePath(int pathIndex, int asset, double* path) const {
    const double* randoms = randomGen->getRandomNumbers()[asset];

    int numPathSteps = simulationDates->value().size();
    // populate paths field
    const DoubleArray& forward = drifts[asset]->value().forward;
    const DoubleArray& logForward = drifts[asset]->value().logForward;
    const ILocalVolGrid& pathVols = *localVolGrids[asset]->value();
    double eqFXCorr = drifts[asset]->value().eqFXCorr;
    const DoubleArray& fxSqrtVar = drifts[asset]->value().fxSqrtVar;
    const BoolArray& assetObserved = observed[asset];

    double logS = log(spotsAtStart[asset]->value() / forward[0]);
    for (int s = 0, p = 0; s < numPathSteps; ++s) {
        if (s) {
            double assetSqrtVar = pathVols.interpLocVar(s-1, logS);
            // double assetSqrtVar = pathVols.interpLocVar(s-1, logS + logForward[s-1]);

            // add quanto adjustment
            if (eqFXCorr != 0.)
                logS += -eqFXCorr * fxSqrtVar[s-1] * assetSqrtVar;

            logS += assetSqrtVar * randoms[s-1]
                  - 0.5 * assetSqrtVar * assetSqrtVar;
        }

        if (assetObserved[s]) {
            path[p] = forward[s] * exp(logS);
            ++p;
        }
    }
}

/** Returns the number of bytes used for random number storage per
    path. Do not invoke if there are no sim dates in future */
int MCPathConfigLV::randomStoragePerPath(IMCProduct* product) const{
    // Otherwise have to worry about how many random numbers we use
    // This is a bit of a pain as currently there is no easy way to get hold
    // of the number of dates unless we build the entire path generator

    // If you think about it's particularly silly to have to do this in the
    // "dependency graph" implementation, but it's all too complicated for
    // me to want to unravel

    smartPtr<MCPathConfigLV> pathConfig(copy(this)); // to avoid const problems
    MCPathGeneratorSP pastPathGenerator(pathConfig->pastPathGenerator(product));
    product->pathGenUpdated(pastPathGenerator.get());
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

    int numDates = pathConfig->cache->simulationDates->value().size();

    if (!numDates){
        throw ModelException("MCPathConfigLV::storagePerPath", "Internal "
                             "error - no simulation dates");
    }
    // we store both correlated and uncorrelated numbers but only for
    // every other path. Hence no times by 2.
    return sizeof(double) * product->getNumAssets() * numDates;
}


/***************************************** 
MCPathBaseLW class
*****************************************/

class MCPathBaseLW: virtual public MCPathBase {

    P<LocalVolCache> cache; // $unregistered
    vector<DoubleArray>        productPaths; // [iAsset][iStep]
    RefLevelDataSP             refData;     // Ref level data
    MCProductTimelineSP        timeline;    // Product timeline

public:
    MCPathBaseLW(P<LocalVolCache>          cache,
                 const MCPathGeneratorSP&  pastPathGenerator,
                 const IMCProduct*         prod);

    /** Returns the MAX(1,d_0)*MAX(1,d_1)*...*MAX(1,d_n) where d_i is
        the drift between simulation date i and simulation date i+1.
        (With d_0 = simulation start date). Range of dates is restricted to
        those in the future so the past path generator will
        return 1.0 */
    double maxDriftProduct(int iAsset) const{
        return cache->localVolGrids[iAsset]->value()->maxDrift();
    }

    /** Obtains timeline object from base */
    MCProductTimelineConstSP getTimeline() const {
        return timeline;
    }

    /** Returns number of assets */
    int NbSimAssets() const;

    /** Returns the reference level for iAsset, iPath */
    double& refLevel(int iAsset, int iPath);
    
    /** Returns the reference level for iAsset, iPath */
    // const double& refLevel(int iAsset, int iPath) const;

    /** Returns the reflevel path */
    IRefLevel::IMCPathSP& refLevelPath() const;

    /** Returns the path for iAsset, iPath */
    const double* Path(int iAsset, int iPath) const;

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

    /** Simulates random numbers */
    virtual void drawRandomNumbers(int pathIdx);

    void generatePath(int pathIdx, int iAsset, int iPath);
};

////////////////////////////////////////////////
// MCPathBaseLW class member functions
////////////////////////////////////////////////
MCPathBaseLW::MCPathBaseLW(P<LocalVolCache>          cache,
                           const MCPathGeneratorSP&  pastPathGenerator,
                           const IMCProduct*         prod):
    cache(cache),
    productPaths(cache->numFactors()) {

    TRACE_METHOD;

    int numAssets = cache->numFactors();

    try{
        // Obtain product timeline

        timeline = getProductTimeline(prod, pastPathGenerator);
        cache->futureObservationDates->setValue(timeline->futureDates, cache->tweaked);
        cache->numPastDates->setValue(timeline->totalNumPastDates, cache->tweaked);
        //ASSERT(cache->futureObservationDates->value() == timeline->futureDates);

        // Obtain market data
        IntArray nbPaths(numAssets, 1);
        refData = getRefData(timeline, nbPaths, prod->getMultiFactors(),
                             prod, pastPathGenerator);

        // Paths
        int iAsset;
        for (iAsset=0; iAsset < numAssets; iAsset++) {
            productPaths[iAsset] = DoubleArray(timeline->totalNumSteps);
        }

        // get the vol 'start dates' (includes historic dates)
        DateTimeArray cliquetDates(prod->getVolStartDates(pastPathGenerator.get()));
        if (cliquetDates.size() != 1){
            throw ModelException(__FUNCTION__, "Cliquets not supported by Local Vol MC");
        }

        // Recall, cliquetDates might be in the past
        cache->volStart->setValue(timeline->today > cliquetDates.front() ?
                                      timeline->today : cliquetDates.front(),
                                  cache->tweaked);
        if (!cache->volStart->value().equals(prod->getEffectiveSimStartDate())){
            // the call to setSimStartDateToFirstRefDate should means this
            // never happens
            throw ModelException(__FUNCTION__, 
                                 "RefLevel sim start date != vol start date");
        }

        // Why is this apparently different from the state vars version?
        // Dunno.

        // loop backwards through futureDates to assign isSampleDate flag
        // resistant to possible duplicate future dates

        const DateTimeArray& simulationDates = cache->simulationDates->value();

        BoolArray& observed = cache->observed[0];
        observed.resize(simulationDates.size());

        int j = observed.size()-1;
        for(int iSearch = 0; iSearch < timeline->numFutSteps; iSearch++)
        {
            DateTime futDate = (&(timeline->futureDates.back()))[-iSearch];
            while( j >= 0 && simulationDates[j] > futDate )
                j--;
            if( j < 0 || simulationDates[j] != futDate )
                throw ModelException(__FUNCTION__, "Internal error. Future date not in simDatesLV");
            
            observed[j] = true;
        }

        for (int i = 1; i < numAssets; ++i) {
            cache->observed[i] = cache->observed[0];
        }

    } catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

/** Returns number of assets */
int MCPathBaseLW::NbSimAssets() const {
    return cache->numFactors();
}

/** Returns the reference level for iAsset, iPath */
double& MCPathBaseLW::refLevel(int iAsset, int iPath) {
    return refData->refLevels[iAsset][0];
}
    
/** Returns the reference level for iAsset, iPath */
//const double& MCPathBaseLW::refLevel(int iAsset, int iPath) const {
//    return refData->refLevels[iAsset][0];
//}

/** Returns the reflevel path */
IRefLevel::IMCPathSP& MCPathBaseLW::refLevelPath() const {
    return refData->refLevelPath;
}

/** Returns the path for iAsset, iPath */
const double* MCPathBaseLW::Path(int iAsset, int iPath) const {
    return &productPaths[iAsset][0];
}

/** Returns the paths */
double* MCPathBaseLW::Path(int iAsset, int iPath) {
    return &productPaths[iAsset][0];
}

void MCPathBaseLW::drawRandomNumbers(int pathIdx) {
    // Delegate to MCRandom class
    cache->randomGen->generate(pathIdx);
}

/** Generate path for specified path in simulation (pathIdx), for
    specified asset (iAsset) and specified vol interp (iPath) */
void MCPathBaseLW::generatePath(int pathIdx, int iAsset, int iPath){
    cache->generatePath(pathIdx, iAsset,
                        &productPaths[iAsset][timeline->numPastDates]);
}

//*******************************************************************
//* State variables
//*******************************************************************


/** Referee class that distributes simulation to components
    e.g. spot etc. */
class MCPathConfigLV::Gen: virtual public MCPathGenerator, // For backward compatibility
                           virtual public IStateVariableGen::IStateGen {
public:
    /** Constructor */
    Gen(P<LocalVolCache>         cache, 
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient,
        bool                     cachingRequested);

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

    virtual int getPathIndex() const;
    
    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen);

private:
    StateVarDBase                    svDBase;        //!< Collection of Generators + statevars
    MCPathConfigLV::PathGenSpotSP    pathGenSpot;    //!< Spot generator
    int                              nowPathIdx;     //!< Current path idx
};

/** Spot path generator for LV process using state variable approach */
class MCPathConfigLV::PathGenSpot: virtual public MCPathGen {

    P<LocalVolCache> cache;
    bool                       simHasPast;      //!< Whether product has past
    vector<BoolArray>          observed;        //!< Simulate spot   at simTimeline     [iAsset][iStep]
    vector<DoubleArray>        spotProdPaths;   //!< Spot price at asset specific dates [iAsset][iStep]
    vector<int>                spotProdOffset;  //!< Starting point for future path per asset

public:
    
    /** Constructor */
    PathGenSpot(P<LocalVolCache> cache,
                const PastPathGenSpotSP& pastPathGenSpot,
                const MCProductClient*   prodClient,
                const SVGenSpotArray&       spotGenArray,
                bool                     cachingRequested,
                StateVarDBase&           svDBase):
        cache(cache),
        simHasPast(pastPathGenSpot->hasPast()),
        observed(cache->numFactors()),
        spotProdPaths(cache->numFactors()),
        spotProdOffset(cache->numFactors())
    {
        try {
            int numAssets = cache->numFactors();
            DateTime today = prodClient->getToday();

            DateTimeArraySP spotGenDates = MCPath::getAllDates(spotGenArray);
            DateTimeArray futDates = today.getFutureDates(*spotGenDates.get());
            if (today >= spotGenDates->front()) {
                futDates.insert(futDates.begin(), today);
            }
            cache->futureObservationDates->setValue(futDates, cache->tweaked);
            cache->numPastDates->setValue(
                spotGenDates->size() - futDates.size(), cache->tweaked);
            cache->volStart->setValue(today > spotGenDates->front() ?
                                          today : spotGenDates->front(),
                                      cache->tweaked);

            // Create spot paths per asset
            DateTimeArrayArray spotDatesPerAsset(numAssets);
            vector<const double*> spotPtrs(numAssets);
            vector<int> spotBeginInd(numAssets);
            vector<int> spotEndInd(numAssets);

            const IPastValues* pastValues = prodClient->getMCPastValues();
            for (int asset = 0; asset < numAssets; ++asset) {
                // 1) SPOT PATHS
                // Populate past values
                DateTimeArraySP spotAssetDates = 
                    MCPath::getAllDates(spotGenArray, asset);
                spotProdPaths[asset]   = DoubleArray(spotAssetDates->size());
                DateTimeArray assetPastDates = today.getPastDates(*spotAssetDates);
                DoubleArray assetPastValues = 
                    pastValues->getPastValues(assetPastDates, asset, today);
                spotProdOffset[asset] = assetPastValues.size();
                int s;
                for (s = 0; s < assetPastValues.size(); s++) {
                    spotProdPaths[asset][s] = assetPastValues[s];
                }

                // Create spot mappings
                spotDatesPerAsset[asset] = *spotAssetDates;
                spotPtrs[asset] = &spotProdPaths[asset][0];
                spotBeginInd[asset] = assetPastDates.size();
                spotEndInd[asset]   = spotAssetDates->size();

                const DateTimeArray& simulationDates =
                        cache->simulationDates->value();

                // 2) SPOT FLAGS - filters sim dates from sample dates
                // does this per-asset
                BoolArray& observed = cache->observed[asset];
                observed = BoolArray(simulationDates.size(), false);
                // Be careful not to use iStep = 0 because today might be
                // an asset date but it has been dealt with in the past
                for (int i = cache->volStart->value() > today ? 0 : 1;
                     i < simulationDates.size(); ++i) {
                    observed[i] = count(spotAssetDates->begin(), 
                                        spotAssetDates->end(), 
                                        simulationDates[i]) > 0;
                }
            }

            // Create SVGenSpot::IStateVars and put them in database

            vector<double> maxDrifts(numAssets);
            for (int a = 0; a < numAssets; ++a)
                maxDrifts[a] = cache->localVolGrids[a]->value()->maxDrift();

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
        } catch (exception& e){
            throw ModelException(e, __FUNCTION__);
        }
    }

    /** Simulates paths for all assets. Part of the MCPathGen IFace. */
    virtual void generatePath(int pathIdx) {
        // Draw random numbers
        cache->randomGen->generate(pathIdx);
        
        for (int asset = 0; asset < cache->numFactors(); ++asset) {
            generatePath(pathIdx, asset);
        }
    }

    /** Part of the MCPathGen IFace. */
    virtual bool doingPast() const {
        return false;
    }

    bool hasPast() const {
        return simHasPast;
    }

    double maxDriftProduct(int asset) const {
        return cache->localVolGrids[asset]->value()->maxDrift();
    }

    virtual const IMultiFactors* getMultiFactors() const {
        return cache->multiFactors->value();
    }

protected:
    /** Generate path for specified path in simulation (pathIdx), for
        specified asset */
    void generatePath(int pathIdx, int asset) {
        cache->generatePath(pathIdx, asset,
                            &spotProdPaths[asset][spotProdOffset[asset]]);
    }
};

MCPathConfigLV::Gen::Gen(P<LocalVolCache>         cache, 
                         const MCPathGeneratorSP& pastPathGenerator,
                         const MCProductClient*   prodClient,
                         bool                     cachingRequested):
    nowPathIdx(0) {

    try {
        // Collect state variables from product and categorize them
        StateVariableCollectorSP svCollector(new StateVariableCollector());
        prodClient->collectStateVars(svCollector);
        IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

        // * Spot requests
        SVGenSpotArray spotGenArray = filterStateVars<SVGenSpot>(stateVarGenArray);
        if(!spotGenArray.size()) {
            throw ModelException("No spot paths specified.");
        }

        // * Create discount factor state variables
        SVGenDiscFactorArray discFactors = 
            filterStateVars<SVGenDiscFactor>(stateVarGenArray);
        unsigned int iVar;
        for(iVar = 0; iVar < discFactors.size(); iVar++) {
            IStateVariableSP sv(discFactors[iVar]->
                             determinsticSV(false /* not doing past */));
            svDBase.append(discFactors[iVar], sv);
        }

        // * Create expected discount factor (ZCBs) state variables
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

        // Create a spot path generator for the SVGenSpot
        PastPathGen* pastPathGen = dynamic_cast<PastPathGen*>(pastPathGenerator.get());
        if(!pastPathGen) {
            throw ModelException("Past path generator is not of PastPathGen type.");
        }
        pathGenSpot = MCPathConfigLV::PathGenSpotSP(new 
            MCPathConfigLV::PathGenSpot(cache, pastPathGen->getPathGenSpot(),
            prodClient, spotGenArray, cachingRequested, svDBase));

    } catch(exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

int MCPathConfigLV::Gen::NbSimAssets() const {
    throw ModelException(__FUNCTION__, "Method is retired for StateVars");
}

const double* MCPathConfigLV::Gen::Path(int iAsset, int iPath) const {
    throw ModelException(__FUNCTION__, "Method is retired for StateVars");
}; 

double MCPathConfigLV::Gen::refLevel(int iAsset, int iPath) const {
    throw ModelException(__FUNCTION__, "Method is retired for StateVars");
}

double MCPathConfigLV::Gen::maxDriftProduct(int iAsset) const {
    try {
        return pathGenSpot->maxDriftProduct(iAsset);
    } catch(exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

int MCPathConfigLV::Gen::begin(int iAsset) const {
    throw ModelException(__FUNCTION__, "Method is retired for StateVars");
}

int MCPathConfigLV::Gen::end(int iAsset) const{
    throw ModelException(__FUNCTION__, "Method is retired for StateVars");
}

// Live methods
bool MCPathConfigLV::Gen::hasPast() const {
    return pathGenSpot->hasPast();
}

bool MCPathConfigLV::Gen::doingPast() const {
    return false;
}

void MCPathConfigLV::Gen::generatePath(int pathIdx) {
    nowPathIdx = pathIdx;
    pathGenSpot->generatePath(pathIdx);
}

int MCPathConfigLV::Gen::getPathIndex() const {
    return nowPathIdx;
}

/** Returns the state variable corresponding to generator.
    Part of the IStateVariableGen::IStateGen IFace */
IStateVariableSP MCPathConfigLV::Gen::create(const IStateVariableGen* svGen) {

    try {
        return svDBase.find(svGen);
    } catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

/** Creates a future path generator */
MCPathGeneratorSP MCPathConfigLV::makePathGenerator(
    bool                               cachingRequested,
    int                                numPaths,
    const MCPathGeneratorSP&           pastPathGenerator,
    const IMCProduct*                  prod,
    Control*                           control, 
    Results*                           results,
    DateTimeArray&                     simDates){

    static const string method = "MCPathConfigLV::makePathGenerator";
    try { 
        if (!cache) {
            TRACE("New cache");
            cache.reset(new LocalVolCache(prod->getNumAssets()));
        }

        cache->pathConfig->setValue(this, false, false);

        const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);

        cache->tweaked = control && !control->isPricing();

        // Allow some nodes to be updated depending on tweak methodology
        switch (tweakLocalVols) {
            size_t f;
            case TWEAK_LOCAL_VOL:
                // Unlock timeline during Theta
                if (!cache->tweaked) {
                    cache->simulationDates->unlockValue();
                }
                // Unlock LocalVol always
                for(f = 0; f < cache->localVolGrids.size(); f++) {
                    // cache->localVolGrids[f]->unlockValue();
                }
                break;
            case FIX_LOCAL_VOL_EXCEPT_THETA:
                // Unlock timeline and LocalVol during Theta
                if (!cache->tweaked) {
                    cache->simulationDates->unlockValue();
                    for(f = 0; f < cache->localVolGrids.size(); f++) {
                        cache->localVolGrids[f]->unlockValue();
                    }
                }
                break;
            case FIX_LOCAL_VOL:
                // Never unlock timeline
                throw ModelException(method,
                                    BoxedEnum<LVTweaks>::toString(FIX_LOCAL_VOL)+
                                    " is not supported");
                break;
            default:
                throw ModelException("Unrecognized tweak methodology");
                break;
        }
        
        cache->today->setValue(prod->getToday(), cache->tweaked);
        cache->multiFactors->setValue(prod->getMultiFactors(), false, false);

        if (cache->tweaked) {
            if (prod->getMultiFactors()->NbAssets() == 1) {
                for (size_t a = 0; a < cache->factors.size(); ++a) {
                    cache->factors[a]->setValue(a, true, true);
                }
            }

            // Always bin correlations for multiple factors (a less
            // conservative policy would be easy but perhaps a bit brittle)

            cache->correls->pretendChanged();

            if (!cacheLocalVols ||
                !prod->getMultiFactors()->getSensitiveAssets(
                    control->getCurrentSensitivity().get(), true).empty()) {
                // Tweaks affecting ALL assets
                // BE CAREFUL: this includes Correlations but also Theta, ParamSens etc.
                for (int a = 0; a < int(cache->factors.size()); ++a) {
                    cache->factors[a]->setValue(a, true, true);
                }
            }
            else {
                // Tweaks affecting asset by name
                IntArray sensitiveAssets(
                    prod->getMultiFactors()->getSensitiveAssets(
                        control->getCurrentSensitivity().get(), false));

                for (int a = 0; a < int(cache->factors.size()); ++a) {
                    int s = sensitiveAssets.size() - 1;
                    for (; s >= 0 && sensitiveAssets[s] != a; --s);
                    cache->factors[a]->setValue(a, s >= 0, s >= 0);
                }
            }
        }

        cache->tweaked = control && !control->isPricing();

        MCPathGeneratorSP it;
        if (prodClient) {
            it.reset(new MCPathConfigLV::Gen(
                cache, pastPathGenerator, prodClient, cachingRequested));
        }
        else {
            it = MCPathBase::createPathGenerator(
                pastPathGenerator, this, prod,
                refCountPtr<MCPathBaseLW>(
                    new MCPathBaseLW(cache, pastPathGenerator, prod)));
        }

        if (control && control->isPricing()){
            results->storeScalarGreek(cache->simulationDates->value().size() - 1,
                                    Results::DEBUG_PACKET, 
                                    OutputName::SP("LV_SIM_STEPS_USED"));
        }

        // Update some nodes depending on tweak methodology
        switch (tweakLocalVols) {
            size_t f;
            case TWEAK_LOCAL_VOL:
                // Lock timeline if not Theta
                if (!cache->tweaked) {
                    cache->simulationDates->lockValue();
                }
                break;
            case FIX_LOCAL_VOL_EXCEPT_THETA:
                // Lock timeline and LocalVol if not Theta
                if (!cache->tweaked) {
                    cache->simulationDates->lockValue();
                    for(f = 0; f < cache->localVolGrids.size(); f++) {
                        cache->localVolGrids[f]->lockValue();
                    }
                }
                break;
            case FIX_LOCAL_VOL:
                // Lock timeline and LocalVol always
                throw ModelException(method,
                                    BoxedEnum<LVTweaks>::toString(FIX_LOCAL_VOL)+
                                    " is not supported");
                break;
            default:
                throw ModelException("Unrecognized tweak methodology");
                break;
        }
        
        cache->randomGen = cache->newMCRandom();

        return it;
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** MC LV model */
class MonteCarloLVDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloLVDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloLVDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloLV);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloLV(){
        return new MonteCarloLVDefault();
    }
};

CClassConstSP const MonteCarloLVDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloLVDefault", typeid(MonteCarloLVDefault), load);

bool MCPathConfigLVLoad(){
    return MCPathConfigLV::TYPE != 0;
}

DRLIB_END_NAMESPACE
