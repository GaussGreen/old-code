//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRMGen.cpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MCPathConfigSRMGen.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/SRMICE.hpp"
#include "edginc/MemoryProfiler.hpp"
#include "edginc/QMCRNGManager.hpp"
#include "edginc/SRMCreditHJMDiffuse.hpp"
#include "edginc/Maths.hpp"
#include <functional>
#include <typeinfo>

DRLIB_BEGIN_NAMESPACE

SrmCorrData::SrmCorrData(MCPathConfigSRM*   pathConfig) :
    corrObjArray(pathConfig->corrObjArray),
    eqeqCorrObjArray(pathConfig->eqeqCorrObjArray),
    betaCorrArray(pathConfig->betaCorrArray),
    eqeqCorrTermObjArray(pathConfig->eqeqCorrTermObjArray),
    strictCorr(pathConfig->strictCorr)
{
}

MCPathConfigSRMGen::SrmEqVolData::SrmEqVolData(
    const DoubleArrayArray&  assetSpotVols,
    const DoubleArrayArray&  fwdAssetCompVars,
    const DoubleArrayArray&  totalAssetCompVars,
    const DoubleArrayArray&  irAssetIntegral,
    const DoubleArray&       irIrIntegral) :
        assetSpotVols(assetSpotVols),
        fwdAssetCompVars(fwdAssetCompVars),
        totalAssetCompVars(totalAssetCompVars),
        irAssetIntegral(irAssetIntegral),
        irIrIntegral(irIrIntegral)
{}

void MCPathConfigSRMGen::advance() {
    for (size_t i = 0; i < advanceSvSet.size(); ++i)
        advanceSvSet[i]->advance();
}

void MCPathConfigSRMGen::reset() {
    for (size_t i = 0; i < advanceSvSet.size(); ++i)
        advanceSvSet[i]->reset();
}

void MCPathConfigSRMGen::generatePath(int _pathIdx)
{
    pathIdx = _pathIdx;

    rngMgr->seekToPath(pathIdx);
    // loop over all diffusible assets

    pathWeight = fractionByPath[pathIdx];
    for(size_t i=0; i<allAssets.size(); ++i)
    {
        allAssets[i]->generatePath(rngMgr, strataByPath[pathIdx]);
    }

    if (!icePerCcy.empty()){
        for (unsigned int i = 0; i < icePerCcy.size(); ++i){
            // then price our swaptions (a type of control variate)
            // ice => iterative calibration engine
            icePerCcy[i].updateAllPrices();
        }
        if (pathIdx+1 == iceCalibIdx){
            for (unsigned int i = 0; i < icePerCcy.size(); ++i){
                icePerCcy[i].recalibrate(pathIdx);
            }
            // then set the next point to recalibrate
            iceCalibIdx += iceCalibFreq;
        }
    }

    reset();
}

    //// constructor - most of the work is done at construction
MCPathConfigSRMGen::MCPathConfigSRMGen(
    int                      numSimPaths,
    MCPathConfigSRM*         mcPathConfig,
    const MCPathGeneratorSP& pastPathGenerator,
    const MCProductClient*   prodClient,
    DateTimeArray&           simDates, // (O)
    Control*                 control,
    Results*                 results):
        control(control),
        results(results),
        pathIdx(-1),
        iceCalibIdx(numSimPaths/2),  // start half way through
        havePast(pastPathGenerator->hasPast()),
        simDates(simDates),
        betaCorrMax(mcPathConfig->getBetaCorrMax()),
        allAssets(),
        pathWeight(1.0),
        strataByPath(numSimPaths, QMCStrataConstSP()), // by default these are empty strata records
        fractionByPath(numSimPaths, 1.0)
{
    static const string method = "MCPathConfigSRMGen::MCPathConfigSRMGen";
    // eg, 40,000 paths, with numICERuns=2, then run ICE at 20,000 and
    // 30,000. With offCycleICE, this becomes 24,000, and 32,000
    iceCalibFreq = numSimPaths/(2*mcPathConfig->numICERuns +
                                (mcPathConfig->offCycleICE? 1: 0));



    if (mcPathConfig->offCycleICE){
        iceCalibIdx += iceCalibFreq/2;
    }
    try {
        const DateTime& today = prodClient->getToday();
        QMCStratificationRecordArray stratification = mcPathConfig->getStratification();

        /*******************************************************************/
        // Phase I: creating SV object - sorted collections of data relevant for diffusion
        /*******************************************************************/

        sv = SVSP(new SV(mcPathConfig, prodClient)); // build sv->byAssetPosition

        allAssets.reserve(sv->byAssetPosition.size());

        QMCGenDiffusibleIR* pDomIRGen = dynamic_cast<QMCGenDiffusibleIR*>(sv->byAssetPosition[0].get());
        ASSERT(pDomIRGen); // should come first in the vector
        ASSERT(pDomIRGen->name == sv->domISOCode); // sanity check

        // reserve space for fx/eq spot sv -- legacy of MultiAsset object
        vector<SVPathSP> tmp(sv->numFX + sv->numEq, SVPathSP(   ));
        assetPaths.resize(sv->spots.size(), tmp);

         // Setting pointers to the domestic IR object (necessary for init of FX, etc)
        for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
           sv->byAssetPosition[i]->setDomesticIRGen(pDomIRGen);
        
        /*******************************************************************/
        // Phase II: Collecting all the assets in the order of their diffusion,
        // create the SV from SVGens for all the assets
        /*******************************************************************/

        // Creating all the empty-shell IQMCDiffusible******* assets
        // and attaching them to the allAssets vector
        for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
        {
            sv->byAssetPosition[i]->initializeDiffusibleAssets();

            allAssets.push_back( sv->byAssetPosition[i]->getDiffusibleAsset().get() );
        }
        
        /** create state variables (dates are being collected by assets) */
        // also alert SVs to moment matching if it is turned on
        for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
        {
            sv->byAssetPosition[i]->createStateVars(svDBase, assetPaths, mcPathConfig->getMomentMatchingFlag());
        }


        // create all the requested Path Weight SV
        // for the time being - there is one and only -- but can be easily extended
        SVPathWeightSP pwSP(new SVQmcPathWeightTimeIndependent(this));
        if (!stratification.empty())
            QLIB_VERIFY(!sv->pathWeightSVGen.empty(),
                "Monte-Carlo stratified sampling cannot be used for pricers "
                "that do not define PathWeighting State Variables. Please remove the stratification.");

        for(size_t i=0; i<sv->pathWeightSVGen.size(); ++i)
            svDBase.append(sv->pathWeightSVGen[i], pwSP);

        // Propagate MaxMaturity date between assets
//         for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
//         {
//             sv->byAssetPosition[i]->updateAssetMaxDiffusion();
//         }

        /*******************************************************************/
        // Phase III: collecting all the requested dates from all the assets
        //            creating the master timeline (=allDates)
        /*******************************************************************/

        //TODO: might use vector<> and set_union instead - it will be faster by a factor of ln(N_T)
        set<DateTime> wholeDateSet;

        // (1) collecting every date requested by the SV (aka AssetDates)
        for(size_t i=0; i<sv->byAssetPosition.size(); ++i)
        {
            DateTimeArray assetDates = sv->byAssetPosition[i]->getAssetDates();
            wholeDateSet.insert(assetDates.begin(), assetDates.end());
        }
        if (wholeDateSet.empty())
            throw ModelException(method, "Product does not specify any dates");

        wholeDateSet.insert(today);
        
        DateTime lastRequested = *wholeDateSet.rbegin();

        // (2) collecting every 'critical' date  - usually implied volatility points
        for(size_t i=0; i<sv->byAssetPosition.size(); ++i)
        {
            DateTimeArray assetDates = sv->byAssetPosition[i]->
                getCriticalDates(mcPathConfig, today, lastRequested);
            wholeDateSet.insert(assetDates.begin(), assetDates.end());
        }

        // (3) filling in the timeline with more dates to have a reasonably dense line
        DateTimeArray allDates(wholeDateSet.size());
        std::copy(wholeDateSet.begin(), wholeDateSet.end(), allDates.begin());

        simDates = *SRMUtil::fillInTimeLine(today,
                                            allDates,
                                            mcPathConfig->timePointsPerYear);

        assert (simDates[0] == today);

        allDates = DateTime::merge(allDates, simDates);
        // Avoid copying timeline all over
        DateTimeArrayConstSP simDatesSP(new DateTimeArray(simDates)); // FIXME move simDates to be pointer itself
        /** build a Rates util for each ccy, an FX util for each fx,
            a Credit util for each credit and an Equity util for each equity:
            createRatesUtil => calibration of IR spot vols
            createFXUtil/createEquityUtil => bootstrapping of FX/EQ spot vols 
            Moment matching flag is passed to asset Util objects.
            */

        for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
            sv->byAssetPosition[i]->createUtil(mcPathConfig, today, simDatesSP);

        /*******************************************************************/
        // Phase IV: initializing the correlation structure "dependency"
        /*******************************************************************/

        /** now do some clever stuff with the correlations (incl EQEQ corr adjustment)
            build dependence using appropriate dependenceMaker **/
        // in order to have easier access to correlations (needed for CR)
        // TO DO:  Pass in hash of corrs and other mcPathConfig parameters

        corrData = SrmCorrDataSP(new SrmCorrData(mcPathConfig));

        SRMRatesHJMUtil* baseRatesUtil( dynamic_cast<SRMRatesHJMUtil*>(
            sv->byAssetPosition[0]->getRatesUtil().get() ) );

        if ( baseRatesUtil != 0 )
            populateSrmEqVolData(baseRatesUtil); // requires sv->orderedCcys[0]->ratesHJMUtil to be initialized

        DependenceMakerSP depMaker = mcPathConfig->getDependenceMaker();
        dependence = depMaker->createDependence(this);  // first place where assets' randomIdx is used.

        /** then ask for correlations (needed by SRMRates & SRMCredit)
            important: no EQEQ corrs needed and thus its safe to ask for timepoint 0 */
        // Do not ask for correlations as a matrix!

        /** do preparation if ICE adjustment is requested */
        if (mcPathConfig->numICERuns > 0 && !mcPathConfig->deICE)
        {
           for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
            {
                // using ICE, build up swaptions and the state variables they need.
                sv->byAssetPosition[i]->pushbackICERecord(svDBase, icePerCcy);

                // Now we need a sanity check by verifying that no new points were created
                // on the timeline as the result of adding ICE objects
                // if this assertion fails - somehow a new requested date has been
                // added by ICE swaption
                DateTimeArray dta = sv->byAssetPosition[i]->getDiffusibleAsset()->getAssetDates();
                if(!DateTime::isSubset(allDates, dta))
                    throw ModelException(method,
                        "logical error: ICE added few points into timeline which must be already frozen.");
            }
        }
        

        /*******************************************************************/
        // Phase V: getting all the diffusible asset objects properly "set"
        /*******************************************************************/
        for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
        {
            sv->byAssetPosition[i]->setDiffusibleAsset(today, mcPathConfig, sv->pastValues, dependence);
        }

        // rebuild the SV for asset spots (re: Multiassets)
        for (size_t i = 0; i < sv->spots.size(); i++){
            // to do: fix maxDrfits
            vector<double> maxDrifts(assetPaths[i].size(), 100000);
            IStateVariableSP spotSV(new SVGenSpot::PathByAsset(assetPaths[i],
                                                            maxDrifts, false));
                svDBase.append(sv->spots[i], spotSV);
        }
        
        /*******************************************************************/
        // Phase VI: finalizing, preparing, etc...
        /*******************************************************************/

        DateTimeArrayConstSP allDatesSP(new DateTimeArray(allDates));
        // finalize the timeline
        for(size_t i=0; i<allAssets.size(); ++i)
        {

            // this is ugly: TO DO - find a better way
            if (sv->byAssetPosition[i]->isRates())
            {
                StringArray::const_iterator found = std::find(mcPathConfig->rateDriverTargets.begin(), mcPathConfig->rateDriverTargets.end(), sv->byAssetPosition[i]->name);
                if (found != mcPathConfig->rateDriverTargets.end())
                {
                    SRMCreditHJMDiffuse* asset = dynamic_cast<SRMCreditHJMDiffuse*>(allAssets[i]);
                    if (asset != 0)
                        asset->setDriver(mcPathConfig->rateDrivers[found - mcPathConfig->rateDriverTargets.begin()]);
                }
            }
            else if (sv->byAssetPosition[i]->isCredit())
            {
                StringArray::const_iterator found = std::find(mcPathConfig->creditDriverTargets.begin(), mcPathConfig->creditDriverTargets.end(), sv->byAssetPosition[i]->name);
                if (found != mcPathConfig->creditDriverTargets.end())
                {
                    SRMRatesHJMDiffuse* asset = dynamic_cast<SRMRatesHJMDiffuse*>(allAssets[i]);
                    if (asset != 0)
                        asset->setDriver(mcPathConfig->creditDrivers[found - mcPathConfig->creditDriverTargets.begin()]);
                }
            }
            else if (sv->byAssetPosition[i]->isEnrg())
            {
                StringArray::const_iterator found = std::find(mcPathConfig->energyDriverTargets.begin(), mcPathConfig->energyDriverTargets.end(), sv->byAssetPosition[i]->name);
                if (found != mcPathConfig->creditDriverTargets.end())
                {
                    if (dynamic_cast<SRMEnergyOilDiffuse*>(allAssets[i]) != 0)
                        dynamic_cast<SRMEnergyOilDiffuse*>(allAssets[i])->setDriver(mcPathConfig->energyDrivers[found - mcPathConfig->energyDriverTargets.begin()]);
                    else if (dynamic_cast<SRMEnergyGasDiffuse*>(allAssets[i]) != 0)
                        dynamic_cast<SRMEnergyGasDiffuse*>(allAssets[i])->setDriver(mcPathConfig->energyDrivers[found - mcPathConfig->energyDriverTargets.begin()]);
                }
            }

            allAssets[i]->finalize(allDatesSP);
        }

        // prepare all the state variables
        for(StateVarDBase::DBase::iterator it=svDBase.svDBase.begin(); it !=svDBase.svDBase.end(); ++it)
            it->second->prepare(mcPathConfig->getMomentMatchingFlag());


        /*******************************************************************/
        // Phase VII: stratified sampling
        /*******************************************************************/
        // first of all we should normalize the weights and calculate the 
        // integer number of paths to be run for each strata, 
        // i.e. what strata each path will be running

        if (!stratification.empty())
        {
            // below also is verification that the stratas form a non-intersecting 
            // complete cover of all the allowed states

            double total = 0.0;
            for(int j=0; j<stratification.size(); ++j)
                total += stratification[j]->getWeight();
            QLIB_VERIFY(total > 0, "Stratification cannot be done -- stratas have weight 0.0");

            int numUsed = 0;
            double totalProbability = 0;

            for(int j=0; j<stratification.size(); ++j)
            {

                // (1) verification of independence of stratas
                for(int k=j+1; k<stratification.size(); ++k)
                    QLIB_VERIFY(stratification[j]->getStrata()->notIntersecting(stratification[k]->getStrata()),
                        "The stratas " + Format::toString(j) + " and "+Format::toString(k) + " have a non-empty intersect. Correct the issue.");


                int numStrataPaths;

                double realProbability = 1.0;
                for(size_t i = 0; i<sv->byAssetPosition.size(); ++i)
                    realProbability *= sv->byAssetPosition[i]->
                        getDiffusibleAsset()->getStrataProbability(stratification[j]->getStrata());

                // (2) we will check the the probabilities sum up to 100%
                totalProbability += realProbability;

                if (j+1 < stratification.size()) // at least one path is going to happen for each strata

                    numStrataPaths = max(int(numSimPaths * stratification[j]->getWeight() / total + 0.5), 1);

                else
                {
                    numStrataPaths = numSimPaths - numUsed; // special treatment for the last slot
                    QLIB_VERIFY(numStrataPaths > 0, 
                        "Due to rounding, the last strata had to get a negative number of paths -- please conside moving the heaviest strata to the last position");
                }

                for(int i=numUsed; i<numUsed+numStrataPaths; ++i)
                {
                    strataByPath[i] = stratification[j]->getStrata();
                    fractionByPath[i] = numSimPaths * realProbability/ numStrataPaths;  // Radon-Nikodym
                }
                numUsed += numStrataPaths;
            }
            // we might want to see if isZero gives enough leeway for roundup/down errors
            QLIB_VERIFY(Maths::isZero(totalProbability - 1.0), 
                "The stratification array does not cover all the possible "
                "scenarios -- please add more stratas to complete the coverage of all the possibilities.");
        }


        /*******************************************************************/
        // Phase VI: Random Number generators
        /*******************************************************************/
        initializeRandGen(numSimPaths, mcPathConfig);

        /*******************************************************************/
        // Phase V: do moment matching
        /*******************************************************************/
        if (mcPathConfig->getMomentMatchingFlag())
        {
            enforceFirstMomentMatching(numSimPaths);
            // re-initialize random number generators
            initializeRandGen(numSimPaths, mcPathConfig);
        }

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

void MCPathConfigSRMGen::initializeRandGen(int numSimPaths, MCPathConfigSRM* mcPathConfig)
{
    ///// this piece is really hacky ...... //////////
    int numFutSteps = simDates.size()-1, numPastDates = 0;
    // Initialize random number generator. For matchLegacySRM3 we need to
    // generate the randoms in a different order which essentially means
    // we need a cache
    MCPathConfig::RandomCacheSP myRandomCache(
        mcPathConfig->matchLegacySRM3 ?
        new MCPathConfig::RandomCache(
        numSimPaths/2,
        false, /* correlated */
        true) : // uncorrelated
    new MCPathConfig::RandomCache());
    myRandomCache->configure(dependence->getNumAssets(), numFutSteps);
    randomGen = MCRandomSP(
        new MCRandom(
        0,
        dependence,
        mcPathConfig->getRandomGenerator(),
        myRandomCache,
        mcPathConfig->carefulRandoms(),  // major problems!
        numFutSteps,
        dependence->getNumAssets(), // bad bad fudge! TODO
        numPastDates));
    if (mcPathConfig->matchLegacySRM3){
        int runSize = numSimPaths/(2 * mcPathConfig->numICERuns +
            (mcPathConfig->offCycleICE? 1: 0));
        randomGen->blockFillByDate(runSize/2);
    }

    rngMgr = IQMCRNGManagerSP(new QMCRNGManager(randomGen, dependence->getNumAssets()));

    // finally update our ICE SV's
    for (unsigned int idx = 0; idx < icePerCcy.size(); idx++){
        icePerCcy[idx].pathGenUpdated(this);
    }
    // now that svDBase is populated, store all advanceables
    svDBase.filterAdvanceables(advanceSvSet);
}



void MCPathConfigSRMGen::enforceFirstMomentMatching(int numSimPaths)
{
     size_t i = 0;
     vector<IQMCStateVariableBase*> allIQMCStateVars;
 
     for (StateVarDBase::DBase::iterator it=svDBase.svDBase.begin(); it != svDBase.svDBase.end(); ++it)
     {
         IQMCStateVariableBase* ptr = dynamic_cast<IQMCStateVariableBase*>(it->second.get());
         if (ptr) 
         {
             allIQMCStateVars.push_back(ptr);
             continue;
         } 

         MCPath::PathByAsset* pathPtr = dynamic_cast<MCPath::PathByAsset*>(it->second.get());
         if (pathPtr)
         {
             for(int i=0; i<pathPtr->numAssets(); ++i)
             {
                 IQMCStateVariableBase* ptr = dynamic_cast<IQMCStateVariableBase*>(pathPtr->getSV(i).get());
                 if (ptr)
                    allIQMCStateVars.push_back(ptr);
             }
             continue;
         }
     }

     for (i = 0; i < allIQMCStateVars.size(); ++i)
     {
         allIQMCStateVars[i]->resetMMCorrection();
     }

    
   
     for (int p = 0; p < numSimPaths; ++p)
     {
         generatePath(p);
         for (i = 0; i < allIQMCStateVars.size(); ++i)
         {
             allIQMCStateVars[i]->accumulateMMCorrection();
         }
         
     }
     for (i = 0; i < allIQMCStateVars.size(); ++i)
     {
         allIQMCStateVars[i]->setMMCorrection();
     }
    
     
}


/** To get Beta correlation structure.
 *  Important:  beta correlation can presently work only between single-factor objects!
 */
SparseDoubleMatrixSP MCPathConfigSRMGen::getBetaCorrelations(
    const DependenceMakerGaussSrm* dependenceMaker) const
{
    static const string method = "MCPathConfigSRMGen::getBetaCorrelations";
    try {
        SparseDoubleMatrixSP betas;
        //  retrieve off-diagonal betas and map then to the indices properly

        // The lookup is going to be reverse compared to createMasterSparseCollection()
        // we will go over the array of beta correlations and search the
        // diffusible assets for all the pairs we met there


        // Let's see if this part is required at all:
        if (corrData->betaCorrArray.empty())
        {

            // DEBUG: trying to see quickly whether the order is correct:
            // imagine p1 corresponds to T1 asset, p1 = 0
            // p2 corresponds to T2 asset, p2 = 3
            // as this experiment shows, correct way of initializing is
            //  push_back(tier1, tier2, beta)
//            int numFactors = sv->numAllFactors();
//            if (numFactors > 3)
//            {
//               SparseCollectionSP sparseCollection(new SparseCollection(numFactors, numFactors));
//               sparseCollection->push_back(0, 3, 0.5555);
//               sparseCollection->push_back(1, 3, 0.3333);
//
//               betas = SparseDoubleMatrixSP(new SparseDoubleMatrix(*sparseCollection));
//            }

            return betas;
        }

        // (1) let's prepare for the quick search -- build a map of
        //      asset name -> position in sv->byAssetPosition


        int numAssets = sv->byAssetPosition.size();
        int numFactors = sv->numAllFactors();
        map<string, int> assetMapper;
        for (int r = 0; r < numAssets; r++){
            assert(sv->byAssetPosition[r]->isDefined());
            string name = sv->byAssetPosition[r]->name;
            // no asset name should exist twice the list
            assert(assetMapper.find(name)==assetMapper.end());
            assetMapper[name] = r;
        }

        BetaCorrelationArray &betaCorrs = corrData->betaCorrArray;
        SparseCollectionSP sparseCollection(new SparseCollection(numFactors, numFactors));

        for(int i=0; i<betaCorrs.size(); ++i)
        {
            BetaCorrelationSP betaCorr = betaCorrs[i];
            string name1 = betaCorr->getTier1Name();
            string name2 = betaCorr->getTier2Name();
            double corr = betaCorr->getBetaCorrelation();

            // TODO: make sure that order (name1, name2) is correct
            if (assetMapper.count(name1) && assetMapper.count(name2)) // if both are present
            {
                // make sure that there is no multifactor objects given for beta-corr,
                // as this is not implemented yet
                int pos1 = assetMapper[name1];
                int pos2 = assetMapper[name2];


                if (sv->byAssetPosition[pos2]->numFactors() >1)
                    throw ModelException(method,
                        "beta correlation is supported only for single-factor Tier2 objects for now alas "
                        +name2+" is multifactor");


                int tier1idx = sv->byAssetPosition[pos1]->getRandomIdx();
                int tier2idx = sv->byAssetPosition[pos2]->getRandomIdx();

                // Expand betas if needed for multi-factor tier1 assets.
                vector<double> tmpCorrs( sv->byAssetPosition[pos1]->expandTier2Betas( corr ) );

                for(size_t i=0; i<tmpCorrs.size(); ++i)
                    sparseCollection->push_back(tier1idx + i, tier2idx, tmpCorrs[i]);

            }
        }
        betas = SparseDoubleMatrixSP(new SparseDoubleMatrix(*sparseCollection));

        return betas;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** DEPENDENCE MAKER GAUSS SRM && NO CORR MAPPING */
vector<SparseDoubleMatrixSP> MCPathConfigSRMGen::createSparseGaussMatrixArray(
                                    const DependenceMakerGaussSrm* dependenceMaker) const {
    static const string method = "MCPathConfigSRMGen::createSparseGaussMatrixArray";
    try {        
        /** create masterSparseCollection */
        SparseCollectionSP masterSparseCollection = createMasterSparseCollection();
        
        /** convert masterSparseCollection into sparse matrix */
        SparseDoubleMatrixSP masterSparseCorrMatrix =
            SparseDoubleMatrixSP(new SparseDoubleMatrix(*masterSparseCollection));

        vector<SparseDoubleMatrixSP> sqrtCorrMatrixArray(1);
        sqrtCorrMatrixArray[0] = // length one, since neither corrTS nor corrMapping
            checkAndExtendSparseCorrMatrix(masterSparseCorrMatrix,
                                           dependenceMaker->getEigenValueFloor(),
                                           dependenceMaker->getMaxSqError(),
                                           simDates.front(), // helper date for error msg
                                           simDates.back()); // helper date for error msg
        return sqrtCorrMatrixArray; 
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** DEPENDENCE MAKER GAUSS SRM && YES CORR MAPPING */
vector<SparseDoubleMatrixSP> MCPathConfigSRMGen::createSparseGaussMatrixArray(
    const DependenceMakerGaussSrm*  dependenceMaker,
    const IntArray&                 fwdCorrIndexArray,
    const DateTimeArray&            fwdCorrDatesArray) const {
    static const string method = "MCPathConfigSRMGen::createSparseGaussMatrixArray";
    try {
                
        /** get raw eqeq correlations */
        CDoubleMatrixSP eqEqCorrelations = getOrderedEqEqCorrelations();        

        /* create vector of (constant) correlation matrices */
        DoubleMatrixArraySP fwdEqEqCorrelations(new DoubleMatrixArray(fwdCorrDatesArray.size()));
        for (int iStep = 0; iStep < fwdCorrDatesArray.size(); iStep++) {
            (*fwdEqEqCorrelations)[iStep] = CDoubleMatrixSP(new DoubleMatrix(*eqEqCorrelations));
        }

        /** map fwd correlations into spot correlations */
        computeAdjEqEqCorrs(fwdEqEqCorrelations, 
                            fwdCorrIndexArray, 
                            fwdCorrDatesArray);        

        /** create array of sparse correlation matrices */
        return replaceOrderedEqEqCorrelations(eqEqCorrelations,
                                              fwdEqEqCorrelations,
                                              dependenceMaker->getEigenValueFloor(),
                                              dependenceMaker->getMaxSqError()); 
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** DEPENDENCE MAKER GAUSS TERM SRM && NO/YES CORR MAPPING */
DoubleMatrix MCPathConfigSRMGen::getFwdVarAtDates() const {
    static const string method = "MCPathConfigSRMGen::getFwdVarAtDates";
    try {
        /** SrmEqVolData.fwdAssetCompVars is array of length nbAssets with arrays of length nbDates
            fwdVarAtDates is DoubleMatrix with cols=nbAssets and rows=nbDates */
        return CDoubleMatrix(thisSrmEqVolData->fwdAssetCompVars);
    } catch (exception& e){
        throw ModelException(e, method);
    }

}

/** DEPENDENCE MAKER GAUSS TERM SRM && NO/YES CORR MAPPING */
void MCPathConfigSRMGen::getCorrelationData(
    CorrelationCommonArray& corrObjArray,
    CorrelationTermArray&   corrTermObjArray) const {
    static const string method = "MCPathConfigSRMGen::getCorrelationData";
    try {
        corrObjArray.resize(corrData->eqeqCorrObjArray.size());
        corrTermObjArray.resize(corrData->eqeqCorrTermObjArray.size());
        int pos = 0;
        for (int iAsset=0; iAsset<sv->numEq; iAsset++) {
            string name1 = sv->orderedEquities[iAsset]->eqAsset->getTrueName();
#if 0
            /*
             *  Note that if post-calibration (aka low-level) correlation may be input,
             *  it is possible that we have non-trivial correlation for, e. g., 
             *  (USD, USD) as a matrix
             */
            for (int jAsset=iAsset; jAsset<sv->numEq; jAsset++)
            {
                string name2 = sv->orderedEquities[jAsset]->eqAsset->getTrueName();
                int index;
                try
                {
                    index = CorrelationCommon::findCorrArrayIndex(
                        corrData->eqeqCorrObjArray, name1, name2);
                }
                catch (ModelException& e)
                {
                    continue;   // ignore if not found
                }
                corrObjArray[pos] = corrData->eqeqCorrObjArray[index];
                corrTermObjArray[pos] = corrData->eqeqCorrTermObjArray[index];
                ++ pos;
            }
#else   
            /*
             *  Yet, the above seems to break a regression test on CorrTerm by an amount of order 10^-4
             */
            for (int jAsset=iAsset + 1; jAsset<sv->numEq; jAsset++)
            {
                string name2 = sv->orderedEquities[jAsset]->eqAsset->getTrueName();
                int index;
                try
                {
                    index = CorrelationCommon::findCorrArrayIndex(
                        corrData->eqeqCorrObjArray, name1, name2);
                }
                catch (ModelException& /*e*/)
                {
                    continue;   // ignore if not found
                }
                corrObjArray[pos] = corrData->eqeqCorrObjArray[index];
                corrTermObjArray[pos] = corrData->eqeqCorrTermObjArray[index];
                ++ pos;
            }
#endif
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

    /** DEPENDENCE MAKER GAUSS TERM SRM && NO/YES CORR MAPPING */
vector<SparseDoubleMatrixSP> MCPathConfigSRMGen::createSparseGaussTermMatrixArray(
    const DependenceMakerGaussSrm*  dependenceMaker,
    DoubleMatrixArraySP                 fwdCorrelations,
    const IntArray&                     fwdCorrIndexArray,
    const DateTimeArray&                fwdCorrDatesArray) const {
    static const string method = "MCPathConfigSRMGen::createSparseGaussTermMatrixArray";
    try {

       /** NO CORR MAPPING */ 
        if (!(dependenceMaker->doCorrMapping())) {

            return replaceOrderedEqEqCorrelations(getOrderedEqEqCorrelations(),
                                                  fwdCorrelations,
                                                  dependenceMaker->getEigenValueFloor(),
                                                  dependenceMaker->getMaxSqError());
        /** YES CORR MAPPING */
        } else {            
            /** do corr mapping */
            computeAdjEqEqCorrs(fwdCorrelations, 
                                fwdCorrIndexArray,
                                fwdCorrDatesArray);

            /** create array of sparse correlation matrices */
            return replaceOrderedEqEqCorrelations(getOrderedEqEqCorrelations(),
                                                  fwdCorrelations,
                                                  dependenceMaker->getEigenValueFloor(),
                                                  dependenceMaker->getMaxSqError());
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}


/** creation of master sparse collection - note that GaussTerm
    does not work with multifactor IR.
    => CALLED ONLY ONCE, since highly inefficient with all the lookups */
SparseCollectionSP MCPathConfigSRMGen::createMasterSparseCollection() const {
    static const string method = "MCPathConfigSRMGen::createMasterSparseCollection";
    try {
        int numAssets = sv->byAssetPosition.size();
        int totalNumFactors = sv->numAllFactors();

        SparseCollectionSP sparseCollection(new SparseCollection);
        for (int p = 0; p < totalNumFactors; p++){
            sparseCollection->push_back(p,p,1.0); // ensure diagonal is 1
        }
        for (int r = 0; r < numAssets; ++r){

            // special treatment for assets not using brownians -- just skip them!
            if (!sv->byAssetPosition[r]->numFactors()) continue; // do nothing

            QLIB_VERIFY(sv->byAssetPosition[r]->isDefined(),
                "Asset at the position "+Format::toString(r)+ 
                " with name " + sv->byAssetPosition[r]->name + 
                " which is required for correlation -- is not defined. Seems like internal ERROR. ");

            string name1 = sv->byAssetPosition[r]->name;
            int numRowsToFill = sv->byAssetPosition[r]->numFactors();

#if 0
            /*
             *  Note that if post-calibration (aka low-level) correlation may be input,
             *  it is possible that we have non-trivial correlation for, e. g., 
             *  (USD, USD) as a matrix
             */
            for (int c = 0; c <= r; c++)
#else
            /*
             *  Yet, the above seems to break a regression test on CorrTerm by an amount of order 10^-4
             */
            for (int c = 0; c < r; c++)
#endif
            {
                // special treatment for assets not using brownians -- just skip them!
                if (!sv->byAssetPosition[c]->numFactors()) continue; // do nothing

                QLIB_VERIFY(sv->byAssetPosition[c]->isDefined(),
                    "Asset at the position "+Format::toString(c)+ 
                    " with name " + sv->byAssetPosition[c]->name + 
                    " which is required for correlation -- is not defined. Seems like internal ERROR. ");

                string name2 = sv->byAssetPosition[c]->name;

                int numColsToFill = sv->byAssetPosition[c]->numFactors();
                CorrelationCommonSP correlation = CorrelationCommon::lookupCorrelation(corrData->corrObjArray,
                    name1, name2,
                    corrData->strictCorr);
                if (correlation.get() == 0)
                        continue;

                static CClassConstSP FactorCorrelationClazz = CClass::forName("FactorCorrelation");
                if (correlation.get() && FactorCorrelationClazz->isAssignableFrom(correlation->getClass())) 
                {
                    // factor correlations found, skip factor calibration:
                    FactorCorrelationSP factorCorr = FactorCorrelationSP::dynamicCast(correlation);
                    DoubleMatrix matrix = factorCorr->getCorrelation();
                    if (correlation->getAsset1Name() != name1) 
                        matrix.transpose();

                    if (numRowsToFill != matrix.numRows() || numColsToFill != matrix.numCols())
                        throw ModelException(method, "size of factor correlation matrix and number of factors mismatch"); 

                    int rowidx = sv->byAssetPosition[r]->getRandomIdx();
                    int colidx = sv->byAssetPosition[c]->getRandomIdx();

                    for (int m=0; m<numRowsToFill; ++m)
                    {
                        // here we need to be careful not including diagonal twice!
                        if (name1 == name2)
                        {
                            ASSERT(Maths::isZero(matrix[m][m]-1.0));
                            for (int n=m+1; n<numColsToFill; ++n)
                            {
                                if (matrix[n][m] != 0.0) // TODO:  Should use Maths::isZero(...)
                                {
                                    sparseCollection->push_back(rowidx+m, colidx+n, matrix[n][m]);
                                    sparseCollection->push_back(colidx+n, rowidx+m, matrix[n][m]);
                                    QLIB_VERIFY(fabs(matrix[n][m]) <= 1.0, "ERROR: Correlation > 1");
                                }
                            }
                        }
                        else
                        {
                            for (int n=0; n<numColsToFill; ++n)
                            {
                                if (matrix[n][m] != 0.0)
                                {
                                    QLIB_VERIFY(rowidx+m != colidx + n, "Writing on the diagonal! But 1 was already put there!");
                                    sparseCollection->push_back(rowidx+m, colidx+n, matrix[n][m]);
                                    sparseCollection->push_back(colidx+n, rowidx+m, matrix[n][m]);
                                    QLIB_VERIFY(fabs(matrix[n][m]) <= 1.0, "ERROR: Correlation > 1");

                                }
                           }

                        }
                   }
                   continue;
                }

                /** look up correlation (may prove to be too inefficient (for non-adj & adj corr)
                    if no adj eqeq corr can be found, then take market corr
                    this applies to all nonEQEQ corr and to EQEQ corr, if no mapping */
                // TODO:  this seems to be a bug which causes expansion of diagonal correlations to 
                // be skipped resulting in 0 off-diagonal factor correlations.
                if (name1 == name2)
                    continue;   // it will be raw correlation of 1., no need to include

                double corr = correlation.get() ? CorrelationSP::dynamicCast(correlation)->getCorrelation() : 0;

                vector<double> tmpCorrs(
                    sv->byAssetPosition[r]->expandAssetAssetCorrs(sv->byAssetPosition[c].get(), corr));

                int rowidx = sv->byAssetPosition[r]->getRandomIdx();
                int colidx = sv->byAssetPosition[c]->getRandomIdx();

                if (tmpCorrs.size())    //may be zero: if rates are deterministic then there are zero factors.
                {
                    /* fill in the extended corr matrix */
                    int x = 0;
                    for (int m = 0; m < numRowsToFill; m++){
                        for (int n = 0; n < numColsToFill; n++, x++)
                        if (tmpCorrs[x] != 0.0) {
                            QLIB_VERIFY(rowidx+m != colidx+n, "Writing on the diagonal! But 1 was already put there!");
                            
                            sparseCollection->push_back(rowidx+m, colidx+n, tmpCorrs[x]);
                            sparseCollection->push_back(colidx+n, rowidx+m, tmpCorrs[x]);
                            QLIB_VERIFY(fabs(tmpCorrs[x]) <= 1.0, "ERROR: Correlation > 1");

                        }
                    }
                }
            }
        }

        /** if output request, the check pos definiteness and report potential modification */
        if (control) {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::MOD_CORR_MATRIX_SQ_ERROR);
            if(request) {
                SparseDoubleMatrixSP sparseCorrMatrix(new SparseDoubleMatrix(*sparseCollection));
                double sqError = 0.0;
                SparseDoubleMatrix sparseCorrMatrixTmp =
                    sparseCorrMatrix->symmToCorrel(&sqError, CorrelationTerm::EIGEN_VALUE_FLOOR);
                results->storeRequestResult(request,sqError);
            }
        }
        return sparseCollection;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

SparseDoubleMatrixSP MCPathConfigSRMGen::checkAndExtendSparseCorrMatrix(
                                                SparseDoubleMatrixSP    inputSparseCorrMatrix,
                                                double                  eigenValueFloor,
                                                double                  maxSqError,
                                                DateTime                startDate,
                                                DateTime                endDate) const {
    static const string method = "MCPathConfigSRMGen::checkAndExtendSparseCorrMatrix";
    try {

        /** check that matrix without extension for IR factors is (almost) pos definite */
        double sqError = 0.0;

        SparseDoubleMatrix outputSparseCorrMatrix =
            inputSparseCorrMatrix->symmToCorrel(&sqError, eigenValueFloor);

        if( Maths::isPositive(sqError - maxSqError)) {
            throw ModelException(method,
                "Average distance to positive definite modification of correlation matrix is "
                    + Format::toString(sqrt(sqError)) + ".\n"
                    + "Average distance should not exceed the square of "
                    + Format::toString(sqrt(maxSqError)) + ". \n"
                    + "The error occurs between " + startDate.toString() 
                    + " and " + endDate.toString() + ".");
        }

        /** now we can for sure do the square root */
        outputSparseCorrMatrix = outputSparseCorrMatrix.computeSquareRoot();

        /** extend matrix for correlation between the different factors
            convert the transformation matrix to be based on the IR exponential factors
            this 2-step sqrt'ing allows to assign weights to the important IR factors */

        for (size_t x=0; x<sv->byAssetPosition.size(); ++x)
        {
            sv->byAssetPosition[x]->expandFactors(outputSparseCorrMatrix);
        }

        /** return modified output sparse matrix */
        return SparseDoubleMatrixSP(new SparseDoubleMatrix(outputSparseCorrMatrix));

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** retrieve raw EqEq corrmatrix from corr obj array via string search */
CDoubleMatrixSP MCPathConfigSRMGen::getOrderedEqEqCorrelations() const {
    static const string method = "MCPathConfigSRMGen::getOrderedEqEqCorrelations";
    try {
        int nbEqAssets = sv->numEq;
        CDoubleMatrix eqEqCorrMatrix(nbEqAssets,nbEqAssets);
        for(int iAsset = 0; iAsset < nbEqAssets; iAsset++) {
            eqEqCorrMatrix[iAsset][iAsset] = 1.0;
            string name1 = sv->orderedEquities[iAsset]->eqAsset->getTrueName();
            // TO DO: may need to consider self correlation if post calibration correlation 
            // is possible here, for Eq-Eq.
            for (int jAsset = iAsset + 1; jAsset < nbEqAssets; jAsset++) {
                string name2 = sv->orderedEquities[jAsset]->eqAsset->getTrueName();
                double eqeqCorr;
                CorrelationCommonSP eqeqCorrSP = CorrelationCommon::find(corrData->eqeqCorrObjArray, name1, name2);
                if (eqeqCorrSP.get() && Correlation::TYPE->isInstance(eqeqCorrSP.get()))
                {
                    CorrelationSP p = CorrelationSP::dynamicCast(eqeqCorrSP);
                    eqeqCorr = p.get() ? p->getCorrelation() : 0;
                     // TO DO: handle the case the correlation is post calibrated
               }
                else
                    eqeqCorr = 0;   // TO DO: what to do if it is post calibrated?
                                    // Presumably Eq-Eq never have multi factor.
                eqEqCorrMatrix[iAsset][jAsset] = eqeqCorr;
                eqEqCorrMatrix[jAsset][iAsset] = eqeqCorr;
            }
        }
        return CDoubleMatrixSP(new CDoubleMatrix(eqEqCorrMatrix));
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** for CorrTS as well as CorrMap */
vector<SparseDoubleMatrixSP> MCPathConfigSRMGen::replaceOrderedEqEqCorrelations(
        CDoubleMatrixSP     eqeqCorrelations,
        DoubleMatrixArraySP adjEqEqCorrelations,
        double              eigenValueFloor,
        double              maxSqError) const {
    static const string method = "MCPathConfigSRMGen::replaceOrderedEqEqCorrelations";
    try {
        SparseCollectionSP masterSparseCollection = createMasterSparseCollection();
        vector<SparseDoubleMatrixSP> modMasterSparseCorrMatrixArray(simDates.size()-1);

        int totalModifications = sv->numEq * (sv->numEq - 1);
        for (int iStep = 0; iStep < adjEqEqCorrelations->size(); iStep++) {
            CDoubleMatrixSP thisAdjEqEqCorrelations = (*adjEqEqCorrelations)[iStep];

            /** for DEBUG purposes */
            //SparseDoubleMatrixSP tmp1a(new SparseDoubleMatrix(*masterSparseCollection));
            //CDoubleMatrixSP tmp1b = tmp1a->toDoubleMatrixSP();

            /** COPY master sparse collection and reserve space */
            SparseCollection currentSparseCollection = *masterSparseCollection;
            currentSparseCollection.reserve(currentSparseCollection.size()+totalModifications);

            for (int iAsset = 0; iAsset < sv->numEq; iAsset++) {
                double* thisCorrVec = (*thisAdjEqEqCorrelations)[iAsset];
                double* thisOldCorrVec = (*eqeqCorrelations)[iAsset];
                int idx1 = sv->orderedEquities[iAsset]->getRandomIdx();
                for (int jAsset = iAsset + 1; jAsset < sv->numEq; jAsset++) {
                    int idx2 = sv->orderedEquities[jAsset]->getRandomIdx();
                    double corrDiff = thisCorrVec[jAsset] - thisOldCorrVec[jAsset];
                    currentSparseCollection.push_back(idx1, idx2, corrDiff);
                    currentSparseCollection.push_back(idx2, idx1, corrDiff);
                }
            }
            /** check, extend & save */
            modMasterSparseCorrMatrixArray[iStep] = checkAndExtendSparseCorrMatrix(
                SparseDoubleMatrixSP(new SparseDoubleMatrix(currentSparseCollection)),
                eigenValueFloor, maxSqError, simDates[iStep], simDates[iStep+1]);

            /** for DEBUG purposes */
            //SparseDoubleMatrixSP tmp2a(new SparseDoubleMatrix(*currentSparseCollection));
            //CDoubleMatrixSP tmp2b = tmp2a->toDoubleMatrixSP();
        }
        return modMasterSparseCorrMatrixArray;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** SrmEqVolData needed for CorrTS and/or CorrMapping */
void MCPathConfigSRMGen::populateSrmEqVolData( SRMRatesHJMUtil* baseRatesUtil ) {
    static const string method = "MCPathConfigSRMGen::populateSrmEqVolData";
    try {
        /** general info */
        int iStep, jStep, nbSteps = simDates.size() - 1; // simDates.size() = nbSteps + 1
        int iAsset, nbEqAssets = sv->numEq;
        unsigned int ll;

        /** allocate some space for the objects we need to save ... */
        DoubleArrayArray    assetSpotVols(nbEqAssets);
        DoubleArrayArray    fwdAssetCompVars(nbEqAssets);
        DoubleArrayArray    irAssetIntegral(nbEqAssets);
        DoubleArray         irIrIntegral(nbSteps);

        /** allocate some space of the objects we do not need to save ... */
        vector< vector<double> >    irAssetCorrs(nbEqAssets);
        DoubleArrayArray            intSqAssetSpotVols(nbEqAssets);  // int_0^t sigma^2(u) du
        DoubleArrayArray            totalAssetCompVars(nbEqAssets);
        for(iAsset = 0; iAsset < nbEqAssets; iAsset++) {
            QMCGenDiffusibleEquitySP eqData = sv->orderedEquities[iAsset];
            vector<double> tmp = eqData->getEquityUtil()->getSpotVol();
            assetSpotVols[iAsset] = DoubleArray(tmp.begin(), tmp.end());
            string name = eqData->eqAsset->getTrueName();

            /** also retrieve corr betw IR and EQ */
            double corrValue;
            CorrelationCommonSP corrSP = CorrelationCommon::lookupCorrelation(corrData->corrObjArray,
                                                                         name,  baseRatesUtil->getCcy(),
                                                                         corrData->strictCorr);
            if (corrSP.get() && Correlation::TYPE->isInstance(corrSP.get()))
            {
                CorrelationSP p = CorrelationSP::dynamicCast(corrSP);
                corrValue = p.get() ? p->getCorrelation() : 0;
            }
            else
                corrValue = 0;

            irAssetCorrs[iAsset] = SRMCorrelation::assetIr(
                    corrValue,
                    *baseRatesUtil,
                    SRMCorrelation::EXPONENTIAL_FACTORS);
            irAssetIntegral[iAsset].resize(nbSteps);
            intSqAssetSpotVols[iAsset].resize(nbSteps);
            fwdAssetCompVars[iAsset].resize(nbSteps);
            totalAssetCompVars[iAsset].resize(nbSteps);
        }

        /** more IR info: ir spot vols and ir fwd rate */
        DoubleArray deltaTime(nbSteps);
        vector<double> tmp = baseRatesUtil->getSwaptionSpotVolsAtSimDates();
        DoubleArray baseIrSpotVol = DoubleArray(tmp.begin(), tmp.end());
        DoubleArray baseIrFwdRateDeltaT(nbSteps);
        IYieldCurveConstSP yc = baseRatesUtil->getDiffYC(); // assumes flat extrapolation
        double Discount1 = yc->pv(simDates[0]);
        for (iStep=0; iStep<nbSteps; iStep++) {
            deltaTime[iStep] = simDates[iStep].yearFrac(simDates[iStep+1]);
            double Discount2 = yc->pv(simDates[iStep+1]);
            baseIrFwdRateDeltaT[iStep] = Discount1 / Discount2 - 1.0;
            /** in order to avoid zero fwdRates */
            if (Maths::isZero(baseIrFwdRateDeltaT[iStep])
                    && iStep>0
                    && Maths::isPositive(deltaTime[iStep-1])) {
                baseIrFwdRateDeltaT[iStep] =
                    baseIrFwdRateDeltaT[iStep-1] / deltaTime[iStep-1] * deltaTime[iStep];
            }
            Discount1 = Discount2;
        }

        int numIrFactors = baseRatesUtil->numFactors();
        int nbIrIrIntegrals = numIrFactors + (numIrFactors - 1) * numIrFactors / 2;

        for (iStep = 0; iStep < nbSteps; iStep++) {
            /** allocate space for IR and each EQ */
            vector<double>              aFactor(numIrFactors);
            vector<double>              irIrIntHelper(nbIrIrIntegrals);
            vector< vector<double> >    irAssetIntHelper(nbEqAssets);
            for(iAsset = 0; iAsset < nbEqAssets; iAsset++) {
                irAssetIntHelper[iAsset].resize(numIrFactors);
            }
            
            ASSERT(baseIrFwdRateDeltaT.size() > iStep);
            ASSERT(baseIrSpotVol.size() > iStep);
            ASSERT(deltaTime.size() > iStep);
            
            /** loop 1: go backwards ... */
            for (jStep = iStep; jStep>=0; jStep--) {
                /** first: pure ir integral*/
                baseRatesUtil->aFactor(baseIrFwdRateDeltaT[jStep],
                                    deltaTime[jStep],
                                    aFactor);
                /** irIrIntHelper += sigma_P^2(jStep,iStep) deltaTime[jStep]
                    attention: we need to take diffs in order to retrieve integrand */
                vector<double>::iterator currentIrIntHelper = irIrIntHelper.begin();
                currentIrIntHelper = baseRatesUtil->factorVariance(
                                        currentIrIntHelper,
                                        aFactor,
                                        baseIrSpotVol[jStep],
                                        deltaTime[jStep]);
                currentIrIntHelper = baseRatesUtil->factorCovariance(
                                        currentIrIntHelper,
                                        aFactor,
                                        baseIrSpotVol[jStep],
                                        deltaTime[jStep]);

                /** second: loop over all equities */
                for (iAsset = 0; iAsset < nbEqAssets; iAsset++) {
                    ASSERT(assetSpotVols[iAsset].size() > jStep);
                    vector<double>::iterator currentirAssetIntHelper =
                        irAssetIntHelper[iAsset].begin();
                    currentirAssetIntHelper =
                        baseRatesUtil->factorFXCovariance(
                            true, // not sure what to choose here ... :S
                            currentirAssetIntHelper,
                            aFactor,
                            irAssetCorrs[iAsset],
                            baseIrSpotVol[jStep]*assetSpotVols[iAsset][jStep],
                            deltaTime[jStep]);
                }
            }

            for (ll = 0; ll < (unsigned int)numIrFactors; ll++) {
                irIrIntegral[iStep] += irIrIntHelper[ll];
                for (iAsset = 0; iAsset < nbEqAssets; iAsset++) {
                    irAssetIntegral[iAsset][iStep] =
                        irAssetIntHelper[iAsset][ll];

                }
            }

            // a few helper vector ...
            for (iAsset = 0; iAsset < nbEqAssets; iAsset++) {
                double totalAssetSpotVar = 0.0;
                for (jStep=0; jStep <= iStep; jStep++) {
                    totalAssetSpotVar +=
                        Maths::square(assetSpotVols[iAsset][jStep]) * deltaTime[jStep];
                }
                totalAssetCompVars[iAsset][iStep] = totalAssetSpotVar + irIrIntegral[iStep]
                    + irAssetIntegral[iAsset][iStep];
                if (iStep == 0) {
                    fwdAssetCompVars[iAsset][iStep] = totalAssetCompVars[iAsset][iStep];
                } else {
                    fwdAssetCompVars[iAsset][iStep] =
                        totalAssetCompVars[iAsset][iStep]
                        - totalAssetCompVars[iAsset][iStep-1];
                }

            }
        } // end of loop trough benchmark dates
        thisSrmEqVolData = SrmEqVolDataSP(new SrmEqVolData(assetSpotVols,
                                                       fwdAssetCompVars,
                                                       totalAssetCompVars,
                                                       irAssetIntegral,
                                                       irIrIntegral));
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** use raw const EqEq or time-dep fwd EqEq corrs and do corr mapping */
void MCPathConfigSRMGen::computeAdjEqEqCorrs(DoubleMatrixArraySP&   fwdEqEqCorrelations,
                                             const IntArray&        fwdCorrIndexArray,
                                             const DateTimeArray&   fwdCorrDatesArray) const {
    static const string method = "MCPathConfigSRMGen::getAdjEqEqCorrs";
    try {
        /** SrmEqVolData holds the following info (for all simstep = simDates.size-1)
                assetSpotVols       sigma(t)
                fwdAssetCompVars    int_0^t compSigma^2(u) du
                irAssetIntegral     2*rho*int_0^t sigma(u) sigma_P(u,t) du
                irIrIntegral        int_0^t sigma_P^2(u,t) du */

        /** dimension checking */
        if (fwdEqEqCorrelations->size() != fwdCorrIndexArray.size()) {
            throw ModelException(method, "Size mitmatch correlation matrix array and fwd index array");
        }
        if (fwdEqEqCorrelations->size() != fwdCorrDatesArray.size()) {
            throw ModelException(method, "Size mitmatch correlation matrix array and fwd dates array");
        }


        int iStep, nbSteps = simDates.size()-1;

        int index;
        double adjCorr;
        for(int iAsset = 0; iAsset < sv->numEq; iAsset++) {
            for (int jAsset = iAsset + 1; jAsset < sv->numEq; jAsset++) {
                index = -1;
                for (int jFwdDate=0; jFwdDate< fwdCorrDatesArray.size(); jFwdDate++) {
                    double deltaT = 0.0;
                    double compVolTerm = 0.0, mixedVolTerm = 0.0, spotVolTerm = 0.0;                    
                    double eqeqCorr = (*(*fwdEqEqCorrelations)[jFwdDate])[iAsset][jAsset];

                    int thisIdx = fwdCorrIndexArray[jFwdDate];
                    if (jFwdDate>0) {
                        int prevIdx = fwdCorrIndexArray[jFwdDate-1];
                        mixedVolTerm = thisSrmEqVolData->irIrIntegral[thisIdx]
                                - thisSrmEqVolData->irIrIntegral[prevIdx]
                            + thisSrmEqVolData->irAssetIntegral[iAsset][thisIdx]/2.0
                                - thisSrmEqVolData->irAssetIntegral[iAsset][prevIdx]/2.0
                            + thisSrmEqVolData->irAssetIntegral[jAsset][thisIdx]/2.0
                                - thisSrmEqVolData->irAssetIntegral[jAsset][prevIdx]/2.0;
                    } else {
                        mixedVolTerm = thisSrmEqVolData->irIrIntegral[thisIdx]
                            + thisSrmEqVolData->irAssetIntegral[iAsset][thisIdx]/2.0
                            + thisSrmEqVolData->irAssetIntegral[jAsset][thisIdx]/2.0;
                    }
                    for (iStep=index+1; iStep<=fwdCorrIndexArray[jFwdDate]; iStep++) {
                        deltaT += simDates[iStep].yearFrac(simDates[iStep+1]);

                        compVolTerm += eqeqCorr *
                            sqrt(thisSrmEqVolData->fwdAssetCompVars[iAsset][iStep]
                                *thisSrmEqVolData->fwdAssetCompVars[jAsset][iStep]);
                        spotVolTerm += thisSrmEqVolData->assetSpotVols[iAsset][iStep] *
                            thisSrmEqVolData->assetSpotVols[jAsset][iStep] * 
                            simDates[iStep].yearFrac(simDates[iStep+1]);
                    } // iStep
                    if ( !Maths::isPositive(deltaT) ) { 
                        adjCorr = eqeqCorr;
                    } else { 
                        adjCorr = (compVolTerm - mixedVolTerm) / spotVolTerm;
                    } 
                    (*(*fwdEqEqCorrelations)[jFwdDate])[iAsset][jAsset] = adjCorr;
                    (*(*fwdEqEqCorrelations)[jFwdDate])[jAsset][iAsset] = adjCorr;
                    index = fwdCorrIndexArray[jFwdDate];
                } // fwdDate
            } // jAsset
        } // iAsset
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
