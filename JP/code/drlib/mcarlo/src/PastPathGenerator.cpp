//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PastPathGenerator.cpp
//
//   Description : Monte Carlo path generator for historic 'pass'
//
//   Date        : Nov 6 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/SVGenDemo.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenExpectedSpot.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVGenPathWeight.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductClient.hpp"

DRLIB_BEGIN_NAMESPACE

/** Class for dealing with past - hopefully we can reuse this for different
    generators. Supports different dates per asset */

void PastPathGenerator::generatePath(int pathIdx){ /* nothing to do */}

int PastPathGenerator::getPathIndex() const {
    return 0; // only ever one of these
}

/** Returns true if this path generator is being used to 'simulate'
        the past */
bool PastPathGenerator::doingPast() const{
    return true;
}
    
/** to do: check that this is the number of simulated assets
        rather than the number of assets that the payoff is expecting */
int  PastPathGenerator::NbSimAssets() const{
    return numAssets;
}

/** Returns 1.0 */
double PastPathGenerator::maxDriftProduct(int iAsset) const{
    return 1.0;
}

const double* PastPathGenerator::Path(int iAsset, int iPath) const{ 
    return paths[iAsset].empty()? 
        0: paths[iAsset][0]; // independent of iPath
}

/** Returns the reference level for given asset and path */
double PastPathGenerator::refLevel(int iAsset, int iPath) const{
    return refLevels[iAsset]; // independent of iPath
}

/** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
int PastPathGenerator::begin(int iAsset) const{
    return 0; // same for everything
}

//// see  begin(int iAsset)
int PastPathGenerator::end(int iAsset) const{
    return paths[iAsset].numRows();
}

/** Returns true if there are historic simulation dates */
bool PastPathGenerator::hasPast() const{
    return simHasPast;
}

//// Constructor
PastPathGenerator::PastPathGenerator(const IMCProduct*  prod):
numAssets(prod->getNumAssets()),
simHasPast(false),
refLevels(numAssets), 
paths(numAssets){
    static const string routine = "PastPathGenerator::PastPathGenerator";
    try{
        // get hold of things we need from IMCProduct
        const IMultiFactors* mAsset = prod->getMultiFactors();
        const SimSeries*   simSeries = prod->getSimSeries();
        const IPastValues* pastValues = prod->getMCPastValues();
        const DateTime&  today = prod->getToday();
        const IRefLevel*  refLevelObj = prod->getRefLevel();
        // before we create the IRefLevel::IMCPath need to decide what
        // values we're going to give for spots at sim start date.
        // Follow 'convention' and set to spot
        DoubleArray spotsAtSimStart(numAssets);
        for (int i = 0; i < numAssets; i++){
            spotsAtSimStart[i] = mAsset->assetGetSpot(i);
        }
        IRefLevel::IMCPathSP refLevelPath(
            refLevelObj->createMCPath(today, spotsAtSimStart, pastValues));
        if(!simHasPast) {
            QLIB_VERIFY(! refLevelObj->getAllDates().empty(), "Calling back() on empty vector");
            if(refLevelObj->getAllDates().back() <= today) {
                simHasPast = true;
            }
        }

        // then store past in our format
        for (int iAsset = 0; iAsset < numAssets; iAsset++){
            DoubleArray past(simSeries->getPastValues(today, iAsset,
                                                      pastValues));
            if (!past.empty()){
                // slightly painful mapping to double matrix
                paths[iAsset] = DoubleMatrix(past);
                if (!simHasPast){
                    simHasPast = true;
                }
            } 
            // sort out refLevels (needed for volInterp)
            refLevels[iAsset] = refLevelPath->refLevel(iAsset, mAsset);
        }
    } 
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


IStateVariableSP PastPathGenerator::create(const IStateVariableGen* svGen) {
    throw ModelException(
        "PastPathGenerator::create", 
        "Path generator does not support state variables");
}


#ifdef STATE_VARIABLES

void PastPathGen::generatePath(int pathIdx){ 
    pathGenSpot->generatePath(pathIdx);
}

void PastPathGen::advance()
{
    pathGenSpot->advance();
}

int PastPathGen::getPathIndex() const {
    return 0; // only ever one of these
}

/** Returns true if this path generator is being used to 'simulate'
        the past */
bool PastPathGen::doingPast() const{
    return true;
}
    
/** to do: check that this is the number of simulated assets
        rather than the number of assets that the payoff is expecting */
int  PastPathGen::NbSimAssets() const{
    static const string routine = "PastPathGen::NbSimAssets";
    throw ModelException(routine, "Method is retired for StateVars");
}

/** Returns 1.0 */
double PastPathGen::maxDriftProduct(int iAsset) const{
    static const string routine = "PastPathGen::maxDriftProduct";
    throw ModelException(routine, "Method is retired for StateVars");
}

const double* PastPathGen::Path(int iAsset, int iPath) const{ 
    static const string routine = "PastPathGen::Path";
    throw ModelException(routine, "Method is retired for StateVars");
}

/** Returns the reference level for given asset and path */
double PastPathGen::refLevel(int iAsset, int iPath) const{
    static const string routine = "PastPathGen::refLevel";
    throw ModelException(routine, "Method is retired for StateVars");
}

/** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
int PastPathGen::begin(int iAsset) const{
    static const string routine = "PastPathGen::begin";
    throw ModelException(routine, "Method is retired for StateVars");
}

//// see  begin(int iAsset)
int PastPathGen::end(int iAsset) const{
    static const string routine = "PastPathGen::end";
    throw ModelException(routine, "Method is retired for StateVars");
}

/** Returns true if there are historic simulation dates */
bool PastPathGen::hasPast() const{
    return (pathGenSpot->hasPast() ||
            pathGenQuadVar->hasPast() ||
            pathGenSqrtAnnualQuadVar->hasPast());
}

/** Constructor */
PastPathGen::PastPathGen(const MCProductClient* prodClient):
    today(prodClient->getToday()) {
    static string const routine = "PastPathGen::PastPathGen";
    try {
        // Collect state variables from product and categorize them
        StateVariableCollectorSP svCollector(new StateVariableCollector());
        prodClient->collectStateVars(svCollector);
        IElemStateVariableGenArray stateVarGenArray(
            svCollector->getElemStateVarGens());

        // 1) Spot requests
        SVGenSpotArray spotGenArray = filterStateVars<SVGenSpot>(stateVarGenArray);

        // 2) MCQuadVar requests
        MCQuadVarArray quadVarGenArray = filterStateVars<MCQuadVar>(stateVarGenArray);

        // 3) MCSqrtAnnualQuadVar requests
        MCSqrtAnnualQuadVarArray sqrtAnnualQuadVarGenArray = 
            filterStateVars<MCSqrtAnnualQuadVar>(stateVarGenArray);

        // 4) Create discount factor state variables
        SVGenDiscFactorArray discFactors(
            filterStateVars<SVGenDiscFactor>(stateVarGenArray));
        unsigned int iVar;
        for(iVar = 0; iVar < discFactors.size(); iVar++) {
            IStateVariableSP sv(discFactors[iVar]->determinsticSV(
                                 true /* doing past */));
            svDBase.append(discFactors[iVar], sv);
        }

        // 5) Create expected discount factor (ZCBs) state variables
        vector<const SVGenExpectedDiscFactor*>  expDiscFactors(
            filterStateVars<SVGenExpectedDiscFactor>(stateVarGenArray));
        for(iVar = 0; iVar < expDiscFactors.size(); iVar++) {
            IStateVariableSP sv(expDiscFactors[iVar]->determinsticSV(
                                 today, true /* doing past */));
            svDBase.append(expDiscFactors[iVar], sv);
        }

        // 6) Create survival discount factor state variables
        SVGenSurvivalDiscFactorArray survivalDiscFactors(
            filterStateVars<SVGenSurvivalDiscFactor>(stateVarGenArray));
        for(iVar = 0; iVar < survivalDiscFactors.size(); iVar++) {
            IStateVariableSP sv(survivalDiscFactors[iVar]->determinsticSV(
                                 true /* doing past */));
            svDBase.append(survivalDiscFactors[iVar], sv);
        }

        // 7) Create expected survival discount factor state variables
        vector<const SVGenExpectedSurvivalDiscFactor*>  expSurvivalDiscFactors(
            filterStateVars<SVGenExpectedSurvivalDiscFactor>(stateVarGenArray));
        for(iVar = 0; iVar < expSurvivalDiscFactors.size(); iVar++) {
            IStateVariableSP sv(expSurvivalDiscFactors[iVar]->determinsticSV(
                                 today, true /* doing past */));
            svDBase.append(expSurvivalDiscFactors[iVar], sv);
        }

        // 8) Create expected spot variables
        vector<const SVGenExpectedSpot*>  expSpot(
            filterStateVars<SVGenExpectedSpot>(stateVarGenArray));
        for(iVar = 0; iVar < expSpot.size(); iVar++) {
            IStateVariableSP sv(expSpot[iVar]->createHistoricZeroSV(today));
            svDBase.append(expSpot[iVar], sv);
        }

        // 9) Create expected energy future SVs
        vector<const SVGenExpectedEnergyFuture*> expEnrgFutureGens(
            filterStateVars<SVGenExpectedEnergyFuture>(stateVarGenArray));
        for (iVar = 0; iVar < expEnrgFutureGens.size(); ++iVar) {
            IStateVariableSP sv(expEnrgFutureGens[iVar]->getPastSV(
                                today, true /* doing past */));
            svDBase.append(expEnrgFutureGens[iVar], sv);
        }

        //10) Create expected basis spread SVs
        vector<const SVGenExpectedBasisFwdSpread*> expBasisSprdGens(
            filterStateVars<SVGenExpectedBasisFwdSpread>(stateVarGenArray));
        for (iVar = 0; iVar < expBasisSprdGens.size(); ++iVar) {
            IStateVariableSP sv(expBasisSprdGens[iVar]->determinsticSV(
                                today, true /* doing past */));
            svDBase.append(expBasisSprdGens[iVar], sv);
        }
        //10) Create credit date of default SVs
        vector<const SVGenDateOfDefault*> dateOfDefaultGens(
            filterStateVars<SVGenDateOfDefault>(stateVarGenArray));
        for (iVar = 0; iVar < dateOfDefaultGens.size(); ++iVar) {
            IStateVariableSP sv(dateOfDefaultGens[iVar]->determinsticSV(
                                true /* doing past */));
            svDBase.append(dateOfDefaultGens[iVar], sv);
        }

        // Create DemoStateVariables
        vector<const SVGenDemo*> jSv(filterStateVars<SVGenDemo>(stateVarGenArray));
        for (iVar = 0; iVar < jSv.size(); iVar++) {
            IStateVariableSP sv(jSv[iVar]->pastNewSV());
            svDBase.append(jSv[iVar], sv);
        }

        // Create PathWeightSV
        vector<const SVGenPathWeight*> pwSv(filterStateVars<SVGenPathWeight>(stateVarGenArray));
        for (iVar = 0; iVar < pwSv.size(); iVar++) {
            IStateVariableSP sv(pwSv[iVar]->determinsticSV(true /* doing past */));
            svDBase.append(pwSv[iVar], sv);
        }

        // Create DemoStateVariables
        for (iVar = 0; iVar < jSv.size(); iVar++) {
            IStateVariableSP sv(jSv[iVar]->pastSV());
            // svDBase.append(jSv[iVar], sv);
        }

        // Ignore remaining state variables as they might refer to future
        // e.g. HitValueBB

        // Create simulation modules
        pathGenSpot = PastPathGenSpotSP(new 
            PastPathGenSpot(prodClient, spotGenArray, svDBase));
        pathGenQuadVar = PastPathGenQuadVarSP(new 
            PastPathGenQuadVar(prodClient, quadVarGenArray,svDBase));
        pathGenSqrtAnnualQuadVar = PastPathGenSqrtAnnualQuadVarSP(new 
            PastPathGenSqrtAnnualQuadVar(prodClient, sqrtAnnualQuadVarGenArray,svDBase));
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** Returns past path generator for spot */
PastPathGenSpotSP PastPathGen::getPathGenSpot() const {
    return pathGenSpot;
}

/** Returns past path generator for QuadVar */
PastPathGenQuadVarSP PastPathGen::getPathGenQuadVar() const {
    return pathGenQuadVar;
}

/** Returns past path generator for SqrtAnnualQuadVar */
PastPathGenSqrtAnnualQuadVarSP PastPathGen::getPathGenSqrtAnnualQuadVar() const {
    return pathGenSqrtAnnualQuadVar;
}

IStateVariableSP PastPathGen::create(const IStateVariableGen* svGen) {
    static const string routine = "PastPathGen::create";

    try {
        return svDBase.find(svGen);
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

///////////////////////////////////////////////////////////////

//// Constructor
PastPathGenSpot::PastPathGenSpot(
        const MCProductClient* prodClient, 
        const SVGenSpotArray&     spotGenArray,
        StateVarDBase&         svDBase):
numAssets(prodClient->getNumAssets()),
simHasPast(false),
paths(numAssets) {    
    static const string routine = "PastPathGenSpot::PastPathGenSpot";
    try{
        if(spotGenArray.size()) {
            // Loop over assets and
            //  i)  aggregate asset timeline from different SVGenSpotGens
            //  ii) populate asset path with past values
            DateTimeArrayArray allDates(numAssets);
            const IPastValues* pastValues = prodClient->getSVGenSpotPastValues();
            const DateTime& today = prodClient->getToday();
            // Obtain path pointers
            vector<const double*> pathPointers(numAssets);
            vector<int> beginIndices(numAssets);
            vector<int> endIndices(numAssets);
            int iAsset;
            for (iAsset = 0; iAsset < numAssets; iAsset++){
                // Get past asset dates and put them in allDates array
                DateTimeArraySP assetDates = MCPath::getAllDates(spotGenArray, iAsset);
                DateTimeArray assetPastDates = today.getPastDates(*assetDates);
                allDates[iAsset] = *assetDates;

                // Get the past values
                DoubleArray assetPastValues = 
                    pastValues->getPastValues(assetPastDates, iAsset, today);

                // Populate asset's past values
                paths[iAsset] = assetPastValues;
                if (!simHasPast && assetPastValues.size() > 0){
                    simHasPast = true;
                }
    
                pathPointers[iAsset] = assetPastDates.size() ? &paths[iAsset][0] : 0;
                beginIndices[iAsset] = 0;
                endIndices[iAsset]   = assetPastDates.size();

            }

            // 5) Create SVGenSpot::IStateVars
            vector<double> maxDrifts(numAssets, 1.0);
            SVGenSpot::IStateVarArray spotStateVarArray(SVGenSpot::createPaths(
                true,
                spotGenArray,
                allDates,
                beginIndices,
                endIndices,
                pathPointers,
                maxDrifts));
    
            // Populate the map
            for(unsigned int iVar = 0; iVar < spotStateVarArray.size(); iVar++) {
                svDBase.append(spotGenArray[iVar], spotStateVarArray[iVar]);
                spotSvSet.push_back(spotStateVarArray[iVar].get());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void PastPathGenSpot::advance()
{
    for (size_t i = 0; i <spotSvSet.size(); ++i)
        spotSvSet[i]->advance();
}

void PastPathGenSpot::generatePath(int pathIdx){ 
    // Nothing to do
}

bool PastPathGenSpot::hasPast() const {
    return simHasPast;
}

bool PastPathGenSpot::doingPast() const {
    return true;
}

//// Constructor
PastPathGenQuadVar::PastPathGenQuadVar(
        const MCProductClient* prodClient, 
        const MCQuadVarArray&  quadVarGenArray,
        StateVarDBase&         svDBase):
numAssets(prodClient->getNumAssets()),
simHasPast(false),
paths(numAssets) {   
    static const string routine = "PastPathGenQuadVar::PastPathGenQuadVar";
    try{
        if(quadVarGenArray.size()) {
            // Loop over assets and
            //  i)  aggregate asset timeline from different MCQuadVarGens
            //  ii) populate asset path with past values
            DateTimeArrayArray allDates(numAssets);
#if 0
            const IPastValues* pastValues = prodClient->getMCPastValues();
#endif
            const DateTime& today = prodClient->getToday();
            // Obtain path pointers
            vector<const double*> pathPointers(numAssets);
            vector<int> beginIndices(numAssets);
            vector<int> endIndices(numAssets);
            int iAsset;
            for (iAsset = 0; iAsset < numAssets; iAsset++){
                // Get past asset dates and put them in allDates array
                DateTimeArraySP assetDates = MCPath::getAllDates(quadVarGenArray, iAsset);
                DateTimeArray assetPastDates = today.getPastDates(*assetDates);
                allDates[iAsset] = *assetDates;

#if 1           //  XXX no past supported for now
                if (assetPastDates.size()){
                    throw ModelException(routine, 
                                         "past quadratic variation values not currently supported");
                }
                DoubleArray assetPastValues;
#else
                // Get the past values
                DoubleArray assetPastValues = 
                    pastValues->getPastValues(assetPastDates, iAsset, today);
#endif

                // Populate asset's past values
                paths[iAsset] = assetPastValues;
                if (!simHasPast && assetPastValues.size() > 0){
                    simHasPast = true;
                }
    
                pathPointers[iAsset] = assetPastDates.size() ? &paths[iAsset][0] : 0;
                beginIndices[iAsset] = 0;
                endIndices[iAsset]   = assetPastDates.size();
            }

            // 5) Create MCQuadVar::IStateVars
            vector<double> maxDrifts(numAssets, 1.0);
            MCQuadVar::IStateVarArray quadVarStateVarArray(MCQuadVar::createPaths(
                true,
                quadVarGenArray,
                allDates,
                beginIndices,
                endIndices,
                pathPointers,
                maxDrifts));

            // Populate the map
            unsigned int iVar;
            for(iVar = 0; iVar < quadVarStateVarArray.size(); iVar++) {
                svDBase.append(quadVarGenArray[iVar], quadVarStateVarArray[iVar]);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


void PastPathGenQuadVar::generatePath(int pathIdx){ 
    // Nothing to do
}

bool PastPathGenQuadVar::hasPast() const {
    return simHasPast;
}

bool PastPathGenQuadVar::doingPast() const {
    return true;
}

//// Constructor
PastPathGenSqrtAnnualQuadVar::PastPathGenSqrtAnnualQuadVar(
        const MCProductClient*          prodClient, 
        const MCSqrtAnnualQuadVarArray& sqrtAnnualQuadVarGenArray,
        StateVarDBase&                  svDBase):
numAssets(prodClient->getNumAssets()),
simHasPast(false),
paths(numAssets) {    
    static const string routine = "PastPathGenSqrtAnnualQuadVar::PastPathGenSqrtAnnualQuadVar";
    try{    
        if(sqrtAnnualQuadVarGenArray.size()) {
            // 4) Loop over assets and
            //    i)  aggregate asset timeline from different MCSqrtAnnualQuadVarGens
            //    ii) populate asset path with past values
            DateTimeArrayArray allDates(numAssets);
            const IPastValues* pastValues = prodClient->getMCSqrtAnnualQuadVarPastValues();
            const DateTime& today = prodClient->getToday();
            // Obtain path pointers
            vector<const double*> pathPointers(numAssets);
            vector<int> beginIndices(numAssets);
            vector<int> endIndices(numAssets);
            int iAsset;
            for (iAsset = 0; iAsset < numAssets; iAsset++){
                // Get past asset dates and put them in allDates array
                DateTimeArraySP assetDates = MCPath::getAllDates(sqrtAnnualQuadVarGenArray, iAsset);
                DateTimeArray assetPastDates = today.getPastDates(*assetDates);
                allDates[iAsset] = *assetDates;

#if 0           //  XXX no past supported for now
                if (assetPastDates.size()){
                    throw ModelException(routine, 
                                         "past quadratic variation values not currently supported");
                }
                DoubleArray assetPastValues;
#else
                // Get the past values
                DoubleArray assetPastValues = 
                    pastValues->getPastValues(assetPastDates, iAsset, today);
#endif

                // Populate asset's past values
                paths[iAsset] = assetPastValues;
                if (!simHasPast && assetPastValues.size() > 0){
                    simHasPast = true;
                }
        
                pathPointers[iAsset] = assetPastDates.size() ? &paths[iAsset][0] : 0;
                beginIndices[iAsset] = 0;
                endIndices[iAsset]   = assetPastDates.size();
            }

            // 5) Create MCSqrtAnnualQuadVar::IStateVars
            vector<double> maxDrifts(numAssets, 1.0);
            MCSqrtAnnualQuadVar::IStateVarArray sqrtAnnualQuadVarStateVarArray(MCSqrtAnnualQuadVar::createPaths(
                true,
                sqrtAnnualQuadVarGenArray,
                allDates,
                beginIndices,
                endIndices,
                pathPointers,
                maxDrifts));
    
            // Populate the map
            unsigned int iVar;
            for(iVar = 0; iVar < sqrtAnnualQuadVarStateVarArray.size(); iVar++) {
                svDBase.append(sqrtAnnualQuadVarGenArray[iVar], sqrtAnnualQuadVarStateVarArray[iVar]);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


void PastPathGenSqrtAnnualQuadVar::generatePath(int pathIdx){ 
    // Nothing to do
}

bool PastPathGenSqrtAnnualQuadVar::hasPast() const {
    return simHasPast;
}

bool PastPathGenSqrtAnnualQuadVar::doingPast() const {
    return true;
}

#endif

DRLIB_END_NAMESPACE
