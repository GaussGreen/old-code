//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRMGenSV.cpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCGenDiffusibleAsset.hpp"

#include "edginc/MemoryProfiler.hpp"
#include "edginc/SparseMatrix.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/CIDParameters.hpp"

USING_DRLIB_NAMESPACE


static /** Extract the dates from a SVGenSpot array */
vector<const DateTimeArray*> getDates(
            int                          assetIdx,
            const vector<const SVGenSpot*>& mcSpot)
{
    vector<const DateTimeArray*> dates(mcSpot.size());
    for (size_t i = 0; i < dates.size(); i++){
        if (assetIdx>=0) {
            dates[i] = &mcSpot[i]->getSimSeries()->getDates(assetIdx);
        } else {
            // !! SNN really want no dates but for the moment hacking round
            // the effort of making the rest of the indexing code safe
            dates[i] = &mcSpot[i]->getSimSeries()->getDates(0);
        }
    }
    return dates;
}


void IQMCGenDiffusibleAsset::createUtil(
    MCPathConfig*          mcPathConfig,
    const DateTime&        today,
    DateTimeArrayConstSP   simDates)
{
    mcPathConfig->createUtil(this, today, simDates );
}

void IQMCGenDiffusibleAsset::setCurrencyTreatment(CurrencyTreatment _ccyTreatment) {
    QLIB_VERIFY(_ccyTreatment == ccyVanilla, "Only vanilla currency treatment is supported by " + name);
    ccyTreatment = _ccyTreatment;
}


DateTimeArray IQMCGenDiffusibleAsset::getCriticalDates(
    MCPathConfig* mcPathConfig,
    const DateTime& start,       // likely to be "today"
    const DateTime& finish) // the latest requested date
{

    return mcPathConfig->getCriticalDates(this,start,finish);
}


void IQMCGenDiffusibleAsset::setDiffusibleAsset(
    const DateTime&     today,              // base date of the run
    MCPathConfig* mcPathConfig,       // to get max/min boundaries, etc
    const IPastValues*  pastValues,         // historic values
    DependenceSP        dependence)    // factorized correlations
{
    mcPathConfig->setDiffusibleAsset(today,this,pastValues,dependence);
}


/************************************************************************/
/*  Interest Rate -related block                                        */
/************************************************************************/

void QMCGenDiffusibleIR::createStateVars(StateVarDBase& svDBase,
                                         vector<vector<SVPathSP> > & /*assetPaths*/,
                                         bool mm)
{
    static const string method("QMCGenDiffusibleIR::createStateVars");
    vector<const SVGenDiscFactor*>&         mcDiscFactors = dfs.first;
    vector<const SVGenExpectedDiscFactor*>& mcExpDiscFactors = dfs.second;
    // then for each SV
    size_t j;
    for (j = 0; j < mcDiscFactors.size(); j++){
        SVDiscFactorSP dfSVSP (getDiffusibleIR()->createDFSV(mcDiscFactors[j]->getDates(), mm));
        svDBase.append(mcDiscFactors[j], dfSVSP);
    }
    for (j = 0; j < mcExpDiscFactors.size(); j++){
        SVExpectedDiscFactorSP edfSVSP(
            getDiffusibleIR()->createExpDiscFactorSV(mcExpDiscFactors[j]->getPVDate(),
                                            mcExpDiscFactors[j]->getDates(),
                                            mcExpDiscFactors[j]->logRequired(),
                                            mcExpDiscFactors[j]->getYieldCurve(), 
                                            mm));
        svDBase.append(mcExpDiscFactors[j], edfSVSP);
    }

}

/************************************************************************/
/*  FX-related block                                                    */
/************************************************************************/

void QMCGenDiffusibleFX::createStateVars(StateVarDBase& svDBase,
                                         vector<vector<SVPathSP> > &assetPaths,
                                         bool mm)
{
    static const string method("QMCGenDiffusibleFX::createStateVars");

    const vector<const SVGenSpot*>&          mcSpot = *multiSpots;
    const vector<const SVGenExpectedSpot*>&  mcExpSpot = fxExpSpots;

    vector<const DateTimeArray*> dates = getDates(assetIdx, mcSpot);

    size_t j;

    for (j = 0; j < mcSpot.size(); j++){

        const DateTimeArray& svDates = *(dates[j]);
        SVSpotFXSP spotSV(getDiffusibleFX()->createSpotFXSV(svDates, mm));
        svDBase.append(mcSpot[j], spotSV);

        if (assetIdx>=0) {
            assetPaths[j][assetIdx] = spotSV;
        }
    }
    for (j = 0; j < mcExpSpot.size(); j++){
        SVExpectedFXSP expSpotSV(
            getDiffusibleFX()->createExpectedFXSV(
                mcExpSpot[j]->getCalcDate(),
                mcExpSpot[j]->getDates(),
                mcExpSpot[j]->logRequired(), mm));
        svDBase.append(mcExpSpot[j], expSpotSV);
    }

}

/************************************************************************/
/*  Equity-related block                                                */
/************************************************************************/
DateTimeArray QMCGenDiffusibleEquity::criticalDates(const DateTime& start,
        const DateTime& end)
{
    return SRMEquityDiffuse::getCriticalDates(eqAsset, start,end);
}


void QMCGenDiffusibleEquity::createStateVars(StateVarDBase& svDBase,
                                              vector<vector<SVPathSP> > &assetPaths,
                                             bool mm)
{
    static const string method("QMCGenDiffusibleEquity::createStateVars");
    
    const vector<const SVGenSpot*>&          mcSpot = *multiSpots;
    const vector<const SVGenExpectedSpot*>&  mcExpSpot = expSpots;

    vector<const DateTimeArray*> dates = getDates(assetIdx, mcSpot);

    size_t j;
    for (j = 0; j < mcSpot.size(); j++){
        const DateTimeArray & measurementDates = * (dates[j]);
        SVSpotEQSP spotSV(getDiffusibleEQ()->createSpotEQSV(measurementDates, mm));
        svDBase.append(mcSpot[j], spotSV);
        if (assetIdx>=0) {
            assetPaths[j][assetIdx] = spotSV;
        }
    }

    for (j = 0; j < mcExpSpot.size(); j++){
        SVExpectedEQSP expSpotSV(
            getDiffusibleEQ()->createExpectedEQSV(
                mcExpSpot[j]->getCalcDate(),
                mcExpSpot[j]->getDates(),
                mcExpSpot[j]->logRequired(), mm));
        svDBase.append(mcExpSpot[j], expSpotSV);
    }

}



/************************************************************************/
/*  Credit-related block                                                    */
/************************************************************************/
//TODO: move to visitor pattern (i.e. visit all of th easset's Gen and create SV

void QMCGenDiffusibleCredit::createStateVars(StateVarDBase& svDBase,
                                              vector<vector<SVPathSP> > & /*assetPaths*/,
                                             bool mm)
{
    static const string method("QMCGenDiffusibleCredit::createStateVars");
    
    vector<const SVGenSurvivalDiscFactor*>&         mcSurvivalDiscFactors = sdf.first;
    vector<const SVGenExpectedSurvivalDiscFactor*>& mcExpSurvivalDiscFactors = sdf.second;
    vector<const SVGenAggregatedSurvivalDiscFactor*>&  mcAggregatedSurvivalDiscFactors = sdf.third;

    vector<const SVGenDateOfDefault*>& mcDateOfDefaults = sdf.dod;

    // then for each SV
    for (size_t j = 0; j < mcSurvivalDiscFactors.size(); j++){
        SVSurvivalDiscFactorSP sdfSP(
            getDiffusibleCredit()->createSDFSV(mcSurvivalDiscFactors[j]->getDates(), mm));
        svDBase.append(mcSurvivalDiscFactors[j], sdfSP);
    }
 
    for (size_t j = 0; j < mcExpSurvivalDiscFactors.size(); j++){
        SVExpSurvDiscFactorSP sdfSP(
            getDiffusibleCredit()->createExpSurvDiscFactorSV(
                mcExpSurvivalDiscFactors[j]->getCalcDate(),
                mcExpSurvivalDiscFactors[j]->getDates(),
                mcExpSurvivalDiscFactors[j]->logRequired(), mm));
        svDBase.append(mcExpSurvivalDiscFactors[j], sdfSP);
    }

    for (size_t j = 0; j < mcAggregatedSurvivalDiscFactors.size(); j++){
        SVAggregatedSurvDiscFactorSP asdfSP(
                getDiffusibleCredit()->createAggregatedSurvDiscFactorSV(
                mcAggregatedSurvivalDiscFactors[j]->sdfDates,
                mcAggregatedSurvivalDiscFactors[j]->sdfIdxSP,
                mcAggregatedSurvivalDiscFactors[j]->esdfRequestedDates,
                mcAggregatedSurvivalDiscFactors[j]->esdfReqIdxSP,
                mcAggregatedSurvivalDiscFactors[j]->esdfForwardDates,
                mcAggregatedSurvivalDiscFactors[j]->esdfForIdxSP,
                mcAggregatedSurvivalDiscFactors[j]->getMaxMatDate(),
                mcAggregatedSurvivalDiscFactors[j]->getMaxCurveDate(),
                mcAggregatedSurvivalDiscFactors[j]->logRequired(), 
                mm));
        svDBase.append(mcAggregatedSurvivalDiscFactors[j], asdfSP);
    }



    if (!mcDateOfDefaults.empty())
    {
        IQMCDiffusibleDefaultableCreditSpread* defaultableAsset = 
            dynamic_cast<IQMCDiffusibleDefaultableCreditSpread*>(getDiffusibleCredit().get());
        QLIB_VERIFY(defaultableAsset != NULL, 
                    "The request for Date Of Default is not implemented in the selected credit diffusion model." );

        SVDateOfDefaultSP dodSP(defaultableAsset->createSVDateOfDefault()); // special! need only one

        for (size_t j = 0; j < mcDateOfDefaults.size(); ++j)
        {
            svDBase.append(mcDateOfDefaults[j], dodSP);
        }
    }
}


/************************************************************************/
/*  Energy block                                                        */
/************************************************************************/

void QMCGenDiffusibleEnergy::createStateVars(
	StateVarDBase& svDBase,
	vector<vector<SVPathSP> > &assetPaths,
    bool mm)
{
	static const string method("SV::QMCGenDiffusibleEnergy::createStateVars");

   	// for each SV gen:
	for (size_t j = 0; j < expFpGens.size(); ++j){
        SVExpEnergyFutureSP expFpSVSP = SVExpEnergyFutureSP(
            getDiffusibleEnergy()->createExpEnergyFutureSV(
				expFpGens[j]->getCalcDate(),
				expFpGens[j]->getDates(),
				expFpGens[j]->logRequired(),
                mm)
			);
		svDBase.append(expFpGens[j], expFpSVSP);
	}
}


void QMCGenDiffusibleEnergy::setTierLevel()
{
    if (pParentQMCGenDiffusibleEnergy) { // check if it has parent
        pParentQMCGenDiffusibleEnergy->setTierLevel();
        tierLevel = pParentQMCGenDiffusibleEnergy->tierLevel + 1;
    }
}

bool DiffEnrgGenSPComparator::operator() (const QMCGenDiffusibleEnergySP& a, const QMCGenDiffusibleEnergySP& b)
{
    return a->tierLevel < b->tierLevel;
}

/************************************************************************/
/*  Basis Spread -related block                                        */
/************************************************************************/

void QMCGenDiffusibleBasisSpread::createStateVars(
    StateVarDBase& svDBase,
    vector<vector<SVPathSP> > & /*assetPaths*/,
    bool mm)
{
    static const string method("QMCGenDiffusibleBasisSpread::createStateVars");

    // then for each SV
    size_t j;
    for (j = 0; j < expSPGens.size(); j++){
        SVExpectedBasisFwdSpreadSP fsSVSP (
            getDiffusibleBasis()->createExpBasisFwdSpreadSV(
                expSPGens[j]->getPVDate(),
                expSPGens[j]->getDates()));
        svDBase.append(expSPGens[j], fsSVSP);
    }
}


void QMCGenPureJump::initializeDiffusibleAssets()
{
    pPureJumpsAsset = QMCPureJumpsSP(new QMCPureJumps());
}

void QMCGenPureJump::setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence)         // factorized correlations
{

    const ICIDParameters* pCIDParameters = mcPathConfig->getCIDParameters();
    QLIB_VERIFY(pCIDParameters != NULL, 
        "Cannot build non-CID model with pure jumps for credit at the moment. Call again later.");

    vector<double> intensities(pCIDParameters->getNumOfJumpFactors());
    for (size_t i= 0; i<intensities.size(); ++i) 
        intensities[i] = pCIDParameters->getJumpFactor(i)->jFreqs;
    pPureJumpsAsset->setQMCPureJumps(today, intensities);
}

/** this one is a generator for jump effects of a jump source on a name */
void QMCGenDiffusibleCreditCIDJumps::createSRMUtil(
    MCPathConfigSRM*       mcPathConfig,
    const DateTime&        today,
    DateTimeArrayConstSP   simDates)
{
    // no util for this type
}

void QMCGenDiffusibleCreditCIDJumps::initializeDiffusibleAssets()
{
    QLIB_VERIFY(pJumpSource.get()!=NULL, "jump source must be initialized before CIDJumps object");
    pAssetCreditCIDJumps = QMCCreditCIDJumpsSP(new QMCCreditCIDJumps(pJumpSource->getPureJumpAsset(), isFullMC));
}

void QMCGenDiffusibleCreditCIDJumps::setSRMDiffusibleAsset(
        const DateTime&        today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,       // to get max/min boundaries, etc
        const IPastValues*     pastValues,         // historic values
        DependenceSP           dependence)         // factorized correlations
{
    const ICIDParameters* pCIDParameters = mcPathConfig->getCIDParameters();

    QLIB_VERIFY(pCIDParameters != NULL, 
        "Cannot initialize CID jumps for credit -- CIDParameters are not defined in MCPathConfig.");

    size_t nJFactors = pCIDParameters->getNumOfJumpFactors();
    vector<double> jumpImpact(nJFactors);   
    vector<double> jumpMeanSize(nJFactors);
    vector<double> jumpMeanReversionSpeed(nJFactors);

    CIDNameRecordConstSP pNameRecord = pCIDParameters->getSingleNameRecord(name);
    for(size_t i=0; i<nJFactors; ++i)
    {
        CIDJumpFactorConstSP pJumpSource = pCIDParameters->getJumpFactor(i);
        jumpImpact[i] = (*pNameRecord->jImpacts)[i];
        jumpMeanSize[i] = (*pNameRecord->jSizes)[i];
        jumpMeanReversionSpeed[i] = pJumpSource->jDecays;
    }

    pAssetCreditCIDJumps->setQMCCreditCIDJumps(
        today, 
        jumpImpact,
        jumpMeanSize,
        jumpMeanReversionSpeed);

}


/** this one is a generator of aggregating asset, which will be able to return 
    everything for a single name */



void QMCGenDiffusibleCreditComposite::createSRMUtil(
    MCPathConfigSRM*       mcPathConfig,
    const DateTime&        today,
    DateTimeArrayConstSP   simDates)
{
    // no particular util is necessary
}

void QMCGenDiffusibleCreditComposite::initializeDiffusibleAssets()
{
    // the components and weights inside this object should be already 
    // populated before this function is called

    QLIB_VERIFY(!components.empty(), "cannot create composite credit as the list of components is empty");
    QLIB_VERIFY(components.size()==weights.size(), "the list of components differs in size from the list of weights");

    vector<IQMCDiffusibleCreditSpreadBaseSP> componentAssets;
    for(size_t i=0; i<components.size(); ++i)
        componentAssets.push_back(components[i]->getDiffusibleCredit());

    pAssetCreditComposite = QMCCreditCompositeSP(new QMCCreditComposite(componentAssets,weights));
}

void QMCGenDiffusibleCreditComposite::setSRMDiffusibleAsset(
        const DateTime&        today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,       // to get max/min boundaries, etc
        const IPastValues*     pastValues,         // historic values
        DependenceSP           dependence)         // factorized correlations
{
    pAssetCreditComposite->setQMCCreditComposite(today, regularRecoveryRate );
}

void QMCGenDiffusibleCreditComposite::addComponent(QMCGenDiffusibleCreditSP component, double weight)
{
    components.push_back(component);
    weights.push_back(weight);
}


void QMCGenDetermRates::initializeDiffusibleAssets()
{
    pAssetRates = QMCDetermRatesSP(new QMCDetermRates());
}

void QMCGenDetermRates::setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence)         // factorized correlations
{
    pAssetRates->setQMCDetermRates(today, ycs.first);
}
 


