//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMGenDiffusibleAsset.cpp
//
//   Description : The implementation of generators for creating QMCDiffusibleAsset classes
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMGenDiffusibleAsset.hpp"
#include "edginc/QMCCreditCIDDiffuse.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/SRMEquityDiffuseEuler.hpp"
#include "edginc/SRMEquityDiffuseMappingMethod.hpp"
#include "edginc/SRMEquityDiffuseSecondOrder.hpp"
#include "edginc/SRMEquityDiffuseStrictSecondOrder.hpp"
#include "edginc/SparseMatrix.hpp"
#include "edginc/SRMEnergyUtil.hpp"
#include "edginc/SRMCorrelation.hpp"

USING_DRLIB_NAMESPACE

/************************************************************************/
/*  Rates SRM block                                                     */
/************************************************************************/
vector<double> SRMGenDiffusibleIR::expandAssetToIRCorrs(
    QMCGenDiffusibleIR* irAsset, 
    double corr)
{
    SRMGenDiffusibleIR* irAssetSRM = dynamic_cast<SRMGenDiffusibleIR*>(irAsset);
    QLIB_VERIFY( irAssetSRM != 0, string() + "Expect SRMGenDiffusibleIR but '"
        + typeid(*irAsset).name() + "' was used");

    //return expandAssetToIRCorrsSRM(irAssetSRM, corr);
    // ir vs ir expansion
    return SRMCorrelation::ExpandAssetAssetCorrs(irAssetSRM, this, corr, SRMCorrelation::PC);
}

vector<double> SRMGenDiffusibleIR::expandAssetToENCorrs(
    QMCGenDiffusibleEnergy* enAsset, 
    double corr)
{
    SRMGenDiffusibleEnergy* enAssetSRM = dynamic_cast<SRMGenDiffusibleEnergy*>(enAsset);
    QLIB_VERIFY( enAssetSRM != 0, string() + "Expect SRMGenDiffusibleEnergy but '"
        + typeid(*enAsset).name() + "' was used");

    // Energy vs ir expansion
    return SRMCorrelation::ExpandAssetAssetCorrs(enAssetSRM, this, corr, SRMCorrelation::PC, SRMCorrelation::PC);
}

vector<double> SRMGenDiffusibleIR::expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to ir
}
vector<double> SRMGenDiffusibleIR::expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to ir
}
vector<double> SRMGenDiffusibleIR::expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to ir
}
vector<double> SRMGenDiffusibleIR::expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to ir
}

/*
vector<double> SRMGenDiffusibleIR::expandAssetToIRCorrsSRM(
    SRMGenDiffusibleIR* irAsset, 
    double corr)
{
    return SRMCorrelation::ExpandAssetAssetCorrs(irAsset, this, corr, SRMCorrelation::PC);
}
*/

/************************************************************************/
/*  Rates HJM-related block                                             */
/************************************************************************/

void SRMGenDiffusibleIRHJM::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                    const DateTime&        today,
                                    DateTimeArrayConstSP   simDates)
{

    CVolProcessedSP processedVol(
                                mcPathConfig->getProcessedIRVol(name, 
                                                                ycs.second));
    // now we can build the SRMRatesHJMUtil object for this ccy
    const string& isoCode = name; 
    ratesHJMUtil = SRMRatesHJMUtilSP(
                    new SRMRatesHJMUtil(today,
                                mcPathConfig->getNumIRFactors(isoCode),
                                mcPathConfig->getIRModelParams(isoCode),
                                mcPathConfig->getIRSmileParams(isoCode),
                                processedVol, 
                                ycs.first, //discYC
                                ycs.second, // diffYC
                                mcPathConfig->skipNegIRVols(isoCode),
                                mcPathConfig->getFlatVolIr(),
                                mcPathConfig->getChoiceCutoff(),
                                mcPathConfig->getCutoffValue(),
                                mcPathConfig->getCorrelationSwapStart(isoCode),
                                mcPathConfig->getCorrelationSwapMat(isoCode),
                                mcPathConfig->getCorrelationSwapDCC(isoCode),
                                mcPathConfig->getCorrelationSwapFreq(isoCode)));

    ratesHJMUtil->setTimeLine(simDates);
    ratesHJMUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleIRHJM::initializeDiffusibleAssets()
{
    static const string method("SRMGenDiffusibleIRHJM::initializeDiffusibleIR()");

    switch(numFactors())
    {
    case 1: pAssetHJMRates = SRMRatesHJMDiffuseSP(new SRMRatesHJM1F());
            break;
    case 2: pAssetHJMRates = SRMRatesHJMDiffuseSP(new SRMRatesHJM2F());
            break;
    case 3: pAssetHJMRates = SRMRatesHJMDiffuseSP(new SRMRatesHJM3F());
            break;
    default: throw ModelException(method, 
                "Number of IR factors > 3 not supported");
    }
}

void SRMGenDiffusibleIRHJM::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  /*pastValues*/,     // historic values
                DependenceSP        dependence)         // factorized correlations
{
    
    double irFxCorr[]={0,0,0};

    if (pQMCGenDiffusibleFX) {  
        // there is FX asset too, so
        // look up the extended correlations of IR vs FX
        for(int nnf = 0; nnf<nFactors; ++nnf)
            irFxCorr[nnf] = dependence->getCorrelationIJ(pQMCGenDiffusibleFX->getRandomIdx(), 
                                                         randomIdx+nnf, 0);
    }

    // IR pastValues always return a vector of zeroes. This simple logic has been moved to
    // SRMRatesCreate for the time being until anything meaningful is needed
    /*    vector<double> pastValues(getDiscFactorSVPast(pastPathGenerator,
                                                        mcPathConfig,
                                                        mcDiscFactors,
                                                        today,
                                                        cmptSVSIndexes,
                                                        cmptSVSBegin,
                                                        svsIndexes,
                                                        svsBegin));*/

    // set all the values of pAssetHJMRates
    pAssetHJMRates->setSRMRatesHJMDiffuse(
        randomIdx, 
        today, 
        ratesHJMUtil, 
        true,   // saveDomLnMoney: was "i == 0 && sv->mapIRGen.size() > 1,"
                // Could be better by doing only equity ccys and dom
        true,   // saveSigmaR: for the time being save 
                // it always; but we would need it only for credit 
        mcPathConfig->getNumSDMaxEffRateIR(),
        mcPathConfig->getNumSDMinEffRateIR(),
        vector<double>(),//TODO: pastValues,
        irFxCorr,
		mcPathConfig->calibrateAgainstSwaptionVols(name));
    ratesHJMUtil = SRMRatesHJMUtilSP(); //FIXME: doublecheck
}

void SRMGenDiffusibleIRHJM::pushbackICERecord(StateVarDBase& svdbICE, vector<ICE>& vICE)
{
    vICE.push_back(ICE(ratesHJMUtil, pAssetHJMRates));
    // collect the state variables
    StateVariableCollectorSP svCollector(new StateVariableCollector());
    vICE.back().collectStateVars(svCollector);
    IElemStateVariableGenArray stateVarGenArray(svCollector->getElemStateVarGens());

    // get the discount factor SVs
    vector<const SVGenDiscFactor*> icedfs(
        filterPointers<const IElemStateVariableGen*, 
                       const SVGenDiscFactor*>(stateVarGenArray));

    // and add them to the array
    for(size_t i=0; i<icedfs.size(); ++i)
    {
        SVDiscFactorSP dfSVSP (getDiffusibleIR()->createDFSV(icedfs[i]->getDates(), false));
        svdbICE.append(icedfs[i], dfSVSP);
    }

    // get the expected discount factor SVs
    vector<const SVGenExpectedDiscFactor*> iceexpDfs(
        filterPointers<const IElemStateVariableGen*, 
                       const SVGenExpectedDiscFactor*>(stateVarGenArray));

    for (size_t j = 0; j < iceexpDfs.size(); j++)
    {
        SVExpectedDiscFactorSP edfSVSP(
            getDiffusibleIR()->createExpDiscFactorSV(iceexpDfs[j]->getPVDate(),
                                            iceexpDfs[j]->getDates(),
                                            iceexpDfs[j]->logRequired(),
                                            iceexpDfs[j]->getYieldCurve(), false));
        svdbICE.append(iceexpDfs[j], edfSVSP);
    }
}

void SRMGenDiffusibleIRHJM::expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const
{
    DoubleMatrix assetTMtx(ratesHJMUtil->get3A());
    int randomIndex = getRandomIdx();
    outputSparseCorrMatrix.selfMultiplyToSmallmtx(assetTMtx, randomIndex, randomIndex);
}

vector<double> SRMGenDiffusibleIRHJM::expandTier2Betas(double beta) const
{
    vector<double> tmpCorrs;
    if ( numFactors() > 1 ) {// multifactor Tier1
        // split the one-factor vs ir correlation
        tmpCorrs = SRMCorrelation::assetIr(
                beta,
                *ratesHJMUtil,
                SRMCorrelation::EXPONENTIAL_FACTORS );
    } else {
        tmpCorrs = vector<double>(1,beta);
    }
    return tmpCorrs;
}

DateTimeArray SRMGenDiffusibleIRHJM::getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
                                                   const DateTime& start, 
                                                   const DateTime& finish)
{
    if (mcPathConfig->calibrateAgainstSwaptionVols(name))
    {
        return CriticalDateCollector::collectVolDates(
                ycs.first, IRVolBase::TYPE, start, finish);
        // ICE: will need to store bmDates or similar. Note that
        // there may be a slight discrepancy here between bmDates
        // and dates returned from getBMDetails for FIX style
        // calibration (latter might be curtailed)
        // expDates = DateTime::merge(expDates, bmDates);
    }
    return DateTimeArray();
}

/************************************************************************/
/*  Rates Libor-related block                                             */
/************************************************************************/

void SRMGenDiffusibleIRLibor::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                      const DateTime&        today,
                                      DateTimeArrayConstSP   simDates)
{
    CVolProcessedSP      processedVol(
                                mcPathConfig->getProcessedIRVol(name, 
                                                                ycs.second)); // why second ?
    // now we can build the SRMRatesLiborUtil object for this ccy
    const string& isoCode = name; 
    ratesLiborUtil = SRMRatesLiborUtilSP(
                    new SRMRatesLiborUtil(today,
                                mcPathConfig->getNumIRFactors(isoCode),
                                mcPathConfig->getIRModelParams(isoCode),
                                mcPathConfig->getIRSmileParams(isoCode),
                                processedVol, 
                                ycs.first, //discYC
                                ycs.second, // diffYC
                                mcPathConfig->skipNegIRVols(isoCode),
                                mcPathConfig->getFlatVolIr(),
                                mcPathConfig->getChoiceCutoff(),
                                mcPathConfig->getCutoffValue(),
                                mcPathConfig->getCorrelationSwapStart(isoCode),
                                mcPathConfig->getCorrelationSwapMat(isoCode),
                                mcPathConfig->getCorrelationSwapDCC(isoCode),
                                mcPathConfig->getCorrelationSwapFreq(isoCode)));

    ratesLiborUtil->setTimeLine(simDates);
    ratesLiborUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleIRLibor::expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const
{
    DoubleMatrix assetTMtx(ratesLiborUtil->get3A()); // This appears to be not doing anything at this time.
    int randomIndex = getRandomIdx();
    outputSparseCorrMatrix.selfMultiplyToSmallmtx(assetTMtx, randomIndex, randomIndex);
}

vector<double> SRMGenDiffusibleIRLibor::expandTier2Betas(double beta) const
{
    vector<double> tmpCorrs;
    if ( numFactors() > 1 ) {// multifactor Tier1
        // split the one-factor vs ir correlation
        tmpCorrs = SRMCorrelation::assetIr(
            beta,
            *ratesLiborUtil,
            SRMCorrelation::EXPONENTIAL_FACTORS );
    } else {
        tmpCorrs = vector<double>(1,beta);
    }
    return tmpCorrs;
}

void SRMGenDiffusibleIRLibor::initializeDiffusibleAssets()
{
    pAssetLiborRates = SRMRatesLiborDiffuseSP(new SRMRatesLiborDiffuse());
}

void SRMGenDiffusibleIRLibor::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  /*pastValues*/,     // historic values
                DependenceSP        dependence)         // factorized correlations
{
    
    vector<double> irFxCorr(nFactors);

    if (pQMCGenDiffusibleFX) {  
        // there is FX asset too, so
        // look up the extended correlations of IR vs FX
        for(int nnf = 0; nnf<nFactors; ++nnf)
            irFxCorr[nnf] = dependence->getCorrelationIJ(pQMCGenDiffusibleFX->getRandomIdx(), 
                                                         randomIdx+nnf, 0);
    }

    const string& isoCode = name; 
	const string calibStyle = mcPathConfig->getCalibrationStyle(isoCode);
	const string calibMaturity = mcPathConfig->getCalibrationMaturity(isoCode);
	const string calibMaturityCMS = mcPathConfig->getCalibrationMaturityCMS(isoCode);

    // set all the values of pAssetLiborRates
    pAssetLiborRates->setSRMRatesLiborDiffuse(
        randomIdx, 
        today, 
        ratesLiborUtil, 
        true,   // saveDomLnMoney: was "i == 0 && sv->mapIRGen.size() > 1,"
                // Could be better by doing only equity ccys and dom
        true,   // saveSigmaR: for the time being save 
                // it always; but we would need it only for credit 
        mcPathConfig->getNumSDMaxEffRateIR(),
        mcPathConfig->getNumSDMinEffRateIR(),
        vector<double>(),//TODO: pastValues,
        irFxCorr,  
		mcPathConfig->calibrateAgainstSwaptionVols(name),
		calibStyle,
		calibMaturity,
		calibMaturityCMS
		);
    ratesLiborUtil = SRMRatesLiborUtilSP(0); //FIXME: doublecheck
}

void SRMGenDiffusibleIRLibor::pushbackICERecord(StateVarDBase& svdbICE, vector<ICE>& vICE)
{
    vICE.push_back(ICE(ratesLiborUtil, pAssetLiborRates));
    // collect the state variables
    StateVariableCollectorSP svCollector(new StateVariableCollector());
    vICE.back().collectStateVars(svCollector);
    IElemStateVariableGenArray stateVarGenArray(svCollector->getElemStateVarGens());

    // get the discount factor SVs
    vector<const SVGenDiscFactor*> icedfs(
        filterPointers<const IElemStateVariableGen*, 
                       const SVGenDiscFactor*>(stateVarGenArray));

    // and add them to the array
    for(size_t i=0; i<icedfs.size(); ++i)
    {
        SVDiscFactorSP dfSVSP (getDiffusibleIR()->createDFSV(icedfs[i]->getDates(), false));
        svdbICE.append(icedfs[i], dfSVSP);
    }

    // get the expected discount factor SVs
    vector<const SVGenExpectedDiscFactor*> iceexpDfs(
        filterPointers<const IElemStateVariableGen*, 
                       const SVGenExpectedDiscFactor*>(stateVarGenArray));

    for (size_t j = 0; j < iceexpDfs.size(); j++)
    {
        SVExpectedDiscFactorSP edfSVSP(
            getDiffusibleIR()->createExpDiscFactorSV(iceexpDfs[j]->getPVDate(),
                                            iceexpDfs[j]->getDates(),
                                            iceexpDfs[j]->logRequired(),
                                            iceexpDfs[j]->getYieldCurve(), false));
        svdbICE.append(iceexpDfs[j], edfSVSP);
    }
}

DateTimeArray SRMGenDiffusibleIRLibor::getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
                                                     const DateTime& start, 
                                                     const DateTime& finish)
{
    if (mcPathConfig->calibrateAgainstSwaptionVols(name))
    {
        return CriticalDateCollector::collectVolDates(
                ycs.first, IRVolBase::TYPE, start, finish);
        // ICE: will need to store bmDates or similar. Note that
        // there may be a slight discrepancy here between bmDates
        // and dates returned from getBMDetails for FIX style
        // calibration (latter might be curtailed)
        // expDates = DateTime::merge(expDates, bmDates);
    }
    return DateTimeArray();
}


/************************************************************************/
/*  Rates Deterministic-related block                                             */
/************************************************************************/

void SRMGenDiffusibleIRDeterm::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                       const DateTime&        today,
                                       DateTimeArrayConstSP   simDates)
{
    CVolProcessedSP processedVol(
        mcPathConfig->getProcessedIRVol(name, 
        ycs.second)); // why second ?
    // now we can build the SRMRatesLiborUtil object for this ccy
    const string& isoCode = name; 

    ratesDetUtil = SRMRatesDetermUtilSP(
        new SRMRatesDetermUtil(today,
        mcPathConfig->getNumIRFactors(isoCode),
        mcPathConfig->getIRModelParams(isoCode),
        mcPathConfig->getIRSmileParams(isoCode),
        processedVol, 
        ycs.first, //discYC
        ycs.second, // diffYC
        mcPathConfig->skipNegIRVols(isoCode),
        mcPathConfig->getFlatVolIr(),
        mcPathConfig->getChoiceCutoff(),
        mcPathConfig->getCutoffValue(),
        mcPathConfig->getCorrelationSwapStart(isoCode),
        mcPathConfig->getCorrelationSwapMat(isoCode),
        mcPathConfig->getCorrelationSwapDCC(isoCode),
        mcPathConfig->getCorrelationSwapFreq(isoCode)));

    ratesDetUtil->setTimeLine(simDates);
}

void SRMGenDiffusibleIRDeterm::initializeDiffusibleAssets()
{
    pAssetDetRates = SRMRatesDetermDiffuseSP(new SRMRatesDetermDiffuse());
}

void SRMGenDiffusibleIRDeterm::setSRMDiffusibleAsset(
    const DateTime&     today,              // base date of the run
    const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
    const IPastValues*  ,     // historic values
    DependenceSP        dependence)         // factorized correlations
{
    vector<double> irFxCorr(nFactors);

    if (pQMCGenDiffusibleFX) {  
        // there is FX asset too, so
        // look up the extended correlations of IR vs FX
        for(int nnf = 0; nnf<nFactors; ++nnf)
            irFxCorr[nnf] = dependence->getCorrelationIJ(pQMCGenDiffusibleFX->getRandomIdx(), 
            randomIdx+nnf, 0);
    }

    const string& isoCode = name; 
    const string calibStyle = mcPathConfig->getCalibrationStyle(isoCode);
    const string calibMaturity = mcPathConfig->getCalibrationMaturity(isoCode);
    const string calibMaturityCMS = mcPathConfig->getCalibrationMaturityCMS(isoCode);


    // set all the values of pAssetLiborRates
    pAssetDetRates->setSRMRatesDetermDiffuse(ratesDetUtil, today);

    ratesDetUtil = SRMRatesDetermUtilSP(0); //FIXME: doublecheck
}


DateTimeArray SRMGenDiffusibleIRDeterm::getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
                                                      const DateTime& start, 
                                                      const DateTime& finish)
{
    if (mcPathConfig->calibrateAgainstSwaptionVols(name))
    {
        return CriticalDateCollector::collectVolDates(ycs.first, IRVolBase::TYPE, start, finish);
        // ICE: will need to store bmDates or similar. Note that
        // there may be a slight discrepancy here between bmDates
        // and dates returned from getBMDetails for FIX style
        // calibration (latter might be curtailed)
        // expDates = DateTime::merge(expDates, bmDates);
    }
    return DateTimeArray();
}


vector<double> SRMGenDiffusibleIRDeterm::expandTier2Betas(double beta) const
{
    return vector<double>(0);
}

/************************************************************************/
/*  SRM FX-related block                                                    */
/************************************************************************/

void SRMGenDiffusibleFX::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                    const DateTime&        today,
                                    DateTimeArrayConstSP   simDates)
{
    // must have done domestic ccy before calling this
    assert(pDomQMCGenDiffusibleIR);
    assert(pForQMCGenDiffusibleIR);
    
    const string& forCcy = ccy;
        
    // look up the relevant model parameters
    const string& bootStrapMode = mcPathConfig->getFXVolBootStrapMode(forCcy);
    double fxCutOffLevel = mcPathConfig->getFXCutOffLevel(forCcy);
    const string& fx = name; //fxAsset->getName();

    double corrFxFor;
    CorrelationCommonSP corrFxForSP = CorrelationCommon::lookupCorrelation(
        mcPathConfig->getCorrObjArray(), 
        fx, 
        forCcy, 
        mcPathConfig->isStrictCorr());
    if (corrFxForSP.get() && Correlation::TYPE->isInstance(corrFxForSP.get()))
    {
        CorrelationSP p = CorrelationSP::dynamicCast(corrFxForSP);
        corrFxFor = p.get() ? p->getCorrelation() : 0;
    }
    else 
        corrFxFor = 0;

    double corrFxDom ;
    CorrelationCommonSP corrFxDomSP = CorrelationCommon::lookupCorrelation(
        mcPathConfig->getCorrObjArray(), 
        fx,
        pDomQMCGenDiffusibleIR->name, 
        mcPathConfig->isStrictCorr());
    if (corrFxDomSP.get() && Correlation::TYPE->isInstance(corrFxDomSP.get()))
    {
        CorrelationSP p = CorrelationSP::dynamicCast(corrFxDomSP);
        corrFxDom = p.get() ? p->getCorrelation() : 0;
    }
    else
        corrFxDom = 0;

    double corrDomFor;
    CorrelationCommonSP corrDomForSP = CorrelationCommon::lookupCorrelation(
        mcPathConfig->getCorrObjArray(), 
        pDomQMCGenDiffusibleIR->name, 
        forCcy, 
        mcPathConfig->isStrictCorr());
    if (corrDomForSP.get() && Correlation::TYPE->isInstance(corrDomForSP.get()))
    {
        CorrelationSP p = CorrelationSP::dynamicCast(corrDomForSP);
        corrDomFor = p.get() ? p->getCorrelation() : 0;
    }
    else
        corrDomFor = 0;

    string fxVolType = mcPathConfig->getFXVolType(forCcy);
    SRMRatesUtilSP forUtil( dynamic_cast<SRMRatesUtil*>(pForQMCGenDiffusibleIR->getRatesUtil().get()) );
    SRMRatesUtilSP domUtil( dynamic_cast<SRMRatesUtil*>(pDomQMCGenDiffusibleIR->getRatesUtil().get()) );
    QLIB_VERIFY( forUtil.get() != 0 && domUtil.get() != 0, "SRM FX requires SRMFXUtil objects");
    fxUtil = SRMFXUtilSP(
                        new SRMFXUtil(
                            forUtil, // foreign
                            domUtil, // domestic
                            fxAsset,
                            corrFxFor, 
                            corrFxDom, 
                            corrDomFor,
                            bootStrapMode, 
                            fxCutOffLevel,
                            fxVolType) );

    fxUtil->setTimeLine(simDates);
    fxUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());

}


void SRMGenDiffusibleFX::initializeDiffusibleAssets() 
{
    assert(pDomQMCGenDiffusibleIR);            // should be set
    assert(pDomQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already
    assert(pForQMCGenDiffusibleIR);                 // should be set
    assert(pForQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already

    // For now we only allow the Rates HJM model or deterministic rates
    SRMGenDiffusibleIRHJM* pForSRMGenDiffusibleIRHJM = dynamic_cast<SRMGenDiffusibleIRHJM*>(pForQMCGenDiffusibleIR);
    SRMGenDiffusibleIRDeterm* pForSRMGenDiffusibleIRDet = dynamic_cast<SRMGenDiffusibleIRDeterm*>(pForQMCGenDiffusibleIR);

    if (!pForSRMGenDiffusibleIRHJM && !pForSRMGenDiffusibleIRDet)
        throw ModelException("SRMGenDiffusibleFX::initializeDiffusibleAssets()", (string) "the FX " + name + 
        " cannot use SRM model when underlying foreign currency "+pForQMCGenDiffusibleIR->name+" is not using it");
    SRMGenDiffusibleIRHJM* pDomSRMGenDiffusibleIRHJM = dynamic_cast<SRMGenDiffusibleIRHJM*>(pDomQMCGenDiffusibleIR);
    SRMGenDiffusibleIRDeterm* pDomSRMGenDiffusibleIRDet = dynamic_cast<SRMGenDiffusibleIRDeterm*>(pDomQMCGenDiffusibleIR);

    if (!pDomSRMGenDiffusibleIRHJM && !pDomSRMGenDiffusibleIRDet)
        throw ModelException("SRMGenDiffusibleFX::initializeDiffusibleAssets()", (string) "the FX " + name + 
        " cannot use SRM model when domestic currency "+pDomQMCGenDiffusibleIR->name+" is not using it");

    pAssetFX = SRMFXDiffuseSP(new SRMFXDiffuse(pDomQMCGenDiffusibleIR->getDiffusibleIR(), pForQMCGenDiffusibleIR->getDiffusibleIR() ));
}


void SRMGenDiffusibleFX::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  pastValues,         // historic values
                DependenceSP        dependence)         // factorized correlations
{
    
    DoubleArray past;
    
    if (assetIdx >= 0)
    {
        DateTimeArray pastDates = today.getPastDates( pAssetFX->getAssetDates() );
        past =  pastValues->getPastValues(pastDates, assetIdx, today);
    }
    
    pAssetFX->setSRMFX(
        randomIdx,
        today,
        fxUtil,
        vector<double>(past.begin(), past.end())); 
}

DateTimeArray SRMGenDiffusibleFX::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
                                                   const DateTime& start, 
                                                   const DateTime& finish)
{
    return CriticalDateCollector::collectVolDates(
                fxAsset, FXVolBase::TYPE, start, finish);
}

vector<double> SRMGenDiffusibleFX::expandAssetToIRCorrs(
    QMCGenDiffusibleIR* irAsset, 
    double corr)
{
    SRMGenDiffusibleIR* irAssetSRM = dynamic_cast<SRMGenDiffusibleIR*>(irAsset);
    QLIB_VERIFY(irAssetSRM != 0, string() + "Expect type SRMGenDiffusibleIR but type '"
        + typeid(*irAsset).name() + "' was used.");

    return SRMCorrelation::ExpandAssetAssetCorrs(irAssetSRM, corr, SRMCorrelation::PC); // 1-factor to ir
}
vector<double> SRMGenDiffusibleFX::expandAssetToENCorrs(
    QMCGenDiffusibleEnergy* enAsset, 
    double corr)
{
    SRMGenDiffusibleEnergy* enAssetSRM = dynamic_cast<SRMGenDiffusibleEnergy*>(enAsset);
    QLIB_VERIFY(enAssetSRM != 0, string() + "Expect type SRMGenDiffusibleEnergy but type '"
        + typeid(*enAsset).name() + "' was used.");

    // TODO: reverse case!
    return SRMCorrelation::ExpandAssetAssetCorrs(enAssetSRM, corr, SRMCorrelation::PC); // 1-factor to energy
}
vector<double> SRMGenDiffusibleFX::expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleFX::expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleFX::expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleFX::expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}


/************************************************************************/
/*  SRM Equity-related block                                                */
/************************************************************************/

void SRMGenDiffusibleEquity::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                        const DateTime&        today,
                                        DateTimeArrayConstSP   simDates)
{
    // must have done domestic ccy before calling this
    assert(pDomQMCGenDiffusibleIR);
    assert(pQMCGenDiffusibleIR); 

    // look up the relevant model parameters
    const string& bootStrapMode = mcPathConfig->getEqVolBootStrapMode(assetIdx);
    double eqCutOffLevel = mcPathConfig->getEqCutOffLevel(assetIdx);
    const string& eqName = eqAsset->getTrueName();

    double corrEqIr;
    CorrelationCommonSP corrEqIrSP = CorrelationCommon::lookupCorrelation(mcPathConfig->getCorrObjArray(), 
        eqName, 
        pDomQMCGenDiffusibleIR->name, //TODO: is this correct?
        mcPathConfig->isStrictCorr());
    if (corrEqIrSP.get() && Correlation::TYPE->isInstance(corrEqIrSP.get()))
    {
        CorrelationSP p = CorrelationSP::dynamicCast(corrEqIrSP);
        corrEqIr = p.get() ? p->getCorrelation() : 0;
    }
    else
        corrEqIr = 0;
        
    if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX != NULL) {
        //If the base ccy of the eq has an associated FX, then this means it is 
        //a 'foreign' currency, so vanilla ccyTreatment doesn't make sense.
        QLIB_VERIFY(ccyTreatment != ccyVanilla, "Inconsistent currency treatment.");
    }

    // TO DO: we need to ensure that SRMEQUtil below ignores the corr values if they are already calibrated. 
    double corrEqFX = 0.0;
    if (ccyTreatment != ccyVanilla) {
        // Gather info for quanto adjustment
        const string& fxName = pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->name;   
        CorrelationCommonSP corrEqFXSP = CorrelationCommon::lookupCorrelation(
            mcPathConfig->getCorrObjArray(), 
            fxName, 
            eqName, 
            mcPathConfig->isStrictCorr());
        if (corrEqFXSP.get() && Correlation::TYPE->isInstance(corrEqFXSP.get()))
        {
            CorrelationSP p = CorrelationSP::dynamicCast(corrEqFXSP);
            corrEqFX = p.get() ? p ->getCorrelation() : 0;
        }
        else { 
            corrEqFX = 0;
        }

        QLIB_VERIFY(mcPathConfig->getEqProtAdjMethod() == "DEFAULT", 
            "QMCGenDiffusibleEquity::createSRMUtil"
            "eqProtAdjMethod must be DEFAULT for now");
    }

    SRMRatesUtilSP ratesUtil( dynamic_cast<SRMRatesUtil*>(pQMCGenDiffusibleIR->getRatesUtil().get()) );
    QLIB_VERIFY(!!ratesUtil, "SRM Equity requires a SRMRatesUtil object.");

    equityUtil = 
        SRMEquityUtilSP(
            new SRMEquityUtil(
                ratesUtil, 
                eqAsset,
                corrEqIr,
                ccyTreatment,
                corrEqFX, 
                bootStrapMode, 
                eqCutOffLevel));

    equityUtil->setTimeLine(simDates);   
    equityUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());

}

void SRMGenDiffusibleEquity::initializeDiffusibleAssets()
{
    assert(pQMCGenDiffusibleIR);                 // should be set
    assert(pQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already
    
    string method = "SRMGenDiffusibleEquity::initializeDiffusibleAssets()";

    // For now we only allow the Rates HJM model or deterministic rates
    SRMGenDiffusibleIRHJM* pSRMGenDiffusibleIRHJM = dynamic_cast<SRMGenDiffusibleIRHJM*>(pQMCGenDiffusibleIR);
    SRMGenDiffusibleIRDeterm* pSRMGenDiffusibleIRDet = dynamic_cast<SRMGenDiffusibleIRDeterm*>(pQMCGenDiffusibleIR);

    if (!pSRMGenDiffusibleIRHJM && !pSRMGenDiffusibleIRDet)
        throw ModelException(method, (string) "the equity " + name + 
        " requires an HJM or deterministic rates model type");

    switch(diffusionStyle)
    {
    case eqDEFAULT:
    case eqEULER:
        pAssetEquity = SRMEquityDiffuseSP(new SRMEquityDiffuseEuler(pQMCGenDiffusibleIR->getDiffusibleIR()));
        break;
    case eqSECOND_ORDER:
        pAssetEquity = SRMEquityDiffuseSP(new SRMEquityDiffuseSecondOrder(pQMCGenDiffusibleIR->getDiffusibleIR()));
        break;
    case eqSTRICT_SECOND_ORDER:
        pAssetEquity = SRMEquityDiffuseSP(new SRMEquityDiffuseStrictSecondOrder(pQMCGenDiffusibleIR->getDiffusibleIR()));
        break;
        /*
    case DiffusionStyleEnum::MAPPING_METHOD:
        pAssetEquity = SRMEquityDiffuseSP(new SRMEquityDiffuseMappingMethod(pQMCGenDiffusibleIR->getDiffusibleIR()));
        break;
        */
    default:
        throw ModelException(method, Format::toString("Unsupported diffusion style %i", (int)diffusionStyle));
    }    
}

void SRMGenDiffusibleEquity::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  pastValues,         // historic values
                DependenceSP        dependence)         // factorized correlations
{
    DoubleArray past;

    if (assetIdx >= 0)
    {
        DateTimeArray pastDates = today.getPastDates( pAssetEquity->getAssetDates() );
        past =  pastValues->getPastValues(pastDates, assetIdx, today);
    }
 
    pAssetEquity->setSRMEquityDiffuse(randomIdx,
                                        today,
                                        equityUtil,
                                        vector<double>(past.begin(), past.end()),
                                        mcPathConfig->getTimePointsPerYear()); 
}

DateTimeArray SRMGenDiffusibleEquity::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
                                                       const DateTime& start, 
                                                       const DateTime& finish)
{
    // Smile dates, and dividend dates
    return SRMEquityDiffuse::getCriticalDates(eqAsset, start, finish);
}

vector<double> SRMGenDiffusibleEquity::expandAssetToIRCorrs(
    QMCGenDiffusibleIR* irAsset, 
    double corr)
{
    SRMGenDiffusibleIR* irAssetSRM = dynamic_cast<SRMGenDiffusibleIR*>(irAsset);
    QLIB_VERIFY(irAssetSRM != 0, string() + "Expect type SRMGenDiffusibleIR but type '"
        + typeid(*irAsset).name() + "' was used.");

    return SRMCorrelation::ExpandAssetAssetCorrs(irAssetSRM, corr, SRMCorrelation::PC); // 1-factor to ir
}
vector<double> SRMGenDiffusibleEquity::expandAssetToENCorrs(
    QMCGenDiffusibleEnergy* enAsset, 
    double corr)
{
    SRMGenDiffusibleEnergy* enAssetSRM = dynamic_cast<SRMGenDiffusibleEnergy*>(enAsset);
    QLIB_VERIFY(enAssetSRM != 0, string() + "Expect type SRMGenDiffusibleEnergy but type '"
        + typeid(*enAsset).name() + "' was used.");

    // TODO: reverse case!
    return SRMCorrelation::ExpandAssetAssetCorrs(enAssetSRM, corr, SRMCorrelation::PC); // 1-factor to energy
}
vector<double> SRMGenDiffusibleEquity::expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleEquity::expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleEquity::expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleEquity::expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}


/************************************************************************/
/*  SRM Credit block                                                    */
/************************************************************************/
vector<double> SRMGenDiffusibleCredit::expandAssetToIRCorrs(
    QMCGenDiffusibleIR* irAsset, 
    double corr)
{
    SRMGenDiffusibleIR* irAssetSRM = dynamic_cast<SRMGenDiffusibleIR*>(irAsset);
    QLIB_VERIFY( irAssetSRM != 0, string() + "Expect SRMGenDiffusibleIR but '"
        + typeid(*irAsset).name() + "' was used");

    return SRMCorrelation::ExpandAssetAssetCorrs(irAssetSRM, corr, SRMCorrelation::PC);
}

vector<double> SRMGenDiffusibleCredit::expandAssetToENCorrs(
    QMCGenDiffusibleEnergy* enAsset, 
    double corr)
{
    SRMGenDiffusibleEnergy* enAssetSRM = dynamic_cast<SRMGenDiffusibleEnergy*>(enAsset);
    QLIB_VERIFY( enAssetSRM != 0, string() + "Expect SRMGenDiffusibleEnergy but '"
        + typeid(*enAsset).name() + "' was used");

    return SRMCorrelation::ExpandAssetAssetCorrs(enAssetSRM, corr, SRMCorrelation::PC);
}

vector<double> SRMGenDiffusibleCredit::expandAssetToCRCorrs(
    QMCGenDiffusibleCredit* crAsset, 
    double corr)
{
    return vector<double>(1, corr); // 1 factor vs 1 factor, no need to expand corr
}

vector<double> SRMGenDiffusibleCredit::expandAssetToEQCorrs(
    QMCGenDiffusibleEquity* eqAsset, 
    double corr)
{
    return vector<double>(1, corr); // 1 factor vs 1 factor, no need to expand corr
}

vector<double> SRMGenDiffusibleCredit::expandAssetToFXCorrs(
    QMCGenDiffusibleFX* fxAsset, 
    double corr)
{
    return vector<double>(1, corr); // 1 factor vs 1 factor, no need to expand corr
}

vector<double> SRMGenDiffusibleCredit::expandAssetToBSCorrs(
    QMCGenDiffusibleBasisSpread* bsAsset, 
    double corr)
{
    return vector<double>(1, corr); // 1 factor vs 1 factor, no need to expand corr
}


/************************************************************************/
/*  Credit HJM-related block                                            */
/************************************************************************/

void SRMGenDiffusibleCreditHJM::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                        const DateTime&        today,
                                        DateTimeArrayConstSP   simDates)
{
    creditHJMUtil = SRMCreditHJMUtilSP(new SRMCreditHJMUtil(today, mcPathConfig->getCRSmileParams(name), cds));

    creditHJMUtil->setTimeLine(simDates);
    creditHJMUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleCreditHJM::initializeDiffusibleAssets()
{
    static const string method("QMCGenDiffusibleCredit::initializeDiffusibleAssets()");

    assert(pQMCGenDiffusibleIR);                 // should be set
    assert(pQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already

    SRMGenDiffusibleIRHJM* pSRMGenDiffusibleIRHJM = dynamic_cast<SRMGenDiffusibleIRHJM*>(pQMCGenDiffusibleIR);
    if (!pSRMGenDiffusibleIRHJM)
        throw ModelException(method, (string) "the credit spread " + name + 
            " cannot use " + MCPathConfigSRM::CREDIT_MODEL_TYPE_HJM + " model when underlying currency " + 
            pQMCGenDiffusibleIR->name+" is not using " + MCPathConfigSRM::RATES_MODEL_TYPE_HJM + " model");

    const int numIRFactors = pQMCGenDiffusibleIR->numFactors();

    switch(numIRFactors)
    {
    case 1: pAssetCreditHJM = SRMCreditHJMDiffuseSP(new SRMCreditHJM1F(pSRMGenDiffusibleIRHJM->pAssetHJMRates));
            break;
    case 2: pAssetCreditHJM = SRMCreditHJMDiffuseSP(new SRMCreditHJM2F(pSRMGenDiffusibleIRHJM->pAssetHJMRates));
            break;
    case 3: pAssetCreditHJM = SRMCreditHJMDiffuseSP(new SRMCreditHJM3F(pSRMGenDiffusibleIRHJM->pAssetHJMRates));
            break;
    default: throw ModelException(method, 
                "Number of IR factors > 3 not supported");
    }
}

void SRMGenDiffusibleCreditHJM::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  pastValues,         // historic values
                DependenceSP        dependence)         // factorized correlations
{

    ///*    vector<double> pastValues(
    //        getSurvivalDiscFactorSVPast(
    //            pastPathGenerator,
    //            mcPathConfig,
    //            mcSurvivalDiscFactors,
    //            today,
    //            cmptSVSIndexes,
    //            cmptSVSBegin,
    //            svsIndexes,
    //            svsBegin));*/
    double NbSigmasMax = mcPathConfig->getNumSDMaxEffRateCR();
    double NbSigmasMin = mcPathConfig->getNumSDMinEffRateCR();

    double corrCRFX = 0.0;
    if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX) {
        corrCRFX = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->getRandomIdx(), 
                                                randomIdx, 0);
    }

    int nFactors = pQMCGenDiffusibleIR->numFactors();
    vector<double> corrCRIR(nFactors);
    for (int j = 0; j < nFactors; j++) {
        corrCRIR[j] = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->getRandomIdx()+j, 
                                                   randomIdx, 0);
    }

    pAssetCreditHJM->setSRMCreditHJMDiffuse(
        randomIdx,
        today,
        creditHJMUtil,
        NbSigmasMax,
        NbSigmasMin,
        corrCRIR,
        corrCRFX,
        vector<double>() /*pastValues*/);
    creditHJMUtil = SRMCreditHJMUtilSP(); // creditHJMUtil is no longer needed

}

DateTimeArray SRMGenDiffusibleCreditHJM::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
                                                       const DateTime& /*start*/, 
                                                       const DateTime& /*finish*/)
{
    /* Credit contribution -- nthg so far, since there is no calibration
        native Srm3 as of March 2005 does not have anything in this regard */
    return DateTimeArray();
}


/************************************************************************/
/*  Credit CIR-related block                                            */
/************************************************************************/

void SRMGenDiffusibleCreditCIR::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                        const DateTime&        today,
                                        DateTimeArrayConstSP   simDates)
{
    creditCIRUtil = SRMCreditCIRUtilSP(new SRMCreditCIRUtil(today, mcPathConfig->getCRSmileParams(name), cds));
    creditCIRUtil->setTimeLine(simDates);
    creditCIRUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleCreditCIR::initializeDiffusibleAssets()
{
    assert(pQMCGenDiffusibleIR);                 // should be set
    assert(pQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already

    pAssetCreditCIR = SRMCreditCIRDiffuseSP(new SRMCreditCIRDiffuse(pQMCGenDiffusibleIR->getDiffusibleIR(), true));
 }

void SRMGenDiffusibleCreditCIR::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  pastValues,         // historic values
                DependenceSP        dependence)         // factorized correlations
{
    ///*    vector<double> pastValues(
    //        getSurvivalDiscFactorSVPast(
    //            pastPathGenerator,
    //            mcPathConfig,
    //            mcSurvivalDiscFactors,
    //            today,
    //            cmptSVSIndexes,
    //            cmptSVSBegin,
    //            svsIndexes,
    //            svsBegin));*/

    double corrCRFX = 0.0;
    if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX) {
        corrCRFX = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->getRandomIdx(), 
                                                randomIdx, 0);
    }

    pAssetCreditCIR->setSRMCreditCIRDiffuse(
        randomIdx,
        today,
        creditCIRUtil,
        corrCRFX,
        vector<double>() /*pastValues*/);

    creditCIRUtil = SRMCreditCIRUtilSP(0); // creditCIRUtil is no longer needed
}

DateTimeArray SRMGenDiffusibleCreditCIR::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
                                                       const DateTime& /*start*/, 
                                                       const DateTime& /*finish*/)
{
    /* Credit contribution -- nthg so far, since there is no calibration */
    return DateTimeArray();
}

/************************************************************************/
/*  Credit Libor-related block                                          */
/************************************************************************/

void SRMGenDiffusibleCreditLibor::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                          const DateTime&        today,
                                          DateTimeArrayConstSP   simDates)
{
    creditLiborUtil = SRMCreditLiborUtilSP(new SRMCreditLiborUtil(today, mcPathConfig->getCRSmileParams(name), cds));
    creditLiborUtil->setTimeLine(simDates);
    creditLiborUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleCreditLibor::initializeDiffusibleAssets()
{
    static const string method("SRMGenDiffusibleCreditLibor::initializeDiffusibleAssets()");

    assert(pQMCGenDiffusibleIR);						    // should be set
    assert(pQMCGenDiffusibleIR->getDiffusibleIR().get());	// should be initialized already

    pAssetCreditLibor = SRMCreditLiborDiffuseSP(new SRMCreditLiborDiffuse(pQMCGenDiffusibleIR->getDiffusibleIR()));
}

void SRMGenDiffusibleCreditLibor::setSRMDiffusibleAsset(
    const DateTime&     today,              // base date of the run
    const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
    const IPastValues*  pastValues,         // historic values
    DependenceSP        dependence)         // factorized correlations
{
    ///*    vector<double> pastValues(
    //        getSurvivalDiscFactorSVPast(
    //            pastPathGenerator,
    //            mcPathConfig,
    //            mcSurvivalDiscFactors,
    //            today,
    //            cmptSVSIndexes,
    //            cmptSVSBegin,
    //            svsIndexes,
    //            svsBegin));*/

    double corrCRFX = 0.0;
    if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX) {
        corrCRFX = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->getRandomIdx(), 
            randomIdx, 0);
    }

    pAssetCreditLibor->setSRMCreditLiborDiffuse(     // is this correct?
        randomIdx,
        today,
        creditLiborUtil,
        corrCRFX,
        vector<double>());
    creditLiborUtil = SRMCreditLiborUtilSP(0); // creditLiborUtil is no longer needed
}

DateTimeArray SRMGenDiffusibleCreditLibor::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
                                                         const DateTime& /*start*/, 
                                                         const DateTime& /*finish*/)
{
    /* Credit contribution -- nthg so far, since there is no calibration */
    return DateTimeArray();
}

/************************************************************************/
/*  SRM Energy block                                                    */
/************************************************************************/
/**
 * 
 * @param mcPathConfig 
 * @param today 
 * @param simDates 
 */
void SRMGenDiffusibleEnergy::createSRMUtil(
	MCPathConfigSRM * mcPathConfig,
	const DateTime & today,
	DateTimeArrayConstSP simDates
	)
{
    if (isTier1) {
        enrgUtil = SRMEnergyUtilCreate(
		    modelType, 
		    nFactors,
		    today, 
		    futureCurve,
		    mcPathConfig->getEnergyCorrInstrStart(name),
		    mcPathConfig->getEnergyCorrInstrMaturity(name)
		    );
    }
    else { // tier 2: create generic energy util object
        enrgUtil = SRMEnergyUtilBaseSP(
            new SRMEnergyUtilTier2(
                nFactors,
                today,
                futureCurve,
                mcPathConfig->getEnergyCorrInstrStart(name), // note: to be decide if we really need this
                mcPathConfig->getEnergyCorrInstrMaturity(name)  // note: to be decide if we really need this
            )
        );
    }

    enrgUtil->setTimeLine(simDates);
    // moment matching is not implemented to energy futrues yet
    //enrgUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleEnergy::expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const
{
    DoubleMatrix assetTMtx(enrgUtil->getMatrixAByPCA());
    int randomIndex = getRandomIdx();
    outputSparseCorrMatrix.selfMultiplyToSmallmtx(assetTMtx, randomIndex, randomIndex);
}

vector<double> SRMGenDiffusibleEnergy::expandTier2Betas(double beta) const
{
    vector<double> tmpCorrs;
    if ( numFactors() > 1 ) {// multifactor Tier1
        // split the one-factor vs ir correlation
        tmpCorrs = SRMCorrelation::enrgToAsset(
            beta,
            SRMCorrelation::EXPONENTIAL_FACTORS,
            *enrgUtil.get());
    } else {
        tmpCorrs = vector<double>(1,beta);
    }
    return tmpCorrs;
}

void SRMGenDiffusibleEnergy::initializeDiffusibleAssets()
{
	static const string method("SV::SRMGenDiffusibleEnergy::initializeDiffusibleAssets()");

	assert(pQMCGenDiffusibleIR);                 // should be set
	assert(pQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already
	SRMGenDiffusibleIRHJM* pSRMGenDiffusibleIRHJM = dynamic_cast<SRMGenDiffusibleIRHJM*>(pQMCGenDiffusibleIR);
	if (!pSRMGenDiffusibleIRHJM)
		throw ModelException(method, (string) "the energy future " + name + 
		    " cannot use SRM model when underlying currency "+pQMCGenDiffusibleIR->name+" is not using it");

    if (isTier1) {
        pAssetEnrg = SRMEnergyDiffuseCreate(modelType, pSRMGenDiffusibleIRHJM->pAssetHJMRates);
    }
    else { // is tier 2
        //pAssetEnrg = QMCEnergyDiffuseSP(new SRMEnergyTier2Diffuse(pSRMGenDiffusibleIRHJM->pAssetHJMRates));
        assert(pParentQMCGenDiffusibleEnergy); // should be set
        assert(pParentQMCGenDiffusibleEnergy->getDiffusibleEnergy().get()); // should be initialized already
        SRMGenDiffusibleEnergy* pParentSRMGenDiffusibleEnergy = dynamic_cast<SRMGenDiffusibleEnergy*>(pParentQMCGenDiffusibleEnergy);
        pAssetEnrg = QMCEnergyDiffuseSP(
            new SRMEnergyTier2Diffuse(
                pSRMGenDiffusibleIRHJM->pAssetHJMRates,
                pParentSRMGenDiffusibleEnergy->pAssetEnrg
            )
        );
    }
}

void SRMGenDiffusibleEnergy::setSRMDiffusibleAsset(
	const DateTime&     today,              // base date of the run
	const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
	const IPastValues*  pastValues,         // historic values
	DependenceSP        dependence)         // factorized correlations
{
	// setup energy-fx correlation vector for quanto adjustment, 
	// only if necessary
	vector<double> corrEnrgFx(nFactors);
	if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX) {
		for (int i = 0; i < nFactors; ++i) {
			corrEnrgFx[i] = dependence->getCorrelationIJ(
				pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->getRandomIdx(), // fx index
				randomIdx + i, // enrg index
				0);
		}
	}

	pAssetEnrg->setQMCEnergyDiffuse(
		futureCurve,
		enrgUtil,
		randomIdx,
		today,
		corrEnrgFx,
		vector<double>() // past values
		);

/*    // Now ready to set up the link to parent energy for tier 2 assets
    if (!isTier1) { 
        assert(pParentQMCGenDiffusibleEnergy); // should be set
        assert(pParentQMCGenDiffusibleEnergy->getDiffusibleEnergy().get()); // should be initialized already
       	SRMGenDiffusibleEnergy* pParentSRMGenDiffusibleEnergy = dynamic_cast<SRMGenDiffusibleEnergy*>(pParentQMCGenDiffusibleEnergy);
        pAssetEnrg->setParent(pParentSRMGenDiffusibleEnergy->pAssetEnrg); 
    }*/
}

DateTimeArray SRMGenDiffusibleEnergy::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
													   const DateTime& start, 
													   const DateTime& finish)
{
	// Energy contribution
	/*
	// Note: energy uses fix-style calibration: the volatility
	// benchmark dates are fixed and we will add them directly
	// to the timeline.
	const DateTimeArray & volDatesAll = futureCurve->getVolBenchmarkDates();
	DateTimeArray volDates = DateTime::getInclusiveDates(start, finish, volDatesAll);
	
	return volDates;
	*/

	/* Energy contribution -- nthg so far, since there is no calibration
	native Srm3 as of March 2006 does not have anything in this regard */
	return DateTimeArray();
}

vector<double> SRMGenDiffusibleEnergy::expandAssetToIRCorrs(
    QMCGenDiffusibleIR* irAsset, 
    double corr)
{
    SRMGenDiffusibleIR* irAssetSRM = dynamic_cast<SRMGenDiffusibleIR*>(irAsset);
    QLIB_VERIFY( irAssetSRM != 0, string() + "Expect SRMGenDiffusibleIR but '"
        + typeid(*irAsset).name() + "' was used");

    // reverse ordering case here: ir is primary asset here
    // ir to energy:
    return SRMCorrelation::ExpandAssetAssetCorrs(irAssetSRM, this, corr, SRMCorrelation::PC, SRMCorrelation::PC);
}

vector<double> SRMGenDiffusibleEnergy::expandAssetToENCorrs(
    QMCGenDiffusibleEnergy* enAsset, 
    double corr)
{
    SRMGenDiffusibleEnergy* enAssetSRM = dynamic_cast<SRMGenDiffusibleEnergy*>(enAsset);
    QLIB_VERIFY( enAssetSRM != 0, string() + "Expect SRMGenDiffusibleEnergy but '"
        + typeid(*enAsset).name() + "' was used");

    return SRMCorrelation::ExpandAssetAssetCorrs(enAssetSRM, this, corr, SRMCorrelation::PC);
}

vector<double> SRMGenDiffusibleEnergy::expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to energy
}
vector<double> SRMGenDiffusibleEnergy::expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to energy
}
vector<double> SRMGenDiffusibleEnergy::expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to energy
}
vector<double> SRMGenDiffusibleEnergy::expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) {
    return SRMCorrelation::ExpandAssetAssetCorrs(this, corr, SRMCorrelation::PC); // 1-factor to energy
}


/************************************************************************/
/*  SRM Basis-related block                                             */
/************************************************************************/

void SRMGenDiffusibleBasisSpread::createSRMUtil(MCPathConfigSRM*      mcPathConfig,
                                            const DateTime&        today,
                                            DateTimeArrayConstSP   simDates)
{
    basisUtil = SRMBasisSpreadUtilSP(new SRMBasisSpreadUtil(today, basisCurve));                                 
    basisUtil->setTimeLine(simDates);
    // momentmatching is not implement for BasisSpread.
    // basisUtil->setMomentMatchingFlag(mcPathConfig->getMomentMatchingFlag());
}

void SRMGenDiffusibleBasisSpread::initializeDiffusibleAssets()
{
    static const string method("QMCGenDiffusibleCredit::initializeDiffusibleCredit()");

    assert(pQMCGenDiffusibleIR);                 // should be set
    assert(pQMCGenDiffusibleIR->getDiffusibleIR().get()); // should be initialized already
    assert(basisCurve.get());

    SRMGenDiffusibleIRHJM* pSRMGenDiffusibleIRHJM = dynamic_cast<SRMGenDiffusibleIRHJM*>(pQMCGenDiffusibleIR);
    if (!pSRMGenDiffusibleIRHJM)
        throw ModelException(method, (string) "the basis spread " + name + 
            " cannot use HJM model when underlying currency "+pQMCGenDiffusibleIR->name+" is not using it");

    SRMRatesHJMDiffuse* pDiffusibleHJM = dynamic_cast<SRMRatesHJMDiffuse*>(pSRMGenDiffusibleIRHJM->pAssetHJMRates.get());
    if (!pDiffusibleHJM)
        throw ModelException(method, (string) "the basis spread " + name + 
            " cannot use HJM model when underlying currency "+pQMCGenDiffusibleIR->name+" is not using it");

    SRMRatesHJMDiffuseSP spDiffusibleHJM(pDiffusibleHJM);

    const int numIRFactors = pQMCGenDiffusibleIR->numFactors();
    switch(numIRFactors)
    {
    case 1: pAsset = SRMBasisSpreadHJMDiffuseSP(new SRMBasisSpreadHJM1F(spDiffusibleHJM, basisCurve));
            break;
    case 2: pAsset = SRMBasisSpreadHJMDiffuseSP(new SRMBasisSpreadHJM2F(spDiffusibleHJM, basisCurve));
            break;
    case 3: pAsset = SRMBasisSpreadHJMDiffuseSP(new SRMBasisSpreadHJM3F(spDiffusibleHJM, basisCurve));
            break;
    default: throw ModelException(method, 
                "Number of IR factors > 3 not supported");
    }
}

void SRMGenDiffusibleBasisSpread::setSRMDiffusibleAsset(
                const DateTime&     today,              // base date of the run
                const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
                const IPastValues*  pastValues,         // historic values
                DependenceSP        dependence)         // factorized correlations
{

    double NbSigmasMax = mcPathConfig->getNumSDMaxEffRateSP();
    double NbSigmasMin = mcPathConfig->getNumSDMinEffRateSP();

    double corrSPFX = 0.0;
    if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX) {
        corrSPFX = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->getRandomIdx(), 
                                                randomIdx, 0);
    }

    int nFactors = pQMCGenDiffusibleIR->numFactors();
    vector<double> corrSPIR(nFactors);
    for (int j = 0; j < nFactors; j++) {
        corrSPIR[j] = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->getRandomIdx()+j, 
                                                   randomIdx, 0);
    }

    pAsset->setSRMBasisSpreadHJMDiffuse(
        randomIdx,
        today,
        basisUtil,
        NbSigmasMax,
        NbSigmasMin,
        corrSPIR,
        corrSPFX);
    basisUtil = SRMBasisSpreadUtilSP(0); // crUtil is no longer needed?

}

DateTimeArray SRMGenDiffusibleBasisSpread::getSRMCriticalDates(const MCPathConfigSRM* /*mcPathConfig*/,
                                                       const DateTime& /*start*/, 
                                                       const DateTime& /*finish*/)
{
    /* Spread contribution -- nthg so far, since there is no calibration
        native Srm3 as of March 2005 does not have anything in this regard */
    return DateTimeArray();
}

vector<double> SRMGenDiffusibleBasisSpread::expandAssetToIRCorrs(
    QMCGenDiffusibleIR* irAsset, 
    double corr)
{
    SRMGenDiffusibleIR* irAssetSRM = dynamic_cast<SRMGenDiffusibleIR*>(irAsset);
    QLIB_VERIFY(irAssetSRM != 0, string() + "Expect type SRMGenDiffusibleIR but type '"
        + typeid(*irAsset).name() + "' was used.");

    return SRMCorrelation::ExpandAssetAssetCorrs(irAssetSRM, corr, SRMCorrelation::PC); // 1-factor to ir
}
vector<double> SRMGenDiffusibleBasisSpread::expandAssetToENCorrs(
    QMCGenDiffusibleEnergy* enAsset, 
    double corr)
{
    SRMGenDiffusibleEnergy* enAssetSRM = dynamic_cast<SRMGenDiffusibleEnergy*>(enAsset);
    QLIB_VERIFY(enAssetSRM != 0, string() + "Expect type SRMGenDiffusibleEnergy but type '"
        + typeid(*enAsset).name() + "' was used.");

    // TODO: reverse case!
    return SRMCorrelation::ExpandAssetAssetCorrs(enAssetSRM, corr, SRMCorrelation::PC); // 1-factor to energy
}
vector<double> SRMGenDiffusibleBasisSpread::expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleBasisSpread::expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleBasisSpread::expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}
vector<double> SRMGenDiffusibleBasisSpread::expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr)
{
    return vector<double>(1, corr); // 1-factor to 1-factor
}


/************************************************************************/
/* QMC Credit CID Block (TODO: Why is QMC stuff this here?)             */
/************************************************************************/
// CID Generators implemented

/** this one is almost the same as the generator of CIR, with few distinctions */
void QMCGenDiffusibleCreditCID::initializeDiffusibleAssets()
{
    QLIB_VERIFY(pQMCGenDiffusibleIR!=NULL, "IR cannot be NULL");                 // should be set
    QLIB_VERIFY(pQMCGenDiffusibleIR->getDiffusibleIR().get()!=NULL, "IR should be initialized"); // should be initialized already

    pAssetCreditCIR = QMCCreditCIDDiffuseSP(new QMCCreditCIDDiffuse(pQMCGenDiffusibleIR->getDiffusibleIR(), isFullMC));
}
void QMCGenDiffusibleCreditCID::createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                                        const DateTime&        today,
                                        DateTimeArrayConstSP   simDates)
{
    // no util here - just direct construction of this object out of CIDparameters
}

void QMCGenDiffusibleCreditCID::setSRMDiffusibleAsset(
        const DateTime&        today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,       // to get max/min boundaries, etc
        const IPastValues*     pastValues,         // historic values
        DependenceSP           dependence)         // factorized correlations
{
    // This is specific to re-using CIR in Ian's implementation -- if we use different one, 
    // then no correlation exists!

    // The Nicolas' extention does not provide for non-zero correlation

    //double corrCRFX = 0.0;
    //if (pQMCGenDiffusibleIR->pQMCGenDiffusibleFX) {
    //    corrCRFX = dependence->getCorrelationIJ(pQMCGenDiffusibleIR->pQMCGenDiffusibleFX->getRandomIdx(), 
    //                                            randomIdx, 0);
    // }

    const ICIDParameters* pCIDParameters = mcPathConfig->getCIDParameters();

    QLIB_VERIFY(pCIDParameters != NULL, 
        "Cannot initialize CIR diffusion for CID-modeled credit -- CIDParameters are not defined in MCPathConfig.");

    CIDCommonRecordConstSP pCIRFactor;
    if (isCommonFactor)
        pCIRFactor = pCIDParameters->getCommonFactor();
    else
        pCIRFactor = pCIDParameters->getSingleNameRecord(name);


    double sigma = pCIRFactor->sigma;
    double meanRev = pCIRFactor->kappa;
    ExpiryArrayConstSP endsForTheta = pCIRFactor->endDatesForTheta;
    DateTimeArray  datesForTheta(endsForTheta->size());
    vector<double> logSurvivalProbs(endsForTheta->size());
    for(size_t i=0; i<logSurvivalProbs.size(); ++i)
    {
        datesForTheta[i] = (*endsForTheta)[i]->toDate(today);
        logSurvivalProbs[i] = log(pCIRFactor->calcProperSurvProbability(today, datesForTheta[i]));
    }

    // calibration of thetas happens entirely inside the diffuse object
    static_cast<QMCCreditCIDDiffuse*>(pAssetCreditCIR.get())->setQMCCreditCIDDiffuse(
        randomIdx,
        today,
        vector<double>(), // historic probabilities -- make no sense in credit context
        sigma,
        meanRev,
        datesForTheta,
        logSurvivalProbs
        );

}

void SRMGenDiffusibleEquity::setCurrencyTreatment(CurrencyTreatment _ccyTreatment)
{
    ccyTreatment = _ccyTreatment;    
}


