//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SimPathIInstrument.cpp
//
//   Description : Pseudo-instrument for RM applications
//
//   Author      : Afshin Bayrooti
//
//   Date        : October 2005
//
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/StateVariableCollector.hpp"
#include "edginc/SimpathIInstrument.hpp"

#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"

#include "edginc/ISimPathICollection.hpp"

DRLIB_BEGIN_NAMESPACE

void SimpathICurrencyCollection::initAssetInfo(RM_Assets::MCWriterDataHolder& writerData)
{
    AssetData ad(name, maxDiffDate.getDate());
    AssetData_Currency adc(maxCurveMaturity.getDate(), vOffsets);
    ptrAssetInfo = AssetInfoWritable_CurrencySP(new AssetInfoWritable_Currency(ad, adc, yc, fx) );
    writerData.push_back(ptrAssetInfo.get());
}

void SimpathICurrencyCollection::initSV(IStateVariableGen::IStateGen* pathGen)
{
    // TODO:  return this line:  fxSV = fxGen->getSpotSV(pathGen);
    // which is done outside for the time being -- because of inseparable fx...

    dfSV = dfGen->getSVDiscFactor(pathGen);

    //just in case
    vExpDFSV.clear();

    for(size_t i = 0; i<vExpDFGen.size(); ++i)
        vExpDFSV.push_back(vExpDFGen[i]->getSVExpectedDiscFactor(pathGen));
}

void SimpathICurrencyCollection::updateDiffusedAssetInfo(int iStep, const DateTime& date)
{

    ptrAssetInfo->set_cur_date(date.getDate());
    //FIXME: I don't understand, iStep is on timeLine, so we should convert iStep to an index ?!
	//A: iStep is not on global timeline; it is on timeLine i.e. dates for which diffusion was requested.
    // do nothing beyond the asset maturity
    if (iStep >= dfSV->path().end()-dfSV->path().begin())
        return;

    ptrAssetInfo->m_ssSpotFX = fxPosition == -1 ? 1.0 : (fxSV->path(fxPosition))[iStep];
    ptrAssetInfo->m_ssDiscount =  (dfSV->path())[iStep];

    size_t npoints =  ptrAssetInfo->get_npoints(ptrAssetInfo->get_cur_date());  //vExpDFSV[iStep]->getDates().size();
    for(size_t i=0; i<npoints; ++i)
        ptrAssetInfo->m_ssLnFutureDF[i] = (vExpDFSV[iStep]->path())[i];
}

void SimpathICurrencyCollection::appendSVGen(IStateVariableCollectorSP svCollector)
{
  // TODO: verify this fxGen been appended many times does not create any issues
    svCollector->append(fxGen.get());
    svCollector->append(dfGen.get());
    for(size_t i = 0; i<vExpDFGen.size(); ++i)
        svCollector->append(vExpDFGen[i].get());
}

void SimpathICreditCollection::initAssetInfo(RM_Assets::MCWriterDataHolder& writerData)
{
    AssetData ad(name, maxDiffDate.getDate());
    AssetData_Credit adc(maxCurveMaturity.getDate(), vOffsets);
    ptrAssetInfo= AssetInfoWritable_CreditSP(new AssetInfoWritable_Credit(ad, adc, cds));
    writerData.push_back(ptrAssetInfo.get());
}

// if we use AggregatedGen we don't have sdf or ESDF generators
// so we shoild be initializing only the stuff we have
void SimpathICreditCollection::initSV(IStateVariableGen::IStateGen* pathGen)
{
    if(sdfGen.get())
    sdfSV = sdfGen->getSVSurvivalDiscFactor(pathGen);

    //just in case
    vExpSDFSV.clear();

    for(size_t i = 0; i<vExpSDFGen.size(); ++i)
        if (vExpSDFGen[i].get())
        vExpSDFSV.push_back(vExpSDFGen[i]->getSVExpectedSurvivalDiscFactor(pathGen));

    if (aggSDFGen.get())
        aggSDFSV = aggSDFGen->getIQSVGenAggregatedSurvivalDiscFactorSV(pathGen);
}
#if 0
void SimpathICreditCollection::updateDiffusedAssetInfo(int iStep, const DateTime& date)
{
    ptrAssetInfo->set_cur_date(date.getDate());

    // do nothing beyond the asset maturity
    if (iStep >= sdfSV->end()-sdfSV->begin()) //FIXME : doesn't work with AggregatedGen
        return;

    ptrAssetInfo->m_ssLnNDP = (sdfSV->path())[iStep];

    size_t npoints =  ptrAssetInfo->get_npoints(ptrAssetInfo->get_cur_date());  //vExpDFSV[iStep]->getDates().size();
    for(size_t i=0; i<npoints; ++i)
        ptrAssetInfo->m_ssLnFutureSpread[i] = (vExpSDFSV[iStep]->path())[i];

}
#endif

// With SimpathI we bounded to work with aggregated state variables, so we should be more clever:
// basically we need to reconstruct again all dates we want to print.

void SimpathICreditCollection::updateDiffusedAssetInfo(int iStep, const DateTime& date)
{
    IQMCDiffusibleCreditSpreadBase * cr = aggSDFSV->getAsset();
    
    ptrAssetInfo->valid = false;

    if (date > cr->getDiffusionBound()->getMaxDiffDate())
        return;
    SpotIdx sdfIdx = cr->getSpotIndex( date);
    ASSERT(iStep == sdfIdx); // FIXME: sanity check: we're moveing aling sdfDates for this SV, so iStep should match index
    
    ptrAssetInfo->set_cur_date(date.getDate());
    
/*    if (sdfIdx == SpotIdx::npos) // FIXME: hack, assumes we cannot have ESDF without SDF
        return; // no stop here
    */
    ptrAssetInfo->m_ssNDP = cr->getSurvivalDiscFactor(sdfIdx);

    FwdIdx startIdx = cr->getForwardForwardIndex( date);
    DateTimeArray fwdDates = SimpathIInstrument::calcForwardDates(date, vOffsets, cr->getDiffusionBound()->getMaxCurveMat());

    ptrAssetInfo->m_ssLnExpSDF.clear();
    size_t npoints =  ptrAssetInfo->get_npoints(ptrAssetInfo->get_cur_date());  //vExpDFSV[iStep]->getDates().size();

    QLIB_VERIFY(size_t(npoints) == fwdDates.size(), "ERROR: Number of requested dates != number of calculated dates");
    for(int i=0; i < fwdDates.size(); ++i) {
        FwdIdx idx = cr->getForwardForwardIndex(fwdDates[i]);
        ptrAssetInfo->m_ssLnExpSDF.push_back(cr->getLnExpectedSurvivalDiscFactor( startIdx, idx));
    }
    ptrAssetInfo->valid = true;
}

void SimpathICreditCollection::appendSVGen(IStateVariableCollectorSP svCollector)
{
    if (sdfGen.get())
        svCollector->append(sdfGen.get()); // FIXME: remove when we move to AggregatedSV
    for(size_t i = 0; i != vExpSDFGen.size(); ++i)
        if (vExpSDFGen[i].get())
            svCollector->append(vExpSDFGen[i].get()); // FIXME: remove when we move to AggregatedSV
    if (aggSDFGen.get())
        svCollector->append(aggSDFGen.get());
}

void SimpathIEquityCollection::initAssetInfo(RM_Assets::MCWriterDataHolder& writerData)
{
    AssetData ad(name, maxDiffDate.getDate());
    AssetData_Equity adc/*(eq, "")*/;
    ptrAssetInfo = AssetInfoWritable_EquitySP(new AssetInfoWritable_Equity(ad, adc, eq));
    writerData.push_back(ptrAssetInfo.get());
}

void SimpathIEquityCollection::initSV(IStateVariableGen::IStateGen* pathGen)
{
    eqSV = eqGen->getSpotSV(pathGen);
}

void SimpathIEquityCollection::updateDiffusedAssetInfo(int iStep, const DateTime& date)
{
    ptrAssetInfo->set_cur_date(date.getDate());

    // do nothing beyond the asset maturity
    if (date > maxDiffDate)
    {
        ptrAssetInfo->m_ssSpotEQ = 0.0;
        return;
    }

    // TODO: verify that access here!
    ptrAssetInfo->m_ssSpotEQ = eqPosition == -1 ? 0.0 : (eqSV->path(eqPosition))[iStep];
}

void SimpathIEquityCollection::appendSVGen(IStateVariableCollectorSP svCollector)
{
  // TODO: verify this fxGen been appended many times does not create any issues
    svCollector->append(eqGen.get());
}

//////////////////////////// Energy stuffs here: ///////////////////////////////

void SimpathIEnergyCollection::initAssetInfo(RM_Assets::MCWriterDataHolder& writerData)
{
    AssetData ad(name, maxDiffDate.getDate());
    //AssetData_Energy ade(maxCurveMaturity.getDate(), vDates);
	AssetData_Energy ade(vDates);
    ptrAssetInfo = AssetInfoWritable_EnergySP(new AssetInfoWritable_Energy(ad, ade, enrg));
    writerData.push_back(ptrAssetInfo.get());
}

void SimpathIEnergyCollection::initSV(IStateVariableGen::IStateGen* pathGen)
{
    //just in case
    vExpFpSV.clear();

    for(size_t i = 0; i < vExpFpGen.size(); ++i)
        vExpFpSV.push_back(vExpFpGen[i]->getSVExpEnergyFuture(pathGen));
}

void SimpathIEnergyCollection::updateDiffusedAssetInfo(int iStep, const DateTime& date)
{
    ptrAssetInfo->set_cur_date(date.getDate());

	if (date > maxDiffDate)
		return;
	/*
    size_t npoints =  ptrAssetInfo->get_npoints(ptrAssetInfo->get_cur_date());  //vExpFpSV[iStep]->getDates().size();
    for(size_t i = 0; i < npoints; ++i)
        ptrAssetInfo->m_ssLnFuturePrices[i] = (vExpFpSV[iStep]->path())[i];
		*/
    size_t start = ptrAssetInfo->get_startingpoint(ptrAssetInfo->get_cur_date());  //vExpFpSV[iStep]->getDates().size();
    for(size_t i = start; i < ptrAssetInfo->m_dates.size(); ++i)
    {
        if (ptrAssetInfo->m_dates[i] < maxCurveMaturity.getDate() /* && ptrAssetInfo->m_dates[i] > date.getDate() */)
            ptrAssetInfo->m_ssLnFuturePrices[i] = vExpFpSV[iStep]->getExpFuturePrice(i - start); // *vExpFpSV[iStep] is the same as vExpFpSV[iStep]->path()
        else
            ptrAssetInfo->m_ssLnFuturePrices[i] = 1.e99;    // outside of max curve mat, should not be used
    }

}

void SimpathIEnergyCollection::appendSVGen(IStateVariableCollectorSP svCollector)
{
    // only expected SV gen for energy
    for(size_t i = 0; i < vExpFpGen.size(); ++i)
        svCollector->append(vExpFpGen[i].get());
}

/////////////////////////// Basis /////////////////////////////////////////

void SimpathIBasisCollection::initAssetInfo(RM_Assets::MCWriterDataHolder& writerData)
{
    AssetData ad(name, maxDiffDate.getDate());
    AssetData_Basis adc(maxCurveMaturity.getDate(), vOffsets);
    ptrAssetInfo = AssetInfoWritable_BasisSP(new AssetInfoWritable_Basis(ad, adc, basis) );
    writerData.push_back(ptrAssetInfo.get());
}

void SimpathIBasisCollection::initSV(IStateVariableGen::IStateGen* pathGen)
{
    // TODO:  return this line:  fxSV = fxGen->getSpotSV(pathGen);
    // which is done outside for the time being -- because of inseparable fx...

    //just in case
    vExpBasisFwdSprdSV.clear();

    for(size_t i = 0; i<vExpBasisFwdSprdGen.size(); ++i)
        vExpBasisFwdSprdSV.push_back(vExpBasisFwdSprdGen[i]->getSVExpectedBasisFwdSpread(pathGen));
}

void SimpathIBasisCollection::updateDiffusedAssetInfo(int iStep, const DateTime& date)
{
    ptrAssetInfo->set_cur_date(date.getDate());

    // do nothing beyond the asset maturity
    //if (iStep >= dfSV->end()-dfSV->begin())
    if (date > maxDiffDate)
        return;

    size_t npoints =  ptrAssetInfo->get_npoints(ptrAssetInfo->get_cur_date());  //vExpDFSV[iStep]->getDates().size();
    for(size_t i=0; i<npoints; ++i)
        ptrAssetInfo->m_ssFutureSpread[i] = (vExpBasisFwdSprdSV[iStep]->path())[i];
}

void SimpathIBasisCollection::appendSVGen(IStateVariableCollectorSP svCollector)
{
    // only expected SV gen for basis
    for(size_t i = 0; i<vExpBasisFwdSprdGen.size(); ++i)
        svCollector->append(vExpBasisFwdSprdGen[i].get());
}


DRLIB_END_NAMESPACE

