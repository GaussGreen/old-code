//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditDiffuse.cpp
//
//   Description : Base class for SRM credit diffusion
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditDiffuse.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/MemoryProfiler.hpp"
#include "edginc/Format.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/LinearInterpolator.hpp"

#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** constructor */
SRMCreditDiffuse::SRMCreditDiffuse(IQMCDiffusibleInterestRateSP pRatesDiffuse) :
    // meaningful initialization
    irAsset(pRatesDiffuse),
    sigmaFX(pRatesDiffuse.get() ? pRatesDiffuse->getSigmaFX() : NULL),
    // default initialization with 'bad' values
    randomIndex(-999),
    crFxCorr(-999.),
    diffBound(pRatesDiffuse->getDiffusionBound()) // dependent asset
{}

/** destructor */
SRMCreditDiffuse::~SRMCreditDiffuse() {
}


/** add all dates that will be needed later */
void SRMCreditDiffuse::addAggregatedDates(
    DateTimeArrayConstSP _sdfDates,
    DateTimeArrayConstSP _esdfRequestedDates,
    DateTimeArrayConstSP _esdfForwardDates,
    const DateTime& maxDiffDate,
    const DateTime& maxCurveMat) 
{
    // SRMCredit is an asset that worth the effort of saving memory and using AggregatedState variables 
    // (because in Sampras we have 5,000 of them per run)
    
    // QMCCreditDiffuse::addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates, maxDiffDate, maxCurveMat);
    
    getTimeLogic()->addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates);
    getTimeLogic()->addAggregatedDates(DateTimeArrayConstSP(), DateTimeArrayConstSP(), _esdfRequestedDates);
    
    irAsset->addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates, maxDiffDate, maxCurveMat);
    
    getDiffusionBound()->updateMaxDiffDate(maxDiffDate);
    getDiffusionBound()->updateCurveMat(maxCurveMat);
}

/** return underlying IR asset */
IQMCDiffusibleInterestRateSP SRMCreditDiffuse::getUnderlyingIRAsset() {
    return irAsset;
}

/** returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
QMCHelperBoundedDiffusion* SRMCreditDiffuse::getDiffusionBound()  {
    return &diffBound;
}

double SRMCreditDiffuse::getWholeTimelineLogSurvProb(size_t i) {
    if (wholeTimelineLogSurvProb.size()<=i)
    {
        QLIB_VERIFY(isWholeTimelineSurvivalRequested(), 
            "this function cannot be called unless the request parameter for the whole timeline surv prob was set.");
        QLIB_VERIFY(!wholeTimelineLogSurvProb.empty(),
            "the particular diffusion model of Credit did not initialize the array to hold the whole timeline surv prob.");
        QLIB_VERIFY(wholeTimelineLogSurvProb.size()>i,
            "the requested point is outside of the generated array of the whole timeline surv prob.");


    }
    return wholeTimelineLogSurvProb.at(i);
}


// For each forward-forward date find its index in IR asset
void SRMCreditDiffuse::calcRemapToIRAssetIdx(const DateTimeArray& esdfForwardDates)
{
        // aux mapping to map between CR fwd dates to IR fwd dates
    this->fwdIdx2IRedfIdx.resize(esdfForwardDates.size());
    for(int i=0; i < esdfForwardDates.size(); ++i)
        fwdIdx2IRedfIdx[i] = irAsset->getForwardForwardIndex(esdfForwardDates[i]);

}

/** returns IR forward index */
FwdIdx SRMCreditDiffuse::getIRfwdIdx(FwdIdx idx) const {
    return fwdIdx2IRedfIdx[idx];
}

DRLIB_END_NAMESPACE
