//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMLossCalculator.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CCMLossUnit.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CcmOnlyParameters.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

/** Retrieves CCM stlye engine params from Credit asset. Fails if they are
    not there or are of wrong type */
CcmOnlyParametersConstSP CCMLossCalculatorBase::getCCMEngineParameters(
    SingleCreditAssetConstSP asset){
    CreditEngineParametersConstSP params(
        asset->getEngineParams(CcmOnlyParameters::TYPE));
    return CcmOnlyParametersConstSP::dynamicCast(params); 
}

CCMLossCalculatorBase::~CCMLossCalculatorBase(){}

/** Constructor - takes in full timeline to allow for optimisations. If
    computeCondCurve if false then lossConditional in loss() below will be
    set to zero. */
CCMLossCalculatorBase::CCMLossCalculatorBase(
        const DateTimeArray&           timeline, /* (I) */
        CreditTrancheLossConfigConstSP tranche,  /* (I) */
        CounterPartyCreditConstSP      cpty):    /* (I) */
    // parent handles survivalProb and counterPartyProb
    CreditMetricsLossCalculatorBase(timeline, tranche, cpty)
{
    try {
        // compute probabilities from 'senior curve' (CCM specific)
        int numNames = tranche->numInnerLossConfigs(); // for ease
        // reserve some space
        floorSurvivalProb.resize(timeline.size());
        for (int t = 0; t < timeline.size(); t++) {
            floorSurvivalProb[t].resize(numNames);
        }
        const DateTime& today = tranche->getToday();
        for (int i = 0; i < numNames; i++) {
            if (tranche->nameDefaulted(i)) {
                for (int t = 0; t < timeline.size(); t++){
                    floorSurvivalProb[t][i] = 0.0; // override any input
                }
            } 
            else {
#if 0
                // this may be the right thing to do (for forward
                // starting tranches)
                const DateTime& protectionStart = 
                    tranche->nameProtectionStartDate(i, today);
#else
                const DateTime& protectionStart = today;
#endif
                const DateTime& protectionEnd =
                    today.max(tranche->nameProtectionEndDate(i,
                                                             timeline.back()));
                SingleCreditAssetConstSP asset(tranche->nameAsset(i));
                CcmOnlyParametersConstSP ccmParams(getCCMEngineParameters(asset));
                for (int t = 0; t < timeline.size(); t++) {
                    floorSurvivalProb[t][i] =
                        ccmParams->dependenceSurvivalProba(
                            protectionStart, timeline[t].min(protectionEnd));
                }
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "CCMLossCalculatorBase::CCMLossCalculatorBase");
    }
}


DRLIB_END_NAMESPACE
