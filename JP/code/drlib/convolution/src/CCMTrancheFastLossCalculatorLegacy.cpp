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
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CCMTrancheFastLossCalculatorLegacy.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CCMTrancheCalculatorLegacy.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

/** Builds calculator skipping all variant data ie data that changes
    across time */
CCMTrancheFastLossCalculatorLegacy::CCMTrancheFastLossCalculatorLegacy(
    bool                           doCreditMetrics,
    bool                           doRFL,
    const bool                     recoverNotional,
    CreditTrancheLossConfigConstSP tranche,
    CounterPartyCreditConstSP      counterParty)
{
    pastLoss = tranche->portfolioLoss();

    if (doRFL) {
        CCMTrancheCalculatorLegacy::createBasketInfoRFL(tranche, basketInfo);
        cpty = CCMTrancheCalculatorLegacy::createCptyInfoRFL(counterParty);
    }
    else {
        CCMTrancheCalculatorLegacy::createBasketInfo(tranche, basketInfo);
        cpty = CCMTrancheCalculatorLegacy::createCptyInfo(counterParty);
    }

    CCMTrancheCalculatorLegacy::populateBasketInfoCM(tranche, recoverNotional, basketInfo);
    CCMTrancheCalculatorLegacy::populateCptyInfoCM(counterParty, cpty);

    if (!doCreditMetrics){
        // do additional stuff for CCM
        CCMTrancheCalculatorLegacy::populateBasketInfoCCM(tranche, basketInfo);
        CCMTrancheCalculatorLegacy::populateCptyInfoCCM(counterParty, *cpty);
    }
    if (doRFL) {
        // do additional stuff for RFL
        CCMTrancheCalculatorLegacy::populateBasketInfoRFL(tranche, basketInfo);
        CCMTrancheCalculatorLegacy::populateCptyInfoRFL(counterParty, *cpty);    
    }
}


/** Same as comments as setupCM above */
void CCMTrancheFastLossCalculatorLegacy::setupCM(
    const DoubleArray&        survivalProb,      /* name survival proba */
    double                    counterPartyProb,  /* name survival proba */
    const DoubleArray&        betaOverride)      // optional per name override
{
    bool useOverride = !betaOverride.empty();
    for (int i = 0; i < basketInfo.size(); i++){
        basketInfo[i]->survival = survivalProb[i];
        if (useOverride){
            basketInfo[i]->beta = betaOverride[i];
        }
    }
    if (cpty.get()){
        cpty->survival      = counterPartyProb;
    }
}

/** Same as comments as setupCM above */
void CCMTrancheFastLossCalculatorLegacy::setupCCM(
    const DoubleArray&        survivalProb,      /* name survival proba */
    const DoubleArray&        floorSurvivalProb, /* from senior curve */
    double                    counterPartyProb, /* name survival proba */
    const DoubleArray&        betaOverride)  // optional per name override
{
    bool useOverride = !betaOverride.empty();
    // just copy over probabilities
    for (int i = 0; i < basketInfo.size(); i++){
        basketInfo[i]->survival = survivalProb[i];
        basketInfo[i]->pdep = floorSurvivalProb[i];
        if (useOverride){
            basketInfo[i]->beta = betaOverride[i];
        }
    }
    if (cpty.get()){
        cpty->survival      = counterPartyProb;
    }
}

/** creates a copy */
CCMTrancheFastLossCalculatorLegacy::CCMTrancheFastLossCalculatorLegacy(
    CCMTrancheFastLossCalculatorLegacy* copy)
{
    // this is a copy and paste of the non fast version - needs 
    // complete refactoring
    basketInfo = *CCMConvolution::NameParamArraySP((
        dynamic_cast<CCMConvolution::NameParamArray*>(copy->basketInfo.clone())));

    if (copy->cpty.get() != 0) {
        cpty    = CCMConvolution::NameParamSP( copy->cpty.clone());
    }
    pastLoss = copy->pastLoss;
}

/** 
 * function calculating expected loss given strikes 
 * can be invoked directly or as a callback by london floor adjustment function
 */
void CCMTrancheFastLossCalculatorLegacy::loss(
        double        K1,    /* (I) lower strike      */
        double        K2,    /* (I) upper strike      */
        double       &L,     /* (O) tranche loss amt  */
        double       &Lcond) /* (O) tranche loss amt cond on cpty surviving */
        const 
{
    CCMConvolution conv(true, //use fast
                        basketInfo,
                        0, 
                        0);

    //CCMConvolution conv(true, basketInfo,cpty.get(),0,0);

    conv.calcFastTrancheLoss(cpty, pastLoss, K1, K2);
    L = conv.getExpLoss();
    Lcond = conv.getExpLossCond();
}


DRLIB_END_NAMESPACE
