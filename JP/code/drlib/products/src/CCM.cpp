//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CCM.cpp
//
//   Description : Simple closed form model for CCM
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CcmOnlyParameters.hpp"
#include "edginc/CCM.hpp"

DRLIB_BEGIN_NAMESPACE

CCM::CCM(const CClassConstSP& clazz): CreditMetricsModel(clazz)
{}

/** Destructor */
CCM::~CCM() 
{}

/** Same as CreditMetricsModel::createLossCalculator() method but
    returns one which is derived from
    CreditMetricsLossCalculatorBase. This method won't make sense
    for certain derived types - or rather its meaning changes.
    The implementation here creates an appropriate CCM loss
    calculator */
CreditMetricsLossCalculatorBase* CCM::createLossCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    // use utility method in parent
    return createCCMLossCalculator(timeline, tranche, cpty);
}


CreditMetricsLossCalculatorBase* CCM::createRecoveredNotionalCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    // use utility method in parent
    return createCCMRecoveredNotionalCalculator(timeline, tranche, cpty);
}


/** Indicates whether this model supports stochastic recovery rates
    AND there are any engine parameters for any names specifying so. */
const bool CCM::hasStochasticRecoveries(
    CreditTrancheLossConfigConstSP tranche) const 
{
    return CCM::ccmHasStochasticRecoveries(tranche);
}


/** Static method that indicates whether the tranche contains names with
    engine parameters specifying stochastic recovery rates */
const bool CCM::ccmHasStochasticRecoveries(
    CreditTrancheLossConfigConstSP tranche) 
{
    bool hasStochastic = false;
    int numNames = tranche->numInnerLossConfigs();

    for (int i=0; i < numNames && !hasStochastic; ++i) {
        ICreditLossConfigConstSP innerLossCfg = tranche->getInnerLossConfig(i);
        // first loss config with stochastic parameters causes the loop to end
        CreditEngineParametersConstSP engineParams = 
            innerLossCfg->getEngineParams(CcmOnlyParameters::TYPE);

        CcmOnlyParametersConstSP ccmParams = 
            CcmOnlyParametersConstSP::dynamicCast(engineParams);

        hasStochastic = ccmParams->hasStochasticRecovery();
    }
    return hasStochastic;
}


class CCMHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CCM, clazz);
        SUPERCLASS(CreditMetricsModel);
        EMPTY_SHELL_METHOD(defaultCCM);
    }
    
    static IObject* defaultCCM(){
        return new CCM(CCM::TYPE);
    }
};

CClassConstSP const CCM::TYPE = CClass::registerClassLoadMethod(
    "CCM", typeid(CCM), CCMHelper::load);

bool CCMLoad() {
    return (CCM::TYPE != 0);
}

DRLIB_END_NAMESPACE
