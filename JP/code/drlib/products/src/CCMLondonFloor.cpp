//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMLondonFloor.cpp
//
//   Description : London floor adjustment Algorithm
//
//   Date        : Dec 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/LondonFloorLossCalculator.hpp" 
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CCM.hpp" // for static CCM::ccmHasStochasticRecoveries
#include "edginc/CCMLondonFloor.hpp"

DRLIB_BEGIN_NAMESPACE

CCMLondonFloor::~CCMLondonFloor(){}

/** Overridden to apply a 'London Floor' adjustment */
CreditMetricsLossCalculatorBase* CCMLondonFloor::createLossCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    // use utility method in parent
    return createCCMLossCalculator(timeline, tranche, cpty);
}

/** Overridden to apply a 'London Floor' adjustment */
CreditMetricsLossCalculatorBase* CCMLondonFloor::createRecoveredNotionalCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    // use utility method in parent
    return createCCMRecoveredNotionalCalculator(timeline, tranche, cpty);
}


/** Indicates whether this model supports stochastic recovery rates
    AND there are any engine parameters for any names specifying so. */
const bool CCMLondonFloor::hasStochasticRecoveries(
    CreditTrancheLossConfigConstSP tranche) const
{
    return CCM::ccmHasStochasticRecoveries(tranche);
}

/** Invoked when Class is 'loaded' */
void CCMLondonFloor::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CCMLondonFloor, clazz);
    SUPERCLASS(CreditMetricsLondonFloor);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

/** Private constructor */
CCMLondonFloor::CCMLondonFloor(const CClassConstSP& clazz) : CreditMetricsLondonFloor(clazz) {}

/** Default constructor */
IObject* CCMLondonFloor::defaultConstructor(){
    return new CCMLondonFloor(CCMLondonFloor::TYPE);
}

CClassConstSP const CCMLondonFloor::TYPE = CClass::registerClassLoadMethod(
    "CCMLondonFloor", typeid(CCMLondonFloor), CCMLondonFloor::load);

bool CCMLondonFloorLoad() {
    return (CCMLondonFloor::TYPE != 0);
}
DRLIB_END_NAMESPACE

