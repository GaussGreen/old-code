//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CCMRFL.cpp
//
//   Description : Composite Copula Model with RFL
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMRFL.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/CCMRFL.hpp"
#include "edginc/CCMRFLLossCalculator.hpp"
#include "edginc/CCMRFLFastLossCalculator.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CCM.hpp" // for static CCM::ccmHasStochasticRecoveries
#include "edginc/RFLParameters.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 3e-15

/* Flag used to get an error message when beta is tweaked outside [-1,1] */
#define BETA_TWEAK_DEBUG 0

/** Useful function to check range of a field */
void inline checkRange(string fieldName, 
                       double fieldValue, 
                       double minValue, 
                       double maxValue) 
{
    if (fieldValue < minValue || fieldValue > maxValue) {
        throw ModelException(fieldName +
                             " (=" + Format::toString(fieldValue) +
                             ") outside [" + Format::toString(minValue) +
                             "," + Format::toString(maxValue) + "]");
    }
}

/** Destructor */
CCMRFL::~CCMRFL() 
{}

/** Private constructor */
CCMRFL::CCMRFL(const CClassConstSP& clazz) : 
    CreditMetricsRFL(clazz)
{}

/** Same as CreditMetricsModel::createLossCalculator() method but
    returns one which is derived from
    CreditMetricsLossCalculatorBase. This method won't make sense
    for certain derived types - or rather its meaning changes.
    The implementation here creates an appropriate CCM loss
    calculator */
CreditMetricsLossCalculatorBase* CCMRFL::createLossCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CCMRFLFastLossCalculator(timeline, 
                                                tranche, 
                                                false, //recoverNotional
                                                cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CCMRFLLossCalculator(timeline, 
                                            tranche, 
                                            lossUnitToUse,
                                            false, //recoverNotional
                                            cpty);
        }
    } 
    catch (exception& e){
        throw ModelException(e, "CCMRFL::createLossCalculatorBase");
    }
}

CreditMetricsLossCalculatorBase* CCMRFL::createRecoveredNotionalCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CCMRFLFastLossCalculator(timeline, 
                                                tranche, 
                                                true, //recoverNotional
                                                cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CCMRFLLossCalculator(timeline, 
                                            tranche, 
                                            lossUnitToUse,
                                            true, //recoverNotional
                                            cpty);
        }
    } 
    catch (exception& e){
        throw ModelException(e, "CCMRFL::createRecoveredNotionalCalculator");
    }
}


/** Indicates whether this model supports stochastic recovery rates
    AND there are any engine parameters for any names specifying so. */
const bool CCMRFL::hasStochasticRecoveries(
    CreditTrancheLossConfigConstSP tranche) const 
{
    return CCM::ccmHasStochasticRecoveries(tranche);
}


/** Default constructor */
IObject* CCMRFL::defaultConstructor() {
    return new CCMRFL(CCMRFL::TYPE);
}

/** Invoked when Class is 'loaded' */
void CCMRFL::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditMetricsRFL, clazz);
    SUPERCLASS(CreditMetricsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
}

CClassConstSP const CCMRFL::TYPE = CClass::registerClassLoadMethod(
    "CCMRFL", typeid(CCMRFL), CCMRFL::load);

bool CCMRFLLoad() {
    return (CCMRFL::TYPE != 0);
}
DRLIB_END_NAMESPACE
