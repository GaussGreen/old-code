//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CreditMetricsRFL.cpp
//
//   Description : Credit Metrics model with RFL
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CreditMetricsRFL.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/CreditMetricsRFLLossCalculator.hpp"
#include "edginc/CreditMetricsRFLFastLossCalculator.hpp"
#include "edginc/RFLParameters.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

///** Called immediately after object constructed */
//void CreditMetricsRFL::validatePop2Object()
//{
//    CreditMetricsModel::validatePop2Object();
//}

/** Destructor */
CreditMetricsRFL::~CreditMetricsRFL() {
}

/** Private constructor */
CreditMetricsRFL::CreditMetricsRFL(const CClassConstSP& clazz) : 
    CreditMetricsModel(clazz)
{}


/** Same as CreditMetricsModel::createLossCalculator() method but
    returns one which is derived from
    CreditMetricsLossCalculatorBase. This method won't make sense
    for certain derived types - or rather its meaning changes.
    The implementation here creates an appropriate CCM loss
    calculator */
CreditMetricsLossCalculatorBase* CreditMetricsRFL::createLossCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CreditMetricsRFLFastLossCalculator(timeline, 
                                                          tranche, 
                                                          false, //recoverNotional
                                                          cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CreditMetricsRFLLossCalculator(timeline, 
                                                      tranche, 
                                                      lossUnitToUse,
                                                      false, //recoverNotional
                                                      cpty);
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "CreditMetricsRFL::createLossCalculatorBase");
    }
}

CreditMetricsLossCalculatorBase* CreditMetricsRFL::createRecoveredNotionalCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    try {
        if (useFastConvolution(tranche)) {
            return new CreditMetricsRFLFastLossCalculator(timeline, 
                                                          tranche, 
                                                          true, //recoverNotional
                                                          cpty);
        } 
        else {
            double lossUnitToUse = calculateLossUnit(tranche);
            return new CreditMetricsRFLLossCalculator(timeline, 
                                                      tranche, 
                                                      lossUnitToUse,
                                                      true, //recoverNotional
                                                      cpty);
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "CreditMetricsModel::"
                             "createRecoveredNotionalCalculatorBase");
    }
}

/** Default constructor */
IObject* CreditMetricsRFL::defaultConstructor(){
    return new CreditMetricsRFL(TYPE);
}

/** Invoked when Class is 'loaded' */
void CreditMetricsRFL::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditMetricsRFL, clazz);
    SUPERCLASS(CreditMetricsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const CreditMetricsRFL::TYPE = CClass::registerClassLoadMethod(
    "CreditMetricsRFL", typeid(CreditMetricsRFL), CreditMetricsRFL::load);

bool CreditMetricsRFLLoad() {
    return (CreditMetricsRFL::TYPE != 0);
}
DRLIB_END_NAMESPACE
