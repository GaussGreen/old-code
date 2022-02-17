//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditMetricsLondonFloor.cpp
//
//   Description : Credit Metrics with London floor adjustment Algorithm
//
//   Date        : April 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CreditMetricsLondonFloor.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/LondonFloorLossCalculator.hpp"
#include "edginc/TrancheLossCalculator.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 1e-10

/** Overridden to apply a 'London Floor' adjustment */
ITrancheLossCalculator* CreditMetricsLondonFloor::createLossCalculator(
    const DateTimeArray&           timeline,   /* (I) */
    CreditTrancheLossConfigConstSP tranche,    /* (I) */
    CounterPartyCreditConstSP      cpty) const /* (I) */
{
    // create original
    ITrancheLossCalculatorConstSP cmCalculator(
        CreditMetricsModel::createLossCalculator(timeline, tranche, cpty));

    // and then our adjustment
    double portfolioNotional = tranche->portfolioNotional(); 
    return new LondonFloorLossCalculator(
        cmCalculator, 
        timeline, 
        tranche->getToday(),
        !!cpty,
        londonFloor.getSP(),
        portfolioNotional,
        londonFloorEqStrike * portfolioNotional, 
        londonFloorSeniorStrike * portfolioNotional);
}

CreditMetricsLondonFloor::~CreditMetricsLondonFloor()
{}


/** Invoked by the containing model after the instrument data is fetched
    ie after CInstrument::GetMarket is invoked. Overrides
    CreditMetricsModel in order to get london floor curve */
void CreditMetricsLondonFloor::postFetchMarketData(IModel*            model,
                                                   MarketDataConstSP market)
{
    //now fetch the london floor curve : now may be adjusted according to the
    //index map treatment specification
    ICDSParSpreads::getMarketData(model, market.get(), londonFloor);
    CreditMetricsModel::postFetchMarketData(model, market);
}


/** Called immediately after object constructed */
void CreditMetricsLondonFloor::validatePop2Object() {
    CreditMetricsModel::validatePop2Object();
    try{
        if (londonFloorEqStrike < 0.0)
            throw ModelException ("London Floor equity strike cannot be "
                                  "less than 0");
        if (londonFloorSeniorStrike > 1.0)
            throw ModelException ("London Floor senior strike cannot be "
                                  "greater than 1");
        if (londonFloorEqStrike > londonFloorSeniorStrike)
            throw ModelException ("London Floor equity strike must not be "
                                  "greater than senior strike");
    } catch (exception& e){
        throw ModelException(e, "CreditMetricsLondonFloor::validatePop2Object");
    }
}

/** Invoked when Class is 'loaded' */
void CreditMetricsLondonFloor::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditMetricsLondonFloor, clazz);
    SUPERCLASS(CreditMetricsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
    FIELD(londonFloor,
                 "London floor curve");
    FIELD(londonFloorEqStrike,
                 "Equity strike for the London floor");
    FIELD(londonFloorSeniorStrike,
                 "Senior strike for the London floor");
}

/** Private constructor */
CreditMetricsLondonFloor::CreditMetricsLondonFloor(const CClassConstSP& clazz) : 
    CreditMetricsModel(clazz), londonFloorEqStrike(0.0), londonFloorSeniorStrike(0.0) {}

/** Default constructor */
IObject* CreditMetricsLondonFloor::defaultConstructor(){
    return new CreditMetricsLondonFloor(CreditMetricsLondonFloor::TYPE);
}

CClassConstSP const CreditMetricsLondonFloor::TYPE = 
CClass::registerClassLoadMethod(
    "CreditMetricsLondonFloor", typeid(CreditMetricsLondonFloor), load);

bool CreditMetricsLondonFloorLoad() {
    return (CreditMetricsLondonFloor::TYPE != 0);
}
DRLIB_END_NAMESPACE

