//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Defines interface to Scenario when pricing 
//                 IInstrumentCollection
//
//   Author      : Mark A Robson
//
//   Date        : 10 Feb 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MultiScenarioInterface.hpp"
#include "edginc/Addin.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE
MultiScenarioInterface::~MultiScenarioInterface(){}


MultiScenarioInterface::MultiScenarioInterface(
    ScenarioSP              scenario,
    IModelSP                model, 
    IInstrumentCollectionSP insts,
    CControlSP              ctrl,
    CMarketDataSP           market):
    CObject(TYPE), scenario(scenario), model(model), insts(insts),
    ctrl(ctrl), market(market) {}

// for reflection
MultiScenarioInterface::MultiScenarioInterface(): CObject(TYPE) {}

/** Runs 'regression test' */
IObjectSP MultiScenarioInterface::runTest() const{
    // don't want to write out an input file
    ctrl->switchOffWriteToFile();
    return go();
}

IObjectSP MultiScenarioInterface::go() const{
    CMarketDataSP marketData(market);
    if (!marketData){
        // review in conjunction with avoiding writing all market data out
        marketData.reset(new MarketData());
    }
    CControlSP control = ctrl->applyCommandLineOptions();
    return model->go(insts, scenario, control, marketData);
}

// EdrAction 
IObjectSP MultiScenarioInterface::run(){
    return go();
}

/** Invoked when Class is 'loaded' */
void MultiScenarioInterface::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MultiScenarioInterface, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(scenario, "scenario");
    FIELD(model, "model");
    FIELD(insts, "instruments");
    FIELD(ctrl, "ctrl");
    FIELD(market, "market");
    FIELD_MAKE_OPTIONAL(market);
    // registration for addin function
    Addin::registerObjectMethod(
        "MULTI_SCENARIO",
        Addin::RISK,
        "Calculates price and greeks for a scenario on an "
        "IInstrumentCollection",
        true,
        Addin::returnHandle,
        &MultiScenarioInterface::run);
}

IObject* MultiScenarioInterface::defaultConstructor(){
    return new MultiScenarioInterface();
}

CClassConstSP const MultiScenarioInterface::TYPE = 
CClass::registerClassLoadMethod(
    "MultiScenarioInterface", typeid(MultiScenarioInterface), load);
   

DRLIB_END_NAMESPACE

