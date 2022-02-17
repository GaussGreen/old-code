//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ScenarioInterface.cpp
//
//   Description : Defines interface to Scenario
//
//   Author      : Andrew J Swain
//
//   Date        : 25 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ScenarioInterface.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

ScenarioInterface::~ScenarioInterface(){}


ScenarioInterface::ScenarioInterface(ScenarioSP    scenario,
                                     IModelSP      model, 
                                     CInstrumentSP inst,
                                     CControlSP    ctrl,
                                     CMarketDataSP market) :
    CObject(TYPE), scenario(scenario), model(model), inst(inst), 
    ctrl(ctrl), market(market) {}

// for reflection
ScenarioInterface::ScenarioInterface(): CObject(TYPE){}

/** Runs 'regression test' */
IObjectSP ScenarioInterface::runTest() const{
    // don't want to write out an input file
    ctrl->switchOffWriteToFile();
    return go();
}

IObjectSP ScenarioInterface::go() const{
    CMarketDataSP marketData(market);
    if (!marketData){
        // review in conjunction with avoiding writing all market data out
        marketData = CMarketDataSP(new MarketData());
    }
    CControlSP control = ctrl->applyCommandLineOptions();
    return model->go(inst, scenario, control, marketData);
}

// EdrAction 
IObjectSP ScenarioInterface::run(){
    return go();
}

/** Invoked when Class is 'loaded' */
void ScenarioInterface::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ScenarioInterface, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultInterface);
    FIELD(scenario, "scenario");
    FIELD(model, "model");
    FIELD(inst, "instrument");
    FIELD(ctrl, "ctrl");
    FIELD(market, "market");
    FIELD_MAKE_OPTIONAL(market);
    // registration for addin function
    Addin::registerObjectMethod("SCENARIO",
                                Addin::RISK,
                                "Calculates price and greeks for a scenario",
                                true,
                                Addin::returnHandle,
                                &ScenarioInterface::run);
}

IObject* ScenarioInterface::defaultInterface(){
    return new ScenarioInterface();
}

CClassConstSP const ScenarioInterface::TYPE = CClass::registerClassLoadMethod(
    "ScenarioInterface", typeid(ScenarioInterface), load);

DRLIB_END_NAMESPACE

